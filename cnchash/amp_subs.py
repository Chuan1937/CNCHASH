"""
S/P amplitude ratio subroutines for focal mechanism determination.

Implements FOCALAMP_MC algorithm (Hardebeck and Shearer, 2003) for finding
focal mechanisms using both P-wave polarities and S/P amplitude ratios.

Optimized with numba JIT compilation.
"""

import math

import numpy as np
from numba import njit

from .core import get_rotation_grid
from .utils import (
    DEG_TO_RAD,
    PI,
    fp_coord_vectors_to_angles,
)

# Table size for amplitude lookup
NTAB = 180

# Cache for amplitude tables
_AMP_TABLE_CACHE = {}


@njit(cache=True)
def _setup_amplitude_tables():
    """
    Set up lookup tables for P and S amplitudes.

    Returns
    -------
    tuple
        (amptable, phitable, thetable)
        - amptable: P and S amplitudes for each (theta, phi)
        - phitable: phi angle lookup
        - thetable: theta angle lookup
    """
    amptable = np.zeros((2, NTAB, 2 * NTAB), dtype=np.float64)
    phitable = np.zeros((2 * NTAB + 1, 2 * NTAB + 1), dtype=np.float64)
    thetable = np.zeros(2 * NTAB + 1, dtype=np.float64)

    astep = 1.0 / float(NTAB)

    # Setup theta and phi tables
    for i in range(2 * NTAB + 1):
        bbb3 = -1.0 + float(i) * astep
        thetable[i] = math.acos(bbb3)
        for j in range(2 * NTAB + 1):
            bbb1 = -1.0 + float(j) * astep
            phitable[i, j] = math.atan2(bbb3, bbb1)
            if phitable[i, j] < 0.0:
                phitable[i, j] += 2.0 * PI

    # Setup amplitude tables
    # amptable[0, :, :] = P amplitude = |sin(2*theta)*cos(phi)|
    # amptable[1, :, :] = S amplitude = sqrt(s1^2 + s2^2)
    #   where s1 = cos(2*theta)*cos(phi)
    #         s2 = -cos(theta)*sin(phi)
    for i in range(2 * NTAB):
        phi = float(i) * PI * astep
        for j in range(NTAB):
            theta = float(j) * PI * astep
            amptable[0, j, i] = abs(math.sin(2.0 * theta) * math.cos(phi))
            s1 = math.cos(2.0 * theta) * math.cos(phi)
            s2 = -math.cos(theta) * math.sin(phi)
            amptable[1, j, i] = math.sqrt(s1 * s1 + s2 * s2)

    return amptable, phitable, thetable


def get_amplitude_tables():
    """
    Get or create the amplitude lookup tables.

    Results are cached to avoid recomputation.

    Returns
    -------
    tuple
        (amptable, phitable, thetable)
    """
    if "tables" not in _AMP_TABLE_CACHE:
        _AMP_TABLE_CACHE["tables"] = _setup_amplitude_tables()
    return _AMP_TABLE_CACHE["tables"]


@njit(cache=True)
def _compute_sp_ratio(p_amp, s_amp):
    """
    Compute S/P amplitude ratio (log10 scale).

    Uses the formula from Hardebeck and Shearer (2003):
    sp_ratio = log10(4.9 * s_amp / p_amp)

    Parameters
    ----------
    p_amp : float
        P-wave amplitude
    s_amp : float
        S-wave amplitude

    Returns
    -------
    float
        log10(S/P ratio), clipped to reasonable range
    """
    if p_amp == 0.0:
        return 4.0  # Upper limit
    elif s_amp == 0.0:
        return -2.0  # Lower limit
    else:
        return math.log10(4.9 * s_amp / p_amp)


@njit(cache=True)
def _find_best_mechanisms_with_amp(
    p_azi_mc,
    p_the_mc,
    sp_amp,
    p_pol,
    b1,
    b2,
    b3,
    nrot,
    npsta,
    nmc,
    nextra,
    ntotal,
    qextra,
    qtotal,
    amptable,
    phitable,
    thetable,
):
    """
    Find acceptable focal mechanisms using both polarities and S/P ratios.

    Parameters
    ----------
    p_azi_mc : ndarray, shape (npsta, nmc)
        Azimuths for each observation and trial
    p_the_mc : ndarray, shape (npsta, nmc)
        Takeoff angles for each observation and trial
    sp_amp : ndarray, shape (npsta,)
        S/P amplitude ratios (log10 scale)
    p_pol : ndarray, shape (npsta,)
        Polarity observations (+1, -1, or 0 for no polarity)
    b1, b2, b3 : ndarray
        Rotation arrays
    nrot : int
        Number of rotations
    npsta : int
        Number of observations
    nmc : int
        Number of Monte Carlo trials
    nextra : int
        Additional polarity misfits allowed above minimum
    ntotal : int
        Total allowed polarity misfits
    qextra : float
        Additional amplitude misfit allowed above minimum
    qtotal : float
        Total allowed amplitude misfit
    amptable, phitable, thetable : ndarray
        Amplitude lookup tables

    Returns
    -------
    ndarray, shape (nrot,)
        Boolean array indicating which rotations meet criteria
    """
    irotgood = np.zeros(nrot, dtype=np.int32)

    astep = 1.0 / float(NTAB)
    rad = DEG_TO_RAD

    for im in range(nmc):
        # Convert to Cartesian
        p_a1 = np.empty(npsta, dtype=np.float64)
        p_a2 = np.empty(npsta, dtype=np.float64)
        p_a3 = np.empty(npsta, dtype=np.float64)

        for i in range(npsta):
            theta = p_the_mc[i, im]
            phi = p_azi_mc[i, im]
            theta_rad = theta * rad
            phi_rad = phi * rad
            p_a3[i] = -1.0 * math.cos(theta_rad)
            p_a1[i] = math.sin(theta_rad) * math.cos(phi_rad)
            p_a2[i] = math.sin(theta_rad) * math.sin(phi_rad)

        # Find misfit for each rotation
        nmis = np.zeros(nrot, dtype=np.int32)
        qmis = np.zeros(nrot, dtype=np.float64)

        nmis0min = 100000
        qmis0min = 1.0e5

        for irot in range(nrot):
            for ista in range(npsta):
                # Project ray direction onto fault frame
                p_b1 = (
                    b1[0, irot] * p_a1[ista] + b1[1, irot] * p_a2[ista] + b1[2, irot] * p_a3[ista]
                )
                p_b3 = (
                    b3[0, irot] * p_a1[ista] + b3[1, irot] * p_a2[ista] + b3[2, irot] * p_a3[ista]
                )

                # S/P amplitude ratio misfit
                if sp_amp[ista] != 0.0:
                    # Project onto plane perpendicular to fault normal
                    p_proj1 = p_a1[ista] - p_b3 * b3[0, irot]
                    p_proj2 = p_a2[ista] - p_b3 * b3[1, irot]
                    p_proj3 = p_a3[ista] - p_b3 * b3[2, irot]
                    plen = math.sqrt(p_proj1 * p_proj1 + p_proj2 * p_proj2 + p_proj3 * p_proj3)

                    if plen > 0:
                        p_proj1 /= plen
                        p_proj2 /= plen
                        p_proj3 /= plen

                    pp_b1 = b1[0, irot] * p_proj1 + b1[1, irot] * p_proj2 + b1[2, irot] * p_proj3
                    pp_b2 = b2[0, irot] * p_proj1 + b2[1, irot] * p_proj2 + b2[2, irot] * p_proj3

                    # Look up theta and phi from tables
                    i_idx = int(round((p_b3 + 1.0) / astep)) + 1
                    if i_idx < 1:
                        i_idx = 1
                    if i_idx > 2 * NTAB + 1:
                        i_idx = 2 * NTAB + 1
                    theta = thetable[i_idx - 1]

                    i_idx = int(round((pp_b2 + 1.0) / astep)) + 1
                    j_idx = int(round((pp_b1 + 1.0) / astep)) + 1
                    if i_idx < 1:
                        i_idx = 1
                    if i_idx > 2 * NTAB + 1:
                        i_idx = 2 * NTAB + 1
                    if j_idx < 1:
                        j_idx = 1
                    if j_idx > 2 * NTAB + 1:
                        j_idx = 2 * NTAB + 1
                    phi = phitable[i_idx - 1, j_idx - 1]

                    # Look up amplitude
                    i_idx = int(round(phi / (PI * astep))) + 1
                    if i_idx > 2 * NTAB:
                        i_idx = 1
                    j_idx = int(round(theta / (PI * astep))) + 1
                    if j_idx > NTAB:
                        j_idx = NTAB

                    p_amp = amptable[0, j_idx - 1, i_idx - 1]
                    s_amp = amptable[1, j_idx - 1, i_idx - 1]

                    sp_ratio = _compute_sp_ratio(p_amp, s_amp)
                    qmis[irot] += abs(sp_amp[ista] - sp_ratio)

                # Polarity misfit
                if p_pol[ista] != 0:
                    prod = p_b1 * p_b3
                    ipol = -1
                    if prod > 0.0:
                        ipol = 1
                    if ipol != p_pol[ista]:
                        nmis[irot] += 1

            if nmis[irot] < nmis0min:
                nmis0min = nmis[irot]
            if qmis[irot] < qmis0min:
                qmis0min = qmis[irot]

        # Set acceptance criteria
        nmis0max = ntotal
        if nmis0max < nmis0min + nextra:
            nmis0max = nmis0min + nextra
        qmis0max = qtotal
        if qmis0max < qmis0min + qextra:
            qmis0max = qmis0min + qextra

        # Find rotations meeting criteria
        nadd = 0
        for irot in range(nrot):
            if nmis[irot] <= nmis0max and qmis[irot] <= qmis0max:
                irotgood[irot] = 1
                nadd += 1

        # If none meet criteria, loosen amplitude criteria
        if nadd == 0:
            qmis0min = 1.0e5
            for irot in range(nrot):
                if nmis[irot] <= nmis0max and qmis[irot] < qmis0min:
                    qmis0min = qmis[irot]
            qmis0max = qtotal
            if qmis0max < qmis0min + qextra:
                qmis0max = qmis0min + qextra

            for irot in range(nrot):
                if nmis[irot] <= nmis0max and qmis[irot] <= qmis0max:
                    irotgood[irot] = 1

    return irotgood


def focalamp_mc(
    p_azi_mc, p_the_mc, sp_amp, p_pol, npsta, nmc, dang, maxout, nextra, ntotal, qextra, qtotal
):
    """
    Perform grid search to find acceptable focal mechanisms using
    both P-wave polarities and S/P amplitude ratios.

    This implements the FOCALAMP_MC algorithm from HASH v1.2.

    Parameters
    ----------
    p_azi_mc : ndarray, shape (npsta, nmc)
        Azimuth to station from event (degrees, East of North)
    p_the_mc : ndarray, shape (npsta, nmc)
        Takeoff angle (degrees, from vertical, up=0, <90 upgoing, >90 downgoing)
    sp_amp : ndarray, shape (npsta,)
        S/P amplitude ratios (log10 scale, 0 if no ratio)
    p_pol : ndarray, shape (npsta,)
        Polarity observations: 1=up, -1=down, 0=no polarity
    npsta : int
        Number of observations
    nmc : int
        Number of Monte Carlo trials
    dang : float
        Desired angle spacing for grid search (degrees)
    maxout : int
        Maximum number of fault planes to return
    nextra : int
        Additional polarity misfits allowed above minimum
    ntotal : int
        Total number of allowed polarity misfits
    qextra : float
        Additional amplitude misfit allowed above minimum
    qtotal : float
        Total allowed amplitude misfit

    Returns
    -------
    dict
        Dictionary containing:
        - 'nf': number of fault planes found
        - 'strike': array of strike angles
        - 'dip': array of dip angles
        - 'rake': array of rake angles
        - 'faults': array of fault normal vectors (3, nf)
        - 'slips': array of slip vectors (3, nf)
    """
    # Get rotation grid
    b1, b2, b3, nrot = get_rotation_grid(dang)

    # Get amplitude tables
    amptable, phitable, thetable = get_amplitude_tables()

    # Find good rotations
    irotgood = _find_best_mechanisms_with_amp(
        p_azi_mc,
        p_the_mc,
        sp_amp,
        p_pol,
        b1,
        b2,
        b3,
        nrot,
        npsta,
        nmc,
        nextra,
        ntotal,
        qextra,
        qtotal,
        amptable,
        phitable,
        thetable,
    )

    # Count good rotations
    nfault = np.sum(irotgood)

    # Get indices of good rotations
    good_indices = np.where(irotgood > 0)[0]

    # Select output solutions
    nf = min(nfault, maxout)

    if nfault <= maxout:
        selected_indices = good_indices
    else:
        np.random.shuffle(good_indices)
        selected_indices = good_indices[:maxout]
        nf = len(selected_indices)

    # Convert to strike, dip, rake
    strike = np.zeros(nf, dtype=np.float64)
    dip = np.zeros(nf, dtype=np.float64)
    rake = np.zeros(nf, dtype=np.float64)
    faults = np.zeros((3, nf), dtype=np.float64)
    slips = np.zeros((3, nf), dtype=np.float64)

    for i, idx in enumerate(selected_indices):
        faultnorm = b3[:, idx].copy()
        slip = b1[:, idx].copy()

        faults[0, i] = faultnorm[0]
        faults[1, i] = faultnorm[1]
        faults[2, i] = faultnorm[2]

        slips[0, i] = slip[0]
        slips[1, i] = slip[1]
        slips[2, i] = slip[2]

        s, d, r = fp_coord_vectors_to_angles(faultnorm, slip)
        strike[i] = s
        dip[i] = d
        rake[i] = r

    return {
        "nf": nf,
        "strike": strike,
        "dip": dip,
        "rake": rake,
        "faults": faults,
        "slips": slips,
    }


@njit(cache=True)
def get_misf_amp_numba(npol, p_azi, p_the, sp_ratio, p_pol, strike, dip, rake):
    """
    Calculate misfit for a given mechanism using polarities and S/P ratios.

    Parameters
    ----------
    npol : int
        Number of observations
    p_azi : ndarray
        Azimuths (degrees)
    p_the : ndarray
        Takeoff angles (degrees)
    sp_ratio : ndarray
        S/P amplitude ratios (log10 scale, 0 if no ratio)
    p_pol : ndarray
        Polarity observations (+1, -1, or 0)
    strike, dip, rake : float
        Mechanism parameters (degrees)

    Returns
    -------
    tuple
        (mfrac, mavg, stdr) - polarity misfit fraction, average S/P misfit,
        and station distribution ratio
    """
    rad = DEG_TO_RAD

    # Convert strike, dip, rake to radians
    s_rad = strike * rad
    d_rad = dip * rad
    r_rad = rake * rad

    # Moment tensor components
    sin_d = math.sin(d_rad)
    cos_d = math.cos(d_rad)
    sin_r = math.sin(r_rad)
    cos_r = math.cos(r_rad)
    sin_s = math.sin(s_rad)
    cos_s = math.cos(s_rad)
    sin_2s = math.sin(2.0 * s_rad)
    sin_2d = math.sin(2.0 * d_rad)
    cos_2d = math.cos(2.0 * d_rad)

    M11 = -sin_d * cos_r * sin_2s - sin_2d * sin_r * sin_s * sin_s
    M22 = sin_d * cos_r * sin_2s - sin_2d * sin_r * cos_s * cos_s
    M33 = sin_2d * sin_r
    M12 = sin_d * cos_r * math.cos(2.0 * s_rad) + 0.5 * sin_2d * sin_r * math.sin(2.0 * s_rad)
    M13 = -cos_d * cos_r * cos_s - cos_2d * sin_r * sin_s
    M23 = -cos_d * cos_r * sin_s + cos_2d * sin_r * cos_s

    # Get fault normal and slip vectors
    fn3 = -cos_d
    fn1 = -sin_d * sin_s
    fn2 = sin_d * cos_s

    sl1 = cos_r * cos_s + cos_d * sin_r * sin_s
    sl2 = cos_r * sin_s - cos_d * sin_r * cos_s
    sl3 = -sin_r * sin_d

    # Auxiliary plane (b2 = cross product of fn and sl)
    b21 = fn2 * sl3 - fn3 * sl2
    b22 = fn3 * sl1 - fn1 * sl3
    b23 = fn1 * sl2 - fn2 * sl1

    mfrac = 0.0
    qcount = 0.0
    scount = 0.0
    mavg = 0.0
    acount = 0.0

    for k in range(npol):
        # Convert to Cartesian
        theta = p_the[k] * rad
        phi = p_azi[k] * rad

        p_a1 = math.sin(theta) * math.cos(phi)
        p_a2 = math.sin(theta) * math.sin(phi)
        p_a3 = -math.cos(theta)

        # Project onto mechanism
        _p_b1 = sl1 * p_a1 + sl2 * p_a2 + sl3 * p_a3  # noqa: F841 - kept for Fortran consistency
        p_b3 = fn1 * p_a1 + fn2 * p_a2 + fn3 * p_a3

        # Project onto plane perpendicular to fault normal
        p_proj1 = p_a1 - p_b3 * fn1
        p_proj2 = p_a2 - p_b3 * fn2
        p_proj3 = p_a3 - p_b3 * fn3

        plen = math.sqrt(p_proj1**2 + p_proj2**2 + p_proj3**2)
        if plen > 0:
            p_proj1 /= plen
            p_proj2 /= plen
            p_proj3 /= plen

        pp_b1 = sl1 * p_proj1 + sl2 * p_proj2 + sl3 * p_proj3
        pp_b2 = b21 * p_proj1 + b22 * p_proj2 + b23 * p_proj3

        phi_ang = math.atan2(pp_b2, pp_b1)
        theta_ang = math.acos(max(-1.0, min(1.0, p_b3)))

        p_amp = abs(math.sin(2.0 * theta_ang) * math.cos(phi_ang))
        wt = math.sqrt(p_amp)

        # Polarity misfit
        if p_pol[k] != 0:
            a1 = math.sin(theta) * math.cos(phi)
            a2 = math.sin(theta) * math.sin(phi)
            a3 = -math.cos(theta)

            b1_val = M11 * a1 + M12 * a2 + M13 * a3
            b2_val = M12 * a1 + M22 * a2 + M23 * a3
            b3_val = M13 * a1 + M23 * a2 + M33 * a3

            dot_val = a1 * b1_val + a2 * b2_val + a3 * b3_val

            if dot_val < 0:
                pol = -1
            else:
                pol = 1

            if pol * p_pol[k] < 0:
                mfrac += wt

            qcount += wt
            scount += 1.0

        # S/P amplitude ratio misfit
        if sp_ratio[k] != 0.0:
            s1 = math.cos(2.0 * theta_ang) * math.cos(phi_ang)
            s2 = -math.cos(theta_ang) * math.sin(phi_ang)
            s_amp = math.sqrt(s1 * s1 + s2 * s2)

            if p_amp > 0:
                sp_rat = math.log10(4.9 * s_amp / p_amp)
                mavg += abs(sp_ratio[k] - sp_rat)

            acount += 1.0
            scount += 1.0

    if qcount > 0:
        mfrac /= qcount
    if acount > 0:
        mavg /= acount
    if scount > 0:
        stdr = qcount / scount
    else:
        stdr = 0.0

    return mfrac, mavg, stdr


def get_misf_amp(npol, p_azi, p_the, sp_ratio, p_pol, strike, dip, rake):
    """
    Find the percent of misfit polarities and average S/P misfit for a given mechanism.

    Parameters
    ----------
    npol : int
        Number of observations
    p_azi : ndarray, shape (npol,)
        Azimuths (degrees, East of North)
    p_the : ndarray, shape (npol,)
        Takeoff angles (degrees)
    sp_ratio : ndarray, shape (npol,)
        S/P amplitude ratios (log10 scale, 0 if no ratio)
    p_pol : ndarray, shape (npol,)
        Polarity observations: 1=up, -1=down, 0=no polarity
    strike, dip, rake : float
        Mechanism parameters (degrees)

    Returns
    -------
    tuple
        (mfrac, mavg, stdr) - weighted fraction misfit polarities,
        average S/P misfit (log10), and station distribution ratio
    """
    return get_misf_amp_numba(npol, p_azi, p_the, sp_ratio, p_pol, strike, dip, rake)


# Export functions
__all__ = [
    "focalamp_mc",
    "get_misf_amp",
    "get_amplitude_tables",
]
