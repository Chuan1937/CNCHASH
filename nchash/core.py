"""
Core HASH algorithms for focal mechanism determination.

Implements the grid search algorithm (FOCALMC) and related functions
for finding acceptable focal mechanisms from polarity data.
Optimized with numba JIT compilation.
"""

import numpy as np
from numba import njit, float64, int32, int8, boolean
import math

from .utils import (
    cross_product,
    cross_product_array,
    to_cartesian_array,
    fp_coord_vectors_to_angles,
    normal_distribution_random_numba,
    normalize_vectors_array,
    dot_product,
    DEG_TO_RAD,
    RAD_TO_DEG,
    PI,
)

# Constants for grid rotation
# These define the grid resolution for searching focal mechanisms
# Higher values = finer grid = more accurate but slower
_GRID_CACHE = {}
_GRID_DANG_OLD = None


@njit(cache=True)
def _setup_rotation_grid(dang):
    """
    Set up the rotation grid for focal mechanism search.

    This exactly matches the Fortran FOCALMC grid setup (lines 60-110).

    Parameters
    ----------
    dang : float
        Grid angle spacing in degrees

    Returns
    -------
    tuple
        (b1, b2, b3, nrot) where b1, b2, b3 are rotation arrays
        and nrot is the number of rotations
    """
    # Maximum number of rotations (adjust based on dang)
    # Match Fortran: ithe=0,int(90.1/dang) and izeta=0,int(179.9/dang)
    max_the = int(90.1 / dang)
    max_zeta = int(179.9 / dang)
    max_rot = (max_the + 1) * (360 + 1) * (max_zeta + 1)

    # Pre-allocate arrays
    b1 = np.zeros((3, max_rot), dtype=np.float64)
    b2 = np.zeros((3, max_rot), dtype=np.float64)
    b3 = np.zeros((3, max_rot), dtype=np.float64)

    irot = 0

    for ithe in range(max_the + 1):
        the = float(ithe) * dang
        rthe = the * DEG_TO_RAD
        costhe = np.cos(rthe)
        sinthe = np.sin(rthe)

        # Match Fortran: fnumang=360./dang; numphi=nint(fnumang*sin(rthe))
        fnumang = 360.0 / dang
        numphi = int(np.round(fnumang * np.sin(rthe)))

        if numphi != 0:
            dphi = 360.0 / float(numphi)
        else:
            dphi = 10000.0

        # Match Fortran: do 4 iphi=0,int(359.9/dphi)
        for iphi in range(int(359.9 / dphi) + 1):
            phi = float(iphi) * dphi
            rphi = phi * DEG_TO_RAD
            cosphi = np.cos(rphi)
            sinphi = np.sin(rphi)

            # bb3 is the fault normal direction (matches Fortran lines 79-81)
            bb3_3 = costhe
            bb3_1 = sinthe * cosphi
            bb3_2 = sinthe * sinphi

            # bb1 is the slip direction (matches Fortran lines 82-84)
            bb1_3 = -sinthe
            bb1_1 = costhe * cosphi
            bb1_2 = costhe * sinphi

            # bb2 is orthogonal (cross product bb3 x bb1, matches Fortran line 85)
            bb2_1 = bb3_2 * bb1_3 - bb3_3 * bb1_2
            bb2_2 = bb3_3 * bb1_1 - bb3_1 * bb1_3
            bb2_3 = bb3_1 * bb1_2 - bb3_2 * bb1_1

            for izeta in range(max_zeta + 1):
                zeta = float(izeta) * dang
                rzeta = zeta * DEG_TO_RAD
                coszeta = np.cos(rzeta)
                sinzeta = np.sin(rzeta)

                if irot >= max_rot:
                    break

                # Store b3 (fault normal) - matches Fortran lines 96-98
                b3[0, irot] = bb3_1
                b3[1, irot] = bb3_2
                b3[2, irot] = bb3_3

                # Store b1 (slip, rotated by zeta) - matches Fortran lines 99-101
                b1[0, irot] = bb1_1 * coszeta + bb2_1 * sinzeta
                b1[1, irot] = bb1_2 * coszeta + bb2_2 * sinzeta
                b1[2, irot] = bb1_3 * coszeta + bb2_3 * sinzeta

                # Store b2 (orthogonal, rotated by zeta) - matches Fortran lines 102-104
                b2[0, irot] = bb2_1 * coszeta - bb1_1 * sinzeta
                b2[1, irot] = bb2_2 * coszeta - bb1_2 * sinzeta
                b2[2, irot] = bb2_3 * coszeta - bb1_3 * sinzeta

                irot += 1

    return b1[:, :irot], b2[:, :irot], b3[:, :irot], irot


def get_rotation_grid(dang):
    """
    Get or create the rotation grid for a given angle spacing.

    Results are cached to avoid recomputation.

    Parameters
    ----------
    dang : float
        Grid angle spacing in degrees

    Returns
    -------
    tuple
        (b1, b2, b3, nrot)
    """
    global _GRID_CACHE, _GRID_DANG_OLD

    cache_key = round(dang, 6)

    if cache_key in _GRID_CACHE:
        return _GRID_CACHE[cache_key]

    # Generate new grid
    b1, b2, b3, nrot = _setup_rotation_grid(dang)
    _GRID_CACHE[cache_key] = (b1, b2, b3, nrot)
    _GRID_DANG_OLD = dang

    return b1, b2, b3, nrot


@njit(cache=True)
def _count_misfits_single_trial(
    p_a1, p_a2, p_a3,
    p_pol, p_qual,
    b1, b2, b3, nrot,
    npol
):
    """
    Count polarity misfits for all rotations for a single trial.

    Parameters
    ----------
    p_a1, p_a2, p_a3 : ndarray
        Cartesian coordinates of ray directions
    p_pol : ndarray
        Polarity observations (+1 or -1)
    p_qual : ndarray
        Quality flags (0=impulsive, 1=emergent)
    b1, b2, b3 : ndarray
        Rotation arrays
    nrot : int
        Number of rotations
    npol : int
        Number of polarity observations

    Returns
    -------
    tuple
        (fit0, fit, nmiss01min) where:
        - fit0[i] = impulsive misfits for rotation i
        - fit[i] = total misfits for rotation i
        - nmiss01min[j] = minimum total misfits with j impulsive misfits
    """
    fit0 = np.zeros(nrot, dtype=np.int32)
    fit = np.zeros(nrot, dtype=np.int32)
    nmiss01min = np.full(npol + 1, 999, dtype=np.int32)

    for irot in range(nrot):
        nmiss = 0
        nmiss0 = 0

        b1_0, b1_1, b1_2 = b1[0, irot], b1[1, irot], b1[2, irot]
        b3_0, b3_1, b3_2 = b3[0, irot], b3[1, irot], b3[2, irot]

        for ista in range(npol):
            # Project ray direction onto fault frame
            p_b1 = b1_0 * p_a1[ista] + b1_1 * p_a2[ista] + b1_2 * p_a3[ista]
            p_b3 = b3_0 * p_a1[ista] + b3_1 * p_a2[ista] + b3_2 * p_a3[ista]

            # Predicted polarity from product
            prod = p_b1 * p_b3
            ipol = -1
            if prod > 0.0:
                ipol = 1

            # Check if prediction matches observation
            if ipol != p_pol[ista]:
                nmiss += 1
                if p_qual[ista] == 0:
                    nmiss0 += 1

        fit0[irot] = nmiss0
        fit[irot] = nmiss

        # Track minimums
        if nmiss < nmiss01min[nmiss0]:
            nmiss01min[nmiss0] = nmiss

    return fit0, fit, nmiss01min


@njit(cache=True)
def _find_best_mechanisms_all_trials(
    p_azi_mc, p_the_mc,
    p_pol, p_qual,
    b1, b2, b3, nrot,
    npol, nmc,
    nextra, ntotal
):
    """
    Find acceptable focal mechanisms across all Monte Carlo trials.

    This exactly matches the Fortran FOCALMC algorithm.

    Parameters
    ----------
    p_azi_mc : ndarray, shape (npol, nmc)
        Azimuths for each observation and trial
    p_the_mc : ndarray, shape (npol, nmc)
        Takeoff angles for each observation and trial
    p_pol : ndarray, shape (npol,)
        Polarity observations
    p_qual : ndarray, shape (npol,)
        Quality flags
    b1, b2, b3 : ndarray
        Rotation arrays
    nrot : int
        Number of rotations
    npol : int
        Number of polarity observations
    nmc : int
        Number of Monte Carlo trials
    nextra : int
        Additional misfits allowed above minimum
    ntotal : int
        Total allowed misfits

    Returns
    -------
    ndarray, shape (nrot,)
        Boolean array indicating which rotations meet criteria
    """
    irotgood = np.zeros(nrot, dtype=np.int32)

    for im in range(nmc):
        # Convert to Cartesian
        p_a1 = np.empty(npol, dtype=np.float64)
        p_a2 = np.empty(npol, dtype=np.float64)
        p_a3 = np.empty(npol, dtype=np.float64)

        for i in range(npol):
            theta = p_the_mc[i, im]
            phi = p_azi_mc[i, im]
            theta_rad = theta * DEG_TO_RAD
            phi_rad = phi * DEG_TO_RAD
            p_a3[i] = -1.0 * math.cos(theta_rad)
            p_a1[i] = math.sin(theta_rad) * math.cos(phi_rad)
            p_a2[i] = math.sin(theta_rad) * math.sin(phi_rad)

        # Count misfits
        fit0, fit, nmiss01min = _count_misfits_single_trial(
            p_a1, p_a2, p_a3, p_pol, p_qual, b1, b2, b3, nrot, npol
        )

        # Find minimums - exactly matching Fortran logic
        nmiss0min = 999
        nmissmin = 999
        for irot in range(nrot):
            if fit0[irot] < nmiss0min:
                nmiss0min = fit0[irot]
            if fit[irot] < nmissmin:
                nmissmin = fit[irot]

        # Choose fit criteria - exactly matching Fortran
        if nmiss0min == 0:
            nmiss0max = ntotal
            nmissmax = ntotal
        else:
            nmiss0max = ntotal
            nmissmax = npol

        if nmiss0max < nmiss0min + nextra:
            nmiss0max = nmiss0min + nextra
        # CRITICAL: Use nmiss01min(nmiss0min) not nmissmin - matches Fortran line 169-171
        if nmissmax < nmiss01min[nmiss0min] + nextra:
            nmissmax = nmiss01min[nmiss0min] + nextra

        # Mark good rotations
        for irot in range(nrot):
            if fit0[irot] <= nmiss0max and fit[irot] <= nmissmax:
                irotgood[irot] = 1

    return irotgood


def focalmc(p_azi_mc, p_the_mc, p_pol, p_qual, npol, nmc,
            dang, maxout, nextra, ntotal):
    """
    Perform grid search to find acceptable focal mechanisms.

    This implements the core FOCALMC algorithm from HASH v1.2.
    Uses a grid search over the focal sphere with Monte Carlo
    sampling for uncertainty.

    Parameters
    ----------
    p_azi_mc : ndarray, shape (npol, nmc)
        Azimuth to station from event (degrees, East of North)
    p_the_mc : ndarray, shape (npol, nmc)
        Takeoff angle (degrees, from vertical, up=0, <90 upgoing, >90 downgoing)
    p_pol : ndarray, shape (npol,)
        First motion polarity: 1=up, -1=down
    p_qual : ndarray, shape (npol,)
        Quality: 0=impulsive, 1=emergent
    npol : int
        Number of first motions
    nmc : int
        Number of trials
    dang : float
        Desired angle spacing for grid search (degrees)
    maxout : int
        Maximum number of fault planes to return
    nextra : int
        Additional misfits allowed above minimum
    ntotal : int
        Total number of allowed misfits

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

    # Find good rotations across all trials
    irotgood = _find_best_mechanisms_all_trials(
        p_azi_mc, p_the_mc, p_pol, p_qual,
        b1, b2, b3, nrot, npol, nmc, nextra, ntotal
    )

    # Count good rotations
    nfault = np.sum(irotgood)

    # Get indices of good rotations
    good_indices = np.where(irotgood > 0)[0]

    # Select output solutions
    nf = min(nfault, maxout)

    if nfault <= maxout:
        # Return all good rotations
        selected_indices = good_indices
    else:
        # Random selection
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
        'nf': nf,
        'strike': strike,
        'dip': dip,
        'rake': rake,
        'faults': faults,
        'slips': slips,
    }


@njit(cache=True)
def get_misfit_numba(
    npol, p_azi, p_the, p_pol, p_qual,
    strike, dip, rake
):
    """
    Calculate the fraction of misfit polarities for a given mechanism.

    Parameters
    ----------
    npol : int
        Number of polarity observations
    p_azi : ndarray
        Azimuths (degrees)
    p_the : ndarray
        Takeoff angles (degrees)
    p_pol : ndarray
        Polarity observations
    p_qual : ndarray
        Quality flags
    strike, dip, rake : float
        Mechanism parameters (degrees)

    Returns
    -------
    tuple
        (mfrac, stdr) - weighted fraction misfit and station distribution ratio
    """
    rad = DEG_TO_RAD

    # Convert strike, dip, rake to radians
    s_rad = strike * rad
    d_rad = dip * rad
    r_rad = rake * rad

    # Compute moment tensor (simplified for polarity calculation)
    sin_d = math.sin(d_rad)
    cos_d = math.cos(d_rad)
    sin_r = math.sin(r_rad)
    cos_r = math.cos(r_rad)
    sin_s = math.sin(s_rad)
    cos_s = math.cos(s_rad)
    sin_2s = math.sin(2.0 * s_rad)
    sin_2d = math.sin(2.0 * d_rad)
    cos_2s = math.cos(2.0 * s_rad)
    cos_2d = math.cos(2.0 * d_rad)

    # Moment tensor components
    M11 = -sin_d * cos_r * sin_2s - sin_2d * sin_r * sin_s * sin_s
    M22 = sin_d * cos_r * sin_2s - sin_2d * sin_r * cos_s * cos_s
    M33 = sin_2d * sin_r
    M12 = sin_d * cos_r * cos_2s + 0.5 * sin_2d * sin_r * math.sin(2.0 * s_rad)
    M13 = -cos_d * cos_r * cos_s - cos_2d * sin_r * sin_s
    M23 = -cos_d * cos_r * sin_s + cos_2d * sin_r * cos_s

    # Get fault normal and slip vectors
    fn3 = -cos_d
    fn1 = -sin_d * sin_s
    fn2 = sin_d * cos_s

    sl1 = cos_r * cos_s + cos_d * sin_r * sin_s
    sl2 = cos_r * sin_s - cos_d * sin_r * cos_s
    sl3 = -sin_r * sin_d

    # Get auxiliary plane (b2 = cross product of fn and sl)
    b21 = fn2 * sl3 - fn3 * sl2
    b22 = fn3 * sl1 - fn1 * sl3
    b23 = fn1 * sl2 - fn2 * sl1

    mfrac = 0.0
    qcount = 0.0
    scount = 0.0

    for k in range(npol):
        # Convert to Cartesian
        theta = p_the[k] * rad
        phi = p_azi[k] * rad

        p_a1 = math.sin(theta) * math.cos(phi)
        p_a2 = math.sin(theta) * math.sin(phi)
        p_a3 = -math.cos(theta)

        # Project onto mechanism
        p_b1 = sl1 * p_a1 + sl2 * p_a2 + sl3 * p_a3
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

        # Calculate predicted polarity from moment tensor
        azi = p_azi[k] * rad
        toff = p_the[k] * rad

        a1 = math.sin(toff) * math.cos(azi)
        a2 = math.sin(toff) * math.sin(azi)
        a3 = -math.cos(toff)

        b1 = M11 * a1 + M12 * a2 + M13 * a3
        b2 = M12 * a1 + M22 * a2 + M23 * a3
        b3 = M13 * a1 + M23 * a2 + M33 * a3

        dot_val = a1 * b1 + a2 * b2 + a3 * b3

        if dot_val < 0:
            pol = -1
        else:
            pol = 1

        if p_qual[k] == 0:
            wo = 1.0
        else:
            wo = 0.5

        if pol * p_pol[k] < 0:
            mfrac += wt * wo

        qcount += wt * wo
        scount += wo

    if qcount > 0:
        mfrac /= qcount
    if scount > 0:
        stdr = qcount / scount
    else:
        stdr = 0.0

    return mfrac, stdr


def get_misfit(npol, p_azi, p_the, p_pol, p_qual, strike, dip, rake):
    """
    Find the percent of misfit polarities for a given mechanism.

    Parameters
    ----------
    npol : int
        Number of polarity observations
    p_azi : ndarray, shape (npol,)
        Azimuths (degrees, East of North)
    p_the : ndarray, shape (npol,)
        Takeoff angles (degrees)
    p_pol : ndarray, shape (npol,)
        Polarity observations: 1=up, -1=down
    p_qual : ndarray, shape (npol,)
        Quality: 0=impulsive, 1=emergent
    strike, dip, rake : float
        Mechanism parameters (degrees)

    Returns
    -------
    tuple
        (mfrac, stdr) - weighted fraction misfit and station distribution ratio
    """
    return get_misfit_numba(
        npol, p_azi, p_the, p_pol, p_qual,
        strike, dip, rake
    )


@njit(cache=True)
def get_gap_numba(npol, p_azi, p_the):
    """
    Find maximum azimuthal and takeoff angle gaps.

    Parameters
    ----------
    npol : int
        Number of polarity observations
    p_azi : ndarray
        Azimuths (degrees)
    p_the : ndarray
        Takeoff angles (degrees)

    Returns
    -------
    tuple
        (magap, mpgap) - maximum azimuthal and takeoff angle gaps
    """
    # Convert takeoff angles > 90
    p2_azi = np.empty(npol, dtype=np.float64)
    p2_the = np.empty(npol, dtype=np.float64)

    for k in range(npol):
        if p_the[k] > 90.0:
            p2_the[k] = 180.0 - p_the[k]
            p2_azi[k] = p_azi[k] - 180.0
            if p2_azi[k] < 0.0:
                p2_azi[k] += 360.0
        else:
            p2_the[k] = p_the[k]
            p2_azi[k] = p_azi[k]

    # Sort arrays
    sorted_azi = np.sort(p2_azi)
    sorted_the = np.sort(p2_the)

    # Find gaps
    magap = 0.0
    mpgap = 0.0

    for k in range(1, npol):
        azi_gap = sorted_azi[k] - sorted_azi[k - 1]
        the_gap = sorted_the[k] - sorted_the[k - 1]

        if azi_gap > magap:
            magap = azi_gap
        if the_gap > mpgap:
            mpgap = the_gap

    # Check wraparound
    azi_wrap = sorted_azi[0] - sorted_azi[npol - 1] + 360.0
    if azi_wrap > magap:
        magap = azi_wrap

    # Check edges for takeoff
    the_edge1 = 90.0 - sorted_the[npol - 1]
    if the_edge1 > mpgap:
        mpgap = the_edge1

    the_edge0 = sorted_the[0]
    if the_edge0 > mpgap:
        mpgap = the_edge0

    return magap, mpgap


def get_gap(npol, p_azi, p_the):
    """
    Find maximum azimuthal and takeoff angle gaps.

    Parameters
    ----------
    npol : int
        Number of polarity observations
    p_azi : ndarray, shape (npol,)
        Azimuths (degrees)
    p_the : ndarray, shape (npol,)
        Takeoff angles (degrees)

    Returns
    -------
    tuple
        (magap, mpgap) - maximum azimuthal and takeoff angle gaps (degrees)
    """
    return get_gap_numba(npol, p_azi, p_the)


# Export all functions
__all__ = [
    "focalmc",
    "get_misfit",
    "get_gap",
    "get_rotation_grid",
]
