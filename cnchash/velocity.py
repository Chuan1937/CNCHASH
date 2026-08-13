"""
Velocity model and ray tracing for takeoff angle calculations.

Implements 1D ray tracing through layered velocity models
and table-based interpolation for takeoff angles.
"""

import math

import numpy as np
from numba import njit

from .utils import RAD_TO_DEG

# Velocity table cache
_VELOCITY_TABLES = {}
_TABLE_PARAMS = {}


def _default_table_params():
    """Get default table generation parameters."""
    return {
        "del1": 0.0,  # Minimum distance (km)
        "del2": 120.0,  # Maximum distance (km)
        "del3": 1.0,  # Distance step (km)
        "dep1": 0.0,  # Minimum depth (km)
        "dep2": 35.0,  # Maximum depth (km)
        "dep3": 1.0,  # Depth step (km)
        "pmin": 0.0,  # Minimum ray parameter
        "nump": 1000,  # Number of ray parameters
    }


@njit(cache=True)
def layer_trace(p, h, utop, ubot, imth):
    """
    Calculate travel time and range offset for ray tracing through a single layer.

    Parameters
    ----------
    p : float
        Horizontal slowness (s/km)
    h : float
        Layer thickness (km)
    utop, ubot : float
        Slowness at top and bottom of layer (s/km)
    imth : int
        Interpolation method:
        1: v(z) = 1/sqrt(a - 2*b*z)
        2: v(z) = a - b*z
        3: v(z) = a*exp(-b*z)

    Returns
    -------
    tuple
        (dx, dt, irtr) - range offset (km), travel time (s), return code
        irtr: -1=zero thickness, 0=ray turned above, 1=ray passed through,
              2=ray turned within layer
    """
    # Check for zero thickness layer
    if h == 0.0:
        return 0.0, 0.0, -1

    # Check for complex vertical slowness
    u = utop
    y = u - p
    if y <= 0.0:
        return 0.0, 0.0, 0

    q = y * (u + p)
    qs = math.sqrt(q)

    # Special function at top of layer
    if imth == 2:
        y_val = u + qs
        if p != 0.0:
            y_val = y_val / p
        qr = math.log(y_val)
    elif imth == 3:
        qr = math.atan2(qs, p)
    else:
        qr = 0.0  # Not used for imth=1

    # Determine b parameter
    if imth == 1:
        b = -(utop**2 - ubot**2) / (2.0 * h)
    elif imth == 2:
        vtop = 1.0 / utop if utop > 0 else 0.0
        vbot = 1.0 / ubot if ubot > 0 else 0.0
        b = -(vtop - vbot) / h
    else:  # imth == 3
        if utop > 0 and ubot > 0:
            b = -math.log(ubot / utop) / h
        else:
            b = 1.0 / h  # Default

    # Constant velocity layer
    if b == 0.0:
        b = 1.0 / h
        etau = qs
        ex_val = p / qs if qs > 0 else 0.0
        irtr = 1
    else:
        # Integral at upper limit (1/b factor omitted until end)
        if imth == 1:
            etau = -q * qs / 3.0
            ex_val = -qs * p
        elif imth == 2:
            ex_val = qs / u
            etau = qr - ex_val
            if p != 0.0:
                ex_val = ex_val / p
        else:  # imth == 3
            etau = qs - p * qr
            ex_val = qr

    # Check lower limit for turning point
    u_bot = ubot
    if u_bot <= p:
        irtr = 2
        dx = ex_val / b
        dtau = etau / b
        dt = dtau + p * dx
        return dx, dt, irtr

    irtr = 1
    q_bot = (u_bot - p) * (u_bot + p)
    qs_bot = math.sqrt(q_bot)

    if imth == 1:
        etau = etau + q_bot * qs_bot / 3.0
        ex_val = ex_val + qs_bot * p
    elif imth == 2:
        y_val = u_bot + qs_bot
        z_val = qs_bot / u_bot
        etau = etau + z_val
        if p != 0.0:
            y_val = y_val / p
            z_val = z_val / p
        qr_bot = math.log(y_val)
        etau = etau - qr_bot
        ex_val = ex_val - z_val
    else:  # imth == 3
        qr_bot = math.atan2(qs_bot, p)
        etau = etau - qs_bot + p * qr_bot
        ex_val = ex_val - qr_bot

    dx = ex_val / b
    dtau = etau / b
    dt = dtau + p * dx

    return dx, dt, irtr


def make_table_from_model(depth, velocity, params=None):
    """
    Create a takeoff angle table from a velocity model.

    Parameters
    ----------
    depth : ndarray
        Depth values (km)
    velocity : ndarray
        P-wave velocity values (km/s)
    params : dict, optional
        Table generation parameters

    Returns
    -------
    dict
        Table data with keys: 'table', 'delttab', 'deptab', 'ndel', 'ndep'
    """
    if params is None:
        params = _default_table_params()

    del1 = params.get("del1", 0.0)
    del2 = params.get("del2", 120.0)
    del3 = params.get("del3", 1.0)
    dep1 = params.get("dep1", 0.0)
    dep2 = params.get("dep2", 35.0)
    dep3 = params.get("dep3", 1.0)
    pmin = params.get("pmin", 0.0)
    nump = params.get("nump", 1000)

    # Set up depth grid
    ndep = int((dep2 - dep1) / dep3) + 1
    deptab = dep1 + dep3 * np.arange(ndep)

    # Set up distance grid
    ndel = int((del2 - del1) / del3) + 1
    delttab = del1 + del3 * np.arange(ndel)

    # Convert to slowness
    npts = len(depth)
    z = np.array(depth, dtype=np.float64)
    alpha = np.array(velocity, dtype=np.float64)
    slow = 1.0 / alpha

    # Extend model to include all table depths as layer boundaries,
    # following the original MK_TABLE insertion order: for each model
    # point from bottom to top, insert every table depth that falls
    # inside that layer (deepest table depth first).
    npts = len(depth)
    z = list(depth)
    alpha = list(velocity)
    npts_old = npts
    for i in range(npts_old - 1, 0, -1):
        for idep in range(ndep - 1, -1, -1):
            if z[i - 1] <= deptab[idep] - 0.1 and z[i] >= deptab[idep] + 0.1:
                dep = deptab[idep]
                z_top = z[i]
                z_bot = z[i - 1]
                a_top = alpha[i]
                a_bot = alpha[i - 1]
                frac = (dep - z_bot) / (z_top - z_bot)
                z.insert(i, dep)
                alpha.insert(i, a_bot + frac * (a_top - a_bot))
                npts += 1
    z = np.array(z, dtype=np.float64)
    alpha = np.array(alpha, dtype=np.float64)
    slow = 1.0 / alpha
    pmax = slow[0]
    plongcut = slow[-1]
    pstep = (pmax - pmin) / float(nump)

    # Ray tracing (2D correction arrays, mirroring the original MK_TABLE:
    # depxcor[np, idep] keeps the offset/travel time at each table depth
    # for every ray parameter).
    npmax = int((pmax + pstep / 2.0 - pmin) / pstep) + 1
    deltab = np.full(npmax, -999.0, dtype=np.float64)
    tttab = np.full(npmax, -999.0, dtype=np.float64)
    depxcor = np.full((npmax, ndep), -999.0, dtype=np.float64)
    deptcor = np.full((npmax, ndep), -999.0, dtype=np.float64)
    depucor = np.full((npmax, ndep), -999.0, dtype=np.float64)

    for np_idx in range(npmax):
        p = pmin + pstep * float(np_idx)

        x = 0.0
        t = 0.0
        imth = 3  # Exponential interpolation

        for idep in range(ndep):
            if deptab[idep] == 0.0:
                depxcor[np_idx, idep] = 0.0
                deptcor[np_idx, idep] = 0.0
                depucor[np_idx, idep] = slow[0]

        # Trace through layers
        for i in range(len(z) - 1):
            if z[i] >= 9999.0:
                deltab[np_idx] = -999.0
                tttab[np_idx] = -999.0
                break

            h = z[i + 1] - z[i]
            if h == 0.0:
                continue

            dx, dt, irtr = layer_trace(p, h, slow[i], slow[i + 1], imth)
            x += dx
            t += dt

            if irtr == 0 or irtr == 2:
                break

            # Check if we're at a table depth
            for idep in range(ndep):
                if abs(z[i + 1] - deptab[idep]) < 0.1:
                    depxcor[np_idx, idep] = x
                    deptcor[np_idx, idep] = t
                    depucor[np_idx, idep] = slow[i + 1]

        deltab[np_idx] = 2.0 * x
        tttab[np_idx] = 2.0 * t

    # Create takeoff angle table
    table = np.full((ndel, ndep), 0.0, dtype=np.float64)

    for idep in range(ndep):
        # Collect branches (original MK_TABLE logic: upgoing rays from the
        # source to the surface, then downgoing rays from the surface).
        xsave = []
        tsave = []
        psave = []
        usave = []

        if deptab[idep] == 0.0:
            i2 = npmax - 1
        else:
            i2 = npmax
            xold = -999.0
            for i in range(npmax):
                x2 = depxcor[i, idep]
                if x2 == -999.0:
                    continue
                if x2 <= xold:
                    i2 = i
                    break
                t2 = deptcor[i, idep]
                xsave.append(x2)
                tsave.append(t2)
                psave.append(-(pmin + pstep * float(i)))
                usave.append(depucor[i, idep])
                xold = x2
            i2 = i2 - 1

        # Downgoing rays (reflected) - from the last upgoing index down.
        for i in range(i2, -1, -1):
            if depxcor[i, idep] == -999.0:
                continue
            if deltab[i] == -999.0:
                continue
            x2 = deltab[i] - depxcor[i, idep]
            t2 = tttab[i] - deptcor[i, idep]
            xsave.append(x2)
            tsave.append(t2)
            psave.append(pmin + pstep * float(i))
            usave.append(depucor[i, idep])

        ncount = len(xsave)

        # Interpolate to grid. The angle stays 0.0 where no valid ray
        # exists (the original MK_TABLE/GET_TTS convention: 0 = invalid),
        # while the travel-time sentinel of 999.0 decides the minimum.
        for idel in range(ndel):
            del_dist = delttab[idel]
            best_t = 999.0

            for i in range(1, ncount):
                x1 = xsave[i - 1]
                x2 = xsave[i]

                if x1 > del_dist or x2 < del_dist:
                    continue

                if psave[i] > 0.0 and psave[i] < plongcut:
                    continue

                frac = (del_dist - x1) / (x2 - x1) if x2 != x1 else 0.0
                t1 = tsave[i - 1] + frac * (tsave[i] - tsave[i - 1])

                if t1 < best_t:
                    best_t = t1
                    # Calculate takeoff angle
                    scr1 = psave[i] / usave[i]
                    angle = math.asin(max(-1.0, min(1.0, scr1))) * RAD_TO_DEG

                    if angle < 0.0:
                        angle = -angle
                    else:
                        angle = 180.0 - angle

                    table[idel, idep] = angle

    # Set zero distance
    if delttab[0] == 0.0:
        for idep in range(ndep):
            table[0, idep] = 0.0

    return {
        "table": table,
        "delttab": delttab,
        "deptab": deptab,
        "ndel": ndel,
        "ndep": ndep,
    }


def make_table(vmodel_files, params=None):
    """
    Create takeoff angle tables from velocity model files.

    Parameters
    ----------
    vmodel_files : str or list
        Path(s) to velocity model file(s)
    params : dict, optional
        Table generation parameters

    Returns
    -------
    int
        Number of tables created
    """
    global _VELOCITY_TABLES, _TABLE_PARAMS

    if isinstance(vmodel_files, str):
        vmodel_files = [vmodel_files]

    if params is None:
        params = _default_table_params()

    _TABLE_PARAMS.update(params)

    ntab = len(vmodel_files)

    for itab, vmodel_file in enumerate(vmodel_files):
        # Read velocity model
        depth = []
        velocity = []

        with open(vmodel_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        depth.append(float(parts[0]))
                        velocity.append(float(parts[1]))
                    except ValueError:
                        continue

        # Create table
        table_data = make_table_from_model(depth, velocity, params)
        _VELOCITY_TABLES[itab] = table_data

    return ntab


def get_tts(ip, del_dist, qdep):
    """
    Get takeoff angle from velocity model table.

    Parameters
    ----------
    ip : int
        Index of velocity model (0-based)
    del_dist : float
        Distance (km)
    qdep : float
        Earthquake depth (km)

    Returns
    -------
    tuple
        (tt, iflag) - takeoff angle (degrees) and flag
        iflag: -1=outside depth range, 0=interpolation, 1=extrapolation
    """
    global _VELOCITY_TABLES

    if ip not in _VELOCITY_TABLES:
        return 999.0, -1

    table = _VELOCITY_TABLES[ip]
    delttab = table["delttab"]
    deptab = table["deptab"]
    tt = table["table"]
    ndel = table["ndel"]
    ndep = table["ndep"]

    # Check depth range
    if qdep < deptab[0] or qdep > deptab[-1]:
        return 999.0, -1

    # Find depth indices
    id1, id2 = 0, ndep - 1
    for id in range(1, ndep):
        if deptab[id] >= qdep:
            id2 = id
            id1 = id - 1
            break

    # Find distance indices
    ix1, ix2 = 0, ndel - 1
    for ix in range(1, ndel):
        if delttab[ix] >= del_dist:
            ix2 = ix
            ix1 = ix - 1
            break

    # Check for valid table values
    if (
        tt[ix1, id1] == 0.0
        or tt[ix1, id2] == 0.0
        or tt[ix2, id1] == 0.0
        or tt[ix2, id2] == 0.0
        or delttab[ix2] < del_dist
    ):
        # Need extrapolation
        iflag = 1

        # Find nearest valid values
        xoffmin1 = 999.0
        xoffmin2 = 999.0
        ixbest1 = -1
        ixbest2 = -1

        for ix in range(1, ndel):
            if tt[ix - 1, id1] != 0.0 and tt[ix, id1] != 0.0:
                xoff = abs((delttab[ix - 1] + delttab[ix]) / 2.0 - del_dist)
                if xoff < xoffmin1:
                    xoffmin1 = xoff
                    ixbest1 = ix

            if tt[ix - 1, id2] != 0.0 and tt[ix, id2] != 0.0:
                xoff = abs((delttab[ix - 1] + delttab[ix]) / 2.0 - del_dist)
                if xoff < xoffmin2:
                    xoffmin2 = xoff
                    ixbest2 = ix

        if ixbest1 < 0 or ixbest2 < 0:
            return 999.0, -1

        # Extrapolate
        xfrac1 = (del_dist - delttab[ixbest1 - 1]) / (delttab[ixbest1] - delttab[ixbest1 - 1])
        tt1 = tt[ixbest1 - 1, id1] + xfrac1 * (tt[ixbest1, id1] - tt[ixbest1 - 1, id1])

        xfrac2 = (del_dist - delttab[ixbest2 - 1]) / (delttab[ixbest2] - delttab[ixbest2 - 1])
        tt2 = tt[ixbest2 - 1, id2] + xfrac2 * (tt[ixbest2, id2] - tt[ixbest2 - 1, id2])

        dfrac = (qdep - deptab[id1]) / (deptab[id2] - deptab[id1])
        tt_final = tt1 + dfrac * (tt2 - tt1)

        return tt_final, 1

    # Bilinear interpolation
    iflag = 0
    xfrac = (del_dist - delttab[ix1]) / (delttab[ix2] - delttab[ix1])
    t1 = tt[ix1, id1] + xfrac * (tt[ix2, id1] - tt[ix1, id1])
    t2 = tt[ix1, id2] + xfrac * (tt[ix2, id2] - tt[ix1, id2])

    dfrac = (qdep - deptab[id1]) / (deptab[id2] - deptab[id1])
    tt_final = t1 + dfrac * (t2 - t1)

    return tt_final, iflag


# Export all functions
__all__ = [
    "make_table",
    "make_table_from_model",
    "get_tts",
    "_default_table_params",
]
