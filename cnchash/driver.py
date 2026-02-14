"""
Main driver program for CNCHASH (HASH v1.2 in Python).

This module provides high-level functions for running the HASH
algorithm on earthquake polarity data.
"""

import math

import numpy as np

from . import amp_subs, core, io, uncertainty, velocity
from .utils import (
    DEG_TO_RAD,
    RAD_TO_DEG,
    reset_random_seed,
)

# Global state for random number generator
_RANDOM_SEED = 314159


def run_hash(
    p_azi,
    p_the,
    p_pol,
    p_qual,
    dang=5.0,
    nmc=30,
    maxout=500,
    badfrac=0.1,
    cangle=45.0,
    prob_max=0.1,
    npolmin=8,
    max_agap=90.0,
    max_pgap=60.0,
):
    """
    Run HASH algorithm on polarity data.

    Parameters
    ----------
    p_azi : ndarray, shape (npol,) or (npol, nmc)
        Azimuths (degrees, East of North)
    p_the : ndarray, shape (npol,) or (npol, nmc)
        Takeoff angles (degrees)
    p_pol : ndarray, shape (npol,)
        Polarity observations: 1=up, -1=down
    p_qual : ndarray, shape (npol,)
        Quality: 0=impulsive, 1=emergent
    dang : float
        Grid angle for focal mechanism search (degrees)
    nmc : int
        Number of Monte Carlo trials
    maxout : int
        Maximum number of acceptable mechanisms
    badfrac : float
        Assumed fraction of bad polarities
    cangle : float
        Angle cutoff for mechanism probability (degrees)
    prob_max : float
        Probability threshold for multiple solutions
    npolmin : int
        Minimum number of polarities
    max_agap : float
        Maximum azimuthal gap (degrees)
    max_pgap : float
        Maximum takeoff angle gap (degrees)

    Returns
    -------
    dict
        Result containing:
        - 'success': bool, whether inversion succeeded
        - 'strike_avg', 'dip_avg', 'rake_avg': preferred mechanism
        - 'quality': quality rating (A-D, E, F)
        - 'prob': probability
        - 'var_est': estimated variance
        - 'mfrac': misfit fraction
        - 'stdr': station distribution ratio
        - 'nmult': number of solutions
        - 'nout2': number of acceptable mechanisms
        - And optionally multiple solutions if nmult > 1
    """
    global _RANDOM_SEED
    reset_random_seed(_RANDOM_SEED)

    npol = len(p_pol)

    # Convert to arrays if needed
    p_azi = np.asarray(p_azi, dtype=np.float64)
    p_the = np.asarray(p_the, dtype=np.float64)
    p_pol = np.asarray(p_pol, dtype=np.int32)
    p_qual = np.asarray(p_qual, dtype=np.int32)

    # Ensure shape is (npol, nmc)
    if p_azi.ndim == 1:
        p_azi = p_azi.reshape(-1, 1)
        p_the = p_the.reshape(-1, 1)
        if nmc > 1:
            # Need to add uncertainty - for now just duplicate
            p_azi = np.repeat(p_azi, nmc, axis=1)
            p_the = np.repeat(p_the, nmc, axis=1)

    # Check minimum polarities
    if npol < npolmin:
        return {
            "success": False,
            "strike_avg": 999.0,
            "dip_avg": 99.0,
            "rake_avg": 999.0,
            "var_est": [99.0, 99.0],
            "mfrac": 0.99,
            "quality": "F",
            "prob": 0.0,
            "stdr": 0.0,
            "nmult": 0,
            "nout2": 0,
            "nout1": 0,
        }

    # Check gaps
    magap, mpgap = core.get_gap(npol, p_azi[:, 0], p_the[:, 0])

    if magap > max_agap or mpgap > max_pgap:
        return {
            "success": False,
            "strike_avg": 999.0,
            "dip_avg": 99.0,
            "rake_avg": 999.0,
            "var_est": [99.0, 99.0],
            "mfrac": 0.99,
            "quality": "E",
            "prob": 0.0,
            "stdr": 0.0,
            "nmult": 0,
            "nout2": 0,
            "nout1": 0,
        }

    # Calculate acceptance criteria
    nmismax = max(int(npol * badfrac), 2)
    nextra = max(int(npol * badfrac * 0.5), 2)

    # Find acceptable mechanisms
    result = core.focalmc(p_azi, p_the, p_pol, p_qual, npol, nmc, dang, maxout, nextra, nmismax)

    nf = result["nf"]

    if nf == 0:
        return {
            "success": False,
            "strike_avg": 999.0,
            "dip_avg": 99.0,
            "rake_avg": 999.0,
            "var_est": [99.0, 99.0],
            "mfrac": 0.99,
            "quality": "F",
            "prob": 0.0,
            "stdr": 0.0,
            "nmult": 0,
            "nout2": 0,
            "nout1": 0,
        }

    # Find preferred mechanism
    prob_result = uncertainty.mech_prob(nf, result["faults"], result["slips"], cangle, prob_max)

    nsltn = prob_result["nsltn"]
    strike_avg = prob_result["strike_avg"]
    dip_avg = prob_result["dip_avg"]
    rake_avg = prob_result["rake_avg"]
    prob = prob_result["prob"]
    rms_diff = prob_result["rms_diff"]

    # Calculate misfit for preferred solution
    if nsltn > 0:
        mfrac, stdr = core.get_misfit(
            npol, p_azi[:, 0], p_the[:, 0], p_pol, p_qual, strike_avg[0], dip_avg[0], rake_avg[0]
        )
    else:
        mfrac = 0.99
        stdr = 0.0

    # Determine quality rating
    quality = "D"
    var_avg = (rms_diff[0, 0] + rms_diff[1, 0]) / 2.0 if nsltn > 0 else 99.0

    if nsltn > 0:
        if prob[0] > 0.8 and var_avg <= 25.0 and mfrac <= 0.15 and stdr >= 0.5:
            quality = "A"
        elif prob[0] > 0.6 and var_avg <= 35.0 and mfrac <= 0.2 and stdr >= 0.4:
            quality = "B"
        elif prob[0] > 0.5 and var_avg <= 45.0 and mfrac <= 0.3 and stdr >= 0.3:
            quality = "C"
        else:
            quality = "D"

    # Build result
    output = {
        "success": True,
        "strike_avg": strike_avg[0] if nsltn > 0 else result["strike"][0],
        "dip_avg": dip_avg[0] if nsltn > 0 else result["dip"][0],
        "rake_avg": rake_avg[0] if nsltn > 0 else result["rake"][0],
        "var_est": [rms_diff[0, 0], rms_diff[1, 0]] if nsltn > 0 else [99.0, 99.0],
        "mfrac": mfrac,
        "quality": quality,
        "prob": prob[0] if nsltn > 0 else 0.0,
        "stdr": stdr,
        "nmult": nsltn,
        "nout2": nf,
        "nout1": min(nf, maxout),
    }

    # Add multiple solutions if present
    if nsltn > 1:
        output["strike_avg"] = strike_avg[:nsltn]
        output["dip_avg"] = dip_avg[:nsltn]
        output["rake_avg"] = rake_avg[:nsltn]
        output["prob"] = prob[:nsltn]
        output["rms_diff"] = rms_diff[:, :nsltn]
        output["quality"] = [quality] * nsltn  # Simplified - could compute per solution

    # Also store raw results
    output["faults"] = result["faults"]
    output["slips"] = result["slips"]
    output["strike"] = result["strike"]
    output["dip"] = result["dip"]
    output["rake"] = result["rake"]

    return output


def run_hash_with_amp(
    p_azi,
    p_the,
    p_pol,
    sp_amp,
    dang=5.0,
    nmc=30,
    maxout=500,
    nmismax=None,
    qbadfac=0.3,
    cangle=45.0,
    prob_max=0.1,
    npolmin=8,
    max_agap=90.0,
    max_pgap=60.0,
):
    """
    Run HASH algorithm with both P-wave polarities and S/P amplitude ratios.

    This implements the FOCALAMP_MC algorithm from Hardebeck and Shearer (2003).

    Parameters
    ----------
    p_azi : ndarray, shape (npsta,) or (npsta, nmc)
        Azimuths (degrees, East of North)
    p_the : ndarray, shape (npsta,) or (npsta, nmc)
        Takeoff angles (degrees)
    p_pol : ndarray, shape (npsta,)
        Polarity observations: 1=up, -1=down, 0=no polarity
    sp_amp : ndarray, shape (npsta,)
        S/P amplitude ratios (log10 scale, 0 if no ratio)
    dang : float
        Grid angle for focal mechanism search (degrees)
    nmc : int
        Number of Monte Carlo trials
    maxout : int
        Maximum number of acceptable mechanisms
    nmismax : int or None
        Maximum number of polarity misfits (default: max(npsta*0.1, 2))
    qbadfac : float
        Assumed noise in amplitude ratios, log10 (0.3 = factor of 2)
    cangle : float
        Angle cutoff for mechanism probability (degrees)
    prob_max : float
        Probability threshold for multiple solutions
    npolmin : int
        Minimum number of polarities
    max_agap : float
        Maximum azimuthal gap (degrees)
    max_pgap : float
        Maximum takeoff angle gap (degrees)

    Returns
    -------
    dict
        Result containing:
        - 'success': bool, whether inversion succeeded
        - 'strike_avg', 'dip_avg', 'rake_avg': preferred mechanism
        - 'quality': quality rating (A-D, E, F)
        - 'prob': probability
        - 'var_est': estimated variance
        - 'mfrac': polarity misfit fraction
        - 'mavg': average S/P amplitude misfit (log10)
        - 'stdr': station distribution ratio
        - 'nmult': number of solutions
        - 'nout2': number of acceptable mechanisms
        - And optionally multiple solutions if nmult > 1
    """
    global _RANDOM_SEED
    reset_random_seed(_RANDOM_SEED)

    npsta = len(p_pol)

    # Count valid polarities and S/P ratios
    npol = np.sum(p_pol != 0)
    nspr = np.sum(sp_amp != 0.0)

    # Convert to arrays if needed
    p_azi = np.asarray(p_azi, dtype=np.float64)
    p_the = np.asarray(p_the, dtype=np.float64)
    p_pol = np.asarray(p_pol, dtype=np.int32)
    sp_amp = np.asarray(sp_amp, dtype=np.float64)

    # Ensure shape is (npsta, nmc)
    if p_azi.ndim == 1:
        p_azi = p_azi.reshape(-1, 1)
        p_the = p_the.reshape(-1, 1)
        if nmc > 1:
            p_azi = np.repeat(p_azi, nmc, axis=1)
            p_the = np.repeat(p_the, nmc, axis=1)

    # Check minimum polarities
    if npol < npolmin:
        return {
            "success": False,
            "strike_avg": 999.0,
            "dip_avg": 99.0,
            "rake_avg": 999.0,
            "var_est": [99.0, 99.0],
            "mfrac": 0.99,
            "mavg": 0.99,
            "quality": "F",
            "prob": 0.0,
            "stdr": 0.0,
            "nmult": 0,
            "nout2": 0,
            "nout1": 0,
        }

    # Check gaps
    magap, mpgap = core.get_gap(npsta, p_azi[:, 0], p_the[:, 0])

    if magap > max_agap or mpgap > max_pgap:
        return {
            "success": False,
            "strike_avg": 999.0,
            "dip_avg": 99.0,
            "rake_avg": 999.0,
            "var_est": [99.0, 99.0],
            "mfrac": 0.99,
            "mavg": 0.99,
            "quality": "E",
            "prob": 0.0,
            "stdr": 0.0,
            "nmult": 0,
            "nout2": 0,
            "nout1": 0,
        }

    # Calculate acceptance criteria
    if nmismax is None:
        nmismax = max(int(npol * 0.1), 2)
    nextra = max(int(npol * 0.05), 0)

    # S/P amplitude ratio criteria
    qmismax = nspr * qbadfac
    qextra = max(nspr * qbadfac * 0.5, 2.0)

    # Find acceptable mechanisms using both polarities and S/P ratios
    result = amp_subs.focalamp_mc(
        p_azi, p_the, sp_amp, p_pol, npsta, nmc, dang, maxout, nextra, nmismax, qextra, qmismax
    )

    nf = result["nf"]

    if nf == 0:
        return {
            "success": False,
            "strike_avg": 999.0,
            "dip_avg": 99.0,
            "rake_avg": 999.0,
            "var_est": [99.0, 99.0],
            "mfrac": 0.99,
            "mavg": 0.99,
            "quality": "F",
            "prob": 0.0,
            "stdr": 0.0,
            "nmult": 0,
            "nout2": 0,
            "nout1": 0,
        }

    # Find preferred mechanism
    prob_result = uncertainty.mech_prob(nf, result["faults"], result["slips"], cangle, prob_max)

    nsltn = prob_result["nsltn"]
    strike_avg = prob_result["strike_avg"]
    dip_avg = prob_result["dip_avg"]
    rake_avg = prob_result["rake_avg"]
    prob = prob_result["prob"]
    rms_diff = prob_result["rms_diff"]

    # Calculate misfit for preferred solution
    if nsltn > 0:
        mfrac, mavg, stdr = amp_subs.get_misf_amp(
            npsta, p_azi[:, 0], p_the[:, 0], sp_amp, p_pol, strike_avg[0], dip_avg[0], rake_avg[0]
        )
    else:
        mfrac = 0.99
        mavg = 0.99
        stdr = 0.0

    # Determine quality rating (using both polarity and amplitude criteria)
    quality = "D"
    var_avg = (rms_diff[0, 0] + rms_diff[1, 0]) / 2.0 if nsltn > 0 else 99.0

    if nsltn > 0:
        if prob[0] > 0.8 and var_avg <= 25.0 and mfrac <= 0.15 and mavg <= 0.2 and stdr >= 0.5:
            quality = "A"
        elif prob[0] > 0.6 and var_avg <= 35.0 and mfrac <= 0.2 and mavg <= 0.3 and stdr >= 0.4:
            quality = "B"
        elif prob[0] > 0.5 and var_avg <= 45.0 and mfrac <= 0.3 and mavg <= 0.4 and stdr >= 0.3:
            quality = "C"
        else:
            quality = "D"

    # Build result
    output = {
        "success": True,
        "strike_avg": strike_avg[0] if nsltn > 0 else result["strike"][0],
        "dip_avg": dip_avg[0] if nsltn > 0 else result["dip"][0],
        "rake_avg": rake_avg[0] if nsltn > 0 else result["rake"][0],
        "var_est": [rms_diff[0, 0], rms_diff[1, 0]] if nsltn > 0 else [99.0, 99.0],
        "mfrac": mfrac,
        "mavg": mavg,
        "quality": quality,
        "prob": prob[0] if nsltn > 0 else 0.0,
        "stdr": stdr,
        "nmult": nsltn,
        "nout2": nf,
        "nout1": min(nf, maxout),
        "npol": npol,
        "nspr": nspr,
    }

    # Add multiple solutions if present
    if nsltn > 1:
        output["strike_avg"] = strike_avg[:nsltn]
        output["dip_avg"] = dip_avg[:nsltn]
        output["rake_avg"] = rake_avg[:nsltn]
        output["prob"] = prob[:nsltn]
        output["rms_diff"] = rms_diff[:, :nsltn]
        output["quality"] = [quality] * nsltn

    # Also store raw results
    output["faults"] = result["faults"]
    output["slips"] = result["slips"]
    output["strike"] = result["strike"]
    output["dip"] = result["dip"]
    output["rake"] = result["rake"]

    return output


def run_hash_from_file(input_file):
    """
    Run HASH from an input file.

    Parameters
    ----------
    input_file : str
        Path to HASH input file

    Returns
    -------
    list
        List of results, one per event
    """
    # Read input parameters
    params = io.read_hash_input_file(input_file)

    # Read phase data
    events = io.read_phase_file(params["phasefile"])

    # Read station data
    stations = io.read_station_file(params.get("station_file", ""))
    if not stations:
        # Use fallback - try to find station file
        import os

        base_dir = os.path.dirname(input_file)
        for st_file in ["scsn.stations", "stations.txt"]:
            st_path = os.path.join(base_dir, st_file)
            if os.path.exists(st_path):
                stations = io.read_station_file(st_path)
                break

    # Read polarity reversal file
    reversals = io.read_polarity_reversal_file(params["polfile"])

    # Read velocity model if needed
    vmodel_file = params.get("vmodel_file", "")
    if vmodel_file and os.path.exists(vmodel_file):
        depth, vel = io.read_velocity_model(vmodel_file)
        velocity.make_table_from_model(depth, vel)

    results = []

    for event in events:
        # Process event
        result = process_event(event, stations, reversals, params)
        results.append(result)

    # Write output files
    if results:
        io.write_mechanism_output(params["outfile1"], events, results)

        for event, result in zip(events, results, strict=True):
            if result.get("success") and result.get("nout2", 0) > 0:
                io.write_acceptable_planes(params["outfile2"], event, result)
                break  # Only write first successful event for example

    return results


def process_event(event, stations, reversals, params):
    """
    Process a single event.

    Parameters
    ----------
    event : dict
        Event data
    stations : dict
        Station locations (includes '_byname' index for fast lookup)
    reversals : dict
        Polarity reversal data
    params : dict
        HASH parameters

    Returns
    -------
    dict
        Result from run_hash
    """
    # Extract event location
    ev_lat = event["lat"]
    ev_lon = event["lon"]
    ev_dep = event["depth"]
    ev_year = event["year"]
    ev_month = event["month"]
    ev_day = event["day"]

    # Convert event time to integer format for polarity check
    ev_date_int = ev_year * 10000 + ev_month * 100 + ev_day

    # Get fast lookup index
    stations_by_name = stations.get("_byname", {})

    # Process station data
    p_azi = []
    p_the = []
    p_pol = []
    p_qual = []

    for sta in event["stations"]:
        sta_name = sta["name"]

        # Fast lookup: try by-name index first (O(1))
        if sta_name in stations_by_name:
            sta_lat, sta_lon, _ = stations_by_name[sta_name]
        else:
            # Fallback to full key lookup
            sta_net = sta.get("network", "CI")
            sta_comp = sta.get("component", "HHZ")

            key = (sta_name, sta_comp, sta_net)
            if key not in stations:
                # Try with different component variations
                found = False
                for comp in ["HHZ", "VHZ", "EHZ", "BHZ", "HH", "VH", "EH", "BH"]:
                    key = (sta_name, comp, sta_net)
                    if key in stations:
                        found = True
                        break
                if not found:
                    continue

            sta_lat, sta_lon, _ = stations[key]

        # Calculate azimuth and distance
        dx = (sta_lon - ev_lon) * 111.2 * math.cos(ev_lat * DEG_TO_RAD)
        dy = (sta_lat - ev_lat) * 111.2
        distance = math.sqrt(dx**2 + dy**2)

        azi = 90.0 - math.atan2(dy, dx) * RAD_TO_DEG
        if azi < 0.0:
            azi += 360.0

        # Check distance limit
        if distance > params["delmax"]:
            continue

        # Get polarity
        pol_char = sta["polarity"].upper()
        if pol_char in ("U", "+"):
            pol = 1
        elif pol_char in ("D", "-"):
            pol = -1
        else:
            continue

        # Check polarity reversal
        if sta_name in reversals:
            for start_date, end_date in reversals[sta_name]:
                if start_date <= ev_date_int <= end_date:
                    pol = -pol
                    break

        # Get quality from onset
        onset = sta.get("onset", "I").upper()
        qual = 0 if onset == "I" else 1

        # For this version, use a simple takeoff angle approximation
        # (In full version, would use velocity model table)
        the = 180.0 - math.degrees(math.atan2(distance, ev_dep)) if ev_dep > 0 else 90.0

        p_azi.append(azi)
        p_the.append(the)
        p_pol.append(pol)
        p_qual.append(qual)

    if len(p_pol) == 0:
        return {
            "success": False,
            "quality": "F",
            "strike_avg": 999.0,
            "dip_avg": 99.0,
            "rake_avg": 999.0,
        }

    # Add Monte Carlo perturbation using vectorized numpy
    nmc = params["nmc"]
    npol = len(p_azi)

    # Vectorized: broadcast base values and add random perturbation
    p_azi_base = np.array(p_azi, dtype=np.float64).reshape(-1, 1)
    p_the_base = np.array(p_the, dtype=np.float64).reshape(-1, 1)

    # Generate all random perturbations at once (Box-Muller transform)
    np.random.seed(_RANDOM_SEED)
    u1 = np.random.random((npol, nmc - 1))
    u2 = np.random.random((npol, nmc - 1))
    z = np.sqrt(-2.0 * np.log(u1)) * np.cos(2.0 * np.pi * u2) * 5.0

    # Combine: first column is original, rest are perturbed
    p_azi_arr = np.hstack([p_azi_base, p_azi_base + z])
    p_the_arr = np.hstack(
        [p_the_base, p_the_base + np.sqrt(-2.0 * np.log(u1)) * np.sin(2.0 * np.pi * u2) * 5.0]
    )

    # Run HASH
    result = run_hash(
        p_azi_arr,
        p_the_arr,
        np.array(p_pol),
        np.array(p_qual),
        dang=params["dang"],
        nmc=nmc,
        maxout=params["maxout"],
        badfrac=params["badfrac"],
        cangle=params["cangle"],
        prob_max=params["prob_max"],
        npolmin=params["npolmin"],
        max_agap=params["max_agap"],
        max_pgap=params["max_pgap"],
    )

    # Add event info to result
    result["npol"] = len(p_pol)

    return result


# Export all functions
__all__ = [
    "run_hash",
    "run_hash_with_amp",
    "run_hash_from_file",
    "process_event",
]
