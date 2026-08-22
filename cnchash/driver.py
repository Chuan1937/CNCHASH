"""
Main driver program for CNCHASH (HASH v1.2 in Python).

This module provides high-level functions for running the HASH
algorithm on earthquake polarity data. Computation is dispatched to a
backend: the Modern Fortran/OpenMP native backend (default when the
compiled library is available) or the pure Python/Numba reference
backend.
"""

import math
import os
import warnings

import numpy as np

from . import backend, io
from .takeoff import TakeoffRangeError, TakeoffTable
from .utils import (
    DEG_TO_RAD,
    RAD_TO_DEG,
)

# Global state for random number generator
_RANDOM_SEED = 314159


def _resolve_backend(backend_name):
    """Resolve a backend name ('auto' -> best available)."""
    return backend.get_backend(backend_name)


def _validate_polarity_inputs(p_azi, p_the, p_pol, p_qual):
    """Validate array shapes before crossing into the native ABI.

    Wrong sizes are memory hazards in the assumed-size C interface, so
    they are rejected here with clear Python errors instead.
    """
    p_azi = np.asarray(p_azi, dtype=np.float64)
    p_the = np.asarray(p_the, dtype=np.float64)
    p_pol = np.asarray(p_pol, dtype=np.int32)
    p_qual = np.asarray(p_qual, dtype=np.int32)

    if p_pol.ndim != 1 or p_qual.ndim != 1:
        raise ValueError("p_pol and p_qual must be 1-D arrays")
    if p_pol.size != p_qual.size:
        raise ValueError(
            f"p_pol and p_qual must have the same length ({p_pol.size} vs {p_qual.size})"
        )
    if p_azi.ndim not in (1, 2) or p_the.ndim != p_azi.ndim:
        raise ValueError("p_azi and p_the must both be 1-D or both be 2-D")
    if p_azi.shape != p_the.shape:
        raise ValueError(
            f"p_azi and p_the must have the same shape ({p_azi.shape} vs {p_the.shape})"
        )
    if p_azi.shape[0] != p_pol.size:
        raise ValueError(
            f"p_azi first dimension ({p_azi.shape[0]}) must match p_pol length ({p_pol.size})"
        )
    return p_azi, p_the, p_pol, p_qual


def _validate_amp_inputs(p_azi, p_the, p_pol, sp_amp):
    """Validate amplitude inputs; sp_amp is 1-D and matches p_pol."""
    p_azi, p_the, p_pol, _ = _validate_polarity_inputs(p_azi, p_the, p_pol, np.zeros_like(p_pol))
    sp_amp = np.asarray(sp_amp, dtype=np.float64)
    if sp_amp.ndim != 1:
        raise ValueError("sp_amp must be a 1-D array")
    if sp_amp.size != p_pol.size:
        raise ValueError(f"sp_amp length ({sp_amp.size}) must match p_pol length ({p_pol.size})")
    return p_azi, p_the, p_pol, sp_amp


def run_hash_event_payload(payload):
    """
    Solve one event from a serialized payload.

    This helper is process-safe and pickle-friendly for use with
    ProcessPoolExecutor (spawn context) in notebooks and scripts.

    Parameters
    ----------
    payload : tuple
        (p_azi_mc, p_the_mc, p_pol, p_qual, dang, nmc)

    Returns
    -------
    dict
        Result dictionary from run_hash.
    """
    p_azi_mc, p_the_mc, p_pol, p_qual, dang, nmc = payload
    return run_hash(p_azi_mc, p_the_mc, p_pol, p_qual, dang=dang, nmc=nmc, num_threads=1)


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
    num_threads=None,
    backend="auto",
    selection=0,
):
    """Run HASH algorithm on polarity data.

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
    num_threads : int or None
        Thread count for the backend (Numba threads or OpenMP threads).
    backend : str
        "auto" (default) or "fortran".
    selection : int
        0 = deterministic solution selection (default), 1 = HASH-style
        random selection (reproduces the original Fortran behavior when
        more than maxout mechanisms are found).

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
    solver = _resolve_backend(backend)
    if num_threads is not None:
        solver.set_num_threads(num_threads)

    p_azi, p_the, p_pol, p_qual = _validate_polarity_inputs(p_azi, p_the, p_pol, p_qual)
    npol = p_pol.size

    if p_azi.ndim == 1:
        p_azi = p_azi.reshape(-1, 1)
        p_the = p_the.reshape(-1, 1)
        if nmc > 1:
            p_azi = np.repeat(p_azi, nmc, axis=1)
            p_the = np.repeat(p_the, nmc, axis=1)

    return solver.run_event(
        p_azi,
        p_the,
        p_pol,
        p_qual,
        npol,
        p_azi.shape[1],
        dang,
        maxout,
        badfrac,
        cangle,
        prob_max,
        npolmin,
        max_agap,
        max_pgap,
        selection,
    )


def run_hash_with_amp(
    p_azi,
    p_the,
    p_pol,
    sp_amp,
    dang=5.0,
    nmc=30,
    maxout=500,
    badfrac=0.1,
    nmismax=None,
    qbadfac=0.3,
    cangle=45.0,
    prob_max=0.1,
    npolmin=8,
    max_agap=90.0,
    max_pgap=60.0,
    num_threads=None,
    backend="auto",
    selection=0,
):
    """Run HASH algorithm with both P-wave polarities and S/P amplitude ratios.

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
    num_threads : int or None
        Thread count for the backend (Numba threads or OpenMP threads).
    backend : str
        "auto" (default) or "fortran".
    selection : int
        0 = deterministic solution selection (default), 1 = HASH-style
        random selection.

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
    solver = _resolve_backend(backend)
    if num_threads is not None:
        solver.set_num_threads(num_threads)

    p_azi, p_the, p_pol, sp_amp = _validate_amp_inputs(p_azi, p_the, p_pol, sp_amp)
    npsta = p_pol.size

    if p_azi.ndim == 1:
        p_azi = p_azi.reshape(-1, 1)
        p_the = p_the.reshape(-1, 1)
        if nmc > 1:
            p_azi = np.repeat(p_azi, nmc, axis=1)
            p_the = np.repeat(p_the, nmc, axis=1)

    result = solver.run_event_amp(
        p_azi,
        p_the,
        p_pol,
        sp_amp,
        npsta,
        p_azi.shape[1],
        dang,
        maxout,
        badfrac=badfrac,
        qbadfac=qbadfac,
        cangle=cangle,
        prob_max=prob_max,
        npolmin=npolmin,
        max_agap=max_agap,
        max_pgap=max_pgap,
        selection=selection,
        nmismax=nmismax,
    )
    return result


def run_hash_batch(
    events,
    dang=5.0,
    nmc=30,
    maxout=500,
    badfrac=0.1,
    cangle=45.0,
    prob_max=0.1,
    npolmin=8,
    max_agap=90.0,
    max_pgap=60.0,
    num_threads=None,
    backend="auto",
    selection=0,
):
    """Run HASH on many events in one call.

    Parameters
    ----------
    events : iterable of tuples
        Each tuple is (p_azi, p_the, p_pol, p_qual) with the same meaning
        as in :func:`run_hash`. Arrays may be 1D (duplicated across MC
        trials) or 2D (npol, nmc).
    dang : float
        Grid angle for focal mechanism search (degrees)
    nmc : int
        Number of Monte Carlo trials (used when inputs are 1D)
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
    num_threads : int or None
        Thread count for the backend.
    backend : str
        "auto" (default) or "fortran".
    selection : int
        Solution selection mode (0 = deterministic, 1 = HASH random).

    Returns
    -------
    list of dict
        One result per event (same schema as :func:`run_hash`).
    """
    solver = _resolve_backend(backend)
    if num_threads is not None:
        solver.set_num_threads(num_threads)

    prepared = []
    for p_azi, p_the, p_pol, p_qual in events:
        p_azi, p_the, p_pol, p_qual = _validate_polarity_inputs(p_azi, p_the, p_pol, p_qual)
        if p_azi.ndim == 1:
            p_azi = p_azi.reshape(-1, 1)
            p_the = p_the.reshape(-1, 1)
            if nmc > 1:
                p_azi = np.repeat(p_azi, nmc, axis=1)
                p_the = np.repeat(p_the, nmc, axis=1)
        prepared.append((p_azi, p_the, p_pol, p_qual))

    return solver.run_batch(
        prepared,
        dang,
        maxout,
        badfrac,
        cangle,
        prob_max,
        npolmin,
        max_agap,
        max_pgap,
        selection,
    )


def run_hash_batch_with_amp(
    events,
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
    num_threads=None,
    backend="auto",
    selection=0,
):
    """Run HASH with S/P amplitude ratios on many events in one call.

    Parameters
    ----------
    events : iterable of tuples
        Each tuple is (p_azi, p_the, p_pol, sp_amp) with the same meaning
        as in :func:`run_hash_with_amp`.

    Returns
    -------
    list of dict
        One result per event.
    """
    solver = _resolve_backend(backend)
    if num_threads is not None:
        solver.set_num_threads(num_threads)

    prepared = []
    for p_azi, p_the, p_pol, sp_amp in events:
        p_azi, p_the, p_pol, sp_amp = _validate_amp_inputs(p_azi, p_the, p_pol, sp_amp)
        if p_azi.ndim == 1:
            p_azi = p_azi.reshape(-1, 1)
            p_the = p_the.reshape(-1, 1)
            if nmc > 1:
                p_azi = np.repeat(p_azi, nmc, axis=1)
                p_the = np.repeat(p_the, nmc, axis=1)
        prepared.append((p_azi, p_the, p_pol, sp_amp))

    return solver.run_batch_amp(
        prepared,
        dang,
        maxout,
        badfrac=0.1,
        qbadfac=qbadfac,
        cangle=cangle,
        prob_max=prob_max,
        npolmin=npolmin,
        max_agap=max_agap,
        max_pgap=max_pgap,
        selection=selection,
        nmismax=nmismax,
    )


def run_hash_from_file(input_file, velocity_models=None):
    """
    Run HASH from an input file.

    Parameters
    ----------
    input_file : str
        Path to HASH input file
    velocity_models : list, optional
        Velocity models as a list of ``(depth, velocity)`` tuples or
        TakeoffTable instances. When None, velocity model file names
        listed at the end of the input file (lines starting with
        ``vz.``) are loaded automatically.

    Returns
    -------
    list
        List of results, one per event
    """
    # Read input parameters
    params = io.read_hash_input_file(input_file)

    # Read phase data
    events = io.read_phase_file(params["phasefile"])

    # Read station data (formats 1 and 5 have no station file)
    stations = {}
    if params.get("station_file"):
        stations = io.read_station_file(params["station_file"])
    if not stations:
        # Use fallback - try to find station file
        base_dir = os.path.dirname(input_file)
        for st_file in ["scsn.stations", "stations.txt"]:
            st_path = os.path.join(base_dir, st_file)
            if os.path.exists(st_path):
                stations = io.read_station_file(st_path)
                break

    # Read polarity reversal file
    reversals = io.read_polarity_reversal_file(params["polfile"])

    # Velocity models: prefer the explicit argument, otherwise load the
    # model files listed in the input file (relative to the input file).
    if velocity_models is None:
        model_files = params.get("velocity_models") or []
        if model_files:
            base_dir = os.path.dirname(os.path.abspath(input_file))
            velocity_models = []
            for model_file in model_files:
                path = (
                    model_file if os.path.isabs(model_file) else os.path.join(base_dir, model_file)
                )
                depth, velocity = io.read_velocity_model(path)
                velocity_models.append((depth, velocity))

    # Optional per-event auxiliary data (drivers 3 and 5), resolved
    # relative to the input file directory.
    base_dir = os.path.dirname(os.path.abspath(input_file))

    amp_data = None
    if params.get("amp_file"):
        amp_path = params["amp_file"]
        if not os.path.isabs(amp_path):
            amp_path = os.path.join(base_dir, amp_path)
        amp_data = io.read_amp_file(amp_path)

    statcor = None
    if params.get("statcor_file"):
        cor_path = params["statcor_file"]
        if not os.path.isabs(cor_path):
            cor_path = os.path.join(base_dir, cor_path)
        statcor = io.read_statcor_file(cor_path)

    simul = None
    if params.get("simul_file"):
        simul_path = params["simul_file"]
        if not os.path.isabs(simul_path):
            simul_path = os.path.join(base_dir, simul_path)
        simul = io.read_simul_takeoff_file(simul_path)

    results = []

    for event in events:
        # Process event
        result = process_event(
            event,
            stations,
            reversals,
            params,
            velocity_models=velocity_models,
            simul=simul,
            amp=amp_data,
            statcor=statcor,
        )
        results.append(result)

    # Write output files
    if results:
        io.write_mechanism_output(params["outfile1"], events, results)

        if params.get("outfile2"):
            for event, result in zip(events, results, strict=True):
                if result.get("success") and result.get("nout2", 0) > 0:
                    io.write_acceptable_planes(params["outfile2"], event, result)
                    break  # Only write first successful event for example

    return results


def _resolve_takeoff_tables(velocity_models):
    """Convert user-supplied velocity models into TakeoffTable instances."""
    tables = []
    for model in velocity_models:
        if isinstance(model, TakeoffTable):
            tables.append(model)
        else:
            depth, velocity = model
            tables.append(TakeoffTable(depth, velocity))
    return tables


def _trial_depths_and_models(ev_dep, sez, nmc, ntab, rng):
    """Trial source depths and velocity-model indices (hash_driver2.f).

    Trial 1 uses the unperturbed depth and model 1; every later trial
    perturbs the depth by ``sez * N(0, 1)`` and rotates the model as
    ``index(nm) = mod(nm, ntab) + 1`` (1-based), exactly like the
    original Fortran driver.
    """
    qdep2 = np.full(nmc, ev_dep, dtype=np.float64)
    if nmc > 1:
        qdep2[1:] = ev_dep + sez * rng.standard_normal(nmc - 1)
    index = np.zeros(nmc, dtype=np.int64)
    for nm in range(2, nmc + 1):
        index[nm - 1] = (nm % ntab + 1) - 1
    return qdep2, index


def _takeoff_matrix(tables, distances, ev_dep, sez, nmc, rng):
    """Per-station, per-trial takeoff angles from velocity-model tables."""
    ntab = len(tables)
    qdep2, index = _trial_depths_and_models(ev_dep, sez, nmc, ntab, rng)
    distances = np.asarray(distances, dtype=np.float64)
    p_the_mc = np.empty((distances.size, nmc), dtype=np.float64)
    for trial in range(nmc):
        p_the_mc[:, trial] = tables[index[trial]].takeoff_batch(distances, qdep2[trial])
    if np.isnan(p_the_mc).any():
        raise TakeoffRangeError(
            "takeoff lookup failed for some station/depth combinations; "
            "extend the velocity-table range (params del2/dep2) or check "
            "the event depths"
        )
    return p_the_mc


def _lookup_station_location(stations, sta):
    """(lat, lon, elev) for a station dict, or None if unknown."""
    sta_name = sta["name"]
    stations_by_name = stations.get("_byname", {})
    if sta_name in stations_by_name:
        return stations_by_name[sta_name]
    sta_net = sta.get("network", "CI")
    sta_comp = sta.get("component", "HHZ")
    key = (sta_name, sta_comp, sta_net)
    if key in stations:
        return stations[key]
    for comp in ["HHZ", "VHZ", "EHZ", "BHZ", "HH", "VH", "EH", "BH"]:
        key = (sta_name, comp, sta_net)
        if key in stations:
            return stations[key]
    # Original GETSTAT_TRI treats V and E components interchangeably.
    if sta_comp and sta_comp[0] in ("V", "E"):
        swapped = ("E" if sta_comp[0] == "V" else "V") + sta_comp[1:]
        key = (sta_name, swapped, sta_net)
        if key in stations:
            return stations[key]
    return None


def _station_distance_azimuth(station_location, ev_lat, ev_lon):
    """(distance km, azimuth deg) from the event to a station."""
    sta_lat, sta_lon, _ = station_location
    dx = (sta_lon - ev_lon) * 111.2 * math.cos(ev_lat * DEG_TO_RAD)
    dy = (sta_lat - ev_lat) * 111.2
    distance = math.sqrt(dx**2 + dy**2)
    azi = 90.0 - math.atan2(dy, dx) * RAD_TO_DEG
    if azi < 0.0:
        azi += 360.0
    return distance, azi


def _polarity_and_qual(sta, reversals, ev_date_int):
    """(polarity, quality) for a station dict; None polarity = skip."""
    pol_char = sta["polarity"].upper()
    if pol_char in ("U", "+"):
        pol = 1
    elif pol_char in ("D", "-"):
        pol = -1
    else:
        return None, None
    if sta["name"] in reversals:
        for start_date, end_date in reversals[sta["name"]]:
            if start_date <= ev_date_int <= end_date:
                pol = -pol
                break
    onset = sta.get("onset", "I").upper()
    qual = 0 if onset == "I" else 1
    return pol, qual


def _file_geometry_mc(p_azi, p_the, sazi, sthe, nmc, rng):
    """Per-station file azimuth/takeoff with per-station uncertainties.

    Trial 1 uses the file values; later trials perturb with sazi/sthe
    (original hash_driver1.f / hash_driver5.f logic).
    """
    azi = np.asarray(p_azi, dtype=np.float64).reshape(-1, 1)
    the = np.asarray(p_the, dtype=np.float64).reshape(-1, 1)
    p_azi_arr = np.repeat(azi, nmc, axis=1)
    p_the_arr = np.repeat(the, nmc, axis=1)
    if nmc > 1:
        sazi = np.asarray(sazi, dtype=np.float64).reshape(-1, 1)
        sthe = np.asarray(sthe, dtype=np.float64).reshape(-1, 1)
        z_azi = rng.standard_normal((azi.shape[0], nmc - 1))
        z_the = rng.standard_normal((azi.shape[0], nmc - 1))
        p_azi_arr[:, 1:] += sazi * z_azi
        p_the_arr[:, 1:] += sthe * z_the
    return p_azi_arr, p_the_arr


def process_event(
    event,
    stations,
    reversals,
    params,
    velocity_models=None,
    simul=None,
    amp=None,
    statcor=None,
):
    """
    Process a single event.

    Parameters
    ----------
    event : dict
        Event data (lat, lon, depth, year, month, day, stations).
        An optional ``sez`` key gives the vertical location error (km)
        used for the per-trial depth perturbation (default 1.0 km,
        like the original HASH driver).
    stations : dict
        Station locations (includes '_byname' index for fast lookup)
    reversals : dict
        Polarity reversal data
    params : dict
        HASH parameters; may carry ``sez`` and ``ratmin`` fallbacks.
    velocity_models : list, optional
        Velocity models as ``(depth, velocity)`` tuples or TakeoffTable
        instances. When given, takeoff angles come from the velocity
        tables with per-trial depth perturbation and velocity-model
        rotation (hash_driver2.f logic).
    simul : dict, optional
        SIMULPS takeoff data from read_simul_takeoff_file; when the
        event id is present, file azimuths/takeoffs with per-station
        uncertainties are used (hash_driver5.f logic).
    amp : dict, optional
        Amplitude records from read_amp_file; when the event id is
        present, S/P amplitude ratios are merged and the FOCALAMP
        pipeline is run (hash_driver3.f logic).
    statcor : dict, optional
        Station log-amplitude corrections from read_statcor_file
        (hash_driver3.f GET_COR).

    Without velocity models, SIMULPS data, or file-provided takeoff
    angles in the station records, a homogeneous geometric
    approximation is used with a warning.

    Returns
    -------
    dict
        Result from run_hash or run_hash_with_amp.
    """
    ev_lat = event["lat"]
    ev_lon = event["lon"]
    ev_dep = event["depth"]
    ev_date_int = event["year"] * 10000 + event["month"] * 100 + event["day"]
    nmc = params["nmc"]
    delmax = params["delmax"]
    amp_mode = amp is not None and str(event["id"]) in amp

    # ---- collect polarity stations ---------------------------------------
    p_azi = []
    p_pol = []
    p_qual = []
    distances = []
    file_azi = []
    file_the = []
    file_sazi = []
    file_sthe = []

    use_simul = simul is not None and str(event["id"]) in simul
    if use_simul:
        # hash_driver5.f: the SIMULPS file defines the station list and
        # geometry; polarities are matched by event id and station name
        # (first three characters).
        simul_stations = simul[str(event["id"])]
        phase_stations = {s["name"][:3]: s for s in event["stations"]}
        for name, geo in simul_stations.items():
            matched = phase_stations.get(name[:3])
            if matched is None:
                continue
            pol, qual = _polarity_and_qual(matched, reversals, ev_date_int)
            if pol is None or qual > 1:
                continue
            if geo["distance"] > delmax:
                continue
            p_azi.append(geo["azi"])
            p_pol.append(pol)
            p_qual.append(qual)
            distances.append(geo["distance"])
            file_azi.append(geo["azi"])
            file_the.append(geo["the"])
            file_sazi.append(geo["sazi"])
            file_sthe.append(geo["sthe"])
    else:
        for sta in event["stations"]:
            pol, qual = _polarity_and_qual(sta, reversals, ev_date_int)
            if pol is None:
                continue

            # hash_driver3.f uses impulsive picks only.
            if amp_mode and qual != 0:
                continue

            if "the" in sta and "azi" in sta:
                # File-provided geometry (hash_driver1.f input format).
                distance = sta.get("distance")
                if distance is None or distance > delmax:
                    continue
                p_azi.append(sta["azi"])
                file_azi.append(sta["azi"])
                file_the.append(sta["the"])
                file_sazi.append(sta.get("sazi", 0.0))
                file_sthe.append(sta.get("sthe", 0.0))
            else:
                location = _lookup_station_location(stations, sta)
                if location is None:
                    continue
                distance, azi = _station_distance_azimuth(location, ev_lat, ev_lon)
                if distance > delmax:
                    continue
                file_azi.append(None)
                file_the.append(None)
                file_sazi.append(None)
                file_sthe.append(None)
                p_azi.append(azi)

            p_pol.append(pol)
            p_qual.append(qual)
            distances.append(distance)

    # ---- merge S/P amplitude stations (hash_driver3.f) -------------------
    sp_amp = None
    if amp_mode:
        ratmin = params.get("ratmin", 3.0)
        sp_amp_values = [0.0] * len(p_pol)
        for record in amp[str(event["id"])]:
            key = (record["name"], record["component"], record["network"])
            qcor = (statcor or {}).get(key)
            if qcor is None or record["p_amp"] == 0.0:
                continue
            if abs(record["p_amp"]) / record["p_noise"] < ratmin:
                continue
            if record["s_amp"] / record["s_noise"] < ratmin:
                continue
            location = _lookup_station_location(stations, record)
            if location is None:
                continue
            distance, azi = _station_distance_azimuth(location, ev_lat, ev_lon)
            if distance > delmax:
                continue
            p_azi.append(azi)
            p_pol.append(0)
            p_qual.append(0)
            distances.append(distance)
            file_azi.append(None)
            file_the.append(None)
            file_sazi.append(None)
            file_sthe.append(None)
            sp_amp_values.append(math.log10(record["s_amp"] / abs(record["p_amp"])) - qcor)
        sp_amp = np.array(sp_amp_values, dtype=np.float64)

    if len(p_pol) == 0:
        return {
            "success": False,
            "quality": "F",
            "strike_avg": 999.0,
            "dip_avg": 99.0,
            "rake_avg": 999.0,
        }

    nsta = len(p_pol)
    distances_arr = np.asarray(distances, dtype=np.float64)
    rng = np.random.default_rng(_RANDOM_SEED)

    # ---- per-trial azimuth/takeoff geometry ------------------------------
    if velocity_models:
        # hash_driver2.f logic: takeoff angles from velocity tables,
        # per-trial source-depth perturbation and velocity-model
        # rotation; azimuths stay fixed across trials.
        tables = _resolve_takeoff_tables(velocity_models)
        sez = float(event.get("sez") or params.get("sez") or 1.0)
        p_the_mc = _takeoff_matrix(tables, distances_arr, ev_dep, sez, nmc, rng)
        p_azi_arr = np.repeat(np.asarray(p_azi).reshape(-1, 1), nmc, axis=1)
        p_the_arr = p_the_mc
    else:
        # No velocity models: file-provided geometry (hash_driver1.f /
        # hash_driver5.f) with per-station uncertainties where
        # available, geometric approximation (legacy 5 deg perturbation)
        # otherwise.
        has_file = [a is not None for a in file_azi]
        if any(has_file) and not all(has_file):
            warnings.warn(
                "mixing file-provided and computed station geometry",
                UserWarning,
                stacklevel=2,
            )
        if not all(has_file):
            warnings.warn(
                "no velocity model or file takeoff angles provided: using "
                "the homogeneous geometric takeoff approximation "
                "(degrees(atan2(distance, depth))); pass "
                "velocity_models=[(depth, velocity), ...] for velocity "
                "table takeoffs (hash_driver2.f)",
                UserWarning,
                stacklevel=2,
            )
        p_the_base = np.empty(nsta, dtype=np.float64)
        sazi_eff = np.empty(nsta, dtype=np.float64)
        sthe_eff = np.empty(nsta, dtype=np.float64)
        for i in range(nsta):
            if has_file[i]:
                p_the_base[i] = file_the[i]
                sazi_eff[i] = file_sazi[i]
                sthe_eff[i] = file_sthe[i]
            else:
                p_the_base[i] = math.degrees(math.atan2(distances_arr[i], ev_dep))
                sazi_eff[i] = 5.0
                sthe_eff[i] = 5.0
        p_azi_arr, p_the_arr = _file_geometry_mc(p_azi, p_the_base, sazi_eff, sthe_eff, nmc, rng)

    # ---- run the appropriate pipeline ------------------------------------
    common = {
        "dang": params["dang"],
        "nmc": nmc,
        "maxout": params["maxout"],
        "badfrac": params["badfrac"],
        "cangle": params["cangle"],
        "prob_max": params["prob_max"],
        "npolmin": params["npolmin"],
        "max_agap": params["max_agap"],
        "max_pgap": params["max_pgap"],
    }
    if sp_amp is not None:
        nppl = sum(1 for pol in p_pol if pol != 0)
        result = run_hash_with_amp(
            p_azi_arr,
            p_the_arr,
            np.array(p_pol),
            sp_amp,
            qbadfac=params.get("qbadfac", 0.3),
            nmismax=max(int(nppl * params["badfrac"]), 2),
            **common,
        )
    else:
        result = run_hash(p_azi_arr, p_the_arr, np.array(p_pol), np.array(p_qual), **common)

    # Add event info to result
    result["npol"] = len(p_pol)

    return result


# Export all functions
__all__ = [
    "run_hash",
    "run_hash_event_payload",
    "run_hash_with_amp",
    "run_hash_from_file",
    "process_event",
]
