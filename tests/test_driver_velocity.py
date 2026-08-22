"""process_event velocity-model integration (hash_driver2.f logic)."""

import math

import numpy as np
import pytest

from cnchash import TakeoffTable, driver
from cnchash import backend as backend_mod

UNIFORM_MODEL = (np.array([0.0, 35.0]), np.array([5.5, 5.5]))
GRADIENT_MODEL = (np.array([0.0, 5.0, 35.0]), np.array([5.0, 8.0, 8.0]))


def _backend_or_skip():
    try:
        return backend_mod.get_backend("fortran")
    except RuntimeError:
        pytest.skip("fortran backend not available")


@pytest.fixture(scope="module")
def uniform_table():
    _backend_or_skip()
    return TakeoffTable(*UNIFORM_MODEL)


@pytest.fixture(scope="module")
def gradient_table():
    _backend_or_skip()
    return TakeoffTable(*GRADIENT_MODEL)


def _make_event(distances, azimuths, depth=10.0):
    """Synthetic event whose stations sit at given (distance, azimuth)."""
    ev_lat, ev_lon = 34.0, -118.0
    stations = []
    for dist, azi in zip(distances, azimuths, strict=True):
        dx = dist * math.sin(math.radians(azi))
        dy = dist * math.cos(math.radians(azi))
        sta_lat = ev_lat + dy / 111.2
        sta_lon = ev_lon + dx / (111.2 * math.cos(math.radians(ev_lat)))
        stations.append(
            {
                "name": f"S{azi:.0f}_{dist:.0f}",
                "network": "CI",
                "component": "HHZ",
                "onset": "I",
                "polarity": "U",
            }
        )
        _BYNAME[f"S{azi:.0f}_{dist:.0f}"] = (sta_lat, sta_lon, 0.0)
    event = {
        "year": 2020,
        "month": 1,
        "day": 2,
        "hour": 3,
        "min": 4,
        "sec": 5.0,
        "lat": ev_lat,
        "lon": ev_lon,
        "depth": depth,
        "mag": 2.0,
        "id": "1",
        "stations": stations,
    }
    return event


_BYNAME = {}


@pytest.fixture(autouse=True)
def _stations():
    _BYNAME.clear()
    return {"_byname": _BYNAME}


def _params(**overrides):
    params = {
        "delmax": 120.0,
        "nmc": 6,
        "dang": 5.0,
        "maxout": 500,
        "badfrac": 0.1,
        "cangle": 45.0,
        "prob_max": 0.1,
        "npolmin": 8,
        "max_agap": 90.0,
        "max_pgap": 60.0,
    }
    params.update(overrides)
    return params


def _capture_run(monkeypatch):
    captured = {}

    def fake_run_hash(p_azi_arr, p_the_arr, p_pol, p_qual, **kwargs):
        captured["p_azi"] = np.asarray(p_azi_arr).copy()
        captured["p_the"] = np.asarray(p_the_arr).copy()
        captured["kwargs"] = kwargs
        return {"success": True, "quality": "B", "nout2": 1}

    monkeypatch.setattr(driver, "run_hash", fake_run_hash)
    return captured


def _expected_trial_depths(depth, sez, nmc):
    rng = np.random.default_rng(driver._RANDOM_SEED)
    qdep2 = np.full(nmc, depth)
    if nmc > 1:
        qdep2[1:] = depth + sez * rng.standard_normal(nmc - 1)
    return qdep2


def test_velocity_branch_matches_driver2_semantics(monkeypatch, uniform_table, gradient_table):
    event = _make_event([10.0, 20.0, 30.0, 40.0], [0.0, 90.0, 180.0, 270.0])
    params = _params(nmc=6)
    captured = _capture_run(monkeypatch)

    tables = [uniform_table, gradient_table]
    result = driver.process_event(event, {"_byname": _BYNAME}, {}, params, velocity_models=tables)

    assert result["success"] is True
    p_the = captured["p_the"]
    p_azi = captured["p_azi"]
    assert p_the.shape == (4, 6)
    assert p_azi.shape == (4, 6)

    # Azimuths stay fixed across trials (original HASH behavior).
    for col in range(6):
        np.testing.assert_allclose(p_azi[:, col], p_azi[:, 0])

    # Trial 1: unperturbed depth, first velocity model.
    expected0 = uniform_table.takeoff_batch([10.0, 20.0, 30.0, 40.0], 10.0)
    np.testing.assert_allclose(p_the[:, 0], expected0, atol=0.1)

    # Model rotation: index(nm) = mod(nm, ntab) + 1 (1-based), trial 1 = 1.
    sez = 1.0  # event/params carry no sez -> default 1.0 km
    qdep2 = _expected_trial_depths(10.0, sez, 6)
    for trial in range(1, 6):
        nm = trial + 1
        model = tables[(nm % 2 + 1) - 1]
        expected = model.takeoff_batch([10.0, 20.0, 30.0, 40.0], qdep2[trial])
        np.testing.assert_allclose(p_the[:, trial], expected, atol=0.1)


def test_velocity_branch_depth_perturbation_and_sez(monkeypatch, uniform_table):
    event = _make_event([10.0, 20.0], [0.0, 90.0])
    event["sez"] = 3.0
    params = _params(nmc=6)
    captured = _capture_run(monkeypatch)

    driver.process_event(event, {"_byname": _BYNAME}, {}, params, velocity_models=[uniform_table])

    p_the = captured["p_the"]
    qdep2 = _expected_trial_depths(10.0, 3.0, 6)
    for trial in range(6):
        expected = uniform_table.takeoff_batch([10.0, 20.0], qdep2[trial])
        np.testing.assert_allclose(p_the[:, trial], expected, atol=0.1)

    # Perturbed trials differ from the unperturbed one.
    assert not np.allclose(p_the[:, 0], p_the[:, 1], atol=0.1)


def test_velocity_branch_deterministic(monkeypatch, uniform_table):
    event = _make_event([10.0, 20.0], [0.0, 90.0])
    params = _params(nmc=6)

    captured1 = _capture_run(monkeypatch)
    driver.process_event(event, {"_byname": _BYNAME}, {}, params, velocity_models=[uniform_table])

    captured2 = _capture_run(monkeypatch)
    driver.process_event(event, {"_byname": _BYNAME}, {}, params, velocity_models=[uniform_table])

    np.testing.assert_array_equal(captured1["p_the"], captured2["p_the"])


def test_geometric_fallback_warns_and_uses_hash_convention(monkeypatch):
    event = _make_event([10.0, 20.0], [0.0, 90.0], depth=10.0)
    params = _params(nmc=3)
    captured = _capture_run(monkeypatch)

    with pytest.warns(UserWarning, match="velocity"):
        driver.process_event(event, {"_byname": _BYNAME}, {}, params, velocity_models=None)

    p_the = captured["p_the"]
    # HASH convention: 0 = up, so the homogeneous straight-ray takeoff
    # is degrees(atan2(distance, depth)).
    np.testing.assert_allclose(
        p_the[:, 0],
        [np.degrees(np.arctan2(10.0, 10.0)), np.degrees(np.arctan2(20.0, 10.0))],
        atol=1e-9,
    )


def test_no_model_same_as_empty_list(monkeypatch):
    event = _make_event([10.0], [0.0])
    params = _params(nmc=3)

    captured = _capture_run(monkeypatch)
    with pytest.warns(UserWarning):
        driver.process_event(event, {"_byname": _BYNAME}, {}, params, velocity_models=[])
    empty = captured["p_the"].copy()

    captured = _capture_run(monkeypatch)
    with pytest.warns(UserWarning):
        driver.process_event(event, {"_byname": _BYNAME}, {}, params, velocity_models=None)
    none = captured["p_the"].copy()

    np.testing.assert_array_equal(empty, none)


def test_out_of_range_trial_depth_raises(monkeypatch, uniform_table):
    # Depth beyond the default table range (dep2=35) must fail loudly.
    event = _make_event([10.0], [0.0], depth=40.0)
    params = _params(nmc=3)
    _capture_run(monkeypatch)

    with pytest.raises(driver.TakeoffRangeError):
        driver.process_event(
            event, {"_byname": _BYNAME}, {}, params, velocity_models=[uniform_table]
        )


def test_run_hash_from_file_loads_velocity_models(monkeypatch, tmp_path):
    _backend_or_skip()
    (tmp_path / "vz.a").write_text("0.0 5.0\n20.0 6.0\n")
    (tmp_path / "vz.b").write_text("0.0 5.5\n20.0 5.5\n")
    inp = tmp_path / "example.inp"
    inp.write_text(
        "scsn.stations\n"
        "scsn.reverse\n"
        "north2.phase\n"
        "out1.txt\n"
        "out2.txt\n"
        "8\n90\n60\n5\n30\n500\n0.1\n120\n45\n0.25\n5\n"
        "vz.a\n"
        "vz.b\n"
    )

    from cnchash import io as io_mod

    monkeypatch.setattr(io_mod, "read_phase_file", lambda path: [{"dummy": True}])
    monkeypatch.setattr(io_mod, "read_station_file", lambda path: {"_byname": {}})
    monkeypatch.setattr(io_mod, "read_polarity_reversal_file", lambda path: {})
    monkeypatch.setattr(io_mod, "write_mechanism_output", lambda *a, **k: None)
    monkeypatch.setattr(io_mod, "write_acceptable_planes", lambda *a, **k: None)

    calls = []
    monkeypatch.setattr(
        driver,
        "process_event",
        lambda event, stations, reversals, params, velocity_models=None, **kwargs: (
            calls.append(velocity_models) or {"success": True, "nout2": 1}
        ),
    )

    driver.run_hash_from_file(str(inp))

    assert len(calls) == 1
    models = calls[0]
    assert models is not None and len(models) == 2
    for (depth, velocity), expected in zip(
        models,
        [
            (np.array([0.0, 20.0]), np.array([5.0, 6.0])),
            (np.array([0.0, 20.0]), np.array([5.5, 5.5])),
        ],
        strict=True,
    ):
        np.testing.assert_array_equal(depth, expected[0])
        np.testing.assert_array_equal(velocity, expected[1])


def test_input_file_params_positional_order(tmp_path):
    from cnchash import io as io_mod

    inp = tmp_path / "example2.inp"
    inp.write_text(
        "scsn.stations\n"
        "scsn.reverse\n"
        "north2.phase\n"
        "out1.txt\n"
        "out2.txt\n"
        "8\n90\n60\n5\n30\n500\n0.1\n120\n45\n0.25\n5\n"
        "vz.socal\n"
        "vz.north\n"
    )
    params = io_mod.read_hash_input_file(str(inp))
    # hash_driver2.f parameter order: npolmin, max_agap, max_pgap, dang,
    # nmc, maxout, badfrac, delmax, cangle, prob_max.
    assert params["npolmin"] == 8
    assert params["max_agap"] == 90.0
    assert params["max_pgap"] == 60.0
    assert params["dang"] == 5.0
    assert params["nmc"] == 30
    assert params["maxout"] == 500
    assert params["badfrac"] == 0.1
    assert params["delmax"] == 120.0
    assert params["cangle"] == 45.0
    assert params["prob_max"] == 0.25
    assert params["velocity_models"] == ["vz.socal", "vz.north"]
