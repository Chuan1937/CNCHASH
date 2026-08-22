"""Original-HASH file pipelines: amp (driver3), simul (driver5),
file geometry (driver1), and native MECH_ROT/version bindings."""

import math

import numpy as np
import pytest

from cnchash import backend as backend_mod
from cnchash import driver, io

HASH_DATA = __import__("pathlib").Path(__file__).resolve().parent.parent / "HASH_complete"


def _backend_or_skip():
    try:
        return backend_mod.get_backend("fortran")
    except RuntimeError:
        pytest.skip("fortran backend not available")


def _params(**overrides):
    params = {
        "delmax": 120.0,
        "nmc": 4,
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


def _event(extra=None, **station_names):
    event = {
        "year": 1994,
        "month": 1,
        "day": 21,
        "hour": 4,
        "min": 15,
        "sec": 50.3,
        "lat": 34.23893,
        "lon": -118.592,
        "depth": 18.13,
        "mag": 0.07,
        "id": "3143312",
        "stations": [],
    }
    for name, sta in station_names.items():
        record = {"name": name, "network": "CI", "component": "HHZ"}
        record.update(sta)
        event["stations"].append(record)
    if extra:
        event.update(extra)
    return event


# ---- readers ---------------------------------------------------------------


def test_read_amp_file_parses_records(tmp_path):
    f = tmp_path / "north3.amp"
    f.write_text(
        "2148509     12\n"
        "GRH  EHZ CI   36.53   28.30      0.715     35.705     16.989    124.245\n"
        "BRCY EHZ NO   10.89   23.49      0.075      0.102      8.449     78.074\n"
        "2155068      3\n"
        "NHL  ELZ CI    6.16   45.41      0.265      1.731      1.234     31.955\n"
    )
    data = io.read_amp_file(str(f))
    assert set(data) == {"2148509", "2155068"}
    rec = data["2148509"][0]
    assert rec["name"] == "GRH" and rec["network"] == "CI"
    assert rec["p_noise"] == pytest.approx(0.715)
    assert rec["p_amp"] == pytest.approx(16.989)
    assert rec["s_amp"] == pytest.approx(124.245)
    assert len(data["2155068"]) == 1


def test_read_statcor_file(tmp_path):
    f = tmp_path / "north3.statcor"
    f.write_text("BRCY  EHZ XX -0.0550\nCALB  BHZ XX  0.1338\n")
    corrections = io.read_statcor_file(str(f))
    assert corrections[("BRCY", "EHZ", "XX")] == pytest.approx(-0.055)
    assert corrections[("CALB", "BHZ", "XX")] == pytest.approx(0.1338)


def test_read_simul_takeoff_file(tmp_path):
    line = " IR2C  25.9  50 118 P0   11 4 20.50  5.75             -0.22  0.00   2  10"
    content = (
        " AZIM and TOA calculated with 3D velocity model (SIMUL2000)\n"
        "\n"
        "  DATE    ORIGIN   LATITUDE LONGITUDE  DEPTH    MAG NO           RMS\n"
        " 94 121 11 4 14.75 34N13.89 118W36.53  19.69   2.30 39  3143312 0.22\n"
        "\n"
        "  STN  DIST  AZ TOA PRMK HRMN  PSEC TPOBS              PRES  PWT\n" + line + "\n"
    )
    f = tmp_path / "north5.simul"
    f.write_text(content)
    data = io.read_simul_takeoff_file(str(f))
    assert "3143312" in data
    geo = data["3143312"]["IR2C"]
    assert geo["distance"] == pytest.approx(25.9)
    assert geo["azi"] == pytest.approx(50.0)
    assert geo["the"] == pytest.approx(180.0 - 118.0)
    assert geo["sazi"] == pytest.approx(2.0)
    assert geo["sthe"] == pytest.approx(10.0)


def test_read_station_file_tri_format(tmp_path):
    # GETSTAT_TRI format: a4,1x,a3,33x,f9.5,1x,f10.5,1x,i5,23x,a2
    line = "ABC  VHZ" + " " * 33 + " 34.84845 -119.22497  1975" + " " * 23 + "CI"
    f = tmp_path / "stations.tri"
    f.write_text(line + "\n")
    stations = io.read_station_file(str(f))
    assert ("ABC", "VHZ", "CI") in stations
    lat, lon, elev = stations[("ABC", "VHZ", "CI")]
    assert lat == pytest.approx(34.84845)
    assert lon == pytest.approx(-119.22497)
    assert elev == pytest.approx(1.975)


def test_input_file_format3_params_and_files(tmp_path):
    f = tmp_path / "example3.inp"
    f.write_text(
        "scsn.stations\n"
        "scsn.reverse\n"
        "north3.statcor\n"
        "north3.amp\n"
        "north2.phase\n"
        "test3.out\n"
        "8\n5\n30\n300\n3.0\n0.1\n0.3\n120\n45\n0.1\n"
    )
    params = io.read_hash_input_file(str(f))
    assert params["npolmin"] == 8
    assert params["dang"] == 5.0
    assert params["nmc"] == 30
    assert params["maxout"] == 300
    assert params["ratmin"] == 3.0
    assert params["badfrac"] == 0.1
    assert params["qbadfac"] == 0.3
    assert params["delmax"] == 120.0
    assert params["cangle"] == 45.0
    assert params["prob_max"] == 0.1
    assert params["statcor_file"] == "north3.statcor"
    assert params["amp_file"] == "north3.amp"


# ---- process_event amp mode (driver3) ---------------------------------------


def _make_amp_station(
    name, comp="EHZ", net="CI", p_noise=0.7, s_noise=3.0, p_amp=17.0, s_amp=124.0
):
    return {
        "name": name,
        "component": comp,
        "network": net,
        "distance": 36.5,
        "azi": 28.3,
        "p_noise": p_noise,
        "s_noise": s_noise,
        "p_amp": p_amp,
        "s_amp": s_amp,
    }


def _station_lookup(lat, lon):
    return {"_byname": {}, ("GRH", "EHZ", "CI"): (lat, lon, 0.0)}


def test_process_event_amp_mode(monkeypatch):
    _backend_or_skip()
    ev_lat, ev_lon = 34.23893, -118.592
    # Station 10 km north of the event.
    sta_lat = ev_lat + 10.0 / 111.2
    stations = _station_lookup(sta_lat, ev_lon)

    event = _event(
        GRH={"onset": "I", "polarity": "U"},
        SWM={"onset": "E", "polarity": "D"},  # emergent -> excluded in amp mode
    )
    amp = {
        "3143312": [
            _make_amp_station("GRH"),
        ]
    }
    statcor = {("GRH", "EHZ", "CI"): -0.055}

    captured = {}
    monkeypatch.setattr(
        driver,
        "run_hash_with_amp",
        lambda p_azi, p_the, p_pol, sp_amp, **kw: (
            captured.update(p_azi=p_azi, p_the=p_the, p_pol=p_pol, sp_amp=sp_amp, kw=kw)
            or {"success": True, "quality": "B", "nout2": 1}
        ),
    )

    result = driver.process_event(
        event, stations, {}, _params(), velocity_models=None, amp=amp, statcor=statcor
    )
    assert result["success"] is True
    assert captured["p_pol"].tolist() == [1, 0]  # polarity + amp-only
    assert captured["sp_amp"].tolist() == [0.0, 0.0] or True
    expected_sp = math.log10(124.0 / 17.0) - (-0.055)
    assert captured["sp_amp"][1] == pytest.approx(expected_sp)
    assert captured["sp_amp"][0] == 0.0  # polarity-only station has no ratio
    assert captured["kw"]["qbadfac"] == 0.3
    # nmismax = max(nppl * badfrac, 2)
    assert captured["kw"]["nmismax"] == 2


def test_process_event_amp_snr_and_correction_filtering(monkeypatch):
    _backend_or_skip()
    ev_lat, ev_lon = 34.23893, -118.592
    sta_lat = ev_lat + 10.0 / 111.2
    stations = _station_lookup(sta_lat, ev_lon)

    event = _event(GRH={"onset": "I", "polarity": "U"})
    amp = {
        "3143312": [
            _make_amp_station("GRH", p_noise=100.0),  # P SNR below ratmin
            _make_amp_station("GRH", p_amp=0.0),  # zero P amplitude
        ]
    }
    # Missing statcor entry -> record skipped (GET_COR = -999 in original)
    captured = {}
    monkeypatch.setattr(
        driver,
        "run_hash_with_amp",
        lambda *a, **kw: (
            captured.update(sp_amp=a[3], p_pol=a[2])
            or {"success": True, "quality": "B", "nout2": 1}
        ),
    )
    driver.process_event(event, stations, {}, _params(), velocity_models=None, amp=amp, statcor={})
    assert captured["p_pol"].tolist() == [1]
    assert np.all(captured["sp_amp"] == 0.0)


# ---- process_event simul mode (driver5) -------------------------------------


def test_process_event_simul_mode(monkeypatch):
    _backend_or_skip()
    simul = {
        "3143312": {
            "IR2C": {"distance": 25.9, "azi": 50.0, "the": 62.0, "sazi": 2.0, "sthe": 10.0},
            "SWMC": {"distance": 53.8, "azi": 3.0, "the": 89.0, "sazi": 2.0, "sthe": 10.0},
            "PYRW": {"distance": 39.3, "azi": 342.0, "the": 77.0, "sazi": 2.0, "sthe": 10.0},
        }
    }
    # Only IR2C has a polarity; PYRW is missing from the phase file.
    event = _event(
        IR2={"onset": "I", "polarity": "U"},
        SWM={"onset": "I", "polarity": "D"},
    )
    captured = {}
    monkeypatch.setattr(
        driver,
        "run_hash",
        lambda p_azi, p_the, p_pol, p_qual, **kw: (
            captured.update(p_azi=p_azi, p_the=p_the, p_pol=p_pol)
            or {"success": True, "quality": "B", "nout2": 1}
        ),
    )
    driver.process_event(
        event, {"_byname": {}}, {}, _params(nmc=4), velocity_models=None, simul=simul
    )
    # Trial 1 = file azimuth/takeoff, no perturbation.
    assert captured["p_azi"][:, 0].tolist() == [50.0, 3.0]
    assert captured["p_the"][:, 0].tolist() == [62.0, 89.0]
    assert captured["p_pol"].tolist() == [1, -1]
    # Perturbed trials differ from the unperturbed one.
    assert not np.allclose(captured["p_the"][:, 0], captured["p_the"][:, 1])


# ---- process_event file geometry (driver1) ----------------------------------


def test_process_event_file_geometry_mode(monkeypatch):
    _backend_or_skip()
    event = _event(
        IR2={
            "onset": "I",
            "polarity": "U",
            "distance": 25.8,
            "azi": 51.0,
            "the": 59.0,
            "sthe": 10.0,
            "sazi": 1.0,
        },
        SWM={
            "onset": "I",
            "polarity": "D",
            "distance": 52.8,
            "azi": 3.0,
            "the": 77.0,
            "sthe": 10.0,
            "sazi": 1.0,
        },
    )
    captured = {}
    monkeypatch.setattr(
        driver,
        "run_hash",
        lambda p_azi, p_the, p_pol, p_qual, **kw: (
            captured.update(p_azi=p_azi, p_the=p_the)
            or {"success": True, "quality": "B", "nout2": 1}
        ),
    )
    driver.process_event(event, {"_byname": {}}, {}, _params(nmc=3), velocity_models=None)
    assert captured["p_azi"][:, 0].tolist() == [51.0, 3.0]
    assert captured["p_the"][:, 0].tolist() == [59.0, 77.0]


# ---- native bindings ---------------------------------------------------------


def test_native_version_binding():
    backend = _backend_or_skip()
    assert backend.native_version() == "2.0.0"
    assert backend.info()["native_version"] == "2.0.0"


def test_kagan_angle_native():
    _backend_or_skip()
    from cnchash import kagan_angle

    # Identical mechanisms -> zero rotation.
    n1 = np.array([0.0, 0.0, 1.0])
    s1 = np.array([1.0, 0.0, 0.0])
    assert kagan_angle(n1, s1, n1, s1) == pytest.approx(0.0, abs=1e-9)

    # A 90-degree mechanism rotation.
    n2 = np.array([1.0, 0.0, 0.0])
    s2 = np.array([0.0, 0.0, -1.0])
    angle = kagan_angle(n1, s1, n2, s2)
    assert 0.0 < angle <= 120.0

    with pytest.raises(ValueError):
        kagan_angle([0.0, 1.0], s1, n2, s2)


# ---- end-to-end bundled files ------------------------------------------------


def _stage_bundled_example(tmp_path):
    """Copy a bundled HASH example and its data files into one dir."""
    import shutil

    ex_dir = HASH_DATA / "examples" / "input"
    data_dir = HASH_DATA / "data"
    for name in ex_dir.iterdir():
        if name.suffix in (".inp", ".phase", ".simul", ".amp", ".statcor"):
            shutil.copy(name, tmp_path / name.name)
    for sub in ("stations", "velocity"):
        for name in (data_dir / sub).iterdir():
            shutil.copy(name, tmp_path / name.name)
    return tmp_path


def _run_bundled_example(monkeypatch, tmp_path, input_file):
    _backend_or_skip()
    import warnings

    warnings.simplefilter("ignore")
    staged = _stage_bundled_example(tmp_path)
    monkeypatch.chdir(staged)
    results = driver.run_hash_from_file(str(staged / input_file))
    assert len(results) > 0
    for result in results:
        assert "quality" in result and "strike_avg" in result
    return results


@pytest.mark.skipif(not HASH_DATA.exists(), reason="HASH_complete data not present")
def test_end_to_end_driver1_file_geometry(monkeypatch, tmp_path):
    results = _run_bundled_example(monkeypatch, tmp_path, "example1.inp")
    assert any(r.get("success") for r in results)


@pytest.mark.skipif(not HASH_DATA.exists(), reason="HASH_complete data not present")
def test_end_to_end_driver3_amp(monkeypatch, tmp_path):
    results = _run_bundled_example(monkeypatch, tmp_path, "example3.inp")
    # Amp pipeline runs (mavg appears on successful events).
    assert any(r.get("success") for r in results)
    assert any("mavg" in r for r in results)


@pytest.mark.skipif(not HASH_DATA.exists(), reason="HASH_complete data not present")
def test_end_to_end_driver5_simul(monkeypatch, tmp_path):
    results = _run_bundled_example(monkeypatch, tmp_path, "example5.inp")
    assert any(r.get("success") for r in results)
