"""Tests for the public API: backend selection, thread control, batch."""

import numpy as np
import pytest
from conftest import make_amp_event, make_synthetic_event

from cnchash import (
    available_backends,
    get_backend_info,
    get_num_threads,
    run_hash,
    run_hash_batch,
    run_hash_batch_with_amp,
    run_hash_with_amp,
    set_num_threads,
)


def test_available_backends_contains_numba():
    assert "numba" in available_backends()


def test_backend_info_shape():
    info = get_backend_info()
    assert info["backend"] in ("fortran", "numba")
    assert info["precision"] == "float64"
    assert info["threads"] >= 1


def test_set_get_num_threads():
    n0 = get_num_threads()
    set_num_threads(max(1, n0 - 1) if n0 > 1 else 1)
    assert get_num_threads() >= 1
    set_num_threads(n0)


def test_unknown_backend_raises():
    with pytest.raises(ValueError):
        run_hash(*make_synthetic_event(nsta=12), backend="nonexistent")


def test_run_hash_basic():
    az, the, pol, qual = make_synthetic_event(nsta=20)
    result = run_hash(az, the, pol, qual, nmc=30)
    assert result["success"] is True
    assert result["nf"] > 0
    assert result["nmult"] >= 1
    q = result["quality"]
    assert all(x in "ABCD" for x in (q if isinstance(q, list) else [q]))
    assert 0.0 <= result["strike_avg"] < 360.0
    assert 0.0 <= result["dip_avg"] <= 90.0
    assert -180.0 <= result["rake_avg"] <= 180.0
    assert len(result["strike"]) == result["nf"]


def test_run_hash_backends_agree():
    az, the, pol, qual = make_synthetic_event(nsta=25)
    results = {}
    for be in available_backends():
        results[be] = run_hash(az, the, pol, qual, nmc=30, backend=be, num_threads=1)
    keys = list(results)
    base = results[keys[0]]
    for be in keys[1:]:
        r = results[be]
        assert r["success"] == base["success"]
        assert r["nf"] == base["nf"]
        assert r["nmult"] == base["nmult"]
        assert r["quality"] == base["quality"]
        assert abs(r["mfrac"] - base["mfrac"]) < 1e-3


def test_run_hash_min_polarities_fails():
    az, the, pol, qual = make_synthetic_event(nsta=5)
    result = run_hash(az, the, pol, qual)
    assert result["success"] is False
    assert result["quality"] == "F"


def test_run_hash_with_amp():
    az, the, pol, sp = make_amp_event(nsta=22)
    result = run_hash_with_amp(az, the, pol, sp, nmc=30)
    assert result["success"] is True
    assert result["nf"] > 0
    assert 0.0 <= result["mavg"] < 5.0


def test_run_hash_batch_matches_scalar():
    events = [make_synthetic_event(nsta=int(n), seed=100 + i) for i, n in enumerate((15, 20, 25))]
    batch = run_hash_batch(events, nmc=30, num_threads=1)
    assert len(batch) == 3
    for (az, the, pol, qual), r in zip(events, batch, strict=True):
        scalar = run_hash(az, the, pol, qual, nmc=30, num_threads=1)
        assert r["nf"] == scalar["nf"]
        assert np.isclose(np.atleast_1d(r["strike_avg"])[0], np.atleast_1d(scalar["strike_avg"])[0], atol=1e-6)
        assert np.isclose(r["mfrac"], scalar["mfrac"], atol=1e-12)


def test_run_hash_batch_with_amp():
    events = [make_amp_event(nsta=18, seed=200 + i) for i in range(3)]
    batch = run_hash_batch_with_amp(events, nmc=30, num_threads=1)
    assert len(batch) == 3
    for (az, the, pol, sp), r in zip(events, batch, strict=True):
        scalar = run_hash_with_amp(az, the, pol, sp, nmc=30, num_threads=1)
        assert r["nf"] == scalar["nf"]
        assert np.isclose(np.atleast_1d(r["strike_avg"])[0], np.atleast_1d(scalar["strike_avg"])[0], atol=1e-6)
        assert np.isclose(r["mavg"], scalar["mavg"], atol=1e-12)


def test_backend_enumeration_consistent():
    """Every registered backend must implement the full contract."""
    az, the, pol, qual = make_synthetic_event(nsta=12)
    for be in available_backends():
        solver = __import__("cnchash").backend.get_backend(be)
        gap = solver.get_gap(len(pol), az, the)
        assert len(gap) == 2
        mfrac, stdr = solver.get_misfit(len(pol), az, the, pol, qual, 45.0, 60.0, 90.0)
        assert 0.0 <= mfrac <= 1.0
        assert 0.0 <= stdr <= 1.0
