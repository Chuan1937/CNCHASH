"""Parity tests: the Fortran backend must match the Numba reference.

Tolerances:
- acceptable mechanism sets: exact (deterministic grid search)
- preferred mechanism angles: within 1 degree (cluster averaging is
  sensitive to float round-off when many mechanisms are near-identical)
- misfit statistics: within 1e-3
"""

import numpy as np
import pytest
from conftest import make_amp_event, make_synthetic_event

from cnchash import backend as backend_mod

ANGLE_TOL = 1.0
MISFIT_TOL = 1e-3


@pytest.fixture(scope="module")
def fortran_backend():
    if "fortran" not in backend_mod.available_backends():
        pytest.skip("fortran backend not available")
    return backend_mod.get_backend("fortran")


@pytest.fixture(scope="module")
def numba_backend():
    return backend_mod.get_backend("numba")


def _scalar(v):
    return float(np.atleast_1d(v)[0])


def _assert_event_parity(rn, rf):
    assert rn["success"] == rf["success"]
    assert rn["nf"] == rf["nf"]
    assert rn["nmult"] == rf["nmult"]
    assert rn["quality"] == rf["quality"]
    assert abs(_scalar(rn["strike_avg"]) - _scalar(rf["strike_avg"])) < ANGLE_TOL
    assert abs(_scalar(rn["dip_avg"]) - _scalar(rf["dip_avg"])) < ANGLE_TOL
    assert abs(_scalar(rn["rake_avg"]) - _scalar(rf["rake_avg"])) < ANGLE_TOL
    assert abs(rn["mfrac"] - rf["mfrac"]) < MISFIT_TOL
    assert abs(rn["stdr"] - rf["stdr"]) < MISFIT_TOL


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_focalmc_mechanism_sets_identical(numba_backend, fortran_backend, seed):
    """The grid search must return identical acceptable mechanisms.

    The fault/slip vectors are compared exactly. Strike/dip/rake angles
    are compared with a tolerance because FPCOOR's strike is numerically
    unstable for near-vertical faults (|fn(3)| within 1e-7 of 1), where
    tiny grid round-off differences can switch the branch.
    """
    az, the, pol, qual = make_synthetic_event(nsta=25, seed=seed)
    azi_mc = az.reshape(-1, 1)
    the_mc = the.reshape(-1, 1)
    n = len(pol)

    rn = numba_backend.run_event(azi_mc, the_mc, pol, qual, n, 1, 5.0, 500, 0.1, 45.0, 0.1, 8, 90, 60, 0)
    rf = fortran_backend.run_event(azi_mc, the_mc, pol, qual, n, 1, 5.0, 500, 0.1, 45.0, 0.1, 8, 90, 60, 0)

    assert rn["nf"] == rf["nf"]
    assert np.allclose(rn["faults"], rf["faults"], atol=1e-12)
    assert np.allclose(rn["slips"], rf["slips"], atol=1e-12)

    # Compare angles with a tolerance; near-vertical faults get extra slack
    # on strike (undefined for horizontal mechanisms).
    diff_s = np.abs(rn["strike"] - rf["strike"]) % 360.0
    diff_s = np.minimum(diff_s, 360.0 - diff_s)
    diff_d = np.abs(rn["dip"] - rf["dip"])
    diff_r = np.abs(rn["rake"] - rf["rake"]) % 360.0
    diff_r = np.minimum(diff_r, 360.0 - diff_r)
    vertical = np.abs(np.abs(rn["faults"][2]) - 1.0) < 1e-6
    tol = np.where(vertical, 10.0, 1e-6)
    assert np.all(diff_s < tol)
    assert np.all(diff_d < 1e-6)
    assert np.all(diff_r < np.where(vertical, 10.0, 1e-6))


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_full_pipeline_parity(numba_backend, fortran_backend, seed):
    az, the, pol, qual = make_synthetic_event(nsta=25, seed=seed)
    rn = numba_backend.run_event(az.reshape(-1, 1), the.reshape(-1, 1), pol, qual,
                                 len(pol), 1, 5.0, 500, 0.1, 45.0, 0.1, 8, 90, 60, 0)
    rf = fortran_backend.run_event(az.reshape(-1, 1), the.reshape(-1, 1), pol, qual,
                                   len(pol), 1, 5.0, 500, 0.1, 45.0, 0.1, 8, 90, 60, 0)
    _assert_event_parity(rn, rf)


@pytest.mark.parametrize("nmc", [1, 30])
def test_full_pipeline_parity_multitrial(numba_backend, fortran_backend, nmc):
    az, the, pol, qual = make_synthetic_event(nsta=20, seed=9)
    azi_mc = np.repeat(az.reshape(-1, 1), nmc, axis=1)
    the_mc = np.repeat(the.reshape(-1, 1), nmc, axis=1)
    rn = numba_backend.run_event(azi_mc, the_mc, pol, qual, len(pol), nmc,
                                 5.0, 500, 0.1, 45.0, 0.1, 8, 90, 60, 0)
    rf = fortran_backend.run_event(azi_mc, the_mc, pol, qual, len(pol), nmc,
                                   5.0, 500, 0.1, 45.0, 0.1, 8, 90, 60, 0)
    _assert_event_parity(rn, rf)


def test_amp_parity(numba_backend, fortran_backend):
    az, the, pol, sp = make_amp_event(nsta=22, seed=11)
    rn = numba_backend.run_event_amp(az.reshape(-1, 1), the.reshape(-1, 1), pol, sp,
                                     len(pol), 1, 5.0, 500, 0.1, 0.3, 45.0, 0.1, 8, 90, 60, 0)
    rf = fortran_backend.run_event_amp(az.reshape(-1, 1), the.reshape(-1, 1), pol, sp,
                                       len(pol), 1, 5.0, 500, 0.1, 0.3, 45.0, 0.1, 8, 90, 60, 0)
    assert rn["success"] == rf["success"]
    assert rn["nf"] == rf["nf"]
    assert rn["nmult"] == rf["nmult"]
    assert abs(_scalar(rn["strike_avg"]) - _scalar(rf["strike_avg"])) < ANGLE_TOL
    assert abs(rn["mfrac"] - rf["mfrac"]) < MISFIT_TOL
    assert abs(rn["mavg"] - rf["mavg"]) < MISFIT_TOL
    assert abs(rn["stdr"] - rf["stdr"]) < MISFIT_TOL


def test_gap_parity(numba_backend, fortran_backend):
    az, the, pol, qual = make_synthetic_event(nsta=30, seed=5)
    n = len(pol)
    gnum = numba_backend.get_gap(n, az, the)
    gfor = fortran_backend.get_gap(n, az, the)
    assert np.allclose(gnum, gfor, atol=1e-9)


def test_mech_rot_parity(numba_backend, fortran_backend):
    import ctypes

    from cnchash import uncertainty

    lib = fortran_backend._lib
    lib.cnchash_mech_rot.argtypes = [
        np.ctypeslib.ndpointer(np.float64),
        np.ctypeslib.ndpointer(np.float64),
        np.ctypeslib.ndpointer(np.float64),
        np.ctypeslib.ndpointer(np.float64),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.cnchash_mech_rot.restype = None

    rng = np.random.default_rng(0)
    for _ in range(50):
        n1 = rng.normal(size=3)
        s1 = rng.normal(size=3)
        n2 = rng.normal(size=3)
        s2 = rng.normal(size=3)
        n1 /= np.linalg.norm(n1)
        s1 /= np.linalg.norm(s1)
        n2 /= np.linalg.norm(n2)
        s2 /= np.linalg.norm(s2)
        rn = uncertainty.mech_rot(n1, s1, n2.copy(), s2.copy())
        ro = ctypes.c_double()
        lib.cnchash_mech_rot(n1, s1, n2, s2, ctypes.byref(ro))
        assert abs(rn - ro.value) < 1e-9
