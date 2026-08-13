"""Shared fixtures and helpers for the CNCHASH test suite."""

import numpy as np
import pytest

from cnchash import backend as backend_mod

BACKEND = "fortran"


def make_synthetic_event(nsta=25, seed=42, mechanism=(45.0, 60.0, 90.0)):
    """Build a realistic synthetic event with a known mechanism.

    Returns (az, the, pol, qual) with ~90% clean polarities predicted
    from the moment tensor of ``mechanism``.
    """
    rng = np.random.default_rng(seed)
    strike, dip, rake = mechanism
    sr, dr, rr = np.deg2rad([strike, dip, rake])

    az = rng.uniform(5.0, 355.0, nsta)
    the = rng.uniform(30.0, 90.0, nsta)

    # Moment tensor prediction
    M11 = -np.sin(dr) * np.cos(rr) * np.sin(2 * sr) - np.sin(2 * dr) * np.sin(rr) * np.sin(sr) ** 2
    M22 = np.sin(dr) * np.cos(rr) * np.sin(2 * sr) - np.sin(2 * dr) * np.sin(rr) * np.cos(sr) ** 2
    M33 = np.sin(2 * dr) * np.sin(rr)
    M12 = np.sin(dr) * np.cos(rr) * np.cos(2 * sr) + 0.5 * np.sin(2 * dr) * np.sin(rr) * np.sin(2 * sr)
    M13 = -np.cos(dr) * np.cos(rr) * np.cos(sr) - np.cos(2 * dr) * np.sin(rr) * np.sin(sr)
    M23 = -np.cos(dr) * np.cos(rr) * np.sin(sr) + np.cos(2 * dr) * np.sin(rr) * np.cos(sr)
    M = np.array([[M11, M12, M13], [M12, M22, M23], [M13, M23, M33]])

    pol = np.zeros(nsta, dtype=np.int32)
    for k in range(nsta):
        tr, pr = np.deg2rad(the[k]), np.deg2rad(az[k])
        a = np.array([np.sin(tr) * np.cos(pr), np.sin(tr) * np.sin(pr), -np.cos(tr)])
        pol[k] = -1 if (a @ M @ a) < 0 else 1

    # ~10% bad picks
    nbad = max(1, int(nsta * 0.1))
    bad = rng.choice(nsta, size=nbad, replace=False)
    pol[bad] = -pol[bad]

    qual = (rng.random(nsta) < 0.2).astype(np.int32)
    return az, the, pol, qual


def make_amp_event(nsta=25, seed=7):
    """Synthetic event with S/P amplitude ratios on ~half the stations."""
    az, the, pol, qual = make_synthetic_event(nsta=nsta, seed=seed)
    rng = np.random.default_rng(seed + 1)
    sp = np.zeros(nsta, dtype=np.float64)
    idx = rng.choice(nsta, size=max(1, nsta // 2), replace=False)
    sp[idx] = rng.uniform(-0.4, 0.4, len(idx))
    return az, the, pol, sp


@pytest.fixture(scope="session")
def backend():
    """The native backend (skips the module when libhashhp is missing)."""
    if not backend_mod.get_backend("fortran").available:
        pytest.skip("fortran backend not available")
    return backend_mod.get_backend("fortran")
