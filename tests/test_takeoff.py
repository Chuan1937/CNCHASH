"""Public takeoff-angle API over the native velocity tables."""

import numpy as np
import pytest

from cnchash import (
    TakeoffRangeError,
    TakeoffTable,
    compute_takeoff_angles,
)
from cnchash import backend as backend_mod

UNIFORM_MODEL = (np.array([0.0, 35.0]), np.array([5.5, 5.5]))


def _backend_or_skip():
    try:
        return backend_mod.get_backend("fortran")
    except RuntimeError:
        pytest.skip("fortran backend not available")


@pytest.fixture(scope="module")
def uniform_table():
    _backend_or_skip()
    depth, velocity = UNIFORM_MODEL
    return TakeoffTable(depth, velocity)


def test_uniform_medium_matches_analytic(uniform_table):
    # Homogeneous medium: straight rays, takeoff = atan2(distance, depth)
    # in degrees (HASH convention, 0 = straight up). The native table
    # interpolates linearly in distance between sampled rays, which
    # introduces up to ~0.2 deg of error near-horizontal takeoffs.
    for dist, dep in [(10.0, 10.0), (20.0, 20.0), (10.0, 2.0), (1.0, 10.0)]:
        expected = np.degrees(np.arctan2(dist, dep))
        assert uniform_table.takeoff(dist, dep) == pytest.approx(expected, abs=0.2)


def test_takeoff_depth_out_of_range_raises(uniform_table):
    with pytest.raises(TakeoffRangeError):
        uniform_table.takeoff(10.0, 100.0)


def test_takeoff_distance_beyond_table_extrapolates(uniform_table):
    # Beyond del2 the native GET_TTS extrapolates (iflag=1) and returns
    # a usable angle instead of raising.
    tt = uniform_table.takeoff(150.0, 10.0)
    assert 0.0 < tt < 180.0


def test_takeoff_batch_broadcasting_and_nan(uniform_table):
    out = uniform_table.takeoff_batch([10.0, 20.0], [10.0, 20.0])
    assert out.shape == (2,)
    assert out[0] == pytest.approx(45.0, abs=0.1)
    assert out[1] == pytest.approx(45.0, abs=0.1)

    out = uniform_table.takeoff_batch(np.array([10.0, 20.0]).reshape(2, 1), [10.0, 100.0, 20.0])
    assert out.shape == (2, 3)
    assert np.isnan(out[:, 1]).all()
    assert out[0, 0] == pytest.approx(45.0, abs=0.1)
    assert out[1, 2] == pytest.approx(45.0, abs=0.1)


def test_compute_takeoff_angles_convenience():
    depth, velocity = UNIFORM_MODEL
    scalar = compute_takeoff_angles(depth, velocity, 10.0, 10.0)
    assert isinstance(scalar, float)
    assert scalar == pytest.approx(45.0, abs=0.1)

    arr = compute_takeoff_angles(depth, velocity, [10.0, 20.0], [10.0, 20.0])
    assert arr.shape == (2,)
    assert arr[0] == pytest.approx(45.0, abs=0.1)


def test_table_cache_reuses_handle():
    _backend_or_skip()
    depth, velocity = UNIFORM_MODEL
    t1 = TakeoffTable(depth, velocity)
    t2 = TakeoffTable(depth.copy(), velocity.copy())
    assert t1.ip == t2.ip

    gradient = TakeoffTable(np.array([0.0, 5.0, 35.0]), np.array([5.0, 8.0, 8.0]))
    assert gradient.ip != t1.ip


def test_table_parameter_grid_reflected_in_properties():
    _backend_or_skip()
    depth, velocity = UNIFORM_MODEL
    table = TakeoffTable(depth, velocity)
    assert table.delttab[0] == 0.0
    assert table.delttab[-1] == 120.0
    assert table.deptab[-1] == 35.0


def test_validation_errors():
    _backend_or_skip()
    with pytest.raises(ValueError):
        TakeoffTable(np.array([0.0, 10.0]), np.array([5.5]))  # mismatched
    with pytest.raises(ValueError):
        TakeoffTable(np.array([5.5]), np.array([5.5]))  # one node
    with pytest.raises(ValueError):
        TakeoffTable(np.array([10.0, 0.0]), np.array([5.5, 5.5]))  # decreasing
    with pytest.raises(ValueError):
        TakeoffTable(np.array([0.0, 10.0]), np.array([5.5, 0.0]))  # non-positive
