"""Native velocity-model takeoff tables.

The table builder is compared against the original HASH MK_TABLE using
the bundled Southern California model (vz.socal) with the original
vel.inc parameters. The original values were extracted with the
compiled HASH_complete MK_TABLE and stored here as a golden reference.
"""

import numpy as np
import pytest

from cnchash import backend as backend_mod

# Golden values from the original HASH_complete MK_TABLE for vz.socal
# with vel.inc parameters (del 0-200 km step 2, dep 0-39 km step 3,
# pmin=0, nump=9000). Format: (del_idx, dep_idx, angle_deg).
GOLDEN = [
    (2, 1, 93.5221481),
    (3, 1, 96.8902469),
    (2, 9, 5.0415606),
    (3, 9, 10.0141449),
    (4, 9, 14.8521449),
    (27, 5, 84.0734024),
    (51, 5, 92.0754998),  # (del=100 km, dep=12 km)
]


def _build_socal_table():
    model = np.loadtxt(
        "/home/yuan/code/CNCHASH/HASH_complete/data/velocity/vz.socal",
        ndmin=2,
    )
    depth, vel = model[:, 0], model[:, 1]
    params = {
        "del1": 0.0,
        "del2": 200.0,
        "del3": 2.0,
        "dep1": 0.0,
        "dep2": 39.0,
        "dep3": 3.0,
        "pmin": 0.0,
        "nump": 9000,
    }
    solver = backend_mod.get_backend("fortran")
    return solver.build_velocity_table(depth, vel, params)


def test_velocity_table_matches_original():
    if not backend_mod.get_backend("fortran").available:
        pytest.skip("fortran backend not available")
    table = _build_socal_table()
    assert table["ndel"] == 101
    assert table["ndep"] == 14
    assert np.count_nonzero(table["table"]) > 1000
    for idel, idep, expected in GOLDEN:
        assert abs(table["table"][idel - 1, idep - 1] - expected) < 0.1


def test_velocity_get_tts_interpolation():
    if not backend_mod.get_backend("fortran").available:
        pytest.skip("fortran backend not available")
    solver = backend_mod.get_backend("fortran")
    _build_socal_table()
    # Inside-range interpolation returns flag 0 with a sane angle.
    tt, iflag = solver.get_tts(1, 10.0, 6.0)
    assert iflag == 0
    assert 0.0 < tt < 180.0
    # Out-of-depth-range returns flag -1.
    tt, iflag = solver.get_tts(1, 10.0, 100.0)
    assert iflag == -1
