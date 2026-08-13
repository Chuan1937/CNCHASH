"""Velocity-model takeoff tables: Fortran backend vs Numba reference."""

import numpy as np
import pytest

from cnchash import backend as backend_mod

# A simple 3-layer crustal model (depth km, Vp km/s)
SIMPLE_MODEL = np.array(
    [
        [0.0, 4.5],
        [5.0, 5.8],
        [15.0, 6.4],
        [30.0, 7.1],
    ]
)

PARAMS = {
    "del1": 0.0,
    "del2": 60.0,
    "del3": 2.0,
    "dep1": 0.0,
    "dep2": 20.0,
    "dep3": 4.0,
    "pmin": 0.0,
    "nump": 200,
}


@pytest.fixture(scope="module")
def velocity_backends():
    solvers = {}
    for name in backend_mod.available_backends():
        solvers[name] = backend_mod.get_backend(name)
    return solvers


def test_velocity_table_parity(velocity_backends):
    depth = SIMPLE_MODEL[:, 0]
    vel = SIMPLE_MODEL[:, 1]
    tables = {}
    for name, solver in velocity_backends.items():
        tables[name] = solver.build_velocity_table(depth, vel, PARAMS)

    names = list(tables)
    base = tables[names[0]]
    for name in names[1:]:
        t = tables[name]
        assert t["ndel"] == base["ndel"] == 31
        assert t["ndep"] == base["ndep"] == 6
        assert np.allclose(t["delttab"], base["delttab"], atol=1e-12)
        assert np.allclose(t["deptab"], base["deptab"], atol=1e-12)
        # takeoff angles: nonzero where defined, equal where both defined
        valid = (t["table"] != 0.0) & (base["table"] != 0.0)
        assert valid.any()
        diff = np.abs(t["table"][valid] - base["table"][valid])
        assert diff.max() < 1e-6


def test_velocity_tts_parity(velocity_backends):
    depth = SIMPLE_MODEL[:, 0]
    vel = SIMPLE_MODEL[:, 1]
    for solver in velocity_backends.values():
        solver.build_velocity_table(depth, vel, PARAMS)

    probes = [(5.0, 5.0), (30.0, 10.0), (45.0, 16.0), (12.0, 4.0)]
    results = {name: [] for name in velocity_backends}
    for name, solver in velocity_backends.items():
        for del_dist, qdep in probes:
            tt, flag = solver.get_tts(1, del_dist, qdep)
            results[name].append((tt, flag))

    names = list(results)
    base = results[names[0]]
    for name in names[1:]:
        for i, (del_dist, qdep) in enumerate(probes):
            b_tt, b_flag = base[i]
            t_tt, t_flag = results[name][i]
            assert t_flag == b_flag
            if b_flag == 0:
                assert abs(t_tt - b_tt) < 1e-6, f"{del_dist} km, {qdep} km depth"
