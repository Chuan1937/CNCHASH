"""Thread-count determinism for the native backend."""

import numpy as np
import pytest
from conftest import make_synthetic_event

from cnchash import backend as backend_mod
from cnchash import run_hash


@pytest.mark.skipif(
    not backend_mod.get_backend("fortran").available, reason="fortran backend unavailable"
)
@pytest.mark.parametrize("threads", [1, 2, 4])
def test_fortran_thread_determinism(threads):
    az, the, pol, qual = make_synthetic_event(nsta=22, seed=3)
    azi_mc = np.repeat(az.reshape(-1, 1), 30, axis=1)
    the_mc = np.repeat(the.reshape(-1, 1), 30, axis=1)

    res = run_hash(azi_mc, the_mc, pol, qual, nmc=30, backend="fortran", num_threads=threads)
    ref = run_hash(azi_mc, the_mc, pol, qual, nmc=30, backend="fortran", num_threads=1)

    assert res["nf"] == ref["nf"]
    assert res["nmult"] == ref["nmult"]
    assert res["quality"] == ref["quality"]
    assert np.allclose(res["strike"], ref["strike"], atol=1e-9)
    assert np.allclose(res["faults"], ref["faults"], atol=1e-12)
    assert abs(res["mfrac"] - ref["mfrac"]) < 1e-12
    assert abs(res["stdr"] - ref["stdr"]) < 1e-12
