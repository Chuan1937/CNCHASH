"""
NCHASH - Python implementation of HASH v1.2 for earthquake focal mechanism inversion.

This package provides a pure Python implementation of the HASH algorithm for
determining earthquake focal mechanisms from polarities, with performance
optimizations using numba and numpy.

Author: Xingchen He
License: BSD 3-Clause
"""

__version__ = "1.0.0"
__author__ = "Xingchen He"

from .amp_subs import (
    focalamp_mc,
    get_misf_amp,
)
from .core import (
    focalmc,
    get_gap,
    get_misfit,
)
from .driver import (
    run_hash,
    run_hash_with_amp,
)
from .io import (
    read_phase_file,
    read_station_file,
    read_velocity_model,
    write_mechanism_output,
)
from .uncertainty import (
    mech_avg,
    mech_prob,
    mech_rot,
)
from .utils import (
    cross_product,
    fp_coord,
    normal_distribution_random,
    strike_dip_rake_to_vectors,
    to_cartesian,
    vectors_to_strike_dip_rake,
)
from .velocity import (
    get_tts,
    make_table,
)

__all__ = [
    "focalmc",
    "focalamp_mc",
    "mech_prob",
    "mech_avg",
    "mech_rot",
    "get_misfit",
    "get_misf_amp",
    "get_gap",
    "make_table",
    "get_tts",
    "cross_product",
    "to_cartesian",
    "fp_coord",
    "normal_distribution_random",
    "strike_dip_rake_to_vectors",
    "vectors_to_strike_dip_rake",
    "read_phase_file",
    "read_station_file",
    "read_velocity_model",
    "write_mechanism_output",
    "run_hash",
    "run_hash_with_amp",
]
