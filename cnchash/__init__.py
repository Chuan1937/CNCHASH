"""
CNCHASH - high-performance HASH for earthquake focal mechanism inversion.

CNCHASH provides the HASH algorithm (Hardebeck and Shearer, 2002, 2003)
for determining earthquake focal mechanisms from P-wave polarities and
optional S/P amplitude ratios, with a Modern Fortran/OpenMP native
backend (libhashhp) and a pure-Python front-end.

Author: He XingChen
License: BSD 3-Clause
"""

__version__ = "2.0.0"
__author__ = "He XingChen"

from .backend import (
    available_backends,
    get_backend_info,
    get_num_threads,
    set_num_threads,
)
from .driver import (
    run_hash,
    run_hash_batch,
    run_hash_batch_with_amp,
    run_hash_event_payload,
    run_hash_from_file,
    run_hash_with_amp,
)
from .io import (
    read_phase_file,
    read_station_file,
    read_velocity_model,
    write_mechanism_output,
)
from .utils import (
    cross_product,
    fp_coord,
    normal_distribution_random,
    strike_dip_rake_to_vectors,
    to_cartesian,
    vectors_to_strike_dip_rake,
)

__all__ = [
    # Main functions
    "run_hash",
    "run_hash_event_payload",
    "run_hash_with_amp",
    "run_hash_batch",
    "run_hash_batch_with_amp",
    "run_hash_from_file",
    # Backend control
    "available_backends",
    "get_backend_info",
    "set_num_threads",
    "get_num_threads",
    # Utilities
    "cross_product",
    "to_cartesian",
    "fp_coord",
    "normal_distribution_random",
    "strike_dip_rake_to_vectors",
    "vectors_to_strike_dip_rake",
    # I/O
    "read_phase_file",
    "read_station_file",
    "read_velocity_model",
    "write_mechanism_output",
]
