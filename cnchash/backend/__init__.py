"""Backend abstraction for CNCHASH.

A single backend implements the computation contract:

- ``fortran``: Modern Fortran/OpenMP native backend (libhashhp)

The public API uses it via ``backend="auto"`` (the default) or
``backend="fortran"``. Without the compiled library the package raises
a clear error explaining how to build the native backend.
"""

from .base import HashBackend
from .fortran_backend import FortranBackend

_BACKENDS = {}
_FORCE_BACKEND = None

_BACKENDS["fortran"] = FortranBackend()


def get_backend(name=None):
    """Get a backend instance by name (``"auto"`` resolves to fortran)."""
    if name is None or name == "auto":
        name = _FORCE_BACKEND or "fortran"
    if name not in _BACKENDS:
        raise ValueError(f"Unknown backend {name!r}. Available: {available_backends()}")
    backend = _BACKENDS[name]
    if not backend.available:
        raise RuntimeError(
            "The fortran backend is not available: libhashhp was not found. "
            "Build it with `cmake -S . -B build && cmake --build build` "
            "(or set CNCHASH_HASHHP_LIB to its path)."
        )
    return backend


def available_backends():
    """Names of the backends that can be used right now."""
    return [name for name, backend in _BACKENDS.items() if backend.available]


def get_backend_info():
    """Detailed information about the active (default) backend."""
    backend = get_backend("auto")
    return backend.info()


def set_num_threads(num_threads):
    """Set the thread count used by the active backend."""
    return get_backend("auto").set_num_threads(num_threads)


def get_num_threads():
    """Current thread count of the active backend."""
    return get_backend("auto").get_num_threads()


def _set_forced_backend(name):
    """Force a backend for the whole process (test/CI support)."""
    global _FORCE_BACKEND
    if name not in (None, "auto", "fortran"):
        raise ValueError(f"Unknown backend {name!r}")
    _FORCE_BACKEND = name


__all__ = [
    "HashBackend",
    "FortranBackend",
    "get_backend",
    "available_backends",
    "get_backend_info",
    "set_num_threads",
    "get_num_threads",
]
