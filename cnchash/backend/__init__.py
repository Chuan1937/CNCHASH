"""Backend abstraction for CNCHASH.

Two backends implement the same computation contract:

- ``numba``: pure Python/Numba reference implementation (always available)
- ``fortran``: Modern Fortran/OpenMP native backend (optional; used when
  the compiled library is available)

The public API picks a backend with ``backend="auto" | "numba" | "fortran"``
(see :func:`cnchash.run_hash`).
"""

from .base import HashBackend
from .fortran_backend import FortranBackend
from .numba_backend import NumbaBackend

_BACKENDS = {}
_DEFAULT_BACKEND = None
_FORCE_BACKEND = None


def _register_backend(name, backend):
    """Register a backend instance under its name."""
    _BACKENDS[name] = backend


_register_backend("numba", NumbaBackend())
_fortran_backend = FortranBackend()
if _fortran_backend.available:
    _register_backend("fortran", _fortran_backend)


def get_backend(name=None):
    """Get a backend instance by name (``"auto"`` resolves to the best
    available: fortran if present, otherwise numba)."""
    if name is None or name == "auto":
        if _FORCE_BACKEND is not None:
            name = _FORCE_BACKEND
        elif "fortran" in _BACKENDS:
            name = "fortran"
        else:
            name = "numba"
    if name not in _BACKENDS:
        raise ValueError(
            f"Unknown backend {name!r}. Available: {available_backends()}"
        )
    return _BACKENDS[name]


def available_backends():
    """Names of the backends that can be used right now."""
    return list(_BACKENDS.keys())


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
    if name not in (None, "auto", "numba", "fortran"):
        raise ValueError(f"Unknown backend {name!r}")
    _FORCE_BACKEND = name


__all__ = [
    "HashBackend",
    "NumbaBackend",
    "FortranBackend",
    "get_backend",
    "available_backends",
    "get_backend_info",
    "set_num_threads",
    "get_num_threads",
]
