"""
Takeoff-angle computation from 1D velocity models.

This module is the Python counterpart of the original HASH MK_TABLE +
GET_TTS workflow used by hash_driver2.f: a takeoff-angle table is built
from a 1D P-wave velocity model by the native Fortran/OpenMP core, and
lookups interpolate the takeoff angle for any (distance, source depth)
pair.

Takeoff angle convention (identical to HASH): degrees from the
vertical, 0 = straight up, 90 = horizontal, 180 = straight down.
"""

import numpy as np

from . import backend as backend_mod

# Default table-generation parameters (same as the native backend and
# the original HASH vel.inc defaults): distance 0-120 km step 1 km,
# depth 0-35 km step 1 km, ray parameter 0-pmax in 1000 steps.
DEFAULT_TABLE_PARAMS = {
    "del1": 0.0,
    "del2": 120.0,
    "del3": 1.0,
    "dep1": 0.0,
    "dep2": 35.0,
    "dep3": 1.0,
    "pmin": 0.0,
    "nump": 1000,
}

# Cache of built tables keyed by model content so that the same
# velocity model reuses the backend table handle (ip) instead of
# rebuilding it. The backend keeps every built table alive, so the
# handles stay valid for the lifetime of the process.
_TABLE_CACHE = {}


def _cache_key(depth, velocity, params, backend_name):
    return (
        backend_name,
        depth.tobytes(),
        velocity.tobytes(),
        tuple(sorted(params.items())),
    )


class TakeoffRangeError(ValueError):
    """Takeoff lookup outside the table's usable depth/distance range."""


class TakeoffTable:
    """Takeoff-angle table built from a 1D velocity model.

    Parameters
    ----------
    depth : array_like
        Depth nodes of the velocity model (km), increasing.
    velocity : array_like
        P-wave velocity at each depth node (km/s).
    params : dict, optional
        Table-generation parameters. Keys: del1/del2/del3 (distance
        range and step, km), dep1/dep2/dep3 (depth range and step, km),
        pmin (minimum ray parameter, s/km), nump (number of ray
        parameters). Defaults to DEFAULT_TABLE_PARAMS.
    backend : str
        Backend name ("auto" or "fortran").
    """

    def __init__(self, depth, velocity, params=None, backend="auto"):
        # Raises RuntimeError with a clear build hint when libhashhp is
        # missing.
        self._backend = backend_mod.get_backend(backend)
        self.depth = np.asarray(depth, dtype=np.float64)
        self.velocity = np.asarray(velocity, dtype=np.float64)
        if self.depth.ndim != 1 or self.velocity.shape != self.depth.shape:
            raise ValueError("depth and velocity must be 1-D arrays of the same length")
        if self.depth.size < 2:
            raise ValueError("velocity model requires at least two depth nodes")
        if np.any(np.diff(self.depth) < 0):
            raise ValueError("velocity model depths must be non-decreasing")
        if np.any(self.velocity <= 0):
            raise ValueError("velocities must be positive")
        self.params = dict(DEFAULT_TABLE_PARAMS if params is None else params)

        key = _cache_key(self.depth, self.velocity, self.params, self._backend.name)
        cached = _TABLE_CACHE.get(key)
        if cached is not None:
            ip, table_data = cached
        else:
            table_data = self._backend.build_velocity_table(self.depth, self.velocity, self.params)
            ip = table_data["ip"]
            _TABLE_CACHE[key] = (ip, table_data)
        self._ip = ip
        self._table_data = table_data

    @property
    def delttab(self):
        """Distance nodes of the table (km)."""
        return self._table_data["delttab"]

    @property
    def deptab(self):
        """Source depth nodes of the table (km)."""
        return self._table_data["deptab"]

    @property
    def ip(self):
        """Backend table handle (1-based index into the backend's tables)."""
        return self._ip

    def takeoff(self, distance, source_depth):
        """Takeoff angle (degrees) for a single distance/depth pair.

        Raises TakeoffRangeError when the lookup falls outside the
        table's usable range (source depth outside deptab, or no ray
        reaching the distance). An extrapolated value (beyond the
        outermost distance coverage) is returned as-is.
        """
        tt, iflag = self._backend.get_tts(self._ip, float(distance), float(source_depth))
        if iflag < 0:
            raise TakeoffRangeError(
                f"takeoff lookup failed for distance={distance} km, "
                f"depth={source_depth} km (table depth range "
                f"{self.deptab[0]:g}-{self.deptab[-1]:g} km); extend the "
                f"table range via params (dep1/dep2/del1/del2)"
            )
        return float(tt)

    def takeoff_batch(self, distances, source_depths):
        """Vectorized takeoff angles; inputs broadcast to a common shape.

        Lookups outside the usable range yield NaN instead of raising.
        """
        dist = np.asarray(distances, dtype=np.float64)
        dep = np.asarray(source_depths, dtype=np.float64)
        shape = np.broadcast_shapes(dist.shape, dep.shape)
        bdist = np.broadcast_to(dist, shape)
        bdep = np.broadcast_to(dep, shape)
        out = np.empty(shape, dtype=np.float64)
        for idx in np.ndindex(shape):
            tt, iflag = self._backend.get_tts(self._ip, float(bdist[idx]), float(bdep[idx]))
            out[idx] = tt if iflag >= 0 else np.nan
        return out

    def __repr__(self):
        return (
            f"TakeoffTable(ip={self._ip}, "
            f"model_points={self.depth.size}, "
            f"del={self.delttab[0]:g}-{self.delttab[-1]:g} km, "
            f"dep={self.deptab[0]:g}-{self.deptab[-1]:g} km)"
        )


def compute_takeoff_angles(depth, velocity, distances, source_depths, params=None, backend="auto"):
    """Build a table from a 1D velocity model and look up takeoff angles.

    Parameters
    ----------
    depth, velocity : array_like
        1D velocity model (km, km/s), see TakeoffTable.
    distances, source_depths : array_like
        Source-station distances (km) and source depths (km); scalar
        inputs are broadcast against each other.
    params : dict, optional
        Table-generation parameters, see TakeoffTable.
    backend : str
        Backend name.

    Returns
    -------
    ndarray or float
        Takeoff angles in degrees. Scalar for scalar scalar inputs,
        otherwise an array with the broadcast shape. NaN where the
        lookup falls outside the table range.
    """
    table = TakeoffTable(depth, velocity, params=params, backend=backend)
    dist = np.asarray(distances, dtype=np.float64)
    dep = np.asarray(source_depths, dtype=np.float64)
    if dist.ndim == 0 and dep.ndim == 0:
        return table.takeoff(float(dist), float(dep))
    return table.takeoff_batch(dist, dep)
