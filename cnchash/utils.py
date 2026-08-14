"""
Utility functions for CNCHASH.

Coordinate transformations, vector operations, and the fault-plane
coordinate conversion used by the Python front-end.
"""

import math

import numpy as np

# Constants
PI = 3.14159265358979323846
DEG_TO_RAD = PI / 180.0
RAD_TO_DEG = 180.0 / PI


def cross_product_array(v1, v2):
    """
    Compute cross product for arrays of 3D vectors.

    Parameters
    ----------
    v1 : ndarray, shape (n, 3)
        First set of vectors
    v2 : ndarray, shape (n, 3)
        Second set of vectors

    Returns
    -------
    ndarray, shape (n, 3)
        Cross products v1 x v2
    """
    result = np.empty_like(v1)
    result[:, 0] = v1[:, 1] * v2[:, 2] - v1[:, 2] * v2[:, 1]
    result[:, 1] = v1[:, 2] * v2[:, 0] - v1[:, 0] * v2[:, 2]
    result[:, 2] = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
    return result


def cross_product(v1, v2):
    """
    Compute cross product of two 3D vectors.

    Parameters
    ----------
    v1 : ndarray, shape (3,)
        First vector
    v2 : ndarray, shape (3,)
        Second vector

    Returns
    -------
    ndarray, shape (3,)
        Cross product v1 x v2
    """
    result = np.empty(3, dtype=np.float64)
    result[0] = v1[1] * v2[2] - v1[2] * v2[1]
    result[1] = v1[2] * v2[0] - v1[0] * v2[2]
    result[2] = v1[0] * v2[1] - v1[1] * v2[0]
    return result


def to_cartesian_array(theta, phi, r):
    """
    Transform spherical coordinates to Cartesian for arrays.

    Parameters
    ----------
    theta : ndarray, shape (n,)
        Takeoff angles (degrees, from vertical, up=0, <90 upgoing, >90 downgoing)
    phi : ndarray, shape (n,)
        Azimuths (degrees, East of North)
    r : float
        Radius

    Returns
    -------
    ndarray, shape (n, 3)
        Cartesian coordinates (x, y, z)
        Uses coordinate system: x=north, y=east, z=down
    """
    theta_rad = theta * DEG_TO_RAD
    phi_rad = phi * DEG_TO_RAD
    result = np.empty((theta.shape[0], 3), dtype=np.float64)
    result[:, 2] = -r * np.cos(theta_rad)  # z (down)
    result[:, 0] = r * np.sin(theta_rad) * np.cos(phi_rad)  # x (north)
    result[:, 1] = r * np.sin(theta_rad) * np.sin(phi_rad)  # y (east)
    return result


def to_cartesian(theta, phi, r):
    """
    Transform spherical coordinates to Cartesian.

    Parameters
    ----------
    theta : float
        Takeoff angle (degrees, from vertical, up=0, <90 upgoing, >90 downgoing)
    phi : float
        Azimuth (degrees, East of North)
    r : float
        Radius

    Returns
    -------
    tuple of float
        (x, y, z) Cartesian coordinates
        Uses coordinate system: x=north, y=east, z=down
    """
    theta_rad = theta * DEG_TO_RAD
    phi_rad = phi * DEG_TO_RAD
    z = -r * math.cos(theta_rad)  # z (down)
    x = r * math.sin(theta_rad) * math.cos(phi_rad)  # x (north)
    y = r * math.sin(theta_rad) * math.sin(phi_rad)  # y (east)
    return x, y, z


def fp_coord_vectors_to_angles(faultnorm, slip):
    """
    Convert fault normal and slip vectors to strike, dip, and rake.

    Reference: Aki and Richards, p. 115
    Uses (x,y,z) coordinate system with x=north, y=east, z=down

    Parameters
    ----------
    faultnorm : ndarray, shape (3,)
        Fault normal vector
    slip : ndarray, shape (3,)
        Slip vector

    Returns
    -------
    tuple of float
        (strike, dip, rake) in degrees
    """
    fnorm = np.asarray(faultnorm, dtype=np.float64)
    slip_vec = np.asarray(slip, dtype=np.float64)

    # Check for horizontal fault (undefined strike)
    if 1.0 - abs(fnorm[2]) <= 1e-7:
        # Horizontal fault (matches FPCOOR idir=2)
        dip = 0.0
        phi = math.atan2(-slip_vec[0], slip_vec[1])
        clam = math.cos(phi) * slip_vec[0] + math.sin(phi) * slip_vec[1]
        slam = math.sin(phi) * slip_vec[0] - math.cos(phi) * slip_vec[1]
        lam = math.atan2(slam, clam)
    else:
        # Normal case
        phi = math.atan2(-fnorm[0], fnorm[1])
        a = math.sqrt(fnorm[0] ** 2 + fnorm[1] ** 2)
        del_angle = math.atan2(a, -fnorm[2])
        clam = math.cos(phi) * slip_vec[0] + math.sin(phi) * slip_vec[1]
        slam = -slip_vec[2] / math.sin(del_angle)
        lam = math.atan2(slam, clam)

        if del_angle > 0.5 * PI:
            del_angle = PI - del_angle
            phi = phi + PI
            lam = -lam

    strike = phi * RAD_TO_DEG
    if strike < 0.0:
        strike += 360.0

    dip = del_angle * RAD_TO_DEG
    rake = lam * RAD_TO_DEG

    if rake <= -180.0:
        rake += 360.0
    if rake > 180.0:
        rake -= 360.0

    return strike, dip, rake


def fp_coord_angles_to_vectors(strike, dip, rake):
    """
    Convert strike, dip, and rake to fault normal and slip vectors.

    Reference: Aki and Richards, p. 115
    Uses (x,y,z) coordinate system with x=north, y=east, z=down

    Parameters
    ----------
    strike : float
        Strike angle (degrees)
    dip : float
        Dip angle (degrees)
    rake : float
        Rake/slip angle (degrees)

    Returns
    -------
    tuple of ndarray
        (faultnorm, slip) as 3-element arrays
    """
    phi = strike * DEG_TO_RAD
    del_angle = dip * DEG_TO_RAD
    lam = rake * DEG_TO_RAD

    faultnorm = np.empty(3, dtype=np.float64)
    slip_vec = np.empty(3, dtype=np.float64)

    faultnorm[0] = -math.sin(del_angle) * math.sin(phi)
    faultnorm[1] = math.sin(del_angle) * math.cos(phi)
    faultnorm[2] = -math.cos(del_angle)

    slip_vec[0] = math.cos(lam) * math.cos(phi) + math.cos(del_angle) * math.sin(lam) * math.sin(
        phi
    )
    slip_vec[1] = math.cos(lam) * math.sin(phi) - math.cos(del_angle) * math.sin(lam) * math.cos(
        phi
    )
    slip_vec[2] = -math.sin(lam) * math.sin(del_angle)

    return faultnorm, slip_vec


def fp_coord(strike, dip, rake, faultnorm=None, slip=None, idir=1):
    """
    Convert between fault plane coordinates and strike/dip/rake.

    idir = 1: compute fnorm, slip from strike, dip, rake
    idir = 2: compute strike, dip, rake from fnorm, slip

    Reference: Aki and Richards, p. 115

    Parameters
    ----------
    strike : float or None
        Strike angle (degrees)
    dip : float or None
        Dip angle (degrees)
    rake : float or None
        Rake angle (degrees)
    faultnorm : array-like or None
        Fault normal vector
    slip : array-like or None
        Slip vector
    idir : int
        Direction of conversion (1 or 2)

    Returns
    -------
    tuple or dict
        Depending on idir
    """
    if idir == 1:
        faultnorm, slip_vec = fp_coord_angles_to_vectors(strike, dip, rake)
        return faultnorm, slip_vec
    else:
        strike_out, dip_out, rake_out = fp_coord_vectors_to_angles(faultnorm, slip)
        return strike_out, dip_out, rake_out


def strike_dip_rake_to_vectors(strike, dip, rake):
    """
    Convert strike, dip, rake to fault normal and slip vectors.

    Wrapper for compatibility with higher-level code.
    """
    return fp_coord_angles_to_vectors(strike, dip, rake)


def vectors_to_strike_dip_rake(faultnorm, slip):
    """
    Convert fault normal and slip vectors to strike, dip, rake.

    Wrapper for compatibility with higher-level code.
    """
    return fp_coord_vectors_to_angles(faultnorm, slip)


# Random number generator state
_random_state = {"jran": 314159, "initialized": False}


def _normal_distribution_random_step(jran):
    """One step of the Numerical Recipes normal deviate generator."""
    im = 120050  # overflow at 2**28
    ia = 2311
    ic = 25367
    fran = 0.0
    for _ in range(12):
        jran = (jran * ia + ic) % im
        fran += float(jran) / float(im)
    fran = fran - 6.0
    return fran, jran


def normal_distribution_random():
    """
    Generate normally-distributed random number.

    Uses a saved random state for reproducibility.

    Returns
    -------
    float
        Normally-distributed random number with mean=0, std=1
    """
    global _random_state
    if not _random_state["initialized"]:
        _random_state["jran"] = 314159
        _random_state["initialized"] = True

    fran, _random_state["jran"] = _normal_distribution_random_step(_random_state["jran"])
    return fran


def reset_random_seed(seed=314159):
    """
    Reset the random number generator seed.

    Parameters
    ----------
    seed : int
        Seed value (default: 314159 to match Fortran)
    """
    global _random_state
    _random_state["jran"] = seed
    _random_state["initialized"] = True


def normalize_vector(v):
    """
    Normalize a vector to unit length.

    Parameters
    ----------
    v : ndarray, shape (3,)
        Input vector

    Returns
    -------
    ndarray, shape (3,)
        Normalized vector
    """
    norm = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    if norm > 0:
        return v / norm
    return v


def normalize_vectors_array(v):
    """
    Normalize array of vectors to unit length.

    Parameters
    ----------
    v : ndarray, shape (n, 3)
        Input vectors

    Returns
    -------
    ndarray, shape (n, 3)
        Normalized vectors
    """
    result = np.empty_like(v)
    for i in range(v.shape[0]):
        norm = math.sqrt(v[i, 0] ** 2 + v[i, 1] ** 2 + v[i, 2] ** 2)
        if norm > 0:
            result[i, 0] = v[i, 0] / norm
            result[i, 1] = v[i, 1] / norm
            result[i, 2] = v[i, 2] / norm
        else:
            result[i] = v[i]
    return result


def dot_product(v1, v2):
    """Compute dot product of two vectors."""
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


def vector_length(v):
    """Compute length of a vector."""
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


# Export all utility functions
__all__ = [
    "cross_product",
    "cross_product_array",
    "to_cartesian",
    "to_cartesian_array",
    "fp_coord",
    "fp_coord_vectors_to_angles",
    "fp_coord_angles_to_vectors",
    "strike_dip_rake_to_vectors",
    "vectors_to_strike_dip_rake",
    "normal_distribution_random",
    "reset_random_seed",
    "normalize_vector",
    "normalize_vectors_array",
    "dot_product",
    "vector_length",
]
