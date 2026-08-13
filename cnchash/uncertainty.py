"""
Uncertainty analysis for focal mechanisms.

Implements functions for determining the average focal mechanism,
rotation angles between mechanisms, and probability estimates.
Optimized with numba JIT compilation.
"""

import math

import numpy as np
from numba import njit

from .utils import (
    DEG_TO_RAD,
    PI,
    RAD_TO_DEG,
    fp_coord_vectors_to_angles,
)

# Maximum number of mechanisms to process
NMAX0 = 500


@njit(cache=True)
def mech_rotation_angle_numba(norm1, slip1, norm2, slip2):
    """
    Find minimum rotation angle between two mechanisms.

    Tries 4 different combinations of nodal planes and returns
    the minimum rotation angle.

    Parameters
    ----------
    norm1, slip1 : ndarray, shape (3,)
        First mechanism's fault normal and slip vector
    norm2, slip2 : ndarray, shape (3,)
        Second mechanism's fault normal and slip vector

    Returns
    -------
    float
        Minimum rotation angle (degrees)
    """
    # B vectors are cross products (orthogonal to both)
    B1_0 = norm1[1] * slip1[2] - norm1[2] * slip1[1]
    B1_1 = norm1[2] * slip1[0] - norm1[0] * slip1[2]
    B1_2 = norm1[0] * slip1[1] - norm1[1] * slip1[0]

    # Try 4 combinations
    min_rotation = 180.0

    for iteration in range(4):
        # Set up norm2_temp and slip2_temp based on iteration
        if iteration < 2:
            norm2_0, norm2_1, norm2_2 = norm2[0], norm2[1], norm2[2]
            slip2_0, slip2_1, slip2_2 = slip2[0], slip2[1], slip2[2]
        else:
            # Swap norm and slip
            norm2_0, norm2_1, norm2_2 = slip2[0], slip2[1], slip2[2]
            slip2_0, slip2_1, slip2_2 = norm2[0], norm2[1], norm2[2]

        if iteration in (1, 3):
            # Negate both
            norm2_0, norm2_1, norm2_2 = -norm2_0, -norm2_1, -norm2_2
            slip2_0, slip2_1, slip2_2 = -slip2_0, -slip2_1, -slip2_2

        # B2 for this combination
        B2_0 = norm2_1 * slip2_2 - norm2_2 * slip2_1
        B2_1 = norm2_2 * slip2_0 - norm2_0 * slip2_2
        B2_2 = norm2_0 * slip2_1 - norm2_1 * slip2_0

        # Angles between corresponding vectors
        dot_n = norm1[0] * norm2_0 + norm1[1] * norm2_1 + norm1[2] * norm2_2
        dot_s = slip1[0] * slip2_0 + slip1[1] * slip2_1 + slip1[2] * slip2_2
        dot_b = B1_0 * B2_0 + B1_1 * B2_1 + B1_2 * B2_2

        # Clamp to valid range
        dot_n = max(-1.0, min(1.0, dot_n))
        dot_s = max(-1.0, min(1.0, dot_s))
        dot_b = max(-1.0, min(1.0, dot_b))

        phi_n = math.acos(dot_n)
        phi_s = math.acos(dot_s)
        phi_b = math.acos(dot_b)

        # Check for very close mechanisms
        if phi_n < 1e-4 and phi_s < 1e-4 and phi_b < 1e-4:
            return 0.0
        elif phi_n < 1e-4:
            rotation = phi_s * RAD_TO_DEG
        elif phi_s < 1e-4:
            rotation = phi_b * RAD_TO_DEG
        elif phi_b < 1e-4:
            rotation = phi_n * RAD_TO_DEG
        else:
            # Find rotation axis from difference vectors
            n1_0 = norm1[0] - norm2_0
            n1_1 = norm1[1] - norm2_1
            n1_2 = norm1[2] - norm2_2

            n2_0 = slip1[0] - slip2_0
            n2_1 = slip1[1] - slip2_1
            n2_2 = slip1[2] - slip2_2

            n3_0 = B1_0 - B2_0
            n3_1 = B1_1 - B2_1
            n3_2 = B1_2 - B2_2

            # Normalize difference vectors
            l1 = math.sqrt(n1_0**2 + n1_1**2 + n1_2**2)
            l2 = math.sqrt(n2_0**2 + n2_1**2 + n2_2**2)
            l3 = math.sqrt(n3_0**2 + n3_1**2 + n3_2**2)

            if l1 > 0:
                n1_0 /= l1
                n1_1 /= l1
                n1_2 /= l1
            if l2 > 0:
                n2_0 /= l2
                n2_1 /= l2
                n2_2 /= l2
            if l3 > 0:
                n3_0 /= l3
                n3_1 /= l3
                n3_2 /= l3

            # Dot products between difference vectors
            q12 = n1_0 * n2_0 + n1_1 * n2_1 + n1_2 * n2_2
            q13 = n1_0 * n3_0 + n1_1 * n3_1 + n1_2 * n3_2
            q23 = n2_0 * n3_0 + n2_1 * n3_1 + n2_2 * n3_2

            # Find two vectors that aren't parallel
            iout = -1
            for i, qdot in enumerate((q23, q13, q12)):
                if qdot > 0.9999:
                    iout = i

            # Use smallest vector if none are nearly parallel
            if iout == -1:
                # Manual argmin for numba compatibility
                min_len = l1
                iout = 2
                if l2 < min_len:
                    min_len = l2
                    iout = 1
                if l3 < min_len:
                    min_len = l3
                    iout = 0

            # Select two vectors for cross product
            if iout == 0:
                v1_0, v1_1, v1_2 = n1_0, n1_1, n1_2
                v2_0, v2_1, v2_2 = n2_0, n2_1, n2_2
            elif iout == 1:
                v1_0, v1_1, v1_2 = n1_0, n1_1, n1_2
                v2_0, v2_1, v2_2 = n3_0, n3_1, n3_2
            else:
                v1_0, v1_1, v1_2 = n2_0, n2_1, n2_2
                v2_0, v2_1, v2_2 = n3_0, n3_1, n3_2

            # Rotation axis from cross product
            R_0 = v1_1 * v2_2 - v1_2 * v2_1
            R_1 = v1_2 * v2_0 - v1_0 * v2_2
            R_2 = v1_0 * v2_1 - v1_1 * v2_0

            R_len = math.sqrt(R_0**2 + R_1**2 + R_2**2)
            if R_len > 0:
                R_0 /= R_len
                R_1 /= R_len
                R_2 /= R_len

            # Find angle using vector furthest from rotation axis
            t1 = math.acos(max(-1.0, min(1.0, norm1[0] * R_0 + norm1[1] * R_1 + norm1[2] * R_2)))
            t2 = math.acos(max(-1.0, min(1.0, slip1[0] * R_0 + slip1[1] * R_1 + slip1[2] * R_2)))
            t3 = math.acos(max(-1.0, min(1.0, B1_0 * R_0 + B1_1 * R_1 + B1_2 * R_2)))

            # Use angle furthest from 90 degrees - manual argmin for numba
            d1 = abs(t1 - PI / 2.0)
            d2 = abs(t2 - PI / 2.0)
            d3 = abs(t3 - PI / 2.0)
            min_diff = d1
            iuse = 0
            if d2 < min_diff:
                min_diff = d2
                iuse = 1
            if d3 < min_diff:
                min_diff = d3
                iuse = 2

            thetas = [t1, t2, t3]
            phis = [phi_n, phi_s, phi_b]

            cos_val = math.cos(phis[iuse]) - math.cos(thetas[iuse]) * math.cos(thetas[iuse])
            sin_val = math.sin(thetas[iuse]) * math.sin(thetas[iuse])

            if abs(sin_val) > 1e-10:
                cos_rot = max(-1.0, min(1.0, cos_val / sin_val))
                rotation = math.acos(cos_rot) * RAD_TO_DEG
            else:
                rotation = 0.0

        if rotation < min_rotation:
            min_rotation = rotation

    return min_rotation


@njit(cache=True)
def _mech_rot_match_numba(n1, s1, n2, s2):
    """
    Minimum rotation angle between two mechanisms, trying the 4 nodal
    plane combinations. As in the original MECH_ROT, mechanism B (n2, s2)
    is mutated to the best-matching orientation.

    Returns the rotation angle in degrees.
    """
    B1 = np.empty(3)
    B1[0] = n1[1] * s1[2] - n1[2] * s1[1]
    B1[1] = n1[2] * s1[0] - n1[0] * s1[2]
    B1[2] = n1[0] * s1[1] - n1[1] * s1[0]

    rotemp = np.zeros(4)
    phi = np.zeros(3)
    theta = np.zeros(3)

    for iteration in range(4):
        n2t = np.empty(3)
        s2t = np.empty(3)
        if iteration < 2:
            n2t = n2.copy()
            s2t = s2.copy()
        else:
            n2t = s2.copy()
            s2t = n2.copy()
        if iteration == 1 or iteration == 3:
            n2t = -n2t
            s2t = -s2t

        B2 = np.empty(3)
        B2[0] = n2t[1] * s2t[2] - n2t[2] * s2t[1]
        B2[1] = n2t[2] * s2t[0] - n2t[0] * s2t[2]
        B2[2] = n2t[0] * s2t[1] - n2t[1] * s2t[0]

        d = n1[0] * n2t[0] + n1[1] * n2t[1] + n1[2] * n2t[2]
        phi[0] = math.acos(max(-1.0, min(1.0, d)))
        d = s1[0] * s2t[0] + s1[1] * s2t[1] + s1[2] * s2t[2]
        phi[1] = math.acos(max(-1.0, min(1.0, d)))
        d = B1[0] * B2[0] + B1[1] * B2[1] + B1[2] * B2[2]
        phi[2] = math.acos(max(-1.0, min(1.0, d)))

        if phi[0] < 1e-4 and phi[1] < 1e-4 and phi[2] < 1e-4:
            rotemp[iteration] = 0.0
        elif phi[0] < 1e-4:
            rotemp[iteration] = RAD_TO_DEG * phi[1]
        elif phi[1] < 1e-4:
            rotemp[iteration] = RAD_TO_DEG * phi[2]
        elif phi[2] < 1e-4:
            rotemp[iteration] = RAD_TO_DEG * phi[0]
        else:
            dn = np.empty((3, 3))
            dn[:, 0] = n1 - n2t
            dn[:, 1] = s1 - s2t
            dn[:, 2] = B1 - B2
            scale = np.empty(3)
            for j in range(3):
                scale[j] = math.sqrt(dn[0, j] ** 2 + dn[1, j] ** 2 + dn[2, j] ** 2)
                if scale[j] > 0.0:
                    dn[:, j] = dn[:, j] / scale[j]
            qdot = np.empty(3)
            qdot[2] = dn[0, 0] * dn[0, 1] + dn[1, 0] * dn[1, 1] + dn[2, 0] * dn[2, 1]
            qdot[1] = dn[0, 0] * dn[0, 2] + dn[1, 0] * dn[1, 2] + dn[2, 0] * dn[2, 2]
            qdot[0] = dn[0, 1] * dn[0, 2] + dn[1, 1] * dn[1, 2] + dn[2, 1] * dn[2, 2]
            iout = -1
            for i in range(3):
                if qdot[i] > 0.9999:
                    iout = i
            if iout == -1:
                qmins = 10000.0
                for i in range(3):
                    if scale[i] < qmins:
                        qmins = scale[i]
                        iout = i
            n1v = np.empty(3)
            n2v = np.empty(3)
            k = 0
            for j in range(3):
                if j != iout:
                    if k == 0:
                        n1v = dn[:, j].copy()
                        k = 1
                    else:
                        n2v = dn[:, j].copy()
            R = np.empty(3)
            R[0] = n1v[1] * n2v[2] - n1v[2] * n2v[1]
            R[1] = n1v[2] * n2v[0] - n1v[0] * n2v[2]
            R[2] = n1v[0] * n2v[1] - n1v[1] * n2v[0]
            rnorm = math.sqrt(R[0] ** 2 + R[1] ** 2 + R[2] ** 2)
            if rnorm > 0.0:
                R = R / rnorm
            d = n1[0] * R[0] + n1[1] * R[1] + n1[2] * R[2]
            theta[0] = math.acos(max(-1.0, min(1.0, d)))
            d = s1[0] * R[0] + s1[1] * R[1] + s1[2] * R[2]
            theta[1] = math.acos(max(-1.0, min(1.0, d)))
            d = B1[0] * R[0] + B1[1] * R[1] + B1[2] * R[2]
            theta[2] = math.acos(max(-1.0, min(1.0, d)))
            qmindif = 1000.0
            iuse = 0
            for i in range(3):
                dval = abs(theta[i] - PI / 2.0)
                if dval < qmindif:
                    qmindif = dval
                    iuse = i
            cval = math.cos(phi[iuse]) - math.cos(theta[iuse]) * math.cos(theta[iuse])
            sval = math.sin(theta[iuse]) * math.sin(theta[iuse])
            if abs(sval) > 1e-10:
                rotemp[iteration] = RAD_TO_DEG * math.acos(max(-1.0, min(1.0, cval / sval)))
            else:
                rotemp[iteration] = RAD_TO_DEG * phi[iuse]

    rota = 180.0
    irot = 0
    for iteration in range(4):
        if abs(rotemp[iteration]) < rota:
            rota = abs(rotemp[iteration])
            irot = iteration
    if irot >= 2:
        # Mutate mechanism B in place (numba: element-wise writes only).
        tmp0, tmp1, tmp2 = n2[0], n2[1], n2[2]
        n2[0] = s2[0]
        n2[1] = s2[1]
        n2[2] = s2[2]
        s2[0] = tmp0
        s2[1] = tmp1
        s2[2] = tmp2
    if irot == 1 or irot == 3:
        n2[0] = -n2[0]
        n2[1] = -n2[1]
        n2[2] = -n2[2]
        s2[0] = -s2[0]
        s2[1] = -s2[1]
        s2[2] = -s2[2]
    return rota


@njit(cache=True)
def _mech_avg_faithful_numba(nf, norm1, norm2):
    """
    Average focal mechanism (faithful MECH_AVG): plane matching via
    minimum rotation, weighted orthogonalization with maxmisf = 0.01 deg.
    """
    if nf <= 1:
        return norm1[:, 0].copy(), norm2[:, 0].copy()

    norm1_avg = norm1[:, 0].copy()
    norm2_avg = norm2[:, 0].copy()
    ref1 = norm1[:, 0].copy()
    ref2 = norm2[:, 0].copy()

    for i in range(1, nf):
        t1 = norm1[:, i].copy()
        t2 = norm2[:, i].copy()
        _rota = _mech_rot_match_numba(ref1, ref2, t1, t2)
        norm1_avg = norm1_avg + t1
        norm2_avg = norm2_avg + t2

    l1 = math.sqrt(norm1_avg[0] ** 2 + norm1_avg[1] ** 2 + norm1_avg[2] ** 2)
    l2 = math.sqrt(norm2_avg[0] ** 2 + norm2_avg[1] ** 2 + norm2_avg[2] ** 2)
    if l1 > 0.0:
        norm1_avg = norm1_avg / l1
    if l2 > 0.0:
        norm2_avg = norm2_avg / l2

    avang1 = 0.0
    avang2 = 0.0
    for i in range(nf):
        d11 = norm1[0, i] * norm1_avg[0] + norm1[1, i] * norm1_avg[1] + norm1[2, i] * norm1_avg[2]
        d22 = norm2[0, i] * norm2_avg[0] + norm2[1, i] * norm2_avg[1] + norm2[2, i] * norm2_avg[2]
        a11 = math.acos(max(-1.0, min(1.0, d11)))
        a22 = math.acos(max(-1.0, min(1.0, d22)))
        avang1 += a11 * a11
        avang2 += a22 * a22
    avang1 = math.sqrt(avang1 / nf)
    avang2 = math.sqrt(avang2 / nf)

    if avang1 + avang2 < 0.0001:
        return norm1_avg, norm2_avg

    fract1 = avang1 / (avang1 + avang2)
    for _ in range(100):
        dot1 = norm1_avg[0] * norm2_avg[0] + norm1_avg[1] * norm2_avg[1] + norm1_avg[2] * norm2_avg[2]
        misf = 90.0 - math.acos(max(-1.0, min(1.0, dot1))) * RAD_TO_DEG
        if abs(misf) <= 0.01:
            break
        theta1 = misf * fract1 * DEG_TO_RAD
        theta2 = misf * (1.0 - fract1) * DEG_TO_RAD
        for j in range(3):
            temp = norm1_avg[j]
            norm1_avg[j] = norm1_avg[j] - norm2_avg[j] * math.sin(theta1)
            norm2_avg[j] = norm2_avg[j] - temp * math.sin(theta2)
        l1 = math.sqrt(norm1_avg[0] ** 2 + norm1_avg[1] ** 2 + norm1_avg[2] ** 2)
        l2 = math.sqrt(norm2_avg[0] ** 2 + norm2_avg[1] ** 2 + norm2_avg[2] ** 2)
        if l1 > 0.0:
            norm1_avg = norm1_avg / l1
        if l2 > 0.0:
            norm2_avg = norm2_avg / l2

    return norm1_avg, norm2_avg


def mech_rot(norm1, slip1, norm2, slip2):
    """
    Find the minimum rotation angle between two mechanisms.

    Parameters
    ----------
    norm1, slip1 : array-like, shape (3,)
        First mechanism's fault normal and slip vector
    norm2, slip2 : array-like, shape (3,)
        Second mechanism's fault normal and slip vector

    Returns
    -------
    float
        Minimum rotation angle (degrees)
    """
    return mech_rotation_angle_numba(
        np.asarray(norm1, dtype=np.float64),
        np.asarray(slip1, dtype=np.float64),
        np.asarray(norm2, dtype=np.float64),
        np.asarray(slip2, dtype=np.float64),
    )


@njit(cache=True)
def mech_average_numba(nf, norm1, norm2):
    """
    Determine the average focal mechanism of a set of mechanisms.

    Parameters
    ----------
    nf : int
        Number of fault planes
    norm1 : ndarray, shape (3, nf)
        Normal vectors to fault planes
    norm2 : ndarray, shape (3, nf)
        Slip vectors

    Returns
    -------
    tuple
        (norm1_avg, norm2_avg) - average normal and slip vectors
    """
    if nf <= 1:
        norm1_avg = norm1[:, 0].copy()
        norm2_avg = norm2[:, 0].copy()
        return norm1_avg, norm2_avg

    # Initialize with first mechanism
    norm1_avg = np.zeros(3, dtype=np.float64)
    norm2_avg = np.zeros(3, dtype=np.float64)

    ref1 = norm1[:, 0].copy()
    ref2 = norm2[:, 0].copy()

    # Accumulate vectors (after matching)
    for i in range(nf):
        # Get current mechanism vectors
        n1 = norm1[:, i].copy()
        n2 = norm2[:, i].copy()

        # Match to reference (try 4 combinations)
        min_rot = 180.0
        best_n1 = n1.copy()
        best_n2 = n2.copy()

        for iter in range(4):
            if iter < 2:
                t1 = n1.copy()
                t2 = n2.copy()
            else:
                t1 = n2.copy()
                t2 = n1.copy()

            if iter in (1, 3):
                t1 = -t1
                t2 = -t2

            rot = mech_rotation_angle_numba(ref1, ref2, t1, t2)
            if rot < min_rot:
                min_rot = rot
                best_n1 = t1.copy()
                best_n2 = t2.copy()

        # Add to average
        norm1_avg += best_n1
        norm2_avg += best_n2

    # Normalize
    l1 = math.sqrt(norm1_avg[0] ** 2 + norm1_avg[1] ** 2 + norm1_avg[2] ** 2)
    l2 = math.sqrt(norm2_avg[0] ** 2 + norm2_avg[1] ** 2 + norm2_avg[2] ** 2)

    if l1 > 0:
        norm1_avg /= l1
    if l2 > 0:
        norm2_avg /= l2

    # Make orthogonal (iterative adjustment)
    for _ in range(100):
        dot = (
            norm1_avg[0] * norm2_avg[0] + norm1_avg[1] * norm2_avg[1] + norm1_avg[2] * norm2_avg[2]
        )
        misf = 90.0 - math.acos(max(-1.0, min(1.0, dot))) * RAD_TO_DEG

        if abs(misf) <= 0.01:
            break

        # Adjust both vectors
        theta1 = misf * 0.5 * DEG_TO_RAD
        theta2 = misf * 0.5 * DEG_TO_RAD

        # Store old values
        n1_0, n1_1, n1_2 = norm1_avg[0], norm1_avg[1], norm1_avg[2]
        n2_0, n2_1, n2_2 = norm2_avg[0], norm2_avg[1], norm2_avg[2]

        # Adjust
        norm1_avg[0] = n1_0 - n2_0 * math.sin(theta1)
        norm1_avg[1] = n1_1 - n2_1 * math.sin(theta1)
        norm1_avg[2] = n1_2 - n2_2 * math.sin(theta1)

        norm2_avg[0] = n2_0 - n1_0 * math.sin(theta2)
        norm2_avg[1] = n2_1 - n1_1 * math.sin(theta2)
        norm2_avg[2] = n2_2 - n1_2 * math.sin(theta2)

        # Renormalize
        l1 = math.sqrt(norm1_avg[0] ** 2 + norm1_avg[1] ** 2 + norm1_avg[2] ** 2)
        l2 = math.sqrt(norm2_avg[0] ** 2 + norm2_avg[1] ** 2 + norm2_avg[2] ** 2)

        if l1 > 0:
            norm1_avg /= l1
        if l2 > 0:
            norm2_avg /= l2

    return norm1_avg, norm2_avg


def mech_avg(nf, norm1, norm2):
    """
    Determine the average focal mechanism of a set of mechanisms.

    Parameters
    ----------
    nf : int
        Number of fault planes
    norm1 : ndarray, shape (3, nf) or (nf, 3)
        Normal vectors to fault planes
    norm2 : ndarray, shape (3, nf) or (nf, 3)
        Slip vectors

    Returns
    -------
    tuple
        (norm1_avg, norm2_avg) - average normal and slip vectors
    """
    # Ensure correct shape (3, nf)
    if norm1.shape[0] != 3:
        norm1 = norm1.T
        norm2 = norm2.T

    return _mech_avg_faithful_numba(nf, norm1.astype(np.float64), norm2.astype(np.float64))


@njit(cache=True)
def mech_probability_numba(nf, norm1in, norm2in, cangle, prob_max):
    """
    Determine average focal mechanism and check for multiple solutions.

    Faithful port of the original MECH_PROB: outlier removal uses the
    full minimum-rotation search, RMS differences are computed after
    matching planes, and the average is computed with the weighted
    orthogonalization from MECH_AVG.

    Parameters
    ----------
    nf : int
        Number of fault planes
    norm1in : ndarray, shape (3, nf)
        Normal vectors to fault planes
    norm2in : ndarray, shape (3, nf)
        Slip vectors
    cangle : float
        Cutoff angle (degrees)
    prob_max : float
        Cutoff percent for multiples

    Returns
    -------
    tuple
        (nsltn, str_avg, dip_avg, rak_avg, prob, rms_diff)
    """
    if nf <= 1:
        s, d, r = fp_coord_vectors_to_angles(norm1in[:, 0], norm2in[:, 0])

        str_avg = np.zeros(5, dtype=np.float64)
        dip_avg = np.zeros(5, dtype=np.float64)
        rak_avg = np.zeros(5, dtype=np.float64)
        prob = np.zeros(5, dtype=np.float64)
        rms_diff = np.zeros((2, 5), dtype=np.float64)

        str_avg[0] = s
        dip_avg[0] = d
        rak_avg[0] = r
        prob[0] = 1.0
        return 1, str_avg, dip_avg, rak_avg, prob, rms_diff

    # Output arrays
    str_avg = np.zeros(5, dtype=np.float64)
    dip_avg = np.zeros(5, dtype=np.float64)
    rak_avg = np.zeros(5, dtype=np.float64)
    prob = np.zeros(5, dtype=np.float64)
    rms_diff = np.zeros((2, 5), dtype=np.float64)

    norm1 = norm1in.copy()
    norm2 = norm2in.copy()
    nfault = nf
    nc = nf
    nsltn = 5

    for imult in range(5):
        if nc < 1:
            nsltn = imult
            break

        # Repeatedly remove the mechanism with the largest angular
        # difference from the running average, until all are within cangle.
        for _icount in range(nf):
            norm1_avg, norm2_avg = _mech_avg_faithful_numba(nc, norm1, norm2)
            maxrot = 0.0
            imax = 0
            for i in range(nc):
                t1 = norm1[:, i].copy()
                t2 = norm2[:, i].copy()
                rot_angle = _mech_rot_match_numba(norm1_avg, norm2_avg, t1, t2)
                if abs(rot_angle) > maxrot:
                    maxrot = abs(rot_angle)
                    imax = i
            if maxrot <= cangle:
                break
            # Move the outlier to the end.
            t1 = norm1[:, imax].copy()
            t2 = norm2[:, imax].copy()
            for j in range(imax, nc - 1):
                norm1[:, j] = norm1[:, j + 1]
                norm2[:, j] = norm2[:, j + 1]
            norm1[:, nc - 1] = t1
            norm2[:, nc - 1] = t2
            nc = nc - 1

        prob[imult] = float(nc) / float(nfault)

        if imult > 0 and prob[imult] < prob_max:
            nsltn = imult
            break

        # Set up for the next round: move outliers to the front.
        for j in range(nfault - nc):
            norm1[:, j] = norm1[:, j + nc]
            norm2[:, j] = norm2[:, j + nc]
        nc = nfault - nc

        # RMS angular differences against the average, after matching planes.
        rms1 = 0.0
        rms2 = 0.0
        for i in range(nfault):
            t1 = norm1in[:, i].copy()
            t2 = norm2in[:, i].copy()
            _rot_angle = _mech_rot_match_numba(norm1_avg, norm2_avg, t1, t2)
            d11 = t1[0] * norm1_avg[0] + t1[1] * norm1_avg[1] + t1[2] * norm1_avg[2]
            d22 = t2[0] * norm2_avg[0] + t2[1] * norm2_avg[1] + t2[2] * norm2_avg[2]
            d11 = max(-1.0, min(1.0, d11))
            d22 = max(-1.0, min(1.0, d22))
            a11 = math.acos(d11)
            a22 = math.acos(d22)
            rms1 += a11 * a11
            rms2 += a22 * a22

        rms_diff[0, imult] = RAD_TO_DEG * math.sqrt(rms1 / nfault)
        rms_diff[1, imult] = RAD_TO_DEG * math.sqrt(rms2 / nfault)

        # Convert to strike, dip, rake.
        s, d, r = fp_coord_vectors_to_angles(norm1_avg, norm2_avg)
        str_avg[imult] = s
        dip_avg[imult] = d
        rak_avg[imult] = r

    return nsltn, str_avg, dip_avg, rak_avg, prob, rms_diff


def mech_prob(nf, norm1, norm2, cangle, prob_max):
    """
    Determine average focal mechanism and check for multiple solutions.

    Parameters
    ----------
    nf : int
        Number of fault planes
    norm1 : ndarray, shape (3, nf) or (nf, 3)
        Normal vectors to fault planes
    norm2 : ndarray, shape (3, nf) or (nf, 3)
        Slip vectors
    cangle : float
        Cutoff angle (degrees)
    prob_max : float
        Cutoff percent for multiples (0-1)

    Returns
    -------
    dict
        Dictionary containing:
        - 'nsltn': number of solutions (up to 5)
        - 'strike_avg': array of strike angles
        - 'dip_avg': array of dip angles
        - 'rake_avg': array of rake angles
        - 'prob': probability for each solution
        - 'rms_diff': RMS difference for each solution
    """
    # Ensure correct shape (3, nf)
    if norm1.shape[0] != 3:
        norm1 = norm1.T
        norm2 = norm2.T

    nsltn, str_avg, dip_avg, rak_avg, prob, rms_diff = mech_probability_numba(
        nf, norm1.astype(np.float64), norm2.astype(np.float64), cangle, prob_max
    )

    # Trim to actual number of solutions
    return {
        "nsltn": nsltn,
        "strike_avg": str_avg[:nsltn],
        "dip_avg": dip_avg[:nsltn],
        "rake_avg": rak_avg[:nsltn],
        "prob": prob[:nsltn],
        "rms_diff": rms_diff[:, :nsltn],
    }


# Export all functions
__all__ = [
    "mech_rot",
    "mech_avg",
    "mech_prob",
]
