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

    return mech_average_numba(nf, norm1.astype(np.float64), norm2.astype(np.float64))


@njit(cache=True)
def mech_probability_numba(nf, norm1in, norm2in, cangle, prob_max):
    """
    Determine average focal mechanism and check for multiple solutions.
    Optimized version with reduced memory allocations.

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

    # Use indices instead of copying arrays
    # active[i] = 1 if mechanism i is still active
    active = np.ones(nf, dtype=np.int32)
    nsltn = 0

    for imult in range(5):
        # Count active mechanisms
        nc = 0
        for i in range(nf):
            if active[i] == 1:
                nc += 1

        if nc < 1:
            nsltn = imult
            break

        # Collect active mechanism indices
        active_idx = np.zeros(nc, dtype=np.int32)
        idx = 0
        for i in range(nf):
            if active[i] == 1:
                active_idx[idx] = i
                idx += 1

        # Iteratively find average and remove outliers
        for _icount in range(nf):
            # Compute average directly without extra function call
            norm1_avg = np.zeros(3, dtype=np.float64)
            norm2_avg = np.zeros(3, dtype=np.float64)
            ref1 = np.zeros(3, dtype=np.float64)
            ref2 = np.zeros(3, dtype=np.float64)

            # Get first active as reference
            ref1[0] = norm1in[0, active_idx[0]]
            ref1[1] = norm1in[1, active_idx[0]]
            ref1[2] = norm1in[2, active_idx[0]]
            ref2[0] = norm2in[0, active_idx[0]]
            ref2[1] = norm2in[1, active_idx[0]]
            ref2[2] = norm2in[2, active_idx[0]]

            # Accumulate matched vectors
            for k in range(nc):
                idx_k = active_idx[k]
                n1_0 = norm1in[0, idx_k]
                n1_1 = norm1in[1, idx_k]
                n1_2 = norm1in[2, idx_k]
                n2_0 = norm2in[0, idx_k]
                n2_1 = norm2in[1, idx_k]
                n2_2 = norm2in[2, idx_k]

                # Find best matching orientation (4 combinations)
                best_rot = 180.0
                best_n1_0, best_n1_1, best_n1_2 = n1_0, n1_1, n1_2
                best_n2_0, best_n2_1, best_n2_2 = n2_0, n2_1, n2_2

                for iteration in range(4):
                    if iteration < 2:
                        t1_0, t1_1, t1_2 = n1_0, n1_1, n1_2
                        t2_0, t2_1, t2_2 = n2_0, n2_1, n2_2
                    else:
                        t1_0, t1_1, t1_2 = n2_0, n2_1, n2_2
                        t2_0, t2_1, t2_2 = n1_0, n1_1, n1_2

                    if iteration in (1, 3):
                        t1_0, t1_1, t1_2 = -t1_0, -t1_1, -t1_2
                        t2_0, t2_1, t2_2 = -t2_0, -t2_1, -t2_2

                    # Quick angle estimate using dot products
                    dot_n = ref1[0] * t1_0 + ref1[1] * t1_1 + ref1[2] * t1_2
                    dot_s = ref2[0] * t2_0 + ref2[1] * t2_1 + ref2[2] * t2_2
                    dot_n = max(-1.0, min(1.0, dot_n))
                    dot_s = max(-1.0, min(1.0, dot_s))
                    rot_est = 0.5 * (math.acos(dot_n) + math.acos(dot_s)) * RAD_TO_DEG

                    if rot_est < best_rot:
                        best_rot = rot_est
                        best_n1_0, best_n1_1, best_n1_2 = t1_0, t1_1, t1_2
                        best_n2_0, best_n2_1, best_n2_2 = t2_0, t2_1, t2_2

                norm1_avg[0] += best_n1_0
                norm1_avg[1] += best_n1_1
                norm1_avg[2] += best_n1_2
                norm2_avg[0] += best_n2_0
                norm2_avg[1] += best_n2_1
                norm2_avg[2] += best_n2_2

            # Normalize
            l1 = math.sqrt(norm1_avg[0] ** 2 + norm1_avg[1] ** 2 + norm1_avg[2] ** 2)
            l2 = math.sqrt(norm2_avg[0] ** 2 + norm2_avg[1] ** 2 + norm2_avg[2] ** 2)
            if l1 > 0:
                norm1_avg[0] /= l1
                norm1_avg[1] /= l1
                norm1_avg[2] /= l1
            if l2 > 0:
                norm2_avg[0] /= l2
                norm2_avg[1] /= l2
                norm2_avg[2] /= l2

            # Make orthogonal
            for _ in range(10):
                dot = (
                    norm1_avg[0] * norm2_avg[0]
                    + norm1_avg[1] * norm2_avg[1]
                    + norm1_avg[2] * norm2_avg[2]
                )
                misf = 90.0 - math.acos(max(-1.0, min(1.0, dot))) * RAD_TO_DEG
                if abs(misf) <= 0.1:
                    break
                theta = misf * 0.5 * DEG_TO_RAD
                n1_0, n1_1, n1_2 = norm1_avg[0], norm1_avg[1], norm1_avg[2]
                n2_0, n2_1, n2_2 = norm2_avg[0], norm2_avg[1], norm2_avg[2]
                norm1_avg[0] = n1_0 - n2_0 * math.sin(theta)
                norm1_avg[1] = n1_1 - n2_1 * math.sin(theta)
                norm1_avg[2] = n1_2 - n2_2 * math.sin(theta)
                norm2_avg[0] = n2_0 - n1_0 * math.sin(theta)
                norm2_avg[1] = n2_1 - n1_1 * math.sin(theta)
                norm2_avg[2] = n2_2 - n1_2 * math.sin(theta)
                l1 = math.sqrt(norm1_avg[0] ** 2 + norm1_avg[1] ** 2 + norm1_avg[2] ** 2)
                l2 = math.sqrt(norm2_avg[0] ** 2 + norm2_avg[1] ** 2 + norm2_avg[2] ** 2)
                if l1 > 0:
                    norm1_avg[0] /= l1
                    norm1_avg[1] /= l1
                    norm1_avg[2] /= l1
                if l2 > 0:
                    norm2_avg[0] /= l2
                    norm2_avg[1] /= l2
                    norm2_avg[2] /= l2

            # Find max rotation from average (simplified)
            maxrot = 0.0
            imax = 0
            for k in range(nc):
                idx_k = active_idx[k]
                # Quick rotation estimate
                dot_n = (
                    norm1_avg[0] * norm1in[0, idx_k]
                    + norm1_avg[1] * norm1in[1, idx_k]
                    + norm1_avg[2] * norm1in[2, idx_k]
                )
                dot_s = (
                    norm2_avg[0] * norm2in[0, idx_k]
                    + norm2_avg[1] * norm2in[1, idx_k]
                    + norm2_avg[2] * norm2in[2, idx_k]
                )
                dot_n = abs(dot_n)
                dot_s = abs(dot_s)
                dot_n = max(-1.0, min(1.0, dot_n))
                dot_s = max(-1.0, min(1.0, dot_s))
                rot = 0.5 * (math.acos(dot_n) + math.acos(dot_s)) * RAD_TO_DEG
                if rot > maxrot:
                    maxrot = rot
                    imax = k

            if maxrot <= cangle:
                break

            # Remove outlier
            active[active_idx[imax]] = 0
            nc -= 1
            if nc < 1:
                break

            # Rebuild active indices
            idx = 0
            for i in range(nf):
                if active[i] == 1:
                    active_idx[idx] = i
                    idx += 1

        # Count remaining
        nc = 0
        for i in range(nf):
            if active[i] == 1:
                nc += 1

        prob[imult] = float(nc) / float(nf)

        if imult > 0 and prob[imult] < prob_max:
            nsltn = imult
            break

        # Calculate RMS differences
        rms1 = 0.0
        rms2 = 0.0
        for i in range(nf):
            d11 = abs(
                norm1in[0, i] * norm1_avg[0]
                + norm1in[1, i] * norm1_avg[1]
                + norm1in[2, i] * norm1_avg[2]
            )
            d22 = abs(
                norm2in[0, i] * norm2_avg[0]
                + norm2in[1, i] * norm2_avg[1]
                + norm2in[2, i] * norm2_avg[2]
            )
            d11 = max(-1.0, min(1.0, d11))
            d22 = max(-1.0, min(1.0, d22))
            a11 = math.acos(d11)
            a22 = math.acos(d22)
            rms1 += a11 * a11
            rms2 += a22 * a22

        rms_diff[0, imult] = RAD_TO_DEG * math.sqrt(rms1 / nf)
        rms_diff[1, imult] = RAD_TO_DEG * math.sqrt(rms2 / nf)

        # Convert to strike, dip, rake
        s, d, r = fp_coord_vectors_to_angles(norm1_avg, norm2_avg)
        str_avg[imult] = s
        dip_avg[imult] = d
        rak_avg[imult] = r

        # Mark current cluster as done, continue with remaining
        for i in range(nf):
            if active[i] == 1:
                active[i] = 2  # Mark as processed

        # Check if any unprocessed remain
        has_remaining = False
        for i in range(nf):
            if active[i] == 0:
                active[i] = 1
                has_remaining = True

        if not has_remaining:
            nsltn = imult + 1
            break

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
