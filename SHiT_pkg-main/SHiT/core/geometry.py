"""
mdanalysis.core.geometry
========================
Shared geometry helpers: distances, angles, dihedrals, PBC wrapping.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Pure-Python helpers (no numpy dependency, used by bond-order analysis)
# ---------------------------------------------------------------------------

def distance_3d(p: Tuple, q: Tuple) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p, q)))


def angle_3d(p1: Tuple, p2: Tuple, p3: Tuple) -> float:
    """Angle (degrees) at vertex p2 for the triplet p1-p2-p3."""
    v1 = tuple(p1[k] - p2[k] for k in range(3))
    v2 = tuple(p3[k] - p2[k] for k in range(3))
    dot = sum(v1[k] * v2[k] for k in range(3))
    n1 = math.sqrt(sum(x ** 2 for x in v1))
    n2 = math.sqrt(sum(x ** 2 for x in v2))
    if n1 == 0 or n2 == 0:
        return 0.0
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (n1 * n2)))))


def _cross(a: Tuple, b: Tuple) -> Tuple:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Tuple, b: Tuple) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def dihedral_3d(p1: Tuple, p2: Tuple, p3: Tuple, p4: Tuple) -> float:
    """Dihedral angle (degrees) for p1-p2-p3-p4."""
    b0 = tuple(p2[k] - p1[k] for k in range(3))
    b1 = tuple(p3[k] - p2[k] for k in range(3))
    b2 = tuple(p4[k] - p3[k] for k in range(3))
    b1_len = math.sqrt(sum(x ** 2 for x in b1))
    if b1_len == 0:
        return 0.0
    b1n = tuple(x / b1_len for x in b1)
    v = _cross(b0, b1n)
    w = _cross(b2, b1n)
    return math.degrees(math.atan2(_dot(_cross(b1n, v), w), _dot(v, w)))


# ---------------------------------------------------------------------------
# NumPy-based helpers (used by RDF, water analysis, etc.)
# ---------------------------------------------------------------------------

def minimum_image(delta: np.ndarray, box: np.ndarray, pbc: np.ndarray) -> np.ndarray:
    """Apply minimum-image convention.

    Parameters
    ----------
    delta : (..., 3) displacement vectors
    box   : (3,) box lengths
    pbc   : (3,) boolean flags
    """
    delta = delta.copy()
    for i in range(3):
        if pbc[i] and box[i] > 0:
            delta[..., i] -= box[i] * np.round(delta[..., i] / box[i])
    return delta


def dist_pbc(pos1: np.ndarray, pos2: np.ndarray,
             box: np.ndarray, pbc: np.ndarray) -> float:
    """Scalar distance between two positions with PBC."""
    delta = minimum_image(pos2 - pos1, box, pbc)
    return float(np.linalg.norm(delta))


def angle_pbc(pos1: np.ndarray, pos2: np.ndarray, pos3: np.ndarray,
              box: np.ndarray, pbc: np.ndarray) -> float:
    """Angle (degrees) at pos2 with PBC."""
    v1 = minimum_image(pos1 - pos2, box, pbc)
    v2 = minimum_image(pos3 - pos2, box, pbc)
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-10 or n2 < 1e-10:
        return 0.0
    cos_a = np.clip(np.dot(v1 / n1, v2 / n2), -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_a)))


def dist_matrix_pbc(pos1: np.ndarray, pos2: np.ndarray,
                    box: np.ndarray, pbc: np.ndarray) -> np.ndarray:
    """Vectorised (n1, n2) distance matrix with PBC."""
    delta = pos2[np.newaxis, :, :] - pos1[:, np.newaxis, :]  # (n1, n2, 3)
    delta = minimum_image(delta, box, pbc)
    return np.sqrt(np.sum(delta ** 2, axis=2))


def wrap_coords(coords: np.ndarray, box_origin: np.ndarray,
                box_dims: np.ndarray, pbc: np.ndarray) -> np.ndarray:
    """Wrap coordinates into the primary simulation box."""
    wrapped = coords.copy()
    for i in range(3):
        if pbc[i] and box_dims[i] > 0:
            wrapped[:, i] = (
                box_origin[i]
                + (wrapped[:, i] - box_origin[i]) % box_dims[i]
            )
    return wrapped
