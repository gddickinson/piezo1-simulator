"""Numerical kernels for the PNP solver, kept apart from the physics.

Both of these exist because of a specific numerical failure, and both are pure
linear algebra with no notion of ions in them — so they are testable on their
own and do not belong in the middle of the transport equations.
"""

from __future__ import annotations

import numpy as np

__all__ = ["_bernoulli", "_solve_with_dirichlet"]


def _bernoulli(x: np.ndarray) -> np.ndarray:
    """B(x) = x / (exp(x) - 1), continuous at 0.

    Written piecewise because the naive form is 0/0 at x = 0 and loses all its
    significant figures nearby — which is the ordinary case here, since most of
    the pore has almost no field across a single grid step.
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-8
    out[small] = 1.0 - x[small] / 2.0
    # The asymptotes are handled explicitly: exp(x) overflows above ~709, and
    # the resulting inf propagates into the matrix as a nan rather than as the
    # zero it should be.
    large = x > 700.0
    out[large] = 0.0
    negative = x < -700.0
    out[negative] = -x[negative]
    rest = ~(small | large | negative)
    out[rest] = x[rest] / np.expm1(x[rest])
    return out


def _solve_with_dirichlet(matrix: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve after row-scaling, because the raw system is numerically singular.

    The interior equations are built from ``D * A / h`` with ``A`` a few square
    angstroms and ``h`` one angstrom, so their coefficients are around 1e-18,
    while a Dirichlet row is 1. That is a condition number near 1e18 and LAPACK
    reports the matrix as singular — not because it is, but because the two
    kinds of row are eighteen orders of magnitude apart. The interior rows are
    homogeneous, so dividing each by its own largest coefficient changes
    nothing about the solution and brings every row to order 1.
    """
    scale = np.max(np.abs(matrix), axis=1)
    scale[scale == 0.0] = 1.0
    return np.linalg.solve(matrix / scale[:, None], rhs / scale)
