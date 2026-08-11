"""Numerical kernels for the PNP solver, kept apart from the physics.

Two of these exist because of a specific numerical failure — :func:`_bernoulli`
and :func:`_solve_with_dirichlet` — and are pure linear algebra with no notion
of what an ion is. The rest are the discretisation itself and the coefficients
that go into it, moved here when Round 81 added fixed charge and
``permeation.py`` reached the project's length limit. The split is along a real
seam: everything here answers *how the equations are turned into a matrix*,
while ``permeation.py`` answers *which equations, and what the answer means*.
"""

from __future__ import annotations

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["_bernoulli", "_solve_with_dirichlet", "_nernst_planck",
           "_ohmic_potential", "_donnan_potential", "_face_conductance",
           "_neutrality_step", "_charge_diagnostics", "F_FARADAY", "R_GAS"]

#: Faraday constant, C/mol. Definitional since the 2019 SI redefinition.
F_FARADAY = 96485.33212
#: Molar gas constant, J/(mol K). Definitional since 2019.
R_GAS = 8.314462618

#: Largest exponent argument passed to ``exp`` in the Donnan solve. 40 thermal
#: voltages is 1 V at body temperature — far outside anything physical here —
#: and clipping there keeps a runaway iterate finite instead of nan.
_EXP_CLIP = 40.0


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


def _nernst_planck(z, area, potential, valence, diffusivity, thermal,
                   bath_left, bath_right):
    """Steady concentration profile and the constant flux it carries.

    Scharfetter-Gummel: the flux between two nodes is exact for a linear
    potential and constant coefficients across the interval, so the scheme does
    not go unstable when drift dominates diffusion. A centred difference does,
    and shows it as negative concentrations rather than as an error. That
    robustness is what lets ``potential`` carry the fixed-charge field, whose
    swings across one grid step are large where a carboxylate sits.

    ``potential`` is the **total** electrostatic potential seen by the ion —
    applied plus fixed-charge — and the two bath values are the concentrations
    imposed at the two ends, which differ when the baths do.

    **The drift term's sign was inverted here until Round 81.** The two
    Bernoulli factors were attached to the wrong nodes, so a cation drifted
    *up* the potential gradient. Nothing caught it for fifty rounds because
    every current this project had computed was between identical baths, where
    the concentration term is zero and reversing the field only reverses the
    current — and the sign was then discarded by ``pore_ohm =
    abs(voltage / current)``. It became visible the moment the two baths
    differed, as a pore that grew *more* anion-selective the more negative
    charge it was given. ``test_selectivity`` pins the direction on the case
    where no discretisation is involved at all: at a field weak enough that
    Scharfetter-Gummel reduces to a centred difference, the flux must equal
    ``-D A z c phi' / phi_T``, sign included.
    """
    n = len(z)
    h = np.diff(z)
    delta = valence * np.diff(potential) / thermal               # (n-1,)

    face_area = 0.5 * (area[:-1] + area[1:])
    a = diffusivity * face_area / h                              # (n-1,)

    lower = np.zeros(n)
    diag = np.zeros(n)
    upper = np.zeros(n)
    rhs = np.zeros(n)

    # J_k = a_k (c_k B(delta_k) - c_{k+1} B(-delta_k)): the upstream node is
    # weighted by B(delta), the downstream one by B(-delta).
    b_here = _bernoulli(delta)          # multiplies c_k
    b_next = _bernoulli(-delta)         # multiplies c_{k+1}

    for k in range(1, n - 1):
        left, right = a[k - 1], a[k]
        lower[k] = left * b_here[k - 1]
        diag[k] = -(left * b_next[k - 1] + right * b_here[k])
        upper[k] = right * b_next[k]

    diag[0] = diag[-1] = 1.0
    rhs[0], rhs[-1] = bath_left, bath_right

    matrix = np.zeros((n, n))
    matrix[np.arange(n), np.arange(n)] = diag
    matrix[np.arange(1, n), np.arange(n - 1)] = lower[1:]
    matrix[np.arange(n - 1), np.arange(1, n)] = upper[:-1]

    concentration = _solve_with_dirichlet(matrix, rhs)
    flux = a[0] * (concentration[0] * b_here[0] - concentration[1] * b_next[0])
    return concentration, float(flux)


def _ohmic_potential(z, face_conductance, v_left, v_right):
    """Potential from ohmic current continuity, for a pore with no fixed charge.

    Where the double layers do not overlap, the pore interior is neutral and
    Poisson reduces to :math:`\\nabla\\cdot(\\sigma A \\nabla\\phi) = 0`, with
    the local conductivity set by the local concentrations. That is well-posed:
    a Laplace problem with positive coefficients, so the Gummel loop around it
    converges where the full Poisson coupling diverges.

    It closes the system only while the concentrations are uniform, which they
    are exactly when there is no fixed charge and the two baths match — the
    condition under which every number this project published before Round 81
    was computed. As soon as either fails, the closure is
    :func:`_neutrality_step` instead: a charged pore's potential is set by
    electroneutrality, not by Ohm's law, and this operator cannot see the fixed
    charge at all.
    """
    n = len(z)
    matrix = np.zeros((n, n))
    rhs = np.zeros(n)
    for k in range(1, n - 1):
        left, right = face_conductance[k - 1], face_conductance[k]
        matrix[k, k - 1] = left
        matrix[k, k] = -(left + right)
        matrix[k, k + 1] = right
    matrix[0, 0] = matrix[-1, -1] = 1.0
    rhs[0], rhs[-1] = v_left, v_right
    return _solve_with_dirichlet(matrix, rhs)


def _neutrality_step(valences, concentrations, fixed, thermal):
    """How far the potential must move to restore local electroneutrality, V.

    The closure for a pore that carries charge. At each interior node it asks:
    if the ions here responded to a potential shift by Boltzmann, what shift
    would make :math:`\\sum_i z_i c_i + X` vanish? That is the same equation
    :func:`_donnan_potential` solves, with the *current iterate's*
    concentrations in place of the bath, so the two share their arithmetic.

    The ends return zero because they are pinned: the applied voltage and the
    boundary Donnan jump fix them, and their concentrations already satisfy
    electroneutrality by construction.

    **Why not current continuity.** Adding a diffusion-current source term to
    :func:`_ohmic_potential` is the obvious extension and it does not work
    here. The source has to be discretised centrally while the flux is
    discretised Scharfetter-Gummel, and where a carboxylate sits the potential
    changes by four thermal voltages across one angstrom — at which point the
    two discretisations disagree by more than the term itself. Measured: the
    electroneutrality residual reached 43 times the ionic content and the
    Gummel loop did not converge in 200 iterations. Electroneutrality imposed
    directly is both the physically correct closure in this limit and the one
    that converges.
    """
    step = _donnan_potential(valences, np.asarray(concentrations, dtype=float),
                             fixed, thermal)
    step[0] = step[-1] = 0.0
    return step


def _donnan_potential(valences, bath, fixed, thermal, max_iterations=60,
                      tol=1e-12):
    """Local Donnan potential, in volts, that a fixed charge imposes.

    Solves, at each node independently,

    .. math:: \\sum_i z_i c_i^\\text{bath}
              e^{-z_i \\psi / \\phi_T} + X = 0

    which is local electroneutrality between the bath ions, partitioned by
    Boltzmann, and the fixed charge ``X``. The left side is strictly decreasing
    in psi, so the root is unique and Newton from zero reaches it monotonically
    once the step is capped — which is what the clip below does.

    **Zero fixed charge returns exactly zero**, not something near it: the
    residual at psi = 0 is the bath's own electroneutrality, which is zero by
    construction, so Newton terminates on its first test and every partition
    factor is exactly ``exp(0) = 1``.
    """
    valences = np.asarray(valences, dtype=float)
    bath = np.asarray(bath, dtype=float)                     # (n_species, n)
    fixed = np.asarray(fixed, dtype=float)
    psi = np.zeros_like(fixed)
    scale = np.maximum(np.abs(fixed),
                       (np.abs(valences)[:, None] * bath).sum(axis=0))

    def residual(p):
        arg = np.clip(-valences[:, None] * p[None, :] / thermal,
                      -_EXP_CLIP, _EXP_CLIP)
        partitioned = bath * np.exp(arg)
        return (valences[:, None] * partitioned).sum(axis=0) + fixed, partitioned

    for _ in range(max_iterations):
        value, partitioned = residual(psi)
        if np.all(np.abs(value) <= tol * np.maximum(scale, 1.0)):
            break
        derivative = -(valences[:, None] ** 2 * partitioned).sum(axis=0) / thermal
        step = np.where(derivative != 0.0, value / derivative, 0.0)
        # A Newton step on a monotone, convex-in-exponential residual can
        # overshoot when the fixed charge dwarfs the bath; capping it at a few
        # thermal voltages keeps it monotone without changing the root.
        psi = psi - np.clip(step, -4.0 * thermal, 4.0 * thermal)
    return psi


def _face_conductance(z, areas, concentrations, species, temperature):
    """Ohmic conductance of each interval, S — the operator's coefficients."""
    h = np.diff(z)
    out = np.zeros(len(z) - 1)
    for s in species:
        local = ((s.valence ** 2 * F_FARADAY ** 2 * s.diffusivity)
                 / (R_GAS * temperature)) * concentrations[s.name] * areas[s.name]
        out += 0.5 * (local[:-1] + local[1:]) / h
    return out


def _charge_diagnostics(concentrations, species, fixed, psi) -> dict:
    """What the fixed charge did, and whether the answer is still physical.

    Two numbers decide that. Local electroneutrality is imposed at the mouths
    and ought to hold inside; the residual says whether it does. And the
    counterion concentration the Donnan partition demands can exceed anything
    ions could physically be packed to — at which point the continuum model has
    stopped describing a solution, and the result is reported with that said
    rather than quietly used.
    """
    net = sum(s.valence * concentrations[s.name] for s in species)
    scale = sum(abs(s.valence) * concentrations[s.name] for s in species)
    residual = np.abs(net + (0.0 if fixed is None else fixed)) / np.maximum(scale, 1e-30)
    peak = max(float(np.max(concentrations[s.name])) for s in species) / 1000.0
    ceiling = _P.value("pore_charge.max_concentration")
    return {"electroneutrality_residual": float(np.max(residual)),
            "donnan_potential_mV": [float(np.min(psi) * 1e3),
                                    float(np.max(psi) * 1e3)],
            "peak_in_pore_M": peak,
            "exceeds_packing_limit": bool(peak > ceiling),
            "packing_limit_M": ceiling,
            "net_fixed_charge_mol_per_m3": (0.0 if fixed is None
                                            else float(np.sum(fixed)))}
