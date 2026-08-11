"""Surface electrostatic potential — Guo & MacKinnon 2017, Figure 4c.

Figure 4c is a surface of the trimer coloured by electrostatic potential "in
aqueous solution containing 150 mM NaCl, calculated using APBS", saturating at
±5 k_BT/e. It is the panel behind two of the paper's structural arguments: the
positive-inside distribution on the arms (von Heijne's rule), and the charged
patches that hold the cap onto the extracellular loops.

**This is not APBS, and the difference matters.** APBS solves the Poisson-
Boltzmann equation on a grid with a dielectric boundary: the protein interior
is low-dielectric (ε ≈ 2-4) and the solvent high (ε ≈ 78.5), and the jump
between them is most of the physics. What this module computes is the
**linear-superposition Debye-Hückel** potential — screened Coulomb from each
formal charge through a uniform solvent dielectric:

.. math::
    \\phi(\\mathbf{r}) = \\ell_B \\sum_i z_i
        \\frac{e^{-|\\mathbf{r}-\\mathbf{r}_i| / \\lambda_D}}{|\\mathbf{r}-\\mathbf{r}_i|}

in units of k_BT/e, with the Bjerrum length ℓ_B and Debye length λ_D. It is the
same approximation PyMOL and Chimera offer as "coulombic surface colouring",
and it is chosen here for a reason this project cares about: it has a closed
form, so it can be calibrated against an answer known independently, and a
120k-atom trimer takes seconds rather than a grid solve.

What it gets right is the **sign and the pattern** — which surfaces are acidic,
which basic, and where the patches are. What it gets wrong, systematically:

* it **under-estimates magnitudes**, because a low-dielectric interior focuses
  field lines into the solvent and a uniform ε = 78.5 does not;
* it has **no ion-exclusion layer**, so screening starts at the atom rather
  than a Stern radius out, again reducing magnitudes near the surface;
* it uses **formal charges** on ionisable side chains rather than force-field
  partial charges, so amide and hydroxyl dipoles contribute nothing.

All three push the same way, so a potential from here should be read as a lower
bound on |φ| and never quoted as an APBS number.
:func:`piezo1.physics.electrostatic_patches.compare_conventions` reports the
effect of the ones that can be varied.

Conventions: Angstrom for lengths, k_BT/e for potential (the unit Figure 4c's
colour bar is in), elementary charges for charge.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from ..parameters import PARAMETERS as _P
from ..structure.geometry import sphere_points
from .pore_charge import CHARGE

__all__ = ["ChargeSet", "SurfacePotential", "bjerrum_length", "debye_length",
           "formal_charges", "surface_points", "potential_at",
           "surface_potential", "TERMINAL_CHARGE"]


#: Elementary charge squared over 4 pi eps0, in Joule * metre:
#: ``(1.602176634e-19)^2 / (4 pi * 8.8541878128e-12)``. Used only to form the
#: Bjerrum length; every result downstream is in k_BT/e.
#:
#: This was written as ``2.307077552e-18`` when the module was first drafted —
#: the value in Joule * Angstrom, with the metre-to-Angstrom conversion then
#: applied a second time below. Nothing failed: the Bjerrum length came out at
#: 7.1e10 Angstrom, the Debye length underflowed to zero, every potential was
#: exactly zero, and the truncation-error check reported a flawless 0.000%
#: because it was comparing zero with zero. The single-point-charge
#: calibration is what said no.
_E2_OVER_4PI_EPS0 = 2.307077552e-28   # J * m
_BOLTZMANN = 1.380649e-23             # J/K, definitional since 2019

#: Formal charge on the chain termini. Applied to the N and C atoms of the
#: first and last modelled residue of each chain — which for a cryo-EM model is
#: usually *not* the real terminus, so this is off by default.
TERMINAL_CHARGE = {"n_terminus": 1.0, "c_terminus": -1.0}


def bjerrum_length(temperature: float | None = None,
                   dielectric: float | None = None) -> float:
    """Separation at which two elementary charges interact with 1 k_BT, Å.

    7.14 Å for water at 298.15 K, which is the standard value and what the
    calibration checks against.
    """
    if temperature is None:
        temperature = _P.value("electrostatics.temperature")
    if dielectric is None:
        dielectric = _P.value("electrostatics.dielectric_solvent")
    # e^2 / (4 pi eps0 eps_r k_B T), converted from metres to Angstrom.
    return float(_E2_OVER_4PI_EPS0 * 1e10
                 / (dielectric * _BOLTZMANN * temperature))


def debye_length(ionic_strength: float | None = None,
                 temperature: float | None = None,
                 dielectric: float | None = None) -> float:
    """Debye screening length, Å. 7.86 Å at 150 mM 1:1 salt in water.

    ``kappa^2 = 8 pi l_B N_A I`` for a symmetric 1:1 electrolyte, with the
    ionic strength in mol/L converted to particles per cubic Angstrom.
    Infinite at zero ionic strength, which is returned rather than raising:
    unscreened Coulomb is a legitimate thing to ask for.
    """
    if ionic_strength is None:
        ionic_strength = _P.value("electrostatics.ionic_strength")
    if ionic_strength <= 0:
        return float("inf")
    l_b = bjerrum_length(temperature, dielectric)
    # mol/L -> molecules per Angstrom^3: N_A * 1e-27 per (mol/L).
    number_density = ionic_strength * 6.02214076e23 * 1e-27
    kappa_squared = 8.0 * np.pi * l_b * number_density
    return float(1.0 / np.sqrt(kappa_squared))


# --------------------------------------------------------------------------
# Charges
# --------------------------------------------------------------------------

@dataclass
class ChargeSet:
    """Point charges: positions in Angstrom, charge in elementary units."""

    xyz: np.ndarray
    charge: np.ndarray
    label: tuple[str, ...] = ()

    @property
    def total(self) -> float:
        return float(self.charge.sum())

    def __len__(self) -> int:
        return len(self.charge)


def formal_charges(structure, include_termini: bool = False,
                   histidine: float | None = None) -> ChargeSet:
    """One point charge per ionisable side chain, at its charge centroid.

    Uses the same ``CHARGE`` table as :mod:`piezo1.physics.pore_charge` — Asp
    and Glu at −1, Lys and Arg at +1 — so the two modules cannot disagree
    about what PIEZO1's charge is. The centroid is over the carboxylate or
    guanidinium/amine atoms when they are resolved, and falls back to C-alpha
    when they are not, which matters here: 11ZC, the only open structure, is
    deposited without side chains at all.

    ``histidine`` is zero by default. At pH 7.4 a typical histidine is about
    10% protonated, and assigning it a full +1 would put 3 spurious positive
    charges per protomer on PIEZO1's surface; a fractional value can be passed
    to see how much that matters.
    """
    histidine = (_P.value("electrostatics.histidine_charge")
                 if histidine is None else histidine)
    #: The atoms carrying each formal charge, matching interactions.py.
    groups = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"),
              "LYS": ("NZ",), "ARG": ("NE", "NH1", "NH2")}
    if histidine:
        groups["HIS"] = ("ND1", "NE2")

    st = structure
    protein = st.mask_protein() & (~st.hetero)
    xyz_out: list[np.ndarray] = []
    q_out: list[float] = []
    labels: list[str] = []

    names = st.res_name
    for start, res_index in zip(st.res_first, range(st.n_residues)):
        stop = (st.res_first[res_index + 1] if res_index + 1 < len(st.res_first)
                else st.n_atoms)
        name = str(names[start])
        if name not in groups or not protein[start]:
            continue
        q = CHARGE.get(name, histidine if name == "HIS" else 0.0)
        if q == 0.0:
            continue
        block = slice(start, stop)
        wanted = np.isin(st.atom_name[block], groups[name])
        if wanted.any():
            position = st.xyz[block][wanted].mean(axis=0)
        else:
            alpha = st.atom_name[block] == "CA"
            if not alpha.any():
                continue
            position = st.xyz[block][alpha][0]
        xyz_out.append(position)
        q_out.append(float(q))
        labels.append(f"{st.chain[start]}/{name}{int(st.res_seq[start])}")

    if include_termini:
        for chain in st.chains:
            mask = protein & (st.chain == chain) & (st.atom_name == "CA")
            if not mask.any():
                continue
            order = np.argsort(st.res_seq[mask])
            positions = st.xyz[mask][order]
            for position, q, tag in ((positions[0], TERMINAL_CHARGE["n_terminus"], "Nterm"),
                                     (positions[-1], TERMINAL_CHARGE["c_terminus"], "Cterm")):
                xyz_out.append(position)
                q_out.append(q)
                labels.append(f"{chain}/{tag}")

    return ChargeSet(xyz=np.array(xyz_out, dtype=np.float64).reshape(-1, 3),
                     charge=np.array(q_out, dtype=np.float64),
                     label=tuple(labels))


# --------------------------------------------------------------------------
# The surface
# --------------------------------------------------------------------------

def surface_points(structure, probe: float | None = None,
                   n_points: int | None = None,
                   mask: np.ndarray | None = None
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Solvent-accessible surface as points, with the atom each belongs to.

    The Shrake-Rupley construction, but keeping the accessible sample points
    instead of counting them. The point count per atom is exactly what
    :func:`piezo1.analysis.measure.sasa` divides by, so the areas implied here
    reproduce that function's — which is the calibration, since a surface
    generator that silently kept buried points would still produce a plausible
    coloured picture.
    """
    if probe is None:
        probe = _P.value("sasa.probe_radius")
    if n_points is None:
        n_points = int(_P.value("sasa.n_points_fast"))
    st = structure
    sel = (st.mask_protein() & (~st.hetero)) if mask is None else np.asarray(mask)
    xyz = st.xyz[sel].astype(np.float64)
    radii = st.vdw_radii()[sel].astype(np.float64) + probe
    index = np.flatnonzero(sel)
    if len(xyz) == 0:
        raise ValueError("no atoms selected")

    sphere = sphere_points(n_points)
    tree = cKDTree(xyz)
    max_r = radii.max()
    r2 = radii * radii
    out_points: list[np.ndarray] = []
    out_atom: list[np.ndarray] = []

    for i in range(len(xyz)):
        neighbours = np.asarray(tree.query_ball_point(xyz[i], radii[i] + max_r))
        neighbours = neighbours[neighbours != i]
        candidates = xyz[i] + radii[i] * sphere
        if not len(neighbours):
            accessible = np.ones(len(sphere), bool)
        else:
            v = xyz[i] - xyz[neighbours]
            d2 = ((v * v).sum(axis=1)[None, :] + r2[i]
                  + (2.0 * radii[i]) * (sphere @ v.T))
            accessible = (d2 >= r2[neighbours][None, :]).all(axis=1)
        if accessible.any():
            out_points.append(candidates[accessible])
            out_atom.append(np.full(int(accessible.sum()), index[i]))

    if not out_points:
        return np.zeros((0, 3)), np.zeros(0, dtype=int)
    return np.vstack(out_points), np.concatenate(out_atom)


def potential_at(points: np.ndarray, charges: ChargeSet,
                 debye: float | None = None, bjerrum: float | None = None,
                 cutoff: float | None = None) -> np.ndarray:
    """Screened-Coulomb potential at each point, in k_BT/e.

    The core sum, deliberately kept free of any structure or surface concept so
    it can be driven with a single charge and checked against the closed form.

    ``cutoff`` bounds the neighbour search. Truncation always *reduces* |φ|
    because every omitted term has the sign of its own charge and the omitted
    set is not sign-balanced in general, so the error is bounded rather than
    cancelling; pass ``cutoff=None`` for the exact all-pairs sum, which is what
    the tests compare against.
    """
    if debye is None:
        debye = debye_length()
    if bjerrum is None:
        bjerrum = bjerrum_length()
    pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
    if len(charges) == 0:
        return np.zeros(len(pts))

    phi = np.zeros(len(pts))
    if cutoff is None:
        # Exact: every charge against every point. Fine for the calibration
        # cases and for a patch of a few hundred points.
        delta = pts[:, None, :] - charges.xyz[None, :, :]
        r = np.linalg.norm(delta, axis=-1)
        r = np.maximum(r, 1e-6)
        screen = np.exp(-r / debye) if np.isfinite(debye) else 1.0
        return bjerrum * (charges.charge[None, :] * screen / r).sum(axis=1)

    charge_tree = cKDTree(charges.xyz)
    point_tree = cKDTree(pts)
    pairs = point_tree.sparse_distance_matrix(charge_tree, cutoff,
                                              output_type="coo_matrix")
    r = np.maximum(pairs.data, 1e-6)
    screen = np.exp(-r / debye) if np.isfinite(debye) else np.ones_like(r)
    contribution = bjerrum * charges.charge[pairs.col] * screen / r
    np.add.at(phi, pairs.row, contribution)
    return phi


@dataclass
class SurfacePotential:
    """Figure 4c: potential on the solvent-accessible surface."""

    points: np.ndarray               # (n, 3) surface points, Angstrom
    potential: np.ndarray            # (n,) k_BT/e at each point
    atom_index: np.ndarray           # (n,) which atom each point belongs to
    #: Mean potential per atom, NaN for atoms with no accessible point.
    atom_potential: np.ndarray
    charges: ChargeSet
    debye_length: float
    bjerrum_length: float
    meta: dict = field(default_factory=dict)

    @property
    def scale(self) -> float:
        """Saturation of the colour map, k_BT/e. Figure 4c's is 5."""
        return float(self.meta.get("scale", _P.value("electrostatics.colour_scale")))

    def fraction_saturated(self) -> float:
        """Fraction of the surface beyond the colour scale.

        Worth looking at before reading a picture: a surface that is largely
        saturated is being displayed on the wrong scale, and one that is
        nowhere near it is being displayed on a scale that hides the structure.
        """
        return float(np.mean(np.abs(self.potential) >= self.scale))

    def summary(self) -> str:
        return (f"{len(self.points)} surface points | "
                f"phi {np.percentile(self.potential, 5):+.2f} to "
                f"{np.percentile(self.potential, 95):+.2f} k_BT/e (5-95%), "
                f"mean {self.potential.mean():+.2f} | "
                f"{100 * self.fraction_saturated():.0f}% beyond +-{self.scale:g}")


def surface_potential(structure, mask: np.ndarray | None = None,
                      ionic_strength: float | None = None,
                      histidine: float | None = None,
                      include_termini: bool = False,
                      n_points: int | None = None) -> SurfacePotential:
    """Compute Figure 4c for a structure.

    Charges come from the whole structure; the surface is restricted by
    ``mask``. That asymmetry is deliberate — a patch's potential depends on
    charges outside it, so masking the charges too would report the patch as
    though the rest of the protein were not there.
    """
    charges = formal_charges(structure, include_termini=include_termini,
                             histidine=histidine)
    debye = debye_length(ionic_strength)
    bjerrum = bjerrum_length()
    cutoff = _P.value("electrostatics.max_distance")

    points, atom_index = surface_points(structure, n_points=n_points, mask=mask)
    phi = potential_at(points, charges, debye=debye, bjerrum=bjerrum,
                       cutoff=cutoff)

    atom_phi = np.full(structure.n_atoms, np.nan)
    if len(points):
        totals = np.bincount(atom_index, weights=phi,
                             minlength=structure.n_atoms)
        counts = np.bincount(atom_index, minlength=structure.n_atoms)
        seen = counts > 0
        atom_phi[seen] = totals[seen] / counts[seen]

    return SurfacePotential(
        points=points, potential=phi, atom_index=atom_index,
        atom_potential=atom_phi, charges=charges,
        debye_length=debye, bjerrum_length=bjerrum,
        meta={"n_charges": len(charges), "net_charge": charges.total,
              "ionic_strength_M": (ionic_strength if ionic_strength is not None
                                   else _P.value("electrostatics.ionic_strength")),
              "cutoff_A": cutoff,
              "scale": _P.value("electrostatics.colour_scale"),
              "method": "linear-superposition Debye-Huckel, formal charges, "
                        "uniform solvent dielectric",
              "not_apbs": ("Figure 4c used APBS with a dielectric boundary. "
                           "This omits it and so under-estimates |phi|; read "
                           "the sign and the pattern, not the value."),
              "citation": "dolinsky2004"})


# The questions Figure 4-figure supplement 1 asks of this machinery — the
# charged patches, whether they attract, and how much the open conventions
# move the answer — live in :mod:`piezo1.physics.electrostatic_patches`. Split
# at the project's length limit and along a real seam: everything above is
# *how to compute a screened potential anywhere*, and everything there is a
# specific question about PIEZO1's surface.

