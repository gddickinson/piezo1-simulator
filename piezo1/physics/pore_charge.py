"""Fixed charge on the pore wall, mapped onto the axial slices.

:func:`piezo1.physics.permeation.solve_pnp` has taken a ``fixed_charge``
argument since it was written and no caller ever supplied one. This is what
supplies it: the acidic and basic residues that line the conduction pathway,
turned into a signed charge density at every height the pore profile already
measures.

**Two routes, reported side by side.** ``curated`` uses only the residues this
project's own annotation calls pore-facing — the four selectivity glutamates,
the pore-lining set, the CTD constrictions — which is the literature's answer to
*which residues line the pore*. ``lining`` uses every ionisable residue whose
side chain could reach the lumen, which is the coordinates' answer. They do not
agree, and the disagreement is the interesting part: the curated set is a claim
about function, the geometric set a claim about position, and neither is
evidence for the other.

**Positions come from C-alpha, on every structure, deliberately.** 11ZC is the
only downloaded entry whose pore is open and the only one deposited without side
chains — it has N, CA, C and O and nothing else. So the charged group's own
coordinates do not exist where they are most needed. Measuring the other entries
from their real side chains and 11ZC from its backbone would make the two
incomparable, so every entry is measured the same way: the charge is placed at
its residue's C-alpha height, and whether the residue reaches the lumen at all
is decided by comparing its radial distance against the pore radius plus the
length of its own side chain fully extended. That criterion is deliberately
permissive — a residue it excludes cannot line the pore in *any* rotamer.

**The density is a space-charge density, and that is justified rather than
assumed.** Wall charge is spread uniformly across the lumen cross-section, which
is only legitimate when the screening length exceeds the pore radius. The
permeation module already measures both: 5.7-8.1 A against 3.3 A. The same
double-layer overlap that stops the Poisson coupling converging is what makes
the uniform-potential limit the right one here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.annotations import Annotations, load_annotations
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["ChargedGroup", "PoreCharge", "charged_groups", "map_charge",
           "pore_charge", "cytosolic_end", "CHARGE", "REACH_KEYS", "MODES",
           "AVOGADRO"]

#: Avogadro constant, 1/mol. Definitional since the 2019 SI redefinition; it
#: converts a count of elementary charges into a molar-equivalent density.
AVOGADRO = 6.02214076e23

#: Formal side-chain charge at physiological pH. Histidine is absent on
#: purpose: its pKa is near 6, so it is mostly neutral at pH 7.4 and including
#: it at +1 would be a decision about protonation that no measurement here
#: supports. Termini are absent for the same reason — the deposited chains are
#: fragments, so their ends are construct boundaries rather than real termini.
CHARGE = {"ASP": -1.0, "GLU": -1.0, "LYS": 1.0, "ARG": 1.0}

#: How far each side chain can put its charge from its own C-alpha, as registry
#: keys. Fully extended values, so the reach test errs towards inclusion.
REACH_KEYS = {"ASP": "pore_charge.reach_asp", "GLU": "pore_charge.reach_glu",
              "LYS": "pore_charge.reach_lys", "ARG": "pore_charge.reach_arg"}


def _reach(res_name: str) -> float:
    """This residue's side-chain reach in Angstrom, resolved at call time.

    The four keys are written out rather than looked up through
    :data:`REACH_KEYS`, and the duplication is deliberate.
    ``provenance_chain.resolved_keys`` finds a wired parameter by *scanning*
    for the call rather than running it, so a key reached only through a
    dictionary is reported as read by nothing — which is the chain's worst
    break, because such a parameter still appears in the dialog, still trips
    the non-default banner, and still stops ``verify_claims``, while doing
    nothing. All four of these were caught that way when they were added.
    ``test_pore_charge`` checks the two routes agree, so they cannot drift.
    """
    if res_name == "ASP":
        return _P.value("pore_charge.reach_asp")
    if res_name == "GLU":
        return _P.value("pore_charge.reach_glu")
    if res_name == "LYS":
        return _P.value("pore_charge.reach_lys")
    if res_name == "ARG":
        return _P.value("pore_charge.reach_arg")
    raise KeyError(f"{res_name} carries no side-chain charge")


#: Which residues a route considers. See the module docstring.
MODES = ("curated", "lining")


@dataclass(frozen=True)
class ChargedGroup:
    """One ionisable side chain, placed against the pore it may or may not line."""

    res_seq: int
    res_name: str
    chain: str
    z: float                     # A along the conduction axis
    radial: float                # A from the axis, at the C-alpha
    pore_radius: float           # A, the profile's radius at this height
    charge: float                # elementary charges, signed
    curated: str | None = None   # the annotation group that names it, if any

    @property
    def reach(self) -> float:
        return _reach(self.res_name)

    @property
    def margin(self) -> float:
        """How far the charge must travel from C-alpha to reach the lumen, A."""
        return self.radial - self.pore_radius

    @property
    def lines_pore(self) -> bool:
        return self.margin <= self.reach

    def label(self) -> str:
        return f"{self.res_name}{self.res_seq}{self.chain}"


@dataclass
class PoreCharge:
    """A fixed-charge density profile, and everything needed to disbelieve it."""

    z: np.ndarray                # A, the profile's own slices
    density: np.ndarray          # mol/m^3, signed (rho_fixed / F)
    groups: list[ChargedGroup] = field(default_factory=list)
    mode: str = "curated"
    meta: dict = field(default_factory=dict)

    @property
    def net_charge(self) -> float:
        """Total elementary charges included, signed."""
        return float(sum(g.charge for g in self.groups))

    @property
    def n_groups(self) -> int:
        return len(self.groups)

    @property
    def peak_density(self) -> float:
        """Largest magnitude of the density, mol/m^3."""
        return float(np.max(np.abs(self.density))) if len(self.density) else 0.0

    def residue_summary(self) -> list[dict]:
        """One row per residue number: how many copies line the pore, and where."""
        rows: dict[int, dict] = {}
        for g in self.groups:
            row = rows.setdefault(g.res_seq, {
                "res_seq": g.res_seq, "res_name": g.res_name,
                "curated": g.curated, "charge": g.charge, "copies": 0,
                "z_A": [], "margin_A": []})
            row["copies"] += 1
            row["z_A"].append(round(g.z, 1))
            row["margin_A"].append(round(g.margin, 2))
        return [rows[k] for k in sorted(rows)]

    def summary(self) -> str:
        return (f"{self.n_groups} charged groups ({self.mode}), net "
                f"{self.net_charge:+.0f} e, peak density "
                f"{self.peak_density / 1000.0:.2f} M-equivalent")


def cytosolic_end(structure: Structure, axis) -> int:
    """Which end of the profile is cytosolic: index ``0`` or ``-1``.

    A permeability ratio has a sign, and the sign is which side the dilute bath
    is on, so getting this backwards would report an anion-selective pore as a
    cation-selective one and the number would look entirely reasonable. It is
    therefore measured rather than assumed, by the same rule
    :mod:`piezo1.structure.frame` uses to orient a trimer: PIEZO1's C-terminal
    domain is cytosolic, so the end the last residues sit at is the cytosolic
    one. The canonical frame puts that at -z, and this returns ``0`` when it
    does — but it checks rather than trusting that the structure was framed.
    """
    from ..structure.frame import CTERM_FRACTION

    mask = structure.mask_ca()
    if not mask.any():
        raise ValueError("no C-alpha atoms to orient the axis with")
    xyz = structure.xyz[mask].astype(float)
    order = np.argsort(structure.res_seq[mask])
    n_tail = max(1, int(round(CTERM_FRACTION * len(order))))
    projection = axis.project(xyz)
    tail = float(np.mean(projection[order[-n_tail:]]))
    return 0 if tail < float(np.mean(projection)) else -1


#: The groups whose residues count as pore-facing **for the charge map**.
#:
#: Named rather than selected by ``category == "pore"``, which is what this did
#: until Round 84d. That coupling was invisible and load-bearing: adding four
#: curated groups for Liu et al. 2025's cap gate and spring linker — all of
#: them genuinely pore elements, all of them correctly categorised — took the
#: curated charge set from 6 to 12 and **flipped the measured selectivity from
#: cation- to anion-selective**, because R2279, D2310, E2318 and E2367 sit in
#: the cap rather than on the conduction pathway. An annotation edit must not
#: silently redefine a recorded measurement, so the set it feeds is written
#: down here and adding a group is now a deliberate act.
CURATED_CHARGE_GROUPS = ("hydrophobic_gate", "pore_lining",
                         "selectivity_acidic", "ctd_constriction")


def _curated_map(annotations: Annotations, species: str) -> dict[int, str]:
    """Residue number -> annotation group, for the groups called pore-facing.

    Numbering goes through the annotation's own human/mouse pair, which was
    built from a real alignment. The human-to-mouse offset is not constant, and
    a structure that is mouse must not be read with human numbers.
    """
    key = "mouse" if species == "mouse" else "human"
    out: dict[int, str] = {}
    for group in annotations.residue_groups:
        if group.id not in CURATED_CHARGE_GROUPS:
            continue
        for detail in group.detail:
            number = detail.get(key)
            if number is not None:
                out.setdefault(int(number), group.id)
    return out


def charged_groups(structure: Structure, profile, axis,
                   mode: str = "curated", species: str = "human",
                   annotations: Annotations | None = None
                   ) -> list[ChargedGroup]:
    """Every ionisable residue that reaches the lumen, in the chosen route.

    ``mode='curated'`` keeps only residues the annotation names as pore-facing;
    ``mode='lining'`` keeps every ionisable residue. Both then apply the reach
    test, so ``curated`` is a subset of ``lining`` by construction.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, not {mode!r}")
    curated = _curated_map(annotations or load_annotations(), species)

    mask = structure.mask_ca()
    xyz = structure.xyz[mask].astype(float)
    res_seq = structure.res_seq[mask]
    res_name = structure.res_name[mask]
    chain = structure.chain[mask]
    projection = axis.project(xyz)
    radial = axis.radial(xyz)

    z_profile = np.asarray(profile.z, dtype=float)
    r_profile = np.asarray(profile.radius, dtype=float)
    lo, hi = float(z_profile.min()), float(z_profile.max())

    out: list[ChargedGroup] = []
    for i in range(len(xyz)):
        name = str(res_name[i])
        if name not in CHARGE:
            continue
        number = int(res_seq[i])
        if mode == "curated" and number not in curated:
            continue
        z = float(projection[i])
        if z < lo or z > hi:
            continue
        group = ChargedGroup(
            res_seq=number, res_name=name, chain=str(chain[i]), z=z,
            radial=float(radial[i]),
            pore_radius=float(np.interp(z, z_profile, r_profile)),
            charge=CHARGE[name], curated=curated.get(number))
        if group.lines_pore:
            out.append(group)
    return sorted(out, key=lambda g: g.z)


def map_charge(groups: list[ChargedGroup], profile, mode: str = "curated",
               smoothing: float | None = None) -> PoreCharge:
    """Smear each charge along the axis and divide by the lumen it sits in.

    The Gaussian is not cosmetic. A charge's height is known to about a side
    chain's length — and in 11ZC only from C-alpha — so putting it on one slice
    would claim a precision the coordinates do not have, and would make the
    answer depend on the slice spacing. The kernel is normalised over the grid
    actually used, so the total charge is conserved however coarse that grid is.
    """
    smoothing = (_P.value("pore_charge.smoothing") if smoothing is None
                 else smoothing)
    z = np.asarray(profile.z, dtype=float)
    radius = np.asarray(profile.radius, dtype=float)
    per_length = np.zeros_like(z)                        # charges per A

    for group in groups:
        kernel = np.exp(-0.5 * ((z - group.z) / smoothing) ** 2)
        total = float(np.trapezoid(kernel, z))
        if total <= 0.0:
            continue
        per_length += group.charge * kernel / total

    # Below the ion radius nothing permeates anyway, and dividing by a vanishing
    # area would manufacture an unbounded density out of one carboxylate.
    floor = _P.value("pore.ion_radius")
    area = np.pi * np.maximum(radius, floor) ** 2        # A^2
    # charges/A^3 -> mol/m^3: 1e30 A^3 per m^3, divided by Avogadro.
    density = (per_length / area) * 1e30 / AVOGADRO

    return PoreCharge(
        z=z, density=density, groups=list(groups), mode=mode,
        meta={"smoothing_A": smoothing,
              "n_groups": len(groups),
              "net_charge_e": float(sum(g.charge for g in groups)),
              "radius_floor_A": floor,
              "note": "space-charge density: wall charge spread across the "
                      "lumen, valid because the screening length exceeds the "
                      "pore radius (see permeation.debye_length)"})


def pore_charge(structure: Structure, profile, axis, mode: str = "curated",
                species: str = "human",
                annotations: Annotations | None = None) -> PoreCharge:
    """The whole route in one call: find the charges, place them, report them."""
    groups = charged_groups(structure, profile, axis, mode=mode,
                            species=species, annotations=annotations)
    charge = map_charge(groups, profile, mode=mode)
    charge.meta["species"] = species
    charge.meta["structure"] = structure.name
    return charge
