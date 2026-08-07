"""What the deposited variant structures can and cannot tell us.

Round 34 set out to compare ion permeation across the four deposited PIEZO1
variant structures and read a direction of change against the measured
phenotype. The survey below is what stopped that, and it is worth having as a
standing, checkable fact rather than a note in a log.

Three things, each measured here rather than asserted:

1. **Every deposited human PIEZO1 structure is closed.** The pore is too narrow
   for a cation in all of them, so every conductance is exactly zero and a
   *difference* in conductance cannot be measured. There is no open variant
   structure to compare against an open wild type.

2. **Most of the variants are not in their own structures.** A1988 is unmodelled
   in both entries named for A1988V, and E756 is unmodelled in the E756del
   entry. Only R2456H (8YFG) actually shows its mutation — arginine everywhere
   else, histidine there.

3. **Three of the entries share one model.** 8ZU3 (wild type + MDFIC), 8YFC
   (A1988V + MDFIC) and 9VMX (E756del + MDFIC) have byte-identical protein
   coordinates: 31,839 atoms, the same hash, 0.000 A RMSD. They are separate
   depositions with separate titles and separate files, so this is a fact about
   the depositions rather than about our download — but it means they cannot
   distinguish anything from each other.

Taken together: of four nominal variant structures, **one** carries its variant
and is coordinate-distinct. That is the honest denominator against sixty-eight
curated variants, and it is the same data limit Round 22 met from the other
side — there, not enough phenotyped variants; here, not enough structures.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

from ..core.structure import Structure

__all__ = ["VariantStructure", "StructureSurvey", "survey_variant_structures",
           "coordinate_fingerprint", "VARIANT_ENTRIES"]

#: The deposited human entries this survey covers, and which variant each is
#: named for. ``None`` marks a wild-type control rather than a variant.
VARIANT_ENTRIES = {
    "8YEZ": None,
    "8ZU3": None,
    "8ZU8": ("A1988V", 1988, "A", "V"),
    "8YFC": ("A1988V", 1988, "A", "V"),
    "8YFG": ("R2456H", 2456, "R", "H"),
    "9VMX": ("E756del", 756, "E", "del"),
}


def coordinate_fingerprint(structure: Structure) -> str:
    """A hash of the protein coordinates, for spotting identical depositions.

    Rounded to 1e-3 A — the precision coordinates are deposited at — so that two
    entries differing only in file formatting still hash the same, while any
    real difference in a single atom does not.
    """
    mask = structure.mask_protein() & ~structure.hetero
    return hashlib.sha1(
        np.round(structure.xyz[mask], 3).tobytes()).hexdigest()[:12]


@dataclass
class VariantStructure:
    """One deposited entry, and what it is actually able to show."""

    pdb: str
    variant: str | None
    residue: int | None
    expected_wt: str | None
    expected_mut: str | None
    observed: str | None                 # residue name actually modelled
    fingerprint: str = ""
    n_protein_atoms: int = 0
    bottleneck_A: float = float("nan")
    wetting_score: float = float("nan")
    conductance_pS: float = float("nan")
    mechanisms: list = field(default_factory=list)
    duplicates: tuple = ()

    @property
    def is_control(self) -> bool:
        return self.variant is None

    @property
    def mutation_resolved(self) -> bool:
        """Is the mutated residue modelled at all?

        A structure that does not resolve its own mutation cannot show what the
        mutation does, however good it is elsewhere.
        """
        return self.observed is not None

    @property
    def shows_mutation(self) -> bool:
        """Is the *mutant* residue actually present, rather than the wild type?"""
        from ..core.structure import AA3TO1
        if self.observed is None or self.expected_mut in (None, "del"):
            return False
        return AA3TO1.get(self.observed, "?") == self.expected_mut

    @property
    def informative(self) -> bool:
        """Can this entry say anything about its variant?"""
        return (not self.is_control and self.shows_mutation
                and not self.duplicates)

    def summary(self) -> str:
        if self.is_control:
            return f"{self.pdb}: wild-type control"
        state = ("shows the mutation" if self.shows_mutation
                 else f"residue {self.residue} unmodelled"
                 if not self.mutation_resolved
                 else f"residue {self.residue} modelled as wild type")
        extra = (f"; coordinates identical to {', '.join(self.duplicates)}"
                 if self.duplicates else "")
        return f"{self.pdb} ({self.variant}): {state}{extra}"


@dataclass
class StructureSurvey:
    """The whole table, plus what it means for the comparison Round 34 wanted."""

    entries: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def by_pdb(self, pdb: str) -> VariantStructure | None:
        return next((e for e in self.entries if e.pdb == pdb), None)

    @property
    def variants(self) -> list:
        return [e for e in self.entries if not e.is_control]

    @property
    def informative(self) -> list:
        return [e for e in self.entries if e.informative]

    @property
    def any_conducting(self) -> bool:
        return any(e.conductance_pS > 0 for e in self.entries)

    def duplicate_groups(self) -> list:
        """Entries sharing byte-identical protein coordinates."""
        groups: dict[str, list] = {}
        for entry in self.entries:
            groups.setdefault(entry.fingerprint, []).append(entry.pdb)
        return [sorted(v) for v in groups.values() if len(v) > 1]

    def coverage(self) -> dict:
        """The denominator, stated plainly."""
        from ..core.annotations import load_annotations
        curated = [v for v in load_annotations("human").variants
                   if v.classification in ("GoF", "LoF")]
        directions = {e.variant: _direction(e.variant) for e in self.variants}
        return {
            "curated_variants": len(curated),
            "deposited_variant_entries": len(self.variants),
            "resolve_their_own_mutation": sum(
                1 for e in self.variants if e.shows_mutation),
            "informative": len(self.informative),
            "distinct_variants_shown": len({e.variant for e in self.informative}),
            "directions_available": sorted(
                {d for d in directions.values() if d}),
            "all_closed": not self.any_conducting,
        }

    def summary(self) -> str:
        cov = self.coverage()
        return (f"{cov['deposited_variant_entries']} deposited variant entries; "
                f"{cov['resolve_their_own_mutation']} resolve their own "
                f"mutation; {cov['informative']} are informative. "
                f"All closed: {cov['all_closed']}. "
                f"Directions represented: {cov['directions_available']} "
                f"of ['GoF', 'LoF'] — against {cov['curated_variants']} "
                f"curated variants.")


def _direction(label: str | None) -> str | None:
    if not label:
        return None
    from ..core.annotations import load_annotations
    for variant in load_annotations("human").variants:
        if variant.label == label:
            return variant.classification
    return None


def survey_variant_structures(entries: dict | None = None,
                              step: float | None = None) -> StructureSurvey:
    """Measure every deposited variant entry and report what it can support.

    Runs the same pore, wetting and permeation pipeline the wild type gets, so
    that a difference between structures cannot come from a difference in
    treatment.
    """
    if step is None:
        step = _P.value("pore.step")
    from ..physics.permeation import solve_pnp
    from ..structure.frame import apply_frame, canonical_transform
    from ..structure.pore import pore_profile
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis
    from .hydration import load_grid, predict_wetting

    entries = VARIANT_ENTRIES if entries is None else entries
    from ..io.registry import load_registry
    registry = load_registry()
    grid = load_grid()

    out = []
    for pdb, spec in entries.items():
        record = registry.get(pdb)
        if record is None or not record.available:
            continue
        raw = Structure.from_file(record.path)
        structure = apply_frame(raw, canonical_transform(raw))

        variant = residue = wt = mut = None
        if spec is not None:
            variant, residue, wt, mut = spec

        observed = None
        if residue is not None:
            mask = (structure.mask_protein() & ~structure.hetero
                    & (structure.res_seq == residue))
            if mask.any():
                observed = str(structure.res_name[mask][0])

        blocks, _ = protomer_blocks(structure)
        profile = pore_profile(structure, detect_c3_axis(blocks), step=step)
        wetting = predict_wetting(structure, profile, grid) if grid.available else None
        result = solve_pnp(profile, wetting)

        out.append(VariantStructure(
            pdb=pdb, variant=variant, residue=residue, expected_wt=wt,
            expected_mut=mut, observed=observed,
            fingerprint=coordinate_fingerprint(raw),
            n_protein_atoms=int((raw.mask_protein() & ~raw.hetero).sum()),
            bottleneck_A=float(profile.bottleneck_radius),
            wetting_score=float(wetting.score) if wetting is not None else float("nan"),
            conductance_pS=float(result.conductance_pS),
            mechanisms=list(result.meta.get("mechanisms", []))))

    survey = StructureSurvey(entries=out)
    for group in survey.duplicate_groups():
        for entry in survey.entries:
            if entry.pdb in group:
                entry.duplicates = tuple(p for p in group if p != entry.pdb)
    survey.meta["note"] = (
        "Every deposited human PIEZO1 structure is closed, so no difference in "
        "conductance can be measured between them. Coverage is reported rather "
        "than worked around.")
    return survey
