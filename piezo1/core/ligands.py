"""Curated PIEZO1 modulators, and how much is known about where they bind.

Loads ``resources/ligands.json``. The resource exists mainly to make one thing
hard to forget: **no PIEZO structure with a bound small-molecule modulator has
been deposited**, so every binding site the field talks about is inferred from
mutagenesis, docking or geometry. A pocket drawn on a structure looks exactly
like a pocket that was observed in one, which is why the evidence level travels
with the site rather than sitting in a footnote.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

__all__ = ["Ligand", "LigandSet", "load_ligands", "SITE_EVIDENCE_ORDER"]

#: Strongest first. ``bound_structure`` is listed because its absence is the
#: point; nothing may claim it.
SITE_EVIDENCE_ORDER = ("bound_structure", "mutagenesis", "docking_md",
                       "geometric", "none")


@dataclass(frozen=True)
class Ligand:
    """One modulator, its potency, and the standing of its binding site."""

    key: str
    name: str
    role: str                      # activator | antagonist | inhibitor
    kind: str                      # small_molecule | peptide
    description: str = ""
    pubchem_cid: int | None = None
    inchikey: str | None = None
    chembl_id: str | None = None
    uniprot: str | None = None
    formula: str | None = None
    molecular_weight: float | None = None
    smiles: str | None = None
    potency: dict | None = None
    site_residues: tuple = ()
    site_numbering: str = "human"
    site_evidence: str = "none"
    site_note: str = ""
    site_citation: str | None = None

    @property
    def has_site(self) -> bool:
        return bool(self.site_residues)

    @property
    def site_is_observed(self) -> bool:
        """Always ``False``: no bound structure exists for any of these."""
        return self.site_evidence == "bound_structure"

    def potency_text(self) -> str:
        if not self.potency:
            return "no potency recorded"
        p = self.potency
        return (f"{p['measure']} {p['value']:g} {p['unit']} "
                f"({p['assay']}, {p['citation']})")

    def site_text(self) -> str:
        if not self.has_site:
            return f"no residue-level site — {self.site_note}"
        residues = ", ".join(str(r) for r in self.site_residues)
        return (f"{residues} ({self.site_numbering} numbering), evidence "
                f"'{self.site_evidence}' — INFERRED, not from a bound structure")

    def summary(self) -> str:
        return f"{self.name} ({self.role}): {self.potency_text()}"


@dataclass
class LigandSet:
    """Every curated modulator, plus what the set as a whole can claim."""

    ligands: list = field(default_factory=list)
    evidence_levels: dict = field(default_factory=dict)
    note: str = ""

    def __len__(self) -> int:
        return len(self.ligands)

    def get(self, key: str) -> Ligand | None:
        return next((x for x in self.ligands if x.key == key), None)

    def by_role(self, role: str) -> list:
        return [x for x in self.ligands if x.role == role]

    def with_sites(self) -> list:
        return [x for x in self.ligands if x.has_site]

    @property
    def any_observed_site(self) -> bool:
        return any(x.site_is_observed for x in self.ligands)

    def summary(self) -> str:
        return (f"{len(self.ligands)} curated modulators "
                f"({len(self.by_role('activator'))} activators, "
                f"{len(self.by_role('antagonist')) + len(self.by_role('inhibitor'))} "
                f"blockers); {len(self.with_sites())} carry a residue-level "
                f"site, none from a bound structure")


def load_ligands() -> LigandSet:
    from ..config import RESOURCE_DIR

    path = RESOURCE_DIR / "ligands.json"
    if not path.exists():
        return LigandSet(note="ligands.json not built; "
                              "run python scripts/build_ligands.py")
    data = json.loads(path.read_text())
    ligands = []
    for record in data["ligands"]:
        fields = {k: v for k, v in record.items()
                  if k in Ligand.__dataclass_fields__}
        fields["site_residues"] = tuple(fields.get("site_residues") or ())
        ligands.append(Ligand(**fields))
    return LigandSet(ligands=ligands,
                     evidence_levels=data.get("site_evidence_levels", {}),
                     note=data.get("note", ""))
