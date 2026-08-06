"""What is actually in a deposited structure, and what each part is.

A PIEZO1 coordinate file is not just the channel. Across the 23 entries this
project ships, the files also contain:

* **MDFIC**, a genuine auxiliary subunit — three copies of a 21-residue
  cysteine-rich peptide (``CCESSDCLEICMECCGICFPS``) that inserts into the pore
  and modulates gating (Zhou et al. 2023). Present in six entries.
* **Poly-UNK** — 6B3R carries three 16-residue chains of unassigned density.
* **Lipids and detergent** — PLX, PEE, P5S, L9Q are phospholipids; D12 is
  dodecane. Up to 1,092 atoms in one entry.
* **Glycan** — NAG on PIEZO2 (6KG7).

**Why this module exists rather than a `hetero` flag.** Those extra chains are
*protein*, so a protein mask includes them, and their residue numbers (226–247
for MDFIC) sit inside PIEZO1's own numbering. Nothing in this project has been
wrong because of it — PIEZO1 resolves from residue 570 upward in exactly the
entries that carry MDFIC, so the ranges happen not to overlap — but that is
luck, not design. A selection keyed on residue number alone would silently pool
MDFIC with PIEZO1 the moment an entry resolved further into the N-terminus.

So every atom gets a category, the principal polymer is identified by size
*relative to the largest chain* rather than an absolute cut-off (4RAX is only
227 residues and is the whole structure), and analyses say which categories they
included.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .structure import AA3TO1, Structure

__all__ = ["EntityClass", "EntityMap", "classify", "LIGAND_NAMES",
           "PRINCIPAL_FRACTION"]


class EntityClass:
    """Categories an atom can belong to. Plain strings, so they serialise."""

    PROTOMER = "protomer"
    AUXILIARY = "auxiliary"
    UNKNOWN_POLYMER = "unknown_polymer"
    LIPID = "lipid"
    DETERGENT = "detergent"
    GLYCAN = "glycan"
    ION = "ion"
    WATER = "water"
    OTHER = "other"

    ALL = (PROTOMER, AUXILIARY, UNKNOWN_POLYMER, LIPID, DETERGENT, GLYCAN,
           ION, WATER, OTHER)

    #: What each category is, for the UI and for reports.
    LABELS = {
        PROTOMER: "Channel protomer",
        AUXILIARY: "Auxiliary protein subunit",
        UNKNOWN_POLYMER: "Unassigned polymer (poly-UNK)",
        LIPID: "Lipid",
        DETERGENT: "Detergent",
        GLYCAN: "Glycan",
        ION: "Ion",
        WATER: "Water",
        OTHER: "Other heterogen",
    }


#: Chemical component names, taken from the PDB chemical component dictionary
#: rather than inferred from the three-letter code. Every one of these appears
#: in a structure this project ships.
LIGAND_NAMES = {
    "L9Q": ("lipid", "phosphatidylethanolamine (C41 H80 N O8 P)"),
    "PLX": ("lipid", "phosphatidylcholine derivative (C42 H89 N O8 P)"),
    "PEE": ("lipid", "1,2-dioleoyl-sn-glycero-3-phosphoethanolamine"),
    "P5S": ("lipid", "phosphatidylserine (C42 H82 N O10 P)"),
    "D12": ("detergent", "dodecane"),
    "NAG": ("glycan", "2-acetamido-2-deoxy-beta-D-glucopyranose"),
}

#: Common codes that may appear in other entries. Kept short and honest: an
#: unrecognised heterogen is reported as `other`, not guessed at.
ION_CODES = {"NA", "K", "CA", "MG", "ZN", "CL", "MN", "FE", "CU", "CD", "NI"}
WATER_CODES = {"HOH", "WAT", "DOD"}
GLYCAN_CODES = {"NAG", "BMA", "MAN", "FUC", "GAL", "GLC", "BGC", "SIA"}
LIPID_HINTS = ("PC", "PE", "PS", "PG", "PI", "PA", "CDL", "CLR", "POV", "OLA",
               "OLC", "PLM", "STE", "MYR", "LMT", "LMN")

#: A chain counts as a principal polymer if it has at least this fraction of the
#: longest chain's residues. Relative rather than absolute because 4RAX is a
#: 227-residue domain and is the entire structure, while the 21-residue MDFIC
#: chains sit alongside 1,280-residue protomers.
PRINCIPAL_FRACTION = 0.5


@dataclass
class EntityMap:
    """Per-atom classification of one structure."""

    categories: np.ndarray            # (n_atoms,) of EntityClass strings
    chain_class: dict = field(default_factory=dict)
    residue_names: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def mask(self, *classes: str) -> np.ndarray:
        return np.isin(self.categories, list(classes))

    @property
    def protomer_chains(self) -> list[str]:
        return sorted(c for c, k in self.chain_class.items()
                      if k == EntityClass.PROTOMER)

    @property
    def auxiliary_chains(self) -> list[str]:
        return sorted(c for c, k in self.chain_class.items()
                      if k in (EntityClass.AUXILIARY,
                               EntityClass.UNKNOWN_POLYMER))

    def counts(self) -> dict:
        unique, totals = np.unique(self.categories, return_counts=True)
        return {str(k): int(v) for k, v in zip(unique, totals)}

    def present(self) -> list[str]:
        """Categories actually present, in a stable display order."""
        counts = self.counts()
        return [k for k in EntityClass.ALL if counts.get(k)]

    def summary(self) -> str:
        counts = self.counts()
        parts = [f"{EntityClass.LABELS[k].lower()} {counts[k]:,}"
                 for k in EntityClass.ALL if counts.get(k)]
        text = " · ".join(parts)
        if self.auxiliary_chains:
            text += (f" · auxiliary chains "
                     f"{', '.join(self.auxiliary_chains)} excluded from "
                     f"channel analyses")
        return text


def _classify_residue(name: str) -> str:
    name = str(name).strip().upper()
    if name in WATER_CODES:
        return EntityClass.WATER
    if name in LIGAND_NAMES:
        kind = LIGAND_NAMES[name][0]
        return {"lipid": EntityClass.LIPID,
                "detergent": EntityClass.DETERGENT,
                "glycan": EntityClass.GLYCAN}[kind]
    if name in GLYCAN_CODES:
        return EntityClass.GLYCAN
    if name in ION_CODES:
        return EntityClass.ION
    if any(hint in name for hint in LIPID_HINTS):
        return EntityClass.LIPID
    return EntityClass.OTHER


def classify(structure: Structure) -> EntityMap:
    """Assign every atom a category.

    Polymer chains are split into principal (channel protomer) and auxiliary by
    length relative to the longest chain; heterogens are classified by residue
    code. Anything unrecognised becomes ``other`` rather than being guessed at.
    """
    n = structure.n_atoms
    categories = np.full(n, EntityClass.OTHER, dtype=object)
    chain_class: dict[str, str] = {}

    polymer = structure.mask_protein() & ~structure.hetero
    lengths: dict[str, int] = {}
    for chain in structure.chains:
        in_chain = structure.chain == chain
        lengths[chain] = int((in_chain & polymer & structure.mask_ca()).sum())
    longest = max(lengths.values()) if lengths else 0

    for chain in structure.chains:
        in_chain = structure.chain == chain
        length = lengths[chain]
        if length == 0:
            continue
        if longest and length >= PRINCIPAL_FRACTION * longest:
            kind = EntityClass.PROTOMER
        else:
            names = set(structure.res_name[in_chain & polymer].tolist())
            kind = (EntityClass.UNKNOWN_POLYMER
                    if names and names <= {"UNK"} else EntityClass.AUXILIARY)
        chain_class[str(chain)] = kind
        categories[in_chain & polymer] = kind

    residue_names: dict[str, str] = {}
    hetero = structure.hetero | ~structure.mask_protein()
    for name in np.unique(structure.res_name[hetero]):
        kind = _classify_residue(name)
        categories[hetero & (structure.res_name == name)] = kind
        residue_names[str(name)] = kind

    return EntityMap(
        categories=categories, chain_class=chain_class,
        residue_names=residue_names,
        meta={"longest_chain_ca": longest,
              "chain_lengths": lengths,
              "principal_fraction": PRINCIPAL_FRACTION,
              "n_protomers": sum(1 for k in chain_class.values()
                                 if k == EntityClass.PROTOMER)})
