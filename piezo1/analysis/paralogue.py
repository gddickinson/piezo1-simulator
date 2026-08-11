"""PIEZO2 as the control: how much of this mechanism is PIEZO1, and how much is the fold?

6KG7 is fetched, entity-classified, and then excluded from every ensemble as a
paralogue. That exclusion is right for a PIEZO1 ensemble and wrong as a final
answer, because PIEZO2 is the only thing available that can separate "PIEZO1
does this" from "a PIEZO does this".

**The registry said 6KG7 resolves residues 8-823. It resolves 8-2822.** 1,817
C-alphas per protomer in sixteen segments — *more resolved residues than any
PIEZO1 entry in the catalogue*, which resolve 1,223 to 1,502. The note was
written from the abstract's emphasis on the distal blade and nobody had asked
the file. Corrected in ``structures.json``; measured here so it cannot drift
back.

**Which protein and which numbering are measured, not assumed.** Every
deposited entry is scored, residue by residue, against all four committed
UniProt sequences: 6KG7 matches mouse Piezo2 at 1.000 and everything else below
0.25, and each PIEZO1 entry matches exactly one PIEZO1 reference the same way.
That matters because mouse Piezo2 is 2,822 aa against human PIEZO2's 2,752 and
mouse Piezo1's 2,547 — three lengths, no constant offset, and reading the wrong
one silently shifts every transmembrane helix.

**The naive comparison is a coverage artefact, and saying so is the result.**
Measured directly, PIEZO2's dome looks dramatically different from PIEZO1's:
8.5 nm deep against 4.9, with 462 nm² of excess area against 256. But 6KG7
resolves all 38 transmembrane helices where 7WLT resolves 22, so the two
measurements trace different amounts of blade. Restricted to the helices
resolved in both, PIEZO2 falls inside the PIEZO1 range on every quantity. The
difference was in what each file contains, not in what the proteins are.

**The gating coordinate is a property of the fold.** With the sites
coverage-matched through a real alignment and the protomer correspondence
searched rather than trusted — it is ``(2, 0, 1)``, so chain labels would have
been wrong — PIEZO1's lowest symmetric mode overlaps a single PIEZO2 symmetric
mode at **0.80**, and **0.93** of it lies in PIEZO2's symmetric subspace. The
motion this project identified as the candidate gating coordinate is not
specific to PIEZO1.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .numbering_check import (PIEZO2_REFERENCES, identify_numbering,
                              reference_entry as _reference)

__all__ = ["DomeRow", "ModeComparison", "paralogue_map",
           "tm_index_correspondence", "dome_comparison",
           "mode_comparison", "compare"]


def paralogue_map(piezo1: str = "mouse", piezo2: str = "mouse_piezo2"):
    """The PIEZO1 to PIEZO2 residue-number map, from a real global alignment.

    Not cached to a resource: unlike the human-mouse PIEZO1 map this is used
    by one analysis, takes a tenth of a second, and caching it would create a
    second place for the numbering to be wrong.
    """
    from ..core.sequence import NumberingMap

    return NumberingMap.from_sequences(
        _reference(piezo1)["sequence"], _reference(piezo2)["sequence"],
        piezo1, piezo2)


def tm_index_correspondence(numbering=None, piezo1: str = "mouse",
                            piezo2: str = "mouse_piezo2") -> dict:
    """Does TM *k* of PIEZO1 align onto TM *k* of PIEZO2?

    The coverage matching below pairs transmembrane helices **by index**, which
    is the architectural claim that both proteins have the same 38 helices in
    the same order. That claim is cheap to make and would quietly ruin the
    comparison if it were false, so it is checked against the alignment rather
    than assumed: take each PIEZO1 helix's midpoint, map it, and ask which
    PIEZO2 helix it lands in.
    """
    numbering = numbering or paralogue_map(piezo1, piezo2)
    first = _reference(piezo1)["transmembrane"]
    second = _reference(piezo2)["transmembrane"]
    slack = _P.value("paralogue.tm_boundary_slack")

    agree, rows = 0, []
    for index, helix in enumerate(first, start=1):
        middle = (helix["start"] + helix["end"]) // 2
        mapped = numbering.to_b(middle)
        landed = [] if mapped is None else [
            j for j, other in enumerate(second, start=1)
            if other["start"] - slack <= mapped <= other["end"] + slack]
        ok = landed == [index]
        agree += ok
        rows.append({"tm": index, "piezo1_mid": middle, "piezo2_residue": mapped,
                     "landed_in": landed, "agrees": ok})
    return {"n_helices": len(first), "n_agree": agree,
            "disagree": [r for r in rows if not r["agrees"]], "rows": rows,
            "identity": numbering.identity}


@dataclass
class DomeRow:
    """One dome measurement, with what it was measured over."""

    pdb: str
    reference: str
    n_helices: int
    radius_nm: float
    depth_nm: float
    area_nm2: float
    excess_nm2: float

    def summary(self) -> str:
        return (f"{self.pdb}: {self.n_helices} TM helices, R_c "
                f"{self.radius_nm:.2f} nm, depth {self.depth_nm:.2f} nm, "
                f"excess {self.excess_nm2:.0f} nm^2")


def _dome_row(pdb: str, structure: Structure, reference: str, keep=None) -> DomeRow:
    from ..structure.geometry import measure_dome, tm_surface_points
    from ..structure.protomers import protomer_blocks

    blocks, _ = protomer_blocks(structure)
    points, resolved = tm_surface_points(structure, reference, keep=keep)
    dome = measure_dome(blocks, points)
    return DomeRow(pdb=pdb, reference=reference,
                   n_helices=len(keep) if keep is not None else len(resolved),
                   radius_nm=dome.radius_of_curvature / 10.0,
                   depth_nm=dome.dome_depth / 10.0,
                   area_nm2=dome.dome_area / 100.0,
                   excess_nm2=dome.excess_area / 100.0)


def dome_comparison(piezo1_structure: Structure, piezo2_structure: Structure,
                    piezo1_pdb: str = "", piezo2_pdb: str = "") -> dict:
    """The dome of each, measured naively and then over the shared helices.

    Both rows are returned because the difference between them *is* the
    finding: the naive comparison says the two proteins have very different
    domes and it is measuring how much blade each file resolves.
    """
    from ..structure.geometry import tm_surface_points

    first = identify_numbering(piezo1_structure)
    second = identify_numbering(piezo2_structure)
    _, resolved_1 = tm_surface_points(piezo1_structure, first.reference)
    _, resolved_2 = tm_surface_points(piezo2_structure, second.reference)
    shared = sorted(resolved_1 & resolved_2)

    naive = [_dome_row(piezo1_pdb or piezo1_structure.name, piezo1_structure,
                       first.reference),
             _dome_row(piezo2_pdb or piezo2_structure.name, piezo2_structure,
                       second.reference)]
    matched = [_dome_row(piezo1_pdb or piezo1_structure.name, piezo1_structure,
                         first.reference, keep=set(shared)),
               _dome_row(piezo2_pdb or piezo2_structure.name, piezo2_structure,
                         second.reference, keep=set(shared))]
    return {"piezo1_numbering": first.reference,
            "piezo2_numbering": second.reference,
            "n_helices_piezo1": len(resolved_1),
            "n_helices_piezo2": len(resolved_2),
            "shared_helices": shared, "naive": naive, "coverage_matched": matched,
            "note": ("the naive rows differ mostly in how much blade each "
                     "entry resolves; the coverage-matched rows are the "
                     "comparison")}


@dataclass
class ModeComparison:
    """Whether PIEZO2's normal modes contain PIEZO1's gating coordinate."""

    n_sites: int
    protomer_order: tuple
    superposition_rmsd: float
    piezo1_symmetry: str
    piezo2_symmetry: str
    piezo1_gating_index: int
    best_overlap: float
    best_index: int
    best_symmetry: str
    symmetric_subspace_overlap: float
    shuffled_control: float
    meta: dict = field(default_factory=dict)

    @property
    def beats_control(self) -> bool:
        return self.best_overlap > 3.0 * self.shuffled_control

    def summary(self) -> str:
        return (f"PIEZO1's gating mode overlaps PIEZO2 mode {self.best_index} "
                f"({self.best_symmetry}) at {self.best_overlap:.3f}, "
                f"{self.symmetric_subspace_overlap:.3f} of it inside PIEZO2's "
                f"symmetric subspace; shuffled control "
                f"{self.shuffled_control:.3f}")


def mode_comparison(piezo1_structure: Structure, piezo2_structure: Structure,
                    n_modes: int | None = None, seed: int = 0
                    ) -> ModeComparison:
    """Coverage-matched elastic networks, and the overlap between their modes.

    Three things have to be right before the overlap means anything, and each
    has bitten this project or a neighbour of it:

    * the two site sets must **correspond**, which is what the alignment is
      for — a constant offset is certain to be wrong across a paralogue;
    * the **protomer order** must be searched rather than read off chain
      labels, and here it is not the identity;
    * the two structures must be in a **common frame**, so the PIEZO2 mode
      vectors are rotated by the same superposition that aligns the sites.

    The control shuffles the site correspondence and repeats the calculation:
    if a scrambled pairing gave the same overlap, the number would be saying
    something about elastic networks in general rather than about these two.
    """
    from ..physics.anm import ANM
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import kabsch, match_protomers

    n_modes = int(_P.value("anm.n_modes")) if n_modes is None else n_modes
    numbering = paralogue_map(identify_numbering(piezo1_structure).reference,
                              identify_numbering(piezo2_structure).reference)

    blocks_1, residues_1 = protomer_blocks(piezo1_structure)
    blocks_2, residues_2 = protomer_blocks(piezo2_structure)
    if not blocks_1 or not blocks_2:
        raise ValueError("both structures need three well-resolved protomers")
    residues_1 = np.asarray(residues_1)
    residues_2 = np.asarray(residues_2)
    available = set(residues_2.tolist())

    pairs = [(int(a), numbering.to_b(int(a))) for a in residues_1]
    pairs = [(a, b) for a, b in pairs if b is not None and b in available]
    if len(pairs) < 3:
        raise ValueError("too few corresponding residues to compare")

    index_1 = np.searchsorted(residues_1, [a for a, _ in pairs])
    index_2 = np.searchsorted(residues_2, [b for _, b in pairs])
    first = [block[index_1] for block in blocks_1]
    second = [block[index_2] for block in blocks_2]

    match = match_protomers(first, second)
    second = [second[i] for i in match.order]

    modes_1 = _modes(ANM.from_trimer(first), n_modes)
    modes_2 = _modes(ANM.from_trimer(second), n_modes)

    stacked_1, stacked_2 = np.vstack(first), np.vstack(second)
    rotation, translation, centroid = kabsch(stacked_2, stacked_1)
    fitted = (stacked_2 - centroid) @ rotation.T + translation
    rmsd = float(np.sqrt(((fitted - stacked_1) ** 2).sum(axis=1).mean()))

    symmetric_1 = [i for i, s in enumerate(modes_1.symmetry) if s == "A"]
    gating = symmetric_1[0] if symmetric_1 else 0
    reference_mode = modes_1.vectors[gating]

    overlaps = _overlaps(reference_mode, modes_2, rotation)
    best = int(np.argmax(overlaps))
    symmetric_share = float(np.sqrt(sum(
        overlaps[j] ** 2 for j, s in enumerate(modes_2.symmetry) if s == "A")))

    rng = np.random.default_rng(seed)
    order = rng.permutation(len(pairs))
    shuffled = [block[order] for block in second]
    control = float(np.max(_overlaps(reference_mode,
                                     _modes(ANM.from_trimer(shuffled), n_modes),
                                     rotation)))

    return ModeComparison(
        n_sites=len(pairs), protomer_order=tuple(int(i) for i in match.order),
        superposition_rmsd=rmsd,
        piezo1_symmetry="".join(modes_1.symmetry),
        piezo2_symmetry="".join(modes_2.symmetry),
        piezo1_gating_index=gating, best_overlap=float(overlaps[best]),
        best_index=best, best_symmetry=str(modes_2.symmetry[best]),
        symmetric_subspace_overlap=symmetric_share, shuffled_control=control,
        meta={"n_modes": n_modes, "alignment_identity": numbering.identity,
              "note": "site correspondence from a global alignment, protomer "
                      "order searched, PIEZO2 modes rotated into PIEZO1's frame"})


def _modes(anm, n_modes: int):
    built = anm.build()
    modes = built.calc_modes(n_modes=n_modes)
    built.label_symmetry(modes)
    return modes


def _overlaps(reference: np.ndarray, modes, rotation: np.ndarray) -> np.ndarray:
    """Cosine overlap of one mode with every mode of another network."""
    flat = reference.reshape(-1)
    norm = np.linalg.norm(flat)
    out = np.zeros(modes.n_modes)
    for j in range(modes.n_modes):
        rotated = (modes.vectors[j].reshape(-1, 3) @ rotation.T).reshape(-1)
        scale = norm * np.linalg.norm(rotated)
        out[j] = 0.0 if scale == 0 else abs(float(flat @ rotated) / scale)
    return out


def compare(piezo1_pdb: str = "7WLT", piezo2_pdb: str = "6KG7",
            n_modes: int | None = None) -> dict:
    """The whole comparison, as the report and the CLI serve it."""
    from ..io.registry import load_registry

    registry = load_registry()
    loaded = {}
    for pdb in (piezo1_pdb, piezo2_pdb):
        record = registry.get(pdb)
        if record is None or not record.available:
            return {"error": f"{pdb} is not downloaded"}
        loaded[pdb] = Structure.from_file(record.path)

    first, second = loaded[piezo1_pdb], loaded[piezo2_pdb]
    identity_1 = identify_numbering(first)
    identity_2 = identify_numbering(second)
    if identity_1.is_piezo2 or not identity_2.is_piezo2:
        return {"error": f"expected a PIEZO1 entry and a PIEZO2 entry, got "
                         f"{identity_1.reference} and {identity_2.reference}"}
    for pdb, identity in ((piezo1_pdb, identity_1), (piezo2_pdb, identity_2)):
        if not identity.confident:
            # 6LQI is the live case: annotation applied to it by canonical
            # residue number is wrong for more than half the chain, so a dome
            # measured on it would be tracing the wrong helices.
            return {"error": f"{pdb} is not in canonical numbering, so its "
                             f"annotation cannot be read by residue number "
                             f"({identity.summary()})"}

    modes = mode_comparison(first, second, n_modes=n_modes)
    return {
        "piezo1": {"pdb": piezo1_pdb, "identified": identity_1.summary(),
                   "confident": identity_1.confident},
        "piezo2": {"pdb": piezo2_pdb, "identified": identity_2.summary(),
                   "confident": identity_2.confident},
        "tm_correspondence": {k: v for k, v in
                              tm_index_correspondence().items() if k != "rows"},
        "dome": dome_comparison(first, second, piezo1_pdb, piezo2_pdb),
        "modes": modes,
        "note": ("PIEZO2 is the only control available for the question of "
                 "whether this mechanism is PIEZO1's or the fold's"),
    }
