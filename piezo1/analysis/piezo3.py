"""The third vertebrate PIEZO, as a structure rather than an argument.

The census's central claim is that vertebrates have a third PIEZO gene the
databases largely missed — as old as PIEZO1 and PIEZO2, both duplications on the
jawed-vertebrate stem — that it is transcribed, spliced, tissue-patterned and
under purifying selection at its pore, and that it has kept the human residues
whose mutation causes disease. All of that is sequence evidence. Nobody has
recorded a current from it, and human piezo3 has been the pseudogene
``PIEZO1P2`` since before the primate radiation, so there is no human protein to
record from.

What there *is* is one AlphaFold model of the zebrafish protein. This module
runs it through the machinery this project already has, and the first thing it
does is say what that machinery can and cannot conclude from a predicted
monomer.

**Three cautions, all of them load-bearing.**

1. **It is a monomer, and half of what this project measures needs three.** The
   dome, the pore and the conduction verdict are all properties of a trimer, so
   a trimer has to be *built* — and :func:`piezo1.structure.assembly.borrowed_fraction`
   measures that 85–96% of the resulting departure from planarity is the
   template's arrangement rather than piezo3's. A dome radius measured this way
   is mostly a measurement of the entry it was built on. It is reported with
   that fraction beside it, never alone.
2. **It is a prediction.** The B-factor column is pLDDT, and the census's own
   figure audit found the piezo3 predictions no less confident *where it
   matters* than predictions of the two genes whose structures have been solved
   — which is a statement about the core, not about the blade.
3. **Two UniProt records exist for this gene and they are not in the same
   numbering.** The model is ``A0AB32U1Q1``; the census scored ``A0AC58GFC9``.
   They differ by a single inserted residue around position 2014, so everything
   before it agrees and everything after is off by one. Nothing here joins them
   by arithmetic — :func:`census_to_model` goes through an alignment, and a test
   pins that the offset is 0 before the indel and -1 after.

What the module can say without any of that: the *fold* comparison needs no
assembly at all, and it is the census's own structural result.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from ..config import STRUCTURE_DIR
from ..core.family import load_family_findings
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["MODEL_FILE", "MODEL_ACC", "CENSUS_ACC", "KeptPosition",
           "TemplateFit", "load_model", "census_to_model", "model_to_census",
           "kept_positions", "fold_comparison", "template_survey",
           "best_paralogue_template"]

#: The only structural representation of piezo3 that exists.
MODEL_FILE = "AF-A0AB32U1Q1-F1-model_v6.cif"
MODEL_ACC = "A0AB32U1Q1"
#: The record the census scored. One residue longer, and not interchangeable.
CENSUS_ACC = "A0AC58GFC9"


@dataclass(frozen=True)
class KeptPosition:
    """A human pathogenic pore position, carried onto the piezo3 model."""

    gene: str
    human_resi: int
    human_aa: str
    element: str
    census_resi: int
    census_aa: str
    model_resi: int | None
    model_aa: str | None
    resolved: bool
    distance_to_axis: float | None = None

    @property
    def kept(self) -> bool:
        """Does the model's own sequence carry the human residue here?"""
        return bool(self.model_aa) and self.model_aa == self.human_aa

    @property
    def agrees_with_census(self) -> bool:
        """Does the model record agree with the record the census scored?

        Not a formality. The two UniProt entries differ by an indel, so a
        position mapped wrongly would show up here as a residue mismatch — this
        is the check that the alignment-backed map is doing its job.
        """
        return bool(self.model_aa) and self.model_aa == self.census_aa


@dataclass(frozen=True)
class TemplateFit:
    """One candidate trimer to build piezo3 against, and what it costs."""

    template: str
    protein: str
    identity: float
    n_core: int
    n_corresponding: int
    placement_rmsd: float
    clashes: int
    borrowed_fraction: float | None
    ok: bool

    @property
    def clash_ratio(self) -> float:
        """Clashes relative to what a real deposited trimer scores (3-8)."""
        return self.clashes / 8.0

    def summary(self) -> str:
        borrowed = ("-" if self.borrowed_fraction is None
                    else f"{self.borrowed_fraction:.0%}")
        return (f"{self.template} ({self.protein}): identity {self.identity:.2f}, "
                f"{self.n_core} core residues fitted to {self.placement_rmsd:.2f} A, "
                f"{self.clashes:,} inter-protomer clashes, {borrowed} of the "
                f"non-planarity borrowed from the template")


@lru_cache(maxsize=1)
def load_model() -> Structure:
    """The zebrafish piezo3 AlphaFold monomer."""
    path = STRUCTURE_DIR / MODEL_FILE
    if not path.exists():
        raise FileNotFoundError(
            f"{MODEL_FILE} is not downloaded; run "
            f"`python -m piezo1.io.fetch` or fetch_alphafold('{MODEL_ACC}')")
    return Structure.from_file(path)


@lru_cache(maxsize=1)
def _census_alignment():
    """The map between the two zebrafish piezo3 UniProt records.

    Built from a real global alignment of the two sequences, not from the
    single-residue offset they happen to differ by. The offset is *measured* by
    the test rather than assumed by the code, which is the same rule this
    project applies to human and mouse PIEZO1.
    """
    import json

    from ..config import RESOURCE_DIR
    from ..core.family import load_constraint
    from ..core.sequence import NumberingMap

    model_seq = json.loads(
        (RESOURCE_DIR / "uniprot_zebrafish_piezo3.json").read_text())["sequence"]
    census_seq = load_constraint("piezo3").sequence
    return NumberingMap.from_sequences(census_seq, model_seq,
                                       CENSUS_ACC, MODEL_ACC)


def census_to_model(residue: int) -> int | None:
    """A residue number in the census's record, in the model's record."""
    return _census_alignment().a_to_b.get(int(residue))


def model_to_census(residue: int) -> int | None:
    """The inverse."""
    return _census_alignment().b_to_a.get(int(residue))


def kept_positions(structure: Structure | None = None) -> list[KeptPosition]:
    """The fourteen pathogenic pore positions, located on the piezo3 model.

    The census's finding is that piezo3 carries the identical residue at all
    fourteen. This checks that against the *model's own* sequence — a different
    UniProt record for the same gene — and, where the model resolves the
    position, reports how far it sits from the protein's own long axis, so a
    reader can see that the kept residues are in the pore rather than scattered.
    """
    structure = structure or load_model()
    findings = load_family_findings()
    sequence = _model_sequence()
    axis_distance = _axis_distances(structure)

    out = []
    for p in findings.pathogenic_pore:
        model_resi = census_to_model(p.piezo3_resi)
        model_aa = (sequence[model_resi - 1]
                    if model_resi and 1 <= model_resi <= len(sequence) else None)
        out.append(KeptPosition(
            gene=p.gene, human_resi=p.resi, human_aa=p.aa, element=p.element,
            census_resi=p.piezo3_resi, census_aa=p.piezo3_aa,
            model_resi=model_resi, model_aa=model_aa,
            resolved=model_resi in axis_distance,
            distance_to_axis=axis_distance.get(model_resi)))
    return out


def _model_sequence() -> str:
    import json

    from ..config import RESOURCE_DIR
    return json.loads(
        (RESOURCE_DIR / "uniprot_zebrafish_piezo3.json").read_text())["sequence"]


def _axis_distances(structure: Structure) -> dict[int, float]:
    """Per-residue distance from the monomer's own principal axis.

    A monomer has no three-fold axis — that is the whole difficulty with piezo3
    — so this is the long axis of its own inertia tensor. It is **not** a pore
    radius, which needs a trimer, and it does not separate the kept residues
    from anything: measured on this model they span 17–63 A with no clustering.

    Reported anyway, and labelled, because the alternative is worse. The obvious
    thing to do with fourteen conserved pore residues is to show them sitting
    together near an axis, and on a *monomer* that picture cannot be made
    honestly — a protomer's own long axis runs blade-tip to pore and is not the
    conduction axis. A number that visibly fails to separate them is a clearer
    statement of that than leaving the column out.
    """
    mask = structure.mask_ca()
    xyz = structure.xyz[mask]
    res = structure.res_seq[mask]
    centred = xyz - xyz.mean(axis=0)
    _, _, vh = np.linalg.svd(centred, full_matrices=False)
    axis = vh[0]
    along = centred @ axis
    perp = centred - np.outer(along, axis)
    d = np.sqrt((perp * perp).sum(axis=1))
    return {int(r): float(v) for r, v in zip(res, d)}


def fold_comparison(target_pdb: str = "6B3R"):
    """The census's own structural result, re-run here.

    Superposes the piezo3 monomer on a deposited PIEZO1 protomer **by the pore
    module alone** and measures where the blades land. The census reported
    3.86 A over 448 C-alpha against mouse Piezo1 6B3R; this uses this project's
    own domain boundaries, its own alignment and its own Kabsch, so agreement is
    a reproduction and disagreement is worth chasing.

    Needs no trimer, which is why it is the one piezo3 result carrying none of
    the assembly caveats.
    """
    from .core_periphery import compare

    target = Structure.from_file(STRUCTURE_DIR / f"{target_pdb}.cif")
    return compare(load_model(), target, "piezo3", target_pdb)


def template_survey(candidates: tuple = ("7WLT", "6B3R", "8YEZ", "6KG7", "9ZIS"),
                    structure: Structure | None = None) -> list[TemplateFit]:
    """Build the trimer against each candidate and report what each costs.

    :func:`piezo1.structure.assembly.best_template` picks "same protein first,
    then most residues resolved". piezo3 is nobody's same protein, so that rule
    falls through to whichever PIEZO trimer resolves most of itself — the worm
    PEZO-1 entry 9ZIS — and building on it gives **28%** sequence identity and
    thousands more inter-protomer clashes than building on a PIEZO1 trimer at
    **44%**. The rule is right for the case it was written for and wrong here,
    so this module chooses explicitly and shows the comparison rather than
    overriding it silently.
    """
    from ..io.registry import load_registry
    from ..structure.assembly import assemble_trimer, borrowed_fraction

    structure = structure or load_model()
    registry = {e.pdb: e for e in load_registry().available()}
    fits = []
    for pdb in candidates:
        record = registry.get(pdb)
        if record is None or record.n_protomers < 3:
            continue
        try:
            assembly = assemble_trimer(structure, template=pdb)
        except Exception:                                    # noqa: BLE001
            continue
        placement = assembly.placement_rmsd
        placement = (float(np.mean(placement))
                     if isinstance(placement, (tuple, list, np.ndarray))
                     else float(placement))
        try:
            borrowed = borrowed_fraction(assembly)["borrowed_fraction"]
        except Exception:                                    # noqa: BLE001
            borrowed = None
        fits.append(TemplateFit(
            template=pdb, protein=record.protein or "?",
            identity=float(assembly.alignment_identity),
            n_core=int(assembly.n_core),
            n_corresponding=int(assembly.n_corresponding),
            placement_rmsd=placement, clashes=int(assembly.clashes),
            borrowed_fraction=borrowed, ok=bool(assembly.ok)))
    return fits


def best_paralogue_template(fits: list[TemplateFit] | None = None) -> str | None:
    """The template to build piezo3 on, chosen on identity rather than size.

    piezo3's siblings are PIEZO1 and PIEZO2 — all three duplications sit on the
    jawed-vertebrate stem — and the worm PEZO-1 is an outgroup. Among templates
    that assemble, the highest sequence identity wins, with the clash count as
    the tie-break, because clashes are what an ill-fitting arrangement looks
    like when nothing models the inter-protomer contact.
    """
    fits = fits if fits is not None else template_survey()
    usable = [f for f in fits if f.ok]
    if not usable:
        return None
    return max(usable, key=lambda f: (round(f.identity, 3), -f.clashes)).template
