"""What two paralogues kept of each other, domain by domain, recomputed here.

Split from ``family_constraint.py`` at its length limit and along a real seam.
That module measures constraint *within* one paralogue — how much 174 versions
of PIEZO1 have been willing to change each residue. This one measures identity
*between* paralogues, which is a different question with a different failure
mode: the within-paralogue comparison is always in one numbering, and this one
never is.

The census's headline is that PIEZO1, PIEZO2 and piezo3 agree at 45-48% of
positions overall and at 73-91% across the anchor, the pore helices and the CTD
— with the **CED**, the lid over the pore, the one part of the pore machinery
*below* the whole-protein figure. Recomputed here from this project's own global
alignment over its own domain boundaries, the cap lands at 0.404 between PIEZO1
and PIEZO2 against the census's 0.402.

**Only PIEZO1-framed pairs are answered**, and that refusal is the reason this
module is worth reading. The boundaries are ``domains.json``'s, which are human
PIEZO1 residue numbers, and this project holds annotation for no other PIEZO.
Framing a comparison in PIEZO2 and indexing it with PIEZO1's ranges put the cap
at **0.85 identity where the census measures 0.35**, and the ordering of every
other domain with it — and announced nothing.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..core.annotations import Annotations, load_annotations

__all__ = ["DomainIdentity", "paralogue_asymmetry", "PAIRS"]


@dataclass(frozen=True)
class DomainIdentity:
    """How alike two paralogues are over one domain, on our boundaries."""

    pair: str
    domain: str
    category: str
    n_columns: int
    identity: float
    whole_protein: float

    @property
    def vs_whole(self) -> float:
        return self.identity - self.whole_protein

    @property
    def above_average(self) -> bool:
        return self.identity > self.whole_protein


def paralogue_asymmetry(pair: str = "PIEZO1_vs_PIEZO2",
                        annotations: Annotations | None = None
                        ) -> list[DomainIdentity]:
    """What two paralogues kept of each other, domain by domain, recomputed here.

    The census's headline is that PIEZO1, PIEZO2 and piezo3 agree at 45-48% of
    positions overall and at 73-91% across the anchor, the pore helices and the
    CTD — with the CED, the lid over the pore, the one part of the pore
    machinery *below* the whole-protein figure.

    Recomputed rather than read back, and the difference is real: the identity
    is taken over **our** domain boundaries from **our** global alignment
    (:func:`piezo1.core.sequence.align_global`, BLOSUM62 with the registered
    gap costs), not over the census's bands from its 117-sequence family
    alignment. The two partitions put the anchor 141 residues apart, so
    agreement here is a property of the proteins and not of a boundary choice.

    Identity is counted over **mutually aligned columns only** — a position
    where one sequence has a gap is not a mismatch, it is a position where the
    question does not arise. The census measured that this convention moves
    pairs by up to 51 identity points, which is why it is stated rather than
    left to an aligner default.

    **Only PIEZO1-framed pairs are answered.** The domain boundaries are
    ``domains.json``'s, which are human PIEZO1 residue numbers, and this
    project has annotation for no other PIEZO. Framing a comparison in PIEZO2
    and indexing it with PIEZO1's ranges is the same error
    :data:`piezo1.core.annotations.ANNOTATED_NUMBERINGS` exists to prevent, and
    it does not announce itself: on ``PIEZO2_vs_piezo3`` it put the cap at
    0.85 identity where the census measures 0.35, and the ordering of every
    other domain with it. A pair framed elsewhere raises, and its census row is
    still readable through ``load_family_findings().table('paralogue_identity')``.
    """
    from ..core.sequence import align_global

    left_name, right_name = _pair_references(pair)
    if left_name != "human":
        raise KeyError(
            f"{pair} is framed in {left_name}, and this project holds domain "
            f"boundaries only for PIEZO1. Read the census's own row for it "
            f"instead: load_family_findings().table('paralogue_identity')")
    left, right = _reference_sequence(left_name), _reference_sequence(right_name)
    a, b = align_global(left, right)

    # Walk the alignment once, carrying each column's residue number in the
    # left sequence, so a domain can be read off without re-aligning.
    columns, index = [], 0
    for x, y in zip(a, b):
        if x != "-":
            index += 1
            columns.append((index, x, y))
    matched = [(i, x, y) for i, x, y in columns if y != "-"]
    whole = (sum(1 for _, x, y in matched if x == y) / len(matched)
             if matched else 0.0)

    ann = annotations or load_annotations("human")
    out = []
    for dom in ann.domains:
        if dom.sub_element or dom.start is None or dom.end is None:
            continue
        inside = [(x, y) for i, x, y in matched if dom.start <= i <= dom.end]
        if len(inside) < 10:
            continue
        identity = sum(1 for x, y in inside if x == y) / len(inside)
        out.append(DomainIdentity(pair=pair, domain=dom.id, category=dom.category,
                                  n_columns=len(inside), identity=identity,
                                  whole_protein=whole))
    return out


#: Which committed reference each half of a pair name refers to. The pair names
#: are the census's, so its table and ours can be put side by side by key.
PAIRS = {
    "PIEZO1_vs_PIEZO2": ("human", "human_piezo2"),
    "PIEZO1_vs_piezo3": ("human", "piezo3"),
    "PIEZO2_vs_piezo3": ("human_piezo2", "piezo3"),
}


def _pair_references(pair: str) -> tuple:
    if pair not in PAIRS:
        raise KeyError(f"unknown pair {pair!r}; have {', '.join(PAIRS)}")
    return PAIRS[pair]


def _reference_sequence(name: str) -> str:
    """One reference sequence, whichever resource holds it.

    piezo3 is the odd one out: it is a family reference like the rest, and its
    sequence also arrives with the imported constraint track. The UniProt
    resource is preferred because that is the record the model is keyed on —
    the two differ by an inserted residue, and picking the wrong one would
    shift every identity past position 2014.
    """
    import json

    from ..config import RESOURCE_DIR
    filename = ("uniprot_zebrafish_piezo3.json" if name == "piezo3"
                else f"uniprot_{name}.json")
    return json.loads((RESOURCE_DIR / filename).read_text())["sequence"]
