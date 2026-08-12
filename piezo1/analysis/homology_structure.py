"""The same comparison ``paralogue`` runs, over the whole family it can reach.

``analysis.paralogue`` asks whether this project's central mechanical claim —
that a symmetric low-frequency mode of the elastic network is the gating
coordinate — is a fact about PIEZO1 or about the fold. It answers with PIEZO2,
and it refuses everything else by construction: ``compare`` returns an error
unless it is handed one PIEZO1 entry and one PIEZO2 entry.

That was the right refusal when PIEZO2 was the only other structure in the
catalogue. It is not any more. Deposited coordinates now exist for **four**
PIEZOs — PIEZO1, PIEZO2, C. elegans PEZO-1 and Drosophila PIEZO — and the
question "is this the fold?" gets a far better answer from ~800 Myr of
divergence than from one vertebrate duplication.

**What this adds over ``paralogue``, and what it deliberately does not.**

*It removes the paralogue restriction* and keeps every guard that made the
PIEZO2 comparison believable: sites corresponded through a real alignment
rather than an offset, protomer order searched rather than read off chain
labels, both structures rotated into a common frame, and a shuffled-
correspondence control on every overlap.

*It renumbers first.* Two of the entries it can now reach are not in canonical
numbering — 9W7X is in a dPIEZO isoform's own coordinates, +3 after residue
1570 for 713 residues — and a cross-species alignment map applied to an entry
numbered in something else is wrong everywhere past the shift. ``paralogue``
never had to care, because no PIEZO2 entry has this problem.

*It refuses index-pairing of transmembrane helices across the family, and that
refusal is the most useful thing here.* ``paralogue.tm_index_correspondence``
pairs TM *k* of one protein with TM *k* of the other and confirms it against
the alignment for 37 of 38, which works because PIEZO1 and PIEZO2 both have 38.
Nothing else in the family does: **PEZO-1 has 36, dPIEZO 40, and the plant and
amoebal PIEZOs 35**. Index pairing is not merely unconfirmed there, it is
arithmetically impossible, and a function that returned "2 of 38 agree" would
be reporting a fact about counting as though it were a fact about structure.

*It reports a range rather than a comparison, because a comparison here is a
cherry-pick.* This was found the way these things usually are — by running a
second pair after the first had already been written down. 7WLT against 9W7X
gives a gating-mode overlap of **0.980** with dPIEZO; **8YEZ against the same
9W7X gives 0.189**, and does not clear its own shuffled control. So
``mode_overlap_spread`` runs every combination and returns a range with a
stability verdict. PIEZO2 is the positive control that keeps "not stable" from
being vacuous: six pairs, all 0.80-0.98, all clearing their controls.

*It refuses monomers, which costs the most interesting member.* The dome, the
pore and the elastic network all need three protomers. The only structural
representation of a non-animal PIEZO is an AlphaFold **monomer** of Arabidopsis
PIEZO — Dictyostelium pzoA has not even that — so the question of whether the
dome is a property of the fold rather than of animals **cannot be asked from
structure at all** with what exists. That is a gap in the world rather than in
this code, and it is reported as one.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..config import N_PROTOMERS
from ..core.numbering_check import identify_numbering, renumber
from ..core.structure import Structure
from .paralogue import dome_comparison, mode_comparison, paralogue_map

__all__ = ["StructuralPair", "HomologComparison", "comparable_entries",
           "helix_counts", "index_pairing_valid", "compare_structures",
           "OverlapSpread", "mode_overlap_spread", "MIN_PROTOMERS"]

#: Everything downstream needs a three-fold axis. A monomer is refused rather
#: than measured, which is what keeps an AlphaFold single chain from producing
#: a dome radius that reads like an experimental one.
#:
#: Taken from ``config.N_PROTOMERS`` rather than written as a 3 here: PIEZO
#: being a homotrimer is a fact about the protein that the whole project
#: already shares, and a second copy of it is a second thing to keep right.
MIN_PROTOMERS = N_PROTOMERS


def helix_counts() -> dict:
    """Annotated transmembrane helices per family member.

    The table that decides whether helices may be paired by index. Read from
    the committed UniProt resources, so it cannot disagree with the annotation
    the dome measurement selects its surface points from.
    """
    from .homology import family

    return {m.key: m.n_transmembrane for m in family()}


def index_pairing_valid(reference_a: str, reference_b: str) -> tuple[bool, str]:
    """May TM *k* of one be paired with TM *k* of the other?

    Only when both are annotated with the same number of helices. This is
    checked before the alignment is consulted, because a mismatch in count is
    already decisive: with 36 against 38 there is no index pairing to confirm
    or refute, and running the confirmation anyway produces a low agreement
    score that reads like a structural finding.
    """
    counts = helix_counts()
    first, second = counts.get(reference_a), counts.get(reference_b)
    if first is None or second is None:
        return False, "one of the two has no committed helix annotation"
    if first != second:
        return False, (f"{reference_a} is annotated with {first} transmembrane "
                       f"helices and {reference_b} with {second}, so pairing "
                       f"them by index is arithmetic rather than structure")
    return True, f"both annotated with {first} helices"


@dataclass
class StructuralPair:
    """One entry, prepared for comparison, or the reason it cannot be."""

    pdb: str
    structure: Structure | None
    reference: str = ""
    protein: str = ""
    n_protomers: int = 0
    renumbered: str = ""
    refusal: str = ""

    @property
    def usable(self) -> bool:
        return self.structure is not None and not self.refusal


def _prepare(pdb: str) -> StructuralPair:
    from ..io.registry import load_registry

    record = load_registry().get(pdb.upper())
    if record is None:
        return StructuralPair(pdb=pdb, structure=None,
                              refusal=f"{pdb} is not in the catalogue")
    if not record.available:
        return StructuralPair(pdb=pdb, structure=None,
                              refusal=f"{pdb} is not downloaded — run "
                                      f"python -m piezo1.io.fetch")

    structure = Structure.from_file(record.path)
    fixed, correction = renumber(structure)
    note = "" if not correction.needed else (
        f"renumbered to canonical: {correction.n_residues} residues, "
        f"identity {correction.identity_before:.3f} -> "
        f"{correction.identity_after:.3f}")
    identity = identify_numbering(fixed)

    protomers = sum(1 for chain in fixed.chains
                    if (fixed.mask_ca() & (fixed.chain == chain)).sum() > 300)
    pair = StructuralPair(
        pdb=pdb.upper(), structure=fixed, reference=identity.reference,
        protein=identity.protein, n_protomers=protomers, renumbered=note)

    if protomers < MIN_PROTOMERS:
        pair.refusal = (
            f"{pdb} models {protomers} protomer(s); the dome, the pore and the "
            f"elastic network all need {MIN_PROTOMERS}. This is what rules out "
            f"the plant PIEZO, whose only structural representation is an "
            f"AlphaFold monomer.")
    elif not identity.explained:
        pair.refusal = (f"{pdb} could not be read in any known numbering "
                        f"({identity.summary()})")
    return pair


def comparable_entries() -> dict:
    """Which catalogued entries can enter a structural comparison, by protein.

    **Best-resolved first within each protein**, measured as C-alphas in the
    largest protomer rather than by deposited resolution: a 3.1 A map of half
    the molecule compares worse than a 3.8 A map of all of it, and every
    quantity here is coverage-limited. Sorting these alphabetically would have
    made ``analysis_homology`` pick PEZO-1's 9UOX (1,324 residues) over 9UOY
    (1,919) for no reason but the letter X.

    Computed rather than listed: the answer changed when PEZO-1 and dPIEZO were
    catalogued and will change again, and a written list would be the thing
    that goes stale.
    """
    from ..io.registry import load_registry

    out: dict[str, list] = {}
    for record in load_registry():
        if not record.available or record.state == "predicted":
            continue
        if record.n_protomers < MIN_PROTOMERS:
            continue
        n_ca = max((c["n_ca"] for c in record.protomer_chains), default=0)
        out.setdefault(record.protein, []).append((-n_ca, record.pdb))
    return {protein: [pdb for _n, pdb in sorted(rows)]
            for protein, rows in sorted(out.items())}


@dataclass
class OverlapSpread:
    """The gating-mode overlap over **every** entry pair, not one of them.

    This exists because a single pair is not a measurement here, and finding
    that out was the most useful thing Round 89 did. Comparing 7WLT with 9W7X
    gives an overlap of 0.980 with dPIEZO; comparing 8YEZ with the same 9W7X
    gives **0.189**. Same two proteins, same method, and a number that would
    have gone into a document as "the gating coordinate is preserved in
    Drosophila".

    Reported as a range with the fraction of pairs that clear their own
    shuffled control, so the instability is the headline rather than a
    footnote under whichever pair was run first.
    """

    protein: str
    n_pairs: int
    overlaps: tuple
    controls: tuple
    pairs: tuple
    n_beating_control: int

    @property
    def low(self) -> float:
        return min(self.overlaps) if self.overlaps else float("nan")

    @property
    def high(self) -> float:
        return max(self.overlaps) if self.overlaps else float("nan")

    @property
    def stable(self) -> bool:
        """Every pair clears its control **and** they agree with each other.

        Both clauses are needed. A set of pairs that all beat their controls
        while ranging from 0.18 to 0.98 is not evidence that the mode is
        shared; it is evidence that which mode wins depends on what each entry
        happens to resolve.
        """
        return (self.n_pairs > 0
                and self.n_beating_control == self.n_pairs
                and (self.high - self.low) < 0.3)

    def summary(self) -> str:
        verdict = ("stable" if self.stable else
                   "NOT STABLE, so no single number here is a property of the "
                   "two proteins")
        return (f"{self.protein}: overlap {self.low:.2f}-{self.high:.2f} over "
                f"{self.n_pairs} entry pairs, {self.n_beating_control} of them "
                f"beating their shuffled control — {verdict}")


def mode_overlap_spread(piezo1_entries, partner_entries, protein: str,
                        n_modes: int | None = None) -> OverlapSpread:
    """Run every combination and report the range rather than a representative."""
    overlaps, controls, pairs = [], [], []
    beating = 0
    for a in piezo1_entries:
        for b in partner_entries:
            comparison = compare_structures(a, b, n_modes=n_modes,
                                            with_dome=False)
            if not comparison.ok:
                continue
            modes = comparison.modes
            overlaps.append(float(modes.best_overlap))
            controls.append(float(modes.shuffled_control))
            pairs.append((a, b))
            beating += bool(modes.beats_control)
    return OverlapSpread(protein=protein, n_pairs=len(pairs),
                         overlaps=tuple(overlaps), controls=tuple(controls),
                         pairs=tuple(pairs), n_beating_control=beating)


@dataclass
class HomologComparison:
    """Two homologues, compared as far as the structures allow."""

    a: StructuralPair
    b: StructuralPair
    dome: dict | None = None
    modes: object | None = None
    tm_pairing: dict = field(default_factory=dict)
    alignment_identity: float = float("nan")
    error: str = ""

    @property
    def ok(self) -> bool:
        return not self.error

    def summary(self) -> str:
        if self.error:
            return f"refused: {self.error}"
        lines = [f"{self.a.protein} {self.a.pdb} vs {self.b.protein} "
                 f"{self.b.pdb}; sequence identity "
                 f"{self.alignment_identity:.3f}"]
        for pair in (self.a, self.b):
            if pair.renumbered:
                lines.append(f"  {pair.pdb}: {pair.renumbered}")
        if not self.tm_pairing.get("valid"):
            lines.append(f"  transmembrane helices NOT paired by index — "
                         f"{self.tm_pairing.get('reason', '')}")
        elif "n_agree" in self.tm_pairing:
            lines.append(f"  TM index pairing confirmed by the alignment for "
                         f"{self.tm_pairing['n_agree']} of "
                         f"{self.tm_pairing['n_helices']}")
        if self.modes is not None:
            lines.append(f"  {self.modes.summary()}")
        if self.dome is not None:
            for row in self.dome["coverage_matched"]:
                lines.append(f"  {row.summary()}")
        return "\n".join(lines)


def compare_structures(pdb_a: str, pdb_b: str, n_modes: int | None = None,
                       with_dome: bool = True) -> HomologComparison:
    """Compare any two catalogued homologues.

    No restriction on which proteins they are — the point of the module — but
    every refusal ``paralogue`` made for a reason is kept, and the two it never
    needed (isoform numbering, monomers) are added.
    """
    first, second = _prepare(pdb_a), _prepare(pdb_b)
    if not first.usable or not second.usable:
        return HomologComparison(
            a=first, b=second,
            error="; ".join(p.refusal for p in (first, second) if p.refusal))
    if first.reference == second.reference:
        return HomologComparison(
            a=first, b=second,
            error=(f"both entries are {first.reference}; this compares "
                   f"homologues, and two entries of the same protein are a "
                   f"state comparison — use the overlay or the morph"))

    valid, reason = index_pairing_valid(first.reference, second.reference)
    pairing = {"valid": valid, "reason": reason}
    if valid:
        numbering = paralogue_map(first.reference, second.reference)
        from .paralogue import tm_index_correspondence

        pairing.update({k: v for k, v in tm_index_correspondence(
            numbering, first.reference, second.reference).items()
            if k != "rows"})

    numbering = paralogue_map(first.reference, second.reference)
    comparison = HomologComparison(
        a=first, b=second, tm_pairing=pairing,
        alignment_identity=numbering.identity)
    try:
        comparison.modes = mode_comparison(first.structure, second.structure,
                                           n_modes=n_modes)
    except ValueError as exc:
        comparison.error = f"modes: {exc}"
        return comparison
    if with_dome:
        try:
            comparison.dome = dome_comparison(
                first.structure, second.structure, first.pdb, second.pdb)
        except (ValueError, KeyError) as exc:
            # A dome needs annotated transmembrane helices resolved in both. It
            # is reported as missing rather than dropped, because a comparison
            # silently missing half its content reads as a completed one.
            comparison.dome = None
            comparison.tm_pairing["dome_error"] = str(exc)
    return comparison
