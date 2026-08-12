"""What protein a deposited file is, and whether we can read it by residue number.

It lives in ``core`` beside :mod:`piezo1.core.entities`, which answers the
neighbouring question of *what is in* a deposited file, and because
``structure`` needs it: the full-length graft has to pick the AlphaFold model
matching the entry it is grafting onto, and ``structure`` importing
``analysis`` would point the dependency arrow backwards.

Round 83 built this to keep a PIEZO1 structure from being read with PIEZO2's
transmembrane annotation. It immediately found two things nobody was looking
for, and both are live — this project reads domains, helices, variants and
functional residues by **residue number**, so a file numbered differently from
the reference is silently annotated in the wrong place.

The measurement is simple and has a known answer. A deposited file states each
residue's number *and* its name; if the numbering belongs to a sequence, the
names agree with it at every position. Across the catalogue the right reference
scores 1.000 and every other one below 0.25, so a score in between is a fact
about the file rather than about the method.

**What it found.**

* **6LQI** is deposited in the Piezo1.1 isoform's own continuous numbering
  across its 24-residue deletion — 1.000 before the splice site, 0.058 after,
  1.000 again shifted by +24, for 764 of 1,301 resolved residues.
* **8ZU3, 8YFC, 9VMX and 8YFG** carry residues 767-857 numbered 22 low. That is
  91 residues inside a whole-file identity of 0.932, which passes any sensible
  floor — which is why runs of disagreement are reported and not just totals.

**And one thing it got wrong first.** 3JAC scored 0.623 and looked like a third
case. Every one of its mismatches was a ``UNK`` — the depositor declining to
name a residue, not disagreeing about one. ``AA3TO1`` maps UNK to X, so
membership in it is not the test. Excluding unassigned residues, 3JAC matches
at 1.000 over the 572 it does name.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field

import numpy as np

from ..config import RESOURCE_DIR
from ..parameters import PARAMETERS as _P
from .structure import AA3TO1, Structure

__all__ = ["NumberingIdentity", "SpliceShift", "MismatchBlock",
           "identify_numbering", "piezo1_numbering", "detect_splice", "mismatch_blocks",
           "Renumbering", "canonical_renumbering", "apply_renumbering", "renumber",
           "reference_entry", "REFERENCES", "PIEZO1_REFERENCES",
           "PIEZO2_REFERENCES", "INVERTEBRATE_REFERENCES"]

#: The committed UniProt resources an entry can be scored against.
REFERENCES = ("human", "mouse", "human_piezo2", "mouse_piezo2",
              "worm_piezo", "fly_piezo")
PIEZO1_REFERENCES = ("human", "mouse")
PIEZO2_REFERENCES = ("human_piezo2", "mouse_piezo2")
#: The invertebrate PIEZOs. Neither is a PIEZO1 or a PIEZO2 — the duplication
#: that produced those two is vertebrate — so they are their own category and
#: the generality question they answer is a wider one.
INVERTEBRATE_REFERENCES = ("worm_piezo", "fly_piezo")


def reference_entry(name: str) -> dict:
    """The committed UniProt resource for one reference, as a dict."""
    import json

    return json.loads((RESOURCE_DIR / f"uniprot_{name}.json").read_text())


def _named(residue) -> bool:
    """Whether the file names this residue at all.

    ``AA3TO1`` maps ``UNK`` to ``X``, so membership is not the test — a
    poly-UNK trace would score as 24 mismatching residues rather than as a
    depositor declining to assign them. 3JAC carries 346 unknowns in 918
    C-alphas and matches at 1.000 on the rest; counting them made it 0.623 and
    looked exactly like a numbering error.
    """
    return AA3TO1.get(str(residue), "X") != "X"


@dataclass(frozen=True)
class MismatchBlock:
    """A contiguous run of residues whose names disagree with the sequence.

    The whole-file identity hides these: 8ZU3 scores 0.932, comfortably past
    any sensible floor, and the 7% it is missing is not spread out. It is a
    single block of 91 residues, 767-857, every one of which disagrees — and
    every one of which agrees again if the numbers are read 22 higher. Four of
    the six human entries carry it.

    A localised register error like that is invisible in a total and matters
    exactly as much as a whole-file one for anything that reads annotation by
    residue number inside the block.
    """

    start: int
    end: int
    n_residues: int
    repaired_by: int | None      # shift that restores agreement, if one does
    repaired_identity: float

    def summary(self) -> str:
        fix = ("no single shift repairs it" if self.repaired_by is None
               else f"reading them {self.repaired_by:+d} restores "
                    f"{self.repaired_identity:.3f}")
        return (f"{self.n_residues} residues, {self.start}-{self.end}, "
                f"disagree; {fix}")


@dataclass(frozen=True)
class SpliceShift:
    """A file numbered in a splice isoform's own coordinates, not the canonical.

    Found by the identification refusing 6LQI, the Piezo1.1 isoform. Its
    residue names agree with canonical mouse Piezo1 **exactly** up to the
    splice site and then not at all — because the deposited numbering runs
    continuously across a 24-residue deletion, so every residue after it is 24
    lower than its canonical counterpart.

    That is 764 of 6LQI's 1,301 resolved residues, more than half the chain.
    Any annotation applied to this entry by canonical residue number — a
    transmembrane helix, a domain boundary, a variant — is wrong there.
    """

    breakpoint: int          # last residue that matches without a shift
    offset: int              # add this to a deposited number to get canonical
    identity_before: float
    identity_after: float
    n_after: int

    def summary(self) -> str:
        return (f"numbered in an isoform's own coordinates: canonical up to "
                f"{self.breakpoint}, then {self.offset:+d} for {self.n_after} "
                f"residues (agreement {self.identity_before:.3f} then "
                f"{self.identity_after:.3f} once shifted)")


@dataclass
class NumberingIdentity:
    """Which protein and numbering a deposited file is actually in.

    Read off the coordinates rather than off a label, because the label is
    what would be wrong. The margin is carried because a confident answer and
    a coin flip look identical once only the winner is reported, and
    ``n_unassigned`` because a file that declines to name most of its residues
    can only be identified from the rest.
    """

    reference: str
    identity: float
    n_scored: int
    runner_up: str
    runner_up_identity: float
    scores: dict = field(default_factory=dict)
    #: Set when the file is numbered in a splice isoform's own coordinates.
    splice: SpliceShift | None = None
    #: Residues the file models without naming — ``UNK``. Not scored, because
    #: an unassigned residue is the depositor saying they do not know, not a
    #: disagreement. 3JAC carries 346 of them in 918 C-alphas and matches at
    #: 1.000 on the rest; scoring them as mismatches put it at 0.623 and would
    #: have been read as a numbering error.
    n_unassigned: int = 0
    #: Contiguous runs of disagreement, which a whole-file identity hides.
    blocks: tuple = ()

    @property
    def margin(self) -> float:
        return self.identity - self.runner_up_identity

    @property
    def is_piezo2(self) -> bool:
        return self.reference in PIEZO2_REFERENCES

    @property
    def is_invertebrate(self) -> bool:
        return self.reference in INVERTEBRATE_REFERENCES

    @property
    def protein(self) -> str:
        """PIEZO1, PIEZO2, PEZO-1 or dPIEZO — what to call this on screen."""
        return {"human": "PIEZO1", "mouse": "PIEZO1",
                "human_piezo2": "PIEZO2", "mouse_piezo2": "PIEZO2",
                "worm_piezo": "PEZO-1", "fly_piezo": "dPIEZO",
                }.get(self.reference, "unknown")

    @property
    def confident(self) -> bool:
        return self.identity >= _P.value("paralogue.min_sequence_identity") \
            and self.margin >= _P.value("paralogue.min_identity_margin")

    @property
    def explained(self) -> bool:
        """Confident, or accounted for by a splice shift that is reported."""
        return self.confident or self.splice is not None

    @property
    def clean(self) -> bool:
        """Confident *and* free of localised register errors."""
        return self.confident and not self.blocks

    def summary(self) -> str:
        unassigned = ("" if not self.n_unassigned
                      else f" ({self.n_unassigned} unassigned, not scored)")
        base = (f"{self.reference} at {self.identity:.3f} over "
                f"{self.n_scored} residues{unassigned}, next best "
                f"{self.runner_up} at {self.runner_up_identity:.3f}")
        if self.splice is not None:
            base = f"{base}; {self.splice.summary()}"
        for block in self.blocks:
            base = f"{base}; {block.summary()}"
        return base


def identify_numbering(structure: Structure) -> NumberingIdentity:
    """Score a structure's own residue names against every reference sequence.

    A deposited file states its residue numbers and its residue names; if the
    numbering belongs to a given sequence, the names agree with it at every
    position. That makes this a measurement with a known answer — an entry
    should match exactly one reference near 1.0 and everything else near the
    background rate for shuffled protein — rather than a guess.
    """
    mask = structure.mask_ca()
    if structure.chains:
        mask = mask & (structure.chain == structure.chains[0])
    numbers = structure.res_seq[mask]
    names = structure.res_name[mask]

    assigned = np.array([_named(n) for n in names])
    numbers, names = numbers[assigned], names[assigned]

    scores, counts = {}, {}
    for name in REFERENCES:
        sequence = reference_entry(name)["sequence"]
        matched = scored = 0
        for number, residue in zip(numbers, names):
            index = int(number) - 1
            if 0 <= index < len(sequence):
                scored += 1
                matched += AA3TO1[str(residue)] == sequence[index]
        scores[name] = matched / max(scored, 1)
        counts[name] = scored

    ranked = sorted(scores, key=lambda k: scores[k], reverse=True)
    identity = NumberingIdentity(
        reference=ranked[0], identity=scores[ranked[0]],
        n_scored=counts[ranked[0]], runner_up=ranked[1],
        runner_up_identity=scores[ranked[1]], scores=scores,
        n_unassigned=int((~assigned).sum()))
    sequence = reference_entry(ranked[0])["sequence"]
    if not identity.confident:
        shift = detect_splice(numbers, names, sequence)
        if shift is not None:
            identity.splice = shift
    if identity.splice is None:
        identity.blocks = tuple(mismatch_blocks(numbers, names, sequence))
    return identity


def mismatch_blocks(numbers, names, sequence: str) -> list:
    """Contiguous runs of disagreement, and the shift that would repair each.

    Run whenever the file is not a splice case, because a total cannot show
    where its errors are. A run has to be long enough not to be a point
    mutation — an engineered variant is a real residue change and not a
    numbering fault — which is what the registered minimum length is for.
    """
    minimum = int(_P.value("paralogue.min_mismatch_block"))
    span = int(_P.value("paralogue.max_splice_shift"))
    pairs = sorted(((int(n), str(r)) for n, r in zip(numbers, names)
                    if _named(r)), key=lambda p: p[0])

    blocks, run = [], []
    for number, residue in pairs:
        index = number - 1
        agrees = 0 <= index < len(sequence) and AA3TO1[residue] == sequence[index]
        if agrees:
            if len(run) >= minimum:
                blocks.append(run)
            run = []
        else:
            run.append((number, residue))
    if len(run) >= minimum:
        blocks.append(run)

    out = []
    for run in blocks:
        best, best_shift = 0.0, None
        for shift in range(-span, span + 1):
            if shift == 0:
                continue
            matched = sum(1 for number, residue in run
                          if 0 <= number + shift - 1 < len(sequence)
                          and AA3TO1[residue] == sequence[number + shift - 1])
            if matched / len(run) > best:
                best, best_shift = matched / len(run), shift
        out.append(MismatchBlock(
            start=run[0][0], end=run[-1][0], n_residues=len(run),
            repaired_by=best_shift if best >= 0.9 else None,
            repaired_identity=best))
    return out


def detect_splice(numbers, names, sequence: str) -> SpliceShift | None:
    """Is a file numbered in a splice isoform's own coordinates?

    Only asked when the straight identification fails, and it answers a
    specific shape of failure: perfect agreement up to some residue and none
    after it, restored by one constant shift. Two constants rather than one, so
    it is still not an offset anyone may hard-code — the breakpoint and the
    shift are both properties of the isoform.

    Returns ``None`` unless *both* halves are essentially perfect, because a
    partial rescue would let any badly numbered file be explained away.
    """
    floor = _P.value("paralogue.min_sequence_identity")

    def agreement(pairs, shift):
        matched = scored = 0
        for number, residue in pairs:
            index = int(number) + shift - 1
            if 0 <= index < len(sequence):
                scored += 1
                matched += AA3TO1[str(residue)] == sequence[index]
        return (matched / scored if scored else 0.0), scored

    pairs = sorted(((n, r) for n, r in zip(numbers, names) if _named(r)),
                   key=lambda p: int(p[0]))
    matches = [AA3TO1[str(r)] == sequence[int(n) - 1]
               if 0 <= int(n) - 1 < len(sequence) else False for n, r in pairs]
    if not matches or not matches[0]:
        return None
    cut = next((i for i, ok in enumerate(matches) if not ok), len(matches))
    if cut == len(matches) or cut < 2:
        return None

    before, _ = agreement(pairs[:cut], 0)
    span = int(_P.value("paralogue.max_splice_shift"))
    best, best_shift, n_after = 0.0, 0, 0
    for shift in range(-span, span + 1):
        if shift == 0:
            continue
        score, scored = agreement(pairs[cut:], shift)
        if score > best:
            best, best_shift, n_after = score, shift, scored
    if before < floor or best < floor:
        return None
    return SpliceShift(breakpoint=int(pairs[cut - 1][0]), offset=best_shift,
                       identity_before=before, identity_after=best,
                       n_after=n_after)


def piezo1_numbering(structure) -> str | None:
    """Which numbering a PIEZO1 entry is in, or ``None`` if it is not one.

    ``"human"`` or ``"mouse"`` when :func:`identify_numbering` matches a PIEZO1
    reference and finds no splice shift; ``None`` otherwise — a paralogue, an
    invertebrate PIEZO, or 6LQI, which is PIEZO1 and still cannot be read by
    residue number because it is deposited in the Piezo1.1 isoform's own
    numbering.

    A ``None`` is a refusal and never a default to human. Anything that reads
    annotation by residue number has to ask this first, which is why it lives
    in ``core`` rather than beside its first caller: ``physics`` needs it to
    bound a conduction path and ``analysis`` needs it to locate a gate, and
    the arrow between those two only points one way.
    """
    try:
        identity = identify_numbering(structure)
    except Exception:                            # noqa: BLE001
        return None
    if identity.reference not in PIEZO1_REFERENCES or identity.splice is not None:
        return None
    return identity.reference


@dataclass(frozen=True)
class Renumbering:
    """A correction from a file's own numbering to its reference's.

    Detection already existed — :func:`detect_splice` and
    :func:`mismatch_blocks` both report the shift that would repair what they
    find. Nothing applied it, so five entries were read at residue numbers that
    point at the wrong residue: every transmembrane helix, domain boundary and
    variant inside the affected range.
    """

    reference: str
    #: ``(first, last, shift)`` in the file's own numbering. Add ``shift`` to a
    #: deposited number in that range to get the canonical one.
    shifts: tuple = ()
    n_corrected: int = 0          # ATOMS whose number changed
    n_residues: int = 0           # distinct residues, which is the readable count
    identity_before: float = float("nan")
    identity_after: float = float("nan")
    reason: str = ""

    @property
    def needed(self) -> bool:
        return bool(self.shifts)

    def summary(self) -> str:
        if not self.needed:
            return (f"{self.reference} numbering, no correction needed "
                    f"({self.identity_before:.3f} as deposited)")
        parts = ", ".join(f"{lo}-{hi} read {sh:+d}" for lo, hi, sh in self.shifts)
        return (f"{self.reference} numbering with {self.n_residues} residues "
                f"({self.n_corrected} atoms) corrected ({parts}); identity "
                f"{self.identity_before:.3f} -> {self.identity_after:.3f}. "
                f"The spans are in the FILE's numbering and may cover "
                f"unresolved gaps; only resolved residues move.")


def canonical_renumbering(structure) -> Renumbering:
    """What correction, if any, this entry's residue numbers need.

    Returns an empty correction for a file that is already right — which is the
    case that makes it a measurement rather than a rewriter. 8YEZ resolves the
    same 767-857 region as the four entries that carry the register error and
    needs nothing; if this returned a shift for it, it would be inventing one.
    """
    identity = identify_numbering(structure)
    if identity.reference not in REFERENCES:
        return Renumbering(reference="", reason="no reference matches this file")

    sequence = reference_entry(identity.reference)["sequence"]
    shifts = []
    if identity.splice is not None:
        shifts.append((identity.splice.breakpoint + 1, 10 ** 9,
                       int(identity.splice.offset)))
    else:
        for block in identity.blocks:
            if block.repaired_by is None:
                continue
            shifts.append((int(block.start), int(block.end),
                           int(block.repaired_by)))

    if not shifts:
        return Renumbering(reference=identity.reference,
                           identity_before=identity.identity,
                           identity_after=identity.identity)

    shifts = _consolidate(structure, sequence, shifts)
    corrected = apply_renumbering(structure.res_seq, shifts)
    after = _identity_of(structure, corrected, sequence)
    changed = corrected != np.asarray(structure.res_seq)
    return Renumbering(
        reference=identity.reference, shifts=tuple(shifts),
        n_corrected=int(changed.sum()),
        n_residues=int(len(set(int(v) for v in
                               np.asarray(structure.res_seq)[changed]))),
        identity_before=identity.identity, identity_after=after)


def _consolidate(structure, sequence: str, shifts) -> tuple:
    """Merge and extend the detected spans, keeping only what helps.

    `mismatch_blocks` finds *runs* of disagreement, and a residue that agrees
    by chance at the wrong numbering ends a run. On 8ZU3 that split one
    91-residue register error into three blocks and left 772-787, 789-834 and
    839-857 corrected with the gaps between them untouched — identity 0.969
    where a uniform read of 767-857 gives 1.000.

    So spans sharing a shift are merged across small gaps and then grown
    outward one residue at a time, and **every step is kept only if the
    corrected identity does not fall**. That makes the extension a measurement
    rather than a guess: on an entry that needs no correction there is nothing
    to extend, and on one that does, the boundary lands where the agreement
    stops improving.
    """
    numbers = np.asarray(structure.res_seq).astype(int)
    if not len(numbers):
        return tuple(shifts)
    low_limit, high_limit = int(numbers.min()), int(numbers.max())

    by_shift: dict[int, list] = {}
    for low, high, shift in shifts:
        by_shift.setdefault(int(shift), []).append(
            [int(low), min(int(high), high_limit)])

    out = []
    for shift, spans in by_shift.items():
        spans.sort()
        merged = [spans[0]]
        for low, high in spans[1:]:
            if low - merged[-1][1] <= _MERGE_GAP:
                merged[-1][1] = max(merged[-1][1], high)
            else:
                merged.append([low, high])
        for span in merged:
            out.append([span[0], span[1], shift])

    score = lambda proposal: _identity_of(                       # noqa: E731
        structure, apply_renumbering(numbers, proposal), sequence)
    best = score(out)
    for i, (low, high, shift) in enumerate(out):
        for step, index in ((-1, 0), (1, 1)):
            while True:
                trial = [list(x) for x in out]
                edge = trial[i][index] + step
                if not low_limit <= edge <= high_limit:
                    break
                trial[i][index] = edge
                value = score(trial)
                if value < best:
                    break
                best, out = value, trial
    return tuple(tuple(x) for x in out)


#: Residues of chance agreement tolerated inside one register error before it
#: is treated as two. Not a physical quantity — a run-joining tolerance.
_MERGE_GAP = 10


def apply_renumbering(res_seq, shifts) -> "np.ndarray":
    """Residue numbers with each shift applied over its own range."""
    out = np.asarray(res_seq).astype(int).copy()
    for low, high, shift in shifts:
        inside = (out >= int(low)) & (out <= int(high))
        out[inside] += int(shift)
    return out


def renumber(structure):
    """``(corrected structure, Renumbering)``.

    The structure is returned **unchanged** when no correction is needed, so a
    caller cannot tell a corrected file from an uncorrected one by identity
    alone — it has to read the report, which is the point.
    """
    correction = canonical_renumbering(structure)
    if not correction.needed:
        return structure, correction
    fixed = dataclasses.replace(
        structure, res_seq=apply_renumbering(structure.res_seq,
                                             correction.shifts),
        name=f"{structure.name}(renumbered)")
    return fixed, correction


def _identity_of(structure, numbers, sequence: str) -> float:
    mask = structure.mask_ca()
    if structure.chains:
        mask = mask & (structure.chain == structure.chains[0])
    hit = total = 0
    for number, residue in zip(numbers[mask], structure.res_name[mask]):
        if not _named(residue):
            continue
        index = int(number) - 1
        if 0 <= index < len(sequence):
            total += 1
            hit += AA3TO1[str(residue)] == sequence[index]
    return hit / total if total else float("nan")
