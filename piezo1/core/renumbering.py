"""Applying what the numbering check detects — the rewriter, not the detector.

Split from :mod:`piezo1.core.numbering_check` at the 500-line limit and along
the seam the test layout already recognised (``tests/test_renumbering.py`` has
existed since Round 86). That module asks *what numbering is this file in*; this
one rewrites a file into canonical numbering, and the two are different kinds of
thing: the first is a measurement and the second is a mutation.

**The null is the load-bearing case.** A renumberer that touched a correct file
would corrupt it, so an entry needing no correction is returned **unchanged, by
identity** — 8YEZ resolves the same 767-857 region the four register-error
entries get wrong and must come back untouched, as must 7WLT, 3JAC and 6B3R. A
caller cannot tell a corrected file from an uncorrected one without reading the
report, which is the point.

What it repairs, measured: 6LQI to 1.000 with +24 over 765 residues; 8ZU3, 8YFC
and 9VMX to 1.000; 9W7X — Drosophila PIEZO, an isoform numbering found by the
detector rather than by reading the paper — to 1.000 with +3 over 713; and 8YFG
to **0.999 and not 1.000**, which is correct, because its R2456H substitution is
a real residue change a numbering fix must not absorb.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field

import numpy as np

from .numbering_check import (REFERENCES, identify_numbering, mismatch_blocks,
                              reference_entry, _named)
from .structure import AA3TO1

__all__ = ["Renumbering", "canonical_renumbering", "apply_renumbering",
           "renumber"]


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
