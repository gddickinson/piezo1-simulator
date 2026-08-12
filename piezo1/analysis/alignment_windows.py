"""Is an alignment locally trustworthy here? — the window instrument.

Split out of ``homology_sites`` at the 500-line limit, along a real seam: this
module answers *is the alignment in register at this position*, which is a
question about a pair of sequences and nothing to do with PIEZO. The module
next door asks what the curated PIEZO sites do across the family, and it is
this instrument that decides which of its answers may be believed.

**The problem it exists for.** A global alignment always produces an answer. Ask
it where human PIEZO1's transmembrane gate sits in Arabidopsis PIEZO and it will
name three residues, whether or not that stretch of the two proteins has
anything to do with each other. Report the letters found there and you have
manufactured a conservation result out of an alignment artefact — a failure this
project has been caught by in other guises and guards against by calibration
everywhere else.

**The statistic is the BLOSUM62 window score, not the window identity.** Identity
is the weak statistic below Rost's line — that is the whole argument of
``homology`` — so judging an alignment's reliability with it would be
self-contradicting, and it is also underpowered: at the family's real divergence
a window of identity cannot separate a correct mapping from a chance one at any
usable width.

**The null is measured, not assumed.** For each pair, the same window statistic
is computed across an alignment of the source against a **composition-matched
shuffle** of the target. One shuffled alignment yields thousands of overlapping
windows, which is ample for the mean and spread, and it gives them for that
specific pair rather than from a table — human against mouse and human against
pzoA have different backgrounds and a single fixed threshold would be wrong for
both.

**The width was chosen by measuring power.** See ``homology.site_window``: at
width 31 the test refuses mappings that are visibly right, and 101 is where the
worm, fly and plant clear the threshold while Dictyostelium — correctly — does
not. The cost is stated rather than hidden: at 101 columns this answers *is this
region in register*, not *is this residue*.

Calibrated in ``tests/test_homology.py`` on three cases whose answers are known
by construction, including the one that makes it a measurement rather than a
formality: a single conserved block planted inside otherwise scrambled sequence,
which the window must find and the whole-sequence identity must miss.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.sequence import align_global
from ..parameters import PARAMETERS as _P

__all__ = ["AlignmentWindows", "alignment_windows", "window_identity",
           "window_score", "column_scores", "window_null_distribution"]


_MATRIX: dict | None = None


def _blosum() -> dict:
    global _MATRIX
    if _MATRIX is None:
        from Bio.Align import substitution_matrices

        matrix = substitution_matrices.load("BLOSUM62")
        _MATRIX = {(a, b): float(matrix[a, b])
                   for a in matrix.alphabet for b in matrix.alphabet}
    return _MATRIX


def column_scores(top: str, bottom: str):
    """Per-column BLOSUM62 score and a mask of the columns it is defined on.

    Gapped columns score zero and are excluded from the mask, so a window
    landing in an indel reports the score of the residues it does align rather
    than being punished for the gap. A window with no aligned column at all
    reports the null, which makes it unreliable — correctly.
    """
    matrix = _blosum()
    scores = np.zeros(len(top))
    scored = np.zeros(len(top))
    for k, (a, b) in enumerate(zip(top, bottom)):
        if a == "-" or b == "-":
            continue
        scored[k] = 1.0
        scores[k] = matrix.get((a, b), 0.0)
    return scores, scored


def _window_mean(scores, scored, column: int, window: int) -> float:
    half = window // 2
    lo, hi = max(0, column - half), min(len(scores), column + half + 1)
    n = float(scored[lo:hi].sum())
    return float(scores[lo:hi].sum()) / n if n else 0.0


def window_null_distribution(source_sequence: str, target_sequence: str,
                              seed: int = 0) -> tuple[float, float]:
    """Mean and spread of the window statistic when there is no homology.

    One shuffled alignment gives thousands of overlapping windows, which is
    plenty for the two moments and costs one alignment rather than twenty.
    The shuffle preserves the target's composition, so the null is what this
    aligner would produce on a sequence of the same amino-acid content and no
    shared ancestry.
    """
    window = int(_P.value("homology.site_window"))
    rng = np.random.default_rng(seed)
    letters = np.frombuffer(target_sequence.encode(), dtype="S1")
    shuffled = rng.permutation(letters).tobytes().decode()
    top, bottom = align_global(source_sequence, shuffled)
    scores, scored = column_scores(top, bottom)
    means = np.array([_window_mean(scores, scored, c, window)
                      for c in range(0, len(scores), max(window // 2, 1))])
    means = means[np.isfinite(means)]
    if means.size < 2:
        return 0.0, 1.0
    return float(means.mean()), float(means.std(ddof=1))


@dataclass(frozen=True)
class AlignmentWindows:
    """One source-to-target alignment, ready to be asked about positions."""

    target_of: dict
    column_of: dict
    scores: np.ndarray
    scored: np.ndarray
    identity: np.ndarray
    null_mean: float
    null_sd: float

    def at(self, residue: int) -> tuple[int | None, float, float]:
        """``(target residue, window score, window identity)``."""
        window = int(_P.value("homology.site_window"))
        column = self.column_of.get(int(residue))
        if column is None:
            return None, self.null_mean, 0.0
        half = window // 2
        lo, hi = max(0, column - half), min(len(self.scores), column + half + 1)
        n = float(self.scored[lo:hi].sum())
        return (self.target_of.get(int(residue)),
                _window_mean(self.scores, self.scored, column, window),
                float(self.identity[lo:hi].sum()) / n if n else 0.0)


_WINDOWS: dict = {}


def alignment_windows(target: str, source: str = "human") -> AlignmentWindows:
    """Align source to target once, and keep everything a position needs.

    One alignment per pair rather than one per residue. That is not only speed:
    a per-residue alignment could place two neighbouring curated residues
    inconsistently, and a site set whose members disagree about the register
    would be worse than useless.
    """
    from ..core.numbering_check import reference_entry

    key = (source, target, _P.value("homology.gap_open"),
           _P.value("homology.site_window"))
    if key in _WINDOWS:
        return _WINDOWS[key]

    source_sequence = reference_entry(source)["sequence"]
    target_sequence = reference_entry(target)["sequence"]
    top, bottom = align_global(source_sequence, target_sequence)
    scores, scored = column_scores(top, bottom)
    identity = np.array([1.0 if (a == b and a != "-") else 0.0
                         for a, b in zip(top, bottom)])

    column_of: dict[int, int] = {}
    target_of: dict[int, int | None] = {}
    i = j = 0
    for column, (a, b) in enumerate(zip(top, bottom)):
        if a != "-":
            i += 1
            column_of[i] = column
            target_of[i] = (j + 1) if b != "-" else None
        if b != "-":
            j += 1

    null_mean, null_sd = window_null_distribution(
        source_sequence, target_sequence)
    windows = AlignmentWindows(
        target_of=target_of, column_of=column_of, scores=scores,
        scored=scored, identity=identity,
        null_mean=null_mean, null_sd=null_sd)
    _WINDOWS[key] = windows
    return windows


def map_positions(target: str, residues, source: str = "human"):
    """Map human residue numbers onto ``target``, with the window statistic.

    Returns ``(mapping, window_score, window_identity)``.
    """
    windows = alignment_windows(target, source)
    mapping, scores, identities = {}, {}, {}
    for residue in residues:
        mapped, score, identity = windows.at(int(residue))
        mapping[int(residue)] = mapped
        scores[int(residue)] = score
        identities[int(residue)] = identity
    return mapping, scores, identities


def window_identity(source_sequence: str, target_sequence: str,
                    column: int) -> float:
    """The identity of one aligned window — a calibration entry point.

    Exposed so a test can drive it on sequences whose answer is known by
    construction, which is the only reason to believe it on real ones.
    """
    window = int(_P.value("homology.site_window"))
    top, bottom = align_global(source_sequence, target_sequence)
    half = window // 2
    lo, hi = max(0, column - half), min(len(top), column + half + 1)
    pairs = [(a, b) for a, b in zip(top[lo:hi], bottom[lo:hi])
             if a != "-" and b != "-"]
    return sum(a == b for a, b in pairs) / max(len(pairs), 1)


def window_score(source_sequence: str, target_sequence: str,
                 column: int) -> float:
    """The BLOSUM62 window score at one column — the other calibration entry."""
    top, bottom = align_global(source_sequence, target_sequence)
    scores, scored = column_scores(top, bottom)
    return _window_mean(scores, scored, column,
                        int(_P.value("homology.site_window")))
