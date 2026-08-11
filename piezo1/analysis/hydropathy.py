"""Hydropathy along the chain — Guo & MacKinnon 2017, Figure 3-supplements 1-3.

Those three supplements are a Kyte-Doolittle trace of the whole 2547-residue
mouse Piezo1 sequence, and they carry an argument the rest of the paper depends
on. The structure resolves 24 blade helices in six 4-TM units; the authors
extend that to **nine** units and 38 helices for the full protein — including
twelve N-terminal helices nobody has ever seen — on the strength of the
hydropathy pattern continuing to the N-terminus.

That inference is load-bearing for this project too. It is why
``domains.json`` defines nine THUs, why the full-length model grafts a distal
blade at all, and why ``uniprot_*.json`` records 38 transmembrane segments. So
it is worth being able to recompute rather than cite.

Three things live here:

* :func:`hydropathy_profile` — the Kyte-Doolittle window average, the curve the
  supplements plot.
* :func:`predict_segments` — runs above threshold, i.e. the panel read as a
  prediction rather than a picture.
* :func:`repeat_periodicity` — the *claim*: that the helices come in fours.
  Measured as the loop length between consecutive helices against position
  modulo four, with a shuffled control, because "you can see the repeat" is
  not a measurement and a 4-periodicity is exactly the kind of pattern the eye
  invents.

**What this is not.** A window average over a hydrophobicity scale is a 1982
method and is comfortably beaten by any modern topology predictor. It is here
because it is what the figure is, not because it is the best available answer;
:func:`compare_with_reference` measures how much worse it is than the UniProt
annotation this project actually uses, so nothing downstream is tempted to
prefer it.

Numbering: whatever the reference sequence uses. Every returned position is a
1-based residue number in that sequence's own numbering, and
:class:`HydropathyProfile` records which.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

from ..config import RESOURCE_DIR
from ..parameters import PARAMETERS as _P

__all__ = ["KYTE_DOOLITTLE", "HydropathyProfile", "Segment", "SegmentAgreement",
           "RepeatPeriodicity", "hydropathy_profile", "predict_segments",
           "compare_with_reference", "repeat_periodicity", "load_reference",
           "threshold_scan", "annotated_hydropathy"]


#: Kyte & Doolittle (1982) J Mol Biol 157:105 Table 1 — the hydropathy index.
#: Positive is hydrophobic. Reproduced exactly; a scale with a single wrong
#: entry produces a curve that still looks right, which is why the test checks
#: the published extremes (Ile +4.5, Arg -4.5) and the published mean.
KYTE_DOOLITTLE = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9,
    "R": -4.5,
}

#: Unknown residues score zero — neither hydrophobic nor hydrophilic — rather
#: than being dropped, so a position never shifts the numbering of the curve.
_UNKNOWN = 0.0


@dataclass
class Segment:
    """A run of the smoothed curve above threshold."""

    start: int                 # 1-based, inclusive
    end: int                   # 1-based, inclusive
    peak: float                # highest window average in the run
    mean: float

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def overlaps(self, other: "Segment | tuple[int, int]") -> int:
        """Residues shared with another segment or ``(start, end)`` pair."""
        lo, hi = (other.start, other.end) if isinstance(other, Segment) else other
        return max(0, min(self.end, hi) - max(self.start, lo) + 1)


@dataclass
class HydropathyProfile:
    """The Kyte-Doolittle trace of one sequence."""

    #: 1-based residue number of each point — the window centre.
    position: np.ndarray
    #: Window-averaged hydropathy at that position.
    value: np.ndarray
    #: Per-residue index before smoothing, same length.
    raw: np.ndarray
    window: int
    reference: str
    sequence_length: int
    meta: dict = field(default_factory=dict)

    def at(self, residue: int) -> float:
        """Smoothed hydropathy at a residue number, NaN outside the sequence."""
        idx = int(residue) - 1
        if not (0 <= idx < len(self.value)):
            return float("nan")
        return float(self.value[idx])

    def mean_over(self, start: int, end: int) -> float:
        lo, hi = max(1, int(start)), min(self.sequence_length, int(end))
        if hi < lo:
            return float("nan")
        return float(np.nanmean(self.value[lo - 1:hi]))


def load_reference(reference: str = "mouse") -> dict:
    """One of the committed ``uniprot_*.json`` resources, by short name."""
    path = RESOURCE_DIR / f"uniprot_{reference}.json"
    if not path.exists():
        raise FileNotFoundError(
            f"no committed reference {reference!r} — expected {path.name}")
    return json.loads(path.read_text())


# --------------------------------------------------------------------------
# The curve
# --------------------------------------------------------------------------

def hydropathy_profile(sequence: str | None = None, window: int | None = None,
                       reference: str = "mouse") -> HydropathyProfile:
    """Kyte-Doolittle window average along a sequence.

    Defaults to the committed mouse UniProt sequence, which is the one the
    paper's supplements plot (E2JF22, 2547 aa).

    The window is centred and **truncated** at the termini rather than padded:
    padding with zeros would pull the first and last residues towards neutral
    and invent a hydrophilic N-terminus that the sequence does not have. The
    number of residues actually averaged at each position is recorded in
    ``meta["window_used"]`` so an edge value is never mistaken for a full one.
    """
    if window is None:
        window = int(_P.value("hydropathy.window"))
    if window < 1:
        raise ValueError("window must be at least one residue")
    if sequence is None:
        ref = load_reference(reference)
        sequence = ref["sequence"]
        length = int(ref["length"])
    else:
        length = len(sequence)

    raw = np.array([KYTE_DOOLITTLE.get(aa.upper(), _UNKNOWN) for aa in sequence],
                   dtype=np.float64)
    n = len(raw)
    half = window // 2
    # Cumulative sums give the truncated centred mean in one pass, which
    # matters: this is called across five sequences by the periodicity control.
    cumulative = np.concatenate([[0.0], np.cumsum(raw)])
    lo = np.clip(np.arange(n) - half, 0, n)
    hi = np.clip(np.arange(n) + half + 1, 0, n)
    counts = (hi - lo).astype(np.float64)
    smoothed = (cumulative[hi] - cumulative[lo]) / counts

    return HydropathyProfile(
        position=np.arange(1, n + 1), value=smoothed, raw=raw, window=window,
        reference=reference, sequence_length=length,
        meta={"scale": "Kyte-Doolittle 1982",
              "window_used": counts,
              "n_unknown": int(sum(1 for aa in sequence
                                   if aa.upper() not in KYTE_DOOLITTLE)),
              "citation": "kyte1982"})


def predict_segments(profile: HydropathyProfile, threshold: float | None = None,
                     min_length: int | None = None) -> list[Segment]:
    """Runs above ``threshold`` at least ``min_length`` residues long."""
    if threshold is None:
        threshold = _P.value("hydropathy.tm_threshold")
    if min_length is None:
        min_length = int(_P.value("hydropathy.min_tm_length"))

    above = profile.value >= threshold
    segments: list[Segment] = []
    start: int | None = None
    for i, flag in enumerate(above):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            if i - start >= min_length:
                block = profile.value[start:i]
                segments.append(Segment(start + 1, i, float(block.max()),
                                        float(block.mean())))
            start = None
    if start is not None and len(above) - start >= min_length:
        block = profile.value[start:]
        segments.append(Segment(start + 1, len(above), float(block.max()),
                                float(block.mean())))
    return segments


# --------------------------------------------------------------------------
# Against the annotation the project actually uses
# --------------------------------------------------------------------------

@dataclass
class SegmentAgreement:
    """How a hydropathy call compares with the curated transmembrane table."""

    n_predicted: int
    n_annotated: int
    #: Annotated helices with at least one residue of predicted overlap.
    n_recovered: int
    #: Predicted segments overlapping no annotated helix.
    n_spurious: int
    #: Residue-level Jaccard index between the two sets of segments.
    jaccard: float
    #: Annotated helices found by no predicted segment, by name.
    missed: tuple[str, ...]
    threshold: float
    note: str = ""

    @property
    def recall(self) -> float:
        return self.n_recovered / self.n_annotated if self.n_annotated else 0.0


def compare_with_reference(profile: HydropathyProfile | None = None,
                           reference: str = "mouse",
                           threshold: float | None = None
                           ) -> SegmentAgreement:
    """Hydropathy segments against the UniProt transmembrane table.

    The point is not to validate UniProt — it is to bound what a 1982 window
    average is worth on this protein, so that the paper's inference from it
    (nine 4-TM units, twelve unseen N-terminal helices) is read with the
    accuracy of the method it was made with.
    """
    profile = profile or hydropathy_profile(reference=reference)
    if threshold is None:
        threshold = _P.value("hydropathy.tm_threshold")
    ref = load_reference(reference)
    annotated = sorted(ref["transmembrane"], key=lambda t: t["start"])
    predicted = predict_segments(profile, threshold=threshold)

    covered_pred = np.zeros(profile.sequence_length + 2, dtype=bool)
    for seg in predicted:
        covered_pred[seg.start:seg.end + 1] = True
    covered_ann = np.zeros(profile.sequence_length + 2, dtype=bool)
    for tm in annotated:
        covered_ann[tm["start"]:tm["end"] + 1] = True

    missed = tuple(tm.get("name", f"{tm['start']}-{tm['end']}")
                   for tm in annotated
                   if not covered_pred[tm["start"]:tm["end"] + 1].any())
    recovered = len(annotated) - len(missed)
    spurious = sum(1 for seg in predicted
                   if not covered_ann[seg.start:seg.end + 1].any())
    union = int((covered_pred | covered_ann).sum())
    intersection = int((covered_pred & covered_ann).sum())

    return SegmentAgreement(
        n_predicted=len(predicted), n_annotated=len(annotated),
        n_recovered=recovered, n_spurious=spurious,
        jaccard=intersection / union if union else 0.0,
        missed=missed, threshold=float(threshold),
        note="Kyte-Doolittle window average against the curated UniProt "
             "table; the annotation is the reference, not the other way round")


def annotated_hydropathy(profile: HydropathyProfile | None = None,
                         reference: str = "mouse") -> dict:
    """Where the annotated transmembrane helices actually sit on the scale.

    The measurement that explains the threshold result, and the reason the
    default cut is left at Kyte & Doolittle's published +1.6 rather than tuned
    down until the recall looks good. PIEZO1's transmembrane helices average
    well below the conventional membrane-spanning cut: they are genuinely
    polar helices, which is consistent with the paper's own description of
    "hydrophobic residues on the TM helices flanked by charged amino acids"
    and with the positive-inside rule it invokes (von Heijne 1992).

    Tuning the threshold to this protein would make the agreement number
    circular — it would be chosen to maximise the agreement it then reports —
    so the whole curve is available through :func:`threshold_scan` instead.
    """
    profile = profile or hydropathy_profile(reference=reference)
    ref = load_reference(reference)
    helices = sorted(ref["transmembrane"], key=lambda t: t["start"])
    means = np.array([profile.mean_over(t["start"], t["end"]) for t in helices])
    peaks = np.array([np.nanmax(profile.value[t["start"] - 1:t["end"]])
                      for t in helices])
    default = _P.value("hydropathy.tm_threshold")
    # Everything outside an annotated helix, as the contrast the segment call
    # is really being asked to make.
    inside = np.zeros(profile.sequence_length, dtype=bool)
    for t in helices:
        inside[t["start"] - 1:t["end"]] = True
    return {
        "reference": reference,
        "n_helices": len(helices),
        "mean_window_hydropathy": float(np.nanmean(means)),
        "median_window_hydropathy": float(np.nanmedian(means)),
        "mean_peak_hydropathy": float(np.nanmean(peaks)),
        "fraction_above_default_threshold":
            float(np.mean(means >= default)),
        "default_threshold": float(default),
        "outside_mean": float(np.nanmean(profile.value[~inside])),
        "separation": float(np.nanmean(means) - np.nanmean(profile.value[~inside])),
        "note": ("the helices sit above their surroundings but below the "
                 "conventional membrane-spanning cut, so a threshold call on "
                 "this protein trades recall for nothing"),
    }


def threshold_scan(profile: HydropathyProfile | None = None,
                   reference: str = "mouse",
                   thresholds: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0,
                                                    1.25, 1.5, 1.6, 1.75, 2.0)
                   ) -> list[dict]:
    """Agreement with the annotation across the whole threshold range.

    Reported instead of a single tuned operating point. Recall rises and
    specificity falls monotonically along it, and the reader can see both
    rather than being handed the threshold that flattered the method.
    """
    profile = profile or hydropathy_profile(reference=reference)
    rows = []
    for threshold in thresholds:
        agreement = compare_with_reference(profile, reference=reference,
                                           threshold=threshold)
        rows.append({"threshold": float(threshold),
                     "n_predicted": agreement.n_predicted,
                     "n_recovered": agreement.n_recovered,
                     "recall": agreement.recall,
                     "n_spurious": agreement.n_spurious,
                     "jaccard": agreement.jaccard})
    return rows


# --------------------------------------------------------------------------
# The claim: helices come in fours
# --------------------------------------------------------------------------

@dataclass
class RepeatPeriodicity:
    """Whether the transmembrane helices group into 4-TM units.

    The paper's argument in one number. If the blade is built from repeated
    4-TM units then the loop *between* units should be systematically longer
    than the loops *inside* one, because a unit's four helices are a compact
    left-handed bundle and the units are strung together by extended linkers.
    """

    #: Mean inter-helix loop length at each position in the repeat, 0-3.
    loop_by_phase: tuple[float, ...]
    #: Which phase carries the long loop. For PIEZO1 numbered from TM1 this is
    #: 3, i.e. the loop *after* every fourth helix — the unit boundary. It is
    #: measured rather than assumed because the register depends on where the
    #: reference's TM1 falls, and in PEZO-1 and dPIEZO it falls elsewhere.
    phase: int
    #: The long-loop phase mean, minus the mean of the other three, residues.
    contrast: float
    #: Same contrast over random regroupings of the same loops.
    control_mean: float
    control_sd: float
    #: How many control standard deviations the real contrast sits above.
    z: float
    n_helices: int
    n_units: int
    period: int
    reference: str

    @property
    def supported(self) -> bool:
        """True when the grouping beats its own shuffled control.

        Two standard deviations, one-sided: the alternative being excluded is
        that the loops are long and short at random and the fours were drawn
        onto them afterwards.
        """
        return self.z >= 2.0

    def summary(self) -> str:
        return (f"{self.n_helices} helices in {self.n_units} units of "
                f"{self.period}: the long loop falls at phase {self.phase} and "
                f"exceeds the other three by {self.contrast:.0f} residues, "
                f"control {self.control_mean:.0f} +- {self.control_sd:.0f} "
                f"(z = {self.z:.1f}) — "
                f"{'supported' if self.supported else 'NOT supported'}")


def repeat_periodicity(reference: str = "mouse", period: int = 4,
                       n_shuffles: int = 2000, seed: int = 0,
                       n_helices: int | None = None) -> RepeatPeriodicity:
    """Test the 4-TM repeat against a shuffled control.

    ``n_helices`` restricts to the first N transmembrane segments, which is how
    the pore helices are excluded: PIEZO1's last two (TM37, TM38) are not part
    of a 4-TM unit and including them puts a non-repeat at the end of the
    series. The default drops however many trailing helices do not complete a
    unit — 38 helices leaves 36 in nine units, with TM37-38 excluded, which is
    exactly the paper's arithmetic.

    The control shuffles the **loops**, not the helices: the null being tested
    is that the loops are the lengths they are but grouped arbitrarily, which
    is the only null that isolates the periodicity from the fact that PIEZO1
    simply has some long linkers.

    The statistic is maximised over the four possible registers rather than
    fixed at "the loop after every fourth helix", because where a repeat
    starts depends on where the reference annotation happens to call TM1 —
    and in PEZO-1 and dPIEZO it starts elsewhere, so a fixed register reports
    a real repeat as absent. The control is maximised over registers too. That
    matters: taking the best of four and comparing against an unmaximised null
    would manufacture roughly a standard deviation of significance from
    nothing, which is the failure mode this project has been bitten by before.
    """
    ref = load_reference(reference)
    helices = sorted(ref["transmembrane"], key=lambda t: t["start"])
    if n_helices is None:
        n_helices = (len(helices) // period) * period
    helices = helices[:n_helices]
    if len(helices) < 2 * period:
        raise ValueError(f"need at least {2 * period} helices to see a "
                         f"period-{period} repeat, got {len(helices)}")

    loops = np.array([helices[i + 1]["start"] - helices[i]["end"] - 1
                      for i in range(len(helices) - 1)], dtype=np.float64)
    # Loop i sits after helix i, so its phase is i mod period; phase
    # ``period - 1`` is the loop that leaves a unit.
    phase = np.arange(len(loops)) % period

    def contrast_of(values: np.ndarray) -> tuple[np.ndarray, int, float]:
        """Best register: the phase whose loops most exceed the other three."""
        means = np.array([values[phase == p].mean() for p in range(period)])
        total = means.sum()
        # For each candidate phase, its mean minus the mean of the rest.
        rest = (total - means) / (period - 1)
        gaps = means - rest
        best = int(np.argmax(gaps))
        return means, best, float(gaps[best])

    means, best_phase, observed = contrast_of(loops)

    rng = np.random.default_rng(seed)
    controls = np.empty(n_shuffles)
    for i in range(n_shuffles):
        controls[i] = contrast_of(rng.permutation(loops))[2]
    sd = float(controls.std(ddof=1))

    return RepeatPeriodicity(
        loop_by_phase=tuple(float(m) for m in means), phase=best_phase,
        contrast=observed, control_mean=float(controls.mean()), control_sd=sd,
        z=float((observed - controls.mean()) / sd) if sd > 0 else 0.0,
        n_helices=len(helices), n_units=len(helices) // period, period=period,
        reference=reference)
