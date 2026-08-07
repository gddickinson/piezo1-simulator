"""Harvesting mutant electrophysiology from the open-access corpus.

Every route to this project's central claim has ended in the same place: not
enough phenotyped variants. Round 36 needed roughly 130 and had 34. New
experiments are not available, so the only remaining source is the published
literature the project has already downloaded.

**What this does and does not do.** It extracts *candidate* mutations from the
JATS full texts in ``ref/papers/``, gates them on the wild-type residue, resolves
which numbering system they are in, and attaches the sentence they appeared in
together with any measurement in it. It does **not** assign a direction. Reading
"slowed inactivation" out of prose and calling it gain-of-function is exactly the
kind of automated curation that would put unreviewed labels into the set the
project's blind tests depend on, and those labels are the one thing that must
stay hand-checked.

**The gate is not a formality.** Of 75 raw regex hits, **20 fail the wild-type
check against both human and mouse** — cDNA changes written like protein ones,
and matches from unrelated text. Of the 55 that pass, **34 are mouse-numbered**,
which is the project's standing trap: most functional literature uses mouse and
most disease variants use human, and the offset is not constant. Conversion goes
through :mod:`piezo1.core.sequence`, never by arithmetic.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from xml.etree import ElementTree as ET

__all__ = ["Candidate", "HarvestReport", "harvest", "extract_candidates",
           "SUBSTITUTION", "MEASUREMENT_PATTERNS"]

#: A protein substitution written in one-letter code. Deliberately loose, and
#: gated afterwards: a tighter pattern would miss real variants while still
#: admitting cDNA changes, which share the shape.
SUBSTITUTION = re.compile(r"\b([ACDEFGHIKLMNPQRSTVWY])(\d{3,4})"
                          r"([ACDEFGHIKLMNPQRSTVWY])\b")

#: Quantities that would make a candidate usable, with the units they carry.
#: Presence of one of these is recorded; interpreting it is a human's job.
MEASUREMENT_PATTERNS = {
    "inactivation_tau": re.compile(r"\b(?:tau|τ)\b[^.]{0,60}?(\d+\.?\d*)\s*ms",
                                   re.I),
    "half_activation": re.compile(r"\b(?:P50|P½|T50|half[- ]maximal)\b"
                                  r"[^.]{0,60}?(-?\d+\.?\d*)", re.I),
    "conductance": re.compile(r"(\d+\.?\d*)\s*pS", re.I),
    "current_density": re.compile(r"(-?\d+\.?\d*)\s*pA\s*/\s*pF", re.I),
}


@dataclass
class Candidate:
    """One substitution found in the literature, with its provenance."""

    label: str                      # as written, e.g. "R2456H"
    wt: str
    position: int
    mut: str
    source: str                     # citation key
    numbering: str                  # human | mouse | both | none
    human_label: str | None = None  # after conversion, when possible
    context: str = ""
    measurements: dict = field(default_factory=dict)
    already_curated: bool = False

    @property
    def passes_gate(self) -> bool:
        return self.numbering != "none"

    @property
    def usable(self) -> bool:
        """Gated, mappable to human, not already held, and carries a number."""
        return (self.passes_gate and self.human_label is not None
                and not self.already_curated and bool(self.measurements))

    def summary(self) -> str:
        return (f"{self.label} ({self.numbering}"
                + (f" -> {self.human_label}" if self.human_label
                   and self.human_label != self.label else "")
                + f", {self.source})"
                + (f" {sorted(self.measurements)}" if self.measurements else ""))


@dataclass
class HarvestReport:
    """Everything found, and what survives each gate."""

    candidates: list = field(default_factory=list)
    n_papers: int = 0
    meta: dict = field(default_factory=dict)

    def passing(self) -> list:
        return [c for c in self.candidates if c.passes_gate]

    def by_numbering(self) -> dict:
        out: dict[str, int] = {}
        for c in self.candidates:
            out[c.numbering] = out.get(c.numbering, 0) + 1
        return out

    def new_usable(self) -> list:
        return [c for c in self.candidates if c.usable]

    def summary(self) -> str:
        counts = self.by_numbering()
        return (f"{len(self.candidates)} candidates from {self.n_papers} "
                f"open-access papers; {len(self.passing())} pass the wild-type "
                f"gate ({counts}); {len(self.new_usable())} are new, mappable "
                f"to human and carry a measurement. Directions are NOT assigned "
                f"here — that stays hand-curated.")


def _numbering_of(wt: str, position: int) -> str:
    from ..core.sequence import human_sequence, mouse_sequence

    human, mouse = human_sequence(), mouse_sequence()

    def fits(seq: str) -> bool:
        return 1 <= position <= len(seq) and seq[position - 1] == wt

    in_human, in_mouse = fits(human), fits(mouse)
    if in_human and in_mouse:
        return "both"
    if in_human:
        return "human"
    if in_mouse:
        return "mouse"
    return "none"


def _sentence_around(text: str, index: int, width: int = 240) -> str:
    start = max(0, text.rfind(".", 0, index) + 1)
    end = text.find(".", index)
    end = len(text) if end < 0 else end + 1
    return re.sub(r"\s+", " ", text[start:end]).strip()[:width]


def extract_candidates(path) -> list:
    """Candidates from one JATS document, with the sentence each came from."""
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return []
    body = root.find(".//body")
    if body is None:
        return []
    text = ET.tostring(body, encoding="unicode", method="text")
    source = str(path).split("/")[-1].split("_")[0]

    seen, out = set(), []
    for match in SUBSTITUTION.finditer(text):
        wt, position, mut = match.group(1), int(match.group(2)), match.group(3)
        label = f"{wt}{position}{mut}"
        if label in seen:
            continue
        seen.add(label)
        context = _sentence_around(text, match.start())
        measurements = {}
        for name, pattern in MEASUREMENT_PATTERNS.items():
            found = pattern.search(context)
            if found:
                measurements[name] = float(found.group(1))
        out.append(Candidate(label=label, wt=wt, position=position, mut=mut,
                             source=source,
                             numbering=_numbering_of(wt, position),
                             context=context, measurements=measurements))
    return out


def harvest(paper_dir=None) -> HarvestReport:
    """Scan the downloaded open-access corpus and report what it yields."""
    from ..config import PAPER_DIR
    from ..core.annotations import load_annotations
    from ..core.sequence import mouse_to_human

    directory = paper_dir or PAPER_DIR
    files = sorted(directory.glob("*.xml")) if directory.exists() else []

    curated = {v.label for v in load_annotations("human").variants}
    candidates, papers = [], set()
    for path in files:
        found = extract_candidates(path)
        if found:
            papers.add(path.name.split("_")[0])
        candidates.extend(found)

    # Deduplicate on (label, source): the same variant discussed twice in one
    # paper is one candidate, but two papers reporting it is worth knowing.
    unique = {}
    for candidate in candidates:
        unique.setdefault((candidate.label, candidate.source), candidate)
    candidates = list(unique.values())

    for candidate in candidates:
        if candidate.numbering in ("human", "both"):
            candidate.human_label = candidate.label
        elif candidate.numbering == "mouse":
            human = mouse_to_human(candidate.position)
            if human:
                candidate.human_label = f"{candidate.wt}{human}{candidate.mut}"
        candidate.already_curated = bool(
            candidate.human_label and candidate.human_label in curated)

    return HarvestReport(candidates=candidates, n_papers=len(papers),
                         meta={"n_files": len(files),
                               "note": "candidates for curation, not labels"})
