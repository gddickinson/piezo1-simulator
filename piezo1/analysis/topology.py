"""The topology cartoon — Guo & MacKinnon 2017, Figure 3a.

Figure 3a is the diagram the rest of the paper is read against: 38 helices in a
row, grouped into nine 4-TM units, the first twelve greyed out because they are
not resolved, the CED as a box, the beam as a long bar underneath, and the cuff
elements tucked in beside the pore.

This module builds that diagram as **data** — elements with residue ranges, a
membrane side, a unit index and 2-D layout coordinates — so it can be drawn by
matplotlib in a script, inspected in a test, or served to the GUI, without any
of them owning the layout. Nothing here imports a plotting library.

Two things it does that a hand-drawn figure cannot:

* **It marks what is resolved.** Pass a structure and every element gets a
  ``resolved`` flag and a modelled-residue count, from that entry's actual
  coordinates. Figure 3a greys out TM1-12 for 6B3R; for 7WLT the greyed set is
  different, and the diagram follows rather than being redrawn.
* **It is built from the same annotation everything else uses.** The units come
  from ``domains.json``'s THUs and the helices from the committed UniProt
  table, so a diagram that disagrees with an analysis is a bug in one of them
  rather than two pictures of different things.

Layout convention: ``x`` increases along the sequence, ``y`` is across the
membrane with **positive extracellular** and negative cytoplasmic, matching
Figure 3a's orientation. Coordinates are dimensionless; a renderer scales them.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

from ..config import RESOURCE_DIR

__all__ = ["TopologyElement", "Topology", "build_topology", "MEMBRANE_HALF",
           "ELEMENT_KINDS"]

#: Half-thickness of the drawn membrane in layout units. A transmembrane helix
#: spans ``-MEMBRANE_HALF`` to ``+MEMBRANE_HALF``.
MEMBRANE_HALF = 1.0

#: The kinds a renderer must know how to draw.
ELEMENT_KINDS = ("tm_helix", "loop", "box", "bar", "cuff")


@dataclass
class TopologyElement:
    """One drawable piece of the topology diagram."""

    kind: str
    label: str
    start: int                    # first residue, in the diagram's numbering
    end: int
    #: Layout box: x from ``x0`` to ``x1``, y from ``y0`` to ``y1``.
    x0: float
    x1: float
    y0: float
    y1: float
    #: Which 4-TM unit this belongs to, 1-9; None for the pore module and the
    #: elements that are not part of a unit.
    unit: int | None = None
    #: 1-based transmembrane helix index, for ``tm_helix`` elements.
    helix: int | None = None
    side: str = "membrane"        # extracellular | cytoplasmic | membrane
    #: Whether the loaded structure models this element at all.
    resolved: bool | None = None
    n_modelled: int = 0
    color: str = "#888888"

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def as_dict(self) -> dict:
        return {"kind": self.kind, "label": self.label,
                "start": self.start, "end": self.end,
                "box": [self.x0, self.y0, self.x1, self.y1],
                "unit": self.unit, "helix": self.helix, "side": self.side,
                "resolved": self.resolved, "n_modelled": self.n_modelled,
                "color": self.color}


@dataclass
class Topology:
    """The whole diagram."""

    elements: list[TopologyElement]
    reference: str
    numbering: str
    sequence_length: int
    n_helices: int
    n_units: int
    structure: str | None = None
    meta: dict = field(default_factory=dict)

    def of_kind(self, kind: str) -> list[TopologyElement]:
        return [e for e in self.elements if e.kind == kind]

    @property
    def unresolved_helices(self) -> tuple[int, ...]:
        """Helix indices the loaded structure does not model at all."""
        return tuple(e.helix for e in self.of_kind("tm_helix")
                     if e.helix is not None and e.resolved is False)

    def summary(self) -> str:
        helices = self.of_kind("tm_helix")
        drawn = sum(1 for e in helices if e.resolved is not False)
        return (f"{self.n_helices} transmembrane helices in {self.n_units} "
                f"4-TM units plus the pore module"
                + (f"; {drawn}/{len(helices)} modelled in {self.structure}"
                   if self.structure else ""))

    def as_dict(self) -> dict:
        return {"reference": self.reference, "numbering": self.numbering,
                "sequence_length": self.sequence_length,
                "n_helices": self.n_helices, "n_units": self.n_units,
                "structure": self.structure,
                "unresolved_helices": list(self.unresolved_helices),
                "n_elements": len(self.elements),
                "elements": [e.as_dict() for e in self.elements],
                "summary": self.summary(), "meta": self.meta}


def _modelled_residues(structure, lo: int, hi: int) -> int:
    """C-alpha count in a residue range, taken from one protomer.

    One protomer, not all three: a trimer would treble every count and the
    number is meant to answer "is this element modelled", not "how many copies".
    """
    if structure is None:
        return 0
    chains = [c for c in structure.chains
              if (structure.mask_ca() & (structure.chain == c)).sum() > 300]
    if not chains:
        chains = list(structure.chains[:1])
    mask = (structure.mask_ca() & (structure.chain == chains[0])
            & (structure.res_seq >= lo) & (structure.res_seq <= hi))
    return int(mask.sum())


def build_topology(reference: str = "mouse", structure=None) -> Topology:
    """Assemble Figure 3a's diagram for a reference sequence.

    ``structure`` is optional. With one, every element is marked resolved or
    not against that entry's coordinates — which is what makes the diagram a
    statement about a particular model rather than about the protein in
    general.
    """
    uniprot = json.loads(
        (RESOURCE_DIR / f"uniprot_{reference}.json").read_text())
    domains = json.loads((RESOURCE_DIR / "domains.json").read_text())["domains"]
    key = "mouse" if reference.startswith("mouse") else "human"
    helices = sorted(uniprot["transmembrane"], key=lambda t: t["start"])
    n_helices = len(helices)
    period = 4
    n_units = (n_helices - 2) // period      # the last two are the pore module

    #: Blues for the blade units, warm for the pore module — matching the
    #: rainbow of Figure 3a well enough to be recognisable without pretending
    #: to be it.
    unit_colors = ["#1f3f7a", "#26508f", "#2d61a4", "#3472b9", "#3b83ce",
                   "#4a97dd", "#5faee2", "#79c2e0", "#95d3d8"]

    elements: list[TopologyElement] = []
    x = 0.0
    step, width, unit_gap = 1.0, 0.62, 1.1

    previous: TopologyElement | None = None
    for index, helix in enumerate(helices, start=1):
        unit = (index - 1) // period + 1 if index <= n_units * period else None
        if index > 1:
            x += step + (unit_gap if unit is not None
                         and (index - 1) % period == 0 else 0.0)
        colour = (unit_colors[(unit - 1) % len(unit_colors)] if unit
                  else "#c0392b")
        n_modelled = _modelled_residues(structure, helix["start"], helix["end"])
        element = TopologyElement(
            kind="tm_helix",
            label=helix.get("name", f"TM{index}"),
            start=int(helix["start"]), end=int(helix["end"]),
            x0=x, x1=x + width, y0=-MEMBRANE_HALF, y1=MEMBRANE_HALF,
            unit=unit, helix=index, side="membrane",
            resolved=(None if structure is None
                      else n_modelled >= max(3, (helix["end"] - helix["start"]) // 3)),
            n_modelled=n_modelled, color=colour)
        elements.append(element)

        if previous is not None:
            # A helix crossing the membrane alternates the side its following
            # loop sits on. Taken from the helix parity rather than from the
            # UniProt topology strings, which are absent for some references.
            side = "extracellular" if index % 2 == 0 else "cytoplasmic"
            y = MEMBRANE_HALF if side == "extracellular" else -MEMBRANE_HALF
            lo, hi = previous.end + 1, element.start - 1
            loop_modelled = _modelled_residues(structure, lo, hi) if hi >= lo else 0
            elements.append(TopologyElement(
                kind="loop", label=f"{previous.label}-{element.label}",
                start=lo, end=max(hi, lo),
                x0=previous.x1, x1=element.x0,
                y0=y, y1=y + (0.45 if y > 0 else -0.45),
                unit=previous.unit if previous.unit == element.unit else None,
                side=side,
                resolved=(None if structure is None else loop_modelled > 0),
                n_modelled=loop_modelled, color="#9aa0a6"))
        previous = element

    right = elements[-1].x1

    # The cap sits above the last two helices; the cuff and beam below.
    #
    # Placement is by where the residue range falls in the sequence, which puts
    # the four cuff elements in two tight pairs — the elbow and base 33
    # residues apart, the PE helix and hairpin 22 — and at this scale their
    # labels sat on top of each other and read as "bo|as" and "PE|irp". Each
    # therefore gets a minimum width, is pushed right of whatever precedes it,
    # and alternates between two rows. Their x is an ordering, not a
    # measurement, which is why moving them is legitimate; the residue range
    # each carries is exact and is what the tooltip and the selection use.
    span = max(uniprot["length"], 1)
    lower_rows = (-MEMBRANE_HALF - 0.75, -MEMBRANE_HALF - 1.35)
    min_width, gap_between = 1.6, 0.25
    placed_end = {row: -1e9 for row in lower_rows}

    cuff_specs = [
        ("cap", "box", "CED", None, "#e74c3c"),
        ("beam", "bar", "Beam", None, "#2ecc71"),
        ("elbow", "cuff", "Elbow", 0, "#c9a227"),
        ("base", "cuff", "Base", 1, "#b8860b"),
        ("pore_extension", "cuff", "PE helix", 0, "#d2691e"),
        ("hairpin", "cuff", "Hairpin", 1, "#a0522d"),
    ]
    for domain_id, kind, label, row_index, colour in cuff_specs:
        record = next((d for d in domains if d["id"] == domain_id), None)
        if record is None:
            continue
        lo, hi = record[key]["start"], record[key]["end"]
        if kind == "box":
            x0, x1 = right - 2.4, right - 0.2
            y0, y1 = MEMBRANE_HALF + 0.5, MEMBRANE_HALF + 1.9
        else:
            x0 = right * lo / span
            x1 = max(x0 + min_width, right * hi / span)
            if kind == "bar":
                # The beam is short in sequence and long in space — 66 residues
                # spanning the whole arm — so Figure 3a draws it as a long bar.
                x1 = max(x1, x0 + 4.0)
                y0, y1 = -MEMBRANE_HALF - 2.1, -MEMBRANE_HALF - 1.75
            else:
                row = lower_rows[row_index]
                if x0 < placed_end[row] + gap_between:
                    width = x1 - x0
                    x0 = placed_end[row] + gap_between
                    x1 = x0 + width
                placed_end[row] = x1
                y0, y1 = row, row + 0.42
        n_modelled = _modelled_residues(structure, lo, hi)
        elements.append(TopologyElement(
            kind=kind, label=label, start=lo, end=hi,
            x0=x0, x1=x1, y0=y0, y1=y1, unit=None,
            side=("extracellular" if kind == "box" else "cytoplasmic"),
            resolved=(None if structure is None else n_modelled > 0),
            n_modelled=n_modelled, color=colour))

    return Topology(
        elements=elements, reference=reference, numbering=key,
        sequence_length=int(uniprot["length"]), n_helices=n_helices,
        n_units=n_units,
        structure=(structure.name if structure is not None else None),
        meta={"period": period,
              "pore_helices": [n_helices - 1, n_helices],
              "membrane_half": MEMBRANE_HALF,
              # Over *every* element, not just the helices: the cuff elements
              # are pushed right to stop their labels colliding and can end
              # past the last helix, and a range that stopped there clipped the
              # hairpin off the edge of the widget.
              "x_range": [0.0, max((e.x1 for e in elements), default=right)],
              "note": ("Guo & MacKinnon number the units from the N-terminus "
                       "— their '4-TM unit 6' is TM21-24 — which is the same "
                       "convention domains.json uses for THU1-9."),
              "citation": "guo2017"})


def unit_extent(topology: Topology) -> dict[int, tuple[float, float]]:
    """Layout x-range of each 4-TM unit, for drawing the red boxes of Fig 3b."""
    out: dict[int, list[float]] = {}
    for element in topology.of_kind("tm_helix"):
        if element.unit is None:
            continue
        span = out.setdefault(element.unit, [element.x0, element.x1])
        span[0] = min(span[0], element.x0)
        span[1] = max(span[1], element.x1)
    return {unit: (lo, hi) for unit, (lo, hi) in out.items()}
