"""Named parts of the assembly, and the residues worth seeing on each.

A PIEZO1 trimer is 120,000 atoms and about 300 A across. Almost every question
anyone asks of it is about one part — the pore module, the cap and its gates,
the beam — and showing the whole propeller to answer them buries the answer in
blades. This module names the parts.

**Every component is built from curated annotation, not from ranges written
down here.** The domains come from ``domains.json`` and the residues to
highlight from ``functional_residues.json``, both of which carry provenance and
both of which state their numbering. A component therefore inherits whatever
those files say, including the corrections they have had, and cannot drift away
from the annotation the rest of the project reads.

Qt-free on purpose: the definitions are data, the masks are numpy, and the
controller that draws them is the only part that needs a window. That is what
lets a test check *which residues a component selects* without a display.

The one that matters most is ``pore_module`` — Liu et al. 2025's Figure 2E view:
the outer helix, the cap, the spring linker and the inner helix enclosing the
four vestibules, with the gates picked out. It is the view in which their whole
argument is legible, and it is 4% of the atoms.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["Component", "COMPONENTS", "component_by_key", "component_masks",
           "ComponentSelection"]


@dataclass(frozen=True)
class Component:
    """One named part, as the annotation ids that make it up."""

    key: str
    label: str
    #: What it is for, in one line. Shown on the status bar.
    shows: str
    #: Curated domain ids whose residues form the backbone drawn.
    domains: tuple[str, ...] = ()
    #: Curated residue-group ids picked out on top of that backbone.
    highlight: tuple[str, ...] = ()
    #: What the view must not be read as, when that is a real risk.
    caveat: str = ""

    @property
    def is_whole(self) -> bool:
        return not self.domains


#: The order they appear in the menu. ``whole`` first, because it is the
#: default and the way back.
COMPONENTS: tuple[Component, ...] = (
    Component(
        "whole", "Whole assembly",
        "Everything the entry resolves, which is the default."),
    Component(
        "pore_module", "Pore module (Liu et al. Figure 2E)",
        "The conduction path end to end: outer helix, cap, spring linker, "
        "inner helix and CTD, with all four gates picked out.",
        domains=("outer_helix", "cap", "inner_helix", "ctd"),
        highlight=("cap_constriction", "cap_gate", "spring_linker",
                   "hydrophobic_gate", "ctd_constriction"),
        caveat="the four cavities are enclosed by this module, but ions enter "
               "and leave it SIDEWAYS — see the ion flux pathway"),
    Component(
        "cap_and_gates", "Cap and its gates",
        "The extracellular cap, with the closed top and the two lateral "
        "cap-gate loops that open on the transition to the intermediate state.",
        domains=("cap",),
        highlight=("cap_constriction", "cap_gate", "cap_gate_loops")),
    Component(
        "tm_gate", "Transmembrane gate",
        "The inner and outer helices, with the hydrophobic gate and the "
        "pore-lining residues.",
        domains=("outer_helix", "inner_helix"),
        highlight=("hydrophobic_gate", "pore_lining")),
    Component(
        "spring_and_ih", "Spring linker and inner helix",
        "The linker that compresses on opening, and the helix it pulls.",
        domains=("cap", "inner_helix"),
        highlight=("spring_linker", "hydrophobic_gate")),
    Component(
        "ctd_vestibule", "Cytoplasmic vestibule",
        "The C-terminal domain beneath the pore, with the constriction neck "
        "and the acidic residues the selectivity is attributed to.",
        domains=("ctd", "pore_extension", "hairpin"),
        highlight=("ctd_constriction", "selectivity_acidic"),
        caveat="three of the four curated selectivity glutamates are not "
               "within side-chain reach of the lumen"),
    Component(
        "beam_and_latch", "Beam and lateral plug",
        "The beam that levers the blade against the pore, its coiled coil, "
        "and the spliced segment that forms the intracellular lateral plug "
        "gate — the region Liu et al. delete in their Figure 5G.",
        domains=("beam", "coiled_coil", "splice_1_1")),
    Component(
        "anchor", "Anchor domain",
        "The anchor and its elbow and base helices, with the apex brake.",
        domains=("anchor", "elbow", "base"),
        highlight=("anchor_brake",)),
    Component(
        "blade", "Blade (THU1-THU9)",
        "The nine four-transmembrane units, with the basic clusters that "
        "contact phosphoinositides.",
        domains=("thu1", "thu2", "thu3", "thu4", "thu5", "thu6", "thu7",
                 "thu8", "thu9"),
        highlight=("pip2_cluster", "basic_cluster_2", "basic_cluster_3",
                   "basic_cluster_4"),
        caveat="most entries resolve only THU4-THU9; what is missing here is "
               "missing from the deposition, not hidden"),
    Component(
        "md_construct", "MD construct (Liu et al. Figure 5)",
        "Exactly what they simulated: the pore module, the beam and the "
        "lateral plug gate, and nothing else.",
        domains=("thu9", "anchor", "outer_helix", "cap", "inner_helix", "ctd",
                 "beam", "coiled_coil", "splice_1_1"),
        highlight=("hydrophobic_gate", "ctd_constriction"),
        caveat="a construct, not a structure — the blades were truncated away "
               "for the simulation and are not absent from the real channel"),
)


@dataclass
class ComponentSelection:
    """Which atoms a component selects on one structure."""

    component: Component
    backbone: np.ndarray                      # per-atom bool, what to draw
    highlight: np.ndarray                     # per-atom bool, what to pick out
    residues: tuple[int, ...] = ()            # residue numbers drawn
    missing: tuple[str, ...] = ()             # annotation ids not resolved
    numbering: str = ""
    note: str = ""

    @property
    def n_atoms(self) -> int:
        return int(self.backbone.sum())

    def summary(self) -> str:
        text = (f"{self.component.label}: {self.n_atoms:,} atoms, "
                f"{len(self.residues):,} residues")
        if self.numbering:
            text += f" ({self.numbering} numbering)"
        if self.missing:
            text += f" · not resolved here: {', '.join(self.missing)}"
        if self.component.caveat:
            text += f" · {self.component.caveat}"
        return text


def component_by_key(key: str) -> Component:
    for component in COMPONENTS:
        if component.key == key:
            return component
    raise KeyError(f"no component {key!r}; have "
                   f"{[c.key for c in COMPONENTS]}")


def component_masks(structure, key: str, numbering: str | None = None
                    ) -> ComponentSelection:
    """Which atoms of ``structure`` belong to the named component.

    ``numbering`` is measured from the coordinates when not given, and a
    structure whose numbering cannot be read falls back to the whole assembly
    with the reason recorded — every range here is a residue number, and
    applying human ranges to a mouse entry would select the wrong helices while
    looking entirely plausible.
    """
    from ..core.annotations import load_annotations
    from ..core.numbering_check import piezo1_numbering

    component = component_by_key(key)
    everything = np.ones(structure.n_atoms, dtype=bool)
    nothing = np.zeros(structure.n_atoms, dtype=bool)

    if component.is_whole:
        return ComponentSelection(component, everything, nothing,
                                  residues=tuple(sorted(
                                      set(int(v) for v in structure.res_seq))))

    measured = numbering or piezo1_numbering(structure)
    if measured is None:
        return ComponentSelection(
            component, everything, nothing, numbering="",
            note="numbering not readable from these coordinates — showing the "
                 "whole assembly rather than selecting the wrong residues")

    annotations = load_annotations(measured)
    wanted: set[int] = set()
    missing: list[str] = []
    for domain_id in component.domains:
        domain = _domain(annotations, domain_id)
        if domain is None:
            missing.append(domain_id)
            continue
        wanted.update(range(int(domain.start), int(domain.end) + 1))

    picked: set[int] = set()
    for group_id in component.highlight:
        group = annotations.group(group_id)
        if group is None:
            missing.append(group_id)
            continue
        picked.update(int(r) for r in group.residues)

    backbone = np.isin(structure.res_seq, sorted(wanted))
    # A highlighted residue is drawn whether or not a domain claimed it: the
    # PIP2 lysines sit just outside THU9's range, and a residue picked out as
    # important that is then not drawn is the worst of both.
    highlight = np.isin(structure.res_seq, sorted(picked))
    backbone = backbone | highlight

    present = sorted(set(int(v) for v in structure.res_seq[backbone]))
    return ComponentSelection(
        component=component, backbone=backbone, highlight=highlight,
        residues=tuple(present), missing=tuple(missing), numbering=measured)


def _domain(annotations, domain_id: str):
    for domain in annotations.domains:
        if domain.id == domain_id:
            return domain
    return None
