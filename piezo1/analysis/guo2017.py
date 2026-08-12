"""Every figure panel of Guo & MacKinnon 2017, and what this project can do with it.

Guo & MacKinnon, *Structure-based membrane dome mechanism for Piezo
mechanosensitivity*, eLife 2017;6:e33660 (PMID 29231809, PDB 6B3R, EMD-7042) is
the paper this project's central claim comes from. The dome, the 10.2 nm radius,
the 120 nm² of projected area, the two-state energetics — all of it is Figure 7
and its supplement.

This is the registry of its panels. Each :class:`Panel` records what the
published panel shows, whether it can be reproduced **from coordinates**, and
either the callable that reproduces it or the reason it cannot be reproduced.

**16 of its 31 panels reproduce** from deposited coordinates, 3 have an
analogue that is a different quantity, and 12 need experimental data this
project does not hold. That count is stated in five places — here, the menu
tooltip, the README, ``docs/SCIENCE.md`` and the in-application help — and a
test makes a new panel move all of them or fail.

**Why the refusals are in the registry rather than left out.** Twelve panels
need experimental data this project does not have and never will — a Fourier
shell correlation needs two half-maps, micrographs need the liposomes. Listing
them as ``not_replicable`` with the reason is the difference between covering
the paper and covering its tractable parts while letting a reader assume the
rest. Same for ``analogue``: a projection of an atomic model is *not* a 2D
class average, and filing it as one would be the most misleading thing here.

Statuses:

``replicated``      computable from deposited coordinates, and where the panel
                    states a number, compared against it.
``analogue``        we can compute what the panel is an estimate *of*, but not
                    the panel. Read the caveat before showing it beside the
                    original.
``not_replicable``  needs experimental data this project does not hold.

Usage::

    from piezo1.analysis.guo2017 import PANELS, coverage, replicate
    coverage()                       # how much of the paper is reachable
    replicate("7-S1")                # one panel
    replicate("6a", structure=st)    # on a chosen entry

Numbering: the paper is **mouse** throughout (E2JF22, 2547 aa) and 6B3R is
deposited in mouse numbering. Every residue number in this module and in
:mod:`piezo1.analysis.guo2017_panels` and
:mod:`piezo1.analysis.guo2017_mechanism` is mouse unless it says otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from . import guo2017_mechanism as _mech
from . import guo2017_panels as _panels

__all__ = ["Panel", "PANELS", "PAPER", "STATUSES", "panel_by_key",
           "coverage", "replicate", "replicate_all", "not_replicable"]


#: The paper, for stamping onto anything this module produces.
PAPER = {
    "citation": "guo2017",
    "authors": "Guo YR, MacKinnon R",
    "title": "Structure-based membrane dome mechanism for Piezo mechanosensitivity",
    "journal": "eLife 2017;6:e33660",
    "doi": "10.7554/eLife.33660",
    "pmid": "29231809",
    "pdb": "6B3R",
    "emdb": "EMD-7042",
    "numbering": "mouse (UniProt E2JF22, 2547 aa)",
}

STATUSES = ("replicated", "analogue", "not_replicable")

#: The entry the paper itself deposited. Panels default to it, and a panel run
#: against anything else says so in its provenance.
DEFAULT_ENTRY = "6B3R"


@dataclass(frozen=True)
class Panel:
    """One published panel and this project's relationship to it."""

    key: str                      # "7a", "3-S1", "6-S1c"
    figure: str
    panel: str
    title: str
    #: What the published panel shows, in the paper's own terms.
    shows: str
    status: str
    #: The module that provides the capability, for the navigation map.
    module: str = ""
    #: What reproduces it. None for anything not ``replicated``/``analogue``.
    compute: Callable[..., dict] | None = None
    #: Why it cannot be reproduced, or what the analogue is not. Required for
    #: anything that is not ``replicated`` — enforced by a test, because an
    #: unexplained refusal is indistinguishable from an oversight.
    reason: str = ""
    #: Experimental data that would be needed. Empty when coordinates suffice.
    needs: tuple[str, ...] = ()
    #: Whether running it requires a loaded structure.
    needs_structure: bool = True
    #: Numbers the panel or its caption states, for the comparison to be
    #: against something rather than merely produced.
    published: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.status not in STATUSES:
            raise ValueError(f"{self.key}: bad status {self.status!r}")
        if self.status != "replicated" and not self.reason:
            raise ValueError(f"{self.key}: {self.status} needs a reason")
        if self.status == "not_replicable" and self.compute is not None:
            raise ValueError(f"{self.key}: not_replicable but has a callable")

    @property
    def label(self) -> str:
        return f"Figure {self.figure}{self.panel}".rstrip()


# --------------------------------------------------------------------------
# The registry
# --------------------------------------------------------------------------

PANELS: list[Panel] = [
    # ---- Figure 1: the gravitational analog -------------------------------
    Panel("1", "1", "", "In-plane area expansion under tension",
          "A cartoon of weights on a tensioned membrane, adapted from Ursell "
          "et al. 2008, showing why in-plane expansion lowers the free energy.",
          "not_replicable",
          reason="A conceptual illustration adapted from another paper. It "
                 "depicts no structure and contains no measurement, so there "
                 "is nothing for a structure to reproduce. The physics it "
                 "illustrates is equation 1 and is exercised by 7d.",
          needs=("nothing — it is a drawing",), needs_structure=False),

    # ---- Figure 2: the reconstruction -------------------------------------
    Panel("2ab", "2", "a,b", "2D class averages, top and side",
          "Representative 2D averaged classes from the particle stack, scale "
          "bar 10 nm.",
          "analogue", module="piezo1.analysis.projection",
          compute=_panels.panel_2ab,
          reason="A 2D class average is an average of thousands of real "
                 "particle images. What we compute is the projection of the "
                 "atomic model — the quantity a class average estimates — with "
                 "no CTF, no defocus, no solvent and, decisively for the side "
                 "view, no detergent micelle. Figure 2b's envelope is "
                 "substantially micelle.",
          needs=("micrographs or the particle stack",),
          published={"scale_bar_nm": 10, "pixel_size_A": 1.3}),
    Panel("2cd", "2", "c,d", "Trimer ribbon, top and side",
          "The atomic model as a ribbon diagram from the top and the side, the "
          "three subunits in red, green and blue.",
          "replicated", module="piezo1.render.representations",
          compute=_panels.panel_2cd,
          published={"residues_modelled": 1518, "of_total": 2547}),
    Panel("2e", "2", "e", "Fourier shell correlation",
          "FSC between two half maps after C1 masked refinement, crossing "
          "0.143 at 3.8 A.",
          "not_replicable",
          reason="An FSC is computed between two independently refined half "
                 "maps. This project holds neither half map and no atomic "
                 "model can produce one; a model-to-map FSC would be a "
                 "different quantity and is not what the panel shows.",
          needs=("the two half maps from EMD-7042",),
          published={"resolution_A": 3.8, "threshold": 0.143}),
    Panel("2f", "2", "f", "Local resolution",
          "The density map coloured by local resolution, estimated by Blocres.",
          "not_replicable",
          reason="Local resolution is a property of the reconstruction, not of "
                 "the model. It requires the map and the half maps.",
          needs=("EMD-7042 and its half maps",),
          published={"range_A": [3.2, 10.7]}),
    Panel("2-S1", "2-S1", "", "CryoEM structure determination",
          "The processing workflow: micrographs, 2D classes, ab initio "
          "reconstruction, symmetry expansion, masked refinement.",
          "not_replicable",
          reason="A record of how the reconstruction was done. It is a "
                 "methods diagram over experimental data, with nothing a "
                 "coordinate file could stand in for.",
          needs=("the raw micrographs and the processing pipeline",),
          needs_structure=False),
    Panel("2-S2", "2-S2", "", "Local EM densities, residues 581-1231",
          "Model segments shown inside their density.",
          "not_replicable",
          reason="Needs the sharpened map. The model half of the panel is "
                 "reproducible and the density half is the point of it.",
          needs=("EMD-7042",)),
    Panel("2-S3", "2-S3", "", "Local EM densities, residues 1280-2543",
          "Model segments shown inside their density.",
          "not_replicable",
          reason="Needs the sharpened map from EMD-7042. As with 2-S2 the "
                 "model half of the panel is reproducible and the density "
                 "half is the whole point of showing it.",
          needs=("EMD-7042",)),

    # ---- Figure 3: topology -----------------------------------------------
    Panel("3a", "3", "a", "Topology cartoon",
          "A monomer as a topology diagram: 38 helices in nine 4-TM units, "
          "the CED as a box, the cuff and beam marked, unresolved regions "
          "dotted. Insets: 4-TM unit 6 (TM21-24) and the pore region.",
          "replicated", module="piezo1.analysis.topology",
          compute=_panels.panel_3a, needs_structure=False,
          published={"n_helices": 38, "n_units": 9,
                     "unit_6": "TM21-24", "unresolved": "TM1-12"}),
    Panel("3b", "3", "b", "Monomer from the top, units boxed",
          "A ribbon of one monomer viewed down the axis with each 4-TM unit "
          "boxed and numbered.",
          "replicated", module="piezo1.analysis.topology",
          compute=_panels.panel_3a, needs_structure=False,
          published={"n_units_visible": 6}),
    Panel("3-S1", "3-S1", "", "Hydropathy, residues 1-900",
          "Kyte-Doolittle hydropathy along the sequence, the evidence that the "
          "4-TM pattern continues to the N-terminus.",
          "replicated", module="piezo1.analysis.hydropathy",
          compute=_panels.panel_3_supplements, needs_structure=False,
          published={"scale": "Kyte-Doolittle"}),
    Panel("3-S2", "3-S2", "", "Hydropathy, residues 901-1800",
          "As 3-S1, continued.",
          "replicated", module="piezo1.analysis.hydropathy",
          compute=_panels.panel_3_supplements, needs_structure=False),
    Panel("3-S3", "3-S3", "", "Hydropathy, residues 1801-2547",
          "As 3-S1, continued.",
          "replicated", module="piezo1.analysis.hydropathy",
          compute=_panels.panel_3_supplements, needs_structure=False),

    # ---- Figure 4: the curved micelle and the surface ---------------------
    Panel("4a", "4", "a", "A monomer in a planar membrane",
          "One protomer from the side with grey lines marking the approximate "
          "planar membrane interfaces — the claim that a subunit fits a plane "
          "and the trimer does not.",
          "replicated", module="piezo1.structure.planarity",
          compute=_panels.panel_4a,
          published={"beam_angle_deg": 60, "arms_out_of_plane_deg": 30}),
    Panel("4b", "4", "b", "Trimer in the curved micelle",
          "Ribbon diagrams inside an unsharpened map contoured at 6 sigma, "
          "showing the micelle curved into a dome. Top, side and bottom.",
          "analogue", module="piezo1.structure.micelle",
          compute=_panels.panel_4b,
          reason="The envelope in this panel is detergent micelle density from "
                 "the unsharpened map, which this project does not hold. What "
                 "we build is the surface at a fixed offset outside the "
                 "hydrophobic belt — the band the paper itself describes — so "
                 "its *thickness* is a parameter carrying no information and "
                 "only its *curvature* is a measurement of the protein. A "
                 "construction from coordinates cannot be evidence that "
                 "PIEZO1 bends its surroundings, which is what the published "
                 "panel is.",
          needs=("the unsharpened map from EMD-7042",)),
    Panel("4c", "4", "c", "Electrostatic surface at 150 mM NaCl",
          "Surface coloured by electrostatic potential computed with APBS, "
          "saturating at +-5 k_BT/e. Top, side and bottom.",
          "analogue", module="piezo1.physics.electrostatics",
          compute=_panels.panel_4c,
          reason="Ours is linear-superposition Debye-Huckel through a uniform "
                 "solvent dielectric, not a Poisson-Boltzmann solve with a "
                 "dielectric boundary. It reproduces the sign and the pattern; "
                 "it systematically under-estimates the magnitude, and on "
                 "6B3R nothing reaches the panel's +-5 saturation.",
          needs=("APBS, or a finite-difference PB solver with a dielectric "
                 "boundary",),
          published={"scale_kT_per_e": 5.0, "ionic_strength_mM": 150}),
    Panel("4-S1", "4-S1", "", "Interface between CED and TM loops",
          "The cap's acidic patch (E2257, E2258, D2264) against the loops' "
          "basic patch (R1761, R1762, R1269), connected by hydrogen bonds and "
          "salt bridges in a domain-swapped manner.",
          "replicated", module="piezo1.analysis.interactions",
          compute=_panels.panel_4_supplement,
          published={"pairs": ["E2257-R1762", "D2264-R1761"],
                     "domain_swapped": True}),

    # ---- Figure 5: proteoliposomes ----------------------------------------
    Panel("5a", "5", "a", "Piezo in a POPE:POPG vesicle",
          "A cryo-EM image of a single unilamellar vesicle containing one "
          "channel, with the molecular model scaled and inserted.",
          "not_replicable",
          reason="A micrograph of a reconstituted vesicle. The model half of "
                 "the panel is ours; the image it is inserted into is an "
                 "experiment this project cannot perform. The geometric claim "
                 "the panel makes — that the model conforms to the curvature "
                 "without adjustment — is what 7a measures.",
          needs=("cryo-EM images of proteoliposomes",)),
    Panel("5b", "5", "b", "Protein-free vesicles",
          "Vesicles of POPC:DOPS:cholesterol with no protein, spherical in "
          "projection — the control for panel c.",
          "not_replicable",
          reason="An experimental control: micrographs of vesicles made with "
                 "no protein in them, establishing that this lipid mixture "
                 "gives spherical vesicles on its own. It is the panel that "
                 "makes 5c an attribution rather than an observation, and it "
                 "contains no protein for a structure to reproduce.",
          needs=("cryo-EM images of protein-free vesicles",)),
    Panel("5c", "5", "c", "Piezo in a POPC:DOPS:cholesterol vesicle",
          "The same as 5a in a lipid composition whose vesicles are spherical "
          "without protein, so the local curvature is attributable to Piezo.",
          "not_replicable",
          reason="A micrograph of a reconstituted vesicle, in the lipid "
                 "composition whose empty vesicles are spherical. The "
                 "attribution of curvature to the channel rests on comparing "
                 "it with panel b, and neither image is something coordinates "
                 "can stand in for.",
          needs=("cryo-EM images of proteoliposomes",)),
    Panel("5-S1", "5-S1", "", "Trimer in liposomes of various sizes",
          "A gallery of vesicles across a size range.",
          "not_replicable",
          reason="A gallery of micrographs across a vesicle size range, "
                 "showing the deformation is not an artefact of one curvature. "
                 "Experimental images throughout.",
          needs=("cryo-EM images of proteoliposomes",)),

    # ---- Figure 6: the pore -----------------------------------------------
    Panel("6a", "6", "a", "Ion-conduction path",
          "The pore viewed from the side, the distance from the axis to the "
          "protein surface drawn as grey spheres, the pore-lining C-alpha "
          "trace in yellow and the pore-facing residues as sticks.",
          "replicated", module="piezo1.ui.pore_controller",
          compute=_mech.panel_6ab,
          published={"lining": "TM37, TM38, hairpin and PE helices"}),
    Panel("6b", "6", "b", "Pore radius profile",
          "van der Waals radius against displacement along the pore axis from "
          "the top, with the constricting residues labelled.",
          "replicated", module="piezo1.structure.pore",
          compute=_mech.panel_6ab,
          published={"radii_A": {"E2537": 0.1, "P2536": 0.4, "M2493": 0.3},
                     "software": "HOLE"}),
    Panel("6c", "6", "c", "Side-chain density at the constrictions",
          "A stereo view of the map around the constricting residues, "
          "contoured at 6 sigma.",
          "not_replicable",
          reason="Needs the sharpened map. A stereo pair of the model alone is "
                 "drawable and is not what the panel shows — its whole "
                 "purpose is that the side chains are supported by density.",
          needs=("EMD-7042",)),
    Panel("6-S1ab", "6-S1", "a,b", "Comparison with P2X and ASIC",
          "The pore region beside P2X and acid-sensing ion channels, "
          "establishing the shared trimeric two-TM-plus-extracellular-domain "
          "architecture.",
          "not_replicable",
          reason="Needs P2X and ASIC coordinates. This project's structure "
                 "registry is deliberately PIEZO-specific — the catalogue, the "
                 "numbering checks and the entity classifier all assume a "
                 "PIEZO — and admitting two unrelated channels to make one "
                 "panel would weaken every one of those guards. The comparison "
                 "is a real gap and is recorded as one rather than bodged.",
          needs=("PDB entries for P2X (e.g. 3H9V) and ASIC (e.g. 4NTW)",)),
    Panel("6-S1cd", "6-S1", "c,d", "The cuff: elbow, base, hairpin and PE",
          "The helices unique to Piezo that surround the narrowest part of the "
          "pore — elbow (2116-2142), base (2149-2175), hairpin (2501-2534) and "
          "the pore-extension helix.",
          "replicated", module="piezo1.core.annotations",
          compute=_mech.panel_6_supplement,
          published={"elbow": [2116, 2142], "base": [2149, 2175],
                     "hairpin": [2501, 2534], "numbering": "mouse"}),

    # ---- Figure 7: the dome model -----------------------------------------
    Panel("7a", "7", "a", "Trimer in the idealised dome",
          "A C-alpha trace inside a semi-spherical membrane 3.6 nm thick, "
          "mid-plane radius 10.2 nm, centred 4.0 nm above the projection "
          "plane.",
          "replicated", module="piezo1.physics.dome_idealised",
          compute=_mech.panel_7a, needs_structure=False,
          published={"radius_nm": 10.2, "center_height_nm": 4.0,
                     "thickness_nm": 3.6}),
    Panel("7b", "7", "b", "Beam and cross-helices in the idealised membrane",
          "The trimer with the beam (1300-1365) red and the cross-helices "
          "yellow, from the top and the side.",
          "replicated", module="piezo1.structure.architecture",
          compute=_mech.panel_7b,
          published={"beam": [1300, 1365], "numbering": "mouse"}),
    Panel("7c", "7", "c", "Projected area against flattening",
          "A schematic of the projected area growing as the dome flattens.",
          "replicated", module="piezo1.physics.dome_idealised",
          compute=_mech.panel_7c, needs_structure=False,
          published={"delta_area_nm2": 120}),
    Panel("7d", "7", "d", "Theoretical activation curves",
          "Open probability against tension for (dG_prot + dG_bend) of 20 or "
          "40 k_BT and dA_proj of 20 or 60 nm^2.",
          "replicated", module="piezo1.physics.dome",
          compute=_mech.panel_7d, needs_structure=False,
          published={"delta_g_kT": [20, 40], "delta_area_nm2": [20, 60]}),
    Panel("7-S1", "7-S1", "", "References for area and energy calculations",
          "The arithmetic behind Figure 7: a 400 nm^2 mid-plane surface, a "
          "280 nm^2 projected area, and a bending energy of about 150 k_BT.",
          "replicated", module="piezo1.physics.dome_idealised",
          compute=_mech.panel_7_supplement, needs_structure=False,
          published={"dome_area_nm2": 400, "projected_area_nm2": 280,
                     "delta_area_nm2": 120, "bending_energy_kT": 150}),
]


# --------------------------------------------------------------------------
# Access
# --------------------------------------------------------------------------

def panel_by_key(key: str) -> Panel:
    """One panel by its key, e.g. ``"7-S1"``. Raises with the valid keys."""
    for panel in PANELS:
        if panel.key == key:
            return panel
    raise KeyError(f"no panel {key!r}; known keys: "
                   + ", ".join(p.key for p in PANELS))


def not_replicable() -> list[Panel]:
    """Panels this project cannot reproduce, each with its reason."""
    return [p for p in PANELS if p.status == "not_replicable"]


def coverage() -> dict:
    """How much of the paper is reachable from coordinates."""
    counts = {status: sum(1 for p in PANELS if p.status == status)
              for status in STATUSES}
    figures: dict[str, dict] = {}
    for panel in PANELS:
        entry = figures.setdefault(panel.figure, {s: 0 for s in STATUSES})
        entry[panel.status] += 1
    return {
        "paper": PAPER,
        "n_panels": len(PANELS),
        "by_status": counts,
        "by_figure": figures,
        "replicable_fraction": (counts["replicated"] + counts["analogue"])
                                / len(PANELS),
        "needs": sorted({need for p in PANELS for need in p.needs}),
        "summary": (
            f"{counts['replicated']} of {len(PANELS)} panels reproduce from "
            f"deposited coordinates, {counts['analogue']} have an analogue "
            f"that is not the same quantity, and {counts['not_replicable']} "
            f"need experimental data this project does not hold."),
    }


def replicate(key: str, structure=None, reference: str | None = None,
              **kw) -> dict:
    """Reproduce one panel.

    ``structure`` is loaded from the paper's own entry when not given and the
    panel needs one. ``reference`` defaults to the entry's own numbering rather
    than to mouse, so running a panel against a human structure does not read
    it with mouse numbers — the failure mode this project has a whole module
    guarding against.
    """
    panel = panel_by_key(key)
    if panel.compute is None:
        return {"panel": panel.key, "label": panel.label,
                "status": panel.status, "title": panel.title,
                "shows": panel.shows, "reason": panel.reason,
                "needs": list(panel.needs), "published": dict(panel.published),
                "paper": PAPER}

    if panel.needs_structure and structure is None:
        from ..config import STRUCTURE_DIR
        from ..core.structure import Structure
        path = STRUCTURE_DIR / f"{DEFAULT_ENTRY}.cif"
        if not path.exists():
            raise FileNotFoundError(
                f"{DEFAULT_ENTRY} not downloaded — run "
                f"`python -m piezo1.io.fetch`")
        structure = Structure.from_file(path)

    if reference is None:
        reference = "mouse"
        if structure is not None:
            from ..io.registry import load_registry
            record = load_registry().get(structure.name)
            if record is not None:
                reference = record.numbering_species

    result = panel.compute(structure, reference=reference, **kw)
    return {"panel": panel.key, "label": panel.label, "status": panel.status,
            "title": panel.title, "shows": panel.shows,
            "module": panel.module, "reason": panel.reason,
            "published": dict(panel.published),
            "structure": (structure.name if structure is not None else None),
            "numbering": reference, "paper": PAPER, "result": result}


def replicate_all(structure=None, keys: list[str] | None = None,
                  skip_slow: bool = False) -> dict:
    """Run every replicable panel and collect the results.

    A panel that raises is recorded with its exception rather than aborting the
    run: one missing structure should not hide the twelve panels that need no
    structure at all.
    """
    slow = {"4c", "6a", "6b"}
    wanted = keys or [p.key for p in PANELS]
    out: dict[str, dict] = {}
    for key in wanted:
        panel = panel_by_key(key)
        if skip_slow and key in slow:
            out[key] = {"panel": key, "status": panel.status,
                        "skipped": "slow"}
            continue
        try:
            out[key] = replicate(key, structure=structure)
        except Exception as exc:                      # noqa: BLE001
            out[key] = {"panel": key, "status": panel.status,
                        "error": f"{type(exc).__name__}: {exc}"}
    return {"coverage": coverage(), "panels": out}
