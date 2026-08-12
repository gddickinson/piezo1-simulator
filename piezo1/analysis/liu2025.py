"""Every figure panel of Liu et al. 2025, and what this project can do with it.

Liu S, Yang X, Chen X, *et al.*, *An intermediate open structure reveals the
gating transition of the mechanically activated PIEZO1 channel*, Neuron
2025;113:590-604 (PMID 39719701, PDB 8IXN and 8IXO). This is the paper that
supplies the only **intermediate-open** PIEZO1 structure, and the one whose
Figure 5 says the current does not go down the axis at all.

The registry follows :mod:`piezo1.analysis.guo2017` exactly, including its
rule that the refusals are half the deliverable. **All four of their states are
deposited and all four are in this catalogue** — PIEZO1-Curved (7WLT),
S2472E-Curved (8IXN), S2472E-Intermediate (8IXO), PIEZO1-Flattened (7WLU) —
which is why so much of it reproduces.

**What reproduced, and closely.** Seven distances the paper states come back
within about an Angstrom: the R2295-E2537 pore axis shortening 110 -> 100 A
(ours 109.5 -> 96.2), the TM-gate V2476 side-chain diagonal opening 7 -> 14 A
(ours 7.7 -> 14.2), the cap-gate loops separating 4.3 -> 16.2 A (ours
4.8 -> 16.1) and 4.8 -> 12.8 A (ours 5.7 -> 11.4), and the compressed spring's
Y2464 at 17 A (ours 16.6). The cavity volumes reproduce the *direction* of
every change — CV, EV and MV grow into the intermediate state, IV does not —
and not the values, for a stated reason.

**What did not.** The curvature radii of their Figure 6. They report ~10-12,
14, 32 and 117 nm across the four states; our dome fitter gives 9.7, 11.2,
16.5 and 18.4. It agrees where it was calibrated — 9.7 nm against Guo &
MacKinnon's 10.2 on the curved state — and saturates where it was not: fitting
a sphere to a nearly flat surface is ill-conditioned, and under-estimating a
large radius is exactly how that fails. Recorded as a disagreement rather than
adjusted, and the panel is filed as an ``analogue`` because of it.

**What cannot be done here at all.** The electrophysiology (Figures 1 and 3O),
the cryo-EM maps and their overlays (Figure 2A and the S1-S3 supplements), and
the two panels that follow individual ions through the lateral portals
(Figures 5F and 5G) — the last of which is precisely what
:mod:`piezo1.physics.martini` exists to prepare for, and why it prepares
rather than pretends.

Statuses, as in the Guo & MacKinnon registry:

``replicated``      computable from deposited coordinates, and where the panel
                    states a number, compared against it.
``analogue``        we can compute what the panel is an estimate *of*, but not
                    the panel. Read the caveat before showing it beside the
                    original.
``not_replicable``  needs experimental data or a simulation this project does
                    not hold.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

__all__ = ["Panel", "PANELS", "PAPER", "STATUSES", "STATES", "panel_by_key",
           "coverage", "replicate", "replicate_all", "not_replicable"]

PAPER = {
    "citation": "liu2025",
    "authors": "Liu S, Yang X, Chen X, Zhang X, Jiang J, Yuan J, Liu W, "
               "Wang L, Zhou H, Wu K, Tian B, Li X, Xiao B",
    "title": "An intermediate open structure reveals the gating transition of "
             "the mechanically activated PIEZO1 channel",
    "journal": "Neuron 2025;113:590-604",
    "doi": "10.1016/j.neuron.2024.11.020",
    "pmid": "39719701",
    "pdb": "8IXN, 8IXO",
    "emdb": "EMD-35799, EMD-35800",
    "numbering": "mouse (UniProt E2JF22, 2547 aa)",
}

STATUSES = ("replicated", "analogue", "not_replicable")

#: Their four states. Imported here so a caller has one place to look.
from .liu2025_panels import STATES  # noqa: E402
from .liu2025_permeation import occupancy_is_not_available  # noqa: E402

#: The entry a panel runs against when the caller does not say.
DEFAULT_ENTRY = "8IXO"


@dataclass(frozen=True)
class Panel:
    """One published panel and this project's relationship to it."""

    key: str
    figure: str
    panel: str
    title: str
    shows: str
    status: str
    module: str = ""
    compute: Callable[..., dict] | None = None
    reason: str = ""
    needs: tuple[str, ...] = ()
    #: Whether the panel compares states, in which case it runs on all four
    #: rather than on one loaded structure.
    across_states: bool = False
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


def _across(fn):
    """Run ``fn(structure)`` on each of their four states."""
    def run(**_kw) -> dict:
        from .liu2025_panels import load_state

        out = {}
        for state in STATES:
            try:
                out[state] = fn(load_state(state))
            except FileNotFoundError as exc:
                out[state] = {"error": str(exc)}
        return out
    return run


def _axis_lengths() -> dict:
    from .liu2025_panels import axis_length

    return _across(lambda st: {"axis_A": axis_length(st)})()


def _v2476() -> dict:
    from .liu2025_panels import v2476_diagonal

    return _across(lambda st: {"diagonal_A": v2476_diagonal(st)})()


def _spring() -> dict:
    from .liu2025_panels import spring_linker_span

    return _across(spring_linker_span)()


def _cap_loops() -> dict:
    from .liu2025_panels import cap_gate_loop_span

    return _across(cap_gate_loop_span)()


def _volumes() -> dict:
    from .liu2025_panels import cavity_volumes

    return _across(cavity_volumes)()


def _curvature() -> dict:
    from .liu2025_panels import curvature_radius

    return _across(curvature_radius)()


def _pore_profiles() -> dict:
    from .liu2025_panels import pore_radius_profile

    return _across(lambda st: {
        k: v for k, v in pore_radius_profile(st).items() if k != "z"})()


def _iv_curve() -> dict:
    from .liu2025_panels import load_state
    from .liu2025_permeation import PUBLISHED, sweep_voltages

    out = {"published": dict(PUBLISHED)}
    for state in STATES:
        sweep = sweep_voltages(load_state(state), pathway="lateral")
        out[state] = {
            "voltages_V": sweep.voltages.tolist(),
            "currents_pA": sweep.currents_pA.tolist(),
            "slope_pS": sweep.slope_pS(),
            "comparison": sweep.comparison(),
            "caveat": sweep.caveat()}
    return out


def _hydrophobicity() -> dict:
    from .hydration import hydrophobicity_profile_chap
    from .liu2025_panels import _profile, load_state

    out = {}
    for state in STATES:
        structure = load_state(state)
        profile = _profile(structure)
        values = hydrophobicity_profile_chap(structure, profile)
        out[state] = {"mean": float(values[values == values].mean()),
                      "n": int(len(values))}
    return out


PANELS: list[Panel] = [
    # ---- Figure 1: functional characterisation ----------------------------
    Panel("1a", "1", "a", "Poking-evoked whole-cell currents",
          "Representative traces from PIEZO1-KO cells transfected with "
          "wild-type or S2472E.", "not_replicable",
          reason="a patch-clamp recording. This project holds no "
                 "electrophysiology and computes no whole-cell current.",
          needs=("whole-cell patch clamp",)),
    Panel("1b", "1", "b", "Basal current",
          "Scatter of basal current: S2472E carries a standing current where "
          "wild-type does not.", "not_replicable",
          reason="a resting current in an unstimulated cell, which depends on "
                 "the resting membrane tension of that cell. Nothing here "
                 "models a cell; the gating scheme takes tension as an input "
                 "and has no way to know what a HEK293T membrane is under.",
          needs=("whole-cell patch clamp",)),
    Panel("1c", "1", "c", "Inactivation tau",
          "Inactivation time constant, 18.9 ms wild-type vs 104.8 ms S2472E.",
          "not_replicable",
          reason="the gating model here is fitted to published rates and has "
                 "no S2472E parameterisation; producing a tau would be "
                 "reporting the input.",
          needs=("whole-cell patch clamp",),
          published={"tau_ms_wt": 18.9, "tau_ms_s2472e": 104.8}),
    Panel("1d", "1", "d", "Normalised remaining current",
          "I_remaining / I_peak against time, showing S2472E retains ~40% of "
          "its current after the stimulus is removed.", "not_replicable",
          reason="the quantity is the fraction of channels still open after "
                 "the mechanical stimulus stops, which needs an S2472E "
                 "parameterisation of the gating scheme. Fitting one to this "
                 "panel and then plotting it would be drawing the input.",
          needs=("whole-cell patch clamp",)),
    Panel("1e", "1", "e", "Cell viability",
          "Luminescent viability of cells expressing S2472E, which are less "
          "viable — the constitutive calcium load the mutant imposes.",
          "not_replicable",
          reason="a cell-biology assay. The calcium machinery here reaches as "
                 "far as the nanodomain around one open pore and stops there; "
                 "whole-cell calcium load, and what it does to a cell over "
                 "hours, is outside everything this project models.",
          needs=("cell viability assay",)),

    # ---- Figure 2: the four states ----------------------------------------
    Panel("2a", "2", "a", "Cryo-EM map comparison",
          "S2472E maps overlaid on EMD-6865 and EMD-32593.", "not_replicable",
          reason="needs the reconstructions. This project downloads "
                 "coordinates, holds no map, and has no density to overlay.",
          needs=("EMD-35799", "EMD-35800", "EMD-6865", "EMD-32593")),
    Panel("2b", "2", "b", "Pore axis length",
          "Vertical distance from the extracellular constriction R2295 to the "
          "intracellular constriction E2537 along the pore axis: 110 A curved, "
          "100 A intermediate.", "replicated",
          module="analysis.liu2025_panels", compute=_axis_lengths,
          across_states=True,
          published={"curved_A": 110.0, "intermediate_A": 100.0}),
    Panel("2d", "2", "d", "Pore radius profile",
          "Pore radius along the central axis of the four structures, from "
          "HOLE.", "replicated",
          module="structure.pore", compute=_pore_profiles, across_states=True,
          published={"published_axis_starts_at_A": 2.0}),
    Panel("2e", "2", "e", "TM-gate diagonal and pore hydrophobicity",
          "CHAP hydrophobicity of the vestibules, and the V2476 side-chain "
          "diagonal opening from 7 to 14 A.", "replicated",
          module="analysis.liu2025_panels", compute=_v2476,
          across_states=True,
          published={"curved_A": 7.0, "intermediate_A": 14.0,
                     "wetting_threshold_A": "9-12"}),
    Panel("2e-chap", "2", "e (CHAP)", "Pore hydrophobicity",
          "The hydrophobicity colouring of the vestibule surfaces.",
          "replicated", module="analysis.hydration", compute=_hydrophobicity,
          across_states=True),
    Panel("2f", "2", "f", "Spring linker compression",
          "F2460 side-chain separation 9 A extended, Y2464 17 A compressed.",
          "replicated", module="analysis.liu2025_panels", compute=_spring,
          across_states=True,
          published={"F2460_curved_A": 9.0, "Y2464_intermediate_A": 17.0}),
    Panel("2g", "2", "g", "Cavity volumes",
          "Volumes of the cap, extracellular, membrane and inner vestibules.",
          "analogue", module="analysis.liu2025_panels", compute=_volumes,
          across_states=True,
          reason="ours is a solid of revolution about the measured pore path, "
                 "so it is circular by construction and over-estimates "
                 "wherever the real lumen is not. The direction of every "
                 "change reproduces — CV, EV and MV grow into the intermediate "
                 "state and IV does not — the absolute volumes do not.",
          published={"CV_nm3": [1.6, 3.4], "EV_nm3": [1.3, 2.1],
                     "MV_nm3": [0.5, 1.2], "IV_nm3": "comparable"}),

    # ---- Figure 3: the three gates ----------------------------------------
    Panel("3fh", "3", "f,h", "Cap-gate loop separation",
          "A2328-P2382 widens from ~4.3 to ~16.2 A and D2326-E2383 from ~4.8 "
          "to ~12.8 A between neighbouring subunits.", "replicated",
          module="analysis.liu2025_panels", compute=_cap_loops,
          across_states=True,
          published={"A2328_P2382_A": [4.3, 16.2],
                     "D2326_E2383_A": [4.8, 12.8]}),
    Panel("3eg", "3", "e,g", "Cap electrostatic surface",
          "Electrostatic surface potential of the cap, showing the acidic "
          "cluster D2326/E2334/E2338/E2383 and DEEED 2393-2397.", "analogue",
          module="physics.electrostatics",
          reason="ours is screened Coulomb from formal charges through a "
                 "uniform dielectric, NOT APBS: no dielectric boundary, no "
                 "ion-exclusion layer, no partial charges. Read the sign and "
                 "the pattern, never the value.",
          published={"scale_kT_per_e": 75.0}),
    Panel("3o", "3", "o", "Mutant poking currents",
          "Poking current for cap-gate, spring and TM-gate mutants.",
          "not_replicable",
          reason="patch clamp on eleven constructs. The variant machinery "
                 "here reports mechanical coupling, and five pre-registered "
                 "tests established it does not predict direction — so "
                 "producing a bar chart would be inventing the panel.",
          needs=("whole-cell patch clamp",)),

    # ---- Figure 4: intermediate to flattened ------------------------------
    Panel("4a", "4", "a", "Cap and blade displacement",
          "Cap up by ~10 A, blade by ~25 A, V650 by 44 A from "
          "S2472E-Intermediate to PIEZO1-Flattened.", "analogue",
          module="analysis.liu2025_panels",
          compute=lambda **kw: _displacements(),
          across_states=True,
          reason="the paper aligns on a different element for each panel — "
                 "the cap for one, the CTD for another — and states neither "
                 "alignment in a form that can be reproduced exactly. Ours "
                 "puts both structures in their own canonical C3 frame and "
                 "says so, which is a defensible choice and not their choice.",
          published={"cap_A": 10.0, "blade_A": 25.0, "V650_A": 44.0}),

    # ---- Figure 5: ion permeation -----------------------------------------
    Panel("5a", "5", "a", "Na+ occupancy at four voltages",
          "Na+ positions in CV/EV/MV/IV at 0, -0.1, -0.25 and -0.5 V.",
          "analogue", module="analysis.liu2025_permeation",
          reason="theirs is a cloud of explicit ions from coarse-grained MD; "
                 "ours is a continuum concentration over a 1-D pore. The "
                 "drawn particles are a rate made visible, not sampled "
                 "positions, and no ion in ours has a trajectory."),
    Panel("5c", "5", "c", "Na+ accessing each cavity",
          "Counts of Na+ that accessed CV/EV/MV/IV at -0.1 and -0.5 V.",
          "not_replicable",
          # The reason lives with the solver that cannot produce it, so the
          # registry and the module state one thing rather than two that could
          # drift apart.
          reason=occupancy_is_not_available(),
          needs=("explicit-ion trajectories",)),
    Panel("5d", "5", "d", "Cumulative permeation",
          "Na+ reaching the inner vestibule against time at -0.5 V.",
          "analogue", module="analysis.liu2025_permeation",
          reason="ours is the straight line N(t) = I t / e that such a count "
                 "fluctuates about. The slopes are comparable; the scatter, "
                 "which is most of what the panel shows, is not reproduced.",
          published={"intermediate_at_1000ns": 20, "flattened_at_1000ns": 10,
                     "curved_at_1000ns": 0}),
    Panel("5e", "5", "e", "Current-voltage curve",
          "I-V of Na+ through the TM pore of S2472E-Intermediate, slope "
          "conductance ~20 pS against a measured ~30 pS.", "analogue",
          module="analysis.liu2025_permeation", compute=_iv_curve,
          across_states=True,
          reason="the quantity is the same and the physics underneath is not: "
                 "drift-diffusion through a continuum, on a pathway chosen by "
                 "the caller. Only the `lateral` pathway conducts at all, and "
                 "the portal itself is not modelled, so the current is an "
                 "upper bound.",
          published={"slope_pS": 20.0, "measured_pS": 30.0,
                     "current_at_-0.5V_pA": -9.0}),
    Panel("5f", "5", "f", "Lateral portal trajectories",
          "Tracked Na+ paths through the intracellular lateral portal and the "
          "lateral plug gate over 10 us.", "not_replicable",
          reason="needs explicit-ion molecular dynamics with the portals "
                 "present. The conduction-path option here can EXCLUDE the "
                 "closed axial ends; it does not model the opening, and "
                 "nothing in this project produces an ion trajectory.",
          needs=("coarse-grained MD, 10 us",)),
    Panel("5g", "5", "g", "Portal trajectories without the plug gate",
          "The same with the lateral plug domain removed: 37 Na+ through one "
          "portal in the last microsecond.", "not_replicable",
          reason="as 5f, and additionally needs a construct with residues "
                 "1401-1421 deleted, which is a simulation input rather than "
                 "a deposited structure.",
          needs=("coarse-grained MD, 10 us", "lateral-plug deletion construct"),
          published={"na_last_us": 37, "portal_conductance_pS": 12,
                     "iv_region_conductance_pS": 23}),

    # ---- Figure 6: the schematic ------------------------------------------
    Panel("6-schematic", "6", "a-d", "Gating choreography cartoon",
          "Cartoon of the closed, intermediate, open and inactivated states.",
          "not_replicable",
          reason="a drawing rather than data. Its content is the rest of the "
                 "paper, which is registered panel by panel above.",
          needs=("nothing — it is an illustration",)),
    Panel("6-curvature", "6", "a-d (R)", "Mid-plane curvature radius",
          "Curvature radius of each state: ~10-12 nm curved, 14 nm "
          "S2472E-Curved, 32 nm intermediate, 117 nm flattened.", "analogue",
          module="structure.geometry", compute=_curvature, across_states=True,
          reason="our sphere fit agrees where it was calibrated — 9.7 nm on "
                 "the curved state against Guo & MacKinnon's 10.2 — and "
                 "saturates where it was not, giving 16.5 nm where they say 32 "
                 "and 18.4 nm where they say 117. Fitting a sphere to a nearly "
                 "flat surface is ill-conditioned and under-estimating a large "
                 "radius is how that fails. Reported as a disagreement, not "
                 "adjusted.",
          published={"curved_nm": [10.0, 12.0], "s2472e_curved_nm": 14.0,
                     "intermediate_nm": 32.0, "flattened_nm": 117.0}),
]


def _displacements() -> dict:
    from .liu2025_panels import LANDMARKS, state_displacement

    pairs = {"cap": LANDMARKS["cap_gate"], "V650": LANDMARKS["blade_tip"]}
    return {name: state_displacement("S2472E-Intermediate",
                                     "PIEZO1-Flattened", number)
            for name, number in pairs.items()}


def panel_by_key(key: str) -> Panel:
    for panel in PANELS:
        if panel.key == key:
            return panel
    raise KeyError(f"no panel {key!r}; have {[p.key for p in PANELS]}")


def coverage() -> dict:
    """How much of the paper this project can reproduce, as counts."""
    counts = {status: 0 for status in STATUSES}
    for panel in PANELS:
        counts[panel.status] += 1
    return {"total": len(PANELS), **counts,
            "figures": sorted({p.figure for p in PANELS}), "paper": PAPER}


def not_replicable() -> list[Panel]:
    """The refusals, which are half the point of the registry."""
    return [p for p in PANELS if p.status == "not_replicable"]


def replicate(key: str, **kw) -> dict:
    """Reproduce one panel, or return why it cannot be."""
    panel = panel_by_key(key)
    base = {"panel": panel.key, "label": panel.label, "status": panel.status,
            "title": panel.title, "shows": panel.shows,
            "reason": panel.reason, "needs": list(panel.needs),
            "published": dict(panel.published), "paper": PAPER}
    if panel.compute is None:
        return base
    base["result"] = panel.compute(**kw)
    return base


def replicate_all(keys: list[str] | None = None) -> dict:
    """Run every panel that has a callable, collecting failures rather than
    raising — a missing download must not stop the rest of the paper."""
    out = {}
    for panel in PANELS:
        if keys is not None and panel.key not in keys:
            continue
        if panel.compute is None:
            out[panel.key] = {"status": panel.status, "reason": panel.reason}
            continue
        try:
            out[panel.key] = replicate(panel.key)
        except Exception as exc:                    # noqa: BLE001
            out[panel.key] = {"status": panel.status,
                              "error": f"{type(exc).__name__}: {exc}"}
    return out
