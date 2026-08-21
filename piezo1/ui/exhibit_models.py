"""The simulations an analysis window can hand the user.

A number in a table is one point on a curve, and which point it is depends on
inputs the reader cannot see. These are those curves: the model re-evaluated
while the user moves the input, so the shape of the dependence is visible
rather than asserted. Where the input is a **registered parameter**, the
control says so and starts at the registry's own default, marked on the plot.

Three rules, and the first is the one that would do damage if it were broken:

1. **Nothing here writes to the parameter registry.** Every model function in
   this project takes its inputs explicitly, so a slider passes a value per
   call. Overriding instead would leave the application quoting non-default
   numbers after the window closed — reports would carry the amber banner and
   ``verify_claims`` would refuse to run, for no reason the user could see.
2. **A moved slider is a sensitivity, not a measurement.** ``uncertainty.py``
   keeps those two apart deliberately; the window repeats it on every
   simulation, because a curve the user just produced looks exactly like one
   that was measured.
3. **Fast enough to be moved.** Everything here is closed-form or a small
   solve; the one that reads coordinates builds the pore profile once and then
   sweeps arithmetic over it.

Qt-free: the widget in :mod:`piezo1.ui.explore_window` builds sliders from
:class:`Control` and paints whatever :class:`~piezo1.ui.exhibits.ChartData`
comes back.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .exhibits import ChartData, Reference, Series, empty_chart

__all__ = ["Control", "Simulation", "Context", "SIMULATIONS", "run_simulation"]


@dataclass(frozen=True)
class Control:
    """One thing the user can move.

    ``parameter`` names the registered parameter this control *is*, when there
    is one — the default then comes from the registry rather than from a
    literal, and the panel can show the citation behind it. ``low`` and
    ``high`` are the slider's range and are deliberately narrower than the
    registry's bounds, which are set wide enough to admit any defensible value
    rather than to be swept.
    """

    key: str
    label: str
    low: float
    high: float
    default: float | None = None
    unit: str = ""
    parameter: str = ""
    log: bool = False
    #: Multiplies the registered value to reach the unit the slider is in —
    #: nanomolar where the registry holds molar, and so on. Applied here
    #: rather than baked into a literal default, so the starting position is
    #: resolved at call time and follows an override in the parameters dialog.
    scale: float = 1.0

    def start(self) -> float:
        if self.parameter:
            from ..parameters import PARAMETERS

            return float(PARAMETERS.value(self.parameter)) * self.scale
        if self.default is not None:
            return float(self.default)
        return float(self.low)


@dataclass
class Context:
    """What a simulation is allowed to read: the loaded structure and the
    result already in the window. Never the other way round — a simulation
    cannot change either."""

    structure: object | None = None
    result: dict = field(default_factory=dict)
    species: str = "human"
    _profile: object = None

    def pore_profile(self):
        """The measured radius profile, built once and kept.

        The same route ``analysis_permeation`` takes, so the bottleneck marked
        on the plot is the bottleneck in the table beside it.
        """
        if self._profile is None:
            from ..parameters import PARAMETERS
            from ..structure.pore import pore_profile
            from ..structure.protomers import protomer_blocks
            from ..structure.superpose import detect_c3_axis

            blocks, _ = protomer_blocks(self.structure)
            if not blocks:
                raise ValueError("needs three well-resolved protomers")
            axis = detect_c3_axis(blocks)
            self._profile = pore_profile(self.structure, axis,
                                         step=PARAMETERS.value("pore.step"))
        return self._profile


@dataclass(frozen=True)
class Simulation:
    key: str
    title: str
    what: str
    caveat: str
    controls: tuple[Control, ...]
    run: Callable[[dict, Context], ChartData]
    needs_structure: bool = False


def _p(key: str) -> float:
    from ..parameters import PARAMETERS

    return float(PARAMETERS.value(key))


# --------------------------------------------------------------------------
# Ion permeation: what the radius decides
# --------------------------------------------------------------------------

def _conductance_curve(values: dict, ctx: Context) -> ChartData:
    from ..physics.permeation import IonSpecies, series_conductance

    profile = ctx.pore_profile()
    z = np.asarray(profile.z, dtype=float) * 1e-10
    radius = np.asarray(profile.radius, dtype=float) * 1e-10

    scale = values["diffusion_scale"]
    cation_radius = values["ion_radius"]
    salt = values["bath"]
    species = [
        IonSpecies("K+", +1, _p("permeation.diffusion_cation") * scale,
                   cation_radius, salt),
        IonSpecies("Cl-", -1, _p("permeation.diffusion_anion") * scale,
                   _p("permeation.radius_anion"), salt),
    ]

    floors = np.linspace(max(cation_radius, 0.5), 8.0, 70)
    out = []
    for floor in floors:
        opened = np.maximum(radius, floor * 1e-10)
        result = series_conductance(z, opened, species=species)
        out.append(result["conductance"] * 1e12)

    return ChartData(
        title="Conductance against how far the constriction is opened",
        x_label="every slice opened to at least this radius (A)",
        y_label="conductance (pS)",
        series=[Series("series-resistor formula", list(floors), out,
                       color="#6fb1ff")],
        references=[
            Reference(_p("permeation.published_conductance"),
                      "published single-channel conductance"),
            Reference(float(profile.bottleneck_radius),
                      "this entry's own bottleneck", vertical=True,
                      color="#f2a65a"),
        ],
        note="Opening the constriction is a construction, not a state of the "
             "channel: the profile is the measured one with a floor applied. "
             "Two of the three controls — the in-pore diffusivity and the ion "
             "radius — are the inputs nothing has measured, and are the ones "
             "the window's own caveat says span 16-94 pS between them. The "
             "bath is a condition you would choose in a recording.")


# --------------------------------------------------------------------------
# HaloTag labelling
# --------------------------------------------------------------------------

def _labelling_curve(values: dict, _ctx: Context) -> ChartData:
    from ..analysis.labelling import (LabellingConditions, detectable_fraction,
                                      fully_labelled_fraction,
                                      site_labelled_fraction)

    conditions = LabellingConditions(
        concentration=values["concentration"] * 1e-9,
        k_on=values["k_on"],
        active_fraction=values["active_fraction"],
        name="explored")
    t = np.logspace(0, 4.5, 200)
    p = np.asarray(site_labelled_fraction(t, conditions), dtype=float)
    return ChartData(
        title="Labelling in time",
        x_label="incubation (s, log)", y_label="fraction",
        log_x=True,
        series=[
            Series("one site", list(t), list(p), color="#6fb1ff"),
            Series("all three sites", list(t),
                   list(np.asarray(fully_labelled_fraction(p), dtype=float)),
                   color="#f2a65a"),
            Series("at least one dye", list(t),
                   list(np.asarray(detectable_fraction(p), dtype=float)),
                   color="#7fd18a"),
        ],
        references=[
            Reference(_p("labelling.incubation_time"),
                      "the protocol's incubation", vertical=True),
            Reference(conditions.asymptote,
                      "ceiling set by the reactive fraction"),
        ],
        note="The ceiling is the reactive fraction cubed, and no incubation "
             "time beats it — which is the finding: a dye mixture needs "
             "unreactive tags, not a shorter incubation.")


# --------------------------------------------------------------------------
# The calcium nanodomain
# --------------------------------------------------------------------------

def _calcium_curve(values: dict, ctx: Context) -> ChartData:
    from ..physics.nanodomain import (calcium_at, distance_for_occupancy,
                                      saturation, screening_length)

    current = values["current_pA"] * 1e-12 * values["calcium_fraction"]
    length = screening_length(buffer_concentration=values["buffer_M"])
    kd = values["sensor_kd_uM"] * 1e-6

    r = np.logspace(-9.3, -6.0, 220)                    # 0.5 nm to 1 um
    molar = np.asarray(calcium_at(r, current, length=length), dtype=float)
    occupancy = np.asarray(saturation(molar, kd), dtype=float)
    nm = r * 1e9

    references = [Reference(kd * 1e6, "sensor K_D"),
                  Reference(_p("nanodomain.resting_calcium") * 1e6,
                            "resting calcium")]
    tag = (ctx.result or {}).get("tag_distance_nm")
    if tag:
        references.append(Reference(float(tag), "modelled tag distance",
                                    vertical=True, color="#f2a65a"))
    ninety = distance_for_occupancy(0.9, current, kd=kd)
    if np.isfinite(ninety):
        references.append(Reference(ninety * 1e9, "90% sensor occupancy",
                                    vertical=True, color="#7fd18a"))

    return ChartData(
        title="Free calcium around an open pore",
        x_label="distance from the pore exit (nm, log)",
        y_label="free calcium (uM, log)", log_x=True, log_y=True,
        series=[Series("free calcium", list(nm), list(molar * 1e6),
                       color="#6fb1ff"),
                Series("sensor occupancy", list(nm), list(occupancy),
                       color="#c678dd", axis=1)],
        references=references,
        note="A point source in free solution: the model does not know the "
             "protein is there. Only the product of buffer concentration and "
             "its on-rate enters, so the two are not separately identifiable "
             "from any nanodomain measurement.")


# --------------------------------------------------------------------------
# Modulators
# --------------------------------------------------------------------------

def _dose_response(values: dict, _ctx: Context) -> ChartData:
    from ..core.ligands import load_ligands

    hill = values["hill"]
    concentration = values["concentration_uM"]
    c = np.logspace(-4, 3, 240)
    series, labels = [], []
    for index, item in enumerate(load_ligands().ligands):
        if not item.potency:
            continue
        midpoint = float(item.potency["value"])
        if item.potency.get("unit") == "nM":
            midpoint /= 1000.0
        response = c ** hill / (c ** hill + midpoint ** hill)
        series.append(Series(f"{item.name} ({item.potency['measure']})",
                             list(c), list(response),
                             color=("#6fb1ff", "#f2a65a", "#7fd18a",
                                    "#c678dd")[index % 4]))
        labels.append(f"{item.name} {item.potency['measure']} "
                      f"{midpoint:g} uM ({item.potency['citation']})")
    if not series:
        return empty_chart("no modulator carries a measured potency")
    return ChartData(
        title="One-site curves through the measured potencies",
        x_label="concentration (uM, log)", y_label="fraction of maximum",
        log_x=True, series=series,
        references=[Reference(concentration, "chosen concentration",
                              vertical=True),
                    Reference(0.5, "half-maximal")],
        note="The midpoints are measured; the curves are assumed. The Hill "
             "coefficient is a slider because nothing here fitted one, and "
             "the three measures are not the same quantity — an EC50 for "
             "activation, an IC50 against Yoda1 and a bilayer K_D. "
             + "; ".join(labels))


# --------------------------------------------------------------------------
# The idealised dome: Guo & MacKinnon's Figure 7
# --------------------------------------------------------------------------

def _flattening(values: dict, _ctx: Context) -> ChartData:
    from ..physics.dome_idealised import flattening_series

    tension = values["tension"]
    delta_g_prot = values["delta_g_prot"]
    points = flattening_series()
    released = [p.delta_projected_nm2 for p in points]
    return ChartData(
        title="Flattening the idealised dome at constant membrane area",
        x_label="projected area released (nm^2)", y_label="energy (k_BT)",
        series=[
            Series("bending energy given up", released,
                   [p.delta_bending_kT for p in points], color="#6fb1ff"),
            Series("free energy of the transition", released,
                   [p.free_energy(tension, delta_g_prot) for p in points],
                   color="#f2a65a"),
        ],
        references=[Reference(0.0, "no change from the closed dome")],
        note="The paper's own idealisation — two lengths and closed-form "
             "spherical-cap geometry, with no structure in it. The protein "
             "term is a slider because the paper does not state it: its own "
             "numbers imply about +180 k_BT, which is what makes the "
             "transition cost anything at all.")


def _activation(values: dict, _ctx: Context) -> ChartData:
    from ..physics.dome import DomeModel

    model = DomeModel(delta_area=values["delta_area"],
                      delta_g0=values["delta_g0"])
    tension = np.linspace(0.0, 20.0, 240)
    return ChartData(
        title="Tension against open probability",
        x_label="membrane tension (mN/m)", y_label="open probability",
        series=[Series("two-state Boltzmann", list(tension),
                       list(np.asarray(model.open_probability(tension),
                                       dtype=float)), color="#6fb1ff")],
        references=[Reference(_p("kinetics.t50_measured"),
                              "measured half-activation", vertical=True,
                              color="#f2a65a"),
                    Reference(0.5, "half open")],
        note=(f"Half-activation here is {model.half_activation_mnm:.2f} mN/m. "
              f"The default parameters were fitted to reproduce a measured "
              f"tension response, so passing through the marked value is a "
              f"consistency check and not a prediction."))


# --------------------------------------------------------------------------

SIMULATIONS: dict[str, Simulation] = {
    "pore_conductance": Simulation(
        key="pore_conductance",
        title="What the radius decides",
        what="The measured profile with a radius floor applied, through the "
             "closed-form series-resistor formula.",
        caveat="A sensitivity, not a measurement: the floor is a construction "
               "and the two sliders are unmeasured inputs.",
        controls=(
            Control("diffusion_scale", "In-pore diffusivity, as a fraction of "
                    "bulk", 0.05, 1.0, parameter="permeation.diffusion_scale"),
            Control("ion_radius", "Cation radius", 0.5, 3.0, unit="A",
                    parameter="permeation.radius_cation"),
            Control("bath", "Bath concentration", 0.01, 1.0, unit="M",
                    parameter="permeation.bath_concentration"),
        ),
        run=_conductance_curve, needs_structure=True),

    "labelling_timecourse": Simulation(
        key="labelling_timecourse",
        title="Labelling in time",
        what="Per-site labelling and the whole-channel fractions that follow "
             "from it.",
        caveat="Kinetics imported from halotag_binding_sim; the reactive "
               "fraction and the linker are unverified assumptions.",
        controls=(
            Control("concentration", "Dye concentration", 1.0, 5000.0,
                    unit="nM", parameter="labelling.concentration",
                    scale=1e9, log=True),
            Control("k_on", "Labelling rate constant", 1e5, 1e7,
                    unit="1/(M s)", parameter="labelling.k_on", log=True),
            Control("active_fraction", "Fraction of tags that are reactive",
                    0.3, 1.0, parameter="labelling.active_fraction"),
        ),
        run=_labelling_curve),

    "calcium_profile": Simulation(
        key="calcium_profile",
        title="Calcium against distance",
        what="The buffered-diffusion Green's function, and the sensor "
             "occupancy it implies.",
        caveat="A point source in free solution, carrying a current borrowed "
               "from the one open-like entry in the catalogue.",
        controls=(
            Control("current_pA", "Unitary current", 0.2, 10.0, default=2.4,
                    unit="pA"),
            Control("calcium_fraction", "Calcium share of the current", 0.005,
                    0.3, parameter="nanodomain.calcium_current_fraction"),
            Control("buffer_M", "Mobile buffer", 1e-6, 1e-2,
                    parameter="nanodomain.buffer_concentration", log=True),
            Control("sensor_kd_uM", "Sensor K_D", 0.01, 100.0, unit="uM",
                    parameter="nanodomain.sensor_kd", scale=1e6, log=True),
        ),
        run=_calcium_curve),

    "dose_response": Simulation(
        key="dose_response",
        title="Dose-response from the measured potencies",
        what="A one-site Hill curve through each modulator's measured "
             "midpoint.",
        caveat="The midpoint is measured; the curve is assumed, and the three "
               "measures are not the same quantity.",
        controls=(
            Control("concentration_uM", "Concentration", 1e-3, 1e3,
                    default=10.0, unit="uM", log=True),
            Control("hill", "Hill coefficient (assumed)", 0.5, 3.0,
                    default=1.0),
        ),
        run=_dose_response),

    "dome_flattening": Simulation(
        key="dome_flattening",
        title="Figure 7c as arithmetic",
        what="The idealised dome flattened at constant membrane area.",
        caveat="The paper's idealisation. No structure enters it, and the "
               "protein term is not something the paper states.",
        controls=(
            Control("tension", "Membrane tension", 0.0, 4.0, default=1.0,
                    unit="k_BT/nm^2"),
            Control("delta_g_prot", "Protein term dG_prot", 0.0, 300.0,
                    default=180.0, unit="k_BT"),
        ),
        run=_flattening),

    "dome_activation": Simulation(
        key="dome_activation",
        title="Figure 7d: tension against open probability",
        what="The two-state Boltzmann the dome mechanism rests on.",
        caveat="The defaults were fitted to a measured tension response, so "
               "agreement with the marked half-activation is not a "
               "prediction.",
        controls=(
            Control("delta_area", "Area change on gating", 1.0, 60.0,
                    unit="nm^2", parameter="dome.delta_area"),
            Control("delta_g0", "Intrinsic closed-state bias", 0.0, 30.0,
                    unit="k_BT", parameter="dome.delta_g0"),
        ),
        run=_activation),
}


def run_simulation(key: str, values: dict, context: Context) -> ChartData:
    """Run one simulation, and never raise into the window.

    Same rule as the chart dispatcher: a model that cannot be evaluated at the
    chosen inputs is a panel saying so. Some of these have genuinely invalid
    corners — a pore floor below the ion radius, a buffer that screens to
    nothing — and the user is entitled to walk into one.
    """
    simulation = SIMULATIONS.get(key)
    if simulation is None:
        return empty_chart(f"no simulation named {key!r}")
    if simulation.needs_structure and context.structure is None:
        return empty_chart("this simulation needs the loaded structure's own "
                           "coordinates; load a structure first")
    try:
        return simulation.run(dict(values), context)
    except Exception as exc:                      # noqa: BLE001 — see docstring
        return empty_chart(f"could not run this model: "
                           f"{type(exc).__name__}: {exc}")
