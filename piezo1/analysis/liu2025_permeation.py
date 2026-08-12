"""Figure 5 as far as a continuum model can take it — and no further.

Liu et al. 2025 (Neuron 113:590–604) ran **coarse-grained molecular dynamics**
on a truncated pore module at four transmembrane potentials, counted Na⁺ into
each cavity, tracked individual permeation events, and read a slope conductance
of ~20 pS off the resulting current–voltage curve. This module runs the
project's existing one-dimensional drift-diffusion solver over the *same*
structures at the *same* four voltages.

**That is an analogue, not a replication, and the two differ in kind.** Their
simulation has explicit ions with trajectories; this has a concentration field.
The consequences are specific rather than decorative:

- Their **5C** counts Na⁺ that *accessed* each cavity over a microsecond,
  including ions that entered and went back. A one-dimensional steady state has
  a single flux everywhere along the pore by construction, so it cannot produce
  different counts for CV, EV, MV and IV — that panel is not reproducible here
  and is recorded as such rather than approximated.
- Their **5D** is a sampled cumulative count. Ours is the straight line
  ``N(t) = I·t/e`` that such a count fluctuates about, so the slopes are
  comparable and the scatter is not.
- Their **5E** is the panel this can genuinely be put beside: an I–V curve and
  the slope conductance read off it.
- Their **5F/5G** track ions through the lateral portals. Nothing here has a
  portal in it — see :mod:`piezo1.physics.conduction_path`, which can *exclude*
  the closed ends but does not model the opening. That is the panel the Martini
  scaffold exists for.

**The pathway matters more than the voltage.** On the axial route every entry
is refused, including 8IXO, because PIEZO1's axis is closed at both ends. These
sweeps therefore default to the ``lateral`` pathway — the route the paper
describes — and every result carries which pathway produced it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["VOLTAGES", "CAVITIES", "PUBLISHED", "VoltagePoint",
           "PermeationSweep", "sweep_voltages", "cumulative_permeation",
           "cavity_bounds", "slope_conductance"]

#: The four transmembrane potentials of their Figure 5A, volts. A transcription
#: of somebody else's protocol, not a tunable input of this project — changing
#: them would stop the comparison being with Figure 5.
VOLTAGES = (0.0, -0.1, -0.25, -0.5)

#: The four cavities of their Figure 2G and Figure 5C, each bounded by the two
#: residues *they* name. Given as ``(group, human number)`` pairs and resolved
#: through the curated group's own human/mouse detail, so a mouse entry is read
#: in mouse numbers and the offset is never assumed — it is 16 in the cap and
#: 26 at the inner helix.
#:
#: Bounded by single residues rather than by whole groups. Groups overlap along
#: the axis — ``pore_lining`` spans the entire inner helix and so straddles
#: ``hydrophobic_gate`` — and taking their extents gave a 1 A membrane
#: vestibule, which is a construction artefact rather than a cavity.
#:
#: Their definitions, mouse numbering: MV between I2466 and the TM gate V2476,
#: IV between K2479 and the intracellular constriction neck at E2537.
CAVITIES = (
    ("CV", ("cap_constriction", 2279), ("cap_gate", 2319)),
    ("EV", ("cap_gate", 2319), ("pore_lining", 2440)),
    ("MV", ("pore_lining", 2440), ("hydrophobic_gate", 2450)),
    ("IV", ("hydrophobic_gate", 2454), ("ctd_constriction", 2510)),
)

#: What the paper reports, for comparison — never used as an input.
PUBLISHED = {
    "slope_conductance_pS": 20.0,
    "measured_conductance_pS": 30.0,
    "current_at_-0.5V_pA": -9.0,
    "source": ("Liu et al. Neuron 2025;113:590-604 (PMID 39719701), "
               "Figure 5E; coarse-grained MD of the truncated pore module"),
}


@dataclass(frozen=True)
class VoltagePoint:
    """One solved voltage."""

    voltage: float                      # V
    current_pA: float
    conductance_pS: float
    blocked_by: str = ""

    @property
    def conducts(self) -> bool:
        return not self.blocked_by and self.current_pA != 0.0


@dataclass
class PermeationSweep:
    """An I–V curve over a chosen conduction pathway."""

    points: list[VoltagePoint] = field(default_factory=list)
    pathway: str = "lateral"
    structure: str = ""
    pathway_caveat: str = ""
    refused: str = ""

    @property
    def voltages(self) -> np.ndarray:
        return np.array([p.voltage for p in self.points], dtype=float)

    @property
    def currents_pA(self) -> np.ndarray:
        return np.array([p.current_pA for p in self.points], dtype=float)

    @property
    def conducts(self) -> bool:
        return any(p.conducts for p in self.points)

    def slope_pS(self) -> float:
        """Slope conductance from a least-squares line through the I–V points.

        Read the same way they read Figure 5E — off the curve — rather than
        taken from the solver's own conductance, so the two numbers are
        obtained by the same operation even though the physics underneath is
        not the same.
        """
        return slope_conductance(self.voltages, self.currents_pA)

    def comparison(self) -> str:
        if self.refused:
            return f"no sweep: {self.refused}"
        if not self.conducts:
            return (f"no conduction at any of the four voltages on the "
                    f"{self.pathway} pathway")
        slope = self.slope_pS()
        published = PUBLISHED["slope_conductance_pS"]
        return (f"{slope:.1f} pS slope against their {published:.0f} pS "
                f"({slope / published:.1f}x), measured "
                f"{PUBLISHED['measured_conductance_pS']:.0f} pS")

    def caveat(self) -> str:
        """Never omitted anywhere this is shown."""
        return ("ANALOGUE, not a reproduction: 1-D continuum drift-diffusion "
                "with no explicit ions, no water and no lateral portal. "
                + self.pathway_caveat)


def cavity_bounds(structure, profile) -> dict:
    """The z range of each cavity, from the curated residues that bound it.

    Returns ``{name: (z_low, z_high)}`` in Angstrom along the conduction axis,
    omitting any cavity whose bounding residues are not both resolved — a
    cavity inferred from one end would be a made-up volume.
    """
    from ..analysis.pore_regions import gate_numbering
    from ..core.annotations import load_annotations

    numbering = gate_numbering(structure)
    if numbering is None:
        return {}
    annotations = load_annotations(numbering)

    # Where the landmark *is*, not where it happens to touch the probe. Asking
    # for a lining hit lost the extracellular and membrane vestibules on 8IXO
    # — I2466 is resolved and simply does not contact the widest sphere that
    # fits beside it — and a cavity that vanishes because its boundary is roomy
    # is a bug, not a measurement.
    axis = profile.axis
    if axis is None:
        return {}
    ca = structure.mask_ca()
    z_of: dict[tuple, float] = {}
    for landmark in {m for _n, a, b in CAVITIES for m in (a, b)}:
        number = _resolve(annotations, landmark, numbering)
        if number is None:
            continue
        here = ca & (structure.res_seq == number)
        if here.any():
            z_of[landmark] = float(np.mean(
                axis.project(structure.xyz[here].astype(float))))

    out = {}
    for name, upper, lower in CAVITIES:
        if upper not in z_of or lower not in z_of:
            continue
        a, b = sorted((z_of[upper], z_of[lower]))
        if b > a:
            out[name] = (a, b)
    return out


def _resolve(annotations, landmark: tuple, numbering: str) -> int | None:
    """``(group, human number)`` -> the number this entry uses, or None.

    Goes through the curated group's own detail, which carries both species'
    numbers from a real alignment. Never an arithmetic offset: it is 16 at the
    cap gate and 26 at the inner helix, and both bound cavities here.
    """
    group_id, human = landmark
    group = annotations.group(group_id)
    if group is None:
        return None
    for detail in group.detail:
        if detail.get("human") == human:
            value = detail.get(numbering)
            return int(value) if value is not None else None
    return None


def sweep_voltages(structure, profile=None, pathway: str = "lateral",
                   voltages=VOLTAGES, grid=None) -> PermeationSweep:
    """Solve the pore at each voltage on the chosen conduction pathway.

    The wetting verdict is evaluated on the **truncated** profile, not the full
    axis — that is the whole point of choosing a pathway, and evaluating it on
    the full axis would refuse every entry however the path was cut.
    """
    from ..analysis.hydration import load_grid, predict_wetting
    from ..physics.conduction_path import conduction_path
    from ..physics.permeation import default_species, solve_pnp
    from ..structure.pore import pore_profile
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis

    name = getattr(structure, "name", "")
    if profile is None:
        blocks, _ = protomer_blocks(structure)
        profile = pore_profile(structure, detect_c3_axis(blocks))

    path = conduction_path(structure, profile, pathway)
    sweep = PermeationSweep(pathway=pathway, structure=name,
                            pathway_caveat=path.caveat())
    if path.refused:
        sweep.refused = path.refused
        return sweep

    grid = grid if grid is not None else load_grid()
    if not grid.available:
        sweep.refused = "CHAP grid not downloaded; run python -m piezo1.io.fetch"
        return sweep

    verdict = predict_wetting(structure, path.profile, grid=grid)
    species = default_species()
    for voltage in voltages:
        result = solve_pnp(path.profile, verdict, voltage=float(voltage),
                           species=species)
        sweep.points.append(VoltagePoint(
            voltage=float(voltage),
            current_pA=float(result.current) * 1e12,
            conductance_pS=float(result.conductance) * 1e12,
            blocked_by=str(result.blocked_by or "")))
    return sweep


def slope_conductance(voltages: np.ndarray, currents_pA: np.ndarray) -> float:
    """Least-squares slope of I against V, in picosiemens.

    Through the origin is *not* assumed: a pore with fixed charge has a
    reversal potential, and forcing the fit through zero would hide it.
    """
    voltages = np.asarray(voltages, dtype=float)
    currents = np.asarray(currents_pA, dtype=float)
    if len(voltages) < 2 or np.allclose(voltages, voltages[0]):
        return float("nan")
    slope = np.polyfit(voltages, currents * 1e-12, 1)[0]
    return float(slope * 1e12)


def cumulative_permeation(current_pA: float, times_ns: np.ndarray,
                          valence: int = 1) -> np.ndarray:
    """Their Figure 5D, as the line a Poisson count fluctuates about.

    ``N(t) = |I| t / (z e)``. A steady-state continuum current has no
    fluctuation in it, so this is the *expectation* their scatter is drawn
    from — comparable in slope, and silent about the variance, which is most of
    what their panel shows.
    """
    from ..physics.charge import ion_rate

    times = np.asarray(times_ns, dtype=float) * 1e-9
    return ion_rate(current_pA, valence) * times


def occupancy_is_not_available() -> str:
    """Why their Figure 5C has no analogue here. Recorded, not approximated.

    A one-dimensional steady state carries the same flux through every slice —
    that is what steady state means — so it cannot produce four different
    cavity counts. Their panel counts ions that *entered a cavity*, including
    the ones that turned round, which is a property of trajectories and not of
    a concentration field.
    """
    return ("Figure 5C counts Na+ that accessed each cavity, including ions "
            "that entered and did not proceed. A 1-D steady state has one flux "
            "through every slice by construction and cannot distinguish the "
            "cavities; reporting a number here would be an artefact of the "
            "discretisation, not a measurement.")
