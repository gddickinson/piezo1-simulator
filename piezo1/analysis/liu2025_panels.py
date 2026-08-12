"""The arithmetic behind Liu et al. 2025's figures, on our own coordinates.

Split from the registry at the seam :mod:`piezo1.analysis.guo2017_panels` uses:
:mod:`piezo1.analysis.liu2025` is the *judgement* about which panels this
project can reproduce and is readable without running anything, and this is the
measurement.

**All four of their states are deposited and all four are in our catalogue**,
which is what makes the comparison possible at all: PIEZO1-Curved (7WLT, also
5Z10 and 6B3R), S2472E-Curved (8IXN), S2472E-Intermediate (8IXO) and
PIEZO1-Flattened (7WLU). The S2472E-Flattened map was never modelled — the
paper says its resolution could not be improved — so any panel that needs it is
refused rather than approximated with one of the others.

Every distance here is measured between **C-alpha** atoms unless the function
says otherwise. The paper quotes side-chain distances for the cap-gate loops
and the spring linker, so those functions use side-chain atoms and say which,
and a structure that does not resolve them returns ``None`` rather than falling
back to C-alpha — a backbone distance reported where a side-chain distance was
asked for would be a smaller number that looks like agreement.
"""

from __future__ import annotations

import itertools

import numpy as np

__all__ = ["STATES", "load_state", "residue_atoms", "pair_distance",
           "axis_length", "v2476_diagonal", "spring_linker_span",
           "cap_gate_loop_span", "cavity_volumes", "curvature_radius",
           "pore_radius_profile", "state_displacement"]

#: Their four states and the entry each is deposited as. The S2472E-Flattened
#: map has no model, so it is deliberately absent.
STATES = {
    "PIEZO1-Curved": "7WLT",
    "S2472E-Curved": "8IXN",
    "S2472E-Intermediate": "8IXO",
    "PIEZO1-Flattened": "7WLU",
}

#: Landmarks the paper names, in **mouse** numbering (all four entries above
#: are mouse). Converted per structure through the curated annotation, never by
#: an offset — it is 16 in the cap and 26 at the inner helix.
LANDMARKS = {
    "cap_constriction": 2295,      # R2295, top of the cap
    "cap_gate": 2335,              # Y2335
    "cap_loop1_a": 2328,           # A2328
    "cap_loop1_d": 2326,           # D2326
    "cap_loop2_p": 2382,           # P2382
    "cap_loop2_e": 2383,           # E2383
    "spring_f": 2460,              # F2460
    "spring_y": 2464,              # Y2464
    "ih_top": 2466,                # I2466
    "tm_gate": 2476,               # V2476
    "ih_lower": 2479,              # K2479
    "neck": 2537,                  # E2537, intracellular constriction
    "blade_tip": 650,              # V650
}


def load_state(state: str):
    """Load one of their four states by name."""
    from ..config import STRUCTURE_DIR
    from ..core.structure import Structure

    if state not in STATES:
        raise ValueError(f"unknown state {state!r}; expected one of "
                         f"{', '.join(STATES)}")
    path = STRUCTURE_DIR / f"{STATES[state]}.cif"
    if not path.exists():
        raise FileNotFoundError(f"{STATES[state]} not downloaded — run "
                                f"`python -m piezo1.io.fetch`")
    return Structure.from_file(path)


def residue_atoms(structure, number: int, names: tuple[str, ...] | None = None
                  ) -> dict:
    """Per-chain coordinates of one residue, optionally restricted by atom name.

    Returns ``{chain: (n, 3) array}``. A chain that does not resolve the
    residue, or resolves it without the requested atoms, is simply absent —
    the caller decides whether a two-of-three answer is usable.
    """
    out = {}
    here = structure.res_seq == number
    if names is not None:
        here = here & np.isin(structure.atom_name, list(names))
    for chain in sorted(set(structure.chain[here].tolist())):
        sel = here & (structure.chain == chain)
        if sel.any():
            out[str(chain)] = structure.xyz[sel].astype(float)
    return out


def pair_distance(structure, a: int, b: int,
                  atoms_a: tuple[str, ...] | None = None,
                  atoms_b: tuple[str, ...] | None = None,
                  across_protomers: bool = False) -> float | None:
    """Distance between two residues, in Angstrom, or None if unresolvable.

    ``across_protomers`` measures every inter-chain pair and returns the
    smallest, which is what "the side-chain distance between neighbouring
    subunits" means; otherwise the two residues are taken within each chain and
    the mean over chains is returned.
    """
    first = residue_atoms(structure, a, atoms_a)
    second = residue_atoms(structure, b, atoms_b)
    if not first or not second:
        return None

    if across_protomers:
        best = []
        for ca, cb in itertools.product(first, second):
            if ca == cb:
                continue
            best.append(float(np.linalg.norm(
                first[ca].mean(axis=0) - second[cb].mean(axis=0))))
        return float(min(best)) if best else None

    shared = sorted(set(first) & set(second))
    if not shared:
        return None
    return float(np.mean([np.linalg.norm(first[c].mean(axis=0)
                                         - second[c].mean(axis=0))
                          for c in shared]))


def axis_length(structure) -> float | None:
    """Figure 2B: R2295 to E2537 along the pore axis, Angstrom.

    Their 110 A (curved) and 100 A (intermediate) are measured "along the
    central pore axis", so this projects onto the C3 axis rather than taking a
    straight-line distance — the two differ by however far off-axis the
    residues sit, which for E2537 is several Angstrom.
    """
    axis = _axis(structure)
    top = residue_atoms(structure, LANDMARKS["cap_constriction"], ("CA",))
    bottom = residue_atoms(structure, LANDMARKS["neck"], ("CA",))
    if not top or not bottom:
        return None
    z_top = np.mean([axis.project(v)[0] for v in top.values()])
    z_bottom = np.mean([axis.project(v)[0] for v in bottom.values()])
    return float(abs(z_top - z_bottom))


def v2476_diagonal(structure) -> float | None:
    """Figure 2E: the TM-gate side-chain diagonal, their 7 A -> 14 A.

    Side chains, not C-alpha: the whole claim is that the *side chains* move
    apart to clear the 9-12 A a hydrated Na+ needs. Valine has CG1 and CG2 and
    nothing further out, so those are the atoms.
    """
    per_chain = residue_atoms(structure, LANDMARKS["tm_gate"], ("CG1", "CG2"))
    if len(per_chain) < 3:
        return None
    centres = [v.mean(axis=0) for v in per_chain.values()]
    return float(np.mean([np.linalg.norm(a - b)
                          for a, b in itertools.combinations(centres, 2)]))


def spring_linker_span(structure) -> dict:
    """Figure 2F: the spring's F2460 and Y2464 side-chain separations.

    Their numbers are 9 A at F2460 in the curved state and 17 A at Y2464 in the
    intermediate — the residue facing the vestibule *changes* on compression,
    so both are reported for both states rather than one being picked to match.
    """
    out = {}
    for key, ring in (("F2460", ("CZ", "CE1", "CE2")),
                      ("Y2464", ("OH", "CZ"))):
        per_chain = residue_atoms(structure, LANDMARKS[
            "spring_f" if key.startswith("F") else "spring_y"], ring)
        if len(per_chain) < 3:
            out[key] = None
            continue
        centres = [v.mean(axis=0) for v in per_chain.values()]
        out[key] = float(np.mean([np.linalg.norm(a - b) for a, b
                                  in itertools.combinations(centres, 2)]))
    return out


def cap_gate_loop_span(structure) -> dict:
    """Figures 3F and 3H: the two cap-gate loop pairs, between neighbours.

    Published: A2328-P2382 widens from ~4.3 to ~16.2 A and D2326-E2383 from
    ~4.8 to ~12.8 A between PIEZO1-Curved and S2472E-Intermediate. Measured
    across protomers, because the loops that separate belong to different
    subunits — within one chain the distance is a different quantity entirely.
    """
    return {
        "A2328-P2382": pair_distance(
            structure, LANDMARKS["cap_loop1_a"], LANDMARKS["cap_loop2_p"],
            ("CB", "CA"), ("CG", "CB"), across_protomers=True),
        "D2326-E2383": pair_distance(
            structure, LANDMARKS["cap_loop1_d"], LANDMARKS["cap_loop2_e"],
            ("OD1", "OD2", "CG"), ("OE1", "OE2", "CD"), across_protomers=True),
    }


def pore_radius_profile(structure, pathway: str = "axial") -> dict:
    """Figure 2D: pore radius against position along the axis.

    Their x-axis starts at 2 A, so any constriction narrower than that is off
    the published plot — which matters, because ours are: R2295 pinches every
    entry to about 1 A and the panel cannot show it.
    """
    from ..physics.conduction_path import conduction_path

    profile = _profile(structure)
    path = conduction_path(structure, profile, pathway)
    used = path.profile
    return {"z": np.asarray(used.z, dtype=float),
            "radius": np.asarray(used.radius, dtype=float),
            "pathway": pathway,
            "bottleneck_A": float(np.min(used.radius)),
            "below_published_axis": int(np.sum(np.asarray(used.radius) < 2.0)),
            "caveat": path.caveat()}


def cavity_volumes(structure) -> dict:
    """Figure 2G: the volume of each cavity, cubic Angstrom.

    A solid of revolution about the measured pore path — ``sum(pi r^2 dz)``
    over the slices inside each cavity. That is **not** what they measured:
    theirs comes from the cavity-detection in their own analysis and follows
    the real, non-circular lumen. This is circular by construction, so it is an
    over-estimate wherever the lumen is not, and only the *ratios between
    states* carry information.
    """
    from .liu2025_permeation import cavity_bounds

    profile = _profile(structure)
    bounds = cavity_bounds(structure, profile)
    z = np.asarray(profile.z, dtype=float)
    r = np.asarray(profile.radius, dtype=float)
    step = float(np.median(np.diff(np.sort(z)))) if len(z) > 1 else 0.0
    out = {}
    for name, (lo, hi) in bounds.items():
        inside = (z >= lo) & (z <= hi)
        out[name] = float(np.pi * np.sum(r[inside] ** 2) * step)
    return out


def curvature_radius(structure) -> dict:
    """Figure 6: the mid-plane curvature radius R, nanometres.

    Published: ~10-12 nm for PIEZO1-Curved (7WLT, 5Z10, 6B3R), 14 nm for
    S2472E-Curved, 32 nm for S2472E-Intermediate and 117 nm for
    PIEZO1-Flattened. This is our own dome fit, which was built for Guo &
    MacKinnon's 10.2 nm and is here asked a question it was not tuned on: a
    structure four times flatter.
    """
    from ..structure.geometry import measure_dome, tm_surface_points
    from ..structure.protomers import protomer_blocks

    from .pore_regions import gate_numbering

    numbering = gate_numbering(structure) or "mouse"
    points, resolved = tm_surface_points(structure, numbering)
    blocks, _ = protomer_blocks(structure)
    if len(points) < 4 or not blocks:
        return {"radius_nm": float("nan"),
                "reason": "too few transmembrane helices resolved to fit a dome"}
    dome = measure_dome(blocks, points)
    return {"radius_nm": float(dome.radius_of_curvature) / 10.0,
            "depth_nm": float(dome.dome_depth) / 10.0,
            "n_helices": len(resolved),
            "rmse_A": float(getattr(dome.sphere, "rmsd", float("nan")))}


def state_displacement(state_a: str, state_b: str, residue: int) -> float | None:
    """How far one residue moves between two states, Angstrom.

    Both structures are put in their own **canonical frame** first — C3 axis on
    +z, cytosolic side down, centred — so the comparison is between two
    consistently framed models rather than between two deposition frames that
    sit 30-150 A apart. The paper aligns on a different element for each panel
    (the cap for Figure 3I, the CTD for Figure 3N), so a number from here will
    not match one of theirs exactly; the frame is stated so the difference is
    attributable rather than mysterious.
    """
    from ..structure.frame import standardise

    positions = []
    for state in (state_a, state_b):
        framed, _ = standardise(load_state(state), mode="canonical")
        per_chain = residue_atoms(framed, residue, ("CA",))
        if not per_chain:
            return None
        positions.append(np.mean([v.mean(axis=0)
                                  for v in per_chain.values()], axis=0))
    return float(np.linalg.norm(positions[0] - positions[1]))


def _axis(structure):
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis

    blocks, _ = protomer_blocks(structure)
    return detect_c3_axis(blocks)


def _profile(structure):
    from ..structure.pore import pore_profile

    return pore_profile(structure, _axis(structure))
