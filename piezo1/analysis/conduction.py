"""Whether a pore conducts, with each criterion read off the right profile.

**The rule, stated before the numbers it produces.**

    The hydrophobic-gate criterion is evaluated on the COMPLETE AXIAL PROFILE.
    Rao et al.'s 0.55 cutoff is a sum of distances over a whole channel's
    lining, so it is not comparable between profiles of different length, and
    applying it to a truncated one is invalid on its own terms.

    The steric criterion is evaluated on the CONDUCTION PATH. Whether an ion
    physically fits is a property of the route it takes, and PIEZO1's route
    does not use the closed ends of its own axis.

    A pore conducts when neither criterion refuses it.

That rule is written here first because the thing it replaced was arrived at
badly. Round 84d evaluated *both* halves on whichever profile the pathway
produced. On the axial route that refuses every entry — the cap apex and the
cytoplasmic neck are shut in all of them — and on a lateral route it refuses
almost nothing, because truncating the path also truncates the sum:

===========================  ============  ============
score                        axial         lateral
===========================  ============  ============
7WLT                         1.35          0.13
6B3R                         2.05          0.22
entries above the cutoff     13 of 18      **0 of 18**
===========================  ============  ============

With the chemistry disabled, the verdict rested entirely on whether a residual
radius cleared the water radius — and that residual was the **cap gate** in 14
of 18 entries and the transmembrane gate in none. The Round 84d lateral
conductances are superseded by this module and should not be quoted; they were
that residual, not a gate.

**How this rule was reached, recorded because it matters.** Three compositions
were tried and the third separated the states. That is the wrong order, and
the rule is adopted on its own argument rather than on that agreement: the
cutoff's calibration is the reason the chemistry stays on the full axis, and it
would be the reason whichever way the numbers came out. The agreement is
reported below as a *check*, not as the derivation.

**What the check found.** Across the catalogue, on the lateral route: 1 of 15
curved entries conducts (3JAC, the 4.8 A entry with 346 unnamed residues), and
all three non-curved states do. The ordering of the three matches Liu et al.'s
Figure 5D, where curved passes ~0 Na+ in a microsecond, flattened ~10 and
intermediate ~20 — including the part that would be easy to get wrong, that the
*flattened* state conducts at all. It does, in their data and in ours, because
its transmembrane gate is dilated while its cap gate is shut.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["RULE", "ConductionVerdict", "conduction_verdict"]

RULE = (
    "hydrophobic gate from the complete axial profile, because Rao et al.'s "
    "cutoff is a sum over a whole channel's lining and is not comparable "
    "between profiles of different length; steric occlusion from the "
    "conduction path, because whether an ion fits is a property of the route"
)


@dataclass
class ConductionVerdict:
    """The composed verdict, carrying which profile decided each half."""

    conductive: bool
    #: The full-axis wetting prediction. Its `score` is the published quantity.
    wetting: object = None
    #: The chosen route.
    path: object = None
    #: Narrowest radius on the route, Angstrom — not on the full axis.
    steric_radius: float = float("nan")
    reasons: list = field(default_factory=list)

    @property
    def hydrophobic_gate(self) -> bool:
        return bool(getattr(self.wetting, "hydrophobic_gate", False))

    @property
    def sterically_occluded(self) -> bool:
        from .hydration import WATER_RADIUS_NM

        return bool(np.isfinite(self.steric_radius)
                    and self.steric_radius / 10.0 < WATER_RADIUS_NM)

    def summary(self) -> str:
        """Both halves, and which profile each came from. Never abbreviated."""
        wetting = self.wetting
        score = getattr(wetting, "score", float("nan"))
        text = (f"score {score:.2f} on the full axis, narrowest "
                f"{self.steric_radius / 10:.3f} nm on the "
                f"{getattr(self.path, 'pathway', '?')} route -> "
                f"{'conductive' if self.conductive else 'non-conductive'}")
        if self.reasons:
            text += " (" + " + ".join(self.reasons) + ")"
        return text

    def caveat(self) -> str:
        return RULE + (" · " + self.path.caveat() if self.path is not None
                       else "")


def conduction_verdict(structure, profile, pathway: str = "axial", grid=None
                       ) -> ConductionVerdict:
    """Apply the rule above to one structure.

    ``pathway="axial"`` makes the two profiles the same object, so the answer
    is **bit-identical** to reading `predict_wetting(...).conductive` directly
    — every number recorded before this module existed is reproduced. A test
    asserts that on every downloaded entry rather than on one.
    """
    from ..physics.conduction_path import conduction_path
    from .hydration import load_grid, predict_wetting

    grid = grid if grid is not None else load_grid()
    path = conduction_path(structure, profile, pathway)

    # Chemistry: always the complete axis, whatever route was chosen.
    wetting = predict_wetting(structure, profile, grid=grid)
    # Sterics: the route, which on the axial pathway is the same profile.
    radius = np.asarray(path.profile.radius, dtype=float)
    steric = float(radius.min()) if len(radius) else float("nan")

    verdict = ConductionVerdict(conductive=False, wetting=wetting, path=path,
                                steric_radius=steric)
    if not getattr(wetting, "available", False):
        verdict.reasons.append(
            getattr(wetting, "verdict", "wetting unavailable"))
        return verdict

    if verdict.hydrophobic_gate:
        verdict.reasons.append(
            f"hydrophobic gate: the lining would dewet (score "
            f"{wetting.score:.2f} > cutoff, measured over the whole axis)")
    if verdict.sterically_occluded:
        verdict.reasons.append(
            f"sterically occluded: the route narrows to "
            f"{steric:.2f} A, below the radius of a water molecule")
    verdict.conductive = not verdict.reasons
    return verdict
