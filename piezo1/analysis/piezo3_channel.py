"""Is piezo3 built like a working channel? Asked of coordinates, not of sequence.

Split from :mod:`piezo1.analysis.piezo3` at the length limit and along a real
seam: that module answers questions about the *protomer* — its fold, the
residues it kept, which template to build on — and every one of those needs no
assembly. This module needs a trimer, and a trimer has to be built, which is a
different kind of claim carrying a different caveat.

**Read this before any number below.** The dome, the pore and the conduction
verdict are properties of three protomers arranged around an axis. piezo3 has
one protomer and no arrangement, so the arrangement is taken from a deposited
PIEZO1 or PIEZO2 trimer. :func:`piezo1.structure.assembly.borrowed_fraction`
measures how much of the result that accounts for, and on every template tried
the answer is **85–96%**. A dome radius measured this way is mostly a
measurement of the template.

So the numbers here are not evidence that piezo3 does or does not conduct. What
they *can* do is fail: a protomer that could not be arranged into a channel at
all — one whose pore-lining helices point the wrong way, or whose blades cannot
be brought round an axis without interpenetrating — would show up as an
assembly that does not close, and this one does close. The honest statement is
a negative that survived, not a positive that was demonstrated, and
:attr:`Piezo3Channel.verdict` says exactly that.

The comparison structure measured by the identical route is carried alongside,
because a number with nothing beside it invites the reader to supply their own
scale.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..config import STRUCTURE_DIR
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.assembly import assemble_trimer, borrowed_fraction
from ..structure.geometry import measure_dome, tm_surface_points
from ..structure.pore import pore_profile
from ..structure.protomers import protomer_blocks
from ..structure.superpose import detect_c3_axis

__all__ = ["ChannelMeasurement", "Piezo3Channel", "measure_channel",
           "build_channel", "COMPARISON_ENTRY"]

#: The PIEZO1 entry measured by the identical route, for scale. 7WLT is the
#: curved state this project treats as its reference throughout.
COMPARISON_ENTRY = "7WLT"


@dataclass(frozen=True)
class ChannelMeasurement:
    """What the standard pipeline says about one trimer."""

    label: str
    reference: str
    n_helices: int
    radius_of_curvature_nm: float | None
    dome_depth_nm: float | None
    excess_area_nm2: float | None
    c3_angle_deg: float | None
    c3_rmsd_A: float | None
    bottleneck_radius_A: float | None
    bottleneck_z_A: float | None
    wetting_score: float | None = None
    conducts: bool | None = None
    conduction_reason: str = ""

    def summary(self) -> str:
        def g(v, fmt="{:.2f}"):
            return "-" if v is None else fmt.format(v)
        return (f"{self.label}: dome R {g(self.radius_of_curvature_nm)} nm, "
                f"depth {g(self.dome_depth_nm)} nm, excess area "
                f"{g(self.excess_area_nm2, '{:.0f}')} nm2 over "
                f"{self.n_helices} helices; pore bottleneck "
                f"{g(self.bottleneck_radius_A)} A")


@dataclass(frozen=True)
class Piezo3Channel:
    """The assembled piezo3 channel, and everything that makes it a model."""

    template: str
    template_protein: str
    identity: float
    clashes: int
    borrowed: float | None
    piezo3: ChannelMeasurement
    comparison: ChannelMeasurement | None
    assembly_note: str = ""
    caveats: tuple = ()

    @property
    def closes(self) -> bool:
        """Did the three protomers arrange into a channel with an axis at all?"""
        return (self.piezo3.c3_rmsd_A is not None
                and self.piezo3.bottleneck_radius_A is not None)

    @property
    def verdict(self) -> str:
        if not self.closes:
            return ("the protomer could not be arranged into a channel on this "
                    "template; nothing about piezo3 follows, because the "
                    "template may simply be wrong for it")
        borrowed = ("" if self.borrowed is None
                    else f" {self.borrowed:.0%} of the departure from planarity "
                         f"is the template's, not piezo3's, and")
        return (f"the piezo3 protomer arranges into a closed trimer on "
                f"{self.template} with an axis and a continuous lumen, which a "
                f"protein not built like a channel need not have done.{borrowed}"
                f" the dome radius of {self.piezo3.radius_of_curvature_nm:.1f} nm "
                f"is therefore not evidence about piezo3's own curvature. This "
                f"is a negative that survived, not a positive demonstrated: no "
                f"current has ever been recorded from this protein")


def measure_channel(structure: Structure, reference: str, label: str,
                    with_conduction: bool = True) -> ChannelMeasurement:
    """Run the dome, pore and conduction pipeline on any PIEZO trimer.

    One function for both sides of the comparison. Two copies of "measure the
    dome" is exactly how a comparison ends up measuring how each side was
    defined — the reason :func:`tm_surface_points` exists at all.
    """
    blocks, _ = protomer_blocks(structure)
    if len(blocks) < 3:
        return ChannelMeasurement(label=label, reference=reference, n_helices=0,
                                  radius_of_curvature_nm=None, dome_depth_nm=None,
                                  excess_area_nm2=None, c3_angle_deg=None,
                                  c3_rmsd_A=None, bottleneck_radius_A=None,
                                  bottleneck_z_A=None,
                                  conduction_reason="fewer than three protomers")
    points, resolved = tm_surface_points(structure, reference)
    dome = measure_dome(blocks, points)
    profile = pore_profile(structure, detect_c3_axis(blocks))

    score = conducts = None
    reason = ""
    if with_conduction:
        try:
            from .conduction import conduction_verdict
            verdict = conduction_verdict(structure, profile)
            score = getattr(verdict.wetting, "score", None)
            conducts = bool(verdict.conductive)
            reason = verdict.summary()
        except Exception as exc:                              # noqa: BLE001
            reason = f"conduction verdict unavailable: {type(exc).__name__}: {exc}"

    return ChannelMeasurement(
        label=label, reference=reference, n_helices=len(resolved),
        radius_of_curvature_nm=dome.radius_of_curvature / 10,
        dome_depth_nm=dome.dome_depth / 10,
        excess_area_nm2=dome.excess_area / 100,
        c3_angle_deg=dome.notes.get("c3_angle_deg"),
        c3_rmsd_A=dome.notes.get("c3_rmsd"),
        bottleneck_radius_A=profile.bottleneck_radius,
        bottleneck_z_A=profile.bottleneck_z,
        wetting_score=score, conducts=conducts, conduction_reason=str(reason))


def build_channel(template: str | None = None,
                  with_comparison: bool = True) -> Piezo3Channel:
    """Assemble piezo3 into a trimer and run the pipeline on it."""
    from ..io.registry import load_registry
    from .piezo3 import best_paralogue_template, load_model

    template = template or best_paralogue_template()
    if template is None:
        raise RuntimeError("no trimeric template assembles the piezo3 monomer")
    assembly = assemble_trimer(load_model(), template=template)
    try:
        borrowed = borrowed_fraction(assembly)["borrowed_fraction"]
    except Exception:                                          # noqa: BLE001
        borrowed = None

    measured = measure_channel(assembly.structure, "zebrafish_piezo3",
                               f"piezo3 on {template}")
    comparison = None
    if with_comparison:
        path = STRUCTURE_DIR / f"{COMPARISON_ENTRY}.cif"
        if path.exists():
            comparison = measure_channel(Structure.from_file(path), "mouse",
                                         f"{COMPARISON_ENTRY} (mouse PIEZO1)")

    record = {e.pdb: e for e in load_registry().available()}.get(template)
    return Piezo3Channel(
        template=template, template_protein=(record.protein if record else "?"),
        identity=float(assembly.alignment_identity),
        clashes=int(assembly.clashes), borrowed=borrowed,
        piezo3=measured, comparison=comparison,
        assembly_note=assembly.summary(),
        caveats=(
            "the arrangement of the three protomers is the template's, not "
            "piezo3's; the borrowed fraction below says how much of the answer "
            "that is",
            "the protomer is an AlphaFold prediction — its B-factor column is "
            "pLDDT, and this project has measured that an AlphaFold PIEZO "
            "monomer's blades sit 33-45 A from where cryo-EM puts them even "
            "for the protein the model is of",
            f"{assembly.clashes:,} inter-protomer clashes, against 3-8 for a "
            f"deposited trimer: nothing models the contacts between protomers, "
            f"so interpenetration is what an imperfect arrangement looks like",
            "the transmembrane annotation is from an unreviewed UniProt entry "
            "naming 21 helices where PIEZO1 and PIEZO2 have 38, so the dome "
            "surface is traced from a sparser set of points",
            "no current has ever been recorded from any piezo3",
        ))
