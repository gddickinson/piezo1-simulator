"""Does the machinery survive the family? — curated residues across all ten PIEZOs.

``homology.py`` measures how related the family members are. This asks the
question that being related is only interesting for: **the gate, the
selectivity glutamates, the PIP2 lysines, the Yoda1 pocket — are they still
there in a fly, a plant, an amoeba?** A residue conserved from human to
Dictyostelium has been held by selection across the whole eukaryotic tree, and
a curated site that turns out not to be is a fact about how far this project's
annotation may be transferred.

**The instrument is an alignment, and below the twilight line an alignment is
exactly the thing that stops being trustworthy.** That is not a caveat added
afterwards; it is why the module is shaped the way it is. A position mapped
into Arabidopsis PIEZO through a global alignment lands *somewhere* — the
alignment always produces an answer — and if the surrounding block is no better
than chance, the letter found there is a lottery ticket, not a measurement.
Reporting "the gate is not conserved in plants" off such a mapping would be
manufacturing a finding of exactly the kind this project has been caught by
before.

So every mapped position is gated by :mod:`piezo1.analysis.alignment_windows`,
which asks whether the alignment is in register around it and is where that
instrument's own calibration and its two corrected mistakes are written up. A
position it cannot vouch for is reported as ``unreliable`` — never as
conserved, and never as *not* conserved, which is the reading that would do the
damage.

**What it measures.** Three things, and none of them was the expected answer.

*The cap does not travel.* The three curated cap groups — the constriction, the
gate between neighbouring subunits, and the two loops whose separation opens
the lateral portals — become **unreadable outside the vertebrates**: no window
around them clears background in the worm, the fly, the plant or the amoeba.
That is consistent with the cap being the PIEZO-specific extracellular domain
and it is a limit on this project's own annotation, since every cap-gate number
in ``liu2025_panels`` is quoted from mouse residue numbers.

*The gate erodes gradually rather than being conserved or not.* 3 of 3
identical across the mammalian PIEZO1s, 2 of 3 in PIEZO2, 1 of 3 in the worm,
the fly and the plant — and in the worm the one that survives is V2450, whose
two neighbours become V and L. A hydrophobic gate stays hydrophobic without
staying the same residue, which is what a gate made of a property rather than a
contact would look like.

*One group is universal, and it is not the pore.* The **anchor brake** — human
P2113 and F2114, the apex of the anchor domain — is identical in every member
where it can be read, Dictyostelium included. Nothing else is. The per-group
numbers come from ``report()`` rather than being written down here, because
they move with the registered alignment parameters.
"""

from __future__ import annotations

from dataclasses import dataclass, field


from ..core.annotations import load_annotations
from ..parameters import PARAMETERS as _P
from .alignment_windows import alignment_windows
from .homology import family

__all__ = ["SitePosition", "GroupConservation", "SiteReport", "CONSERVATIVE",
           "map_positions", "conservation_of", "report", "family_lengths"]

#: Substitutions counted as conservative — the BLOSUM62 positive-score
#: groupings, collapsed to the classes a reader of a channel paper means by
#: "the charge is preserved" or "still hydrophobic". Used only to *label* a
#: substitution, never to decide whether a site is conserved: that decision is
#: identity, so a conservative call can never quietly become a match.
CONSERVATIVE = (
    frozenset("AVLIM"),      # aliphatic
    frozenset("FWY"),        # aromatic
    frozenset("ST"),         # small hydroxyl
    frozenset("DE"),         # acidic
    frozenset("KR"),         # basic
    frozenset("NQ"),         # amide
    frozenset("HKR"),        # basic including histidine
)


def _classify(a: str, b: str) -> str:
    if a == b:
        return "identical"
    if any(a in group and b in group for group in CONSERVATIVE):
        return "conservative"
    return "different"


@dataclass(frozen=True)
class SitePosition:
    """One curated residue, looked up in one family member."""

    human_residue: int
    human_aa: str
    target: str
    target_residue: int | None
    target_aa: str | None
    #: Mean BLOSUM62 score per aligned column over the window, and the mean and
    #: spread of the same statistic under a composition-matched shuffle of the
    #: target. Judged per pair rather than against a fixed cut, because the
    #: background depends on both compositions.
    window_score: float
    window_null: float
    window_null_sd: float
    #: Carried for display only. It is what a reader expects to see and it is
    #: the statistic this module exists to argue against relying on.
    window_identity: float

    @property
    def aligned(self) -> bool:
        return self.target_residue is not None

    @property
    def window_z(self) -> float:
        if self.window_null_sd == 0:
            return float("inf")
        return (self.window_score - self.window_null) / self.window_null_sd

    @property
    def reliable(self) -> bool:
        """Is the alignment here distinguishable from chance over the window?

        The whole-sequence statistic cannot answer this: a global alignment of
        two twilight-zone sequences is at background *on average* while still
        being correct in the conserved core, which is the region every curated
        site is in.
        """
        return self.aligned and self.window_z >= _P.value("homology.min_z")

    @property
    def status(self) -> str:
        if not self.aligned:
            return "not aligned"
        if not self.reliable:
            return "unreliable"
        return _classify(self.human_aa, self.target_aa or "")

    def summary(self) -> str:
        found = ("—" if self.target_aa is None
                 else f"{self.target_aa}{self.target_residue}")
        return (f"{self.human_aa}{self.human_residue} -> {found} in "
                f"{self.target} ({self.status}, window score "
                f"{self.window_score:+.2f} vs null {self.window_null:+.2f}"
                f"±{self.window_null_sd:.2f}, z {self.window_z:.1f}; "
                f"identity {self.window_identity:.2f})")



@dataclass
class GroupConservation:
    """One curated residue group, across one family member."""

    group_id: str
    label: str
    category: str
    target: str
    positions: tuple[SitePosition, ...]

    def _count(self, status: str) -> int:
        return sum(1 for p in self.positions if p.status == status)

    @property
    def n_identical(self) -> int:
        return self._count("identical")

    @property
    def n_conservative(self) -> int:
        return self._count("conservative")

    @property
    def n_unreliable(self) -> int:
        return self._count("unreliable") + self._count("not aligned")

    @property
    def n_readable(self) -> int:
        return len(self.positions) - self.n_unreliable

    @property
    def identity(self) -> float | None:
        """Fraction identical **among the positions that can be read**.

        ``None`` rather than 0.0 when nothing can be read. A group whose every
        position falls in an untrustworthy block would otherwise report perfect
        non-conservation, which is the most confident possible way to say
        nothing.
        """
        return (self.n_identical / self.n_readable) if self.n_readable else None

    def summary(self) -> str:
        if self.identity is None:
            return (f"{self.label} in {self.target}: none of "
                    f"{len(self.positions)} positions land in a block the "
                    f"alignment can be trusted in")
        return (f"{self.label} in {self.target}: {self.n_identical}/"
                f"{self.n_readable} identical"
                + (f", {self.n_conservative} conservative"
                   if self.n_conservative else "")
                + (f" ({self.n_unreliable} not readable)"
                   if self.n_unreliable else ""))


@dataclass
class SiteReport:
    """Every curated group against every family member."""

    rows: tuple[GroupConservation, ...]
    background: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def for_group(self, group_id: str) -> list[GroupConservation]:
        return [r for r in self.rows if r.group_id == group_id]

    def for_target(self, target: str) -> list[GroupConservation]:
        return [r for r in self.rows if r.target == target]

    @property
    def group_ids(self) -> list[str]:
        seen: list[str] = []
        for row in self.rows:
            if row.group_id not in seen:
                seen.append(row.group_id)
        return seen

    def universal(self) -> list[str]:
        """Groups identical in **every** member where they can be read at all.

        The strong statement, and the reason for the reliability gate: without
        it this would silently include groups that are "conserved" in the plant
        because a scrambled block happened to offer the right letter.
        """
        out = []
        for group_id in self.group_ids:
            rows = [r for r in self.for_group(group_id) if r.n_readable]
            if rows and all(r.identity == 1.0 for r in rows):
                out.append(group_id)
        return out

    def lost_beyond(self, group: str) -> list[str]:
        """Groups readable and not conserved outside ``group``'s members."""
        keys = {m.key for m in family() if m.group == group}
        out = []
        for group_id in self.group_ids:
            outside = [r for r in self.for_group(group_id)
                       if r.target not in keys and r.n_readable]
            if outside and all((r.identity or 0.0) < 1.0 for r in outside):
                out.append(group_id)
        return out

    def summary(self) -> str:
        universal = self.universal()
        return (f"{len(self.group_ids)} curated groups across "
                f"{len({r.target for r in self.rows})} family members; "
                f"{len(universal)} identical wherever they are readable"
                + (f" ({', '.join(universal)})" if universal else ""))


def conservation_of(group, target: str) -> GroupConservation:
    """One curated group in one family member."""
    sequence = _target_sequence(target)
    windows = alignment_windows(target)
    detail = {d["human"]: d.get("human_aa", "") for d in group.detail}

    positions = []
    for residue in group.residues:
        mapped, score, identity = windows.at(int(residue))
        positions.append(SitePosition(
            human_residue=int(residue),
            human_aa=detail.get(int(residue), _human_aa(int(residue))),
            target=target, target_residue=mapped,
            target_aa=(sequence[mapped - 1] if mapped else None),
            window_score=score, window_null=windows.null_mean,
            window_null_sd=windows.null_sd, window_identity=identity))
    return GroupConservation(group_id=group.id, label=group.label,
                             category=group.category, target=target,
                             positions=tuple(positions))


def _target_sequence(key: str) -> str:
    from ..core.numbering_check import reference_entry

    return reference_entry(key)["sequence"]


def _human_aa(residue: int) -> str:
    from ..core.sequence import human_sequence

    sequence = human_sequence()
    return sequence[residue - 1] if 1 <= residue <= len(sequence) else "?"


def report(targets=None, groups=None) -> SiteReport:
    """The whole table: every curated group in every family member.

    ``human`` is included deliberately, as the positive control. Every position
    must come back identical with window identity 1.0, and if one does not, the
    mapping is broken and nothing else in the table means anything.
    """
    annotations = load_annotations("human")
    wanted = list(groups) if groups is not None else \
        [g.id for g in annotations.residue_groups]
    selected = [g for g in annotations.residue_groups if g.id in wanted]
    members = list(targets) if targets is not None else \
        [m.key for m in family()]

    backgrounds = {}
    rows = []
    for target in members:
        windows = alignment_windows(target)
        backgrounds[target] = (windows.null_mean, windows.null_sd)
        for group in selected:
            rows.append(conservation_of(group, target))
    return SiteReport(
        rows=tuple(rows), background=backgrounds,
        meta={"source": "human", "n_groups": len(selected),
              "window": int(_P.value("homology.site_window")),
              "note": "a position whose local window is at background is "
                      "reported as unreliable, not as non-conserved",
              "lengths": {m.key: m.length for m in family()
                          if m.key in members},
              "helices": {m.key: m.n_transmembrane for m in family()
                          if m.key in members}})


def family_lengths() -> dict:
    """Length and helix count per member — the reason indices cannot transfer."""
    return {m.key: (m.length, m.n_transmembrane) for m in family()}
