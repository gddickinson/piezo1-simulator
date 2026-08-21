"""The census's per-residue constraint, read on this project's own structures.

The ``piezo_genes`` census measured, for every residue of human PIEZO1, human
PIEZO2 and zebrafish piezo3, how much 121–192 orthologues of that gene have been
willing to change it. That is a scalar per residue — the one form of external
result a structural project can genuinely *use* rather than merely quote.

**What this module is for, and what it deliberately is not.** It is not a
re-derivation of the census: the alignments and the statistic are theirs and are
imported whole. It is the join — putting their number on our coordinates, our
domain boundaries and our numbering — and every one of those three is a place
the join can go wrong quietly:

- **Numbering.** The track is in Q92508. A mouse entry read at human numbers is
  off by a non-constant offset, so the conversion goes through
  :mod:`piezo1.core.sequence` and an entry whose numbering we cannot read is
  refused rather than assumed human.
- **Boundaries.** The census partitioned the protein its own way. Repeating its
  per-domain result on *our* ``domains.json`` is a test, not a restatement — the
  two partitions disagree about where the blades stop and the anchor begins, so
  a finding that survives both is a finding about the protein rather than about
  a boundary choice.
- **Unscored residues.** A residue with no score is ``nan``, never zero. Painting
  the two the same colour would show "we did not measure this" as "nothing here
  matters", which is exactly backwards for the blade tips, where coverage is
  worst and the census's own claim is that constraint is *low*.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.annotations import Annotations, load_annotations
from ..core.family import ConstraintTrack, load_constraint, load_family_findings
from ..core.numbering_check import piezo1_numbering
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["DomainConstraint", "ConstraintOnStructure", "StructureRefusal",
           "ConservationCrossCheck", "domain_constraint", "blade_gradient",
           "constraint_on_structure", "selection_track",
           "census_domain_constraint", "compare_with_own_conservation",
           "DomainIdentity", "paralogue_asymmetry", "POREWARD_ORDER"]

#: The census cut its two blade bands at mouse residue 1300, and ended the
#: proximal band where the anchor begins. Both are needed to reproduce its
#: numbers on our sequence, and neither is ours.
CENSUS_CUT_MOUSE = 1300
CENSUS_CUT_HUMAN_FALLBACK = 1305
CENSUS_PROXIMAL_END = 1935

#: The census's own ordering of its domains, pore-ward. Carried so the
#: replication on our boundaries can be reported *against* an expectation
#: instead of being read off after the fact.
POREWARD_ORDER = ("proximal_blades", "distal_blades", "beam", "CED", "anchor",
                  "outer_helix", "CTD", "inner_helix", "pore_linker")


@dataclass(frozen=True)
class DomainConstraint:
    """Constraint summarised over one domain of one partition."""

    domain: str
    category: str
    n_residues: int
    n_scored: int
    mean: float | None
    median: float | None
    vs_whole: float | None

    @property
    def coverage(self) -> float:
        return self.n_scored / self.n_residues if self.n_residues else 0.0


@dataclass(frozen=True)
class StructureRefusal:
    """Why a structure could not be scored. Never a zero-filled array."""

    reason: str
    numbering: str | None = None

    def __bool__(self) -> bool:
        return False


@dataclass(frozen=True)
class ConstraintOnStructure:
    """The track mapped onto one structure's atoms."""

    gene: str
    track: str
    numbering: str
    per_atom: np.ndarray                 # (n_atoms,) nan where unscored
    per_residue: dict = field(default_factory=dict)
    n_atoms_scored: int = 0
    n_residues_scored: int = 0
    n_residues: int = 0
    converted: bool = False
    note: str = ""

    def __bool__(self) -> bool:
        return self.n_atoms_scored > 0

    @property
    def coverage(self) -> float:
        return self.n_residues_scored / self.n_residues if self.n_residues else 0.0


def _finite(values: np.ndarray) -> np.ndarray:
    return values[~np.isnan(values)]


def domain_constraint(track: ConstraintTrack | None = None,
                      annotations: Annotations | None = None,
                      include_sub_elements: bool = False) -> list[DomainConstraint]:
    """The census constraint summarised over **our** domain partition.

    This is the replication of the census's headline ordering — blades least
    constrained, pore machinery most — on boundaries it did not choose. The
    census used ten bands derived from the mouse cryo-EM literature; this uses
    ``domains.json``, which is built from UniProt topology plus Guo & MacKinnon's
    named elements. Where the two agree, the ordering is a property of the
    protein.

    Sub-elements (the cuff's elbow, base, hairpin and PE helix) are excluded by
    default for the same reason :meth:`Annotations.domain_at` skips them: they
    lie *inside* the anchor and the CTD rather than partitioning the chain
    beside them, so including them double-counts those residues.
    """
    track = track or load_constraint("PIEZO1")
    ann = annotations or load_annotations("human")
    whole = float(np.nanmean(track.values))
    out: list[DomainConstraint] = []
    for dom in ann.domains:
        if dom.sub_element and not include_sub_elements:
            continue
        if dom.start is None or dom.end is None:
            continue
        lo, hi = max(1, dom.start), min(track.length, dom.end)
        if hi < lo:
            continue
        values = track.values[lo - 1:hi]
        finite = _finite(values)
        out.append(DomainConstraint(
            domain=dom.id, category=dom.category, n_residues=hi - lo + 1,
            n_scored=int(finite.size),
            mean=float(finite.mean()) if finite.size else None,
            median=float(np.median(finite)) if finite.size else None,
            vs_whole=float(finite.mean() - whole) if finite.size else None))
    return out


def census_domain_constraint(gene: str = "PIEZO1") -> list[dict]:
    """The census's own per-domain numbers, on its own boundaries.

    Read back from the import so a report can put the two partitions side by
    side. Nothing here is recomputed; that is the point of showing it beside
    :func:`domain_constraint`.
    """
    rows = load_family_findings().table("constraint_by_domain")
    return [r for r in rows if r["paralog"] == gene]


def blade_gradient(track: ConstraintTrack | None = None,
                   annotations: Annotations | None = None) -> dict:
    """Is the *distal* blade more constrained than the proximal one?

    The census says yes in all three paralogues, and reports it as the opposite
    of what "peripheral means dispensable" predicts. **It does not survive being
    asked of the blade units**, and the reason is composition rather than
    biology.

    The census's two bands are a single cut of the chain at mouse 1300: whatever
    lies below is "distal blade", whatever lies above up to the anchor is
    "proximal". Both bands therefore contain the long unstructured stretches
    between transmembrane units as well as the units themselves — and they
    contain very different *amounts* of it. Measured here on human PIEZO1:

    * the census bands reproduce (0.649 and 0.558 against their 0.656 and
      0.558), so the import and the numbering join are sound;
    * inter-unit linker scores the same in both bands (0.517 and 0.515), so the
      linker is not more conserved at one end than the other;
    * the proximal band is 77% linker and the distal band 29%;
    * restricted to the THU units themselves the ordering **reverses** —
      THU1–6 at 0.698 against THU7–9 at 0.737.

    So the finding is a statement about how much linker each band happens to
    contain. What is left after the composition is removed is the *opposite*
    gradient, and a weaker one. Reported as a disagreement rather than resolved,
    because the census's bands are the ones its other numbers were computed on.
    """
    track = track or load_constraint("PIEZO1")
    ann = annotations or load_annotations("human")
    thus = {d.id: d for d in ann.domains
            if d.id.startswith("thu") and d.start and d.end}
    # Resolved at call time, so an override in the dialog takes effect on the
    # next call rather than at import. Which THU the split falls at is exactly
    # what the census's distal-versus-proximal finding turns on.
    last_distal = int(_P.value("family.distal_last_thu"))

    in_unit = np.zeros(track.length, dtype=bool)
    distal = np.zeros(track.length, dtype=bool)
    proximal = np.zeros(track.length, dtype=bool)
    per_thu = {}
    for name, dom in thus.items():
        span = slice(dom.start - 1, dom.end)
        in_unit[span] = True
        (distal if int(name[3:]) <= last_distal else proximal)[span] = True
        finite = _finite(track.values[span])
        if finite.size:
            per_thu[name] = float(finite.mean())

    resi = np.arange(1, track.length + 1)
    cut = _census_cut(track)
    band_distal = resi <= cut
    band_proximal = (resi > cut) & (resi <= CENSUS_PROXIMAL_END)

    def summarise(mask):
        finite = _finite(track.values[mask])
        return (float(finite.mean()) if finite.size else None, int(finite.size))

    unit_distal, n_ud = summarise(distal)
    unit_proximal, n_up = summarise(proximal)
    band_d, n_bd = summarise(band_distal)
    band_p, n_bp = summarise(band_proximal)
    link_d, n_ld = summarise(band_distal & ~in_unit)
    link_p, n_lp = summarise(band_proximal & ~in_unit)

    holds_on_bands = (band_d is not None and band_p is not None and band_d > band_p)
    holds_on_units = (unit_distal is not None and unit_proximal is not None
                      and unit_distal > unit_proximal)
    return {
        "gene": track.gene,
        "track": track.track,
        "per_thu": dict(sorted(per_thu.items())),
        "band_cut_human": cut,
        "band_distal": band_d, "n_band_distal": n_bd,
        "band_proximal": band_p, "n_band_proximal": n_bp,
        "unit_distal": unit_distal, "n_unit_distal": n_ud,
        "unit_proximal": unit_proximal, "n_unit_proximal": n_up,
        "linker_distal": link_d, "n_linker_distal": n_ld,
        "linker_proximal": link_p, "n_linker_proximal": n_lp,
        # Fraction of each *band* that is not part of any transmembrane unit.
        # Not 1 - n_unit/n_band: the distal band also contains THU7, which the
        # distal unit mask does not, and using the unit mask here understated
        # the gap that is the whole mechanism.
        "linker_fraction_distal": n_ld / n_bd if n_bd else None,
        "linker_fraction_proximal": n_lp / n_bp if n_bp else None,
        "holds_on_census_bands": holds_on_bands,
        "holds_on_thu_units": holds_on_units,
        "census_distal": 0.6561,
        "census_proximal": 0.5583,
        "verdict": _blade_verdict(holds_on_bands, holds_on_units,
                                  link_d, link_p),
    }


def _census_cut(track: ConstraintTrack) -> int:
    """The census's band boundary, in this track's numbering.

    Their cut is at mouse 1300 and the track is human, so it goes through the
    alignment map. If the map cannot place it the fallback is stated rather than
    silent, because a cut in the wrong place would move every number above.
    """
    if track.numbering != "human":
        return CENSUS_CUT_MOUSE
    from ..core.sequence import mouse_to_human
    mapped = mouse_to_human(CENSUS_CUT_MOUSE)
    return int(mapped) if mapped else CENSUS_CUT_HUMAN_FALLBACK


def _blade_verdict(holds_bands: bool, holds_units: bool,
                   link_d: float | None, link_p: float | None) -> str:
    if holds_bands and not holds_units:
        gap = (abs(link_d - link_p) if link_d is not None and link_p is not None
               else None)
        detail = (f"and inter-unit linker scores the same either side "
                  f"({gap:.3f} apart)" if gap is not None else "")
        return ("boundary-dependent: it holds on the census's chain-cut bands "
                f"and reverses on the transmembrane units {detail}".strip())
    if holds_bands and holds_units:
        return "holds on both partitions"
    if not holds_bands:
        return ("does not reproduce even on the census's own bands - suspect "
                "the import or the numbering join before the finding")
    return "holds on the units but not on the bands"


def constraint_on_structure(structure: Structure,
                            track: ConstraintTrack | None = None
                            ) -> ConstraintOnStructure | StructureRefusal:
    """Map the track onto every atom of a structure, or refuse and say why.

    The refusal is the important half. The track is human PIEZO1's; a mouse
    entry needs the alignment-backed conversion, and anything that is not
    PIEZO1 at all — a PIEZO2 entry, PEZO-1, 6LQI's splice numbering — has no
    business being coloured by it. Returning zeros for those would produce a
    picture indistinguishable from a genuinely unconstrained protein.
    """
    track = track or load_constraint("PIEZO1")
    numbering = piezo1_numbering(structure)
    if numbering is None:
        return StructureRefusal(
            "this entry's numbering is not one annotation can be read in "
            "(a paralogue, a non-vertebrate PIEZO, or a splice numbering); "
            "the PIEZO1 constraint track cannot be placed on it")
    if numbering not in ("human", "mouse"):
        return StructureRefusal(
            f"no constraint track for {numbering} numbering", numbering)

    converted = numbering == "mouse"
    per_atom = np.full(structure.n_atoms, np.nan)
    per_residue: dict[int, float] = {}
    for resi in np.unique(structure.res_seq):
        # res_seq is in the entry's own numbering, which is what `converted`
        # names. A mouse entry is looked up through the alignment map, never by
        # adding an offset.
        value = (_mouse_lookup(track, int(resi)) if converted
                 else track.value(int(resi)))
        if value is None:
            continue
        mask = structure.res_seq == resi
        per_atom[mask] = value
        per_residue[int(resi)] = value
    return ConstraintOnStructure(
        gene=track.gene, track=track.track, numbering=numbering,
        per_atom=per_atom, per_residue=per_residue,
        n_atoms_scored=int(np.count_nonzero(~np.isnan(per_atom))),
        n_residues_scored=len(per_residue),
        n_residues=int(np.unique(structure.res_seq).size),
        converted=converted,
        note=("residue numbers converted mouse->human through the alignment map"
              if converted else "read directly at human PIEZO1 numbering"))


def _mouse_lookup(track: ConstraintTrack, mouse_resi: int) -> float | None:
    """The human-numbered score at a mouse residue number.

    Goes through the alignment map both ways rather than by subtraction: the
    human/mouse offset is not constant, and this project has a standing rule
    against any code that pretends it is.
    """
    from ..core.sequence import mouse_to_human
    human = mouse_to_human(mouse_resi)
    return None if human is None else track.value(human)


def selection_track(gene: str = "PIEZO1") -> list[dict]:
    """The codon-level result, which shares no statistic with the JSD track.

    Carried because agreement between two routes that share nothing is worth
    more than either alone — and because its coverage is visibly worse (70-74%
    of positions carry a usable rate against ~100% for the amino-acid track),
    which is the honest reason it is the confirmation and not the headline.
    """
    rows = load_family_findings().table("selection_by_domain")
    return [r for r in rows if r["paralog"] == gene]


@dataclass(frozen=True)
class ConservationCrossCheck:
    """This project's own conservation profile against the census's constraint.

    Two evolutionary measurements of the same protein that share no data and no
    statistic. Ours is 1 - normalised Shannon entropy over a few dozen fetched
    vertebrate orthologues; the census's is Jensen-Shannon divergence from a
    background distribution over 174 genome-backed PIEZO1 loci, one per genome.

    They should agree, and how much they *do* is the useful number. A high
    correlation would say the cheap local route recovers the expensive one; a
    low one says our column is measuring something narrower, which matters
    because ``features.py`` ships it as ``conservation`` and several rankings in
    this project use it.
    """

    n: int
    pearson: float
    spearman: float
    n_orthologues_census: int
    n_orthologues_ours: int | None
    note: str = ""

    @property
    def agree(self) -> bool:
        """Rank correlation above 0.5 — half the ordering shared."""
        return self.spearman >= 0.5

    def summary(self) -> str:
        ours = "unknown" if self.n_orthologues_ours is None else str(self.n_orthologues_ours)
        return (f"over {self.n} shared positions the two routes correlate at "
                f"rho = {self.spearman:.2f} (Pearson {self.pearson:.2f}); the "
                f"census used {self.n_orthologues_census} orthologues to our "
                f"{ours}")


def compare_with_own_conservation(track: ConstraintTrack | None = None,
                                  profile=None) -> ConservationCrossCheck | None:
    """Correlate the imported constraint with this project's own conservation.

    Returns ``None`` rather than a number when our own profile is unavailable —
    it needs fetched orthologues, and on a fresh clone there are none. A zero
    correlation and an absent input must not be reachable through the same
    value.
    """
    from .conservation import conservation_profile, load_orthologs
    from .fluctuations import pearson, spearman

    track = track or load_constraint("PIEZO1")
    if profile is None:
        try:
            profile = conservation_profile(load_orthologs())
        except Exception:                                     # noqa: BLE001
            return None
    if profile is None:
        return None

    # Joined by residue number, not by position: the profile covers only the
    # residues its alignment reached, and zipping two arrays of different
    # lengths would slide one against the other silently.
    ours, theirs = [], []
    for resi, score, coverage in zip(profile.residues, profile.conservation,
                                     profile.coverage):
        value = track.value(int(resi))
        if value is None or np.isnan(score) or coverage <= 0:
            continue
        ours.append(float(score))
        theirs.append(value)
    if len(ours) < 50:
        return None
    ours, theirs = np.array(ours), np.array(theirs)
    return ConservationCrossCheck(
        n=int(ours.size),
        pearson=float(pearson(ours, theirs)),
        spearman=float(spearman(ours, theirs)),
        n_orthologues_census=track.n_orthologues,
        n_orthologues_ours=getattr(profile, "n_orthologs", None),
        note=("ours is 1 - normalised Shannon entropy over fetched vertebrate "
              "orthologues; the census's is Jensen-Shannon divergence over "
              "genome-backed loci. No data and no statistic in common."))

# Split to ``paralogue_identity`` at the 500-line limit. Re-exported so every
# existing import of ``family_constraint.paralogue_asymmetry`` still works —
# the seam is between measuring *within* a paralogue and *between* two, and a
# caller that only wants a number should not have to know where it moved.
from .paralogue_identity import DomainIdentity, paralogue_asymmetry  # noqa: E402,F401
