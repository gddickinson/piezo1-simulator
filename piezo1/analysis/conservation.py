"""Evolutionary conservation across PIEZO1 orthologs.

Conservation answers a question structure cannot: which residues has selection
refused to let change? Crossed with the variant table it answers a more useful
one — **which strongly constrained positions have no reported variant at all**,
and are therefore candidate functional sites that nobody has yet tested.

The alignment is *reference-anchored*: every ortholog is aligned pairwise to
human PIEZO1 and its residues indexed by the human position they align to. This
is not a true multiple alignment, and for distant homologues that distinction
matters. Here it does not much: PIEZO1 orthologs across vertebrates are
colinear, with no domain rearrangement, and anchoring to human is what keeps
every number in the human numbering the rest of the project uses.

**A caution on interpretation.** A conserved residue is not necessarily a
functional one — it may simply be structurally load-bearing, and buried
positions are conserved for reasons that have nothing to do with mechanism.
Conservation narrows a search; it does not identify a mechanism.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

import numpy as np

from ..config import SEQUENCE_DIR
from ..core.sequence import align_global, human_sequence

__all__ = ["Ortholog", "OrthologSet", "ConservationProfile", "fetch_orthologs",
           "load_orthologs", "conservation_profile", "constrained_positions",
           "rank_candidates", "AMINO_ACIDS"]

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
UA = "piezo1-simulator/0.1 (research use)"


@dataclass
class Ortholog:
    accession: str
    organism: str
    length: int
    sequence: str
    reviewed: bool = False


@dataclass
class OrthologSet:
    members: list[Ortholog] = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.members)

    def organisms(self) -> list[str]:
        return [m.organism for m in self.members]


# --------------------------------------------------------------------------
# Retrieval
# --------------------------------------------------------------------------

def fetch_orthologs(min_length: int = 2000, max_length: int = 3000,
                    taxon: int = 7742, size: int = 100,
                    force: bool = False) -> OrthologSet:
    """Download vertebrate PIEZO1 orthologs from UniProt, one per species.

    ``taxon`` 7742 is Vertebrata. The length window excludes fragments and
    mis-annotated entries, which are common among the unreviewed records and
    would otherwise contribute long gap runs that look like divergence.

    Only one entry per organism is kept — reviewed first, then longest.
    Multiple entries for one species are isoforms or partial records, and
    counting a species twice silently weights it twice in the conservation.
    """
    cache = SEQUENCE_DIR / "piezo1_orthologs.json"
    if cache.exists() and not force:
        return load_orthologs()

    query = (f"gene:PIEZO1 AND taxonomy_id:{taxon} AND "
             f"length:[{min_length} TO {max_length}]")
    url = (f"{UNIPROT}?query={urllib.parse.quote(query)}&format=json&size={size}"
           f"&fields=accession,id,organism_name,length,reviewed,sequence")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=180) as r:
        data = json.loads(r.read())

    best: dict[str, Ortholog] = {}
    for rec in data.get("results", []):
        organism = rec.get("organism", {}).get("scientificName", "?")
        reviewed = rec.get("entryType", "").startswith("UniProtKB reviewed")
        o = Ortholog(accession=rec["primaryAccession"], organism=organism,
                     length=rec["sequence"]["length"],
                     sequence=rec["sequence"]["value"], reviewed=reviewed)
        current = best.get(organism)
        if current is None or (o.reviewed, o.length) > (current.reviewed, current.length):
            best[organism] = o

    members = sorted(best.values(), key=lambda m: (not m.reviewed, m.organism))
    SEQUENCE_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(
        {"query": query, "n_species": len(members),
         "orthologs": [vars(m) for m in members]}, indent=1))
    return OrthologSet(members=members,
                       meta={"query": query, "source": "UniProt"})


def load_orthologs() -> OrthologSet:
    cache = SEQUENCE_DIR / "piezo1_orthologs.json"
    if not cache.exists():
        raise FileNotFoundError(
            f"{cache} missing — run fetch_orthologs() once with network access")
    data = json.loads(cache.read_text())
    return OrthologSet(members=[Ortholog(**o) for o in data["orthologs"]],
                       meta={"query": data.get("query"), "source": "cache"})


# --------------------------------------------------------------------------
# Conservation
# --------------------------------------------------------------------------

@dataclass
class ConservationProfile:
    """Per-residue conservation in human PIEZO1 numbering."""

    residues: np.ndarray          # 1-based human positions
    entropy: np.ndarray           # Shannon entropy, normalised to [0, 1]
    identity: np.ndarray          # fraction of orthologs matching human
    coverage: np.ndarray          # fraction of orthologs aligning at all
    n_orthologs: int = 0
    organisms: tuple[str, ...] = ()
    meta: dict = field(default_factory=dict)

    @property
    def conservation(self) -> np.ndarray:
        """1 − normalised entropy. 1 is invariant, 0 is maximally variable."""
        return 1.0 - self.entropy

    def at(self, residue: int) -> dict | None:
        idx = np.flatnonzero(self.residues == residue)
        if len(idx) == 0:
            return None
        i = int(idx[0])
        return {"residue": residue, "conservation": float(self.conservation[i]),
                "identity": float(self.identity[i]),
                "coverage": float(self.coverage[i])}

    def top_conserved(self, n: int = 20, min_coverage: float = 0.7
                      ) -> list[tuple[int, float]]:
        ok = self.coverage >= min_coverage
        order = np.argsort(np.where(ok, self.conservation, -np.inf))[::-1][:n]
        return [(int(self.residues[i]), float(self.conservation[i])) for i in order]

    def domain_means(self, annotations) -> dict[str, float]:
        out: dict[str, list[float]] = {}
        for res, cons, cov in zip(self.residues, self.conservation, self.coverage):
            if cov < 0.5:
                continue
            d = annotations.domain_at(int(res))
            out.setdefault(d.id if d else "none", []).append(float(cons))
        return {k: float(np.mean(v)) for k, v in out.items() if v}


def conservation_profile(orthologs: OrthologSet,
                         reference: str | None = None,
                         min_coverage: float = 0.5) -> ConservationProfile:
    """Per-position conservation, anchored to human numbering.

    Each ortholog is aligned pairwise to the reference and its residues mapped
    onto reference positions. Gaps are excluded from the frequency count but
    tracked as ``coverage``, so a column that only three species align to is
    not mistaken for an invariant one.
    """
    ref = reference or human_sequence()
    n = len(ref)
    counts = np.zeros((n, len(AMINO_ACIDS)), dtype=np.float64)
    aligned = np.zeros(n, dtype=np.float64)
    identical = np.zeros(n, dtype=np.float64)
    index = {a: i for i, a in enumerate(AMINO_ACIDS)}

    # Skip the reference entry itself, but only once. Skipping every sequence
    # equal to the reference would silently discard genuinely identical
    # orthologs from closely related species, and with a small set that empties
    # the alignment entirely.
    reference_seen = False
    for member in orthologs.members:
        if not reference_seen and member.sequence == ref:
            reference_seen = True
            continue
        gap_ref, gap_other = align_global(ref, member.sequence)
        pos = 0
        for a, b in zip(gap_ref, gap_other):
            if a == "-":
                continue
            if b != "-" and b in index:
                counts[pos, index[b]] += 1.0
                aligned[pos] += 1.0
                identical[pos] += float(a == b)
            pos += 1

    total = np.maximum(counts.sum(axis=1), 1.0)
    freq = counts / total[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        logs = np.where(freq > 0, np.log(freq), 0.0)
    entropy = -(freq * logs).sum(axis=1) / np.log(len(AMINO_ACIDS))
    # A column nothing aligns to is unknown, not conserved.
    n_other = max(len(orthologs.members) - (1 if reference_seen else 0), 1)
    coverage = aligned / n_other
    entropy = np.where(coverage > 0, entropy, 1.0)

    return ConservationProfile(
        residues=np.arange(1, n + 1), entropy=entropy,
        identity=np.divide(identical, np.maximum(aligned, 1.0)),
        coverage=coverage, n_orthologs=n_other,
        organisms=tuple(orthologs.organisms()),
        meta={"min_coverage": min_coverage, "reference_length": n,
              "note": ("reference-anchored pairwise alignment, not a true MSA; "
                       "adequate here because vertebrate PIEZO1 orthologs are "
                       "colinear")})


# --------------------------------------------------------------------------
# Crossing conservation with the variant record
# --------------------------------------------------------------------------

def constrained_positions(profile: ConservationProfile, annotations,
                          conservation_threshold: float = 0.9,
                          min_coverage: float = 0.7,
                          resolved: set[int] | None = None) -> list[dict]:
    """Highly conserved positions with **no reported variant**.

    The intersection worth looking at: selection has held these residues fixed,
    yet no clinical or engineered variant has ever been recorded there. Those
    are candidates for mutagenesis — positions the literature has not tested.

    ``resolved`` optionally restricts to residues present in a structure, since
    an untestable prediction is not much use.
    """
    known = {v.residue for v in annotations.variants if v.residue is not None}
    out = []
    for res, cons, cov in zip(profile.residues, profile.conservation,
                              profile.coverage):
        res = int(res)
        if cov < min_coverage or cons < conservation_threshold:
            continue
        if res in known:
            continue
        if resolved is not None and res not in resolved:
            continue
        domain = annotations.domain_at(res)
        groups = [g.label for g in annotations.residue_groups
                  if res in g.residues]
        out.append({"residue": res, "conservation": float(cons),
                    "coverage": float(cov),
                    "domain": domain.id if domain else None,
                    "domain_name": domain.name if domain else None,
                    "annotated_sites": groups})
    out.sort(key=lambda r: (-r["conservation"], r["residue"]))
    return out


def rank_candidates(candidates: list[dict], features: dict[str, dict],
                    weights: dict[str, float] | None = None) -> list[dict]:
    """Rank constrained positions by additional per-residue evidence.

    Conservation alone is not selective enough to be a hypothesis. Across 62
    vertebrate orthologs a quarter of PIEZO1 is invariant, so "invariant and
    never mutated in the literature" returns hundreds of positions — true, and
    nearly useless.

    The intersection worth looking at is **conserved *and* mechanically coupled
    to the gate**. Neither is specific alone: conservation catches everything
    structurally load-bearing, and mechanical coupling catches everything near
    the pore. Positions scoring highly on both are a much smaller set, and they
    are the ones this project can nominate that a sequence-only method cannot.

    ``features`` maps a name to a ``{residue: value}`` dict — for example the
    PRS gate response and path betweenness from
    :mod:`piezo1.analysis.allostery`. Each is converted to a percentile rank so
    that quantities in different units combine sensibly, then averaged with
    ``weights``.

    This ranks hypotheses for testing. It is not a prediction of function.
    """
    if not candidates:
        return []
    weights = weights or {k: 1.0 for k in features}

    ranks: dict[str, dict[int, float]] = {}
    for name, table in features.items():
        if not table:
            continue
        residues = np.array(sorted(table))
        values = np.array([table[r] for r in residues], dtype=float)
        order = np.argsort(np.argsort(values))
        pct = order / max(len(values) - 1, 1)
        ranks[name] = {int(r): float(p) for r, p in zip(residues, pct)}

    out = []
    for c in candidates:
        row = dict(c)
        total = 0.0
        used = 0.0
        for name, table in ranks.items():
            value = table.get(c["residue"])
            row[f"{name}_percentile"] = value
            if value is not None:
                w = weights.get(name, 1.0)
                total += w * value
                used += w
        # Conservation is already near 1 for every candidate by construction,
        # so it is not re-weighted here; the ranking is on the extra evidence.
        row["combined_score"] = total / used if used else float("nan")
        row["n_features"] = int(used > 0)
        out.append(row)

    out.sort(key=lambda r: (-(r["combined_score"] if np.isfinite(r["combined_score"])
                              else -1), r["residue"]))
    return out
