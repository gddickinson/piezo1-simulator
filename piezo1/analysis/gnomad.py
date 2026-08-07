"""Population constraint from gnomAD, cached and offline-tolerant.

Rounds 22, 34 and 36 all ended in the same place: not enough phenotyped
variants, and no loss-of-function structure. gnomAD is the one source that gives
a phenotype-adjacent signal at **every residue** rather than at 39 — how much
natural selection has removed variation there.

**What it can and cannot say.** Constraint is not a direction. A depleted
position is one where variation is selected against, which is consistent with
either gain- or loss-of-function being harmful. Any test built on it is asking
whether the two classes sit in differently-constrained *positions*, which is a
position-level question and therefore vulnerable to exactly the confound that
killed Round 7 — where 99.8% of the predictor's variance was between-position.
That is why Round 41 pre-registers before testing.

**The gene-level answer arrived first and it is discouraging**, so it is
recorded here rather than buried: PIEZO1 is **not a constrained gene**.
LOEUF 1.10 and pLI ≈ 0 mean loss-of-function variants are not depleted at all,
which is consistent with the common E756del allele; and ``oe_mis`` 1.45 with
``mis_z`` −11.3 means missense variation is *enriched* rather than depleted.
A signal that does not exist gene-wide may still exist regionally, but the prior
is worse than the round assumed.

Follows :mod:`piezo1.analysis.external`: cache to disk, degrade to the cache
when offline, and record the licence. gnomAD data are released under a
permissive licence for research use; the cache is git-ignored like every other
download.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from dataclasses import dataclass, field

import numpy as np

__all__ = ["GnomadClient", "GeneConstraint", "MissenseDensity",
           "missense_density", "GNOMAD_API", "GNOMAD_LICENCE",
           "GNOMAD_CITATION"]

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_LICENCE = ("gnomAD data are made freely available for research use; see "
                  "https://gnomad.broadinstitute.org/policies")
GNOMAD_CITATION = "chen2024gnomad"
_USER_AGENT = "piezo1-simulator/0.1 (research use)"


@dataclass
class GeneConstraint:
    """Gene-level constraint, as gnomAD reports it."""

    gene_id: str
    symbol: str
    pli: float | None = None
    oe_lof: float | None = None
    loeuf: float | None = None          # oe_lof_upper, the recommended metric
    mis_z: float | None = None
    obs_mis: int | None = None
    exp_mis: float | None = None
    oe_mis: float | None = None
    meta: dict = field(default_factory=dict)

    @property
    def lof_intolerant(self) -> bool:
        """LOEUF < 0.35 is the conventional threshold for constraint."""
        return self.loeuf is not None and self.loeuf < 0.35

    @property
    def missense_depleted(self) -> bool:
        """Positive mis_z means fewer missense than expected."""
        return self.mis_z is not None and self.mis_z > 0

    def summary(self) -> str:
        return (f"{self.symbol}: LOEUF {self.loeuf:.2f} "
                f"({'constrained' if self.lof_intolerant else 'LoF-TOLERANT'}), "
                f"pLI {self.pli:.3g}, mis_z {self.mis_z:+.1f}, "
                f"o/e missense {self.oe_mis:.2f} "
                f"({'depleted' if self.missense_depleted else 'ENRICHED'})")


@dataclass
class MissenseDensity:
    """Observed missense variation per residue, and a local rate around it."""

    residues: np.ndarray               # 1..n_residues
    observed: np.ndarray               # count of distinct missense alleles
    allele_count: np.ndarray           # summed allele count
    window: int = 25
    meta: dict = field(default_factory=dict)

    def local_rate(self) -> np.ndarray:
        """Missense alleles per residue in a window centred on each position.

        A single residue carries at most a handful of observed variants, so a
        per-residue count is mostly shot noise. The window is what makes the
        quantity a *regional* constraint estimate — which is the thing the
        literature says is informative — rather than a per-site coin flip.
        """
        kernel = np.ones(2 * self.window + 1)
        padded = np.convolve(self.observed, kernel, mode="same")
        counts = np.convolve(np.ones_like(self.observed), kernel, mode="same")
        return padded / np.maximum(counts, 1)

    def at(self, residue: int) -> float:
        index = int(residue) - 1
        if not 0 <= index < len(self.observed):
            return float("nan")
        return float(self.local_rate()[index])


class GnomadClient:
    """Cached GraphQL client. Never requires the network twice for the same query."""

    def __init__(self, offline: bool = False, cache_dir=None) -> None:
        from ..config import CACHE_DIR

        self.offline = offline
        self.cache_dir = (cache_dir or CACHE_DIR / "gnomad")

    def _query(self, name: str, query: str) -> dict | None:
        path = self.cache_dir / f"{name}.json"
        if path.exists():
            return json.loads(path.read_text())
        if self.offline:
            return None
        request = urllib.request.Request(
            GNOMAD_API, data=json.dumps({"query": query}).encode(),
            headers={"Content-Type": "application/json",
                     "User-Agent": _USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                payload = json.loads(response.read())
        except (urllib.error.URLError, TimeoutError, OSError):
            return None
        if "errors" in payload and payload.get("errors"):
            return None
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload))
        return payload

    def constraint(self, symbol: str = "PIEZO1") -> GeneConstraint | None:
        query = ('{ gene(gene_symbol: "%s", reference_genome: GRCh38) { '
                 'gene_id symbol gnomad_constraint { pLI oe_lof oe_lof_upper '
                 'mis_z exp_mis obs_mis oe_mis } } }' % symbol)
        payload = self._query(f"constraint_{symbol}", query)
        if payload is None:
            return None
        gene = (payload.get("data") or {}).get("gene")
        if not gene or not gene.get("gnomad_constraint"):
            return None
        c = gene["gnomad_constraint"]
        return GeneConstraint(
            gene_id=gene["gene_id"], symbol=gene["symbol"], pli=c.get("pLI"),
            oe_lof=c.get("oe_lof"), loeuf=c.get("oe_lof_upper"),
            mis_z=c.get("mis_z"), obs_mis=c.get("obs_mis"),
            exp_mis=c.get("exp_mis"), oe_mis=c.get("oe_mis"),
            meta={"licence": GNOMAD_LICENCE, "citation": GNOMAD_CITATION})

    def variants(self, symbol: str = "PIEZO1",
                 dataset: str = "gnomad_r4") -> list | None:
        query = ('{ gene(gene_symbol: "%s", reference_genome: GRCh38) { '
                 'variants(dataset: %s) { consequence hgvsp pos '
                 'exome { ac an } genome { ac an } } } }' % (symbol, dataset))
        payload = self._query(f"variants_{symbol}_{dataset}", query)
        if payload is None:
            return None
        gene = (payload.get("data") or {}).get("gene")
        return None if not gene else gene.get("variants")


def _protein_position(hgvsp: str | None) -> int | None:
    """Residue number from an HGVS protein string like ``p.Arg2456His``."""
    if not hgvsp or "p." not in hgvsp:
        return None
    import re

    match = re.search(r"p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}", hgvsp)
    return int(match.group(1)) if match else None


def missense_density(variants: list, n_residues: int = 2521,
                     window: int = 25) -> MissenseDensity:
    """Per-residue observed missense counts, from the variant list.

    Counts *distinct alleles*, and separately sums allele counts. The two answer
    different questions — how many ways a position has been hit, versus how
    often — and a position can be high on one and low on the other.
    """
    observed = np.zeros(n_residues, dtype=float)
    allele = np.zeros(n_residues, dtype=float)
    unplaced = 0
    for variant in variants:
        if variant.get("consequence") != "missense_variant":
            continue
        position = _protein_position(variant.get("hgvsp"))
        if position is None or not 1 <= position <= n_residues:
            unplaced += 1
            continue
        observed[position - 1] += 1
        for source in ("exome", "genome"):
            block = variant.get(source) or {}
            allele[position - 1] += float(block.get("ac") or 0)
    return MissenseDensity(
        residues=np.arange(1, n_residues + 1), observed=observed,
        allele_count=allele, window=window,
        meta={"unplaced": unplaced, "licence": GNOMAD_LICENCE,
              "citation": GNOMAD_CITATION,
              "n_missense": int(observed.sum())})
