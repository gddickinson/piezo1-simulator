"""External variant predictors, via the ProtVar API.

ProtVar (EMBL-EBI) serves AlphaMissense, EVE, ESM-1b, per-position conservation
and **precomputed FoldX ΔΔG** for a UniProt accession and position, under
**CC BY 4.0**. That licence is the reason this route was chosen over running
the tools locally:

* **FoldX** is not redistributable at all — academic use requires a signed
  agreement, and the community Python wrapper has no licence file, meaning all
  rights reserved.
* **SIFT4G** is GPL-3.0, which would force the whole PyQt application to GPL.
* Every tool on ``biosig.lab.uq.edu.au`` (mCSM, DynaMut2 and relatives) carries
  **no licence text whatsoever**, so nothing is granted.
* **VarSite** and **VarMap** are both retired.

Going through ProtVar sidesteps all of that: one CC BY 4.0 source, no local
models, no licence traps. Attribution is recorded in every cached response.

**What these predictors can and cannot say.** AlphaMissense, EVE and ESM-1b all
emit a *pathogenicity* score — a single axis running benign to damaging. That
axis has no room to express **direction**: a variant that opens a channel too
easily and one that stops it opening at all are both "damaging". They are
therefore complementary to, not a replacement for, the mechanical features this
project computes, and Round 22 exists to test whether the combination does
better than either. Nothing here is compared against the variant labels yet.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from ..config import CACHE_DIR

__all__ = ["ProtVarClient", "ExternalScores", "PROTVAR_BASE", "PROTVAR_LICENCE",
           "PROTVAR_CITATION"]

PROTVAR_BASE = "https://www.ebi.ac.uk/ProtVar/api"
PROTVAR_LICENCE = "CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)"
PROTVAR_CITATION = ("ProtVar, EMBL-EBI. https://www.ebi.ac.uk/ProtVar/ — "
                    "serves AlphaMissense, EVE, ESM-1b, conservation and "
                    "precomputed FoldX under CC BY 4.0.")
USER_AGENT = "piezo1-simulator/0.1 (research use)"


@dataclass
class ExternalScores:
    """Predictor output for one substitution."""

    accession: str
    position: int
    wild_type: str | None = None
    mutant: str | None = None

    conservation: float | None = None
    alphamissense: float | None = None
    alphamissense_class: str | None = None
    eve: float | None = None
    eve_class: str | None = None
    esm1b: float | None = None
    foldx_ddg: float | None = None
    plddt: float | None = None
    missense3d: str | None = None

    source: str = "protvar"
    licence: str = PROTVAR_LICENCE
    cached: bool = False
    meta: dict = field(default_factory=dict)

    @property
    def available(self) -> list[str]:
        return [k for k in ("conservation", "alphamissense", "eve", "esm1b",
                            "foldx_ddg") if getattr(self, k) is not None]

    def as_dict(self) -> dict:
        return {
            "accession": self.accession, "position": self.position,
            "wild_type": self.wild_type, "mutant": self.mutant,
            "conservation": self.conservation,
            "alphamissense": self.alphamissense,
            "alphamissense_class": self.alphamissense_class,
            "eve": self.eve, "eve_class": self.eve_class,
            "esm1b": self.esm1b, "foldx_ddg": self.foldx_ddg,
            "plddt": self.plddt, "missense3d": self.missense3d,
            "source": self.source, "licence": self.licence,
        }


class ProtVarClient:
    """Cached, offline-tolerant client for the ProtVar API.

    Every response is written to disk, so a run that has been done once works
    without a network afterwards. When a request fails and nothing is cached,
    the client returns ``None`` rather than raising: a missing external score
    should degrade an analysis, not abort it.
    """

    def __init__(self, cache_dir: Path | None = None, timeout: int = 45,
                 offline: bool = False, delay: float = 0.0) -> None:
        self.cache_dir = Path(cache_dir or (CACHE_DIR / "protvar"))
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.offline = offline
        self.delay = delay
        self.stats = {"hit": 0, "miss": 0, "fail": 0}

    # ------------------------------------------------------------ transport

    def _cache_path(self, key: str) -> Path:
        return self.cache_dir / (key.replace("/", "_").replace("?", "__")
                                 .replace("=", "-") + ".json")

    def _get(self, endpoint: str) -> object | None:
        path = self._cache_path(endpoint)
        if path.exists():
            self.stats["hit"] += 1
            try:
                return json.loads(path.read_text())
            except json.JSONDecodeError:
                path.unlink(missing_ok=True)
        if self.offline:
            self.stats["fail"] += 1
            return None

        url = f"{PROTVAR_BASE}/{endpoint}"
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as response:
                data = json.loads(response.read())
        except (urllib.error.URLError, urllib.error.HTTPError,
                TimeoutError, json.JSONDecodeError):
            self.stats["fail"] += 1
            return None
        self.stats["miss"] += 1
        path.write_text(json.dumps(data))
        if self.delay:
            time.sleep(self.delay)
        return data

    # -------------------------------------------------------------- queries

    def scores(self, accession: str, position: int,
               mutant: str | None = None) -> ExternalScores | None:
        """Predictor scores for one position, optionally one substitution.

        **Pass ``mutant``.** Without it the endpoint returns nineteen entries
        per predictor — one per possible substitution — in an order the payload
        never states, so there is no safe way to tell which score belongs to
        which mutation. The ``mt`` query parameter resolves it properly;
        guessing an alphabetical ordering would be a silent correctness bug.
        """
        endpoint = f"score/{accession}/{position}"
        if mutant:
            endpoint += f"?mt={urllib.parse.quote(mutant)}"
        payload = self._get(endpoint)
        if payload is None:
            return None
        if mutant is None:
            # Only the position-level conservation is unambiguous here.
            payload = [x for x in payload if x.get("type") == "CONSERV"]

        out = ExternalScores(accession=accession, position=int(position),
                             mutant=mutant,
                             cached=self._cache_path(endpoint).exists())
        for entry in payload:
            kind = entry.get("type")
            if kind == "CONSERV":
                out.conservation = entry.get("score")
            elif kind == "AM":
                out.alphamissense = entry.get("amPathogenicity")
                out.alphamissense_class = entry.get("amClass")
            elif kind == "EVE":
                out.eve = entry.get("score")
                out.eve_class = entry.get("eveClass")
            elif kind == "ESM":
                out.esm1b = entry.get("score")
            elif kind == "M3D":
                out.missense3d = entry.get("prediction")
        return out

    def foldx(self, accession: str, position: int) -> dict[str, dict]:
        """Precomputed FoldX ΔΔG for every substitution at a position.

        Keyed by mutant residue. Unlike ``/score``, this endpoint labels each
        entry with its ``mutatedType``, so no disambiguation is needed.
        """
        payload = self._get(f"prediction/foldx/{accession}/{position}")
        if not payload:
            return {}
        return {entry["mutatedType"]: {"ddg": entry.get("foldxDdg"),
                                       "plddt": entry.get("plddt"),
                                       "wild_type": entry.get("wildType")}
                for entry in payload if entry.get("mutatedType")}

    def pockets(self, accession: str, position: int) -> list[dict]:
        """Predicted pockets containing a residue, from the AlphaFold model."""
        return list(self._get(f"prediction/pocket/{accession}/{position}") or [])

    def function(self, accession: str, position: int) -> object | None:
        return self._get(f"function/{accession}/{position}")

    # ------------------------------------------------------------ enrichment

    def annotate(self, accession: str, position: int, wild_type: str | None,
                 mutant: str | None) -> ExternalScores | None:
        """Everything available for one variant, in one object."""
        result = self.scores(accession, position, mutant)
        if result is None:
            return None
        result.wild_type = wild_type
        if mutant:
            entry = self.foldx(accession, position).get(mutant)
            if entry:
                result.foldx_ddg = entry.get("ddg")
                result.plddt = entry.get("plddt")
                if wild_type and entry.get("wild_type") and \
                        entry["wild_type"] != wild_type:
                    result.meta["wild_type_mismatch"] = (
                        f"caller said {wild_type}, ProtVar says "
                        f"{entry['wild_type']} at {position}")
        return result

    def annotate_variants(self, variants, accession: str = "Q92508",
                          verbose: bool = False) -> dict[str, ExternalScores]:
        """Annotate an iterable of :class:`piezo1.core.annotations.Variant`."""
        out: dict[str, ExternalScores] = {}
        for variant in variants:
            if variant.residue is None:
                continue
            wt = variant.wt_aa if variant.wt_aa and len(variant.wt_aa) == 1 else None
            mt = variant.mut_aa if variant.mut_aa and len(variant.mut_aa) == 1 else None
            if mt is None:
                continue
            result = self.annotate(accession, variant.residue, wt, mt)
            if result is not None:
                out[variant.label] = result
            if verbose:
                got = ",".join(result.available) if result else "unavailable"
                print(f"   {variant.label:10s} {got}")
        return out
