#!/usr/bin/env python
"""Expand the phenotyped variant set from ClinVar, behind the usual gates.

Round 20 established that the binding constraint on this project's central
claim is **data, not method**: 42 variants are needed for a large effect and 98
for a medium one, against the 25 that survived Round 7's inclusion criteria.

**The problem ClinVar does not solve on its own.** It reports *pathogenicity*,
not *direction*. A gain-of-function and a loss-of-function variant are both
"Pathogenic", and this project needs to tell them apart. What makes direction
recoverable for PIEZO1 is that the two diseases have opposite, well-established
mechanisms:

* **Dehydrated hereditary stomatocytosis / xerocytosis** — dominant, and the
  mechanism is *slowed inactivation*, i.e. gain of function (Zarychanski 2012;
  Albuisson 2013; Andolfo 2013).
* **Generalised lymphatic dysplasia / lymphatic malformation** — recessive, and
  the mechanism is loss of function (Fotiou 2015).

So a condition can imply a direction. **That is weaker evidence than measuring
the current**, and the difference is recorded per variant rather than averaged
away: every entry carries an ``evidence`` level, and a downstream analysis can
restrict to the measured ones if it wants to.

Gates, all inherited from `build_variants.py`:

* the wild-type residue must match UniProt Q92508 at that position;
* the protein change must parse to a single, unambiguous substitution or a
  clearly typed truncation;
* anything failing either is **reported, not silently dropped**.

Usage::

    python scripts/build_variants_clinvar.py            # fetch and report
    python scripts/build_variants_clinvar.py --write    # also write the resource
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import CACHE_DIR, RESOURCE_DIR  # noqa: E402
from piezo1.core.sequence import human_sequence  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "piezo1-simulator/0.1 (research use)"

#: Condition text → implied direction, with the paper establishing the
#: mechanism. Matched case-insensitively as a substring. Deliberately short:
#: a condition not listed here yields **no direction**, which is the honest
#: default rather than a guess.
CONDITION_DIRECTION = {
    "dehydrated hereditary stomatocytosis": ("GoF", "zarychanski2012"),
    "hereditary xerocytosis": ("GoF", "zarychanski2012"),
    "stomatocytosis": ("GoF", "albuisson2013"),
    "lymphatic malformation": ("LoF", "fotiou2015"),
    "generalized lymphatic dysplasia": ("LoF", "fotiou2015"),
    "lymphedema": ("LoF", "fotiou2015"),
    "hydrops fetalis": ("LoF", "fotiou2015"),
}

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}

PROTEIN_CHANGE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|fs|=)")


def _get(url: str, cache_name: str) -> dict:
    """Fetch with a disk cache, so a rerun needs no network."""
    path = CACHE_DIR / "clinvar" / cache_name
    if path.exists():
        return json.loads(path.read_text())
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        data = json.loads(response.read())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    time.sleep(0.4)
    return data


def fetch_records(limit: int = 2000) -> list[dict]:
    """Pathogenic and likely-pathogenic PIEZO1 variants, with their conditions."""
    term = ('PIEZO1[gene] AND (("pathogenic"[Clinical_significance]) OR '
            '("likely pathogenic"[Clinical_significance]))')
    search = _get(
        f"{EUTILS}/esearch.fcgi?db=clinvar&term={urllib.parse.quote(term)}"
        f"&retmax={limit}&retmode=json", "search.json")
    ids = search["esearchresult"]["idlist"]

    records = []
    for start in range(0, len(ids), 200):
        chunk = ids[start:start + 200]
        summary = _get(
            f"{EUTILS}/esummary.fcgi?db=clinvar&id={','.join(chunk)}"
            f"&retmode=json", f"summary_{start}.json")["result"]
        for uid in summary.get("uids", []):
            entry = summary[uid]
            germline = entry.get("germline_classification") or {}
            traits = [t.get("trait_name", "")
                      for t in germline.get("trait_set") or []]
            records.append({
                "uid": uid,
                "accession": entry.get("accession", ""),
                "title": entry.get("title", ""),
                "protein_change": entry.get("protein_change", ""),
                "significance": germline.get("description", ""),
                "review_status": germline.get("review_status", ""),
                "conditions": traits,
            })
    return records


def parse_change(record: dict) -> dict | None:
    """Extract (wild type, position, mutant) from the HGVS protein change."""
    match = PROTEIN_CHANGE.search(record["title"])
    if not match:
        return None
    wt3, position, mut3 = match.groups()
    wt = THREE_TO_ONE.get(wt3)
    if wt is None:
        return None
    if mut3 == "fs":
        kind, mut = "frameshift", "fs"
    elif mut3 == "Ter":
        kind, mut = "nonsense", "*"
    elif mut3 == "=":
        kind, mut = "synonymous", wt
    else:
        mut = THREE_TO_ONE.get(mut3)
        if mut is None:
            return None
        kind = "missense"
    return {"wt_aa": wt, "residue": int(position), "mut_aa": mut, "kind": kind}


def direction_of(conditions: list[str]) -> tuple[str | None, str, list[str]]:
    """Implied direction, the citation for the mechanism, and any conflict."""
    votes: dict[str, str] = {}
    for condition in conditions:
        lowered = condition.lower()
        for needle, (direction, citation) in CONDITION_DIRECTION.items():
            if needle in lowered:
                votes[direction] = citation
    if not votes:
        return None, "", []
    if len(votes) > 1:
        # A variant reported under both a gain- and a loss-of-function disease.
        # Recorded as ambiguous rather than resolved by preferring one.
        return None, "", sorted(votes)
    direction = next(iter(votes))
    return direction, votes[direction], []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="write piezo1/resources/variants_clinvar.json")
    args = parser.parse_args()

    reference = human_sequence()
    records = fetch_records()
    print(f"fetched {len(records)} pathogenic/likely-pathogenic ClinVar records")

    kept, rejected = [], Counter()
    ambiguous, wt_mismatch = [], []
    for record in records:
        parsed = parse_change(record)
        if parsed is None:
            rejected["unparseable protein change"] += 1
            continue
        if parsed["kind"] == "synonymous":
            rejected["synonymous"] += 1
            continue

        position = parsed["residue"]
        if not (1 <= position <= len(reference)):
            rejected["position outside Q92508"] += 1
            continue
        if reference[position - 1] != parsed["wt_aa"]:
            wt_mismatch.append(
                f"{parsed['wt_aa']}{position}: Q92508 has "
                f"{reference[position - 1]}")
            rejected["wild type disagrees with Q92508"] += 1
            continue

        direction, citation, conflict = direction_of(record["conditions"])
        if conflict:
            ambiguous.append(f"{parsed['wt_aa']}{position}{parsed['mut_aa']}: "
                             f"{'/'.join(conflict)} — {record['conditions']}")
        kept.append({
            **parsed,
            "label": f"{parsed['wt_aa']}{position}{parsed['mut_aa']}",
            "classification": direction,
            "evidence": "disease_mechanism" if direction else "none",
            "mechanism_citation": citation,
            "significance": record["significance"],
            "review_status": record["review_status"],
            "conditions": record["conditions"],
            "accession": record["accession"],
            "source": "clinvar",
            "ambiguous_direction": bool(conflict),
        })

    print(f"\nparsed and wild-type verified: {len(kept)}")
    for reason, count in rejected.most_common():
        print(f"  rejected — {reason}: {count}")
    if wt_mismatch:
        print(f"\n  wild-type mismatches (first 5): {wt_mismatch[:5]}")

    by_kind = Counter(v["kind"] for v in kept)
    directed = [v for v in kept if v["classification"]]
    by_direction = Counter(v["classification"] for v in directed)
    missense_directed = [v for v in directed if v["kind"] == "missense"]

    print(f"\nby consequence: {dict(by_kind)}")
    print(f"with an implied direction: {len(directed)} {dict(by_direction)}")
    print(f"  of those, missense: {len(missense_directed)} "
          f"{dict(Counter(v['classification'] for v in missense_directed))}")
    if ambiguous:
        print(f"\n  AMBIGUOUS — reported under both directions ({len(ambiguous)}):")
        for line in ambiguous[:8]:
            print(f"    {line}")

    if args.write:
        dest = RESOURCE_DIR / "variants_clinvar.json"
        dest.write_text(json.dumps({
            "schema": 1,
            "source": "ClinVar via NCBI E-utilities",
            "note": ("Direction is INFERRED from the disease mechanism, not "
                     "measured. Weaker evidence than electrophysiology; the "
                     "`evidence` field records which."),
            "condition_map": {k: list(v) for k, v in CONDITION_DIRECTION.items()},
            "variants": kept,
        }, indent=1))
        print(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
