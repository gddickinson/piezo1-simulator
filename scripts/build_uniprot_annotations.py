#!/usr/bin/env python
"""Distil the UniProt entry for PIEZO1 into a compact annotation resource.

Reads the cached UniProt JSON in ``ref/sequences/`` and writes
``piezo1/resources/uniprot_<species>.json``, which *is* committed: it is small,
and pinning it means the application's residue annotations cannot silently
change under us when UniProt revises the entry.

Usage::

    python scripts/build_uniprot_annotations.py            # human + mouse
    python scripts/build_uniprot_annotations.py --fetch    # re-download first
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import (FLY_PIEZO_ACC, HUMAN_ACC,  # noqa: E402
                           HUMAN_PIEZO2_ACC, MOUSE_ACC, MOUSE_PIEZO2_ACC,
                           RESOURCE_DIR, SEQUENCE_DIR, WORM_PIEZO_ACC)

#: The PIEZO2 entries are here because 6KG7 is a PIEZO2 structure and the
#: paralogue comparison needs its transmembrane annotation from the same source
#: PIEZO1's comes from — otherwise the two dome measurements would differ by
#: how their membrane surfaces were defined rather than by their shape. Mouse
#: as well as human because 6KG7 is deposited in mouse numbering.
SPECIES = {"human": HUMAN_ACC, "mouse": MOUSE_ACC,
           "human_piezo2": HUMAN_PIEZO2_ACC, "mouse_piezo2": MOUSE_PIEZO2_ACC,
           "worm_piezo": WORM_PIEZO_ACC, "fly_piezo": FLY_PIEZO_ACC}

#: UniProt feature types we keep, mapped to the key used in the output file.
KEEP_RANGES = {
    "Transmembrane": "transmembrane",
    "Topological domain": "topology",
    "Region": "regions",
    "Coiled coil": "coiled_coil",
    "Disulfide bond": "disulfide",
}
KEEP_SITES = {
    "Modified residue": "modified_residues",
    "Glycosylation": "glycosylation",
    "Natural variant": "natural_variants",
    "Mutagenesis": "mutagenesis",
}


def fetch(acc: str, dest: Path) -> None:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    print(f"  fetching {url}")
    with urllib.request.urlopen(url, timeout=120) as r:
        dest.write_bytes(r.read())


def _span(feature: dict) -> tuple[int, int]:
    loc = feature["location"]
    return int(loc["start"]["value"]), int(loc["end"]["value"])


def distil(entry: dict) -> dict:
    out: dict = {
        "accession": entry["primaryAccession"],
        "id": entry.get("uniProtkbId"),
        "organism": entry.get("organism", {}).get("scientificName"),
        "gene": (entry.get("genes") or [{}])[0].get("geneName", {}).get("value"),
        "length": entry["sequence"]["length"],
        "mass_da": entry["sequence"]["molWeight"],
        "sequence": entry["sequence"]["value"],
        "entry_version": entry.get("entryAudit", {}).get("entryVersion"),
        "sequence_version": entry.get("entryAudit", {}).get("sequenceVersion"),
        "last_modified": entry.get("entryAudit", {}).get("lastSequenceUpdateDate"),
    }
    for key in list(KEEP_RANGES.values()) + list(KEEP_SITES.values()):
        out[key] = []

    for f in entry.get("features", []):
        ftype = f["type"]
        start, end = _span(f)
        desc = f.get("description", "")
        if ftype in KEEP_RANGES:
            out[KEEP_RANGES[ftype]].append(
                {"start": start, "end": end, "description": desc}
            )
        elif ftype in KEEP_SITES:
            rec = {"position": start, "end": end, "description": desc}
            alt = f.get("alternativeSequence")
            if alt:
                rec["wt"] = alt.get("originalSequence", "")
                alts = alt.get("alternativeSequences") or [""]
                rec["mut"] = alts[0]
            xrefs = [x["id"] for x in f.get("featureCrossReferences", [])
                     if x.get("database") == "dbSNP"]
            if xrefs:
                rec["dbsnp"] = xrefs[0]
            out[KEEP_SITES[ftype]].append(rec)

    # Number the transmembrane helices in sequence order: TM1..TM38.
    out["transmembrane"].sort(key=lambda d: d["start"])
    for i, tm in enumerate(out["transmembrane"], start=1):
        tm["name"] = f"TM{i}"

    out["n_transmembrane"] = len(out["transmembrane"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fetch", action="store_true", help="re-download from UniProt")
    args = ap.parse_args()

    SEQUENCE_DIR.mkdir(parents=True, exist_ok=True)
    RESOURCE_DIR.mkdir(parents=True, exist_ok=True)

    for species, acc in SPECIES.items():
        src = SEQUENCE_DIR / f"{acc}_{species}_PIEZO1.json"
        if args.fetch or not src.exists():
            fetch(acc, src)
        entry = json.loads(src.read_text())
        data = distil(entry)
        dest = RESOURCE_DIR / f"uniprot_{species}.json"
        dest.write_text(json.dumps(data, indent=1))
        print(f"{species:6s} {acc}: {data['length']} aa, "
              f"{data['n_transmembrane']} TM, "
              f"{len(data['natural_variants'])} natural variants -> {dest.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
