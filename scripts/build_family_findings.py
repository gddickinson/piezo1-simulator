#!/usr/bin/env python
"""Import the ``piezo_genes`` census findings behind a provenance gate.

Writes two committed resources:

``resources/family_findings.json``
    The thirteen statements of :mod:`family_table`, each with the numbers it
    rests on, the file those numbers came from, what this project does with it,
    and what it does not establish.

``resources/family_constraint.json``
    The per-residue evolutionary constraint tracks for human PIEZO1, human
    PIEZO2 and zebrafish piezo3 — the one bulk import, because a per-residue
    scalar is the thing a structural project can actually *use*.

**The gate.** An imported finding is somebody else's result, and the failure
mode is specific: the source project keeps working, corrects a number, and this
copy silently becomes a confident quotation of a superseded value. So every
number in the table declares how to re-read it, and the build:

1. re-reads each ``table`` check from the source TSV and compares it, failing
   on any disagreement past the stated tolerance;
2. requires each ``document`` check to appear literally in the source project's
   own ``FINDINGS.md``, which is how claims never reduced to a table are held to
   the same standard;
3. refuses to write at all if the source project is not on disk — a rebuild
   with no source would otherwise re-stamp the existing resource with a fresh
   date and no verification behind it;
4. records the source project's git commit in both resources, so a reader can
   tell which version of the census this is a copy of;
5. checks the per-residue tracks against the sequence this project already
   holds — every residue's amino acid must match ``uniprot_human.json``, or the
   track is in a numbering we do not think it is in.

Step 5 is the one that matters most, and it is not decoration: the census works
in three numbering systems and this project works in two others. A track joined
by residue number to the wrong sequence would produce a per-residue colouring
that looks entirely plausible and is off by an indel.

Usage::

    python scripts/build_family_findings.py
    python scripts/build_family_findings.py --check          # verify, write nothing
    python scripts/build_family_findings.py --source ../piezo_genes
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from family_table import (CONSTRAINT_GENES, FINDINGS, SOURCE_PROJECT,  # noqa: E402
                          TABLE_IMPORTS)
from family_table_clinical import (CENSUS, EQUIVALENT_POSITIONS,  # noqa: E402
                                   PORE_MODULE_PATHOGENIC)

from piezo1.config import RESOURCE_DIR  # noqa: E402

#: Fractional tolerance on a re-read number. The table quotes the source's own
#: rounding, so this only has to absorb the last printed digit.
TOLERANCE = 5e-4

#: Per-residue columns carried across. ``deep`` is the headline track — one
#: locus per genome across the paralogue's own orthologues — and the other two
#: are kept because a finding that only holds in one layer is worth being able
#: to see fail.
TRACKS = ("deep_jsd", "vert_jsd", "family_jsd")


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _as_number(text: str):
    try:
        return int(text)
    except (TypeError, ValueError):
        try:
            return float(text)
        except (TypeError, ValueError):
            return None


def _match(rows: list[dict], selector: dict) -> dict | None:
    for row in rows:
        if all(row.get(k) == v for k, v in selector.items()):
            return row
    return None


def verify(results: Path, document: str) -> list[str]:
    """Re-read every declared number. Returns a list of failures, empty if none."""
    failures: list[str] = []
    cache: dict[str, list[dict]] = {}
    for finding in FINDINGS:
        key = finding["key"]
        for check in finding["check"]:
            if check[0] == "document":
                if check[1] not in document:
                    failures.append(
                        f"{key}: '{check[1]}' does not appear in the source "
                        f"FINDINGS.md")
                continue
            _, selector, column, expected = check
            src = finding["source"]
            path = results / src
            if not path.exists():
                failures.append(f"{key}: source {src} missing")
                continue
            rows = cache.setdefault(src, read_tsv(path))
            row = _match(rows, selector)
            if row is None:
                failures.append(f"{key}: no row in {src} matching {selector}")
                continue
            got = _as_number(row.get(column, ""))
            if got is None:
                failures.append(f"{key}: {src} column {column!r} is not a number")
                continue
            scale = max(abs(expected), 1.0)
            if abs(got - expected) > TOLERANCE * scale:
                failures.append(
                    f"{key}: {src}[{selector}].{column} is {got}, table says "
                    f"{expected} — the census has moved, or this copy is stale")
    return failures


def _source_commit(source: Path) -> str:
    try:
        out = subprocess.run(["git", "-C", str(source), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=15)
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _human_sequence() -> str:
    data = json.loads((RESOURCE_DIR / "uniprot_human.json").read_text())
    return data["sequence"]


def build_constraint(results: Path) -> dict:
    """The per-residue tracks, checked against the sequence we already hold."""
    human = _human_sequence()
    genes: dict[str, dict] = {}
    problems: list[str] = []
    for spec in CONSTRAINT_GENES:
        path = results / spec["source"]
        if not path.exists():
            problems.append(f"{spec['gene']}: {spec['source']} missing")
            continue
        rows = read_tsv(path)
        length = spec["length"]
        by_resi = {int(r["resi"]): r for r in rows}
        if len(by_resi) != length:
            problems.append(
                f"{spec['gene']}: {len(by_resi)} rows for a {length}-residue "
                f"protein")
        aa = [by_resi[i]["aa"] if i in by_resi else "" for i in range(1, length + 1)]
        if spec["accession"] == "Q92508":
            mismatch = sum(1 for i, a in enumerate(aa) if a and a != human[i])
            if mismatch:
                problems.append(
                    f"{spec['gene']}: {mismatch} residues disagree with "
                    f"uniprot_human.json — the track is not in the numbering "
                    f"it claims")
        tracks = {}
        for track in TRACKS:
            tracks[track] = [
                (round(float(by_resi[i][track]), 4)
                 if i in by_resi and by_resi[i].get(track) not in (None, "", "NA")
                 else None)
                for i in range(1, length + 1)
            ]
        genes[spec["gene"]] = {
            "accession": spec["accession"],
            "numbering": spec["numbering"],
            "length": length,
            "n_orthologues": spec["n_orthologues"],
            "sequence": "".join(aa),
            "reliable": [i in by_resi and by_resi[i].get("deep_reliable") == "True"
                         for i in range(1, length + 1)],
            "in_pore_module": [i in by_resi and by_resi[i].get("in_pore_module") == "True"
                               for i in range(1, length + 1)],
            "census_domain": [by_resi[i]["domain"] if i in by_resi else None
                              for i in range(1, length + 1)],
            **tracks,
        }
    return {"genes": genes, "problems": problems}


def build_tables(results: Path) -> tuple[dict, list[str]]:
    tables, problems = {}, []
    for spec in TABLE_IMPORTS:
        path = results / spec["source"]
        if not path.exists():
            problems.append(f"{spec['key']}: {spec['source']} missing")
            continue
        rows = read_tsv(path)
        kept = [{c: row.get(c) for c in spec["columns"]} for row in rows]
        tables[spec["key"]] = {"note": spec["note"], "source": spec["source"],
                               "rows": kept}
    return tables, problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default=None,
                        help="path to the piezo_genes project")
    parser.add_argument("--check", action="store_true",
                        help="verify against the source and write nothing")
    args = parser.parse_args()

    default = Path(__file__).resolve().parent.parent / SOURCE_PROJECT["default_path"]
    source = (Path(args.source).expanduser() if args.source else default).resolve()
    results = source / SOURCE_PROJECT["results_subdir"]

    if not results.is_dir():
        print(f"REFUSED: the census project is not on disk at {source}.\n"
              f"         Nothing is written. The committed resource stays as it "
              f"is,\n         which is the honest outcome: it cannot be "
              f"re-verified from here.", file=sys.stderr)
        return 2

    findings_md = source / "FINDINGS.md"
    document = findings_md.read_text() if findings_md.exists() else ""
    failures = verify(results, document)
    constraint = build_constraint(results)
    tables, table_problems = build_tables(results)
    failures += constraint["problems"] + table_problems

    if failures:
        print("REFUSED — the import does not verify against its source:",
              file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1

    commit = _source_commit(source)
    provenance = {
        "source_project": SOURCE_PROJECT["name"],
        "source_title": SOURCE_PROJECT["title"],
        "source_commit": commit,
        "source_note": SOURCE_PROJECT["note"],
        "built_by": "scripts/build_family_findings.py",
        "verified": (f"{len(FINDINGS)} findings, "
                     f"{sum(len(f['check']) for f in FINDINGS)} numbers re-read "
                     f"from the source on this build"),
    }

    findings_doc = {
        "provenance": provenance,
        "findings": [
            dict({k: v for k, v in f.items() if k != "check"},
                 n_checks=len(f["check"]))
            for f in FINDINGS
        ],
        "tables": tables,
        "pore_module_pathogenic": [
            {"gene": g, "resi": r, "aa": a, "element": e, "msa_col": c,
             "piezo3_resi": p3, "piezo3_aa": p3aa}
            for g, r, a, e, c, p3, p3aa in PORE_MODULE_PATHOGENIC
        ],
        "equivalent_positions": list(EQUIVALENT_POSITIONS),
        "census": CENSUS,
    }

    constraint_doc = {"provenance": provenance,
                      "tracks": list(TRACKS),
                      "genes": constraint["genes"]}

    if args.check:
        print(f"OK — {len(FINDINGS)} findings verify against {source} "
              f"({commit}); nothing written.")
        return 0

    (RESOURCE_DIR / "family_findings.json").write_text(
        json.dumps(findings_doc, indent=1) + "\n")
    (RESOURCE_DIR / "family_constraint.json").write_text(
        json.dumps(constraint_doc) + "\n")
    print(f"wrote family_findings.json ({len(FINDINGS)} findings, "
          f"{len(tables)} tables) and family_constraint.json "
          f"({len(constraint['genes'])} genes) from {source} ({commit})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
