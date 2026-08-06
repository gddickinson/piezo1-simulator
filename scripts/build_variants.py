#!/usr/bin/env python
"""Promote the researched variant table into a committed resource.

Reads the curated JSON array embedded in ``ref/research/04_variants_disease.md``
and writes ``piezo1/resources/variants.json``.

The important step is the **validation gate**, not the copying. Every entry is
checked against the real Q92508 sequence, its mouse equivalent is computed by
alignment, and — critically — each variant is checked against every downloaded
human PIEZO1 structure to record whether the mutated residue is actually
modelled there. That last check matters because most of the interesting
variants sit in regions no experimental structure resolves, and a viewer that
silently fails to highlight a residue is worse than one that says why.

Usage::

    python scripts/build_variants.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import RESEARCH_DIR, RESOURCE_DIR, STRUCTURE_DIR  # noqa: E402
from piezo1.core.sequence import human_sequence, load_numbering_map, mouse_sequence  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402

SOURCE = RESEARCH_DIR / "04_variants_disease.md"

#: Human PIEZO1 entries whose residue coverage we report against.
HUMAN_STRUCTURES = ["8YEZ", "8ZU3", "8ZU8", "8YFC", "8YFG", "9VMX"]

_RESIDUE_RE = re.compile(r"^([A-Z])(\d+)")


def extract_json_array(text: str) -> list[dict]:
    """Pull the largest JSON array out of a markdown document."""
    best: list[dict] = []
    for match in re.finditer(r"```(?:json)?\s*(\[.*?\])\s*```", text, re.S):
        try:
            data = json.loads(match.group(1))
        except json.JSONDecodeError:
            continue
        if isinstance(data, list) and len(data) > len(best):
            best = data
    if not best:
        # Fall back to a bare array not wrapped in a fence.
        for match in re.finditer(r"(\[\s*\{.*?\}\s*\])", text, re.S):
            try:
                data = json.loads(match.group(1))
            except json.JSONDecodeError:
                continue
            if isinstance(data, list) and len(data) > len(best):
                best = data
    return best


def residue_coverage() -> dict[str, set[int]]:
    """Residues modelled in every protomer of each human structure."""
    out: dict[str, set[int]] = {}
    for pdb in HUMAN_STRUCTURES:
        path = STRUCTURE_DIR / f"{pdb}.cif"
        if not path.exists():
            continue
        st = Structure.from_file(path)
        per_chain = []
        for ch in st.chains:
            m = st.mask_ca() & (st.chain == ch)
            if m.sum() > 500:
                per_chain.append(set(st.res_seq[m].tolist()))
        if per_chain:
            common = set.intersection(*per_chain)
            out[pdb] = common
    return out


def normalise(entry: dict, hs: str, ms: str, nm, coverage: dict) -> dict | None:
    """Validate one raw entry and return the resource record."""
    pos = entry.get("residue_human_Q92508") or entry.get("residue_human")
    wt = entry.get("wt_aa") or entry.get("wt")
    mut = entry.get("mut_aa") or entry.get("mut")

    # Some entries encode the position only in the HGVS string.
    hgvs = entry.get("hgvs_protein") or entry.get("hgvs") or ""
    if pos in (None, "", "UNVERIFIED"):
        m = re.search(r"p\.\(?([A-Za-z]{1,3})(\d+)", str(hgvs))
        pos = int(m.group(2)) if m else None
    try:
        pos = int(pos)
    except (TypeError, ValueError):
        pos = None

    rec = {
        "id": entry.get("id") or hgvs or "unnamed",
        "hgvs_protein": hgvs,
        "residue": pos,
        "wt_aa": wt,
        "mut_aa": mut,
        "domain": entry.get("domain"),
        "classification": entry.get("classification"),
        "phenotype": entry.get("phenotype"),
        "functional_effect": entry.get("functional_effect")
                             or entry.get("functional_effect_summary"),
        "effect_magnitude": entry.get("effect_magnitude"),
        "pmid": entry.get("pmid"),
        "source_url": entry.get("source_url"),
        "confidence": entry.get("confidence", "unknown"),
        "structure_pdb": entry.get("structure_pdb"),
    }

    if pos is None or not (1 <= pos <= len(hs)):
        rec["sequence_verified"] = False
        rec["note"] = "no usable human residue number"
        rec["modelled_in"] = []
        return rec

    actual = hs[pos - 1]
    rec["human_aa_actual"] = actual
    rec["sequence_verified"] = (wt is None or wt == actual)
    mouse_pos = nm.to_b(pos)
    rec["mouse_residue"] = mouse_pos
    rec["mouse_aa"] = ms[mouse_pos - 1] if mouse_pos else None
    rec["offset"] = (mouse_pos - pos) if mouse_pos else None
    rec["conserved"] = (rec["mouse_aa"] == actual) if mouse_pos else None
    rec["modelled_in"] = sorted(p for p, res in coverage.items() if pos in res)
    return rec


def main() -> int:
    if not SOURCE.exists():
        raise SystemExit(f"{SOURCE} not found — run the variants research first")

    raw = extract_json_array(SOURCE.read_text())
    if not raw:
        raise SystemExit(f"no JSON array found in {SOURCE}")

    hs, ms = human_sequence(), mouse_sequence()
    nm = load_numbering_map()
    coverage = residue_coverage()

    records = [normalise(e, hs, ms, nm, coverage) for e in raw]
    records = [r for r in records if r]
    records.sort(key=lambda r: (r["residue"] or 0))

    bad = [r for r in records if not r.get("sequence_verified")]
    unmodelled = [r for r in records
                  if r["residue"] and not r["modelled_in"]]

    out = {
        "protein": "PIEZO1",
        "numbering": "human Q92508 (2521 aa)",
        "note": ("Every wild-type residue was checked against the Q92508 "
                 "sequence at build time. 'modelled_in' lists the human PIEZO1 "
                 "structures in which the residue is actually resolved in all "
                 "three protomers; an empty list means the viewer cannot show "
                 "it on that structure and must fall back to a predicted model."),
        "structures_checked": sorted(coverage),
        "coverage_summary": {p: len(r) for p, r in sorted(coverage.items())},
        "variants": records,
    }
    dest = RESOURCE_DIR / "variants.json"
    dest.write_text(json.dumps(out, indent=1))

    from collections import Counter
    classes = Counter(r.get("classification") for r in records)
    print(f"wrote {dest}")
    print(f"  {len(records)} variants: " +
          ", ".join(f"{k} {v}" for k, v in classes.most_common()))
    print(f"  sequence-verified: {len(records) - len(bad)}/{len(records)}")
    for r in bad:
        print(f"    ! {r['id']}: declared {r.get('wt_aa')}{r['residue']} but "
              f"Q92508 has {r.get('human_aa_actual')}")
    print(f"\n  residues resolved per human structure: "
          f"{ {p: len(c) for p, c in sorted(coverage.items())} }")
    print(f"  variants NOT modelled in any human structure: "
          f"{len(unmodelled)}/{len(records)}")
    shown = ", ".join(f"{r['wt_aa']}{r['residue']}" for r in unmodelled[:12]
                      if r.get("wt_aa"))
    if shown:
        print(f"    e.g. {shown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
