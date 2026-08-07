#!/usr/bin/env python
"""Author ``piezo1/resources/ligands.json`` behind a provenance gate.

Same shape as ``build_parameters.py`` and ``build_variants.py``: authored
content, validated on the way out, committed as a resource.

**The gate this one enforces.** No PIEZO structure with a bound small-molecule
modulator has ever been deposited, so no binding site here may claim one. The
build:

1. refuses any ``site_evidence`` of ``bound_structure``;
2. **verifies** that claim against the downloaded structures rather than
   trusting it — if a modulator ever appears as a heteroatom in a deposited
   entry, the build fails and this file is out of date;
3. requires every citation to resolve in ``references.json``;
4. fetches chemistry from PubChem and checks the returned InChIKey matches the
   one recorded here, so a wrong CID cannot pass silently;
5. requires a site with residues to carry a citation, and a site without
   residues to say why.

Usage::

    python scripts/build_ligands.py
    python scripts/build_ligands.py --check      # validate without writing
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ligand_table import LIGANDS, SITE_EVIDENCE  # noqa: E402

from piezo1.config import CACHE_DIR, RESOURCE_DIR  # noqa: E402

PUBCHEM = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
           "/property/MolecularFormula,MolecularWeight,CanonicalSMILES,"
           "InChIKey/JSON")
_UA = "piezo1-simulator/0.1 (research use)"

#: Heteroatom codes that are lipid, detergent, glycan or ion — i.e. everything
#: legitimately present in a deposited PIEZO entry. Anything outside this set
#: appearing in a structure would be a candidate bound modulator and is worth a
#: human look, which is why the check reports rather than guesses.
EXPECTED_HETERO = {"D12", "L9Q", "P5S", "PEE", "PLX", "NAG", "CL", "HOH",
                   "NA", "K", "MG", "CA", "ZN", "SO4", "PC1", "POV", "CLR"}


def fetch_chemistry(cid: int) -> dict | None:
    """PubChem properties for one CID, cached to disk."""
    path = CACHE_DIR / "pubchem" / f"cid_{cid}.json"
    if path.exists():
        return json.loads(path.read_text())
    request = urllib.request.Request(PUBCHEM.format(cid=cid),
                                     headers={"User-Agent": _UA})
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
            payload = json.loads(response.read())
    except Exception:
        return None
    properties = payload.get("PropertyTable", {}).get("Properties", [])
    if not properties:
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(properties[0]))
    return properties[0]


def deposited_modulators() -> list:
    """Heteroatom codes in downloaded structures that are not expected.

    Verifies the central claim rather than asserting it. If a Yoda1-bound
    structure is ever deposited and downloaded, this returns its ligand code and
    the build fails — which is the correct outcome, because the resource would
    then be wrong.
    """
    import numpy as np

    from piezo1.core.structure import Structure
    from piezo1.io.registry import load_registry

    found = {}
    for entry in load_registry().entries:
        if not entry.available:
            continue
        structure = Structure.from_file(entry.path)
        mask = structure.hetero
        if not mask.any():
            continue
        for code in np.unique(structure.res_name[mask]):
            if str(code) not in EXPECTED_HETERO:
                found.setdefault(str(code), []).append(entry.pdb)
    return sorted(found.items())


def validate(ligands: list, references: set) -> list:
    problems = []
    for ligand in ligands:
        key = ligand["key"]
        if ligand["site_evidence"] not in SITE_EVIDENCE:
            problems.append(f"{key}: unknown site_evidence "
                            f"{ligand['site_evidence']!r}")
        if ligand["site_evidence"] == "bound_structure":
            problems.append(f"{key}: claims a bound structure; none exists")
        if ligand["site_residues"] and not ligand.get("site_citation"):
            problems.append(f"{key}: names residues without a citation")
        if not ligand["site_residues"] and not ligand.get("site_note"):
            problems.append(f"{key}: no site and no explanation of why")

        potency = ligand.get("potency")
        if potency:
            if potency["citation"] not in references:
                problems.append(f"{key}: potency cites {potency['citation']!r}, "
                                f"which is not in references.json")
            if potency["value"] <= 0:
                problems.append(f"{key}: non-positive potency")
        if ligand.get("site_citation") and ligand["site_citation"] not in references:
            problems.append(f"{key}: site cites {ligand['site_citation']!r}, "
                            f"which is not in references.json")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="validate without writing")
    args = parser.parse_args()

    references = {r["key"] for r in json.loads(
        (RESOURCE_DIR / "references.json").read_text())["references"]}

    problems = validate(LIGANDS, references)

    print("verifying the central claim against the downloaded structures…")
    unexpected = deposited_modulators()
    if unexpected:
        print("  UNEXPECTED heteroatoms found:")
        for code, entries in unexpected:
            print(f"    {code}: {entries}")
        problems.append(
            "unexpected heteroatoms in deposited structures — one may be a "
            "bound modulator, in which case this resource is out of date")
    else:
        print("  no unexpected heteroatoms: no bound modulator in any of the "
              "downloaded entries, as recorded")

    records = []
    for ligand in LIGANDS:
        record = dict(ligand)
        record["site_residues"] = list(ligand["site_residues"])
        cid = ligand.get("pubchem_cid")
        if cid:
            chemistry = fetch_chemistry(cid)
            if chemistry is None:
                problems.append(f"{ligand['key']}: PubChem CID {cid} "
                                f"unreachable and not cached")
            else:
                if chemistry.get("InChIKey") != ligand["inchikey"]:
                    problems.append(
                        f"{ligand['key']}: CID {cid} returns InChIKey "
                        f"{chemistry.get('InChIKey')!r}, recorded "
                        f"{ligand['inchikey']!r}")
                record["formula"] = chemistry.get("MolecularFormula")
                record["molecular_weight"] = chemistry.get("MolecularWeight")
                record["smiles"] = chemistry.get("CanonicalSMILES")
        records.append(record)

    print(f"\n{len(records)} ligands: "
          + ", ".join(f"{r['name']} ({r['role']})" for r in records))
    with_site = [r for r in records if r["site_residues"]]
    print(f"  {len(with_site)} carry a residue-level site, "
          f"all evidence level "
          f"{sorted({r['site_evidence'] for r in with_site}) or 'n/a'}")
    print(f"  {sum(1 for r in records if r.get('potency'))} carry a potency "
          f"with a citation")

    if problems:
        print("\nREFUSING TO WRITE:")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    if args.check:
        print("\nvalidation passed (not written)")
        return 0

    destination = RESOURCE_DIR / "ligands.json"
    destination.write_text(json.dumps({
        "schema": 1,
        "note": ("Binding sites here are INFERRED from mutagenesis, docking or "
                 "geometry. No PIEZO structure with a bound small-molecule "
                 "modulator has been deposited, and the build verifies that "
                 "against the downloaded entries rather than asserting it."),
        "site_evidence_levels": SITE_EVIDENCE,
        "ligands": records,
    }, indent=1))
    print(f"\nwrote {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
