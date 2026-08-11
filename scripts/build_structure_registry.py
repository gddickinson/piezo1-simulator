#!/usr/bin/env python
"""Author the structure registry resource.

Writes ``piezo1/resources/structures.json``: the catalogue the application
offers in its structure chooser. Metadata comes from the RCSB entry API;
conformational state and "recommended for" come from the structural literature
review in ``ref/research/02_structures.md`` and are hand-declared here.

Usage::

    python scripts/build_structure_registry.py            # use cached metadata
    python scripts/build_structure_registry.py --refresh  # re-query RCSB
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import CACHE_DIR, RESOURCE_DIR, STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402

# Hand-declared scientific interpretation, keyed by PDB ID. Anything not listed
# is still included, but with state "unclassified".
#: Downloaded coordinate files that are **not** PIEZO structures and must not
#: join a catalogue documented as one. The build sweeps ``data/structures`` for
#: every ``.cif``, which was harmless until Round 31 downloaded the HaloTag
#: crystal structure for the fusion geometry; the next rebuild — Round 83's —
#: swept it in as an unclassified "unknown" species, and it took two unrelated
#: tests failing to notice. `fusion.load_halotag` reads 6U32 by path and never
#: through the registry, so nothing needs it here.
NOT_A_PIEZO = {
    "6U32",   # HaloTag bound to its tetramethylrhodamine ligand
}

CURATION = {
    "8YEZ": dict(species="human", state="curved", gating="closed",
                 note="Highest-resolution apo human PIEZO1.",
                 recommended_for=["default", "human_reference", "variant_mapping"]),
    "8ZU3": dict(species="human", state="curved", gating="closed",
                 note="Human PIEZO1 in complex with its auxiliary subunit MDFIC.",
                 recommended_for=["mdfic_complex"]),
    "8ZU8": dict(species="human", state="curved", gating="closed",
                 note="Gain-of-function variant A1988V. Note A1988 itself is "
                      "not modelled.",
                 recommended_for=["variant_structure"]),
    "8YFC": dict(species="human", state="curved", gating="closed",
                 note="A1988V with MDFIC.", recommended_for=["variant_structure"]),
    "8YFG": dict(species="human", state="curved", gating="closed",
                 note="R2456H with MDFIC. The only disease variant whose "
                      "mutated residue is actually resolved in its own structure.",
                 recommended_for=["variant_structure"]),
    "9VMX": dict(species="human", state="curved", gating="closed",
                 note="E756del (the African malaria-associated allele) with "
                      "MDFIC. E756 is not modelled.",
                 recommended_for=["variant_structure"]),
    "7WLT": dict(species="mouse", state="curved", gating="closed",
                 note="Curved mPIEZO1 in a lipid bilayer. Curved endpoint of "
                      "the classic flattening pair.",
                 recommended_for=["morph_start", "dome_reference"]),
    "7WLU": dict(species="mouse", state="flattened", gating="closed",
                 note="Flattened mPIEZO1 in a lipid bilayer. Low resolution "
                      "(6.81 A) and the pore remains closed.",
                 recommended_for=["morph_end"]),
    "11YE": dict(species="mouse", state="curved", gating="closed",
                 note="Curved mPiezo1 in native plasma-membrane vesicles.",
                 recommended_for=["morph_start", "dome_reference"]),
    "11ZC": dict(species="mouse", state="flat", gating="open-like",
                 note="Flat mPiezo1 in native plasma-membrane vesicles. The "
                      "only flat endpoint with pore changes consistent with "
                      "conduction. 6.0 A and largely backbone-only.",
                 recommended_for=["morph_end", "open_state"]),
    "8IXN": dict(species="mouse", state="curved", gating="closed",
                 note="S2472E phosphomimetic, curved.",
                 recommended_for=["morph_start"]),
    "8IXO": dict(species="mouse", state="intermediate", gating="intermediate",
                 note="S2472E intermediate. Best-matched resolution pair with "
                      "8IXN and the only pair crossing a gating boundary.",
                 recommended_for=["morph_end", "intermediate_state"]),
    "8IMZ": dict(species="mouse", state="curved", gating="closed",
                 note="mouse Piezo1-MDFIC complex.", recommended_for=["mdfic_complex"]),
    "6B3R": dict(species="mouse", state="curved", gating="closed",
                 note="Classic 2017 structure; largest resolved residue count.",
                 recommended_for=["historical"]),
    "6BPZ": dict(species="mouse", state="curved", gating="closed",
                 note="Guo & MacKinnon; source of the dome model.",
                 recommended_for=["historical", "dome_reference"]),
    "5Z10": dict(species="mouse", state="curved", gating="closed",
                 note="Zhao et al. 2018; the lever-like transduction model.",
                 recommended_for=["historical"]),
    "3JAC": dict(species="mouse", state="curved", gating="closed",
                 note="The first Piezo1 cryo-EM structure (2015). Much of the "
                      "chain is poly-UNK with arbitrary numbering.",
                 recommended_for=["historical"]),
    "6LQI": dict(species="mouse", state="curved", gating="closed",
                 note="The Piezo1.1 splice isoform (delta 1382-1405).",
                 recommended_for=["isoform"]),
    "4RAX": dict(species="mouse", state="fragment", gating="n/a",
                 note="X-ray structure of the cap/CED at 1.45 A - by far the "
                      "highest resolution available for any part of PIEZO1.",
                 recommended_for=["cap_detail"]),
    # The note said "resolves residues 8-823" until Round 83 measured it. It
    # resolves 8-2822 in 16 segments, 1817 C-alphas per protomer - more than
    # any PIEZO1 entry here - and it is in MOUSE Piezo2 numbering (Q8CD54,
    # 2822 aa), not human PIEZO2's 2752. Both facts are checked by
    # tests/test_paralogue.py against the file itself.
    "6KG7": dict(species="mouse", state="curved", gating="closed",
                 note="PIEZO2 (mouse Piezo2, Q8CD54 numbering), the paralogue "
                      "control. Resolves 1,817 residues from 8 to 2,822 in 16 "
                      "segments - more than any PIEZO1 entry - including all "
                      "38 transmembrane helices, so it is both the best "
                      "experimental view of the distal blade and the only "
                      "structure that can separate PIEZO1 from the fold.",
                 recommended_for=["piezo2", "distal_blade", "paralogue"]),
    "9VED": dict(species="mouse", state="curved", gating="closed",
                 note="mouse Piezo1-MDFI complex (2026).",
                 recommended_for=["mdfic_complex"]),
}


def rcsb_metadata(pdb: str) -> dict:
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb}"
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()

    cache = CACHE_DIR / "rcsb_piezo_entries.json"
    meta: dict[str, dict] = {}
    if cache.exists() and not args.refresh:
        for row in json.loads(cache.read_text()):
            meta[row["pdb"]] = row

    entries = []
    for path in sorted(STRUCTURE_DIR.glob("*.cif")):
        pdb = path.stem.upper()
        if pdb.startswith("AF-") or pdb in NOT_A_PIEZO:
            continue
        info = meta.get(pdb)
        if info is None or args.refresh:
            try:
                e = rcsb_metadata(pdb)
            except Exception as exc:
                print(f"  ! {pdb}: {exc}")
                continue
            cite = e.get("rcsb_primary_citation", {})
            res = e.get("rcsb_entry_info", {}).get("resolution_combined") or [None]
            info = {
                "pdb": pdb, "res": res[0],
                "date": e.get("rcsb_accession_info", {}).get("initial_release_date", "")[:10],
                "title": e.get("struct", {}).get("title", ""),
                "journal": cite.get("rcsb_journal_abbrev"), "year": cite.get("year"),
                "pmid": cite.get("pdbx_database_id_pub_med"),
                "doi": cite.get("pdbx_database_id_doi"),
                "emdb": e.get("rcsb_entry_container_identifiers", {}).get("emdb_ids", []),
            }

        st = Structure.from_file(path)
        chains = []
        for ch in st.chains:
            m = st.mask_ca() & (st.chain == ch)
            if m.sum() > 300:
                chains.append({"chain": ch, "n_ca": int(m.sum()),
                               "first": int(st.res_seq[m].min()),
                               "last": int(st.res_seq[m].max())})
        ligands = sorted(set(st.res_name[st.mask_ligands()].tolist()))

        cur = CURATION.get(pdb, dict(species="unknown", state="unclassified",
                                     gating="unknown", note="", recommended_for=[]))
        entries.append({
            "pdb": pdb, "file": path.name,
            "resolution": info.get("res"), "released": info.get("date"),
            "title": info.get("title"), "journal": info.get("journal"),
            "year": info.get("year"), "pmid": info.get("pmid"),
            "doi": info.get("doi"), "emdb": info.get("emdb", []),
            "n_atoms": st.n_atoms, "n_protomers": len(chains),
            "protomer_chains": chains, "ligands": ligands,
            **cur,
        })

    entries.sort(key=lambda e: (e["species"] != "human", e.get("resolution") or 99))
    dest = RESOURCE_DIR / "structures.json"
    dest.write_text(json.dumps({"entries": entries}, indent=1))

    print(f"wrote {dest} with {len(entries)} entries\n")
    print(f"{'PDB':6s} {'species':7s} {'state':12s} {'res':>5s} {'prot':>4s} "
          f"{'range':>12s}  ligands")
    for e in entries:
        rng = (f"{e['protomer_chains'][0]['first']}-{e['protomer_chains'][0]['last']}"
               if e["protomer_chains"] else "-")
        print(f"{e['pdb']:6s} {e['species']:7s} {e['state']:12s} "
              f"{str(e['resolution']):>5s} {e['n_protomers']:>4d} {rng:>12s}  "
              f"{','.join(e['ligands'][:4])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
