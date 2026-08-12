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
from piezo1.core.numbering_check import identify_numbering  # noqa: E402
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
#: AlphaFold models that belong in the chooser, keyed by file stem. They were
#: skipped entirely until now, which meant the only way to look at one was to
#: open the file by hand — and the graft that uses them could not be compared
#: against the thing it grafts.
#:
#: They are MONOMERS. Everything in this project that needs a three-fold axis —
#: the dome, the pore, the elastic network, the paralogue comparison — will
#: refuse them, and that is correct rather than a gap. AlphaFold DB serves one
#: model per accession, so there is exactly one conformation each: for the same
#: fold in different states, use a deposited entry with the Completeness
#: selector instead.
PREDICTED_MODELS = {
    "AF-Q92508-F1-MODEL_V6": dict(
        species="human", state="predicted", gating="unknown",
        note="AlphaFold DB model of human PIEZO1, the whole 2,521 residues. "
             "A MONOMER and a PREDICTION: the B-factor column holds pLDDT, "
             "and anything needing three protomers will refuse it. This is "
             "the model the full-length graft takes its distal blade from.",
        recommended_for=["predicted", "full_length_source"]),
    "AF-E2JF22-F1-MODEL_V6": dict(
        species="mouse", state="predicted", gating="unknown",
        note="AlphaFold DB model of mouse Piezo1, 2,547 residues. Monomer, "
             "prediction, pLDDT in the B-factor column.",
        recommended_for=["predicted", "full_length_source"]),
    "AF-Q0KL00-F1-MODEL_V6": dict(
        species="rat", state="predicted", gating="unknown",
        note="AlphaFold DB model of rat Piezo1, 2,535 residues. The third "
             "mammalian PIEZO1 and the only structural representation it has; "
             "much of the mechanosensitivity electrophysiology is rat.",
        recommended_for=["predicted", "family"]),
    "AF-A0A061ACU2-F1-MODEL_V6": dict(
        species="worm", state="predicted", gating="unknown",
        note="AlphaFold DB model of C. elegans PEZO-1 isoform g, 2,442 "
             "residues. Held beside the deposited 9UOY/9ZIS so the prediction "
             "can be scored against experiment for the same protein — the "
             "control the human and mouse models have never had.",
        recommended_for=["predicted", "family"]),
    "AF-M9MSG8-F1-MODEL_V6": dict(
        species="fly", state="predicted", gating="unknown",
        note="AlphaFold DB model of Drosophila PIEZO, 2,551 residues. In "
             "CANONICAL numbering, where the deposited 9W7X is not.",
        recommended_for=["predicted", "family"]),
    "AF-F4IN58-F1-MODEL_V6": dict(
        species="plant", state="predicted", gating="unknown",
        note="AlphaFold DB model of Arabidopsis PIEZO, 2,462 residues. The "
             "only structural representation of a non-animal PIEZO that "
             "exists — Dictyostelium pzoA has neither a structure nor a "
             "model — and therefore the only way to ask whether the dome is "
             "a property of the fold rather than of animals. A PREDICTION and "
             "a MONOMER: read the pLDDT before believing the blade.",
        recommended_for=["predicted", "family", "generality"]),
}

#: Downloaded coordinate files that do not become catalogue entries, each with
#: the reason. A refusal with a stated reason is a record; a name quietly
#: missing from a glob is a gap nobody can tell from an oversight.
EXCLUDED = {
    "6U32": "HaloTag bound to tetramethylrhodamine — not a PIEZO at all. It is "
            "downloaded because the fusion model needs the tag's own fold.",
    # These two were previously excluded with the wrong reason on file: "a
    # Piezo domain from a distant organism ... cataloguing it would need a
    # seventh reference". Both halves were wrong. The RCSB cross-references
    # them to A0A061ACU2, which this project has held all along — they are
    # C. elegans PEZO-1, the same protein as 9UOY — and adding the reference
    # does not help, because they still score 0.081 and 0.077 against every
    # one of the nine. The numbering is the construct's own, running 14-278
    # for a 291-residue expression fragment, and `canonical_renumbering` finds
    # no shift that repairs it, correctly: recovering canonical numbers from a
    # construct needs an alignment, not an offset.
    "4PKE": "C. elegans PEZO-1 beta-sandwich domain, 211 modelled residues in "
            "the construct's own numbering (14-278, matching no reference "
            "above 0.081). A monomer, so the dome, the pore and the elastic "
            "network all refuse it; and unreadable by residue number, so no "
            "annotation can be applied to it either.",
    "4PKX": "The second crystal form of the same domain, 235 residues, same "
            "construct numbering (0.077).",
}

NOT_A_PIEZO = set(EXCLUDED)

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
    # --- the homologues, added when a search for PIEZO entries found six the
    # --- catalogue was missing. Two are human PIEZO2, which makes the
    # --- paralogue comparison a same-species one for the first time; three are
    # --- invertebrate, which asks the generality question across half a
    # --- billion years rather than across the mammals.
    "9VEE": dict(species="human", state="curved", gating="closed",
                 note="Human PIEZO2 with MDFIC2. The paralogue in the SAME "
                      "species as our PIEZO1 reference, so a PIEZO1-PIEZO2 "
                      "comparison no longer has to cross species too.",
                 recommended_for=["piezo2", "paralogue", "human_piezo2"]),
    "9VEF": dict(species="human", state="curved", gating="closed",
                 note="Human PIEZO2 with MDFIC, the partner complex to 9VEE.",
                 recommended_for=["piezo2", "paralogue"]),
    "9UOY": dict(species="worm", state="curved", gating="closed",
                 note="C. elegans PEZO-1. An invertebrate PIEZO, so neither a "
                      "PIEZO1 nor a PIEZO2 — the duplication that made those "
                      "two is vertebrate. The most distant structure available.",
                 recommended_for=["invertebrate", "generality"]),
    "9ZIS": dict(species="worm", state="curved", gating="closed",
                 note="C. elegans PEZO-1 isoform G, the full-length 2,442 "
                      "residue product — an independent 3.5 A dataset of the "
                      "same isoform as 9UOY. The catalogue's only replicate "
                      "pair, which is what lets an isoform difference be told "
                      "apart from a dataset difference.",
                 recommended_for=["invertebrate", "generality", "isoform",
                                  "replicate"]),
    "9ZIT": dict(species="worm", state="curved", gating="closed",
                 note="C. elegans PEZO-1 isoform K, which begins at residue "
                      "757 and so lacks a third of the blade. Deposited in "
                      "CANONICAL numbering (identity 1.000 over 801-2442), "
                      "not the isoform's own — unlike 6LQI and 9W7X, which "
                      "are not, and unlike 4PKE/4PKX, which are in a "
                      "construct's. Three conventions in one catalogue is why "
                      "the numbering is measured on every load.",
                 recommended_for=["invertebrate", "isoform"]),
    "9UOX": dict(species="worm", state="curved", gating="closed",
                 note="C. elegans PEZO-1 isoform K at 3.8 A, the replicate of "
                      "9ZIT. Also canonical numbering, 808-2437.",
                 recommended_for=["invertebrate", "isoform", "replicate"]),
    "9W7X": dict(species="fly", state="curved", gating="closed",
                 note="Drosophila PIEZO. Deposited in an isoform's own "
                      "numbering, +3 after residue 1570 — found by the "
                      "numbering check, not by reading the paper.",
                 recommended_for=["invertebrate", "generality"]),

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
        if pdb in NOT_A_PIEZO:
            continue
        if pdb.startswith("AF-") and pdb not in PREDICTED_MODELS:
            continue
        info = meta.get(pdb)
        if pdb in PREDICTED_MODELS:
            # AlphaFold models have no RCSB entry, and asking for one 404s and
            # drops them. Their provenance comes from the AlphaFold DB fetch,
            # which already records the version it discovered.
            info = {"pdb": pdb, "res": None, "date": "",
                    "title": "AlphaFold DB predicted model",
                    "journal": "Nucleic Acids Res", "year": 2024,
                    "pmid": None, "doi": "10.1093/nar/gkad1011", "emdb": []}
        elif info is None or args.refresh:
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

        # Which PIEZO this is, **measured** rather than curated: the file's own
        # residue names are scored against every reference sequence, and the
        # right one wins at 1.000 with the runner-up below 0.25. Curating it
        # would have been a second place for the answer to be wrong, and this
        # is the field the structure chooser filters on.
        identity = identify_numbering(st)
        protein = identity.protein
        if not identity.explained:
            print(f"  ! {pdb}: numbering not identified — {identity.summary()}")

        cur = PREDICTED_MODELS.get(pdb) or CURATION.get(
            pdb, dict(species="unknown", state="unclassified",
                      gating="unknown", note="", recommended_for=[]))
        entries.append({
            "pdb": pdb, "file": path.name,
            "resolution": info.get("res"), "released": info.get("date"),
            "title": info.get("title"), "journal": info.get("journal"),
            "year": info.get("year"), "pmid": info.get("pmid"),
            "doi": info.get("doi"), "emdb": info.get("emdb", []),
            "n_atoms": st.n_atoms, "n_protomers": len(chains),
            "protomer_chains": chains, "ligands": ligands,
            "protein": protein, "numbering": identity.reference,
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
