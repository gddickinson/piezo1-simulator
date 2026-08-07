#!/usr/bin/env python
"""Build the project bibliography and download the open-access papers.

Two outputs:

* ``piezo1/resources/references.json`` — committed, machine-readable. Every
  citation the code and docs rely on, with the metadata resolved from Europe
  PMC rather than typed by hand, so author lists and years cannot drift.
* ``docs/REFERENCES.md`` — generated from that JSON, human-readable.

With ``--download`` it also fetches the open-access full texts into
``ref/papers/`` (git-ignored). Only genuinely open-access items are fetched:
the script asks Europe PMC whether a full text is available rather than
guessing at publisher URLs, and records what it could not get.

Usage::

    python scripts/build_references.py
    python scripts/build_references.py --download
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import PAPER_DIR, RESOURCE_DIR  # noqa: E402

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
UA = "piezo1-simulator/0.1 (research use; bibliography builder)"

#: The bibliography seed. ``key`` is the citation key used in code and docs;
#: ``topic`` says what the project uses it for; ``expect`` is a word that MUST
#: appear in the resolved title.
#:
#: That last field exists because of a real failure. Several PMIDs entered from
#: memory resolved cleanly to entirely unrelated papers — one PIEZO1 structure
#: citation came back as a bone-marrow transplantation study, another as a
#: stem-cell reprogramming paper. Europe PMC returns whatever that identifier
#: actually points at; it has no way to know you meant something else. A
#: bibliography with confidently wrong citations is worse than no bibliography,
#: so every entry is now resolved by title and then checked against ``expect``,
#: and anything failing the check is reported rather than written out.
SEED = [
    # --- structure -------------------------------------------------------
    ("ge2015", "PMID:26390154", "First Piezo1 cryo-EM structure (3JAC)", "piezo1"),
    ("kamajaya2014", "PMID:25242456", "First Piezo structures: C. elegans CED (4PKE/4PKX)", "piezo"),
    ("guo2017", "PMID:29231809", "The dome model; PDB 6BPZ", "piezo"),
    ("saotome2018", "PMID:29261642", "Piezo1 structure 6B3R", "piezo1"),
    ("zhao2018", "PMID:29469092", "Lever-like transduction; PDB 5Z10", "piezo1"),
    ("wang2019piezo2", "PMID:31435011", "Piezo2 structure 6KG7", "piezo2"),
    ("geng2020", "PMID:32142647", "Plug-and-latch gating; Piezo1.1 isoform 6LQI", "piezo"),
    ("yang2022", "PMID:35388220", "Curved and flattened mPIEZO1 in bilayer; 7WLT/7WLU", "piezo1"),
    ("zhou2023mdfic", "PMID:37590348", "MDFIC is a PIEZO auxiliary subunit", "piezo"),
    ("vaisey2026", "PMID:42234740", "Lipid cofactor required; force alone insufficient", "piezo1"),
    # --- membrane mechanics ----------------------------------------------
    ("haselwandter2018", "PMID:30480546", "Membrane footprint; R_c 10.2 nm, lambda 14 nm", "piezo"),
    ("haselwandter2022", "PMID:36166476", "Elastic properties and shape of the Piezo dome", "piezo"),
    ("haselwandter2022b", "PMID:36166475", "Quantitative prediction and measurement of the membrane footprint", "piezo"),
    ("dixit2025", "DOI:10.7554/eLife.105138.3", "Nanodome excess area and elasticity", "piezo"),
    ("chong2021", "PMID:33582137", "Full-length model; dome depth 6-7 nm", "piezo1"),
    ("devecchis2021", "PMID:33582135", "MD of opening by membrane tension", "piezo1"),
    # --- electrophysiology and kinetics ----------------------------------
    ("coste2015", "PMID:26008989", "Pore properties dictated by the C-terminal region", "piezo1"),
    ("bae2013", "PMID:23487776", "Xerocytosis mutations slow inactivation", "piezo1"),
    ("lewis2015", "PMID:26646186", "Tension sensitivity; T50 2.7 and 4.7 mN/m", "piezo1"),
    ("cox2016", "PMID:26785635", "Bilayer tension gating; dG0 9.7 kT, dA 8 nm2", "piezo1"),
    ("young2023", "PMID:36795747", "Four-state TENSION model - the one implemented", "mechanotransduction"),
    ("lewis2021", "PMID:34711306", "Clustering does not alter gating", "piezo1"),
    # --- lipids and pharmacology -----------------------------------------
    ("borbiro2015", "PMID:25670203", "Phosphoinositide dependence", "piezo"),
    ("ridone2020", "PMID:32582958", "Cholesterol; P50 shift on depletion", "piezo1"),
    ("shi2020", "PMID:33027663", "Sphingomyelinase disables inactivation", "piezo1"),
    ("romero2019", "PMID:30867417", "Dietary fatty acids tune the mechanical response", "piezo1"),
    ("romero2020", "PMID:32561714", "Margaric acid and PIEZO2", "mechanical"),
    ("buyan2020", "PMID:32949489", "PIP2 and cholesterol sites; K2166-K2169", "piezo1"),
    ("buyan2023", "PMID:35927961", "Lipid redistribution in the curved footprint", "piezo1"),
    ("hashad2025", "PMID:41433068", "PIP2 corrects an endothelial channelopathy", "piezo1"),
    ("syeda2015", "PMID:26001275", "Yoda1", "mechanotransduction"),
    ("botellosmith2019", "PMID:31582801", "Yoda1 mechanism and binding pocket", "piezo1"),
    ("wang2018jedi", "PMID:29610524", "Jedi1/2 and the lever transduction pathway", "piezo1"),
    ("evans2018dooku", "PMID:29498036", "Dooku1 antagonises Yoda1", "yoda1"),
    ("bae2011", "PMID:21696149", "GsMTx4 inhibits Piezo1", "piezo1"),
    ("poole2014", "PMID:24662763", "STOML3 tunes the displacement threshold", "piezo"),
    ("qi2015", "PMID:26443885", "STOML3 membrane stiffening", "stoml3"),
    ("wetzel2017", "PMID:27941788", "STOML3 inhibitors reverse mechanical hypersensitivity", "stoml3"),
    # --- genetics and disease --------------------------------------------
    ("zarychanski2012", "PMID:22529292", "PIEZO1 mutations cause hereditary xerocytosis", "piezo1"),
    ("albuisson2013", "PMID:23695678", "DHS1 mutations; M2225R, R2456H", "piezo1"),
    ("andolfo2013", "PMID:23479567", "DHS1; multiple clinical forms", "piezo1"),
    ("fotiou2015", "PMID:26333996", "PIEZO1 loss of function causes lymphatic dysplasia", "piezo1"),
    ("ma2018malaria", "PMID:29576450", "E756del and malaria resistance", "piezo1"),
    ("karamaticcrew2023", "PMID:36723926", "PIEZO1 carries the Er blood group antigens", "er"),
    # --- methods ----------------------------------------------------------
    ("atilgan2001", "PMID:11159421", "Anisotropic network model", "elastic network"),
    ("bahar2010", "PMID:19785456", "Normal mode analysis of membrane proteins - review", "normal mode"),
    ("smart1996hole", "PMID:9195488", "HOLE pore-radius algorithm", "pore"),
    ("rao2019heuristic", "PMID:31235590",
     "Hydrophobic-gating heuristic: the (hydrophobicity, radius) landscape",
     "hydrophobic"),
    ("klesse2019chap", "PMID:31220459",
     "CHAP - pore annotation; source of the MIT-licensed grid", "chap"),
    ("beckstein2003", "PMID:12740433",
     "Liquid-vapour oscillations of water in hydrophobic nanopores", "water"),
    ("aryal2015", "PMID:25106689", "Hydrophobic gating in ion channels - review",
     "hydrophobic"),
    ("labesse1997", "PMID:9183534", "P-SEA: secondary structure from CA geometry", "secondary structure"),
    ("leguilloux2009", "PMID:19486540",
     "fpocket: alpha-sphere pocket detection - source of the 3.0-5.5 A radii",
     "pocket"),
    ("shrake1973", "PMID:4760134",
     "Shrake-Rupley numerical SASA; source of the 1.4 A water probe", "solvent"),
    ("kabsch1976", "DOI:10.1107/S0567739476001873", "Optimal rotation superposition", "rotation"),
    ("jumper2021", "PMID:34265844", "AlphaFold", "protein structure prediction"),
    ("varadi2024", "PMID:37933859", "AlphaFold DB", "alphafold"),
    # --- external variant predictors, reached through ProtVar (CC BY 4.0) ---
    # Free-text queries rather than remembered PMIDs: six citations in Round 8
    # were resolved from memory to entirely unrelated papers, so every entry
    # here is looked up by title and gated on `expect`.
    ("stephenson2024protvar",
     'TITLE:"ProtVar: mapping and contextualizing human missense variation"',
     "ProtVar API - serves the predictors below under CC BY 4.0", "protvar"),
    ("cheng2023alphamissense",
     'TITLE:"Accurate proteome-wide missense variant effect prediction with AlphaMissense"',
     "AlphaMissense pathogenicity", "alphamissense"),
    ("frazer2021eve",
     'TITLE:"Disease variant prediction with deep generative models of evolutionary data"',
     "EVE - unsupervised variant effect from evolutionary data", "variant"),
    ("brandes2023esm1b",
     'TITLE:"Genome-wide prediction of disease variant effects with a deep protein language model"',
     "ESM-1b variant effects", "language model"),
    ("schymkowitz2005foldx",
     'TITLE:"The FoldX web server: an online force field"',
     "FoldX force field - source of the precomputed ddG", "foldx"),
    # --- HaloTag labelling, imported with the kinetics in Round 32 ----------
    # Looked up by title rather than by remembered PMID, for the reason stated
    # above the ProtVar block.
    ("los2008halotag",
     'TITLE:"HaloTag: a novel protein labeling technology for cell imaging and protein analysis"',
     "HaloTag chemistry; source of the covalent on-rate and its irreversibility",
     "halotag"),
    ("grimm2015jf",
     'TITLE:"A general method to improve fluorophores for live-cell and single-molecule microscopy"',
     "Janelia Fluor dyes including JF646 - cell-permeable, so partition ~ 1",
     "dye"),
    # --- ion permeation, Round 33 ------------------------------------------
    ("coste2010piezo",
     'TITLE:"Piezo1 and Piezo2 are essential components of distinct mechanically activated cation channels"',
     "The original PIEZO1 characterisation; single-channel conductance and "
     "cation non-selectivity", "piezo"),
    ("gnanasambandam2015",
     'TITLE:"Ionic Selectivity and Permeation Properties of Human PIEZO1 Channels"',
     "PIEZO1 selectivity and permeation - the direct target of the PNP model",
     "selectivity"),
    ("hall1975access",
     'TITLE:"Access resistance of a small circular pore"',
     "Access resistance of a circular pore mouth - the term that limits a short "
     "wide pore", "access resistance"),
    # --- population constraint, Round 41 ------------------------------------
    ("chen2024gnomad",
     'TITLE:"A genomic mutational constraint map using variation in 76,156 human genomes"',
     "gnomAD constraint: LOEUF, pLI and missense z-scores", "constraint"),
    # --- calcium nanodomain, Round 35 ---------------------------------------
    ("stern1992",
     'TITLE:"Buffering of calcium in the vicinity of a channel pore"',
     "The steady-state buffered-diffusion Green's function this model uses",
     "calcium"),
    ("naraghi1997",
     'TITLE:"Linearized buffered Ca2+ diffusion in microdomains and its implications for calculation of [Ca2+] at the mouth of a calcium channel"',
     "Linearised buffered diffusion; the screening length and its validity",
     "calcium"),
    ("allbritton1992",
     'TITLE:"Range of messenger action of calcium ion and inositol 1,4,5-trisphosphate"',
     "Cytosolic calcium diffusion coefficient and buffering range", "calcium"),
    ("tsien1980bapta",
     'TITLE:"New calcium indicators and buffers with high selectivity against magnesium and protons: design, synthesis, and properties of prototype structures"',
     "BAPTA - the chelator scaffold of the JF646-BAPTA sensor, and its Kd",
     "calcium"),
    ("bertaccini2025piezo1",
     'TITLE:"Visualizing PIEZO1 Localization and Activity in hiPSC-Derived Single Cells and Organoids with HaloTag Technology"',
     "The tagged-PIEZO1 experiment this labelling model describes; three tags "
     "per channel and the multi-level brightness histogram",
     "halotag"),
]


#: Entries Europe PMC cannot resolve, typed by hand. Kept deliberately small —
#: anything here is a citation whose metadata is not machine-verified.
MANUAL = {
    "kabsch1976": {
        "title": "A solution for the best rotation to relate two sets of vectors",
        "authors": "Kabsch W.",
        "journal": "Acta Crystallographica Section A",
        "year": "1976", "volume": "32", "pages": "922-923",
        "pmid": None, "pmcid": None,
        "doi": "10.1107/S0567739476001873",
        "is_open_access": False, "has_pdf": False,
        "url": "https://doi.org/10.1107/S0567739476001873",
        "note": "Predates PubMed indexing; metadata entered by hand.",
    },
}


def _get(url: str, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _to_epmc_query(seed: str) -> str:
    """Translate a seed identifier into Europe PMC search syntax.

    Europe PMC does **not** understand ``PMID:12345`` — it silently returns
    zero hits rather than erroring, which is the worst possible failure mode
    for a bibliography builder. The PubMed identifier lives in ``EXT_ID``, and
    the source must be pinned to MED or the same ID can match a preprint.
    """
    if seed.startswith("PMID:"):
        return f'EXT_ID:{seed[5:]} AND SRC:MED'
    if seed.startswith("PMCID:"):
        return f'PMCID:{seed[6:]}'
    if seed.startswith("DOI:"):
        return f'DOI:"{seed[4:]}"'
    return seed


def resolve(seed: str) -> dict | None:
    """Resolve one seed identifier to full Europe PMC metadata."""
    for query in (_to_epmc_query(seed), seed.split(":", 1)[-1]):
        url = (f"{EPMC}/search?query={urllib.parse.quote(query)}"
               f"&format=json&resultType=core&pageSize=1")
        try:
            data = json.loads(_get(url))
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            print(f"    ! network error: {exc}")
            return None
        hits = data.get("resultList", {}).get("result", [])
        if hits:
            return hits[0]
    return None


def distil(key: str, topic: str, rec: dict) -> dict:
    ft = rec.get("fullTextIdList", {}).get("fullTextId", [])
    return {
        "key": key,
        "topic": topic,
        "title": rec.get("title", "").rstrip("."),
        "authors": rec.get("authorString", ""),
        "journal": (rec.get("journalInfo", {}) or {}).get("journal", {}).get("title")
                   or rec.get("journalTitle", ""),
        "year": rec.get("pubYear"),
        "volume": (rec.get("journalInfo", {}) or {}).get("volume"),
        "pages": rec.get("pageInfo"),
        "pmid": rec.get("pmid"),
        "pmcid": rec.get("pmcid"),
        "doi": rec.get("doi"),
        "is_open_access": rec.get("isOpenAccess") == "Y",
        "has_pdf": "pdf" in ft,
        "url": (f"https://doi.org/{rec['doi']}" if rec.get("doi") else
                f"https://pubmed.ncbi.nlm.nih.gov/{rec.get('pmid')}/"),
    }


def download_paper(entry: dict) -> tuple[str, str]:
    """Fetch the open-access full text. Returns ``(status, detail)``."""
    pmcid = entry.get("pmcid")
    if not pmcid:
        return "skipped", "no PMC record"
    if not entry.get("is_open_access"):
        return "skipped", "not open access"

    PAPER_DIR.mkdir(parents=True, exist_ok=True)
    stem = f"{entry['key']}_{pmcid}"

    # Prefer the PDF; fall back to the open-access XML full text.
    targets = [
        (f"{EPMC}/{pmcid}/fullTextPDF", PAPER_DIR / f"{stem}.pdf", 20000),
        (f"{EPMC}/{pmcid}/fullTextXML", PAPER_DIR / f"{stem}.xml", 2000),
    ]
    for url, dest, min_size in targets:
        if dest.exists() and dest.stat().st_size > min_size:
            return "have", dest.name
        try:
            blob = _get(url, timeout=180)
        except Exception:
            continue
        if len(blob) >= min_size:
            dest.write_bytes(blob)
            return "downloaded", f"{dest.name} ({len(blob) / 1e6:.1f} MB)"
    return "failed", "no retrievable full text"


def write_markdown(entries: list[dict], dest: Path) -> None:
    groups: dict[str, list[dict]] = {}
    for e in entries:
        groups.setdefault(e.get("section", "Other"), []).append(e)

    lines = [
        "# References",
        "",
        "Every source this project relies on, with what it is used *for*.",
        "Generated by `scripts/build_references.py` from Europe PMC metadata —",
        "do not edit by hand; edit the seed list in that script instead.",
        "",
        "Open-access full texts are downloaded to `ref/papers/` with",
        "`python scripts/build_references.py --download`. That directory is",
        "git-ignored: the papers are other people's copyright, and the",
        "bibliography here is enough to retrieve them.",
        "",
        f"**{len(entries)} references.** "
        f"{sum(1 for e in entries if e['is_open_access'])} are open access.",
        "",
    ]
    for section in ("Structure", "Membrane mechanics",
                    "Electrophysiology and kinetics", "Lipids and pharmacology",
                    "Genetics and disease", "Methods", "Other"):
        rows = groups.get(section)
        if not rows:
            continue
        lines += [f"## {section}", ""]
        for e in sorted(rows, key=lambda x: (x.get("year") or "", x["key"])):
            oa = " · **OA**" if e["is_open_access"] else ""
            ident = []
            if e.get("pmid"):
                ident.append(f"[PMID {e['pmid']}]"
                             f"(https://pubmed.ncbi.nlm.nih.gov/{e['pmid']}/)")
            if e.get("doi"):
                ident.append(f"[doi:{e['doi']}](https://doi.org/{e['doi']})")
            lines.append(
                f"- **`{e['key']}`** — {e['authors']} "
                f"*{e['title']}.* {e['journal']} {e.get('year') or ''}"
                f"{(';' + e['volume']) if e.get('volume') else ''}"
                f"{(':' + e['pages']) if e.get('pages') else ''}. "
                f"{' · '.join(ident)}{oa}  \n"
                f"  <sub>Used for: {e['topic']}</sub>")
        lines.append("")
    dest.write_text("\n".join(lines) + "\n")


SECTION_OF = {
    "ge2015": "Structure", "kamajaya2014": "Structure", "guo2017": "Structure",
    "saotome2018": "Structure", "zhao2018": "Structure",
    "wang2019piezo2": "Structure", "geng2020": "Structure",
    "yang2022": "Structure", "zhou2023mdfic": "Structure",
    "vaisey2026": "Structure",
    "haselwandter2018": "Membrane mechanics",
    "haselwandter2022": "Membrane mechanics", "haselwandter2022b": "Membrane mechanics", "dixit2025": "Membrane mechanics",
    "chong2021": "Membrane mechanics", "devecchis2021": "Membrane mechanics",
    "coste2015": "Electrophysiology and kinetics",
    "bae2013": "Electrophysiology and kinetics",
    "lewis2015": "Electrophysiology and kinetics",
    "cox2016": "Electrophysiology and kinetics",
    "lewis2017": "Electrophysiology and kinetics",
    "young2023": "Electrophysiology and kinetics",
    "lewis2021": "Electrophysiology and kinetics",
    "borbiro2015": "Lipids and pharmacology", "ridone2020": "Lipids and pharmacology",
    "shi2020": "Lipids and pharmacology", "romero2019": "Lipids and pharmacology",
    "romero2020": "Lipids and pharmacology", "buyan2020": "Lipids and pharmacology",
    "buyan2023": "Lipids and pharmacology", "hashad2025": "Lipids and pharmacology",
    "syeda2015": "Lipids and pharmacology",
    "botellosmith2019": "Lipids and pharmacology",
    "wang2018jedi": "Lipids and pharmacology",
    "evans2018dooku": "Lipids and pharmacology", "bae2011": "Lipids and pharmacology",
    "poole2014": "Lipids and pharmacology", "qi2015": "Lipids and pharmacology",
    "wetzel2017": "Lipids and pharmacology",
    "zarychanski2012": "Genetics and disease",
    "albuisson2013": "Genetics and disease", "andolfo2013": "Genetics and disease",
    "fotiou2015": "Genetics and disease", "ma2018malaria": "Genetics and disease",
    "karamaticcrew2023": "Genetics and disease",
    "atilgan2001": "Methods", "yang2009": "Methods", "bahar2010": "Methods",
    "smart1996hole": "Methods", "rao2019heuristic": "Methods",
    "leguilloux2009": "Methods", "shrake1973": "Methods",
    "klesse2019chap": "Methods", "beckstein2003": "Methods",
    "aryal2015": "Methods", "labesse1997": "Methods", "kabsch1976": "Methods",
    "jumper2021": "Methods", "varadi2024": "Methods",
    "stephenson2024protvar": "Methods", "cheng2023alphamissense": "Methods",
    "frazer2021eve": "Methods", "brandes2023esm1b": "Methods",
    "schymkowitz2005foldx": "Methods",
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--download", action="store_true",
                    help="also fetch open-access full texts into ref/papers/")
    ap.add_argument("--refresh", action="store_true",
                    help="re-query Europe PMC even for cached entries")
    args = ap.parse_args()

    dest = RESOURCE_DIR / "references.json"
    cached = {}
    if dest.exists() and not args.refresh:
        cached = {e["key"]: e for e in json.loads(dest.read_text())["references"]}

    entries, unresolved = [], []
    for key, query, topic, expect in SEED:
        if key in cached:
            e = cached[key]
            e["topic"] = topic
            e["section"] = SECTION_OF.get(key, "Other")
            entries.append(e)
            continue
        print(f"  resolving {key} ({query})")
        rec = resolve(query)
        if rec is not None and expect:
            haystack = (rec.get("title", "") + " " +
                        rec.get("abstractText", "") or "").lower()
            if expect.lower() not in haystack:
                print(f"    ! REJECTED: resolved title does not mention "
                      f"{expect!r} -> {rec.get('title', '')[:70]!r}")
                rec = None
        if rec is None:
            if key in MANUAL:
                e = dict(MANUAL[key], key=key, topic=topic,
                         section=SECTION_OF.get(key, "Other"))
                entries.append(e)
                print(f"    (manual entry)")
                continue
            print(f"    ! could not resolve {key}")
            unresolved.append((key, query, topic))
            continue
        e = distil(key, topic, rec)
        e["section"] = SECTION_OF.get(key, "Other")
        entries.append(e)
        time.sleep(0.2)

    entries.sort(key=lambda e: e["key"])
    dest.write_text(json.dumps(
        {"note": ("Bibliography for the PIEZO1 simulator. Metadata resolved "
                  "from Europe PMC, not hand-typed. Regenerate with "
                  "scripts/build_references.py."),
         "n_references": len(entries),
         "unresolved": [{"key": k, "query": q, "topic": t}
                        for k, q, t in unresolved],
         "references": entries}, indent=1))

    write_markdown(entries, Path("docs/REFERENCES.md"))
    oa = sum(1 for e in entries if e["is_open_access"])
    print(f"\n{len(entries)}/{len(SEED)} references resolved, {oa} open access")
    for k, q, *_ in unresolved:
        print(f"  UNRESOLVED  {k}  ({q})")

    if args.download:
        print(f"\ndownloading open-access full texts into {PAPER_DIR}")
        tally: dict[str, int] = {}
        for e in entries:
            status, detail = download_paper(e)
            tally[status] = tally.get(status, 0) + 1
            if status in ("downloaded", "failed"):
                print(f"  {status:11s} {e['key']:20s} {detail}")
        print("  " + ", ".join(f"{v} {k}" for k, v in sorted(tally.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
