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
sys.path.insert(0, str(Path(__file__).resolve().parent))

from piezo1.config import PAPER_DIR, RESOURCE_DIR  # noqa: E402
from reference_seed import SEED  # noqa: E402

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
    "rawicz2000": "Membrane mechanics", "perozo2002": "Membrane mechanics",
    "kyte1982": "Methods", "vonheijne1992": "Methods",
    "pauling1951": "Methods",
    "dolinsky2004": "Methods",
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
            # `MANUAL` was referenced here and never defined, so the *rejection*
            # path — the one that matters — raised NameError instead of
            # reporting the rejection. Latent until Round 84d, when the title
            # gate correctly refused a wrong PMID for Liu et al. 2025
            # (39674176 resolves to a paper about IgG in adipose tissue; the
            # right one is 39719701) and the script died instead of saying so.
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
