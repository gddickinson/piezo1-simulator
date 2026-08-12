"""Cached downloaders for every external resource the project uses.

Nothing downloaded is committed, so this module is what makes a fresh clone
reproducible. Every function is idempotent: it skips a file that already
exists unless ``force=True``.

Run the whole set with::

    python -m piezo1.io.fetch
"""

from __future__ import annotations

import json
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from ..config import (DERIVED_DIR, DICTY_PIEZO_ACC, FLY_PIEZO_ACC, HUMAN_ACC,
                      HUMAN_PIEZO2_ACC, LIGAND_DIR, MOUSE_ACC,
                      MOUSE_PIEZO2_ACC, PLANT_PIEZO_ACC, RAT_ACC, SEQUENCE_DIR,
                      STRUCTURE_DIR, WORM_PIEZO_ACC, ensure_dirs)

__all__ = ["fetch_pdb", "fetch_alphafold", "fetch_uniprot", "fetch_ligand",
           "fetch_all", "fetch_chap_grid", "fetch_cds", "CDS_TRANSCRIPTS",
           "DEFAULT_PDB_IDS", "DEFAULT_LIGANDS", "FAMILY_ACCESSIONS",
           "ALPHAFOLD_ACCESSIONS", "ALPHAFOLD_UNAVAILABLE",
           "CHAP_GRID_URL", "CHAP_LICENCE", "CONTENT_CHECKS"]

#: Every reviewed PIEZO, keyed by the short name its resource file carries.
#:
#: Nine, and that is the whole family — ``reviewed:true AND protein_name:piezo``
#: against UniProt returns exactly this set. Recorded as a list rather than
#: re-queried at runtime because a pinned accession is provenance and a live
#: query is not: the same call next year may return ten, and the difference is
#: something to notice deliberately rather than absorb silently.
FAMILY_ACCESSIONS = {
    "human": HUMAN_ACC, "mouse": MOUSE_ACC, "rat": RAT_ACC,
    "human_piezo2": HUMAN_PIEZO2_ACC, "mouse_piezo2": MOUSE_PIEZO2_ACC,
    "worm_piezo": WORM_PIEZO_ACC, "fly_piezo": FLY_PIEZO_ACC,
    "plant_piezo": PLANT_PIEZO_ACC, "dicty_piezo": DICTY_PIEZO_ACC,
}

#: The family members AlphaFold DB holds a model of the **canonical** sequence
#: for. Five of nine, and the two absences are facts about the database rather
#: than failures of this fetch, so they are named in
#: :data:`ALPHAFOLD_UNAVAILABLE` and not requested.
ALPHAFOLD_ACCESSIONS = (HUMAN_ACC, MOUSE_ACC, RAT_ACC,
                        WORM_PIEZO_ACC, FLY_PIEZO_ACC, PLANT_PIEZO_ACC)

#: What AlphaFold DB does not have, and why. Recorded rather than left as a
#: gap in the list above, because "no model" and "nobody asked" look identical
#: once a name is simply missing.
ALPHAFOLD_UNAVAILABLE = {
    HUMAN_PIEZO2_ACC: ("no model of the canonical 2,752-residue sequence; the "
                       "database holds isoform 2 (2,689 aa) and isoform 3 "
                       "(709 aa) only"),
    DICTY_PIEZO_ACC: ("no model at any length; at 3,080 residues pzoA is past "
                      "the ceiling of the whole-proteome predictions"),
    MOUSE_PIEZO2_ACC: "not requested; the mouse PIEZO2 structure 6KG7 resolves "
                      "all 38 transmembrane helices, so nothing needs one",
}

#: The water free-energy landscape underlying the Rao et al. 2019 hydrophobic
#: gating heuristic, as published in the CHAP repository. 100x100 grid over
#: (normalised Wimley-White hydrophobicity, pore radius in nm) -> kJ/mol,
#: derived from ~600 MD simulations of ~200 channel structures.
CHAP_GRID_URL = ("https://raw.githubusercontent.com/channotation/chap/master/"
                 "scripts/heuristic/heuristic_grid.json")
#: CHAP is MIT licensed, which is why this project can use its published grid
#: directly rather than reconstructing the boundary by eye from a figure.
CHAP_LICENCE = "MIT (Klesse, Rao, Sansom & Tucker)"

USER_AGENT = "piezo1-simulator/0.1 (research use)"
TIMEOUT = 300

#: Structures the application ships a registry entry for.
#:
#: The list is the answer to a structured RCSB query — every entry whose
#: polymer cross-references one of the nine reviewed PIEZO accessions — rather
#: than a full-text search for "piezo", which returns 200 entries of which most
#: are unrelated crystals that merely mention the word.
DEFAULT_PDB_IDS = [
    # human PIEZO1, including three disease variants
    "8YEZ", "8ZU3", "8ZU8", "8YFC", "8YFG", "9VMX",
    # mouse curved/flattened pairs
    "7WLT", "7WLU", "11YE", "11ZC", "8IXN", "8IXO",
    # mouse reference and historical
    "6B3R", "6BPZ", "5Z10", "3JAC", "8IMZ", "6LQI", "4RAX", "9VED",
    # PIEZO2, for comparison and for the distal blade
    "6KG7",
    # C. elegans PEZO-1, deposited as two *isoforms*: g is the full-length
    # 2442-residue product and k starts at 757. Two entries each, which is
    # what makes them worth holding — an isoform pair with a replicate is the
    # only place in the catalogue where "the difference between these two
    # models" can be separated from "the difference between two datasets".
    "9UOY", "9ZIS",      # isoform g
    "9ZIT", "9UOX",      # isoform k
    # The beta-sandwich domain of PEZO-1 at 2.5 A, in two crystal forms. The
    # only atomic-resolution PIEZO coordinates that exist.
    "4PKE", "4PKX",
    # Drosophila PIEZO
    "9W7X",
    # human PIEZO2
    "9VEE", "9VEF",
]

#: PubChem compounds, by name and CID.
DEFAULT_LIGANDS = {
    "Yoda1": 2746822,
    "Yoda2": 170908031,
    "Jedi1": 736516,
    "Jedi2": 2796026,
    "Dooku1": 137321150,
    "Yaddle1": 171378918,
}


@dataclass
class FetchResult:
    path: Path
    downloaded: bool
    size: int
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.error is None


def _check_cif(data: bytes) -> str:
    head = data[:4096].decode("utf-8", "replace")
    if not head.lstrip().startswith("data_"):
        return "not an mmCIF: no data_ block header"
    if b"_atom_site" not in data:
        return "mmCIF has no _atom_site loop, so it carries no coordinates"
    return ""


def _check_fasta(data: bytes) -> str:
    if not data.lstrip().startswith(b">"):
        return "not FASTA: no > header"
    body = b"".join(data.split(b"\n")[1:])
    if not body.strip():
        return "FASTA header with no sequence"
    return ""


def _check_json(data: bytes) -> str:
    try:
        json.loads(data)
    except ValueError as exc:
        return f"not JSON: {exc}"
    return ""


def _check_sdf(data: bytes) -> str:
    if b"$$$$" not in data:
        return "not SDF: no $$$$ record terminator"
    if b"V2000" not in data and b"V3000" not in data:
        return "SDF has no counts line, so it carries no atoms"
    return ""


#: What each kind of download has to look like once it arrives.
#:
#: The size guard is necessary and not sufficient, which the project has learned
#: twice: Round 60 found an Ensembl endpoint returning an HTML error, and Round
#: 65 found two 127-byte error pages **stored as structures** — the file existed,
#: had a plausible name, and every later step treated it as data. An error page
#: from a CDN is comfortably larger than 200 bytes, so only content can catch it.
CONTENT_CHECKS = {
    "cif": _check_cif,
    "fasta": _check_fasta,
    "json": _check_json,
    "sdf": _check_sdf,
}


def _download(url: str, dest: Path, force: bool = False,
              headers: dict | None = None, kind: str | None = None) -> FetchResult:
    """Fetch ``url`` to ``dest``, refusing to write anything that is not ``kind``.

    Nothing is written until the payload has been checked, so a rejected
    download leaves no file behind for a later step to pick up and believe.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force and dest.stat().st_size > 0:
        return FetchResult(dest, False, dest.stat().st_size)
    request_headers = {"User-Agent": USER_AGENT}
    request_headers.update(headers or {})
    req = urllib.request.Request(url, headers=request_headers)
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            data = r.read()
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
        return FetchResult(dest, False, 0, f"{type(exc).__name__}: {exc}")
    if len(data) < 200:
        return FetchResult(dest, False, len(data),
                           f"suspiciously small response ({len(data)} bytes)")
    check = CONTENT_CHECKS.get(kind or "")
    if check is not None:
        problem = check(data)
        if problem:
            return FetchResult(dest, False, len(data),
                               f"{problem} ({len(data)} bytes from {url})")
    dest.write_bytes(data)
    return FetchResult(dest, True, len(data))


# --------------------------------------------------------------------------
# Individual fetchers
# --------------------------------------------------------------------------

def fetch_pdb(pdb_id: str, force: bool = False) -> FetchResult:
    """Download an mmCIF coordinate file from the RCSB."""
    pdb_id = pdb_id.upper()
    return _download(f"https://files.rcsb.org/download/{pdb_id}.cif",
                     STRUCTURE_DIR / f"{pdb_id}.cif", force, kind="cif")


def fetch_alphafold(accession: str, force: bool = False,
                    with_pae: bool = False) -> list[FetchResult]:
    """Download the AlphaFold DB model **for the canonical sequence**.

    The version number is discovered from the API rather than guessed. This is
    not pedantry: the v4 URLs that most documentation still shows now return
    404, and the current model for PIEZO1 is v6.

    Neither is the isoform check. The endpoint returns one entry *per isoform*
    and the canonical is not first — it is not always present at all — and this
    took ``entries[0]``:

    * **Q9H5I5, human PIEZO2**: two entries, isoform 3 (**709 aa**) and isoform
      2 (2,689 aa). AlphaFold DB has no model for the canonical 2,752-residue
      sequence. ``entries[0]`` is the 709-residue one, a quarter of the
      protein, which arrives named ``AF-Q9H5I5-3-F1`` and parses perfectly.
    * **A0A061ACU2, PEZO-1**: **twelve** entries, one per annotated isoform,
      from 1,038 to 2,442 residues. The canonical happens to be first, so this
      was right by luck rather than by construction.

    A model of a different isoform read in canonical numbering is wrong
    everywhere past its first splice difference, and nothing downstream could
    tell: it is a well-formed mmCIF of the right protein. So the entry is
    selected by exact accession — an isoform carries a ``-N`` suffix — and if
    there is no canonical model this **refuses and says what it found**,
    because substituting the nearest isoform is the failure mode, not the
    fallback.
    """
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    dest_hint = STRUCTURE_DIR / f"AF-{accession}.cif"
    try:
        req = urllib.request.Request(api, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            entries = json.loads(r.read())
    except Exception as exc:
        return [FetchResult(dest_hint, False, 0, f"AlphaFold API failed: {exc}")]
    if not entries:
        return [FetchResult(dest_hint, False, 0, "no AlphaFold entry")]

    canonical = [e for e in entries
                 if e.get("uniprotAccession", "") == accession]
    if not canonical:
        offered = ", ".join(
            f"{e.get('uniprotAccession')} ({len(e.get('uniprotSequence', ''))} aa)"
            for e in entries)
        return [FetchResult(
            dest_hint, False, 0,
            f"AlphaFold DB has no model for canonical {accession}; it offers "
            f"only isoforms ({offered}). Refused rather than substituted — an "
            f"isoform model read in canonical numbering is wrong past its "
            f"first splice difference and nothing downstream could tell.")]

    entry = canonical[0]
    out = []
    cif_url = entry["cifUrl"]
    out.append(_download(cif_url, STRUCTURE_DIR / Path(cif_url).name, force,
                         kind="cif"))
    if with_pae and entry.get("paeDocUrl"):
        url = entry["paeDocUrl"]
        out.append(_download(url, STRUCTURE_DIR / Path(url).name, force,
                             kind="json"))
    return out


def fetch_uniprot(accession: str, species: str, force: bool = False) -> list[FetchResult]:
    """Download the UniProt FASTA and full JSON entry."""
    out = [
        _download(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
                  SEQUENCE_DIR / f"{accession}_{species}_PIEZO1.fasta", force,
                  kind="fasta"),
        _download(f"https://rest.uniprot.org/uniprotkb/{accession}.json",
                  SEQUENCE_DIR / f"{accession}_{species}_PIEZO1.json", force,
                  kind="json"),
    ]
    return out


#: Compounds PubChem has no computed 3D conformer for. Their 2D records are
#: still fetched; the 3D request is skipped so it does not read as a failure.
NO_3D_CONFORMER = {170908031}   # Yoda2


def fetch_ligand(name: str, cid: int, force: bool = False) -> list[FetchResult]:
    """Download 2D and, where PubChem has one, 3D SDF records."""
    base = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
    out = [_download(f"{base}/SDF?record_type=2d",
                     LIGAND_DIR / f"{name}_CID{cid}_2d.sdf", force, kind="sdf")]
    if cid not in NO_3D_CONFORMER:
        out.append(_download(f"{base}/SDF?record_type=3d",
                             LIGAND_DIR / f"{name}_CID{cid}_3d.sdf", force,
                             kind="sdf"))
    return out


# --------------------------------------------------------------------------
# Bulk
# --------------------------------------------------------------------------

def fetch_chap_grid(force: bool = False) -> FetchResult:
    """Download the CHAP hydrophobic-gating free-energy grid.

    Kept out of the repository like every other download; regenerable with
    ``python -m piezo1.io.fetch``. Analyses degrade to "unavailable" without
    it rather than failing.
    """
    return _download(CHAP_GRID_URL, DERIVED_DIR / "chap_heuristic_grid.json",
                     force, kind="json")


#: Canonical coding sequences, from the Ensembl transcript UniProt cross-links
#: to. Recorded here rather than looked up each time so the exact transcript is
#: part of the provenance: PIEZO1 has splice isoforms, and 6LQI is one of them.
CDS_TRANSCRIPTS = {
    "human": ("ENST00000301015", "Q92508"),
    "mouse": ("ENSMUST00000156333", "E2JF22"),
}


def fetch_cds(species: str = "human", force: bool = False) -> FetchResult:
    """Download the canonical coding sequence for a species.

    Real DNA, not a back-translation. Back-translating a protein would produce
    a sequence that looks like a gene and is not one — the codon choices would
    be invented, and silent variants could not be represented at all.
    """
    transcript, _accession = CDS_TRANSCRIPTS[species]
    # The content type must be an HTTP header. Ensembl used to honour it as a
    # `;content-type=` query parameter and no longer does — the plain URL now
    # answers 415, and this fetch had been failing silently on any machine
    # without a warm cache. Found by Round 60's empty-clone run, which is the
    # only way a broken download surfaces once the file is already on disk.
    return _download(
        f"https://rest.ensembl.org/sequence/id/{transcript}?type=cds",
        SEQUENCE_DIR / f"{transcript}_{species}_PIEZO1_cds.fasta", force,
        headers={"Content-Type": "text/x-fasta"}, kind="fasta")


def fetch_all(force: bool = False, structures: bool = True,
              sequences: bool = True, ligands: bool = True,
              alphafold: bool = True, pae: bool = False,
              verbose: bool = True) -> list[FetchResult]:
    """Download everything the application needs."""
    ensure_dirs()
    results: list[FetchResult] = []

    def note(label: str, res: FetchResult) -> None:
        if verbose:
            if not res.ok:
                print(f"  FAIL {label}: {res.error}")
            elif res.downloaded:
                print(f"  got  {label}  ({res.size / 1e6:.1f} MB)")
            else:
                print(f"  have {label}")

    if sequences:
        print("UniProt sequences")
        for sp, acc in FAMILY_ACCESSIONS.items():
            for r in fetch_uniprot(acc, sp, force):
                note(f"{acc} {r.path.suffix}", r)
                results.append(r)

    if structures:
        print(f"RCSB structures ({len(DEFAULT_PDB_IDS)})")
        for pdb in DEFAULT_PDB_IDS:
            r = fetch_pdb(pdb, force)
            note(pdb, r)
            results.append(r)

    if alphafold:
        print(f"AlphaFold DB models ({len(ALPHAFOLD_ACCESSIONS)})")
        for acc in ALPHAFOLD_ACCESSIONS:
            # PAE only for the two the full-length graft is built from; the
            # matrices are ~50 MB each and nothing reads the others.
            want_pae = pae and acc in (HUMAN_ACC, MOUSE_ACC)
            for r in fetch_alphafold(acc, force, with_pae=want_pae):
                note(f"AF {acc}", r)
                results.append(r)

    if ligands:
        print(f"PubChem ligands ({len(DEFAULT_LIGANDS)})")
        for name, cid in DEFAULT_LIGANDS.items():
            for r in fetch_ligand(name, cid, force):
                note(f"{name} {r.path.stem[-2:]}", r)
                results.append(r)

    if sequences:
        print("Coding sequences (Ensembl)")
        for species in CDS_TRANSCRIPTS:
            r = fetch_cds(species, force)
            note(f"{species} CDS", r)
            results.append(r)

    print("CHAP hydrophobic-gating grid")
    r = fetch_chap_grid(force)
    note("chap heuristic grid", r)
    results.append(r)

    ok = sum(1 for r in results if r.ok)
    new = sum(1 for r in results if r.downloaded)
    failed = [r for r in results if not r.ok]
    print(f"\n{ok}/{len(results)} resources present ({new} newly downloaded)")
    if failed:
        print(f"{len(failed)} failed:")
        for r in failed:
            print(f"  {r.path.name}: {r.error}")
    return results


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Download all external data.")
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    ap.add_argument("--pae", action="store_true",
                    help="also fetch AlphaFold PAE matrices (large)")
    ap.add_argument("--only", choices=["structures", "sequences", "ligands",
                                       "alphafold"], default=None)
    args = ap.parse_args()
    kw = dict(structures=True, sequences=True, ligands=True, alphafold=True)
    if args.only:
        kw = {k: (k == args.only) for k in kw}
    results = fetch_all(force=args.force, pae=args.pae, **kw)
    return 0 if all(r.ok for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
