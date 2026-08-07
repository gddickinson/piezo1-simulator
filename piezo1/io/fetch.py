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

from ..config import (DERIVED_DIR, HUMAN_ACC, HUMAN_PIEZO2_ACC, LIGAND_DIR,
                      MOUSE_ACC, SEQUENCE_DIR, STRUCTURE_DIR, ensure_dirs)

__all__ = ["fetch_pdb", "fetch_alphafold", "fetch_uniprot", "fetch_ligand",
           "fetch_all", "fetch_chap_grid", "fetch_cds", "CDS_TRANSCRIPTS",
           "DEFAULT_PDB_IDS", "DEFAULT_LIGANDS",
           "CHAP_GRID_URL", "CHAP_LICENCE"]

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
DEFAULT_PDB_IDS = [
    # human PIEZO1, including three disease variants
    "8YEZ", "8ZU3", "8ZU8", "8YFC", "8YFG", "9VMX",
    # mouse curved/flattened pairs
    "7WLT", "7WLU", "11YE", "11ZC", "8IXN", "8IXO",
    # mouse reference and historical
    "6B3R", "6BPZ", "5Z10", "3JAC", "8IMZ", "6LQI", "4RAX", "9VED",
    # PIEZO2, for comparison and for the distal blade
    "6KG7",
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


def _download(url: str, dest: Path, force: bool = False,
              headers: dict | None = None) -> FetchResult:
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
    dest.write_bytes(data)
    return FetchResult(dest, True, len(data))


# --------------------------------------------------------------------------
# Individual fetchers
# --------------------------------------------------------------------------

def fetch_pdb(pdb_id: str, force: bool = False) -> FetchResult:
    """Download an mmCIF coordinate file from the RCSB."""
    pdb_id = pdb_id.upper()
    return _download(f"https://files.rcsb.org/download/{pdb_id}.cif",
                     STRUCTURE_DIR / f"{pdb_id}.cif", force)


def fetch_alphafold(accession: str, force: bool = False,
                    with_pae: bool = False) -> list[FetchResult]:
    """Download the current AlphaFold DB model for a UniProt accession.

    The version number is discovered from the API rather than guessed. This is
    not pedantry: the v4 URLs that most documentation still shows now return
    404, and the current model for PIEZO1 is v6.
    """
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    try:
        req = urllib.request.Request(api, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            entries = json.loads(r.read())
    except Exception as exc:
        return [FetchResult(STRUCTURE_DIR / f"AF-{accession}.cif", False, 0,
                            f"AlphaFold API failed: {exc}")]
    if not entries:
        return [FetchResult(STRUCTURE_DIR / f"AF-{accession}.cif", False, 0,
                            "no AlphaFold entry")]

    entry = entries[0]
    out = []
    cif_url = entry["cifUrl"]
    out.append(_download(cif_url, STRUCTURE_DIR / Path(cif_url).name, force))
    if with_pae and entry.get("paeDocUrl"):
        url = entry["paeDocUrl"]
        out.append(_download(url, STRUCTURE_DIR / Path(url).name, force))
    return out


def fetch_uniprot(accession: str, species: str, force: bool = False) -> list[FetchResult]:
    """Download the UniProt FASTA and full JSON entry."""
    out = [
        _download(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
                  SEQUENCE_DIR / f"{accession}_{species}_PIEZO1.fasta", force),
        _download(f"https://rest.uniprot.org/uniprotkb/{accession}.json",
                  SEQUENCE_DIR / f"{accession}_{species}_PIEZO1.json", force),
    ]
    return out


#: Compounds PubChem has no computed 3D conformer for. Their 2D records are
#: still fetched; the 3D request is skipped so it does not read as a failure.
NO_3D_CONFORMER = {170908031}   # Yoda2


def fetch_ligand(name: str, cid: int, force: bool = False) -> list[FetchResult]:
    """Download 2D and, where PubChem has one, 3D SDF records."""
    base = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
    out = [_download(f"{base}/SDF?record_type=2d",
                     LIGAND_DIR / f"{name}_CID{cid}_2d.sdf", force)]
    if cid not in NO_3D_CONFORMER:
        out.append(_download(f"{base}/SDF?record_type=3d",
                             LIGAND_DIR / f"{name}_CID{cid}_3d.sdf", force))
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
                     force)


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
        headers={"Content-Type": "text/x-fasta"})


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
        for acc, sp in ((HUMAN_ACC, "human"), (MOUSE_ACC, "mouse"),
                        (HUMAN_PIEZO2_ACC, "human_piezo2")):
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
        print("AlphaFold DB models")
        for acc in (HUMAN_ACC, MOUSE_ACC):
            for r in fetch_alphafold(acc, force, with_pae=pae):
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
