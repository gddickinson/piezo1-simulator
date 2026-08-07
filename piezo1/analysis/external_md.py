"""What deposited molecular dynamics can actually contribute here.

Round 42 proposed comparing this project's geometrically-found lipid sites
against occupancies from simulations other people have already run — on the
principle that an independent method agreeing is worth more than a better
version of ours. That principle is right. The data to act on it is not there.

**What was measured, with a working control.** Probing MemProtMD for every one
of the 21 catalogued PIEZO entries returns exactly **one** hit, 3JAC, the 2015
structure. The other twenty — including 7WLT, 7WLU, 6B3R and 8YEZ, which are the
structures this project actually uses — are absent. Canonical MemProtMD entries
(2RH1, 1M0L) return 200 on the same probe, so the absence is about PIEZO rather
than about the request.

**And the one entry cannot answer the question.** 3JAC resolves 918 of 2,547
residues (36%). Of the four curated lipid-associated residue groups it resolves
the PIP2 cluster in full and **none** of the three blade basic clusters — 11 of
the 15 curated lipid-binding residues are simply not in the model that was
simulated.

The other two named sources do not help. Zenodo carries PIEZO1 datasets, but the
ones that exist are microscopy TIFFs and PDFs rather than trajectories; GPCRmd
is GPCR-specific and PIEZO1 is not a GPCR.

So this module does not implement a comparison. It implements the *check*, so
that the conclusion is reproducible and will change by itself when the situation
does — if MemProtMD ever ingests a modern PIEZO structure, the coverage function
says so and the test guarding this null fails.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from dataclasses import dataclass, field

__all__ = ["MD_SOURCES", "MemProtMDCoverage", "memprotmd_coverage",
           "lipid_site_coverage", "assess"]

MEMPROTMD = "https://memprotmd.bioch.ox.ac.uk/_ref/PDB/{pdb}/"
_USER_AGENT = "piezo1-simulator/0.1 (research use)"

#: Entries known to be in MemProtMD, used to prove the probe works. Without a
#: control, "everything is absent" is indistinguishable from "the request is
#: wrong", and this project has already been caught by that class of error.
CONTROLS = ("2rh1", "1m0l")

#: What each source the roadmap named actually provides, measured rather than
#: assumed.
MD_SOURCES = {
    "MemProtMD": {
        "provides": "coarse-grained and atomistic membrane simulations with "
                    "lipid-contact analysis, keyed by PDB entry",
        "piezo_coverage": "1 of 21 catalogued entries (3JAC only)",
        "usable": False,
        "why": "the single entry resolves 36% of the sequence and none of the "
               "three blade basic clusters; and the API is not publicly "
               "documented, so the analysis is browsable but not fetchable",
    },
    "Zenodo": {
        "provides": "author-deposited datasets",
        "piezo_coverage": "PIEZO1 records exist but are microscopy images and "
                          "PDFs, not trajectories",
        "usable": False,
        "why": "no deposited PIEZO1 trajectory or lipid-occupancy table found",
    },
    "GPCRmd": {
        "provides": "deposited GPCR simulations",
        "piezo_coverage": "none",
        "usable": False,
        "why": "PIEZO1 is not a GPCR; the resource is out of scope",
    },
}


@dataclass
class MemProtMDCoverage:
    """Which catalogued structures MemProtMD holds."""

    present: list = field(default_factory=list)
    absent: list = field(default_factory=list)
    controls_ok: bool = False
    checked: bool = False
    note: str = ""

    @property
    def fraction(self) -> float:
        total = len(self.present) + len(self.absent)
        return len(self.present) / total if total else float("nan")

    def summary(self) -> str:
        if not self.checked:
            return f"not checked: {self.note}"
        return (f"MemProtMD holds {len(self.present)} of "
                f"{len(self.present) + len(self.absent)} catalogued PIEZO "
                f"entries ({self.fraction:.0%}): {self.present or 'none'}. "
                f"Probe control {'passed' if self.controls_ok else 'FAILED'}.")


def _reachable(pdb: str, timeout: int = 40) -> bool | None:
    """``True``/``False`` if answered, ``None`` if the network did not respond.

    The distinction matters: a network failure must not be recorded as an
    absence, or an offline run would manufacture this round's conclusion.
    """
    request = urllib.request.Request(MEMPROTMD.format(pdb=pdb.lower()),
                                     headers={"User-Agent": _USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status == 200
    except urllib.error.HTTPError as exc:
        return False if exc.code == 404 else None
    except (urllib.error.URLError, TimeoutError, OSError):
        return None


def memprotmd_coverage(pdbs=None, offline: bool = False,
                       cache_name: str = "memprotmd_coverage"
                       ) -> MemProtMDCoverage:
    """Probe MemProtMD for each catalogued entry, with a control and a cache."""
    from ..config import CACHE_DIR

    path = CACHE_DIR / "external_md" / f"{cache_name}.json"
    if path.exists():
        stored = json.loads(path.read_text())
        return MemProtMDCoverage(present=stored["present"],
                                 absent=stored["absent"],
                                 controls_ok=stored["controls_ok"],
                                 checked=True, note="from cache")
    if offline:
        return MemProtMDCoverage(checked=False, note="offline and not cached")

    if pdbs is None:
        from ..io.registry import load_registry
        pdbs = [e.pdb for e in load_registry().entries]

    controls = [_reachable(p) for p in CONTROLS]
    if not all(c is True for c in controls):
        # Without a working control the probe proves nothing.
        return MemProtMDCoverage(checked=False,
                                 note="probe control failed; not recording "
                                      "absences that may be network failures")

    present, absent = [], []
    for pdb in pdbs:
        answer = _reachable(pdb)
        if answer is None:
            return MemProtMDCoverage(
                checked=False,
                note=f"network did not answer for {pdb}; refusing to record "
                     f"a partial result as coverage")
        (present if answer else absent).append(pdb)

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"present": present, "absent": absent,
                                "controls_ok": True}))
    return MemProtMDCoverage(present=present, absent=absent, controls_ok=True,
                             checked=True)


def lipid_site_coverage(pdb: str = "3JAC") -> dict:
    """Which curated lipid-associated residues a structure actually resolves.

    The question Round 42 wanted answered is about lipid contacts, so a
    simulation of a model that omits the lipid-binding residues cannot answer
    it however good the simulation is.
    """
    from ..core.annotations import load_annotations
    from ..core.sequence import human_to_mouse
    from ..core.structure import Structure
    from ..io.registry import load_registry
    from ..structure.protomers import modelled_residues

    record = load_registry().get(pdb)
    if record is None or not record.available:
        return {"error": f"{pdb} not downloaded"}
    structure = Structure.from_file(record.path)
    resolved = modelled_residues(structure)

    groups = {}
    for group in load_annotations("human").residue_groups:
        if not any(key in group.label.lower()
                   for key in ("pip2", "basic", "lipid")):
            continue
        mouse = [human_to_mouse(r) for r in group.residues]
        have = [m for m in mouse if m and m in resolved]
        groups[group.label] = {"resolved": len(have),
                               "total": len(group.residues),
                               "complete": len(have) == len(group.residues)}
    total = sum(g["total"] for g in groups.values())
    have = sum(g["resolved"] for g in groups.values())
    return {"pdb": pdb, "residues_resolved": len(resolved),
            "groups": groups, "lipid_residues_resolved": have,
            "lipid_residues_total": total,
            "can_address_lipid_contacts": have == total}


def assess(offline: bool = False) -> dict:
    """The whole answer: what exists, and whether it can support the comparison."""
    coverage = memprotmd_coverage(offline=offline)
    sites = lipid_site_coverage("3JAC") if coverage.present else {}
    return {
        "sources": MD_SOURCES,
        "memprotmd": {"summary": coverage.summary(),
                      "present": coverage.present, "checked": coverage.checked,
                      "fraction": coverage.fraction if coverage.checked else None},
        "only_entry_lipid_sites": sites,
        "comparison_possible": bool(
            coverage.checked and coverage.present
            and sites.get("can_address_lipid_contacts")),
        "note": ("Round 42's premise was that deposited MD makes this cheap. "
                 "Measured, it does not: 1 of 21 entries, and that one resolves "
                 "none of the three blade basic clusters."),
    }
