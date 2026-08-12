"""Reproducible analysis reports with full provenance.

A number without its provenance is not a result. Every report here records the
software version, the exact input file and where it came from, every parameter
used, and the versions of the libraries that did the arithmetic — so that a
figure in a talk can be traced back to the run that produced it, months later,
by someone who was not there.

The report is emitted as JSON (machine-readable, diffable) and Markdown
(readable). Both come from the same object, so they cannot disagree.
"""

from __future__ import annotations

import json
import platform
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ..parameters import PARAMETERS as _P

from .. import __version__
from ..core.structure import Structure
from ..io.registry import load_registry
from ..parameters import PARAMETERS
from ..structure.protomers import protomer_blocks

__all__ = ["Provenance", "AnalysisReport", "collect_provenance",
           "build_report", "ANALYSES"]


def _library_versions() -> dict[str, str]:
    out: dict[str, str] = {}
    for name in ("numpy", "scipy", "Bio", "moderngl", "PyQt6"):
        try:
            module = __import__(name)
            out[name] = str(getattr(module, "__version__", "unknown"))
        except Exception:
            out[name] = "not installed"
    return out


@dataclass
class Provenance:
    """Everything needed to reproduce, or to distrust, a number."""

    software_version: str
    timestamp: str
    python: str
    platform: str
    libraries: dict[str, str] = field(default_factory=dict)
    structure: dict = field(default_factory=dict)
    parameters: dict = field(default_factory=dict)
    #: Registry parameters differing from their documented defaults. Empty is
    #: the normal case and means every number below came from the published
    #: parameter set; non-empty means it did not, and the report says so at the
    #: top rather than leaving the reader to assume.
    parameter_overrides: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "software_version": self.software_version,
            "timestamp": self.timestamp, "python": self.python,
            "platform": self.platform, "libraries": self.libraries,
            "structure": self.structure, "parameters": self.parameters,
            "parameter_overrides": self.parameter_overrides,
            "warnings": self.warnings,
        }


def collect_provenance(structure: Structure | None = None,
                       parameters: dict | None = None) -> Provenance:
    prov = Provenance(
        software_version=__version__,
        timestamp=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        python=sys.version.split()[0],
        platform=f"{platform.system()} {platform.machine()}",
        libraries=_library_versions(),
        parameters=dict(parameters or {}),
        parameter_overrides=PARAMETERS.overrides(),
    )
    if prov.parameter_overrides:
        prov.warnings.append(
            "Computed with non-default parameters — "
            + PARAMETERS.override_summary()
            + ". These numbers are NOT comparable with the values in "
              "docs/SCIENCE.md, which were produced at the documented defaults.")
    if structure is not None:
        record = load_registry().get(structure.name)
        prov.structure = {
            "name": structure.name,
            "source_file": structure.source,
            "n_atoms": int(structure.n_atoms),
            "n_residues": int(structure.n_residues),
            "chains": structure.chains,
        }
        if record is not None:
            prov.structure.update({
                "species": record.species, "state": record.state,
                "gating": record.gating, "resolution": record.resolution,
                "citation": record.citation(), "doi": record.doi,
                "numbering": record.numbering_species,
            })
    return prov


@dataclass
class AnalysisReport:
    """A collection of named results sharing one provenance record."""

    provenance: Provenance
    results: dict = field(default_factory=dict)
    title: str = "PIEZO1 analysis report"

    def add(self, name: str, value: dict) -> None:
        self.results[name] = value

    def as_dict(self) -> dict:
        return {"title": self.title, "provenance": self.provenance.as_dict(),
                "results": self.results}

    def to_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.as_dict(), indent=1, default=_encode))
        return path

    def to_markdown(self, path: str | Path | None = None) -> str:
        p = self.provenance
        lines = [f"# {self.title}", ""]
        if p.parameter_overrides:
            # At the top rather than buried in the provenance footer: a reader
            # who skims must not take these numbers for the documented ones.
            lines += ["> **⚠ Non-default parameters.** "
                      + "; ".join(f"`{k}` = {v:g} (default "
                                  f"{PARAMETERS.default(k):g})"
                                  for k, v in sorted(p.parameter_overrides.items()))
                      + ". These results are not comparable with the values in "
                        "`docs/SCIENCE.md`.", ""]
        if p.structure:
            s = p.structure
            lines += [f"**{s.get('name')}** — {s.get('species', '?')}, "
                      f"{s.get('state', '?')}, "
                      f"{s.get('resolution', '?')} Å, "
                      f"{s.get('n_atoms', 0):,} atoms",
                      f"*{s.get('citation', '')}*", ""]
        for name, value in self.results.items():
            lines += [f"## {name.replace('_', ' ').capitalize()}", ""]
            lines += _render(value)
            lines.append("")
        lines += ["---", "", "## Provenance", "",
                  f"- piezo1 {p.software_version}, Python {p.python}, {p.platform}",
                  f"- run {p.timestamp}",
                  "- libraries: " + ", ".join(f"{k} {v}" for k, v in p.libraries.items())]
        if p.structure.get("source_file"):
            lines.append(f"- input: `{p.structure['source_file']}`")
        if p.parameters:
            lines.append("- parameters: " +
                         ", ".join(f"{k}={v}" for k, v in p.parameters.items()))
        lines.append("- parameter registry: " +
                     ("all at documented defaults"
                      if not p.parameter_overrides
                      else PARAMETERS.override_summary()))
        for w in p.warnings:
            lines.append(f"- ⚠ {w}")
        text = "\n".join(lines) + "\n"
        if path is not None:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        return text


def _encode(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)


def _render(value, indent: int = 0) -> list[str]:
    pad = "  " * indent
    if isinstance(value, dict):
        out = []
        for k, v in value.items():
            if isinstance(v, (dict, list)):
                out.append(f"{pad}- **{k}**:")
                out += _render(v, indent + 1)
            else:
                out.append(f"{pad}- **{k}**: {_fmt(v)}")
        return out
    if isinstance(value, list):
        return [f"{pad}- {_fmt(v)}" for v in value[:25]] + \
               ([f"{pad}- … {len(value) - 25} more"] if len(value) > 25 else [])
    return [f"{pad}{_fmt(value)}"]


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


# --------------------------------------------------------------------------
# The analyses
# --------------------------------------------------------------------------

# `_protomer_blocks` lived here until Round 78 and was a duplicate of
# `structure.protomers.protomer_blocks` — identical output on every trimer
# tested, differing only in its sentinel for a non-trimer and in hardcoding the
# 300-C-alpha floor instead of using `well_resolved_chains`. Removing it also
# broke a real import cycle: `report_tags` imported it from here while this
# module imports `report_tags` at the bottom, so `import report_tags` failed in
# a fresh interpreter and only ever worked because something imported `report`
# first.


def analysis_dome(st: Structure, species: str, **kw) -> dict:
    from ..structure.geometry import measure_dome, tm_surface_points
    blocks, _ = protomer_blocks(st)
    if not blocks:
        return {"error": "needs three well-resolved protomers"}
    pts, resolved = tm_surface_points(st, species)
    dome = measure_dome(blocks, pts)
    return {"n_transmembrane_resolved": len(resolved),
            "radius_of_curvature_nm": dome.radius_of_curvature / 10,
            "dome_depth_nm": dome.dome_depth / 10,
            "dome_area_nm2": dome.dome_area / 100,
            "projected_area_nm2": dome.projected_area / 100,
            "excess_area_nm2": dome.excess_area / 100,
            "c3_angle_deg": dome.notes["c3_angle_deg"],
            "c3_rmsd_A": dome.notes["c3_rmsd"],
            "reference": "closed-state R_c = 10.2 nm "
                         "(Haselwandter & MacKinnon 2018)"}


def analysis_pore(st: Structure, species: str, step: float | None = None,
                  **kw) -> dict:
    if step is None:
        step = _P.value("pore.step")
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    blocks, _ = protomer_blocks(st)
    if not blocks:
        return {"error": "needs three well-resolved protomers"}
    prof = pore_profile(st, detect_c3_axis(blocks), step=step)
    return {"bottleneck_radius_A": prof.bottleneck_radius,
            "bottleneck_z_A": prof.bottleneck_z,
            "bottleneck_lining": list(prof.bottleneck_lining()),
            "conductive": prof.is_conductive(),
            "constrictions": [{"z": s.z, "radius": s.radius,
                               "lining": list(s.lining)}
                              for s in prof.constrictions()],
            "note": prof.meta["note"]}


def analysis_hydration(st: Structure, species: str, step: float | None = None,
                       **kw) -> dict:
    """Hydrophobic-gating prediction (Rao et al. 2019) for this structure."""
    if step is None:
        step = _P.value("pore.step")
    from .hydration import load_grid, predict_wetting
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    blocks, _ = protomer_blocks(st)
    if not blocks:
        return {"error": "needs three well-resolved protomers"}
    grid = load_grid()
    if not grid.available:
        return {"error": "CHAP grid not downloaded; run python -m piezo1.io.fetch"}
    prof = pore_profile(st, detect_c3_axis(blocks), step=step)
    pred = predict_wetting(st, prof, grid)
    # Which constriction the verdict is about. "Sterically occluded" reads as a
    # shut gate, and the narrowest point is not at the gate in any entry in the
    # catalogue — see `analysis.pore_regions`.
    from .pore_regions import Bottleneck, describe_bottleneck
    try:
        where = describe_bottleneck(st, prof)
    except Exception as exc:                   # never lose the verdict itself
        where = Bottleneck(reason=f"bottleneck not located: {exc}")
    return {"score": pred.score,
            "cutoff": pred.meta["cutoff"],
            "bottleneck_radius_nm": pred.min_radius / 10.0,
            "bottleneck_region": where.narrowest_region,
            "bottleneck_location": where.sentence(),
            "transmembrane_gate_radius_nm": (
                where.gate.radius / 10.0 if where.gate else None),
            "hydrophobic_gate": pred.hydrophobic_gate,
            "sterically_occluded": pred.sterically_occluded,
            "conductive": pred.conductive,
            "verdict": pred.verdict,
            "dewetted": [{"residue": p.residue, "name": p.name,
                          "radius_nm": p.radius / 10.0,
                          "hydrophobicity": p.hydrophobicity,
                          "energy_kJ_per_mol": p.energy,
                          "distance": p.distance}
                         for p in pred.dewetted[:12]],
            "citation": pred.meta["citation"]}


def analysis_modes(st: Structure, species: str, n_modes: int = 20, **kw) -> dict:
    from ..physics.anm import ANM
    blocks, _ = protomer_blocks(st)
    if not blocks:
        return {"error": "needs three well-resolved protomers"}
    anm = ANM.from_trimer(blocks, cutoff=15.0).build()
    modes = anm.calc_modes(n_modes=n_modes)
    anm.label_symmetry(modes)
    a_modes = [i + 1 for i, s in enumerate(modes.symmetry) if s == "A"]
    return {"n_modes": modes.n_modes, "n_sites": modes.n_sites,
            "n_components": modes.meta["n_components"],
            "symmetric_modes": a_modes,
            "lowest_symmetric_mode": a_modes[0] if a_modes else None,
            "eigenvalues": [float(v) for v in modes.eigenvalues[:10]],
            "collectivity": [modes.collectivity(i) for i in range(min(5, modes.n_modes))],
            "note": "only A (three-fold symmetric) modes can couple to "
                    "isotropic membrane tension at first order"}


def analysis_pockets(st: Structure, species: str, top: int = 5, **kw) -> dict:
    from .pockets import find_pockets, ligand_contact_residues
    pockets = find_pockets(st)
    return {"n_pockets": len(pockets),
            "top": [{"index": p.index, "volume_A3": p.volume,
                     "buriedness": p.buriedness,
                     "n_lining_residues": len(p.residues),
                     "residues": list(p.residues[:20])}
                    for p in pockets[:top]],
            "ligand_contacts": ligand_contact_residues(st)}


def analysis_liu2025(st: Structure, species: str, **kw) -> dict:
    """Liu et al. 2025 panel by panel: what reproduces and what cannot.

    Runs across their four deposited states rather than against the loaded
    structure — the paper's claims are comparisons between states, and a single
    entry can only ever supply one side of one.
    """
    from .liu2025 import PANELS, coverage, replicate_all

    summary = coverage()
    return {
        "paper": summary["paper"],
        "coverage": {k: summary[k] for k in
                     ("total", "replicated", "analogue", "not_replicable")},
        "cannot_replicate": {p.key: p.reason for p in PANELS
                             if p.status == "not_replicable"},
        "analogues": {p.key: p.reason for p in PANELS
                      if p.status == "analogue"},
        "panels": replicate_all(),
    }


def analysis_guo2017(st: Structure, species: str, **kw) -> dict:
    """Guo & MacKinnon 2017 panel by panel: what reproduces and what cannot.

    Returns the coverage summary plus every panel that needs no structure and
    the ones this entry supports, rather than the full ``replicate_all`` dump —
    the report is meant to be read, and the per-panel detail is a CLI or an
    API call away.
    """
    from .guo2017 import PANELS, coverage, replicate

    summary = coverage()
    results, failures = {}, {}
    for panel in PANELS:
        if panel.compute is None:
            continue
        try:
            results[panel.key] = replicate(panel.key, structure=st,
                                           reference=species)["result"]
        except Exception as exc:                          # noqa: BLE001
            failures[panel.key] = f"{type(exc).__name__}: {exc}"
    return {
        "paper": summary["paper"],
        "coverage": summary["summary"],
        "by_status": summary["by_status"],
        "cannot_replicate": {p.key: p.reason for p in PANELS
                             if p.status == "not_replicable"},
        "analogue_not_the_same_quantity": {p.key: p.reason for p in PANELS
                                           if p.status == "analogue"},
        "panels": results,
        "failed": failures,
    }


def analysis_interactions(st: Structure, species: str, **kw) -> dict:
    from .interactions import detect_interactions
    contacts = detect_interactions(st, min_sequence_separation=3)
    return {"counts": contacts.counts(),
            "disulfides": [str(i) for i in contacts.of_kind("disulfide")],
            "caveat": contacts.meta["caveat"]}


#: Name to function. The CLI and the report share this registry, so a new
#: analysis becomes available in both at once.
from .report_tags import (analysis_hybrid,   # noqa: E402
                          analysis_fusion, analysis_labelling,
                          analysis_ligands, analysis_nanodomain,
                          analysis_paired_variant,
                          analysis_permeation, analysis_prediction_record)
from .report_validation import (analysis_fluctuations,   # noqa: E402
                                analysis_homology, analysis_paralogue)

ANALYSES = {
    "dome": analysis_dome,
    "permeation": analysis_permeation,
    "nanodomain": analysis_nanodomain,
    "prediction_record": analysis_prediction_record,
    "ligands": analysis_ligands,
    "paired_variant": analysis_paired_variant,
    "fusion": analysis_fusion,
    "hybrid": analysis_hybrid,
    "labelling": analysis_labelling,
    "pore": analysis_pore,
    "hydration": analysis_hydration,
    "modes": analysis_modes,
    "fluctuations": analysis_fluctuations,
    "paralogue": analysis_paralogue,
    "homology": analysis_homology,
    "pockets": analysis_pockets,
    "interactions": analysis_interactions,
    "guo2017": analysis_guo2017,
    "liu2025": analysis_liu2025,
}


def build_report(structure: Structure, species: str | None = None,
                 analyses: list[str] | None = None,
                 parameters: dict | None = None,
                 title: str | None = None) -> AnalysisReport:
    """Run the named analyses and package them with provenance."""
    record = load_registry().get(structure.name)
    species = species or (record.numbering_species if record else "human")
    wanted = analyses or list(ANALYSES)

    prov = collect_provenance(structure, parameters)
    prov.parameters.setdefault("species_numbering", species)
    report = AnalysisReport(provenance=prov,
                            title=title or f"PIEZO1 analysis — {structure.name}")

    for name in wanted:
        fn = ANALYSES.get(name)
        if fn is None:
            report.provenance.warnings.append(f"unknown analysis {name!r}")
            continue
        try:
            report.add(name, fn(structure, species, **(parameters or {})))
        except Exception as exc:
            report.add(name, {"error": f"{type(exc).__name__}: {exc}"})
            report.provenance.warnings.append(f"{name} failed: {exc}")
    return report
