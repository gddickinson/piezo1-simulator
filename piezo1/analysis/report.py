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

from .. import __version__
from ..core.structure import Structure
from ..io.registry import load_registry
from ..parameters import PARAMETERS

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

def _protomer_blocks(st: Structure):
    chains = []
    for c in st.chains:
        m = st.mask_ca() & (st.chain == c)
        if m.sum() > 300:
            chains.append((st.xyz[m], st.res_seq[m]))
    if len(chains) < 3:
        return None, None
    common = set(chains[0][1].tolist())
    for _, seq in chains[1:3]:
        common &= set(seq.tolist())
    arr = np.array(sorted(common))
    return ([x[np.searchsorted(s, arr)].astype(float) for x, s in chains[:3]],
            arr)


def analysis_dome(st: Structure, species: str, **kw) -> dict:
    import json as _json
    from ..config import RESOURCE_DIR
    from ..structure.geometry import measure_dome
    blocks, _ = _protomer_blocks(st)
    if blocks is None:
        return {"error": "needs three well-resolved protomers"}
    tms = _json.loads((RESOURCE_DIR / f"uniprot_{species}.json").read_text())["transmembrane"]
    pts = []
    for c in st.chains:
        m = st.mask_ca() & (st.chain == c)
        if m.sum() < 300:
            continue
        xyz, seq = st.xyz[m], st.res_seq[m]
        for tm in tms:
            mid = 0.5 * (tm["start"] + tm["end"])
            half = max(2.0, (tm["end"] - tm["start"]) / 6.0)
            sel = (seq >= mid - half) & (seq <= mid + half)
            if sel.sum() >= 3:
                pts.append(xyz[sel].mean(axis=0))
    dome = measure_dome(blocks, np.array(pts))
    return {"radius_of_curvature_nm": dome.radius_of_curvature / 10,
            "dome_depth_nm": dome.dome_depth / 10,
            "dome_area_nm2": dome.dome_area / 100,
            "projected_area_nm2": dome.projected_area / 100,
            "excess_area_nm2": dome.excess_area / 100,
            "c3_angle_deg": dome.notes["c3_angle_deg"],
            "c3_rmsd_A": dome.notes["c3_rmsd"],
            "reference": "closed-state R_c = 10.2 nm "
                         "(Haselwandter & MacKinnon 2018)"}


def analysis_pore(st: Structure, species: str, step: float = 1.5, **kw) -> dict:
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    blocks, _ = _protomer_blocks(st)
    if blocks is None:
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


def analysis_hydration(st: Structure, species: str, step: float = 1.0,
                       **kw) -> dict:
    """Hydrophobic-gating prediction (Rao et al. 2019) for this structure."""
    from .hydration import load_grid, predict_wetting
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    blocks, _ = _protomer_blocks(st)
    if blocks is None:
        return {"error": "needs three well-resolved protomers"}
    grid = load_grid()
    if not grid.available:
        return {"error": "CHAP grid not downloaded; run python -m piezo1.io.fetch"}
    prof = pore_profile(st, detect_c3_axis(blocks), step=step)
    pred = predict_wetting(st, prof, grid)
    return {"score": pred.score,
            "cutoff": pred.meta["cutoff"],
            "bottleneck_radius_nm": pred.min_radius / 10.0,
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
    blocks, _ = _protomer_blocks(st)
    if blocks is None:
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


def analysis_interactions(st: Structure, species: str, **kw) -> dict:
    from .interactions import detect_interactions
    contacts = detect_interactions(st, min_sequence_separation=3)
    return {"counts": contacts.counts(),
            "disulfides": [str(i) for i in contacts.of_kind("disulfide")],
            "caveat": contacts.meta["caveat"]}


#: Name to function. The CLI and the report share this registry, so a new
#: analysis becomes available in both at once.
from .report_tags import (analysis_fusion, analysis_labelling,   # noqa: E402
                          analysis_ligands, analysis_nanodomain,
                          analysis_permeation, analysis_prediction_record)

ANALYSES = {
    "dome": analysis_dome,
    "permeation": analysis_permeation,
    "nanodomain": analysis_nanodomain,
    "prediction_record": analysis_prediction_record,
    "ligands": analysis_ligands,
    "fusion": analysis_fusion,
    "labelling": analysis_labelling,
    "pore": analysis_pore,
    "hydration": analysis_hydration,
    "modes": analysis_modes,
    "pockets": analysis_pockets,
    "interactions": analysis_interactions,
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
