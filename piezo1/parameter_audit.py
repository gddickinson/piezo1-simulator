"""Find numbers in the scientific modules that are not in the registry.

The rule this enforces: **any number a calculation depends on must be a
registered parameter, carrying a unit, bounds and a citation.** A constant
written into a function default is invisible — it cannot be listed, shown to a
user, or traced to a paper — and this project has already had to correct several
numbers that were invisible in exactly that way.

Not every literal is a scientific parameter. A convergence tolerance, an
iteration cap, a random seed and a zero-initialised dataclass field are
implementation details, and pretending otherwise would bury the real parameters
in noise. Those are **exempt, individually and with a stated reason** — the
point is that the exemption is a decision someone made and recorded, not an
oversight.

Run standalone::

    python -m piezo1.parameter_audit
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path

from .config import PROJECT_ROOT
from .parameters import PARAMETERS

__all__ = ["Finding", "audit", "SCANNED_PACKAGES", "EXEMPT"]

#: Packages whose numbers are scientific. `render` and `ui` are presentation and
#: are not scanned; `io` is transport.
SCANNED_PACKAGES = ["piezo1/physics", "piezo1/structure", "piezo1/analysis"]

#: Argument and field names that are implementation details wherever they occur,
#: with the reason. Keyed by name because the reason is the same every time.
EXEMPT_NAMES = {
    "tol": "convergence tolerance, not a physical quantity",
    "tolerance": "convergence tolerance",
    "atol": "absolute tolerance",
    "rtol": "relative tolerance",
    "seed": "random seed; results are reported with it fixed",
    "iterations": "iteration cap",
    "maxiter": "iteration cap",
    "max_iterations": "iteration cap",
    "relaxation": "damping on a fixed-point iteration; changes how fast it "
                  "converges, not what it converges to",
    "calcium": "an experimental condition set by the caller, not a parameter",
    "max_nodes": "solver node cap",
    "max_pairs": "sampling cap, reported when it truncates",
    "max_n": "search ceiling",
    "n_max": "search ceiling; returning it would be reported as no answer",
    "ROUND": "which recorded validation a design analysis refers to; an "
             "identifier, not a quantity — the numbers are read from the record",
    "chunk": "memory blocking; cannot change a result",
    "n_points": "quadrature or plotting density",
    "n_samples": "Monte-Carlo sample count",
    "n_simulations": "Monte-Carlo replicate count",
    "n_resamples": "bootstrap replicate count; changes the precision of an "
                   "interval, never the quantity being estimated",
    "n_pairs": "how many site pairs are sampled; changes the precision of the "
               "correlation, not its expected value",
    "digits": "how many decimals to print; presentation only",
    "level": "derived from stats.alpha rather than set independently",
    "n_bins": "histogram resolution",
    "n_angular": "mesh resolution for rendering",
    "n_frames": "animation length",
    "n_rays": "ray-cast sampling density",
    "n_protomers": "structural fact fixed by the assembly, registered as anm.n_protomers",
    "n_channels": "simulation multiplicity chosen by the caller",
    "n_points_": "plotting density",
    "n": "grid or sample size",
    "size": "page size for a paged API",
    "timeout": "network timeout",
    "delay": "polite rate limit for a public API",
    "coarse": "initial grid resolution before refinement",
    "refine_steps": "refinement iterations",
    "indent": "text formatting",
    "top": "how many results to show",
    "pc": "which principal component to report",
    "scale": "display scaling",
    "amplitude": "display amplitude; explicitly not a physical prediction",
    "min_ca": "chain-length floor, registered as geometry.min_ca_per_protomer",
    "min_common": "coverage floor chosen per analysis",
    "min_length": "sequence-length filter for ortholog search",
    "max_length": "sequence-length filter for ortholog search",
    "min_spheres": "cluster-size floor for reporting",
    "max_pockets": "how many pockets to report",
    "sensitivity": "tunable knob, reported with every prediction",
    "resting": "an experimental condition set by the caller",
    "resting_tension": "an experimental condition set by the caller",
    "duration": "protocol length set by the caller",
    "tension": "an independent variable, not a parameter",
    "lo": "search bracket",
    "hi": "search bracket",
    "ratio": "group-size ratio chosen per design",
    "alpha_": "see stats.alpha",
    "n_modes_": "see anm.n_modes",
}

#: Specific (module, owner, name) triples exempt for reasons of their own.
EXEMPT = {
    ("structure/frame.py", "Frame", "n_atoms_fitted"): "counter, initialised to zero",
    ("structure/fusion.py", "AccessibleVolume", "n_before_clash"):
        "counter, initialised to zero",
    ("physics/permeation.py", None, "F_FARADAY"):
        "SI-definitional since the 2019 redefinition; a unit conversion, not a "
        "measured quantity that could be revised",
    ("physics/permeation.py", None, "R_GAS"):
        "SI-definitional since the 2019 redefinition",
    ("analysis/variant_structures.py", "VariantStructure", "n_protein_atoms"):
        "counter, initialised to zero",
    ("physics/permeation.py", "PermeationResult", "access_ohm"): "computed field",
    ("physics/permeation.py", "PermeationResult", "pore_ohm"): "computed field",
    ("analysis/gnomad.py", None, "window"):
        "half-width of the window over which observed missense variation is "
        "averaged. A smoothing scale for a population-genetics readout, not a "
        "property of PIEZO1; frozen at 25 by "
        "docs/PREREGISTRATION_ROUND41.md and tested at 10 and 50 as "
        "pre-registered secondaries so no result can rest on the choice.",
    ("analysis/harvest.py", "HarvestReport", "n_papers"):
        "counter, initialised to zero",
    ("analysis/harvest.py", "_sentence_around", "width"):
        "how much of the surrounding sentence to keep for a human to read; "
        "presentation only, and it cannot change which candidates are found",
    ("analysis/prediction_confidence.py", None, "n_residues"):
        "length of human PIEZO1, a fact about the sequence rather than a "
        "parameter; the same 2521 core.sequence carries",
    ("analysis/prediction_confidence.py", "saturation", "min_separation"):
        "sequence separation beyond which PAE is read as long-range. A "
        "reporting threshold for a confidence readout, not a physical "
        "quantity; the separation-binned table beside it shows the whole "
        "curve so no conclusion rests on the single cut.",
    ("analysis/prediction_confidence.py", "assess_seam", "seam"):
        "where a hybrid model would be cut, which is set by which residues the "
        "experimental structures resolve rather than chosen",
    ("analysis/gnomad.py", "missense_density", "n_residues"):
        "length of human PIEZO1, a fact about the sequence rather than a "
        "parameter of any calculation; the same 2521 that "
        "core.sequence carries.",
    ("analysis/model_error.py", "pore_convention_error", "uniform_radius"):
        "the radius the *alternative* convention gives every atom. It defines "
        "that convention rather than describing PIEZO1, and the whole point of "
        "the check is that the answer depends on it — which is measured and "
        "reported, not fixed to one value.",
    ("analysis/crosscheck_methods.py", "conservation_by_kmer_anchoring", "k"):
        "seed length for the DP-free anchoring used only as a cross-check. It "
        "is a property of the alternative instrument, not of PIEZO1, and the "
        "pipeline's conservation does not depend on it.",
    ("analysis/crosscheck_methods.py", "conservation_by_kmer_anchoring",
     "window"):
        "how far the seed may slide when anchoring; a search bound for the "
        "cross-check route only",
    ("analysis/labelling.py", "predicted_brightness", "per_dye_intensity"):
        "defines the brightness unit rather than measuring anything: the "
        "histogram is in units of one dye, so 1.0 IS the unit",
    ("analysis/labelling.py", "predicted_brightness", "background"):
        "zero-offset default; a real background is supplied by the caller in "
        "whatever units their camera reports",
    ("structure/frame.py", None, "CTERM_FRACTION"):
        "which slice of residues counts as the cytosolic end when orienting a "
        "structure. A criterion for reading topology off a model, not a "
        "measured quantity — there is no experiment that returns it and no "
        "paper to cite. Its value is justified in the module against all 20 "
        "downloaded entries and pinned by "
        "test_frame.test_no_downloaded_structure_loads_upside_down.",
    ("structure/protomers.py", None, "MIN_CA_PER_PROTOMER"):
        "how many C-alphas a chain needs before it counts as a protomer rather "
        "than a bound peptide. A parsing threshold for deposited files, not a "
        "property of PIEZO1; covered by test_entities.",
    ("physics/kinetics.py", "GatingModel", "conductance_pS"):
        "registered as kinetics.conductance_pS; the field mirrors it",
    ("analysis/measure.py", "SASAResult", "probe"):
        "record of the probe used, not the default that chose it",
    ("analysis/measure.py", "SASAResult", "n_points"):
        "record of the density used",
    ("analysis/pockets.py", "AlphaSpheres", "n_total"): "counter, initialised to zero",
    ("analysis/pockets.py", "Pocket", "buriedness"): "computed field, zero-initialised",
    ("analysis/validation.py", "EffectSize", "n_bootstrap"): "records what was used",
    ("analysis/conservation.py", "ConservationProfile", "n_orthologs"):
        "counter filled in after the alignment; zero means nothing was fetched",
    ("analysis/variant_impact.py", "CouplingScore", "gating_cost_change"): "computed field",
    ("analysis/variant_impact.py", "CouplingScore", "cost_change_normalised"): "computed field",
    ("analysis/variant_impact.py", "CouplingScore", "spring_scale"): "computed field",
    ("analysis/variant_impact.py", "CouplingScore", "n_contacts"): "computed field",
    ("analysis/variant_impact.py", "CouplingScore", "local_strain"): "computed field",
    ("analysis/hydration.py", "LiningPoint", "distance"): "computed field",
    ("analysis/hydration.py", None, "HYDROPHOBICITY_FALLBACK"):
        "CHAP's fallback for residues outside the scale; a definition, not a measurement",
    ("physics/elastica.py", "ElasticaSolution", "r0"): "computed field",
    ("physics/elastica.py", "ElasticaSolution", "slope"): "computed field",
    ("physics/membrane.py", "FootprintSolution", "r0"): "computed field",
    ("physics/membrane.py", "FootprintSolution", "slope"): "computed field",
    ("physics/membrane.py", "FootprintSolution", "energy"): "computed field",
    ("physics/dome.py", "DomeGeometrySummary", "dome_depth"): "computed field",
    ("physics/dome.py", "DomeGeometrySummary", "footprint_radius"): "computed field",
    ("structure/geometry.py", "DomeGeometry", "radius_of_curvature"): "computed field",
    ("structure/geometry.py", "DomeGeometry", "dome_depth"): "computed field",
    ("structure/geometry.py", "DomeGeometry", "footprint_radius"): "computed field",
    ("structure/geometry.py", "DomeGeometry", "dome_area"): "computed field",
    ("structure/geometry.py", "DomeGeometry", "projected_area"): "computed field",
    ("structure/geometry.py", "DomeGeometry", "n_atoms_used"):
        "counter recording how many points survived the sphere-fit trim",
    ("structure/superpose.py", "SymmetryAxis", "order"):
        "the rotational order being tested, an argument of the question",
    ("structure/superpose.py", "SymmetryAxis", "rmsd"): "computed field",
    ("structure/superpose.py", "SymmetryAxis", "angle_deg"): "computed field",
    ("analysis/ensemble.py", "StructureEnsemble", "n_protomers"): "structural fact",
    ("structure/morph.py", None, "max_bond_jump"):
        "a diagnostic threshold for reporting bond distortion, not an input",
    ("structure/pore.py", None, "min_separation"):
        "how far apart two reported constrictions must be; presentation",
    ("physics/kinetics.py", None, "peak_open_probability"): "search bracket",
    ("physics/anm.py", None, "mode"): "display amplitude",
}

#: Registry keys the scanner should recognise for a given (module, name).
MAPPED = {
    ("physics/membrane.py", "kappa"): "membrane.kappa",
    ("physics/membrane.py", "tension"): "membrane.tension",
    ("physics/dome.py", "delta_area"): "dome.delta_area",
    ("physics/dome.py", "delta_g0"): "dome.delta_g0",
    ("physics/anm.py", "cutoff"): "anm.cutoff",
    ("physics/anm.py", "gamma"): "anm.gamma",
    ("physics/anm.py", "d0"): "anm.d0",
    ("physics/anm.py", "n_modes"): "anm.n_modes",
    ("analysis/features.py", "cutoff"): "anm.cutoff",
    ("analysis/features.py", "n_modes"): "anm.n_modes",
    ("analysis/variant_impact.py", "cutoff"): "anm.cutoff",
    ("analysis/variant_impact.py", "d0"): "anm.d0",
    ("analysis/variant_impact.py", "gamma"): "anm.gamma",
    ("structure/pore.py", "step"): "pore.step",
    ("analysis/model_error.py", "n_modes"): "anm.n_modes",
    ("analysis/model_error.py", "cutoff"): "anm.cutoff",
    ("analysis/model_error.py", "step"): "pore.step",
    ("analysis/crosscheck_methods.py", "leash"): "pore.leash",
    ("analysis/crosscheck_methods.py", "search"): "pore.search",
    ("analysis/crosscheck_methods.py", "probe"): "sasa.probe_radius",
    ("structure/pore.py", "leash"): "pore.leash",
    ("structure/pore.py", "search"): "pore.search",
    ("structure/pore.py", "threshold"): "pore.constriction_threshold",
    ("structure/pore.py", "ion_radius"): "pore.ion_radius",
    ("structure/geometry.py", "trim"): "geometry.sphere_trim",
    ("analysis/measure.py", "probe"): "sasa.probe_radius",
    ("analysis/pockets.py", "r_min"): "pockets.r_min",
    ("analysis/pockets.py", "r_max"): "pockets.r_max",
    ("analysis/pockets.py", "min_neighbours"): "pockets.min_neighbours",
    ("analysis/pockets.py", "neighbour_radius"): "pockets.neighbour_radius",
    ("analysis/pockets.py", "cutoff"): "pockets.ligand_cutoff",
    ("analysis/allostery.py", "contact_cutoff"): "allostery.contact_cutoff",
    ("analysis/allostery.py", "min_correlation"): "allostery.min_correlation",
    ("analysis/conservation.py", "taxon"): "conservation.taxon",
    ("analysis/conservation.py", "min_coverage"): "conservation.min_coverage",
    ("analysis/design.py", "alpha"): "stats.alpha",
    ("analysis/design.py", "target_power"): "stats.target_power",
    ("analysis/design.py", "n_permutations"): "stats.n_permutations",
    ("analysis/validation.py", "alpha"): "stats.alpha",
    ("analysis/validation.py", "n_permutations"): "stats.n_permutations",
    ("analysis/validation.py", "n_bootstrap"): "stats.n_bootstrap",
    ("analysis/interactions.py", "min_sequence_separation"):
        "interactions.min_sequence_separation",
    ("analysis/hydration.py", "max_radius"): "hydration.max_radius",
    ("analysis/conservation.py", "conservation_threshold"): "conservation.constrained_threshold",
    ("analysis/features.py", "sasa_points"): "sasa.n_points_fast",
    ("analysis/measure.py", "radius"): "measure.hydrophobicity_radius",
    ("analysis/pockets.py", "reach"): "pockets.buriedness_reach",
    ("analysis/pockets.py", "clearance"): "pockets.buriedness_clearance",
    ("analysis/pockets.py", "cluster_distance"): "pockets.cluster_distance",
    ("analysis/pockets.py", "lining_cutoff"): "pockets.lining_cutoff",
    ("analysis/report.py", "step"): "pore.step",
    ("analysis/report_tags.py", "step"): "pore.step",
    ("analysis/variant_structures.py", "step"): "pore.step",
    ("analysis/report.py", "n_modes"): "anm.n_modes",
}


@dataclass
class Finding:
    """One numeric literal the audit could not account for."""

    module: str
    owner: str | None
    name: str
    value: float

    @property
    def where(self) -> str:
        return f"{self.module}: {self.owner + '.' if self.owner else ''}{self.name}"

    def __str__(self) -> str:
        return f"{self.where} = {self.value!r}"


def _numeric(node) -> bool:
    return (isinstance(node, ast.Constant)
            and isinstance(node.value, (int, float))
            and not isinstance(node.value, bool))


def _scan_module(path: Path, relative: str) -> list[Finding]:
    tree = ast.parse(path.read_text())
    found: list[Finding] = []

    for node in tree.body:
        if isinstance(node, ast.Assign) and _numeric(node.value):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id.isupper():
                    found.append(Finding(relative, None, target.id,
                                         node.value.value))

    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            for statement in node.body:
                if (isinstance(statement, ast.AnnAssign)
                        and statement.value is not None
                        and _numeric(statement.value)):
                    found.append(Finding(relative, node.name,
                                         statement.target.id,
                                         statement.value.value))
        elif isinstance(node, ast.FunctionDef):
            positional = node.args.posonlyargs + node.args.args
            tail = positional[len(positional) - len(node.args.defaults):]
            for argument, default in zip(tail, node.args.defaults):
                if _numeric(default):
                    found.append(Finding(relative, node.name, argument.arg,
                                         default.value))
            for argument, default in zip(node.args.kwonlyargs,
                                         node.args.kw_defaults):
                if default is not None and _numeric(default):
                    found.append(Finding(relative, node.name, argument.arg,
                                         default.value))
    return found


def _accounted(finding: Finding) -> bool:
    """Whether this literal is registered, mapped or explicitly exempt."""
    if finding.name in EXEMPT_NAMES:
        return True
    if (finding.module, finding.owner, finding.name) in EXEMPT:
        return True
    if (finding.module, None, finding.name) in EXEMPT:
        return True
    key = MAPPED.get((finding.module, finding.name))
    if key is not None:
        return key in PARAMETERS
    # A dataclass field or argument named exactly like a registry key.
    return finding.name in PARAMETERS


def audit(root: Path | None = None) -> list[Finding]:
    """Every numeric literal in the scientific packages that is unaccounted for."""
    root = Path(root or PROJECT_ROOT)
    unaccounted: list[Finding] = []
    for package in SCANNED_PACKAGES:
        directory = root / package
        if not directory.exists():
            continue
        for path in sorted(directory.glob("*.py")):
            relative = str(path.relative_to(root / "piezo1"))
            for finding in _scan_module(path, relative):
                if not _accounted(finding):
                    unaccounted.append(finding)
    return unaccounted


def main() -> int:
    findings = audit()
    if not findings:
        print(f"clean: every number in {', '.join(SCANNED_PACKAGES)} is "
              f"registered or explicitly exempt "
              f"({len(PARAMETERS)} registered parameters)")
        return 0
    print(f"{len(findings)} unaccounted numeric literal(s):\n")
    for finding in findings:
        print(f"  {finding}")
    print("\nEach must either be registered in scripts/build_parameters.py "
          "with a unit and a citation, or added to EXEMPT/EXEMPT_NAMES in "
          "piezo1/parameter_audit.py with the reason it is not a scientific "
          "parameter.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
