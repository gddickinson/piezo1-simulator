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
oversight. The exemptions themselves live in
:mod:`piezo1.parameter_audit_exemptions`, split off at the 500-line limit along
the seam this file already had: here is the mechanism, there is the judgement,
and the judgement is the half a reader has to check.

Run standalone::

    python -m piezo1.parameter_audit
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path

from .config import PROJECT_ROOT
from .parameters import PARAMETERS


from .parameter_audit_exemptions import EXEMPT, EXEMPT_NAMES  # noqa: F401
__all__ = ["Finding", "audit", "SCANNED_PACKAGES", "EXEMPT"]

#: Packages whose numbers are scientific. `render` and `ui` are presentation and
#: are not scanned; `io` is transport.
SCANNED_PACKAGES = ["piezo1/physics", "piezo1/structure", "piezo1/analysis"]

#: Argument and field names that are implementation details wherever they occur,
#: with the reason. Keyed by name because the reason is the same every time.
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
