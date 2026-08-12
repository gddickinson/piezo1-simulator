"""Headless command line for batch analysis.

Everything the GUI can compute, scriptable and reproducible::

    python -m piezo1.cli list
    python -m piezo1.cli dome 8YEZ
    python -m piezo1.cli pore 8YEZ --step 1.0
    python -m piezo1.cli modes 8YEZ --n-modes 30
    python -m piezo1.cli pockets 8YEZ --top 5
    python -m piezo1.cli variants --classification GoF
    python -m piezo1.cli conservation --top 20
    python -m piezo1.cli report 8YEZ -o report.md
    python -m piezo1.cli batch --state curved -o results/

Every command accepts ``--json`` to emit machine-readable output instead of a
table, so results can be piped straight into other tools.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from .analysis.report import ANALYSES, build_report
from .config import STRUCTURE_DIR
from .core.annotations import load_annotations
from .core.structure import Structure
from .io.registry import load_registry

__all__ = ["main", "build_parser"]


def _load(pdb: str) -> Structure:
    record = load_registry().get(pdb)
    path = record.path if record and record.available else STRUCTURE_DIR / f"{pdb}.cif"
    if not Path(path).exists():
        raise SystemExit(f"{pdb} not found. Run `python -m piezo1.io.fetch` "
                         f"or use `python -m piezo1.cli list`.")
    return Structure.from_file(path)


def _emit(data, as_json: bool, renderer=None) -> None:
    if as_json:
        print(json.dumps(data, indent=1, default=_default))
    elif renderer is not None:
        renderer(data)
    else:
        print(json.dumps(data, indent=1, default=_default))


def _default(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------

def cmd_list(args) -> int:
    registry = load_registry()
    rows = [r for r in registry
            if (not args.species or r.species == args.species)
            and (not args.state or r.state == args.state)]
    if args.json:
        _emit([{"pdb": r.pdb, "species": r.species, "state": r.state,
                "resolution": r.resolution, "available": r.available,
                "recommended_for": list(r.recommended_for)} for r in rows], True)
        return 0
    print(f"{'PDB':6s} {'species':7s} {'state':13s} {'res':>5s} {'':1s} note")
    for r in rows:
        mark = " " if r.available else "!"
        print(f"{r.pdb:6s} {r.species:7s} {r.state:13s} "
              f"{str(r.resolution):>5s} {mark} {r.note[:60]}")
    missing = [r.pdb for r in rows if not r.available]
    if missing:
        print(f"\n! not downloaded ({len(missing)}): {' '.join(missing)}")
        print("  run: python -m piezo1.io.fetch")
    return 0


def cmd_analysis(args) -> int:
    """Run one named analysis from the shared registry."""
    st = _load(args.structure)
    record = load_registry().get(args.structure)
    species = args.species or (record.numbering_species if record else "human")
    params = {k: v for k, v in vars(args).items()
              if k in ("step", "n_modes", "top") and v is not None}
    result = ANALYSES[args.command](st, species, **params)
    _emit(result, args.json)
    return 0 if "error" not in result else 1


def cmd_liu2025(args) -> int:
    """Liu et al. 2025 panel by panel, including what cannot be reproduced.

    Most panels compare their four deposited states, so ``structure`` is
    accepted for symmetry with the other commands and used by none of them —
    a single entry can only ever supply one side of a comparison.
    """
    from .analysis.liu2025 import (PANELS, coverage, not_replicable,
                                   panel_by_key, replicate)

    counts = {k: coverage()[k] for k in
              ("total", "replicated", "analogue", "not_replicable")}
    if args.panel:
        _emit(replicate(panel_by_key(args.panel).key), args.json)
    elif args.coverage:
        _emit(counts, args.json)
    else:
        _emit({
            "paper": coverage()["paper"],
            "coverage": counts,
            "panels": {p.key: {"label": p.label, "status": p.status,
                               "title": p.title} for p in PANELS},
            "cannot_replicate": {p.key: p.reason for p in not_replicable()},
        }, args.json)
    return 0


def cmd_export(args) -> int:
    """Write a computed per-residue scalar into a PDB B-factor column.

    The point of the command is that the number leaves. Every scalar this
    project computes was reachable only from inside it or from a JSON blob,
    and the standard interoperability route is this file.
    """
    from .core.export import write_scalar_pdb

    st = _load(args.structure)
    scalars = _scalar_sources()
    if args.scalar not in scalars:
        print(f"unknown scalar {args.scalar!r}; choose from "
              f"{', '.join(sorted(scalars))}")
        return 2
    values, default_scale, label = scalars[args.scalar](st, args)
    if not values:
        print(f"{args.scalar} scored nothing on {st.name}")
        return 1
    scale = args.scale if args.scale is not None else default_scale
    report = write_scalar_pdb(st, values, args.out, scale=scale, name=label)
    print(report.summary())
    return 0


def _scalar_sources() -> dict:
    """The per-residue scalars that can be exported, and their natural scale.

    A registry rather than a chain of ifs, so adding one is a line here and the
    CLI help lists it automatically.
    """
    def wetting(st, args):
        from .analysis.hydration import load_grid, predict_wetting
        from .structure.pore import pore_profile
        from .structure.protomers import protomer_blocks
        from .structure.superpose import detect_c3_axis

        blocks, _ = protomer_blocks(st)
        profile = pore_profile(st, detect_c3_axis(blocks))
        prediction = predict_wetting(st, profile, grid=load_grid())
        return ({p.residue: p.energy for p in prediction.points}, 10.0,
                "wetting energy kJ/mol")

    def conservation(st, args):
        from .analysis.conservation import conservation_profile

        profile = conservation_profile()
        return ({int(r): float(v) for r, v in
                 zip(profile.residues, profile.conservation)}, 100.0,
                "conservation")

    def fluctuation(st, args):
        from .analysis.fluctuations import predicted_msf

        values = predicted_msf(st)
        return ({int(r): float(v) for r, v in values.items()}, 1.0,
                "predicted mean-square fluctuation")

    return {"wetting": wetting, "conservation": conservation,
            "fluctuation": fluctuation}


def cmd_guo2017(args) -> int:
    """Replicate Guo & MacKinnon 2017 panel by panel."""
    from .analysis.guo2017 import (PANELS, coverage, panel_by_key, replicate,
                                   replicate_all)

    if args.coverage:
        report = coverage()
        if args.json:
            _emit(report, True)
            return 0
        print(report["summary"])
        print()
        for panel in PANELS:
            mark = {"replicated": "yes", "analogue": "~  ",
                    "not_replicable": "no "}[panel.status]
            print(f"  {mark} {panel.label:22s} {panel.title}")
            if panel.status != "replicated":
                print(f"        {panel.reason.split('.')[0]}.")
        return 0

    structure = None
    needed = ([panel_by_key(args.panel)] if args.panel else PANELS)
    if any(p.needs_structure and p.compute is not None for p in needed):
        structure = _load(args.structure)

    if args.panel:
        result = replicate(args.panel, structure=structure,
                           reference=args.species)
    else:
        result = replicate_all(structure=structure)
    _emit(result, args.json)
    return 0


def cmd_variants(args) -> int:
    ann = load_annotations("human")
    rows = [v for v in ann.variants
            if not args.classification or v.classification == args.classification]
    if args.residue:
        rows = [v for v in rows if v.residue == args.residue]
    if args.json:
        _emit([{"label": v.label, "residue": v.residue,
                "classification": v.classification, "domain": v.domain,
                "phenotype": v.phenotype, "pmid": v.pmid,
                "modelled_in": list(v.modelled_in)} for v in rows], True)
        return 0
    print(f"{'variant':10s} {'class':12s} {'domain':16s} {'resolved in':22s} phenotype")
    for v in rows:
        where = ",".join(v.modelled_in[:3]) or "none"
        print(f"{v.label:10s} {str(v.classification):12s} {str(v.domain)[:16]:16s} "
              f"{where:22s} {str(v.phenotype)[:44]}")
    n_missing = sum(1 for v in rows if not v.modelled_in)
    print(f"\n{len(rows)} variants; {n_missing} resolved in no human structure")
    return 0


def cmd_conservation(args) -> int:
    from .analysis.conservation import conservation_profile, load_orthologs
    try:
        orthologs = load_orthologs()
    except FileNotFoundError as exc:
        raise SystemExit(str(exc))
    profile = conservation_profile(orthologs)
    ann = load_annotations("human")
    if args.json:
        _emit({"n_orthologs": profile.n_orthologs,
               "domain_means": profile.domain_means(ann),
               "top_conserved": profile.top_conserved(args.top)}, True)
        return 0
    print(f"{profile.n_orthologs} orthologs\n")
    print("mean conservation by domain:")
    for k, v in sorted(profile.domain_means(ann).items(), key=lambda kv: -kv[1]):
        print(f"   {k:14s} {v:.3f}")
    return 0


def cmd_report(args) -> int:
    st = _load(args.structure)
    params = {k: v for k, v in vars(args).items()
              if k in ("step", "n_modes", "top") and v is not None}
    report = build_report(st, species=args.species,
                          analyses=args.analyses, parameters=params)
    out = Path(args.output) if args.output else None
    if out is None:
        print(report.to_markdown())
    elif out.suffix == ".json":
        report.to_json(out)
        print(f"wrote {out}")
    else:
        report.to_markdown(out)
        report.to_json(out.with_suffix(".json"))
        print(f"wrote {out} and {out.with_suffix('.json')}")
    return 0


def cmd_batch(args) -> int:
    """Run a report over every matching structure."""
    registry = load_registry()
    rows = [r for r in registry.available()
            if (not args.species or r.species == args.species)
            and (not args.state or r.state == args.state)
            and r.state != "fragment"]
    if not rows:
        raise SystemExit("no structures match")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    summary = []
    for record in rows:
        print(f"  {record.pdb} ...", end=" ", flush=True)
        try:
            st = Structure.from_file(record.path)
            report = build_report(st, analyses=args.analyses)
            report.to_json(out / f"{record.pdb}.json")
            dome = report.results.get("dome", {})
            pore = report.results.get("pore", {})
            summary.append({"pdb": record.pdb, "state": record.state,
                            "radius_of_curvature_nm": dome.get("radius_of_curvature_nm"),
                            "bottleneck_A": pore.get("bottleneck_radius_A"),
                            "conductive": pore.get("conductive")})
            print("ok")
        except Exception as exc:
            print(f"FAILED {type(exc).__name__}: {exc}")
            summary.append({"pdb": record.pdb, "error": str(exc)})
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=_default))

    print(f"\n{'PDB':6s} {'state':13s} {'R_c (nm)':>9s} {'bottleneck':>11s}  conductive")
    for row in summary:
        if "error" in row:
            print(f"{row['pdb']:6s} error: {row['error'][:50]}")
            continue
        rc = row["radius_of_curvature_nm"]
        bn = row["bottleneck_A"]
        print(f"{row['pdb']:6s} {row['state']:13s} "
              f"{rc:9.1f} {bn:11.2f}  {row['conductive']}"
              if rc is not None and bn is not None else f"{row['pdb']:6s} incomplete")
    print(f"\nwrote {len(summary)} reports to {out}")
    return 0


# --------------------------------------------------------------------------
# Parser
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    # --json is offered both before and after the subcommand. Users naturally
    # type `cli dome 8YEZ --json`, and a top-level-only flag rejects that with
    # an unhelpful "unrecognized arguments" error.
    common = argparse.ArgumentParser(add_help=False)
    # SUPPRESS matters here. A subparser sharing a dest writes its default over
    # whatever the parent already parsed, so with a plain default of False,
    # `cli --json dome 8YEZ` silently comes back as json=False. With SUPPRESS
    # the attribute only appears when the flag was actually given, and main()
    # supplies the default.
    common.add_argument("--json", action="store_true",
                        default=argparse.SUPPRESS,
                        help="emit JSON instead of a table")

    ap = argparse.ArgumentParser(
        prog="python -m piezo1.cli", parents=[common],
        description="Headless PIEZO1 structural and physical analysis.",
        epilog="What this project established, and what it could not: "
               "docs/CONCLUSION.md — one page, every number traceable "
               "to the code. The variant-effect prediction it was built "
               "for does not work; five pre-registered tests, five nulls.")
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("list", help="available structures", parents=[common])
    p.add_argument("--species", choices=["human", "mouse"])
    p.add_argument("--state")
    p.set_defaults(func=cmd_list)

    for name, help_text in (
            ("dome", "membrane dome geometry"),
            ("pore", "pore radius profile and bottleneck"),
            ("hydration", "hydrophobic-gating prediction (Rao et al. 2019)"),
            ("modes", "elastic network normal modes with symmetry labels"),
            ("fluctuations", "predicted fluctuation against the deposited "
                             "B-factors, with a burial control"),
            ("paralogue", "PIEZO1 against PIEZO2: is this mechanism the "
                          "fold's rather than PIEZO1's?"),
            ("pockets", "alpha-sphere cavities and ligand contacts"),
            ("interactions", "hydrogen bonds, salt bridges, disulfides"),
            ("hybrid", "full-length model: experimental core plus the "
                       "predicted distal blade, kept apart"),
            ("fusion", "where a C-terminal HaloTag would sit (a model)"),
            ("labelling", "HaloTag labelling kinetics on the three sites"),
            ("permeation", "ion current through the pore (1-D PNP)"),
            ("nanodomain", "calcium at the HaloTag when the channel opens"),
            ("prediction_record", "what a variant prediction is worth"),
            ("ligands", "modulators and the standing of their binding sites"),
            ("paired_variant", "R2456H against wild type, with a control")):
        p = sub.add_parser(name, help=help_text, parents=[common])
        p.add_argument("structure")
        p.add_argument("--species", choices=["human", "mouse"])
        p.add_argument("--step", type=float, help="pore slice spacing, A")
        p.add_argument("--n-modes", dest="n_modes", type=int)
        p.add_argument("--top", type=int)
        p.set_defaults(func=cmd_analysis)

    p = sub.add_parser("guo2017", parents=[common],
                       help="replicate the figures of Guo & MacKinnon 2017, "
                            "the paper the dome model comes from")
    p.add_argument("structure", nargs="?", default="6B3R",
                   help="entry to measure on (default 6B3R, the paper's own)")
    p.add_argument("--panel", help="one panel key, e.g. 7-S1; omit for all")
    p.add_argument("--species", choices=["human", "mouse"])
    p.add_argument("--coverage", action="store_true",
                   help="only report how much of the paper is reachable")
    p.set_defaults(func=cmd_guo2017)

    p = sub.add_parser("export", parents=[common],
                       help="write a computed per-residue scalar into a PDB "
                            "B-factor column, for PyMOL or ChimeraX")
    p.add_argument("structure", nargs="?", default="8YEZ")
    p.add_argument("--scalar", default="wetting",
                   help="which scalar to write (wetting, conservation, "
                        "fluctuation)")
    p.add_argument("--out", default="scalar.pdb", help="output path")
    p.add_argument("--scale", type=float,
                   help="multiply before writing; the B-factor column carries "
                        "two decimals, so a 0-1 score needs scaling to use it")
    p.set_defaults(func=cmd_export)

    p = sub.add_parser("liu2025", parents=[common],
                       help="replicate the figures of Liu et al. 2025, the "
                            "paper the intermediate-open structure comes from")
    p.add_argument("structure", nargs="?", default="8IXO",
                   help="ignored for most panels, which compare all four of "
                        "their states (default 8IXO, the intermediate-open one)")
    p.add_argument("--panel", help="one panel key, e.g. 2b; omit for all")
    p.add_argument("--coverage", action="store_true",
                   help="only report how much of the paper is reachable")
    p.set_defaults(func=cmd_liu2025)

    p = sub.add_parser("variants", parents=[common], help="the curated variant table")
    p.add_argument("--classification")
    p.add_argument("--residue", type=int)
    p.set_defaults(func=cmd_variants)

    p = sub.add_parser("conservation", parents=[common], help="ortholog conservation")
    p.add_argument("--top", type=int, default=20)
    p.set_defaults(func=cmd_conservation)

    p = sub.add_parser("report", parents=[common], help="full report with provenance")
    p.add_argument("structure")
    p.add_argument("-o", "--output")
    p.add_argument("--species", choices=["human", "mouse"])
    p.add_argument("--analyses", nargs="+", choices=list(ANALYSES))
    p.add_argument("--step", type=float)
    p.add_argument("--n-modes", dest="n_modes", type=int)
    p.add_argument("--top", type=int)
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("batch", parents=[common], help="report over every matching structure")
    p.add_argument("-o", "--output", default="data/derived/reports")
    p.add_argument("--species", choices=["human", "mouse"])
    p.add_argument("--state")
    p.add_argument("--analyses", nargs="+", choices=list(ANALYSES),
                   default=["dome", "pore"])
    p.set_defaults(func=cmd_batch)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not hasattr(args, "json"):
        args.json = False
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
