"""The ``homology`` command — the family, rather than one loaded structure.

Split out of ``cli.py`` at the 500-line limit and along a real seam: every
other command there takes ``--structure`` and answers a question about one
deposited entry. These answer questions about the **family**, and two of the
three do not load coordinates at all.

Three subjects, because they are three different amounts of work and a user
should choose deliberately:

* ``--family`` — the ten members, from committed resources. Instant.
* ``--sites`` — every curated functional residue across all ten, gated by
  whether the alignment is in register there. A few seconds.
* the default matrix — every pair with its composition-matched null. **Two
  minutes**, almost all of it computing the nulls, which is the price of not
  printing a bare percentage. Said out loud before it starts.
* ``--structural A B`` — dome, helix correspondence and gating-mode overlap
  for two catalogued entries. Loads coordinates and builds two elastic
  networks.
"""

from __future__ import annotations

__all__ = ["cmd_homology", "add_homology_parser"]


def _emit(payload, as_json: bool) -> None:
    import json

    print(json.dumps(payload, indent=1, default=str) if as_json else payload)


def cmd_homology(args) -> int:
    if args.structural:
        return _structural(args)
    if args.sites:
        return _sites(args)
    if args.family:
        return _family(args)
    if getattr(args, "structure", None):
        return _for_structure(args)
    return _matrix(args)


def _for_structure(args) -> int:
    """One entry against every catalogued homologue — the registry analysis.

    Dispatches to the shared ``ANALYSES`` entry rather than reimplementing it,
    so the CLI, the report and the GUI dialog cannot disagree about what the
    comparison is. Present because every other analysis in the registry is
    reachable as ``piezo1.cli <name> <PDB>`` and a test enforces it; the
    family-level subjects above are the additions, not the default.
    """
    from .analysis.report import ANALYSES
    from .cli import _load
    from .io.registry import load_registry

    structure = _load(args.structure)
    record = load_registry().get(args.structure)
    species = (getattr(args, "species", None)
               or (record.numbering_species if record else "human"))
    result = ANALYSES["homology"](structure, species)
    if args.json:
        _emit(result, True)
        return 0 if "error" not in result else 1
    if "error" in result:
        print(result["error"])
        return 1
    print(f"{result['loaded']['pdb']} is {result['loaded']['protein']}, read "
          f"in {result['loaded']['numbering']} numbering\n")
    print(f"{'partner':9s} {'entries':>7s} {'overlap range':>15s} "
          f"{'ctrl':>6s} {'beats':>7s}  verdict")
    for row in result["homologues"]:
        span = (f"{row['gating_mode_overlap_low']:.2f}-"
                f"{row['gating_mode_overlap_high']:.2f}")
        print(f"{row['protein']:9s} {row['n_entries']:7d} {span:>15s} "
              f"{row['shuffled_control_max']:6.2f} "
              f"{row['n_beating_control']:3d}/{row['n_entries']:<3d}  "
              f"{'stable' if row['stable'] else 'NOT STABLE'}")
    print(f"\n{result['caveat']}\n\n{result['note']}")
    return 0


def _family(args) -> int:
    from .analysis.homology import family

    members = family()
    if args.json:
        _emit([vars(m) for m in members], True)
        return 0
    print(f"{'key':14s} {'accession':12s} {'protein':8s} {'len':>5s} {'TM':>3s}  "
          f"{'group':12s} organism")
    for m in members:
        print(f"{m.key:14s} {m.accession:12s} {m.protein:8s} {m.length:5d} "
              f"{m.n_transmembrane:3d}  {m.group:12s} {m.organism}")
    counts = {m.n_transmembrane for m in members}
    print(f"\n{len(members)} reviewed PIEZOs, {len({m.length for m in members})} "
          f"distinct lengths and {len(counts)} distinct helix counts "
          f"({sorted(counts)}). No two share a numbering system, and the four "
          f"non-vertebrate members do not share the 38-helix architecture the "
          f"domain table is built on — so nothing may transfer a residue "
          f"number or a helix index between them by arithmetic.")
    return 0


def _matrix(args) -> int:
    from .analysis.homology import family_matrix
    from .parameters import PARAMETERS as _P

    replicates = int(_P.value("homology.null_replicates"))
    if not args.json:
        print(f"Aligning every pair and {replicates} composition-matched "
              f"shuffles of each — about two minutes.\n"
              f"The nulls are the cost and the point: below Rost's "
              f"{_P.value('homology.twilight_identity'):.0%} line a bare "
              f"identity is not evidence.\n")
    matrix = family_matrix()
    if args.json:
        _emit({f"{a}|{b}": {
            "identity": r.identity, "identity_z": r.identity_z,
            "local_score": r.local_score, "local_z": r.local_z,
            "in_twilight_zone": r.in_twilight_zone, "verdict": r.verdict}
            for (a, b), r in matrix.relationships.items()}, True)
        return 0

    print(f"{'pair':32s} {'ident':>6s} {'null':>6s} {'z':>6s} "
          f"{'local':>7s} {'null':>5s} {'z':>7s}")
    for (a, b), r in matrix.relationships.items():
        flag = " *" if (r.local_beats_null and not r.identity_beats_null) else ""
        print(f"{a + ' vs ' + b:32s} {r.identity:6.3f} {r.null_identity:6.3f} "
              f"{r.identity_z:6.1f} {r.local_score:7.0f} {r.null_local:5.0f} "
              f"{r.local_z:7.1f}{flag}")
    print(f"\n{matrix.summary()}")
    print("* the identity is indistinguishable from a shuffled sequence of the "
          "same composition while the alignment score is overwhelming — the "
          "pairs a percentage would mislead a reader about.")
    return 0


def _sites(args) -> int:
    from .analysis.homology_sites import report

    result = report()
    targets = [r.target for r in result.rows]
    order = list(dict.fromkeys(targets))
    if args.json:
        _emit({r.group_id + "|" + r.target: {
            "identical": r.n_identical, "readable": r.n_readable,
            "conservative": r.n_conservative, "unreadable": r.n_unreliable,
            "positions": [p.summary() for p in r.positions]}
            for r in result.rows}, True)
        return 0

    print(f"{'curated group':22s}" + "".join(f"{t[:9]:>10s}" for t in order))
    for group_id in result.group_ids:
        cells = []
        for row in result.for_group(group_id):
            cells.append("  -" if not row.n_readable
                         else f"{row.n_identical}/{row.n_readable}")
        print(f"{group_id:22s}" + "".join(f"{c:>10s}" for c in cells))
    print(f"\nidentical / readable. A dash means no position in the group "
          f"lands where the alignment is in register, which is reported as "
          f"'cannot tell' and never as 'not conserved'.")
    print(f"\n{result.summary()}")
    return 0


def _structural(args) -> int:
    from .analysis.homology_structure import compare_structures

    first, second = args.structural
    comparison = compare_structures(first, second, n_modes=args.n_modes)
    if args.json:
        modes = comparison.modes
        _emit({"a": comparison.a.pdb, "b": comparison.b.pdb,
               "error": comparison.error,
               "sequence_identity": comparison.alignment_identity,
               "tm_pairing": comparison.tm_pairing,
               "gating_mode_overlap": modes.best_overlap if modes else None,
               "shuffled_control": modes.shuffled_control if modes else None,
               "dome": ([r.summary() for r in comparison.dome["coverage_matched"]]
                        if comparison.dome else None)}, True)
        return 0 if comparison.ok else 1
    print(comparison.summary())
    if comparison.ok:
        print("\nThe overlap is the column to read. Sequence identity between "
              "these two may be inside the twilight zone, where a percentage "
              "says almost nothing — the mode overlap and its shuffled "
              "control do not depend on it.")
    return 0 if comparison.ok else 1


def add_homology_parser(sub, common) -> None:
    """Register the command. Called from ``cli.build_parser``."""
    p = sub.add_parser(
        "homology", parents=[common],
        help="the PIEZO family: members, sequence relationships, curated "
             "sites across all ten, and structural comparison")
    p.add_argument("structure", nargs="?", default=None,
                   help="a PDB id — compare that entry against every "
                        "catalogued homologue. Omit for the family matrix.")
    p.add_argument("--species", choices=["human", "mouse"], default=None)
    p.add_argument("--family", action="store_true",
                   help="just the ten members and their annotation (instant)")
    p.add_argument("--sites", action="store_true",
                   help="curated functional residues across the family")
    p.add_argument("--structural", nargs=2, metavar=("A", "B"),
                   help="compare two catalogued entries by dome, helix "
                        "correspondence and gating-mode overlap")
    p.add_argument("--n-modes", type=int, default=None,
                   help="elastic-network modes for --structural")
    p.set_defaults(func=cmd_homology)
