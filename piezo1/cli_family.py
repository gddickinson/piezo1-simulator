"""The ``family`` command, split from ``cli.py`` at the length limit.

Split along the same seam ``cli_homology.py`` uses, and for a sharper version of
the same reason. Every other command in ``cli.py`` takes a structure and answers
a question about one deposited entry; ``homology`` answers questions about the
*family*. This one answers a question about **somebody else's project**: what
was imported from the ``piezo_genes`` census, what each statement rests on, and
which module here explores it.

It needs no structure at all, and accepts one anyway — optionally and unused —
because every analysis in the shared registry has to be reachable by the same
command shape and a test enforces it. ``liu2025`` does the same.

Printing is the point rather than a formality. The four structure-taking family
commands beside it all *measure* something; this one is the only place a reader
is told which numbers were measured here and which were not, so it prints the
source, the commit, how many numbers were re-read at build time, and for every
statement its caveat.
"""

from __future__ import annotations

__all__ = ["cmd_family", "add_family_parser"]


def add_family_parser(sub, common) -> None:
    """Register the ``family`` subcommand on an existing subparser set."""
    p = sub.add_parser("family", parents=[common],
                       help="the imported PIEZO-family census: every finding, "
                            "its source, and what is done with it here")
    p.add_argument("structure", nargs="?", default=None,
                   help="accepted for symmetry with the other commands and "
                        "used by none of them")
    p.add_argument("--kind", choices=["sequence", "clinical", "structure",
                                      "census"],
                   help="only findings of one kind")
    p.add_argument("--key", help="one finding by key")
    p.set_defaults(func=cmd_family)


def cmd_family(args) -> int:
    """The imported census: every statement, its source, and what is done here.

    The right first stop before any of the four structure-taking family
    commands: everything it prints is a property of the imported resource —
    which findings came across, how many of their numbers were re-read from the
    source on the last build, and which module in this project explores each.
    """
    from .core.family import load_family_findings

    findings = load_family_findings()
    if args.key:
        finding = findings.by_key(args.key)
        if finding is None:
            print(f"no finding {args.key!r}; have: {', '.join(findings.keys)}")
            return 1
        selected = (finding,)
    elif args.kind:
        selected = findings.by_kind(args.kind)
    else:
        selected = findings.findings

    payload = {
        "source": findings.source,
        "verified": findings.provenance["verified"],
        "census": findings.census,
        "findings": [
            {"key": f.key, "session": f.session, "kind": f.kind,
             "title": f.title, "statement": f.statement,
             "numbers": f.numbers, "source_file": f.source,
             "explored_here": f.here, "caveat": f.caveat}
            for f in selected],
    }
    if args.json:
        from .cli import _emit
        _emit(payload, True)
        return 0

    print(f"Imported from {findings.source}")
    print(f"  {findings.provenance['verified']}\n")
    for f in selected:
        print(f"[{f.session}] {f.title}   ({f.kind})")
        print(f"  {f.statement}")
        if f.numbers:
            print("  numbers: " + ", ".join(f"{k}={v}" for k, v in f.numbers.items()))
        print(f"  here:    {f.here}")
        print(f"  caveat:  {f.caveat}\n")
    return 0
