#!/usr/bin/env python
"""Reproduce the entire working state of the project from a fresh clone.

Aim A5 in `CLAUDE.md`: a clone plus an environment plus this script should
rebuild everything. The steps are ordered by dependency, each is skippable, and
each reports what it did rather than what it intended to do.

The step that matters most is ``--verify``. It recomputes every headline number
the documentation asserts and reports any that has drifted. Documentation is
prose and prose does not fail a test suite; without this, a solver rewrite can
leave `docs/SCIENCE.md` confidently stating a number the code stopped
producing. See :mod:`piezo1.analysis.claims`.

Usage::

    python scripts/reproduce.py                 # everything
    python scripts/reproduce.py --verify        # only check the documented numbers
    python scripts/reproduce.py --quick         # skip the slow analyses
    python scripts/reproduce.py --skip fetch    # data already downloaded
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def run(label: str, command: list[str], optional: bool = False) -> bool:
    """Run one step, echoing its command so a failure can be reproduced by hand."""
    print(f"\n=== {label}")
    print(f"    $ {' '.join(command)}", flush=True)
    start = time.time()
    result = subprocess.run(command, cwd=ROOT)
    seconds = time.time() - start
    if result.returncode == 0:
        print(f"    ok ({seconds:.0f} s)")
        return True
    message = f"    FAILED with code {result.returncode} after {seconds:.0f} s"
    print(message + ("  (optional, continuing)" if optional else ""))
    return optional


STEPS = [
    ("fetch", "Download structures, sequences, ligands and the CHAP grid",
     [sys.executable, "-m", "piezo1.io.fetch"], False),
    ("parameters", "Rebuild the parameter registry (checks every citation)",
     [sys.executable, "scripts/build_parameters.py"], False),
    ("audit", "Check no number in the science modules is unregistered",
     [sys.executable, "-m", "piezo1.parameter_audit"], False),
    ("resources", "Rebuild the curated annotation resources",
     [sys.executable, "scripts/build_uniprot_annotations.py"], True),
    ("domains", "Rebuild domain definitions",
     [sys.executable, "scripts/build_domains.py"], True),
    ("residues", "Rebuild functional-residue groups",
     [sys.executable, "scripts/build_functional_residues.py"], True),
    ("variants", "Rebuild the curated variant table",
     [sys.executable, "scripts/build_variants.py"], True),
    ("registry", "Rebuild the structure registry",
     [sys.executable, "scripts/build_structure_registry.py"], True),
    ("references", "Resolve the bibliography",
     [sys.executable, "scripts/build_references.py"], True),
    ("tests", "Run the test suite",
     [sys.executable, "-m", "pytest"], False),
    ("validation", "Re-run the Round 7 blind test",
     [sys.executable, "scripts/run_validation.py"], True),
    ("validation22", "Re-run the Round 22 exploratory test",
     [sys.executable, "scripts/run_validation_round22.py"], True),
    ("figures", "Regenerate the documentation figures",
     [sys.executable, "scripts/make_figures.py"], True),
    ("screenshots", "Drive the GUI and capture screenshots",
     [sys.executable, "scripts/screenshot_app.py", "--structure", "8YEZ",
      "--analysis"], True),
]


def verify(max_cost: str) -> int:
    """Recompute every documented number and report drift."""
    from piezo1.analysis.claims import claims_by_cost, verify_claims

    claims = claims_by_cost(max_cost)
    print(f"\n=== Verifying {len(claims)} documented numbers "
          f"(cost <= {max_cost})\n")
    results = verify_claims(claims)

    drifted = [r for r in results if not r.ok and not r.error]
    skipped = [r for r in results if r.error]
    print(f"\n  {len(results) - len(drifted) - len(skipped)} reproduced, "
          f"{len(drifted)} drifted, {len(skipped)} skipped")

    if drifted:
        print("\n  DOCUMENTATION DRIFT — these documents now state numbers the "
              "code does not produce:")
        for result in drifted:
            claim = result.claim
            print(f"    {claim.document}: {claim.description}")
            print(f"      documented {claim.expected} {claim.unit}, "
                  f"computed {result.value:.4f}")
            if claim.frozen:
                print("      *** FROZEN RESULT *** — this is a recorded "
                      "finding. Do not edit the document to match; work out "
                      "why the computation changed.")
            elif claim.published:
                print(f"      published reference: {claim.published}")
    if skipped:
        print("\n  Skipped (missing data — run the fetch step):")
        for result in skipped:
            print(f"    {result.claim.key}: {result.error}")
    return 1 if drifted else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--verify", action="store_true",
                        help="only recompute and check the documented numbers")
    parser.add_argument("--quick", action="store_true",
                        help="skip slow steps and slow claims")
    parser.add_argument("--skip", nargs="*", default=[], metavar="STEP",
                        help=f"steps to skip: {', '.join(s[0] for s in STEPS)}")
    parser.add_argument("--only", nargs="*", default=[], metavar="STEP",
                        help="run only these steps")
    args = parser.parse_args()

    max_cost = "medium" if args.quick else "slow"
    if args.verify:
        return verify(max_cost)

    print("PIEZO1 simulator — full reproduction")
    print(f"root: {ROOT}")

    failures = []
    for key, label, command, optional in STEPS:
        if key in args.skip or (args.only and key not in args.only):
            print(f"\n=== {label}\n    skipped")
            continue
        if args.quick and key in ("figures", "screenshots", "references"):
            print(f"\n=== {label}\n    skipped (--quick)")
            continue
        if not run(label, command, optional):
            failures.append(key)

    code = verify(max_cost)

    print("\n" + "=" * 70)
    if failures:
        print(f"FAILED steps: {', '.join(failures)}")
    if code:
        print("Documented numbers have drifted — see above.")
    if not failures and not code:
        print("Reproduced cleanly: every step ran and every documented number "
              "still comes out.")
    return 1 if (failures or code) else 0


if __name__ == "__main__":
    raise SystemExit(main())
