"""Walk a documented number back to everything it came from.

:mod:`claims` answers "is this number still what the code produces?".
This answers the question underneath it: **can the path to the number be
reconstructed at all?** A claim that recomputes correctly is still unprovenanced
if nobody can say which structure file it read, which registered parameters it
consumed, or which commit produced the document it is written in.

The five links, and how each is established:

===================  =========================================================
link                 how it is determined
===================  =========================================================
document             the file is read and searched for a number matching the
                     claim within its own tolerance — **not** assumed
code                 the compute callable's module, qualified name and source
                     line, read off the function object
parameters           **recorded while the claim runs**, by wrapping the single
                     registry read path — not declared, and not guessed from
                     the source text
data                 structure files and resource files opened during the run,
                     with a content hash of each
commit               ``git rev-parse HEAD`` plus whether the tree is dirty
===================  =========================================================

The parameter and data links are *measured*, which is the point. A claim that
silently depends on a parameter nobody listed will show that dependency here,
and a claim that reads no data file at all is one whose number cannot be traced
to any input.

**What a broken link means.** ``ChainTrace.broken`` lists the links that could
not be established. A break is not necessarily a defect — a claim computed from
pure constants legitimately reads no structure file — so each break carries the
reason, and :func:`walk` separates the ones that indicate real drift (a number
absent from the document that claims to state it) from the ones that are simply
facts about the claim.
"""

from __future__ import annotations

import hashlib
import inspect
import re
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass, field

__all__ = ["ChainTrace", "ChainReport", "trace", "walk", "record_sources",
           "git_state", "number_in_document", "resolved_keys",
           "unwired_parameters", "LINKS"]

#: The five links a claim's provenance chain is made of, in the order a reader
#: follows them: from the written number back to the commit that produced it.
LINKS = ("document", "code", "parameters", "data", "commit")

#: Numbers written in prose, including exponents and thousands separators.
#:
#: The leading sign class includes U+2212 MINUS SIGN, because the documents in
#: this project are typeset and write negative values as "−0.211" rather than
#: "-0.211". Matching only the ASCII hyphen reported the Round 22 effect size
#: as missing from its own document when it was written there all along — the
#: checker's first run found a bug in the checker.
#:
#: U+2013 EN DASH is deliberately **not** treated as a minus: the science
#: documents use it for ranges ("2.7–4.7 mN/m"), and accepting it would parse
#: the upper bound of every range as a negative number.
_NUMBER = re.compile(r"[-+\u2212]?\d[\d,]*\.?\d*(?:[eE][-+\u2212]?\d+)?")


@dataclass
class ChainTrace:
    """One claim, and everything reachable from it."""

    key: str
    document: str = ""
    document_exists: bool = False
    stated_in_document: bool = False
    matched_text: str = ""
    code_module: str = ""
    code_qualname: str = ""
    code_line: int | None = None
    parameters: list = field(default_factory=list)
    data_files: list = field(default_factory=list)
    commit: str = ""
    dirty: bool = False
    computed: float | None = None
    error: str = ""
    broken: dict = field(default_factory=dict)

    @property
    def complete(self) -> bool:
        """Every link established. Reported, never assumed."""
        return not self.broken

    def summary(self) -> str:
        state = "complete" if self.complete else f"broken: {', '.join(self.broken)}"
        return (f"{self.key}: {len(self.parameters)} parameters, "
                f"{len(self.data_files)} data files, {state}")


@dataclass
class ChainReport:
    """Every claim's chain, and which breaks are real problems."""

    traces: list = field(default_factory=list)
    commit: str = ""
    dirty: bool = False

    def get(self, key: str) -> ChainTrace | None:
        return next((t for t in self.traces if t.key == key), None)

    @property
    def complete(self) -> list:
        return [t for t in self.traces if t.complete]

    @property
    def drifted(self) -> list:
        """Claims whose number is **not** in the document that states it.

        This is the break that means something is wrong, as opposed to a claim
        that legitimately reads no structure file.
        """
        return [t for t in self.traces
                if t.document_exists and not t.stated_in_document]

    @property
    def unparameterised(self) -> list:
        """Claims that consumed no registered parameter while running."""
        return [t for t in self.traces if not t.parameters and not t.error]

    def explain(self) -> dict:
        """Why the incomplete chains are incomplete, by category.

        The bare fraction reads as failure and is not. A claim computed from a
        frozen validation record legitimately consumes no registered parameter;
        an analytic claim legitimately reads no structure. What would be a real
        break is a number missing from its own document, or code that cannot be
        located — and both of those are zero.
        """
        counts: dict = {}
        for trace in self.traces:
            for link, reason in trace.broken.items():
                key = ("no registered parameter"
                       if "no registered parameter" in reason else
                       "no structure or resource file"
                       if "no structure or resource" in reason else
                       f"{link}: {reason[:40]}")
                counts[key] = counts.get(key, 0) + 1
        return {"complete": len(self.complete), "total": len(self.traces),
                "breaks": counts,
                "document_breaks": len(self.drifted),
                "benign": all(
                    k in ("no registered parameter",
                          "no structure or resource file") for k in counts)}

    def summary(self) -> str:
        return (f"{len(self.complete)}/{len(self.traces)} chains complete; "
                f"{len(self.drifted)} numbers missing from their own document; "
                f"commit {self.commit[:8] or 'unknown'}"
                f"{' (dirty)' if self.dirty else ''}")


# --------------------------------------------------------------------------
# The individual links
# --------------------------------------------------------------------------

def git_state(root=None) -> tuple[str, bool]:
    """``(commit, dirty)``. Empty commit when git cannot answer."""
    from ..config import PROJECT_ROOT

    root = root or PROJECT_ROOT
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root,
                                capture_output=True, text=True, timeout=30)
        status = subprocess.run(["git", "status", "--porcelain"], cwd=root,
                                capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return "", False
    if commit.returncode != 0:
        return "", False
    return commit.stdout.strip(), bool(status.stdout.strip())


def number_in_document(text: str, expected: float,
                       tolerance: float) -> tuple[bool, str]:
    """Is a number matching ``expected`` written anywhere in ``text``?

    Every number in the prose is extracted and compared, rather than searching
    for one formatting of the value — documents round, add units and use
    thousands separators, and a string search would report drift that is really
    a formatting difference.

    A coincidental match is possible in a document full of numbers, so this can
    pass when it should not. It can also fail when it should pass, if a document
    writes a number in a form the pattern does not recognise — which happened on
    the first run, with the Unicode minus sign. Both directions are stated
    because the second was found the hard way.
    """
    window = max(abs(tolerance), abs(expected) * 1e-3, 1e-9)
    for match in _NUMBER.finditer(text):
        raw = (match.group().rstrip(".").replace(",", "")
               .replace("\u2212", "-"))
        try:
            value = float(raw)
        except ValueError:
            continue
        if abs(value - expected) <= window:
            return True, match.group()
    return False, ""


@contextmanager
def record_sources():
    """Record registry reads and file opens for the duration of the block.

    Yields ``{"parameters": set, "files": set}``, filled in as the wrapped code
    runs. Both the registry read path and the two file doors are restored
    afterwards even if the block raises.

    Covers ``ParameterRegistry.value`` (the only way a registered number is
    read), ``Structure.from_file`` (the door for coordinates) and
    ``Path.read_text`` (the door for the curated JSON resources). A claim that
    reached data some other way would show no files, which is itself reported.

    **Caching limits what this can see, and the limit is real.** Most loaders
    here memoise, so only the *first* call reads the file or the registry. Run
    the whole claims registry and ``hydration.score_11zc`` reports zero
    parameters; run it alone and it reports four — the batch run reached a
    warm cache. So an empty result means "read nothing during this call", not
    "depends on nothing", and a chain is most informative when traced with
    ``--key`` on a cold process. ``test_provenance_chain`` pins this.
    """
    import pathlib

    from ..config import DERIVED_DIR, RESOURCE_DIR, STRUCTURE_DIR
    from ..core.structure import Structure
    from ..parameters import ParameterRegistry

    seen = {"parameters": set(), "files": set()}
    tracked = [p for p in (STRUCTURE_DIR, RESOURCE_DIR, DERIVED_DIR) if p]

    real_value = ParameterRegistry.value
    real_from_file = Structure.from_file.__func__
    real_read_text = pathlib.Path.read_text

    def value(self, key):
        seen["parameters"].add(key)
        return real_value(self, key)

    def from_file(cls, path, *args, **kwargs):
        seen["files"].add(str(path))
        return real_from_file(cls, path, *args, **kwargs)

    def read_text(self, *args, **kwargs):
        if any(str(self).startswith(str(root)) for root in tracked):
            seen["files"].add(str(self))
        return real_read_text(self, *args, **kwargs)

    ParameterRegistry.value = value
    Structure.from_file = classmethod(from_file)
    pathlib.Path.read_text = read_text
    try:
        yield seen
    finally:
        ParameterRegistry.value = real_value
        Structure.from_file = classmethod(real_from_file)
        pathlib.Path.read_text = real_read_text


#: How a registered number is read. These are the only two forms in the
#: package, checked rather than assumed — see ``test_provenance_chain``.
_RESOLVE = re.compile(
    r"(?:_P|PARAMETERS)\.(?:value|default)\(\s*[\"\']([\w.]+)[\"\']")


def resolved_keys(package=None) -> set:
    """Every registry key the package actually reads, by scanning for the call.

    Static rather than dynamic because a parameter read only on a branch no
    test exercises is still wired; running the code would under-report it. The
    risk in the other direction — a key named in a string but never reached —
    makes this an optimistic bound, so a key missing from here is definitely
    unwired while a key present here is only probably wired.
    """
    from pathlib import Path

    from ..config import PROJECT_ROOT

    root = Path(package) if package else PROJECT_ROOT / "piezo1"
    keys = set()
    for file in root.rglob("*.py"):
        keys.update(_RESOLVE.findall(file.read_text()))
    return keys


def unwired_parameters() -> list:
    """Registered parameters that no code reads — the chain's worst break.

    A parameter here is **actively misleading**, not merely unused: it appears
    in the parameters dialog with a unit and a citation, an override on it is
    recorded, reports carry the non-default banner because of it, and
    ``verify_claims`` refuses to run against it — while the number it claims to
    control does not move.

    That is a strictly worse failure than an unregistered literal, which is at
    least honestly invisible. The parameter audit cannot see it: the audit
    checks that a literal *is declared* to correspond to a registered
    parameter, and a declaration is not a wire.
    """
    from ..parameters import PARAMETERS

    registered = set(PARAMETERS.as_dict(only_overrides=False))
    return sorted(registered - resolved_keys())


def _digest(path: str) -> str:
    from pathlib import Path

    file = Path(path)
    if not file.exists():
        return ""
    digest = hashlib.sha256()
    with file.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()[:12]


# --------------------------------------------------------------------------
# The walk
# --------------------------------------------------------------------------

def trace(claim, commit: str = "", dirty: bool = False,
          run: bool = True) -> ChainTrace:
    """Establish every link for one claim.

    ``run=False`` skips executing the computation, which leaves the parameter
    and data links unestablished — useful for checking the document link over
    the whole registry without paying for the slow claims.
    """
    from pathlib import Path

    from ..config import PROJECT_ROOT

    out = ChainTrace(key=claim.key, document=claim.document,
                     commit=commit, dirty=dirty)

    # --- code
    function = claim.compute
    out.code_module = getattr(function, "__module__", "")
    out.code_qualname = getattr(function, "__qualname__", "")
    try:
        out.code_line = inspect.getsourcelines(function)[1]
    except (OSError, TypeError):
        out.broken["code"] = "source for the compute callable is not readable"

    # --- document
    path = PROJECT_ROOT / claim.document
    out.document_exists = path.exists()
    if not out.document_exists:
        out.broken["document"] = f"{claim.document} does not exist"
    else:
        found, text = number_in_document(path.read_text(), claim.expected,
                                         claim.tolerance)
        out.stated_in_document, out.matched_text = found, text
        if not found:
            out.broken["document"] = (
                f"{claim.expected} +-{claim.tolerance} {claim.unit} is not "
                f"written anywhere in {claim.document}")

    # --- parameters and data, measured by running it
    if not run:
        out.broken["parameters"] = "not run"
        out.broken["data"] = "not run"
    else:
        try:
            with record_sources() as seen:
                out.computed = float(claim.compute())
            out.parameters = sorted(seen["parameters"])
            files = []
            for raw in seen["files"]:
                name = str(raw)
                if name.startswith(str(PROJECT_ROOT)):
                    name = str(Path(raw).relative_to(PROJECT_ROOT))
                files.append({"path": name, "sha256": _digest(raw)})
            out.data_files = sorted(files, key=lambda d: d["path"])
        except Exception as exc:  # a claim that cannot run has no chain
            out.error = f"{type(exc).__name__}: {exc}"
            out.broken["parameters"] = "the computation did not complete"
            out.broken["data"] = "the computation did not complete"
        if run and not out.error and not out.parameters:
            out.broken["parameters"] = "consumed no registered parameter"
        if run and not out.error and not out.data_files:
            out.broken["data"] = "read no structure or resource file"

    # --- commit
    if not commit:
        out.broken["commit"] = "git could not identify the working tree"

    return out


def walk(claims=None, cost: tuple = ("fast", "medium"),
         run: bool = True) -> ChainReport:
    """Walk every claim's chain and report where it breaks.

    ``cost`` filters by the claim's own runtime tier, because establishing the
    parameter and data links means running the computation. Passing ``None``
    walks all of them.
    """
    from .claims import CLAIMS

    selected = list(claims if claims is not None else CLAIMS)
    if cost is not None:
        selected = [c for c in selected if c.cost in cost]
    commit, dirty = git_state()
    return ChainReport(
        traces=[trace(c, commit=commit, dirty=dirty, run=run)
                for c in selected],
        commit=commit, dirty=dirty)


def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Walk each documented number back to its sources.")
    parser.add_argument("--all", action="store_true",
                        help="include slow claims")
    parser.add_argument("--no-run", action="store_true",
                        help="check documents only; do not recompute")
    parser.add_argument("--key", help="one claim by key")
    args = parser.parse_args()

    from .claims import CLAIMS

    claims = [c for c in CLAIMS if c.key == args.key] if args.key else None
    report = walk(claims=claims, cost=None if (args.all or args.key)
                  else ("fast", "medium"), run=not args.no_run)

    for item in report.traces:
        mark = "ok  " if item.complete else "BREAK"
        print(f"  {mark} {item.key:32s} {item.code_module}.{item.code_qualname}"
              f":{item.code_line}")
        if item.parameters:
            print(f"        parameters: {', '.join(item.parameters)}")
        for data in item.data_files:
            print(f"        data: {data['path']}  sha256:{data['sha256']}")
        for link, reason in item.broken.items():
            print(f"        BROKEN {link}: {reason}")

    dead = unwired_parameters()
    print(f"\n{len(dead)} registered parameters are read by no code:")
    for key in dead:
        print(f"    {key}")

    print("\n" + report.summary())
    if report.drifted:
        print("\nNumbers missing from the document that states them:")
        for item in report.drifted:
            print(f"  {item.key} -> {item.document}")
    return 1 if report.drifted else 0


if __name__ == "__main__":
    raise SystemExit(_main())
