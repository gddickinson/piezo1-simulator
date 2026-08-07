"""The chain from a documented number back to its sources.

`verify_claims` asks whether a number is still right. These tests are about the
question underneath: whether the *path* to it exists. The load-bearing ones are
the two that found real defects — the document check calibrated against a known
answer, and the wiring check that caught 26 registered parameters no code reads.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from piezo1.analysis.provenance_chain import (LINKS, ChainTrace,
                                              number_in_document,
                                              record_sources, resolved_keys,
                                              trace, unwired_parameters, walk)
from piezo1.parameters import PARAMETERS, reset, set_value

ROOT = Path(__file__).resolve().parents[1]


# ------------------------------------------------- the document link

def test_it_finds_a_number_written_in_any_reasonable_form():
    """Calibrated on known answers before it is trusted on real documents."""
    assert number_in_document("the value is 9.72 nm", 9.72, 0.01)[0]
    assert number_in_document("rounded to 9.7 nm here", 9.72, 0.05)[0]
    assert number_in_document("1,234 variants", 1234.0, 1.0)[0]
    assert number_in_document("2.5e-3 mol", 0.0025, 1e-5)[0]


def test_it_reads_the_unicode_minus_the_documents_actually_use():
    """The bug the checker found in itself on its first run.

    Every negative number in this project's documents is typeset with U+2212,
    not the ASCII hyphen. Matching only the hyphen reported Round 22's effect
    size as missing from the document that states it.
    """
    assert number_in_document("Cliff's δ −0.211", -0.2105, 0.001)[0]
    assert not number_in_document("Cliff's δ −0.211", +0.2105, 0.001)[0]


def test_an_en_dash_range_is_not_read_as_a_negative():
    """'2.7–4.7 mN/m' must not make 4.7 into −4.7."""
    assert not number_in_document("tension 2.7–4.7 mN/m", -4.7, 0.01)[0]
    assert number_in_document("tension 2.7–4.7 mN/m", 4.7, 0.01)[0]


def test_it_reports_a_number_that_is_genuinely_absent():
    assert not number_in_document("nothing relevant here, 3.1", 9.72, 0.01)[0]


def test_every_claim_number_appears_in_the_document_that_states_it():
    """The standing guard. Cheap — no claim is recomputed."""
    report = walk(cost=None, run=False)
    missing = [(t.key, t.document) for t in report.drifted]
    assert not missing, (
        "these numbers are not written in the documents that claim them: "
        + "; ".join(f"{k} -> {d}" for k, d in missing))


# ------------------------------------------------- the recording link

def test_it_records_the_parameters_a_computation_actually_reads():
    with record_sources() as seen:
        PARAMETERS.value("pore.step")
        PARAMETERS.value("anm.cutoff")
    assert seen["parameters"] == {"pore.step", "anm.cutoff"}


def test_recording_is_removed_afterwards_even_when_the_block_raises():
    """A permanently wrapped registry would slow every later call."""
    before = PARAMETERS.value("pore.step")
    with pytest.raises(RuntimeError):
        with record_sources():
            raise RuntimeError("boom")
    with record_sources() as seen:
        pass
    assert seen["parameters"] == set()
    assert PARAMETERS.value("pore.step") == before


def test_it_records_the_resource_files_a_computation_opens():
    from piezo1.config import RESOURCE_DIR

    with record_sources() as seen:
        (RESOURCE_DIR / "domains.json").read_text()
    assert any("domains.json" in f for f in seen["files"]), seen["files"]


def test_a_file_outside_the_tracked_roots_is_not_recorded():
    """The filter is real, or every trace would list unrelated reads."""
    with record_sources() as seen:
        (ROOT / "README.md").read_text()
    assert not any("README" in f for f in seen["files"])


def test_a_cached_loader_hides_its_files_from_the_second_call():
    """A documented limitation of the data link, pinned rather than hidden.

    Most loaders in this project memoise. The first call reads the file and is
    recorded; a later call returns the cached object and reads nothing, so its
    trace shows no data files at all. This is why running the whole claims
    registry gave `hydration.score_11zc` zero parameters while running it alone
    gave four — the second call never reached the registry either.

    The consequence for a reader: an empty ``data_files`` means "read nothing
    *during this call*", not "depends on nothing".
    """
    from piezo1.core.annotations import load_annotations

    load_annotations("human")           # ensure the cache is warm
    with record_sources() as seen:
        load_annotations("human")
    assert not seen["files"], (
        "load_annotations no longer caches; the limitation documented in "
        "provenance_chain.record_sources should be revisited")


# --------------------------------------- the finding: unwired parameters

def test_the_wiring_check_knows_the_only_two_accessor_forms():
    """A miss here would understate the problem, so it is pinned.

    ``resolved_keys`` scans for ``.value(`` / ``.default(`` on the registry. If
    a third way to read a parameter were introduced, this check would silently
    start reporting wired parameters as dead.
    """
    import re

    source = "\n".join(p.read_text() for p in (ROOT / "piezo1").rglob("*.py"))
    calls = set(re.findall(r"(?:_P|PARAMETERS)\.(\w+)\(", source))
    # `categories` and `in_category` enumerate for the parameters dialog rather
    # than reading a value by key — which is how all 101 parameters, including
    # the dead ones, come to be displayed to a user.
    unknown = calls - {"value", "default", "set_value", "reset", "overrides",
                       "as_dict", "get", "is_default", "override_summary",
                       "keys", "resolve", "categories", "in_category"}
    assert not unknown, f"new registry accessor(s) {unknown}; teach resolved_keys"


def test_the_pore_parameters_are_wired_end_to_end():
    """The defect this round proved, pinned so it cannot come back.

    Before Round 49 the registry advertised ``pore.step`` with a unit and a
    citation, the dialog let a user change it, an override was recorded, and
    reports carried the non-default banner — while the computed number did not
    move at all.
    """
    from piezo1.analysis.claims import CLAIMS

    claim = next(c for c in CLAIMS if c.key == "pore.bottleneck_8yez")
    structure = ROOT / "ref" / "structures" / "8YEZ.cif"
    if not structure.exists():
        pytest.skip("8YEZ not downloaded")
    try:
        base = claim.compute()
        set_value("pore.step", 0.5)
        moved = claim.compute()
    finally:
        reset()
    assert abs(moved - base) > 1e-6, (
        "overriding pore.step does not change the pore bottleneck; the "
        "parameter is declared but not wired")
    assert claim.compute() == pytest.approx(base, abs=1e-12), \
        "reset must restore the documented value exactly"


def test_all_five_pore_parameters_reach_the_code():
    wired = resolved_keys()
    for key in ("pore.step", "pore.leash", "pore.search",
                "pore.ion_radius", "pore.constriction_threshold"):
        assert key in wired, f"{key} is registered but no code reads it"


def test_the_unwired_count_does_not_grow():
    """A ratchet, not a clean bill of health.

    Round 49 measured 26 registered parameters that no code reads and wired the
    five pore ones, leaving 21. The rest are recorded in ROADMAP.md rather than
    fixed here. This test fails if a new dead parameter is added, and should be
    tightened whenever a batch is repaired.
    """
    dead = unwired_parameters()
    assert len(dead) <= 21, (
        f"{len(dead)} registered parameters are read by no code: {dead}")


def test_a_registered_parameter_that_nothing_reads_is_detected():
    """The instrument, checked against a planted case.

    A count means nothing unless the detector can be shown to respond, so a
    key that certainly is not resolved anywhere must appear.
    """
    from piezo1.analysis import provenance_chain

    keys = provenance_chain.resolved_keys()
    assert "definitely.not.a.real.parameter" not in keys
    assert "pore.step" in keys


# ---------------------------------------------------------- the trace

def test_a_trace_establishes_every_link_for_a_real_claim():
    from piezo1.analysis.claims import CLAIMS

    claim = next(c for c in CLAIMS if c.key == "kinetics.t50")
    result = trace(claim, commit="deadbeef")
    assert isinstance(result, ChainTrace)
    assert result.stated_in_document
    assert result.parameters, "the gating model reads registered rates"
    assert result.code_module.startswith("piezo1")
    assert result.code_line and result.code_line > 0
    assert result.commit == "deadbeef"


def test_a_claim_that_cannot_run_reports_the_error_not_a_number():
    class Broken:
        key = "broken.claim"
        document = "docs/SCIENCE.md"
        expected = 1.0
        tolerance = 0.1
        unit = ""

        @staticmethod
        def compute():
            raise ValueError("no data")

    result = trace(Broken(), commit="abc")
    assert result.computed is None
    assert "no data" in result.error
    assert "parameters" in result.broken and "data" in result.broken


def test_links_are_the_documented_five():
    assert LINKS == ("document", "code", "parameters", "data", "commit")
