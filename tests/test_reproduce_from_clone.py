"""What an empty clone surfaces that a warm cache hides.

Round 60 ran the whole chain from a fresh clone with no `ref/` or `data/`.
Three defects appeared, none of which can be seen on a developer machine
because the files are already on disk:

1. **Eight tests failed instead of skipping.** `conftest.py` states the rule —
   skip when data is absent — and eight tests did not follow it. Three were
   written in recent rounds, by me.
2. **The Ensembl CDS download was broken.** The content type was sent as a
   `;content-type=` query parameter, which Ensembl no longer honours; the plain
   URL now answers 415. Anyone who already had the file never noticed.
3. **`feasibility.assess()` silently lowered its ceiling** from 59 to 34 when
   the open-access corpus is missing — a documented number changing with the
   cache state rather than reporting that it could not be computed.

Measured chain: clone 2 s · empty-clone suite 38 s · fetch 21 s for 43
resources · full suite 209 s · claims 6 s.
"""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


# ------------------------------------------- the download that was broken

def test_the_cds_fetch_sends_the_content_type_as_a_header():
    """A query parameter is not a header, and Ensembl stopped pretending.

    Checked on the source rather than over the network, so the suite still runs
    offline — the network behaviour is what motivated it, not what pins it.
    """
    source = (ROOT / "piezo1" / "io" / "fetch.py").read_text()

    # Comments are stripped first: the fix's own comment explains the old
    # query-parameter form, and a naive scan trips on the explanation rather
    # than the code — the same way Round 55's planted probe found itself.
    code = "\n".join(line.split("#")[0] for line in source.splitlines())
    assert ";content-type=" not in code, (
        "the content type must be an HTTP header; as a query parameter "
        "Ensembl answers 415 and the fetch fails silently on a cold cache")
    assert 'headers={"Content-Type": "text/x-fasta"}' in code


def test_the_downloader_accepts_headers_at_all():
    import inspect

    from piezo1.io.fetch import _download

    assert "headers" in inspect.signature(_download).parameters


# ------------------------------- the number that changed with the cache

def test_the_feasibility_ceiling_reports_when_it_cannot_be_computed():
    """The bug: a missing corpus made the ceiling *smaller*, not *unavailable*.

    59 is a documented number. Without the harvest the same call returned 34
    and said nothing, so a reader on a fresh clone would have drawn a stronger
    conclusion than the data supports — from the same code, silently.
    """
    from piezo1.analysis.feasibility import FeasibilityReport, assess

    empty = FeasibilityReport(meta={"harvest_available": False})
    assert not empty.harvest_available
    assert "not downloaded" in empty.summary()
    assert "fetch" in empty.summary()

    report = assess(n_simulations=200)
    if not report.harvest_available:
        pytest.skip("corpus not downloaded here; the guard is what matters")
    assert report.ceiling_n > report.meta["directional_variants"]


def test_a_report_defaults_to_claiming_availability_only_when_told():
    """A missing key must not silently mean "fine"... except where it does.

    `harvest_available` defaults True so that hand-built reports in tests are
    not all marked unavailable; `assess` always sets it explicitly. This pins
    that `assess` does so, which is what makes the default safe.
    """
    from piezo1.analysis.feasibility import assess

    assert "harvest_available" in assess(n_simulations=100).meta


# --------------------------------- tests must skip, not fail, without data

def test_every_data_dependent_test_declares_how_it_skips():
    """The eight that failed all now name the fetch command in their skip.

    A skip that does not say what to run is a dead end for whoever hits it.
    """
    offenders = []
    for name in ("test_harvest_curation.py", "test_feasibility.py",
                 "test_variant_structures.py", "test_workflow.py",
                 "test_sequence_and_resources.py"):
        text = (ROOT / "tests" / name).read_text()
        if "pytest.skip" in text and "io.fetch" not in text:
            offenders.append(name)
    assert not offenders, (
        f"these skip without saying how to fix it: {offenders}")


def test_the_conftest_rule_is_stated_where_it_can_be_found():
    text = (ROOT / "tests" / "conftest.py").read_text().lower()
    assert "skip" in text and ("download" in text or "fetch" in text)
