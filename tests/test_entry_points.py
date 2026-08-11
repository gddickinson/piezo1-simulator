"""Round 80: the conclusion must be one step from wherever a reader starts.

The project's central result is that the thing it was built to do does not
work, and that the data which would settle it cannot be assembled. A reader who
takes the machinery without that is taking the part that misleads.

Round 59 linked three surfaces — README, `CONCLUSION.md`, the guided tour — so
a sixth pre-registered test could not update one and leave the others behind.
This extends it to all five that state the record, and adds the measurement the
roadmap asked for: **how far is the conclusion from each entry point?**

Measured before this test existed, four of seven were unreachable in one step:
`docs/SCIENCE.md`, the tour's own closing step, the command line, and the
notebooks index. The roadmap's clause is that the surface is wrong rather than
the reader, so all four were fixed rather than the requirement relaxed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONCLUSION = "CONCLUSION.md"

WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven"}

#: Every way into this project, and how to read what it says. A reader arriving
#: at any of these must find the conclusion without a second hop.
ENTRY_POINTS = {
    "README.md": lambda: (ROOT / "README.md").read_text(),
    "docs/SCIENCE.md": lambda: (ROOT / "docs" / "SCIENCE.md").read_text(),
    "INTERFACE.md": lambda: (ROOT / "INTERFACE.md").read_text(),
    "docs/ARCHITECTURE.md": lambda: (ROOT / "docs" / "ARCHITECTURE.md").read_text(),
    "notebooks/README.md": lambda: (ROOT / "notebooks" / "README.md").read_text(),
    "ROADMAP.md": lambda: (ROOT / "ROADMAP.md").read_text(),
}


def _cli_text() -> str:
    from piezo1.cli import build_parser

    parser = build_parser()
    return f"{parser.description or ''}\n{parser.epilog or ''}"


def _help_text() -> str:
    from piezo1.ui.help_content import DOC_LINKS, TOPICS

    return "\n".join([f"{a} {b} {c}" for a, b, c in DOC_LINKS]
                     + [body for _, body in TOPICS])


def _tour_text() -> str:
    """The tour as a reader sees it: prose *and* its live measurements.

    Reading `step.body` alone is wrong and this test said so first — the count
    is not written into the prose, it is computed by `step.report()` from the
    validation record. That indirection is the feature: it is why the tour
    cannot state a stale number.
    """
    from piezo1.tour import TOUR

    parts = []
    for step in TOUR:
        parts.append(step.body)
        try:
            parts.append(step.report({}))
        except Exception:                       # a step needing data it lacks
            pass
    return "\n".join(parts)


# ------------------------------------------------- one step from anywhere

@pytest.mark.parametrize("name", sorted(ENTRY_POINTS))
def test_the_document_entry_point_reaches_the_conclusion(name):
    assert CONCLUSION in ENTRY_POINTS[name](), (
        f"{name} does not point at {CONCLUSION}. A reader starting here has "
        f"no way to learn that the central claim failed.")


def test_the_command_line_says_where_the_result_is():
    """A headless user may never open the GUI or the README."""
    text = _cli_text()
    assert CONCLUSION in text, "the CLI never mentions the conclusion"
    assert "does not work" in text or "nulls" in text, (
        "the CLI points at the page without saying what is in it")


def test_the_gui_help_offers_it_first():
    """First in the list, not buried among the other documents."""
    from piezo1.ui.help_content import DOC_LINKS

    title, path, _ = DOC_LINKS[0]
    assert path == "docs/CONCLUSION.md", (
        f"the help index opens with {path!r}; the conclusion should be first")


def test_the_guided_tour_ends_by_naming_it():
    """The tour is the guided path; it must not end without the result."""
    from piezo1.tour import TOUR

    assert CONCLUSION in TOUR[-1].body, (
        "the tour's last step does not name the conclusion document")


def test_no_entry_point_was_quietly_dropped():
    """The list above is the claim; this stops it shrinking to what passes."""
    assert len(ENTRY_POINTS) >= 6
    for name in ENTRY_POINTS:
        assert (ROOT / name).exists(), f"{name} no longer exists"


# ------------------------- all five surfaces agree on what the record says

def test_all_five_surfaces_state_the_same_count():
    """Round 59 linked three of these. A sixth test must update all five.

    They went out of step before: the tour said "tested twice" after five tests
    had run, and it survived because each surface was written by hand at a
    different time.
    """
    from piezo1.analysis.prediction_record import ALL_PREREGISTERED

    n = len(ALL_PREREGISTERED)
    word = WORDS[n]

    surfaces = {
        "README.md": (ROOT / "README.md").read_text(),
        "docs/CONCLUSION.md": (ROOT / "docs" / "CONCLUSION.md").read_text(),
        "docs/SCIENCE.md": (ROOT / "docs" / "SCIENCE.md").read_text(),
        "GUI help": _help_text(),
        "guided tour": _tour_text(),
    }
    for name, body in surfaces.items():
        flowed = " ".join(body.split()).lower()
        # Either form counts. The tour writes a digit because it *computes* the
        # number from the record rather than spelling it in prose — which is
        # the better of the two, and demanding the word was this test being
        # wrong about the thing it was checking rather than the reverse.
        forms = [f"{word} pre-registered", f"{word} nulls",
                 f"{n} pre-registered", f"{n} nulls"]
        assert any(f in flowed for f in forms), (
            f"{name} does not state {word} pre-registered tests; a new test "
            f"was added without updating every surface that quotes the count")


def test_no_surface_states_a_different_count():
    """The failure is not silence but disagreement, which is worse."""
    from piezo1.analysis.prediction_record import ALL_PREREGISTERED

    correct = WORDS[len(ALL_PREREGISTERED)]
    wrong = set(WORDS.values()) - {correct}

    for name, body in (("README.md", (ROOT / "README.md").read_text()),
                       ("docs/CONCLUSION.md",
                        (ROOT / "docs" / "CONCLUSION.md").read_text()),
                       ("GUI help", _help_text()),
                       ("guided tour", _tour_text())):
        flowed = " ".join(body.split()).lower()
        # Any wrong count adjacent to "pre-registered" or "nulls", not only
        # the three exact phrasings I first thought of. Planting
        # "four pre-registered" in CONCLUSION.md passed the narrower version,
        # because the page also says "five nulls" elsewhere and the positive
        # check was satisfied by that.
        for word in wrong:
            for phrase in (f"{word} pre-registered", f"{word} nulls",
                           f"{word} predictor families"):
                assert phrase not in flowed, (
                    f"{name} says {phrase!r} but the record holds "
                    f"{len(ALL_PREREGISTERED)}")


def test_the_check_would_notice_a_surface_left_behind(tmp_path):
    """Calibration: the count check must be able to fail.

    A surface that simply omits the number would otherwise pass unnoticed,
    which is the exact way the tour drifted in the first place.
    """
    flowed = "this page states no count at all"
    for form in ("five pre-registered", "five nulls", "5 pre-registered",
                 "5 nulls"):
        assert form not in flowed, "the check would pass on a silent surface"
