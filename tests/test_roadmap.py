"""The roadmap split, checked against the file it was split from.

Round 75's own validation clause: *no completed item may be lost — every
measured result recorded in a checkbox must survive the split, and a test
should count them before and after.*

"Before" is not a number typed in from memory. It is the actual pre-split file,
read out of git at the commit that split it, and compared item by item. The
frozen counts below are a second, weaker check that still works where git
history is unavailable — a shallow clone, or a source tarball.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ROADMAP = ROOT / "ROADMAP.md"
ARCHIVE = ROOT / "docs" / "ROADMAP_COMPLETED.md"

#: The commit whose ROADMAP.md was split. Its content is the ground truth for
#: "nothing was lost", and it does not move as later rounds are added.
PRE_SPLIT_COMMIT = "4c1c61c"

#: Measured after the split, with the duplicated Round 68 merged and Round 75
#: itself ticked. The pre-split file had 358 completed and 14 open.
COMPLETED_AT_SPLIT = 362
OPEN_AT_SPLIT = 10


def items(text: str, mark: str) -> list[str]:
    return [line for line in text.splitlines() if line.startswith(mark)]


def pre_split_roadmap() -> str | None:
    """The file as it was before the split, or None if history is unavailable."""
    try:
        done = subprocess.run(
            ["git", "show", f"{PRE_SPLIT_COMMIT}:ROADMAP.md"],
            cwd=ROOT, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):       # pragma: no cover
        return None
    return done.stdout if done.returncode == 0 else None


# ------------------------------------------------------- the split preserved

def test_both_files_exist_and_neither_is_a_stub():
    assert ROADMAP.exists() and ARCHIVE.exists()
    assert len(ARCHIVE.read_text().splitlines()) > 2000, \
        "the completed record is too short to hold 75 rounds"


def test_every_completed_item_survived_the_split():
    """The real check: item by item against the pre-split file in git."""
    original = pre_split_roadmap()
    if original is None:
        pytest.skip(f"no git history for {PRE_SPLIT_COMMIT}")

    before = items(original, "- [x]")
    after = items(ROADMAP.read_text(), "- [x]") + items(ARCHIVE.read_text(), "- [x]")

    missing = sorted(set(before) - set(after))
    assert not missing, (
        f"{len(missing)} completed item(s) lost in the split, first: "
        f"{missing[0][:120]}")
    assert len(after) >= len(before), "an item was dropped as a duplicate"


def test_no_completed_item_is_in_both_files():
    """Duplication would inflate the count and hide a loss underneath it."""
    in_roadmap = set(items(ROADMAP.read_text(), "- [x]"))
    in_archive = set(items(ARCHIVE.read_text(), "- [x]"))
    assert not (in_roadmap & in_archive), sorted(in_roadmap & in_archive)[:3]


def test_every_open_item_survived_and_lives_in_the_roadmap():
    original = pre_split_roadmap()
    if original is None:
        pytest.skip(f"no git history for {PRE_SPLIT_COMMIT}")

    # Round 68 was open in the pre-split file and is closed here, so the open
    # set may shrink by exactly that round — but nothing may vanish untraced.
    before = set(items(original, "- [ ]"))
    after = set(items(ROADMAP.read_text(), "- [ ]"))
    closed = before - after
    assert len(closed) == 4, (
        f"{len(closed)} open items disappeared; only four were closed — the "
        f"duplicated Round 68's two (Round 63 had done them) and Round 75's "
        f"own two, this round: {sorted(closed)[:5]}")
    assert not items(ARCHIVE.read_text(), "- [ ]"), \
        "an unfinished item is filed under completed work"


def test_the_duplicated_round_68_is_merged_and_what_happened_is_recorded():
    """The defect the split surfaced, kept rather than tidied away.

    The file carried two adjacent `Round 68` headings — one ticked as
    superseded by Round 63, one still open asking the same question, because
    Round 67 recorded the supersession by adding a heading instead of ticking
    the original. Silently deleting the open copy would have been the tidy
    option and would have destroyed the only evidence that the roadmap had
    drifted out of step with the work.
    """
    text = ARCHIVE.read_text()
    assert text.count("### Round 68") == 1, "the duplicate heading is back"
    entry = text[text.index("### Round 68"):text.index("### Round 69")]
    assert "Round 63" in entry, "the entry must name the round that did it"
    assert "Round 75" in entry, "and the round that found the duplication"
    assert "- [ ]" not in entry
    assert entry.count("- [x]") == 4, (
        "the supersession note, the original item, its validate clause, and "
        "the record of the duplication")


# ------------------------------- the counts, for a clone without git history

def test_the_completed_count_cannot_fall():
    """A ratchet. Work is only ever finished, so this number only ever grows."""
    total = (len(items(ROADMAP.read_text(), "- [x]"))
             + len(items(ARCHIVE.read_text(), "- [x]")))
    assert total >= COMPLETED_AT_SPLIT, (
        f"{total} completed items, down from {COMPLETED_AT_SPLIT} at the split")


def test_the_open_items_are_all_in_the_live_roadmap():
    open_items = items(ROADMAP.read_text(), "- [ ]")
    assert open_items, "the roadmap lists nothing to do; is Block Q finished?"
    assert len(open_items) <= OPEN_AT_SPLIT + 40, \
        "the roadmap has grown a great deal; check a block was not duplicated"


# ---------------------------------------------------- the files stay useful

def test_the_roadmap_is_short_enough_to_read():
    """The point of the split. 2,702 lines is not a list of what is left."""
    assert len(ROADMAP.read_text().splitlines()) < 400


def test_each_file_points_at_the_other():
    assert "ROADMAP_COMPLETED.md" in ROADMAP.read_text()
    assert "ROADMAP.md" in ARCHIVE.read_text()


def test_the_roadmap_keeps_the_standing_sections():
    """The destination and the per-round checklist are live, not history."""
    text = ROADMAP.read_text()
    for heading in ("The destination", "Standing per-round checklist",
                    "Deliberately not doing"):
        assert heading in text, f"{heading} was filed as completed work"


def test_the_archive_keeps_the_review_sections():
    """The five-round reviews decided what came next; they are the record of
    how the plan changed, and belong with the rounds they reviewed."""
    text = ARCHIVE.read_text()
    assert text.count("## Review after Rounds") >= 4


def test_the_documents_index_offers_both():
    from piezo1.ui.help_content import DOC_LINKS

    paths = {path for _, path, _ in DOC_LINKS}
    assert "ROADMAP.md" in paths
    assert "docs/ROADMAP_COMPLETED.md" in paths


def test_every_shipped_document_link_resolves():
    """Not specific to the split, but this is where a renamed file shows up."""
    from piezo1.ui.help_content import DOC_LINKS

    missing = [path for _, path, _ in DOC_LINKS if not (ROOT / path).exists()]
    assert not missing, f"help links to documents that do not exist: {missing}"
