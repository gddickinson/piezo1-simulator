"""The README, guarded against the two ways it actually went wrong.

It **went stale**. The Status section listed as "Planned" a membrane solver,
gating kinetics, morphing, pore profiling and pocket detection — all of which
had been implemented and shipped, some of them for dozens of rounds. A reader
would have concluded the project did far less than it does.

And it **carried claims with nothing behind them**. Rewriting it, I introduced
two of my own in a single paragraph: that E756del is carried by a third of
people of African ancestry (the project's own data says the gnomAD AFR
frequency is 0.166–0.173, so a sixth), and that it protects against malaria
(`docs/SCIENCE.md` records that as contested, with the odds ratio and the fact
that the mouse work tested a different allele).

So: no "planned" claim about something that exists, every citation traceable to
the verified bibliography, and every link and image resolving. The scientific
numbers in the closing summary have their own guard in `test_conclusion.py`.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"


@pytest.fixture(scope="module")
def text() -> str:
    return README.read_text()


# ------------------------------------------------------- it must not go stale

def test_nothing_implemented_is_described_as_planned(text):
    """The specific failure: a Status section years out of date.

    Each of these is shipped, tested and reachable from the GUI. If the README
    ever calls one of them planned or future work again, this fails.
    """
    shipped = ("membrane", "kinetics", "morph", "pore", "pocket", "permeation",
               "docking")
    lowered = " ".join(text.split()).lower()
    for word in ("planned:", "not yet implemented", "still to come",
                 "coming soon"):
        assert word not in lowered, f"the README claims future work: {word!r}"

    # And the modules it describes must actually import.
    import importlib

    for module in ("piezo1.physics.membrane", "piezo1.physics.kinetics",
                   "piezo1.physics.permeation", "piezo1.structure.morph",
                   "piezo1.structure.pore", "piezo1.structure.hybrid",
                   "piezo1.structure.fusion_pose", "piezo1.analysis.pockets"):
        importlib.import_module(module)
    assert shipped        # the list above is the record of what went stale


def test_every_feature_the_readme_names_has_a_menu_or_a_command(text):
    """A README that advertises what the application cannot do is worse than
    one that is merely out of date."""
    from piezo1.analysis.report import ANALYSES

    for name in ("hybrid", "permeation", "fusion", "pore", "dome"):
        assert name in ANALYSES, f"the README describes {name}, which is absent"
        assert f"piezo1.cli {name}" in text or name in text


# ----------------------------------------------------------- links and images

def test_every_relative_link_and_image_resolves(text):
    broken = []
    for match in re.finditer(r"!?\[[^\]]*\]\(([^)]+)\)", text):
        target = match.group(1)
        if target.startswith(("http", "#", "mailto")):
            continue
        if not (ROOT / target).exists():
            broken.append(target)
    assert not broken, f"README points at files that do not exist: {broken}"


#: Figure stem -> the command that rebuilds it. Screenshots are named from an
#: f-string (``f"app_{name}.png"``), so a plain substring search over the
#: sources misses them; that is a fault in the checker, not in the figures, and
#: it is why this is a stated map rather than a scan. What the scan *did* find
#: was real: `gating_morph_small.gif` was committed with nothing to rebuild it,
#: which is now `make_animations.py --only readme`.
FIGURE_SOURCES = {
    "hero_human_piezo1": "make_figures.py --only hero",
    "curved_vs_flat": "make_figures.py --only dome",
    "domain_key": "make_figures.py --only key",
    "halotag_fold": "make_model_figures.py --only halotag",
    "hybrid_model": "make_model_figures.py --only hybrid",
    "record_nulls": "make_record_figure.py",
    "gating_morph_small": "make_animations.py --only readme",
    "app_modes": "screenshot_app.py --modes",
    "app_pore": "screenshot_app.py --analysis",
}


def test_the_figures_it_shows_are_regenerable(text):
    """A committed figure nobody can rebuild is a fossil.

    `docs/anim/` is git-ignored, so anything the README shows has to be
    committed — which makes it exactly the sort of file that quietly stops
    matching the code. Each one names the command that rebuilds it.
    """
    images = set(re.findall(r"!\[[^\]]*\]\((docs/img/[^)]+)\)", text))
    assert images, "the README shows no figures at all"

    undeclared = [i for i in images if Path(i).stem not in FIGURE_SOURCES]
    assert not undeclared, (
        f"README figures with no recorded way to rebuild them: "
        f"{sorted(undeclared)}. Add the command to FIGURE_SOURCES.")

    for stem, command in FIGURE_SOURCES.items():
        script = ROOT / "scripts" / command.split()[0]
        assert script.exists(), f"{stem} names a script that is gone: {script}"


# -------------------------------------------------------------- the citations

def _reference_section(text: str) -> str:
    return text[text.index("## References"):]


def test_there_is_a_reference_list_and_every_entry_carries_a_doi(text):
    section = _reference_section(text)
    entries = re.findall(r"^\s*\d+\. ", section, re.M)
    dois = re.findall(r"doi\.org/[^\)]+", section)
    assert len(entries) >= 15, f"only {len(entries)} references listed"
    assert len(dois) == len(entries), (
        f"{len(entries)} references but {len(dois)} DOIs; each must be "
        f"resolvable rather than an author-year the reader has to chase")


def test_every_citation_traces_to_the_verified_bibliography(text):
    """`docs/REFERENCES.md` is built behind a title-verification gate against
    Europe PMC. A citation absent from it has not been checked by anything."""
    verified = (ROOT / "docs" / "REFERENCES.md").read_text().lower()
    cited = re.findall(r"doi\.org/([^\)]+)", _reference_section(text))
    missing = [d for d in cited if d.lower() not in verified]
    assert not missing, f"citations not in the verified bibliography: {missing}"


def test_the_reference_list_points_at_the_full_bibliography(text):
    assert "docs/REFERENCES.md" in _reference_section(text)


def test_inline_citations_name_a_paper_in_the_list(text):
    """An "(Author et al. YEAR)" in the prose must be findable at the end."""
    body = text[:text.index("## References")]
    section = _reference_section(text).lower()
    cited = set(re.findall(r"\(([A-Z][a-zA-Z]+)(?: & [A-Z][a-zA-Z]+)?"
                           r"(?: et al\.)? (\d{4})[a-z]?\)", body))
    unresolved = [f"{name} {year}" for name, year in cited
                  if name.lower() not in section]
    assert not unresolved, f"cited in the prose, absent from the list: {unresolved}"


# ------------------------------------------- claims I got wrong when rewriting

def test_the_malaria_claim_is_stated_as_contested(text):
    """`docs/SCIENCE.md` records the association as disputed, and the README
    said it plainly until this test existed."""
    flowed = " ".join(text.split()).lower()
    assert "e756del" in flowed
    assert "contested" in flowed, (
        "the README states the malaria protection without the dispute that "
        "docs/SCIENCE.md records")


def test_the_allele_frequency_matches_the_projects_own_figure(text):
    """0.166-0.173 in gnomAD AFR — a sixth, not a third."""
    science = (ROOT / "docs" / "SCIENCE.md").read_text()
    assert "0.166" in science, "SCIENCE.md no longer states the frequency"
    flowed = " ".join(text.split())
    assert "0.166" in flowed or "sixth" in flowed.lower(), (
        "the README's E756del frequency does not match docs/SCIENCE.md")
    assert "a third of people" not in flowed.lower()


def test_the_conductance_disagreement_is_not_presented_as_agreement(text):
    """41 pS against a published 25-30 pS is the one number that misses, and
    two of its inputs have never been measured. Quoting it as a success would
    be the project's own worst failure mode."""
    flowed = " ".join(text.split())
    assert "16–94 pS" in flowed or "16-94 pS" in flowed, (
        "the README quotes a conductance without the range its unmeasured "
        "inputs span")
