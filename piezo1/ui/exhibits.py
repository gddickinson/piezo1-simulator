"""What an analysis window can *show* you, beside its own table.

Every entry in the Analysis menu opens a window of numbers, and a number is the
*end* of a piece of reasoning. The reasoning is what a reader has to see to
judge it: the figure the result is drawn in, the curve the model traces when
you move the parameter it depends on, the place on the structure it was
measured. Until this existed, all of that was in ``docs/img`` and in scripts,
reachable only by someone who already knew it was there — which is the same gap
Round 34 closed for the analyses themselves, one level up.

So each result window carries an **Explore** button, and each analysis declares
what exploring it means. There are four kinds of exhibit and they differ in
what they are *evidence of*, which is why the kind is shown on every one:

``figure``
    A generated picture from ``docs/img``. Regenerable and sometimes absent —
    a missing figure degrades to the command that rebuilds it, the way a tour
    step degrades to prose, because refusing to open would be worse.
``chart``
    Drawn from **the result already in the window**, never recomputed. The
    picture and the table are then guaranteed to be of one run; recomputing
    would make two, and the project has been bitten by exactly that (see
    ``pore_controller``, which reads the analysis object rather than re-running
    it).
``simulation``
    A model the user drives. The controls are registered parameters wherever
    one exists, passed **per call** — nothing here writes to the registry, so
    exploring cannot leave the application quoting non-default numbers.
``model``
    A button that turns on the corresponding overlay in the 3-D view. It draws
    nothing new: it calls the controller the View menu calls.

Every exhibit states its ``basis`` — whether what you are looking at was
measured here, imported, modelled or reproduced from a paper — and a
``not_this`` line saying what it must not be read as. A picture is more
persuasive than the number it came from, which is the standing hazard the
overlay controllers were written around, and it applies twice over to a curve
the user has just moved a slider to produce.

Qt-free, so the catalogue can be read, tested and audited without a display.
The window is :mod:`piezo1.ui.explore_window`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["Exhibit", "Series", "Reference", "ChartData", "EXHIBITS", "BASES",
           "KINDS", "exhibits_for", "figure_path", "empty_chart", "registry"]

#: What the thing on screen is evidence of. Shown on every exhibit, because
#: "measured on the entry you loaded" and "somebody else's published figure"
#: look identical once they are both pictures in the same window.
BASES = {
    "measured": "measured here, on coordinates this project holds",
    "modelled": "a model — no experiment fixes this",
    "imported": "an external project's result, imported with its provenance",
    "curated": "curated annotation, with the evidence level it carries",
    "published": "reproduced from the paper named beside it",
    "record": "a recorded result, frozen; recomputed only to check it",
}

KINDS = ("figure", "chart", "simulation", "model")


# --------------------------------------------------------------------------
# Chart data. Qt-free so a builder can be tested without a display, and so the
# same description can be drawn by the widget or dumped in a test.
# --------------------------------------------------------------------------

@dataclass
class Series:
    """One drawn set. ``kind`` is line, bar or point."""

    name: str
    x: list[float]
    y: list[float]
    kind: str = "line"
    color: str = ""
    axis: int = 0


@dataclass
class Reference:
    """A value the series has to be read against.

    ``high`` makes it a band rather than a line — a published range is a band,
    and drawing it as a single line would claim a precision the literature does
    not have.
    """

    value: float
    label: str = ""
    high: float | None = None
    vertical: bool = False
    color: str = "#8a919e"


@dataclass
class ChartData:
    """Everything the chart widget needs, and nothing about how to paint it."""

    title: str = ""
    x_label: str = ""
    y_label: str = ""
    series: list[Series] = field(default_factory=list)
    references: list[Reference] = field(default_factory=list)
    #: Bar charts label the x axis with these instead of numbers.
    categories: list[str] = field(default_factory=list)
    log_x: bool = False
    log_y: bool = False
    #: What the picture leaves out, or why there is nothing to draw.
    note: str = ""

    @property
    def empty(self) -> bool:
        return not any(len(s.x) for s in self.series)


def empty_chart(note: str, title: str = "") -> ChartData:
    """A chart with nothing in it and the reason why.

    Builders read a result dict that may legitimately lack the key they want —
    a shut pore has no current, an entry with no partner has no comparison. The
    honest answer is to say so on the panel. Returning a chart of zeros would
    draw a flat line that reads as a measurement of nothing happening.
    """
    return ChartData(title=title, note=note)


# --------------------------------------------------------------------------
# The exhibit itself
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Exhibit:
    """One thing an analysis window can show, beside its numbers."""

    analysis: str
    kind: str
    title: str
    #: One line: what you are looking at.
    what: str
    basis: str
    #: What it must not be read as. Never empty — see the module docstring.
    not_this: str
    #: kind == "figure": file name under docs/img, and the command that makes it.
    figure: str = ""
    rebuild: str = ""
    #: kind == "chart": key in :mod:`piezo1.ui.exhibit_plots`.
    plot: str = ""
    #: kind == "simulation": key in :mod:`piezo1.ui.exhibit_models`.
    simulation: str = ""
    #: kind == "model": key in :mod:`piezo1.ui.explore_window`'s action table.
    action: str = ""

    def figure_file(self):
        """Absolute path to the figure, or ``None`` when it is not built."""
        return figure_path(self.figure) if self.figure else None


def figure_path(name: str):
    """Resolve a ``docs/img`` name, or ``None`` if it has not been generated.

    Kept as a name in the catalogue rather than a path, so the catalogue does
    not depend on where the project is installed — the same reason
    :class:`piezo1.tour.TourStep` keeps one.
    """
    from ..config import PROJECT_ROOT

    path = PROJECT_ROOT / "docs" / "img" / name
    return path if path.exists() else None


def _validate(exhibits: tuple[Exhibit, ...]) -> tuple[Exhibit, ...]:
    """Refuse a malformed catalogue at import, rather than at the click.

    The same rule ``ui.theme`` applies to a missing token: an exhibit with no
    ``not_this`` line, or one naming no source for its content, is a window
    that opens empty in front of a user, and the failure would be invisible
    until someone clicked it.
    """
    required = {"figure": "figure", "chart": "plot", "simulation": "simulation",
                "model": "action"}
    for item in exhibits:
        if item.kind not in KINDS:
            raise ValueError(f"{item.title}: unknown exhibit kind {item.kind!r}")
        if item.basis not in BASES:
            raise ValueError(f"{item.title}: unknown basis {item.basis!r}")
        if not item.not_this.strip():
            raise ValueError(f"{item.title}: no 'must not be read as' line")
        field_name = required[item.kind]
        if not getattr(item, field_name):
            raise ValueError(f"{item.title}: {item.kind} exhibit names no "
                             f"{field_name}")
        if item.kind == "figure" and not item.rebuild:
            raise ValueError(f"{item.title}: no command to rebuild the figure")
    return exhibits


def _catalogue() -> tuple[Exhibit, ...]:
    from .exhibit_catalogue import MEASURED
    from .exhibit_catalogue_family import IMPORTED
    from .exhibit_catalogue_structure import STRUCTURAL

    return _validate(tuple(MEASURED) + tuple(IMPORTED) + tuple(STRUCTURAL))


def _by_analysis() -> dict[str, tuple[Exhibit, ...]]:
    out: dict[str, list[Exhibit]] = {}
    for item in _catalogue():
        out.setdefault(item.analysis, []).append(item)
    return {key: tuple(value) for key, value in out.items()}


_REGISTRY: dict[str, tuple[Exhibit, ...]] | None = None


def registry() -> dict[str, tuple[Exhibit, ...]]:
    """Analysis key -> what exploring it means, built once on first use.

    Built lazily rather than at import, because the catalogues import
    :class:`Exhibit` from here and this module reads them back: assembling at
    import time works only if *this* module is imported first, and
    ``tests/test_imports.py`` imports every module alone in a cleared
    ``sys.modules`` precisely because a cycle is invisible once the package is
    loaded. It found this one.
    """
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = _by_analysis()
    return _REGISTRY


def exhibits_for(analysis: str) -> tuple[Exhibit, ...]:
    """What the window for ``analysis`` can show. Empty is a real answer."""
    return registry().get(analysis, ())


def __getattr__(name: str):
    """``EXHIBITS`` as a module attribute, resolved on first read (PEP 562).

    Assembled from two catalogues, at the seam ``report_family.py`` uses: what
    this project measured, and what it imported from somebody else.
    """
    if name == "EXHIBITS":
        return registry()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
