"""Ways a result can put something on the 3-D view — and one rule for all of them.

The point of this application is the structure, so a result window that can
only print numbers is doing half its job. These are the structural displays an
exhibit may ask for: load a named entry, draw a second one beside it, superpose
one on another (including on the pore module alone), show one named part of the
assembly, mark a set of residues, recolour, or build the gating morph.

**Every one presses a control the user could press.** The Model panel's entry
selector, the Overlay panel's combos and button, the View menu's component
group, the Model panel's colour combo, the Physics panel's morph button. Not
the controllers behind them — otherwise an overlay would be on screen while the
menu that owns it said it was off, and there would be two answers to what is
being drawn. The cost is that this table can drift when a control is renamed,
which is why a test resolves every entry against a real window.

The table of displays itself — which ones are offered and which control each
belongs to — is in :mod:`piezo1.ui.model_action_table`, split off at the
length limit along the seam ``parameter_audit_exemptions`` uses: this is the
mechanism, that is the judgement.

**Residue numbers are the other hazard.** A highlight is a set of residue
numbers, and this project has been bitten three times by numbers read in the
wrong system — most recently in ``analysis/features.py`` (Round 93) and in the
Annotation panel itself. So a highlight resolves through the loaded entry's own
numbering: curated groups come from the annotation already loaded for that
entry, and human-numbered records from the census are converted through
``core.sequence``. An entry whose numbering cannot be read is **refused**, with
the reason, rather than highlighted at whatever sits at those numbers.
"""

from __future__ import annotations

from dataclasses import dataclass

from PyQt6.QtGui import QAction

__all__ = ["ModelActionSpec", "BoundAction", "MODEL_ACTIONS", "ACTION_KINDS",
           "CENSUS_MARKS", "find_menu_action", "find_button",
           "highlight_residues"]


#: Imported positions have their own View-menu entries, so they are reachable
#: without opening a result window — and so an exhibit can press the entry
#: rather than calling the highlight itself.
CENSUS_MARKS = {
    "family:pathogenic_pore": "&Pathogenic pore positions",
    "family:equivalent": "Positions &equivalent to a PIEZO2 disease residue",
}

#: What a model action can do. Each kind names the control it presses.
ACTION_KINDS = ("control", "load", "companion", "overlay", "component",
                "highlight", "colour", "morph")


def __getattr__(name: str):
    """``MODEL_ACTIONS`` lives in :mod:`piezo1.ui.model_action_table`.

    Resolved on first read (PEP 562) rather than imported at the top, because
    the table imports :class:`ModelActionSpec` from here — the same cycle the
    exhibit catalogues have, and the same fix.
    """
    if name == "MODEL_ACTIONS":
        from .model_action_table import MODEL_ACTIONS

        return MODEL_ACTIONS
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


@dataclass(frozen=True)
class ModelActionSpec:
    """Which existing control an exhibit's button should press.

    ``kind`` decides how ``target`` is read: a PDB id for ``load``,
    ``companion`` and ``overlay``; a component key for ``component``; a colour
    name for ``colour``; a source name for ``highlight`` (see
    :func:`highlight_residues`).
    """

    label: str
    kind: str = "control"
    menu: str = ""
    button: str = ""
    target: str = ""
    mode: str = ""


def _plain(text: str) -> str:
    return text.replace("&", "").strip().lower()


def find_menu_action(window, label: str):
    """The menu entry with this text, or ``None``."""
    if window is None or not hasattr(window, "menuBar"):
        return None
    for action in window.menuBar().findChildren(QAction):
        if _plain(action.text()) == _plain(label):
            return action
    return None


def find_button(window, path: str):
    target = window
    for part in path.split("."):
        target = getattr(target, part, None)
        if target is None:
            return None
    return target


def _component_action(window, key: str):
    """The View-menu entry for a component, found by its own label."""
    from ..structure.components import component_by_key

    component = component_by_key(key)
    return find_menu_action(window, component.label) if component else None


def highlight_residues(window, source: str, result: dict | None = None):
    """Residue numbers for a highlight, **in the loaded entry's numbering**.

    Returns ``(residues, label)``, or ``([], reason)`` when the entry's
    numbering cannot carry them. Two sources:

    ``site:<id>``
        A curated residue group, read from the annotation already loaded for
        this entry — so it is mouse 2473/2476/2480 on a mouse entry and human
        2447/2450/2454 on a human one, rather than one of those on both.
    ``family:<name>``
        Positions imported from the census, which are recorded in **human**
        numbering and converted here through ``core.sequence`` — never by
        arithmetic, because the offset is not constant.
    """
    from ..core.numbering_check import piezo1_numbering

    structure = getattr(window, "structure", None)
    if structure is None:
        return [], "load a structure first"

    if source.startswith("site:"):
        annotations = getattr(window, "annotations", None)
        group = annotations.group(source.split(":", 1)[1]) if annotations else None
        if group is None or not group.residues:
            return [], ("no curated residue group of that name in this entry's "
                        "numbering — this protein has none")
        return list(group.residues), group.label

    numbering = piezo1_numbering(structure)
    if numbering is None:
        return [], ("this entry is not PIEZO1 in a numbering the annotation "
                    "can be read at, so these positions are refused rather "
                    "than drawn at whatever sits at those numbers")

    if source.startswith("variant:"):
        annotations = getattr(window, "annotations", None)
        label = source.split(":", 1)[1]
        variant = next((v for v in getattr(annotations, "variants", [])
                        if v.label == label), None)
        if variant is None:
            return [], f"{label} is not in this entry's curated variant set"
        if variant.position is None:
            return [], f"{label} has no position in this entry's numbering"
        return [int(variant.position)], variant.label

    from ..core.family import load_family_findings
    from ..core.sequence import human_to_mouse

    findings = load_family_findings()
    name = source.split(":", 1)[1]
    if name == "pathogenic_pore":
        human = [p.resi for p in findings.pathogenic_pore if p.gene == "PIEZO1"]
        label = "pathogenic pore-module positions (census)"
    elif name == "equivalent":
        human = [e.piezo1 for e in findings.equivalent]
        label = "positions equivalent to a PIEZO2 disease residue"
    else:
        return [], f"no imported position list named {name!r}"

    if numbering == "human":
        return sorted(human), label
    converted = [human_to_mouse(r) for r in human]
    kept = sorted(r for r in converted if r is not None)
    if not kept:
        return [], "none of these positions could be carried into this entry"
    return kept, f"{label}, converted to {numbering} numbering"


class BoundAction:
    """A model action resolved against one window."""

    def __init__(self, spec: ModelActionSpec, window, result=None) -> None:
        self.spec, self.window, self.result = spec, window, result
        self.label = spec.label
        self.action = None
        self.button = None
        if spec.kind == "control":
            self.action = (find_menu_action(window, spec.menu) if spec.menu
                           else None)
            self.button = find_button(window, spec.button) if spec.button else None
        elif spec.kind == "component":
            self.action = _component_action(window, spec.target)
        elif spec.kind == "colour":
            self.button = find_button(window, "structure_panel.color_combo")
        elif spec.kind == "load":
            self.button = find_button(window, "structure_panel.structure_combo")
        elif spec.kind == "overlay":
            self.button = find_button(window, "overlay_panel.load_button")
        elif spec.kind == "morph":
            self.button = find_button(window, "physics_panel.morph_button")
        elif spec.kind == "companion":
            self.button = find_button(window, "add_companion")
        elif spec.kind == "highlight":
            # A highlight presses the control that owns it: a curated group is
            # a row in the Annotation panel's Sites list, a variant is a row in
            # its table, and an imported position list is a View-menu entry.
            # Pressing those also *moves* them, so the user can see where the
            # selection came from and switch it off the same way.
            if spec.target in CENSUS_MARKS:
                self.action = find_menu_action(window, CENSUS_MARKS[spec.target])
            elif spec.target.startswith("site:"):
                self.button = find_button(window, "annotation_panel.site_list")
            elif spec.target.startswith("variant:"):
                self.button = find_button(window,
                                          "annotation_panel.variant_table")

    @property
    def resolved(self) -> bool:
        return self.action is not None or self.button is not None

    # ------------------------------------------------------------- pressing

    def run(self) -> str:
        """Press the control, and say what happened."""
        if not self.resolved:
            return "no control in this window draws that"
        handler = getattr(self, f"_run_{self.spec.kind}", None)
        return handler() if handler else self._run_control()

    def _run_control(self) -> str:
        target = self.action if self.action is not None else self.button
        if not target.isEnabled():
            return (f"'{self.spec.menu or self.spec.button}' is not available "
                    f"yet — it usually needs a structure loaded, or the normal "
                    f"modes computed first.")
        if getattr(target, "isCheckable", None) and target.isCheckable():
            if target.isChecked():
                return "already drawn — see the 3-D view."
            target.setChecked(True)
        else:
            target.trigger()
        return ("drawn. The status line under the 3-D view carries what this "
                "overlay must be read with; switch it off the same way you "
                "would from the View menu.")

    def _run_component(self) -> str:
        if self.action is None:
            return "no component entry of that name in the View menu"
        if self.action.isChecked():
            return "already showing that part — see the 3-D view."
        self.action.setChecked(True)
        return ("showing that part only. It HIDES rather than subsets: every "
                "analysis still runs on the whole assembly, which is what the "
                "status line says.")

    def _run_load(self) -> str:
        panel = getattr(self.window, "structure_panel", None)
        if panel is None:
            return "no structure selector in this window"
        panel.select(self.spec.target)
        record = getattr(self.window, "record", None)
        if record is None or record.pdb.upper() != self.spec.target.upper():
            return (f"{self.spec.target} is not available locally — run "
                    f"python -m piezo1.io.fetch")
        return (f"loaded {self.spec.target}. Every analysis now runs on it, "
                f"and results computed on the previous entry are stale — the "
                f"Analysis panel says so where one is.")

    def _run_companion(self) -> str:
        add = getattr(self.window, "add_companion", None)
        if add is None:
            return "this window cannot draw a second structure"
        add(self.spec.target)
        return (f"{self.spec.target} is drawn alongside, in its own colour and "
                f"the same frame. The analyses still run on the primary "
                f"structure only.")

    def _run_overlay(self) -> str:
        panel = getattr(self.window, "overlay_panel", None)
        if panel is None:
            return "no overlay panel in this window"
        combo = panel.structure_combo
        for index in range(combo.count()):
            if str(combo.itemData(index)).upper() == self.spec.target.upper():
                combo.setCurrentIndex(index)
                break
        else:
            return (f"{self.spec.target} is not in the overlay list — it may "
                    f"not be downloaded, or it is the structure already loaded")
        for index in range(panel.mode_combo.count()):
            if panel.mode_combo.itemData(index) == (self.spec.mode or "protomer"):
                panel.mode_combo.setCurrentIndex(index)
                break
        panel.load_button.click()
        if self.spec.mode == "core":
            return ("superposing on the pore module alone. The blades are "
                    "then a MEASUREMENT rather than part of the fit, and the "
                    "status line reports how far they land apart.")
        return "superposing. The status line reports the fit when it finishes."

    def _run_colour(self) -> str:
        combo = self.button
        for index in range(combo.count()):
            if combo.itemText(index) == self.spec.target:
                combo.setCurrentIndex(index)
                return (f"coloured by {self.spec.target}. The Model panel's "
                        f"own selector moved with it.")
        return f"this build has no '{self.spec.target}' colouring"

    def _run_highlight(self) -> str:
        target = self.spec.target
        if target in CENSUS_MARKS:
            if self.action is None:
                return "no menu entry marks those positions"
            self.action.trigger()
            return self._marked("These positions are the census's, converted "
                                "into this entry's numbering.")
        if target.startswith("site:"):
            return self._press_site(target.split(":", 1)[1])
        if target.startswith("variant:"):
            return self._press_variant(target.split(":", 1)[1])
        return f"no way to mark {target!r}"

    def _press_site(self, group_id: str) -> str:
        from PyQt6.QtCore import Qt

        panel = getattr(self.window, "annotation_panel", None)
        listing = self.button
        if panel is None or listing is None:
            return "no annotation panel in this window"
        for row in range(listing.count()):
            item = listing.item(row)
            if item.data(Qt.ItemDataRole.UserRole) == group_id:
                panel.tabs.setCurrentIndex(1)
                listing.setCurrentRow(row)
                return self._marked("The Annotation panel's Sites tab moved "
                                    "with it, which is where the group and its "
                                    "evidence are.")
        return ("no curated group of that name in this entry's numbering — "
                "this protein has none")

    def _press_variant(self, label: str) -> str:
        panel = getattr(self.window, "annotation_panel", None)
        table = self.button
        if panel is None or table is None:
            return "no annotation panel in this window"
        panel.class_filter.setCurrentIndex(0)
        panel.search.setText("")
        rows = getattr(panel, "_variant_rows", [])
        for row, variant in enumerate(rows):
            if variant.label == label:
                panel.tabs.setCurrentIndex(2)
                table.selectRow(row)
                return self._marked(
                    f"{label} keeps its published name and is marked at its "
                    f"position in this entry's numbering.")
        return f"{label} is not in this entry's curated variant set"

    def _marked(self, note: str) -> str:
        residues = getattr(self.window, "selected_residues", []) or []
        if not residues:
            status = getattr(self.window, "status_label", None)
            return (status.text() if status is not None and status.text()
                    else "nothing was marked — see the status line")
        plural = "residue" if len(residues) == 1 else "residues"
        return (f"marked {len(residues)} {plural} on all three protomers. "
                f"{note}")

    def _run_morph(self) -> str:
        if self.button is None or not self.button.isEnabled():
            return "the morph needs a structure with a curved/flat partner"
        self.button.click()
        return ("building the interpolation between the curved and flattened "
                "endpoints. It is an INTERPOLATION, not a simulated "
                "trajectory; scrub it with the slider in the Physics panel.")
