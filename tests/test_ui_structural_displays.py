"""What a result can put on the 3-D view, and whether it lands where it says.

The point of this application is that PIEZO1's shape is its mechanism, so a
result window that can only print numbers is doing half its job. These cover
the half that draws: loading the entry a result is about, superposing a
partner, showing one component, marking a residue set, recolouring, and the
pore-module-only superposition the family results are made of.

Two of these guard a defect that was live in the application when they were
written, and both are the same shape — **a residue number read in the wrong
numbering, landing on a real, wrong residue**:

* the Annotation panel loaded human annotation once, at construction, and
  never asked again, so on a mouse entry — which is most of the catalogue —
  every domain range, site and variant in it was a human number on mouse
  coordinates;
* and a variant's ``residue`` is its *identity* (R2456H is R2456H in every
  paper), which is not its position: on 7WLT that position is 2482, and 2456
  there is a proline.

Both are calibrated against the wrong answer rather than only the right one.
"""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core.annotations import load_annotations  # noqa: E402
from piezo1.ui.model_actions import (ACTION_KINDS, MODEL_ACTIONS,  # noqa: E402
                                     BoundAction, highlight_residues)


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance()
    if app is None:
        try:
            app = QApplication([])
        except Exception as exc:                       # pragma: no cover
            pytest.skip(f"no Qt platform available: {exc}")
    return app


@pytest.fixture(scope="module")
def window(qapp):
    from piezo1.ui.main_window import MainWindow

    win = MainWindow()
    yield win
    for controller in ("analysis", "physics", "overlay"):
        cleanup = getattr(getattr(win, controller, None), "cleanup", None)
        if cleanup is not None:
            cleanup()
    qapp.processEvents()


# --------------------------------------------------------------------------
# The numbering the panel reads annotation at
# --------------------------------------------------------------------------

def test_the_annotation_panel_reads_the_entry_s_own_numbering(window):
    """The defect: human annotation on a mouse structure.

    Calibrated on the gate, where the two numberings disagree by 26 residues
    and both numbers exist in the file — so the wrong one is not an error, it
    is a different residue.
    """
    panel = window.annotation_panel
    panel.set_species("mouse")
    gate = next(g for g in panel.annotations.residue_groups
                if "hydrophobic gate" in g.label.lower())
    assert tuple(gate.residues) == (2473, 2476, 2480)

    panel.set_species("human")
    gate = next(g for g in panel.annotations.residue_groups
                if "hydrophobic gate" in g.label.lower())
    assert tuple(gate.residues) == (2447, 2450, 2454)


def test_a_protein_with_no_curated_annotation_gets_none_of_it(window):
    """PIEZO2, PEZO-1 and dPIEZO have no curated domains here. Showing
    PIEZO1's would be worse than showing nothing, because it would look
    exactly like an answer."""
    panel = window.annotation_panel
    panel.set_species("mouse_piezo2")
    assert panel.domain_list.count() == 0
    assert panel.site_list.count() == 0
    assert panel.variant_table.rowCount() == 0
    assert "No curated annotation" in panel.info.text()
    panel.set_species("human")
    assert panel.domain_list.count() > 0


def test_a_variant_keeps_its_name_and_moves_its_position():
    """R2456H is R2456H in every document this project holds. Its *position*
    on a mouse entry is 2482 — and 2456 there is a proline, which is why the
    old behaviour never looked broken."""
    human = next(v for v in load_annotations("human").variants
                 if v.label == "R2456H")
    mouse = next(v for v in load_annotations("mouse").variants
                 if v.label == "R2456H")
    assert human.label == mouse.label == "R2456H"
    assert human.residue == mouse.residue == 2456, "the identity must not move"
    assert human.position == 2456 and mouse.position == 2482


def test_the_wrong_number_lands_on_a_real_wrong_residue():
    """The calibration for the test above, on coordinates.

    If the human number simply did not exist in the mouse file this would
    have been caught the first time anyone clicked it.
    """
    from piezo1.core import Structure

    path = STRUCTURE_DIR / "7WLT.cif"
    if not path.exists():
        pytest.skip("7WLT.cif not downloaded — run python -m piezo1.io.fetch")
    structure = Structure.from_file(path)
    variant = next(v for v in load_annotations("mouse").variants
                   if v.label == "R2456H")

    at_position = set(structure.res_name[structure.res_seq == variant.position])
    at_identity = set(structure.res_name[structure.res_seq == variant.residue])
    assert at_position == {"ARG"}, "the variant's own residue is an arginine"
    assert at_identity and at_identity != {"ARG"}, (
        "the human number should land on a different, real residue")


# --------------------------------------------------------------------------
# Every structural display resolves to a control
# --------------------------------------------------------------------------

def test_every_action_names_a_known_kind():
    for key, spec in MODEL_ACTIONS.items():
        assert spec.kind in ACTION_KINDS, (key, spec.kind)
        assert spec.label.strip()


def test_every_structural_display_resolves_to_a_control(window):
    """The one-control rule, checked against a real window.

    If this table drifts — a menu entry renamed, a panel button moved — the
    exhibit would silently do nothing, and the only symptom would be a button
    that appears to work.
    """
    for key, spec in MODEL_ACTIONS.items():
        bound = BoundAction(spec, window, result={})
        assert bound.resolved, f"{key} resolves to no control ({spec})"


def test_every_named_entry_component_and_colour_exists(window):
    """A display naming a structure nobody has, a component that does not
    exist, or a colouring this build has no entry for."""
    from piezo1.structure.components import component_by_key
    from piezo1.ui.panels.structure_panel import COLOR_LABELS

    known_colours = {label for label, _ in COLOR_LABELS}
    catalogued = {r.pdb.upper() for r in window.registry.entries}
    for key, spec in MODEL_ACTIONS.items():
        if spec.kind in ("load", "companion", "overlay"):
            assert spec.target.upper() in catalogued, (key, spec.target)
        elif spec.kind == "component":
            assert component_by_key(spec.target) is not None, (key, spec.target)
        elif spec.kind == "colour":
            assert spec.target in known_colours, (key, spec.target)


# --------------------------------------------------------------------------
# Highlights, on the structure
# --------------------------------------------------------------------------

def _load(window, pdb):
    if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
        pytest.skip(f"{pdb}.cif not downloaded")
    from piezo1.core import Structure
    from piezo1.io.registry import load_registry

    # No GL context offscreen, so `load_structure` returns early; the parts a
    # highlight reads are the structure, the record and the annotation, which
    # is what this sets.
    record = load_registry().get(pdb)
    window.structure = Structure.from_file(record.path)
    window.record = record
    window.annotations = load_annotations(record.numbering_species)
    window.annotation_panel.set_species(record.numbering_species)
    return record


def test_a_curated_group_is_marked_in_the_entrys_own_numbering(window):
    record = _load(window, "7WLT")
    assert record.numbering_species == "mouse"
    residues, label = highlight_residues(window, "site:hydrophobic_gate")
    assert tuple(residues) == (2473, 2476, 2480), label

    _load(window, "8YEZ")
    residues, _ = highlight_residues(window, "site:hydrophobic_gate")
    assert tuple(residues) == (2447, 2450, 2454)


def test_imported_human_positions_are_converted_not_copied(window):
    """The census records positions in human numbering. On a mouse entry they
    are carried across by the alignment map — never by an offset, because the
    offset is not constant."""
    from piezo1.core.family import load_family_findings

    human_positions = sorted(p.resi for p in load_family_findings().pathogenic_pore
                             if p.gene == "PIEZO1")
    _load(window, "8YEZ")
    on_human, _ = highlight_residues(window, "family:pathogenic_pore")
    assert on_human == human_positions

    _load(window, "7WLT")
    on_mouse, label = highlight_residues(window, "family:pathogenic_pore")
    assert on_mouse != human_positions, "nothing was converted"
    assert "mouse" in label
    offsets = {b - a for a, b in zip(human_positions, on_mouse)}
    assert len(offsets) > 1, ("a single offset would mean this is arithmetic "
                              "rather than an alignment")


def test_an_entry_the_annotation_cannot_be_read_at_is_refused(window):
    """A PIEZO2 structure must not be marked at PIEZO1's numbers — the
    reading that would put a 'pathogenic position' on whatever sits there."""
    _load(window, "6KG7")
    residues, reason = highlight_residues(window, "family:pathogenic_pore")
    assert residues == []
    assert "refused" in reason


# --------------------------------------------------------------------------
# The pore-module-only superposition
# --------------------------------------------------------------------------

def test_the_core_only_overlay_is_the_analysis_not_a_second_fit(window):
    """What is drawn has to be the superposition the result window describes,
    so the worker calls `core_fit` rather than fitting again."""
    from piezo1.analysis.core_periphery import core_fit
    from piezo1.core import Structure
    from piezo1.ui.overlay_controller import OverlayWorker

    for pdb in ("7WLT", "6KG7"):
        if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
            pytest.skip(f"{pdb}.cif not downloaded")
    reference = Structure.from_file(STRUCTURE_DIR / "7WLT.cif")
    mobile = Structure.from_file(STRUCTURE_DIR / "6KG7.cif")

    result = OverlayWorker(reference, mobile, "core")._superpose()
    fit = core_fit(mobile, reference, mobile.name, reference.name)
    assert result.rmsd == pytest.approx(fit.comparison.core_rmsd)
    assert result.n_common == fit.comparison.n_core
    assert result.meta["splay_ratio"] == pytest.approx(
        fit.comparison.splay_ratio)
    # Painted on the reference, so keyed by the reference's residue numbers.
    assert set(result.per_residue) == set(fit.deviation_target)


def test_the_core_mode_answers_a_pair_the_others_refuse(window):
    """The reason the mode exists. Matching PIEZO1 to PIEZO2 by residue number
    is refused — correctly, it gives a confident 47.9 Å — and the core fit
    corresponds them through a real alignment instead."""
    from piezo1.core import Structure
    from piezo1.ui.overlay_controller import OverlayController, OverlayWorker

    for pdb in ("7WLT", "6KG7"):
        if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
            pytest.skip(f"{pdb}.cif not downloaded")
    window.structure = Structure.from_file(STRUCTURE_DIR / "7WLT.cif")
    window.record = window.registry.get("7WLT")
    mobile = Structure.from_file(STRUCTURE_DIR / "6KG7.cif")

    controller = OverlayController(window)
    refusal = controller._numbering_refusal("6KG7", mobile)
    assert refusal, "the by-number modes should refuse this pair"

    result = OverlayWorker(window.structure, mobile, "core")._superpose()
    assert result.rmsd < 6.0
    assert "pore module" in result.summary()
    assert result.meta["cross_paralogue"] is True


def test_the_core_mode_reports_the_blades_as_a_measurement(window):
    """The whole construction: the blades are not in the fit, so the status
    line must not read as though they were."""
    from piezo1.core import Structure
    from piezo1.ui.overlay_controller import OverlayWorker

    for pdb in ("7WLT", "7WLU"):
        if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
            pytest.skip(f"{pdb}.cif not downloaded")
    reference = Structure.from_file(STRUCTURE_DIR / "7WLT.cif")
    mobile = Structure.from_file(STRUCTURE_DIR / "7WLU.cif")

    result = OverlayWorker(reference, mobile, "core")._superpose()
    summary = result.summary()
    assert "blades" in summary and "MEASUREMENT" in summary
    # PIEZO1's own gating pair is the extreme case of core-conserved,
    # periphery-free — it is what makes the ratio worth reporting at all.
    assert result.meta["splay_ratio"] > 3.0, result.meta["splay_ratio"]


def test_the_panel_reports_both_numbers_for_a_core_fit(window, qapp):
    """Reporting the core RMSD alone would read as "these two agree to 3.7 Å",
    which is true of the pore module and not of the protein."""
    from piezo1.core import Structure
    from piezo1.ui.overlay_controller import OverlayWorker

    for pdb in ("7WLT", "6KG7"):
        if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
            pytest.skip(f"{pdb}.cif not downloaded")
    reference = Structure.from_file(STRUCTURE_DIR / "7WLT.cif")
    mobile = Structure.from_file(STRUCTURE_DIR / "6KG7.cif")
    result = OverlayWorker(reference, mobile, "core")._superpose()

    window.overlay_panel.set_result(result)
    text = window.overlay_panel.result_label.text()
    assert "Core" in text and "blades" in text and "splay" in text
    assert "measured" in text and "fitted" in text


def test_the_mode_is_offered_in_the_panel_a_user_can_reach(window):
    """A display nobody can trigger without opening a result window would be
    half a feature."""
    from piezo1.ui.panels.overlay_panel import SUPERPOSITION_MODES

    assert "core" in {key for _label, key in SUPERPOSITION_MODES}
    modes = {window.overlay_panel.mode_combo.itemData(i)
             for i in range(window.overlay_panel.mode_combo.count())}
    assert "core" in modes


# --------------------------------------------------------------------------
# Pressing them
# --------------------------------------------------------------------------

def test_showing_a_component_moves_the_menu_the_user_would_have_clicked(window):
    from piezo1.ui.model_actions import BoundAction

    bound = BoundAction(MODEL_ACTIONS["component_pore"], window)
    assert bound.action is not None
    message = bound.run()
    assert bound.action.isChecked()
    assert "HIDES rather than subsets" in message
    BoundAction(MODEL_ACTIONS["component_whole"], window).run()


def test_recolouring_moves_the_model_panels_own_selector(window):
    from piezo1.ui.model_actions import BoundAction

    combo = window.structure_panel.color_combo
    before = combo.currentText()
    message = BoundAction(MODEL_ACTIONS["colour_bfactor"], window).run()
    assert combo.currentText() == "B-factor", message
    assert combo.currentText() != before or before == "B-factor"


def test_an_action_that_needs_something_first_says_so_rather_than_failing(window):
    from piezo1.ui.model_actions import BoundAction

    button = window.physics_panel.fluctuation_button
    button.setEnabled(False)
    assert "not available" in BoundAction(
        MODEL_ACTIONS["colour_fluctuation"], window).run()
