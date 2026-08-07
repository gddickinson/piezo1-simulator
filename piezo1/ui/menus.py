"""Menu bar construction, kept out of the window shell.

Menus are mostly declarative and grow every round; leaving them inline was
pushing `main_window.py` past the project's 500-line limit for reasons that had
nothing to do with what the window *is*.
"""

from __future__ import annotations

from PyQt6.QtCore import QSettings
from PyQt6.QtGui import QAction, QActionGroup, QKeySequence

__all__ = ["build_menus", "SETTINGS_ORG", "SETTINGS_APP"]

SETTINGS_ORG = "piezo1-simulator"
SETTINGS_APP = "PIEZO1 Dynamic Structural Simulator"


def _action(parent, label: str, callback, shortcut: str = "",
            tip: str = "", checkable: bool = False, checked: bool = False):
    action = QAction(label, parent)
    if shortcut:
        action.setShortcut(QKeySequence(shortcut))
    if tip:
        action.setToolTip(tip)
        action.setStatusTip(tip)
    action.setCheckable(checkable)
    if checkable:
        action.setChecked(checked)
        action.toggled.connect(callback)
    else:
        action.triggered.connect(callback)
    parent.addAction(action)
    return action


def build_menus(win) -> None:
    """Build every menu on ``win``. Expects docks and controllers to exist."""
    bar = win.menuBar()
    # setToolTipsVisible is on QMenu, not QMenuBar — each menu enables it
    # individually. The tips carry what each analysis actually computes, so
    # leaving them at Qt's default of hidden would waste them.
    _file_menu(win, bar)
    _view_menu(win, bar)
    _analysis_menu(win, bar)
    _options_menu(win, bar)
    _help_menu(win, bar)


def _file_menu(win, bar) -> None:
    menu = bar.addMenu("&File")
    menu.setToolTipsVisible(True)
    _action(menu, "&Open structure…", win._open_file, "Ctrl+O",
            "Load an mmCIF or PDB file from disk")
    menu.addSeparator()
    _action(menu, "&Save session…", win.session.save, "Ctrl+S",
            "Save what you are looking at — structure, style, camera, "
            "selection. Never results.")
    _action(menu, "&Load session…", win.session.load, "Ctrl+L",
            "Restore a saved view")
    menu.addSeparator()
    _action(menu, "&Export analysis report…", win.session.export_report,
            "Ctrl+E",
            "Write a provenance-stamped report of the analyses that have run, "
            "as Markdown or JSON")
    menu.addSeparator()
    _action(menu, "&Quit", win.close, "Ctrl+Q", "Close the application")


def _view_menu(win, bar) -> None:
    menu = bar.addMenu("&View")
    menu.setToolTipsVisible(True)

    panels = menu.addMenu("&Panels")
    panels.setToolTipsVisible(True)
    for _key, action in win.docks.view_actions():
        panels.addAction(action)
    panels.addSeparator()
    _action(panels, "Show all panels", win.docks.show_all, "",
            "Bring every panel back into view")
    _action(panels, "Hide all panels", win.docks.hide_all, "",
            "Leave only the 3-D viewport")
    _action(panels, "Float all panels", lambda: win.docks.float_all(True), "",
            "Detach every panel into its own window")
    _action(panels, "Dock all panels", lambda: win.docks.float_all(False), "",
            "Return every floating panel into the main window")

    menu.addSeparator()
    _action(menu, "&Reset layout", win._reset_layout, "Ctrl+R",
            "Restore the panel arrangement the application ships with")
    _action(menu, "Reset &window size", win._size_to_screen, "",
            "Resize the window to fit the current screen")
    menu.addSeparator()
    _action(menu, "&Presentation mode (full screen)", win._toggle_fullscreen,
            "F11",
            "Hide the panels and menu so the 3-D view fills the screen")
    _action(menu, "&Display options…", win._show_display_options, "Ctrl+D",
            "Choose what the overlay shows: scale bar, animation clock, "
            "orientation axes and which measured values")
    menu.addSeparator()

    _action(menu, "&Ion flux animation", win.ion_flux.show, "",
            checkable=True, checked=False,
            tip="Animate ions crossing the pore at the rate the computed\n"
                "current sets. A channel passes ~10^7 ions/s, so the stream\n"
                "runs about a MILLIONFOLD slow and the HUD states the factor.\n"
                "A pore the wetting model calls shut shows no ions at all.")

    halotag = menu.addMenu("&HaloTag fusion")
    halotag.setToolTipsVisible(True)
    _action(halotag, "Show modelled &tags", win.fusion.show, "",
            checkable=True, checked=False,
            tip="Draw a HaloTag at each of the three cytosolic C-termini.\n"
                "THIS IS A MODEL: there is no structure of the fusion, so the\n"
                "tag body is drawn as a sphere of its radius of gyration and\n"
                "the linker as a straight seam.")
    _action(halotag, "Show accessible &volume", win.fusion.set_envelope, "",
            checkable=True, checked=False,
            tip="The region the tag centre can occupy without clashing, as a\n"
                "point cloud. Shown so a single sphere is not mistaken for a\n"
                "determined position.")
    _action(halotag, "Show &dyes", win.fusion.set_dyes, "",
            checkable=True, checked=False,
            tip="Draw a dye on each tag the labelling model says is occupied.")

    menu.addSeparator()

    _action(menu, "Show &multiple structures at once", win.set_multi_structure,
            checkable=True, checked=win.multi_structure,
            tip="Keep the current structure on screen when another is loaded, "
                "drawn in its own colour in the same frame. Off by default: "
                "two entries in the same frame sit on top of each other, and a "
                "structure left behind reads as extra density. Analyses always "
                "run on the primary structure, whatever else is drawn.")
    _action(menu, "Remove e&xtra structures", win.clear_companions, "",
            tip="Drop everything except the primary structure")

    menu.addSeparator()

    align = menu.addMenu("Structure &alignment")
    align.setToolTipsVisible(True)
    align_group = QActionGroup(align)
    align_group.setExclusive(True)
    for label, key, tip in (
            ("As deposited", "deposited",
             "Use the coordinate frame from the file. Different entries were "
             "refined in unrelated frames, so they will not overlap."),
            ("Canonical (three-fold axis on z)", "canonical",
             "Put each structure in a frame defined by its own C3 symmetry: "
             "axis vertical, cytosolic side down, centred on the origin. Works "
             "for any trimer, including PIEZO2 and mouse entries."),
            ("Superpose on the loaded structure", "reference",
             "Least-squares fit onto the first structure loaded, over the "
             "C-alphas they share. Maximises overlap, but needs a shared "
             "residue numbering — falls back to canonical across species.")):
        action = _action(align, label,
                         lambda on, k=key: on and win.set_alignment_mode(k),
                         checkable=True,
                         checked=(win.alignment_mode == key), tip=tip)
        align_group.addAction(action)

    menu.addSeparator()

    menu.addSeparator()
    _action(menu, "Reset &camera", win._reset_camera, "",
            "Reframe the model to fill the viewport")
    _action(menu, "&Clear highlight", lambda: win._highlight([], ""), "",
            "Remove the current residue selection")


def _analysis_menu(win, bar) -> None:
    menu = bar.addMenu("&Analysis")
    menu.setToolTipsVisible(True)
    _action(menu, "Measure &dome", win.physics.measure_dome, "",
            "Fit a sphere to the transmembrane surface. Curved structures "
            "should give about 9.7 nm against the published 10.2 nm.")
    _action(menu, "&Pore profile", win.analysis.compute_pore, "",
            "Pore radius along the conduction axis, with the hydrophobicity "
            "trace and a wetting verdict")
    _action(menu, "Find poc&kets", lambda: win.analysis.compute_pockets(10), "",
            "Alpha-sphere cavity detection with a burial filter")
    _action(menu, "&Conservation", win.analysis.compute_conservation, "",
            "Per-residue conservation across vertebrate PIEZO1 orthologs")
    _action(menu, "Coupling to the &gate (PRS)", win.analysis.compute_allostery,
            "", "Perturbation response scanning. Needs normal modes first.")
    menu.addSeparator()
    _action(menu, "&Ion permeation…", win.show_permeation, "",
            "1-D drift-diffusion through the measured pore, gated by the\n"
            "wetting verdict. Gives the unitary conductance and, when the pore\n"
            "is shut, EVERY mechanism shutting it. The in-pore diffusivity and\n"
            "ion radius are unmeasured, so the answer spans 16-94 pS across\n"
            "their plausible ranges against a published 25-30 pS.")
    _action(menu, "&Interactions…", win.show_interactions, "",
            "Hydrogen bonds, salt bridges, hydrophobic contacts, pi-stacking,\n"
            "cation-pi and disulfides, using published geometric criteria.")
    _action(menu, "&Modulators…", win.show_ligands, "",
            "Yoda1, Yoda2, Jedi1/2, Dooku1 and GsMTx4: chemistry, measured\n"
            "potency and what is known about where each binds.\n"
            "NO PIEZO structure with a bound modulator has been deposited, so\n"
            "every site is inferred from mutagenesis, docking or geometry.")
    _action(menu, "Variant &prediction record…", win.show_prediction_record, "",
            "What a variant score from this application is entitled to claim.\n"
            "The central claim - predicting gain- vs loss-of-function from\n"
            "structure - has FAILED three pre-registered tests. Shows all\n"
            "three, the power statement, and what the score may still be\n"
            "used for.")
    _action(menu, "R2456H vs wild &type…", win.show_paired_variant, "",
            "The only variant-versus-wild-type structural comparison this\n"
            "project can make: 8YFG is the one deposited entry that resolves\n"
            "its own mutation. Reported against how much wild-type entries\n"
            "differ among themselves, which is what makes a single pair\n"
            "interpretable. n = 1.")
    _action(menu, "&Variant structures…", win.show_variant_structures, "",
            "What the deposited variant structures can actually support.\n"
            "A null result: every deposited human PIEZO1 structure is closed,\n"
            "only one of four variant entries resolves its own mutation, and\n"
            "three of them share one set of coordinates.")
    menu.addSeparator()
    _action(menu, "HaloTag &labelling…", win.show_labelling, "",
            "Per-site and whole-channel labelling over time, and the\n"
            "1:2:3-dye mixture. Kinetics imported from halotag_binding_sim.")
    _action(menu, "Calcium &nanodomain…", win.show_nanodomain, "",
            "Free calcium where the tag sits when this channel opens, from the\n"
            "buffered-diffusion Green's function. Predicts ~114 uM at 4 nm\n"
            "against a 0.2 uM sensor Kd, so a BAPTA sensor is SATURATED\n"
            "whenever its own channel opens — meaning puncta brightness\n"
            "reports labelling and open probability, not calcium amplitude.")
    _action(menu, "HaloTag &geometry…", win.show_fusion_numbers, "",
            "Where a C-terminal HaloTag would sit: accessible volume, distance\n"
            "to the pore exit and clearance. Draw it with View > HaloTag fusion.")

    menu.addSeparator()
    _action(menu, "&Sequences…", win._show_sequences, "Ctrl+Shift+S",
            "Browse the protein and coding sequences, select onto the model, "
            "and compare sequences with alignment options")


def _options_menu(win, bar) -> None:
    menu = bar.addMenu("&Options")
    menu.setToolTipsVisible(True)
    settings = win.settings

    _action(menu, "Remember &layout on exit",
            win._set_remember_layout, "",
            "Reopen with the panel arrangement and window size you left",
            checkable=True,
            checked=settings.value("options/remember_layout", True, type=bool))
    _action(menu, "Show status-bar &hints", win._set_show_hints, "",
            "The mouse and keyboard reminder in the status bar",
            checkable=True,
            checked=settings.value("options/show_hints", True, type=bool))
    menu.addSeparator()

    focus = menu.addMenu("When something is &selected")
    focus.setToolTipsVisible(True)
    focus_group = QActionGroup(focus)
    focus_group.setExclusive(True)
    mode = win.focus_mode()
    for label, key, tip in (
            ("Keep the view still", "none",
             "Highlight the selection without moving the camera"),
            ("Centre on the selection", "centre",
             "Move the pivot to the selection, keeping the zoom"),
            ("Centre and zoom to the selection", "frame",
             "Move and zoom so the selection fills the viewport, keeping "
             "the current orientation")):
        action = _action(focus, label,
                         lambda on, k=key: on and win._set_focus_mode(k),
                         checkable=True, checked=(mode == key), tip=tip)
        focus_group.addAction(action)

    menu.addSeparator()

    spin = menu.addMenu("&Spin speed")
    spin.setToolTipsVisible(True)
    group = QActionGroup(spin)
    group.setExclusive(True)
    current = settings.value("options/spin_speed", 28.0, type=float)
    for label, value in (("Off", 0.0), ("Slow", 12.0), ("Normal", 28.0),
                         ("Fast", 60.0)):
        action = _action(spin, label,
                         lambda on, v=value: on and win._set_spin_speed(v),
                         checkable=True, checked=abs(current - value) < 1e-6,
                         tip=f"Rotate at {value:.0f} degrees per second")
        group.addAction(action)

    menu.addSeparator()
    _action(menu, "&Parameters…", win._show_parameters, "Ctrl+P",
            "Every number the calculations use, with its default, its unit "
            "and the paper it came from. Editable.")
    menu.addSeparator()
    _action(menu, "&Restore default options", win._reset_options, "",
            "Forget every remembered setting and layout")


def _help_menu(win, bar) -> None:
    menu = bar.addMenu("&Help")
    menu.setToolTipsVisible(True)
    _action(menu, "&Guided tour", win._start_tour, "F2",
            "Walk the mechanism — dome, blades, lever, gate — with every "
            "number measured live rather than quoted")
    _action(menu, "&Feature guide…", win._show_help, "F1",
            "What every panel does, and what the numbers mean")
    _action(menu, "&Keyboard and mouse…",
            lambda: win._show_help("shortcuts"), "",
            "Every shortcut in one table")
    menu.addSeparator()

    docs = menu.addMenu("&Documents")
    docs.setToolTipsVisible(True)
    from .help_content import DOC_LINKS
    for title, path, description in DOC_LINKS:
        _action(docs, title, lambda _=None, p=path: win._open_document(p),
                tip=description)

    menu.addSeparator()
    _action(menu, "What this application will &not do…",
            lambda: win._show_help("Limits"), "",
            "The blind test that failed, what it rules out, and the numbers "
            "that have since been corrected")
    _action(menu, "&About", win._about, "", "Version and provenance")


def make_settings() -> QSettings:
    return QSettings(SETTINGS_ORG, SETTINGS_APP)
