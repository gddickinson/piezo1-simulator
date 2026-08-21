"""Menu bar construction, kept out of the window shell.

Menus are mostly declarative and grow every round; leaving them inline was
pushing `main_window.py` past the project's 500-line limit for reasons that had
nothing to do with what the window *is*.
"""

from __future__ import annotations

from PyQt6.QtCore import QSettings
from PyQt6.QtGui import QAction, QKeySequence

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
    _action(menu, "Export &coloured structure…", win.session.export_scalar, "",
            "Write the value currently colouring the model into a PDB\n"
            "B-factor column, so PyMOL or ChimeraX can colour by it.\n"
            "Residues the analysis did not score go out with occupancy\n"
            "0.00 rather than a zero, so they stay distinguishable.")
    _action(menu, "&Export analysis report…", win.session.export_report,
            "Ctrl+E",
            "Write a provenance-stamped report of the analyses that have run, "
            "as Markdown or JSON")
    menu.addSeparator()
    _action(menu, "&Quit", win.close, "Ctrl+Q", "Close the application")


def build_component_menu_stub(win, menu) -> None:
    from .menus_flux import build_component_menu
    from .menus_styles import build_component_style_menu

    build_component_menu(win, menu, _action)
    build_component_style_menu(win, menu, _action)


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
    menu.addSeparator()
    build_component_menu_stub(win, menu)
    menu.addSeparator()

    _action(menu, "&Ion flux animation", win.ion_flux.show, "",
            checkable=True, checked=False,
            tip="Animate ions crossing the pore at the rate the computed\n"
                "current sets. A channel passes ~10^7 ions/s, so the stream\n"
                "runs about a MILLIONFOLD slow and the HUD states the factor.\n"
                "A pore the wetting model calls shut shows no ions at all -\n"
                "which is 17 of the 19 deposited entries on the default\n"
                "AXIAL pathway, because PIEZO1's axis is closed at both ends.")
    from .menus_flux import build_flux_settings, build_pore_opacity_menu
    build_flux_settings(win, menu, _action)
    build_pore_opacity_menu(win, menu, _action)

    _action(menu, "&Contacts", win.contacts.show, "",
            checkable=True, checked=False,
            tip="Draw the contacts the interaction analysis finds, between the\n"
                "atoms it found them between — same cutoffs, same geometry, no\n"
                "second implementation. Colour by kind, because they are not\n"
                "the same evidence: a disulfide is a covalent bond the map\n"
                "resolved, a salt bridge is two charged groups inside a\n"
                "published cutoff, a hydrophobic contact is two carbons near\n"
                "each other. Hydrophobic contacts are OFF by default; they\n"
                "outnumber the rest several times over and would bury them.\n"
                "Criteria are heavy-atom based: no deposited entry has\n"
                "hydrogens, so a drawn hydrogen bond is geometry rather than\n"
                "an observed proton.")
    _action(menu, "Colour by evolutionary c&onstraint", win.constraint_colour.show, "",
            checkable=True, checked=False,
            tip="How much 174 PIEZO1 orthologues have refused to change each\n"
                "residue, from the piezo_genes census. NOT MEASURED HERE.\n"
                "Scale FIXED at 0-1 so two entries stay comparable; an\n"
                "auto-ranged map would repaint the same protein differently\n"
                "depending on how much blade the entry resolved.\n"
                "An UNSCORED residue is grey, not dark - the blade tips are\n"
                "where coverage is worst and where low constraint is exactly\n"
                "the claim being made, so the two must not share a colour.\n"
                "A mouse entry is read through the alignment map. A PIEZO2 or\n"
                "invertebrate entry is REFUSED rather than coloured by\n"
                "whatever sits at those numbers in PIEZO1.")
    _action(menu, "Colour by &electrostatics", win.electrostatics.show, "",
            checkable=True, checked=False,
            tip="Figure 4c's surface potential, on the SAME FIXED SCALE:\n"
                "red at -5 k_BT/e, white at zero, blue at +5. Fixed, not\n"
                "auto-ranged - an auto-ranged potential map paints an almost\n"
                "neutral protein in full red and blue and cannot be compared\n"
                "with a published surface.\n"
                "NOT APBS. Screened Coulomb from formal charges through a\n"
                "uniform dielectric: no dielectric boundary, no ion-exclusion\n"
                "layer, no partial charges. All three under-estimate the\n"
                "magnitude, and on 6B3R nothing reaches the saturation the\n"
                "published panel visibly reaches. Read the sign and the\n"
                "pattern, not the value.")
    _action(menu, "&Micelle density (modelled)", win.micelle.show, "",
            checkable=True, checked=False,
            tip="Figure 4b's detergent envelope, MODELLED rather than\n"
                "observed. The published panel is the unsharpened cryo-EM map\n"
                "at 6 sigma; this project holds no map, so what is drawn is\n"
                "the surface a fixed distance outside the hydrophobic\n"
                "transmembrane belt. The SHELL THICKNESS is a parameter and\n"
                "carries no information. The CURVATURE is a sphere fitted to\n"
                "the belt atoms themselves and is a measurement of the\n"
                "protein: 9.8 nm on 6B3R against the paper's 10.2 nm\n"
                "idealisation. The status line says which is which.")
    _action(menu, "Planar &membrane (one protomer)",
            win.planar_membrane.show, "", checkable=True, checked=False,
            tip="Figure 4a: a single subunit with the two planar membrane\n"
                "interfaces drawn across it. The paper's point is the\n"
                "CONTRAST — a protomer sits in a plane and the trimer does\n"
                "not — so the status line reports both residuals and the slab\n"
                "thickness each would need against a real 36 A bilayer.\n"
                "Every point set has a best-fit plane, so read the residual\n"
                "rather than the lines.")
    _action(menu, "&Dome surface", win.dome_surface.show, "",
            checkable=True, checked=False,
            tip="Draw the membrane dome that the Physics panel measures.\n"
                "TWO surfaces, in two colours. The BLUE cap is a sphere fitted\n"
                "to the transmembrane helices on screen — a measurement, and\n"
                "the one thing that shows whether it was fitted to the right\n"
                "atoms. The PURPLE skirt outside it is the linearised Helfrich\n"
                "footprint: a solution to an equation with two registered\n"
                "parameters in it, which overestimates this footprint 3.65x at\n"
                "PIEZO1's contact slope. Nothing here resolves it.")
    _action(menu, "&Pore surface", win.pore_surface.show, "",
            checkable=True, checked=False,
            tip="Draw the pore as the probe spheres it was measured with: at\n"
                "each height, the largest sphere that fits without touching\n"
                "an atom. Same profile the Analysis panel plots — read, not\n"
                "recomputed, so the picture and the plot cannot be of\n"
                "different runs. Red is narrower than a bare ion, amber\n"
                "clears that but not the hydrated cut, blue clears both, and\n"
                "the residues lining the narrowest slice are marked.\n"
                "These are the spheres that FIT, not the pore wall — and a\n"
                "radius does not settle whether it conducts, because a wide\n"
                "hydrophobic lumen dewets.")
    _action(menu, "P&ockets", win.pocket_view.show, "",
            checkable=True, checked=False,
            tip="Draw the top-ranked cavities as the alpha spheres the\n"
                "detector found them with. The spheres OVERLAP heavily, so\n"
                "their count is the sampling and the reported volume is\n"
                "their union rather than their sum. Ligands are excluded\n"
                "before detection — a bound lipid fills the very pocket\n"
                "being looked for — so a drawn pocket may sit on top of one.\n"
                "A cavity is geometry, not a binding site: no deposited\n"
                "PIEZO entry holds a bound modulator.")
    _action(menu, "&Allosteric path", win.path.show, "",
            checkable=True, checked=False,
            tip="Draw the cheapest route from the blade to the hydrophobic\n"
                "gate through a contact graph weighted by -log|DCC| — the\n"
                "picture of this project's central mechanical claim. Needs\n"
                "normal modes first. The tube is coloured by each step's own\n"
                "correlation, so the weakest link is visible rather than\n"
                "averaged away. A drawn line reads as unique and is not, so\n"
                "the search is re-run with this route's edges removed and\n"
                "the status line gives what the best alternative costs.")
    _action(menu, "&Calcium nanodomain", win.nanodomain.show, "",
            checkable=True, checked=False,
            tip="Draw the calcium field an open channel makes around its own\n"
                "pore exit, as the two surfaces that matter: where the sensor\n"
                "is 90% occupied and where it is half-occupied at its Kd.\n"
                "The Green's function is spherically symmetric, so these are\n"
                "exactly spheres — nothing is idealised for drawing. The\n"
                "model is a point source in free solution and does not know\n"
                "the protein is there. A shut structure draws NOTHING and\n"
                "says why: no current is borrowed from another entry.")
    _action(menu, "&Full-length model", win.hybrid.show, "",
            checkable=True, checked=False,
            tip="Graft the AlphaFold distal blade onto the experimental core.\n"
                "The experiment resolves ~570-2521; the other 569 residues are\n"
                "a PREDICTION and are coloured by pLDDT so they cannot be read\n"
                "as measured. The seam is marked, and the status line gives the\n"
                "75 A by which the two models disagree away from it.")

    from .menus_styles import build_halotag_menu, build_hybrid_style_menu
    build_hybrid_style_menu(win, menu, _action)
    build_halotag_menu(win, menu, _action)

    menu.addSeparator()

    _action(menu, "Remove e&xtra structures", win.clear_companions, "",
            tip="Drop everything except the primary structure. Whether "
                "loading keeps or replaces what is on screen is an option: "
                "Options → Show multiple structures at once.")

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
    _action(menu, "&Fluctuation vs B-factor…", win.show_fluctuations, "",
            "The standard validation of an elastic network: does the predicted\n"
            "mean-square fluctuation track the deposited B-factor? Reported\n"
            "with a contact-number control that uses no network at all, and\n"
            "with the quality of the column itself, because a grouped or\n"
            "predicted-model column cannot answer the question.")
    _action(menu, "PIEZO&2 comparison…", win.show_paralogue, "",
            "The paralogue control. Dome and gating-mode symmetry,\n"
            "coverage-matched through a real alignment: 6KG7 resolves all\n"
            "38 TM helices where a PIEZO1 entry resolves 22-26.")
    _action(menu, "PIEZO &family comparison…", win.show_homology, "",
            "The same question across the catalogue — PIEZO2, PEZO-1, dPIEZO.\n"
            "Reports a RANGE: one entry pair gives 0.98 and another 0.19.")
    menu.addSeparator()
    _action(menu, "The PIEZO family &census…", win.show_family, "",
            "Everything imported from the piezo_genes census: 13 findings,\n"
            "each with the number it rests on, the file it came from and\n"
            "what this project does with it. NOTHING here is measured on the\n"
            "loaded structure - the four entries below are the ones that\n"
            "measure.")
    _action(menu, "Evolutionary &constraint…", win.show_constraint, "",
            "Per-residue constraint over 174 PIEZO1 orthologues, summarised\n"
            "on THIS project's domain boundaries rather than the census's.\n"
            "The pore machinery comes out most constrained and the blades\n"
            "least - and one census finding, the distal-versus-proximal\n"
            "blade gradient, REVERSES on our boundaries because its bands\n"
            "differ in how much inter-unit linker they contain.")
    _action(menu, "Where &disease sits…", win.show_disease_geography, "",
            "Does pathogenic missense concentrate in the pore module?\n"
            "Re-tested on our own variants against gnomAD population\n"
            "variation. Reported under BOTH domain partitions, because the\n"
            "two disagree by 120 residues and the answer follows.")
    _action(menu, "Core and peri&phery…", win.show_core_periphery, "",
            "Superpose a partner on this entry by the pore module alone,\n"
            "then measure where the blades land. Experimental cross-paralogue\n"
            "pairs splay 0.8-2.5x; an AlphaFold monomer of the SAME protein\n"
            "splays 7.2-9.1x, which is what makes the control worth having.")
    _action(menu, "piezo&3 - the third PIEZO…", win.show_piezo3, "",
            "The paralogue vertebrates have and humans lost: transcribed,\n"
            "spliced and under selection at its pore, with the identical\n"
            "residue at all 14 pathogenic pore positions. Assembled into a\n"
            "trimer on a deposited template and run through the pipeline -\n"
            "with 96% of the resulting shape borrowed from that template.")
    menu.addSeparator()
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
    _action(menu, "&Liu 2025 figures…", win.show_liu2025, "",
            tip="Every panel of the paper the intermediate-open structure\n"
                "(8IXO) comes from, and what this project can do with each.\n"
                "6 reproduce, 7 have an ANALOGUE that is a different\n"
                "quantity, 11 need patch clamp, a cryo-EM map or a\n"
                "molecular-dynamics trajectory we do not hold.\n"
                "The curvature panel DISAGREES and says so.")
    _action(menu, "&Guo && MacKinnon 2017 figures…", win.show_guo2017, "",
            "Replicate the paper the dome model comes from, panel by panel.\n"
            "Sixteen of its thirty-one panels reproduce from coordinates,\n"
            "including every number in Figure 7 and its supplement. Three\n"
            "have only an analogue that is a different quantity, and twelve\n"
            "need the cryo-EM map or the micrographs; each says which.")
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
    _action(menu, "Full-length model &numbers…", win.show_hybrid, "",
            "How much of the full-length model is predicted rather than\n"
            "measured: the grafted range, the fraction clearing pLDDT 70, and\n"
            "the 75 A by which the two models disagree away from the seam.\n"
            "Draw it with View > Full-length model.")
    _action(menu, "HaloTag &geometry…", win.show_fusion_numbers, "",
            "Where a C-terminal HaloTag would sit: accessible volume, distance\n"
            "to the pore exit and clearance. Draw it with View > HaloTag fusion.")

    menu.addSeparator()
    _action(menu, "&Topology diagram…", win._show_topology, "Ctrl+Shift+T",
            "One protomer's membrane topology, after Guo & MacKinnon 2017\n"
            "Figure 3: 38 transmembrane helices in nine 4-TM units, the cap,\n"
            "the beam and the cuff. Helices this entry does not model are\n"
            "drawn dashed rather than dropped, so the numbering cannot shift.\n"
            "Tick a unit to box it as Figure 3b does and select it on the\n"
            "model; shift-click a helix to do the same in one step.")
    _action(menu, "&Sequences…", win._show_sequences, "Ctrl+Shift+S",
            "Browse the protein and coding sequences, select onto the model, "
            "and compare sequences with alignment options")


def _options_menu(win, bar) -> None:
    """Every preference the application remembers, in one menu.

    Built in `menus_options.py`: the rule (Options holds what is remembered
    across sessions, View holds what is shown right now) and the two
    deliberate exceptions are documented there.
    """
    from .menus_options import build_options_menu

    build_options_menu(win, bar, _action)


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
