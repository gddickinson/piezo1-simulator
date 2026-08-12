"""Rendering-style choices for the modelled overlays and extra structures.

Split from :mod:`piezo1.ui.menus` at the project's 500-line limit and along a
real seam: everything here chooses **how a feature is drawn**, never what is
computed or which caveat accompanies it. The primary structure has had a style
selector since the beginning (the Model panel) and the superposition overlay
since Round 44 (the Overlay panel); the HaloTag fold, the full-length graft,
the extra structures and the component highlight were each drawn in exactly
one hard-coded representation, reachable by nobody. These submenus close that
gap.

The rule every entry obeys: restyling moves no caveat. The fold keeps its
UNDETERMINED status line in every style, the graft keeps its pLDDT bands and
seam marker, the highlight keeps its gold. A style that could soften what the
picture admits about itself would not belong in a menu.
"""

from __future__ import annotations

from PyQt6.QtGui import QActionGroup

__all__ = ["build_halotag_menu", "build_companion_style_menu",
           "build_hybrid_style_menu", "build_component_style_menu"]


def _radio_group(win, menu, _action, choices, current, apply, tip) -> None:
    """One exclusive group of style actions, sharing a tooltip."""
    group = QActionGroup(menu)
    group.setExclusive(True)
    for key, label in choices:
        action = _action(
            menu, label,
            (lambda checked, k=key: checked and apply(k)),
            "", checkable=True, checked=(key == current), tip=tip)
        group.addAction(action)


def build_halotag_menu(win, parent, _action) -> None:
    """The HaloTag fusion submenu, with the fold's representation choices.

    Moved here from `menus.py` when the style group joined it: the submenu and
    the styles belong together, and `menus.py` was at the length limit.
    """
    halotag = parent.addMenu("&HaloTag fusion")
    halotag.setToolTipsVisible(True)
    _action(halotag, "Show modelled &tags", win.fusion.show, "",
            checkable=True, checked=False,
            tip="Draw a HaloTag at each of the three cytosolic C-termini.\n"
                "THIS IS A MODEL: there is no structure of the fusion, so the\n"
                "tag body is drawn as a sphere of its radius of gyration and\n"
                "the linker as a straight seam.")
    _action(halotag, "Show tag &structure", win.fusion.set_atoms, "",
            checkable=True, checked=False,
            tip="Draw the tag's real fold — the deposited 6U32 coordinates —\n"
                "instead of the sphere, placed so its N-terminus faces the\n"
                "channel's C-terminus. The POSITION is the model's; the SPIN\n"
                "about the linker is undetermined, so this is one draw of\n"
                "many. Atoms inside the channel are red.")
    _action(halotag, "T&urn tag orientation", win.fusion.rotate_tags, "",
            tip="Rotate the fold about the linker by 10°. Nothing else moves:\n"
                "the free angle is shown rather than asserted, because a\n"
                "drawn fold otherwise reads as a determined pose.")

    from .fusion_controller import FOLD_STYLES

    styles = halotag.addMenu("Tag structure st&yle")
    styles.setToolTipsVisible(True)
    _radio_group(
        win, styles, _action, FOLD_STYLES, win.fusion.fold_style,
        win.fusion.set_fold_style,
        tip="How the fold is drawn once it is shown at all — the same\n"
            "representations the channel offers. PRESENTATION ONLY: every\n"
            "style is the same rigidly placed 6U32 at the same undetermined\n"
            "spin, the contact atoms stay red, the dye keeps its colour and\n"
            "the status line keeps its caveat. A cartoon of the fold is no\n"
            "more a determined pose than a sphere cloud of it.")

    _action(halotag, "Show accessible &volume", win.fusion.set_envelope, "",
            checkable=True, checked=False,
            tip="The region the tag centre can occupy without clashing, as a\n"
                "point cloud. Shown so a single sphere is not mistaken for a\n"
                "determined position.")
    _action(halotag, "Show &dyes", win.fusion.set_dyes, "",
            checkable=True, checked=False,
            tip="Draw a dye on each tag the labelling model says is occupied.")


def build_hybrid_style_menu(win, parent, _action) -> None:
    """How the full-length model is drawn, next to the action that draws it."""
    from .hybrid_controller import HYBRID_STYLES

    menu = parent.addMenu("Full-length model st&yle")
    menu.setToolTipsVisible(True)
    _radio_group(
        win, menu, _action, HYBRID_STYLES, win.hybrid.style,
        win.hybrid.set_style,
        tip="How the full-length model is drawn. Every style keeps the\n"
            "grey-versus-pLDDT colouring, the seam marker and the status\n"
            "line — the things that stop a complete-looking trimer from\n"
            "reading as measured. The ribbon styles trace the C-alphas, so\n"
            "the chain's path shows where the sphere cloud shows bulk.")


def build_companion_style_menu(win, parent, _action) -> None:
    """How extra structures are drawn, next to the toggle that shows them."""
    from .panels.structure_panel import STYLE_LABELS

    menu = parent.addMenu("Extra structures st&yle")
    menu.setToolTipsVisible(True)
    current = win.companion_style()
    choices = [(style.value, label) for label, style in STYLE_LABELS]
    _radio_group(
        win, menu, _action, choices, current.value, win.set_companion_style,
        tip="How the extra structures are drawn when several are shown at\n"
            "once. One choice for all of them: companions are told apart by\n"
            "COLOUR, and a mixture of styles would hand that job to shape\n"
            "as well. Backbone by default so the primary stays the thing\n"
            "being looked at. Analyses run on the primary regardless.")


def build_component_style_menu(win, parent, _action) -> None:
    """How the component's curated residues are highlighted."""
    from .component_controller import HIGHLIGHT_STYLES

    menu = parent.addMenu("Highlighted residues st&yle")
    menu.setToolTipsVisible(True)
    _radio_group(
        win, menu, _action, HIGHLIGHT_STYLES,
        win.components.highlight_style, win.components.set_style,
        tip="How the curated residues of the chosen component are drawn on\n"
            "top of its backbone. All three styles stay gold — the colour is\n"
            "what says these are annotation rather than structure — and none\n"
            "changes WHICH residues are highlighted, which comes from the\n"
            "curated groups alone.")
