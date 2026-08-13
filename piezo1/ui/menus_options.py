"""The Options menu: every preference the application remembers.

Split from :mod:`piezo1.ui.menus` when the preferences were consolidated, and
along a rule rather than a length limit: **Options holds what is remembered
across sessions; View holds what is shown right now.** Before this, four
persisted preferences lived under View — the structure alignment mode, the
multi-structure toggle, the companion style and the display-options dialog —
so the answer to "where do I change how this opens next time" depended on
which menu a feature had happened to be added to.

Two things that look like options deliberately stay where they are. The ion
flux pathway and voltage sit beside the animation they parameterise, because
neither is a preference — each changes the physics of the number on the
status line (see `menus_flux.py`). And the per-feature style submenus stay
beside their feature toggles: how the HaloTag fold is drawn is a choice about
the fold, made while looking at it.
"""

from __future__ import annotations

from PyQt6.QtGui import QActionGroup

__all__ = ["build_options_menu"]


def _radio(win, menu, _action, choices, current, apply, tip) -> None:
    group = QActionGroup(menu)
    group.setExclusive(True)
    for key, label in choices:
        action = _action(
            menu, label,
            (lambda checked, k=key: checked and apply(k)),
            "", checkable=True, checked=(key == current), tip=tip)
        group.addAction(action)


def build_options_menu(win, bar, _action) -> None:
    menu = bar.addMenu("&Options")
    menu.setToolTipsVisible(True)
    settings = win.settings

    # ---------------------------------------------------------- appearance
    from .preferences import BACKGROUNDS
    from .theme import THEMES

    theme = menu.addMenu("Interface &theme")
    theme.setToolTipsVisible(True)
    _radio(win, theme, _action, THEMES, win.ui_theme(), win._set_ui_theme,
           tip="The application chrome: panels, menus, buttons. Dark is the\n"
               "default — depth cueing and rim lighting read best beside it.\n"
               "Light suits a bright office and manuscript work; System\n"
               "hands the styling back to the platform. The 3-D viewport's\n"
               "own background is the option next to this one.")

    background = menu.addMenu("Viewport &background")
    background.setToolTipsVisible(True)
    _radio(win, background, _action,
           [(key, label) for key, label, _c in BACKGROUNDS],
           win.background_key(), win._set_background,
           tip="The colour behind the structure — also the colour the\n"
               "depth-cue fog fades into, so distant atoms recede into the\n"
               "background rather than into a haze of a different colour.\n"
               "White and light grey are for figures on white pages. The\n"
               "scale bar and readouts keep their dark halo and stay\n"
               "legible on every choice.")

    _action(menu, "&Display options…", win._show_display_options, "Ctrl+D",
            "Choose what the overlay shows: scale bar, animation clock, "
            "orientation axes and which measured values")

    menu.addSeparator()

    # ------------------------------------------------------------ behaviour
    _action(menu, "Remember &layout on exit",
            win._set_remember_layout, "",
            "Reopen with the panel arrangement and window size you left",
            checkable=True,
            checked=settings.value("options/remember_layout", True, type=bool))
    _action(menu, "Show status-bar &hints", win._set_show_hints, "",
            "The mouse and keyboard reminder in the status bar",
            checkable=True,
            checked=settings.value("options/show_hints", True, type=bool))

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

    _action(menu, "Show &multiple structures at once", win.set_multi_structure,
            checkable=True, checked=win.multi_structure,
            tip="Keep the current structure on screen when another is loaded, "
                "drawn in its own colour in the same frame. Off by default: "
                "two entries in the same frame sit on top of each other, and a "
                "structure left behind reads as extra density. Analyses always "
                "run on the primary structure, whatever else is drawn.")
    from .menus_styles import build_companion_style_menu
    build_companion_style_menu(win, menu, _action)

    menu.addSeparator()
    _action(menu, "&Parameters…", win._show_parameters, "Ctrl+P",
            "Every number the calculations use, with its default, its unit "
            "and the paper it came from. Editable.")
    menu.addSeparator()
    _action(menu, "&Restore default options", win._reset_options, "",
            "Forget every remembered setting and layout, including the "
            "theme and the viewport background")
