"""What is drawn, and what the ion animation is a picture of.

Split from :mod:`piezo1.ui.menus` at the project's 500-line limit and along a
real seam: every other entry in that module chooses what is *shown*, and these
two choose what is *computed*. Neither is a preference — the pathway decides
whether any structure conducts at all, and the voltage is the axis Liu et al.
2025's Figure 5 is swept along. Both default to what this project already
recorded: the axial route and the registered test voltage.
"""

from __future__ import annotations

from PyQt6.QtGui import QActionGroup

__all__ = ["build_flux_settings", "build_component_menu", "build_pore_opacity_menu", "OPACITIES"]


def build_flux_settings(win, parent, _action) -> None:
    """The two choices that decide what the ion animation is a picture of.

    Both are opt-in and both default to what the project already recorded:
    the axial pathway and the registered test voltage. Kept in a submenu next
    to the animation rather than in Options, because neither is a preference —
    each changes the physics of the number on the status line.
    """
    from ..physics.conduction_path import PATHWAY_LABELS, PATHWAYS

    route = parent.addMenu("Ion flux &pathway")
    route.setToolTipsVisible(True)
    group = QActionGroup(route)
    group.setExclusive(True)
    for pathway in PATHWAYS:
        action = _action(
            route, PATHWAY_LABELS[pathway],
            (lambda checked, p=pathway: checked and win.ion_flux.set_pathway(p)),
            "", checkable=True, checked=(pathway == "axial"),
            tip="Which route the ions may take.\n"
                "AXIAL is bulk-to-bulk down the three-fold axis and is what\n"
                "every recorded number here was computed on. It refuses every\n"
                "deposited entry but two, because PIEZO1's axis is pinched\n"
                "below the water radius at R2295 on top and at the\n"
                "cytoplasmic constriction beneath.\n"
                "The LATERAL options exclude those closed ends, which is the\n"
                "route Liu et al. 2025 report: in through the cap gates, out\n"
                "through the intracellular portals. They do NOT model the\n"
                "portal, so the current they give is an upper bound.")
        group.addAction(action)

    volts = parent.addMenu("Ion flux &voltage")
    volts.setToolTipsVisible(True)
    vgroup = QActionGroup(volts)
    vgroup.setExclusive(True)
    choices = [("Registered default", None)] + [
        (f"{v:+.2f} V (Liu et al. Figure 5)", v)
        for v in (0.0, -0.1, -0.25, -0.5)]
    for label, value in choices:
        action = _action(
            volts, label,
            (lambda checked, v=value: checked and win.ion_flux.set_voltage(v)),
            "", checkable=True, checked=(value is None),
            tip="The transmembrane potential the current is computed at.\n"
                "The four negative values are the ones their Figure 5A\n"
                "sweeps. At 0 V the current is zero and nothing animates,\n"
                "which is the correct picture and not a failure.")
        vgroup.addAction(action)


#: Opacity steps offered for the drawn pore. Coarse on purpose — a slider
#: invites fiddling with a number that carries no information, and the useful
#: question is only whether the lining is visible through the probe spheres.
OPACITIES = (("Opaque", 1.0), ("75%", 0.75), ("50%", 0.5), ("25%", 0.25))


def build_component_menu(win, parent, _action) -> None:
    """Show one named part of the assembly instead of the whole propeller.

    A trimer fills the viewport with blades and almost every question is about
    one part of it. Each entry hides the rest and draws the curated residues of
    that part as ball-and-stick on top — and **hides rather than subsets**, so
    every analysis still runs on the whole assembly. The status line says so on
    every switch, because a selector that silently changed what the dome fit
    measures would be the most confusing possible way to break this.
    """
    from ..structure.components import COMPONENTS

    menu = parent.addMenu("Show &component")
    menu.setToolTipsVisible(True)
    group = QActionGroup(menu)
    group.setExclusive(True)
    for component in COMPONENTS:
        tip = component.shows
        if component.caveat:
            tip += "\n\nCAVEAT: " + component.caveat
        action = _action(
            menu, component.label,
            (lambda checked, k=component.key: checked and win.components.show(k)),
            "", checkable=True, checked=(component.key == "whole"), tip=tip)
        group.addAction(action)


def build_pore_opacity_menu(win, parent, _action) -> None:
    """How solid the drawn pore is.

    The probe spheres are the space left over, and at full opacity they hide
    the residues that bound them — which are the reason anyone looks. Below 1.0
    the batch moves into the scene's blended pass with depth writes off.
    """
    menu = parent.addMenu("Pore surface &opacity")
    menu.setToolTipsVisible(True)
    group = QActionGroup(menu)
    group.setExclusive(True)
    for label, value in OPACITIES:
        action = _action(
            menu, label,
            (lambda checked, v=value: checked and win.pore_surface.set_opacity(v)),
            "", checkable=True, checked=(value == 1.0),
            tip="How solid the drawn probe spheres are.\n"
                "A probe sphere is the space left over, not the wall — at full\n"
                "opacity it hides the lining residues that define it.")
        group.addAction(action)
