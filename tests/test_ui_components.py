"""Showing one part of the assembly, and the three ways that could mislead.

1. **A component could hide atoms from the analyses**, not just from the eye.
   It must not: the pore profile, the dome and the modes all still run on the
   whole trimer, and the status line says so on every switch.
2. **It could hide only some representations.** The cartoon traces are built
   once when the view is made, so setting the residue filter and calling
   `rebuild()` filtered the atoms and left the ribbon drawing the whole chain —
   at the default style that is most of what is on screen, so hiding 97% of the
   atoms changed the picture by a tenth. Measured in pixels here, because that
   is the only assertion that would have caught it.
3. **It could select the wrong residues.** Every range comes from curated
   annotation read in the entry's own numbering, and a structure whose
   numbering cannot be read falls back to the whole assembly with the reason
   rather than selecting plausible-looking nonsense.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from piezo1.config import RenderSettings, STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402
from piezo1.structure.components import (COMPONENTS, component_by_key,  # noqa: E402
                                         component_masks)

SIZE = (320, 240)


def _load(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


@pytest.fixture(scope="module")
def structure():
    return _load("8IXO")


# ------------------------------------------------------- the selection itself

def test_the_whole_assembly_selects_everything(structure):
    selection = component_masks(structure, "whole")
    assert selection.backbone.all()
    assert not selection.highlight.any(), (
        "the default view must not pick anything out")


@pytest.mark.parametrize("component", [c for c in COMPONENTS if not c.is_whole])
def test_every_component_selects_a_real_subset(structure, component):
    selection = component_masks(structure, component.key)
    assert 0 < selection.n_atoms < structure.n_atoms, component.key
    assert selection.numbering == "mouse", "8IXO is a mouse entry"


@pytest.mark.parametrize("component", [c for c in COMPONENTS
                                       if c.highlight and not c.is_whole])
def test_highlighted_residues_are_always_drawn(structure, component):
    """A residue picked out as important that is then not drawn is the worst
    of both. The PIP2 lysines sit just outside THU9's range, which is how this
    came up."""
    selection = component_masks(structure, component.key)
    assert selection.highlight.any(), component.key
    assert (selection.backbone | selection.highlight).sum() == selection.n_atoms


def test_the_pore_module_contains_all_four_gates(structure):
    """The Figure 2E view is only that view if every gate is in it."""
    selection = component_masks(structure, "pore_module")
    from piezo1.core.annotations import load_annotations

    annotations = load_annotations("mouse")
    drawn = set(selection.residues)
    for group_id in ("cap_constriction", "cap_gate", "hydrophobic_gate",
                     "ctd_constriction"):
        residues = set(annotations.group(group_id).residues)
        assert residues & drawn, f"{group_id} is not in the pore module"


def test_the_pore_module_is_a_small_fraction_of_the_trimer(structure):
    """The reason to have it at all."""
    selection = component_masks(structure, "pore_module")
    assert selection.n_atoms < 0.4 * structure.n_atoms


def test_an_unreadable_numbering_falls_back_rather_than_guessing():
    """Every range is a residue number; PIEZO2's gate is not at PIEZO1's.

    A confident selection of the wrong helices looks exactly like a right one,
    so the refusal shows the whole assembly and records why.
    """
    selection = component_masks(_load("6KG7"), "pore_module")
    assert selection.backbone.all(), "should fall back to everything"
    assert not selection.highlight.any()
    assert "numbering" in selection.note


def test_an_unknown_component_is_refused():
    with pytest.raises(KeyError, match="no component"):
        component_by_key("nonsense")


def test_every_component_states_what_it_shows():
    for component in COMPONENTS:
        assert len(component.shows) > 40, component.key


# ------------------------------------------------ it hides, it does not subset

def test_selecting_a_component_does_not_change_what_is_measured(structure):
    """The rule this project runs on: drawing never decides computing.

    The dome is measured before and after a component is selected and must be
    bit-identical — the controller changes the *view*, and the analyses read
    the structure.
    """
    from piezo1.structure.geometry import measure_dome, tm_surface_points
    from piezo1.structure.protomers import protomer_blocks

    points, _ = tm_surface_points(structure, "mouse")
    blocks, _ = protomer_blocks(structure)
    before = measure_dome(blocks, points).radius_of_curvature

    component_masks(structure, "tm_gate")          # the narrowest selection
    after = measure_dome(blocks, points).radius_of_curvature
    assert before == after


def test_the_status_line_says_it_is_hidden_not_removed():
    from piezo1.ui.component_controller import ComponentController
    import inspect

    source = inspect.getsource(ComponentController._announce)
    assert "hidden, not removed" in source


# ------------------------------------------------------ drawn, counted in pixels

@pytest.fixture(scope="module")
def gl():
    moderngl = pytest.importorskip("moderngl")
    try:
        return moderngl.create_standalone_context(require=410)
    except Exception as exc:                              # pragma: no cover
        pytest.skip(f"no OpenGL 4.1 context available: {exc}")


def _scene(gl):
    from piezo1.render.scene import Scene

    scene = Scene(gl, RenderSettings(samples=1))
    scene.resize(*SIZE)
    return scene


def _lit(gl, scene):
    fbo = gl.simple_framebuffer(SIZE)
    fbo.use()
    fbo.clear(0.05, 0.05, 0.07, 1.0)
    scene.render()
    pixels = np.frombuffer(fbo.read(components=3), np.uint8)
    pixels = pixels.reshape(-1, 3).astype(int)
    on = pixels.sum(axis=1) > 60
    return int(on.sum()), (float(pixels[on].mean()) if on.any() else 0.0)


def test_hiding_a_component_really_hides_the_cartoon(gl, structure):
    """The defect this test exists for, measured rather than asserted.

    `traces` is built once in `__post_init__`. Setting `visible_residues` and
    calling `rebuild()` filtered every atom representation and left the ribbon
    untouched — and the ribbon is the default style, so the picture barely
    changed. `set_visible_residues` rebuilds the traces too.
    """
    from piezo1.render.representations import MolecularView, Style

    scene = _scene(gl)
    view = MolecularView(scene, structure, name="model")
    view.style = Style.CARTOON
    view.set_species("mouse")
    view.rebuild()
    scene.camera.frame(structure.xyz.astype(np.float32))
    whole, _ = _lit(gl, scene)

    selection = component_masks(structure, "tm_gate")
    view.set_visible_residues(selection.residues)
    part, _ = _lit(gl, scene)

    assert part < 0.3 * whole, (
        f"the transmembrane gate is 3% of the atoms and drew {part} pixels "
        f"against the whole assembly's {whole} — the cartoon is not being "
        f"filtered")

    view.set_visible_residues(None)
    assert _lit(gl, scene)[0] == whole, "restoring the whole assembly failed"


def test_the_pore_opacity_dims_without_changing_the_geometry(gl):
    """Opacity must blend, not shrink.

    A guard that only counted lit pixels would pass on a batch that had simply
    stopped drawing, so brightness is the assertion and the pixel count is the
    control.
    """
    scene = _scene(gl)
    batch = scene.spheres("pore:probe")
    batch.upload(np.zeros((1, 3), np.float32), np.array([6.0], np.float32),
                 np.array([[0.9, 0.3, 0.3]], np.float32))
    scene.camera.frame(np.array([[-12, -12, -12], [12, 12, 12]], np.float32))

    assert not batch.transparent
    opaque_count, opaque_brightness = _lit(gl, scene)

    batch.alpha = 0.25
    assert batch.transparent, (
        "a batch below full opacity must move into the blended pass, or it "
        "writes depth and hides what it was meant to reveal")
    faded_count, faded_brightness = _lit(gl, scene)

    assert faded_count == opaque_count, "opacity changed the geometry"
    assert faded_brightness < 0.75 * opaque_brightness, (
        f"brightness went {opaque_brightness:.0f} -> {faded_brightness:.0f}; "
        f"the alpha is not reaching the shader")


def test_the_pore_controller_exposes_the_opacity():
    from piezo1.ui.pore_controller import PoreSurfaceController

    assert hasattr(PoreSurfaceController, "set_opacity")


# ------------------------------------------------------------- the colouring

def test_hydrophobicity_puts_the_scale_extremes_at_the_extremes(structure):
    from piezo1.render.colormaps import hydrophobicity_colors

    colors = hydrophobicity_colors(structure)
    ile = colors[structure.res_name == "ILE"][0]
    arg = colors[structure.res_name == "ARG"][0]
    assert ile[0] > ile[2], "isoleucine is the most apolar residue on the scale"
    assert arg[2] > arg[0], "arginine is the most polar"


def test_hydrophobicity_is_fixed_scale_and_not_auto_ranged(structure):
    """Two different structures must be comparable, and comparable with a
    published panel — which auto-ranging destroys."""
    from piezo1.render.colormaps import hydrophobicity_colors

    subset = structure.subset(structure.res_name == "ALA")
    if subset.n_atoms == 0:
        pytest.skip("no alanine")
    everywhere = hydrophobicity_colors(structure)[structure.res_name == "ALA"][0]
    alone = hydrophobicity_colors(subset)[0]
    assert np.allclose(everywhere, alone), (
        "alanine changed colour when nothing else was on screen")


def test_an_unscored_residue_is_not_drawn_as_neutral(structure):
    """"Not scored" must not read as "neither polar nor apolar"."""
    from piezo1.render.colormaps import hydrophobicity_colors

    unknown = structure.res_name == "UNK"
    if not unknown.any():
        pytest.skip("this entry names every residue")
    colors = hydrophobicity_colors(structure)
    gly = colors[structure.res_name == "GLY"][0]
    assert not np.allclose(colors[unknown][0], gly)


def test_hydrophobicity_is_offered_in_the_gui():
    from piezo1.render.representations import ColorBy
    from piezo1.ui.panels.structure_panel import COLOR_LABELS

    assert ColorBy.HYDROPHOBICITY in {value for _label, value in COLOR_LABELS}
