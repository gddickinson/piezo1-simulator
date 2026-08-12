"""Loading an entry this project has no annotation for must not crash, or lie.

``domains.json``, ``variants.json`` and ``functional_residues.json`` are curated
in human and mouse PIEZO1 numbering and in nothing else. The catalogue holds
eight numberings. That gap had two faces and both were live before Round 89b:

* **It crashed.** ``DomainPalette`` indexed ``d[self.species]`` straight into
  the domain records, so loading any PIEZO2, PEZO-1 or dPIEZO entry raised
  ``KeyError`` at the moment the renderer coloured by domain — which is the
  default colouring, so the entry could not be opened at all. It surfaced when
  the plant filter was added because that put a new numbering one click away,
  but 6KG7 had been in the catalogue since Round 83 and behaved identically.

* **It lied, silently, which is worse.** ``Annotations`` built each functional
  residue group with ``r["human"] if species == "human" else r["mouse"]``, so
  **every non-human numbering got mouse PIEZO1 residue numbers**. Selecting the
  hydrophobic gate on a *C. elegans* structure highlighted whatever sits at
  mouse 2473/2476/2480. Nothing raised, the picture looked right, and the
  residues were a different protein's.

The rule is now one line in ``core.annotations``: an unannotated numbering
yields **nothing**, and the caller is told why.
"""

from __future__ import annotations

import pytest

from piezo1.core.annotations import (ANNOTATED_NUMBERINGS, annotation_gap,
                                     is_annotated, load_annotations)
from piezo1.core.structure import Structure
from piezo1.io.registry import load_registry
from piezo1.render import colormaps


def _one_entry_per_numbering():
    seen = {}
    for record in load_registry().available():
        seen.setdefault(record.numbering_species, record)
    if not seen:
        pytest.skip("no structures downloaded — run python -m piezo1.io.fetch")
    return seen


def test_every_catalogued_numbering_can_be_coloured_by_domain():
    """The crash, driven over the real path rather than asserted about.

    ``domain_colors`` is what ``MolecularView.rebuild`` calls, and rebuild is
    what ``load_structure`` calls. If this passes, the entry opens.
    """
    for numbering, record in _one_entry_per_numbering().items():
        structure = Structure.from_file(record.path)
        palette = colormaps.load_domain_palette(numbering)
        colors = colormaps.domain_colors(structure, palette)
        assert colors.shape == (structure.n_atoms, 3), f"{record.pdb} {numbering}"


def test_the_catalogue_really_does_span_numberings_without_annotation():
    """Otherwise the test above passes for want of a case that exercises it."""
    numberings = set(_one_entry_per_numbering())
    assert numberings & set(ANNOTATED_NUMBERINGS), "no annotated entry"
    unannotated = {n for n in numberings if not is_annotated(n)}
    assert len(unannotated) >= 4, unannotated


def test_an_unannotated_numbering_yields_nothing_rather_than_mouse():
    """The silent one. A PEZO-1 gate must be empty, never mouse PIEZO1's.

    Checked against the values it used to return, so the test fails if the
    fallback ever comes back rather than merely if the lists are non-empty.
    """
    mouse = load_annotations("mouse")
    mouse_gate = next(g for g in mouse.residue_groups
                      if g.id == "hydrophobic_gate")
    assert mouse_gate.residues, "the mouse gate is the thing that leaked"

    for numbering in ("worm_piezo", "fly_piezo", "human_piezo2", "plant_piezo",
                      "dicty_piezo", "rat"):
        annotations = load_annotations(numbering)
        assert annotations.residue_groups == [], numbering
        assert annotations.domains == [], numbering
        assert annotations.variants == [], numbering
        # And specifically not the mouse numbers, which is what it did return.
        found = {r for g in annotations.residue_groups for r in g.residues}
        assert not (found & set(mouse_gate.residues)), numbering


def test_the_gap_is_stated_rather_than_left_as_an_empty_list():
    """An empty domain list drawn as uniform grey reads as 'no domains here'.

    That is the opposite of what is true, so the reason has to travel with the
    emptiness — on the object, and from there onto the status line.
    """
    for numbering in ANNOTATED_NUMBERINGS:
        assert annotation_gap(numbering) == ""
        assert not load_annotations(numbering).gap

    gap = annotation_gap("worm_piezo")
    assert gap and "worm_piezo" in gap
    assert "human and mouse" in gap
    assert load_annotations("worm_piezo").meta["unannotated"] == gap


def test_the_status_line_says_so_when_an_unannotated_entry_is_loaded():
    """The guard that keeps the reason from being computed and then dropped."""
    import inspect

    from piezo1.ui import main_window

    source = inspect.getsource(main_window.MainWindow.load_structure)
    assert "annotation_gap(" in source
    assert "NO ANNOTATION FOR THIS PROTEIN" in source


def test_annotated_numberings_are_the_two_the_resources_actually_carry():
    """Read from the files, so adding a curated numbering updates the constant."""
    import json

    from piezo1.config import RESOURCE_DIR

    domains = json.loads((RESOURCE_DIR / "domains.json").read_text())["domains"]
    keyed = {k for d in domains for k in d
             if isinstance(d.get(k), dict) and "start" in d[k]}
    assert keyed == set(ANNOTATED_NUMBERINGS), keyed


def test_piezo1_numbering_gates_on_annotation_not_on_being_piezo1():
    """Rat is PIEZO1 and has no annotation, and those are different questions.

    Round 89 added rat as a third PIEZO1 reference — right, because a rat entry
    should be identified as rat rather than mis-read as mouse — and that
    silently widened ``piezo1_numbering``'s return set to a numbering no
    resource carries. Every caller takes the string straight to
    ``load_annotations``, so the component selector, the conduction path and
    the pore charge map would have got an empty result where they should get a
    refusal.
    """
    from piezo1.core.numbering_check import PIEZO1_REFERENCES, piezo1_numbering

    assert "rat" in PIEZO1_REFERENCES, "the premise: rat is a PIEZO1 reference"
    assert not is_annotated("rat")

    records = {r.pdb: r for r in load_registry().available()}
    rat = records.get("AF-Q0KL00-F1-MODEL_V6")
    if rat is None:
        pytest.skip("rat model not downloaded")
    assert piezo1_numbering(Structure.from_file(rat.path)) is None

    human = records.get("8YEZ")
    if human is not None:
        assert piezo1_numbering(Structure.from_file(human.path)) == "human"


def test_the_component_selector_refuses_an_unannotated_entry():
    """It must say the numbering is unreadable, not report missing domains."""
    from piezo1.structure.components import COMPONENTS, component_masks

    records = {r.pdb: r for r in load_registry().available()}
    rat = records.get("AF-Q0KL00-F1-MODEL_V6")
    if rat is None:
        pytest.skip("rat model not downloaded")
    selection = component_masks(Structure.from_file(rat.path), "pore_module")
    assert selection.numbering == ""
    assert "numbering not readable" in selection.note
