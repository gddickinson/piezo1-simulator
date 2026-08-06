"""Pore-radius profiling.

The scientific checks here are that a closed structure comes out occluded, an
open-like structure comes out conductive, and the constrictions the profiler
finds coincide with the residues independently curated as the gate and the
cytoplasmic constrictions.
"""

import numpy as np
import pytest

from piezo1.core.annotations import load_annotations
from piezo1.structure.pore import pore_profile
from piezo1.structure.superpose import detect_c3_axis
from conftest import protomer_blocks


def _axis(st):
    blocks, _ = protomer_blocks(st)
    return detect_c3_axis(blocks)


@pytest.fixture(scope="module")
def human_profile(human_structure):
    return pore_profile(human_structure, _axis(human_structure), step=1.0)


def test_profile_shape_and_positivity(human_profile):
    p = human_profile
    assert len(p.z) == len(p.radius) == len(p.centers) == len(p.slices)
    assert (p.radius >= 0).all()
    # A leashed probe cannot report a bulk-solvent radius. Unleashed, this
    # escapes the protein entirely and returns thousands of Angstrom.
    assert p.radius.max() < 30.0


def test_probe_stays_on_its_leash(human_profile):
    p = human_profile
    radial = p.axis.radial(p.centers)
    assert radial.max() <= p.meta["leash"] + 1e-6


def test_closed_structure_is_occluded(human_profile):
    """8YEZ is a closed human PIEZO1: nothing should get through."""
    assert human_profile.bottleneck_radius < 1.6
    assert human_profile.is_conductive() is False


def test_flat_structure_is_more_open(human_structure, flat_structure):
    """The flattened state must be measurably wider than the curved one.

    This is the structural correlate of gating, so if it ever reverses,
    something is badly wrong.
    """
    curved = pore_profile(human_structure, _axis(human_structure), step=1.5)
    flat = pore_profile(flat_structure, _axis(flat_structure), step=1.5)
    assert flat.bottleneck_radius > curved.bottleneck_radius


def test_constrictions_match_the_curated_functional_residues(human_profile):
    """The profiler should rediscover the annotated gate and CTD constrictions.

    Neither the profiler nor the annotation knows about the other: the residue
    list comes from the literature, the profile from coordinates alone.
    """
    ann = load_annotations("human")
    gate = set(ann.group("hydrophobic_gate").residues)
    ctd = set(ann.group("ctd_constriction").residues)

    lining = set()
    for s in human_profile.constrictions(threshold=3.2):
        lining |= set(s.lining)
    # Allow +/-2 residues, since the probe touches side chains of neighbours.
    def near(target):
        return any(abs(l - t) <= 2 for l in lining for t in target)

    assert near(gate), f"no constriction near the hydrophobic gate; found {sorted(lining)}"
    assert near(ctd), f"no constriction near the CTD constrictions; found {sorted(lining)}"


def test_bottleneck_lining_is_reported(human_profile):
    assert human_profile.bottleneck_lining()
    assert all(isinstance(r, int) for r in human_profile.bottleneck_lining())


def test_ligands_are_excluded_by_default(human_structure):
    """A lipid in the vestibule must not be reported as a protein constriction.

    The two profiles are given an explicit, identical z range. Left to choose
    their own, they pick it from the atoms near the axis — which differ between
    the two selections — and the bottleneck is a sharp minimum, so even a
    fraction of an Angstrom of grid offset changes the reported value. That is
    a genuine sampling subtlety, not a bug, but it makes an automatic range
    the wrong basis for a comparison.
    """
    ax = _axis(human_structure)
    kw = dict(step=1.0, z_min=-60.0, z_max=60.0)
    with_lig = pore_profile(human_structure, ax, protein_only=False, **kw)
    without = pore_profile(human_structure, ax, protein_only=True, **kw)

    assert without.meta["protein_only"] is True
    assert without.meta["n_atoms"] < with_lig.meta["n_atoms"]
    assert np.allclose(with_lig.z, without.z)
    # On a shared grid, extra atoms can only ever reduce the clearance.
    assert (with_lig.radius <= without.radius + 1e-6).all()
