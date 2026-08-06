"""Geometric and surface measurements."""

import numpy as np
import pytest

from piezo1.analysis.measure import (KYTE_DOOLITTLE, angle, buried_area,
                                     crossing_angle, dihedral, distance,
                                     helix_axis, hydrophobicity_profile,
                                     principal_axes, radius_of_gyration, sasa,
                                     tilt_angle)
from piezo1.structure.pore import pore_profile
from piezo1.structure.superpose import detect_c3_axis
from conftest import protomer_blocks


def test_distance_angle_dihedral_on_known_geometry():
    assert distance([0, 0, 0], [3, 4, 0]) == pytest.approx(5.0)
    assert angle([1, 0, 0], [0, 0, 0], [0, 1, 0]) == pytest.approx(90.0)
    assert angle([1, 0, 0], [0, 0, 0], [-1, 0, 0]) == pytest.approx(180.0)
    # A planar cis arrangement is 0 degrees, trans is 180.
    assert abs(dihedral([1, 1, 0], [0, 1, 0], [0, 0, 0], [1, 0, 0])) == pytest.approx(0.0, abs=1e-6)
    assert abs(dihedral([1, 1, 0], [0, 1, 0], [0, 0, 0], [-1, 0, 0])) == pytest.approx(180.0, abs=1e-6)


def test_dihedral_sign_flips_with_chirality():
    a = dihedral([1, 1, 0], [0, 1, 0], [0, 0, 0], [0.5, 0, 0.87])
    b = dihedral([1, 1, 0], [0, 1, 0], [0, 0, 0], [0.5, 0, -0.87])
    assert a * b < 0


def test_radius_of_gyration_of_a_known_shell():
    """For points on a sphere of radius R, Rg == R."""
    rng = np.random.default_rng(0)
    v = rng.normal(size=(1000, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    # Antipodal pairs put the centroid exactly at the origin, so Rg == R
    # exactly rather than approximately.
    v = np.vstack([v, -v])
    assert radius_of_gyration(v * 12.0) == pytest.approx(12.0, rel=1e-9)


def test_principal_axes_of_an_elongated_cloud():
    rng = np.random.default_rng(1)
    pts = rng.normal(size=(500, 3)) * np.array([10.0, 1.0, 1.0])
    vals, vecs = principal_axes(pts)
    assert vals[0] > vals[1] > vals[2]
    assert abs(abs(vecs[0, 0]) - 1.0) < 0.05      # long axis is x


def test_helix_axis_on_an_ideal_alpha_helix():
    """A synthetic ideal helix: 1.5 A rise, 100 degrees per residue, r 2.3 A."""
    n = 18
    t = np.arange(n)
    ca = np.stack([2.3 * np.cos(np.radians(100 * t)),
                   2.3 * np.sin(np.radians(100 * t)),
                   1.5 * t], axis=1)
    axis, centre = helix_axis(ca)
    assert abs(abs(float(np.dot(axis, [0, 0, 1]))) - 1.0) < 1e-3
    assert axis[2] > 0                              # oriented N to C
    assert centre[2] == pytest.approx(1.5 * (n - 1) / 2, abs=0.1)


def test_tilt_and_crossing_angles():
    assert tilt_angle([0, 0, 1], [0, 0, 1]) == pytest.approx(0.0)
    assert tilt_angle([1, 0, 0], [0, 0, 1]) == pytest.approx(90.0)
    # Tilt is direction-agnostic: an antiparallel helix has the same tilt.
    assert tilt_angle([0, 0, -1], [0, 0, 1]) == pytest.approx(0.0)
    assert crossing_angle([1, 0, 0], [0, 1, 0]) == pytest.approx(90.0)


def test_pore_lining_helix_is_the_least_tilted(human_structure):
    """TM38 lines the pore, so it must run nearly parallel to the axis.

    The comparison is restricted to *pore-proximal* helices, within 25 A of the
    symmetry axis. Some blade helices out at 50-60 A are also near-vertical
    (TM30 is 3.1 degrees), so "least tilted overall" would be false; "least
    tilted of those forming the pore module" is the claim the data supports.
    """
    import json
    from piezo1.config import RESOURCE_DIR
    st = human_structure
    blocks, _ = protomer_blocks(st)
    axis = detect_c3_axis(blocks)
    tms = json.loads((RESOURCE_DIR / "uniprot_human.json").read_text())["transmembrane"]
    tilts = {}
    for tm in tms:
        m = st.mask_ca() & (st.chain == "A") & (st.res_seq >= tm["start"]) \
            & (st.res_seq <= tm["end"])
        if m.sum() >= 6:
            xyz = st.xyz[m].astype(float)
            if axis.radial(xyz).mean() > 25.0:
                continue                      # out in the blade, not the pore
            a, _ = helix_axis(xyz)
            tilts[tm["name"]] = tilt_angle(a, axis.direction)
    assert "TM38" in tilts and len(tilts) >= 2
    assert tilts["TM38"] < 15.0
    assert tilts["TM38"] == min(tilts.values())


def test_sasa_of_an_isolated_atom(human_structure):
    """A lone atom must return the full area of its expanded sphere."""
    single = human_structure.subset(np.arange(1))
    r = sasa(single, probe=1.4, n_points=512)
    expected = 4 * np.pi * (single.vdw_radii()[0] + 1.4) ** 2
    assert r.total == pytest.approx(expected, rel=1e-6)


def test_sasa_is_deterministic_and_sane(human_structure):
    mask = human_structure.mask_protein() & ~human_structure.hetero
    a = sasa(human_structure, n_points=64, mask=mask)
    b = sasa(human_structure, n_points=64, mask=mask)
    assert np.array_equal(a.atom, b.atom)         # golden spiral, not random
    assert (a.atom >= 0).all()
    assert a.total > 10000
    assert len(a.residue) == a.residue_seq.size


def test_buried_area_is_positive_for_contacting_chains(human_structure):
    st = human_structure
    a = (st.chain == "A") & st.mask_protein()
    b = (st.chain == "B") & st.mask_protein()
    area = buried_area(st, a, b, n_points=64)
    assert area > 500        # protomers of a trimer bury a real interface


def test_hydrophobicity_scale_is_the_published_one():
    assert KYTE_DOOLITTLE["ILE"] == 4.5
    assert KYTE_DOOLITTLE["ARG"] == -4.5
    assert len(KYTE_DOOLITTLE) == 20


def test_hydrophobicity_profile_matches_the_pore(human_structure):
    blocks, _ = protomer_blocks(human_structure)
    axis = detect_c3_axis(blocks)
    prof = pore_profile(human_structure, axis, step=3.0)
    h = hydrophobicity_profile(human_structure, prof)
    assert len(h) == len(prof.z)
    finite = h[np.isfinite(h)]
    assert len(finite) > 5
    assert (finite >= -4.5).all() and (finite <= 4.5).all()
