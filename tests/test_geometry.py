"""Dome geometry, including the published-value regression test."""

import json

import numpy as np
import pytest

from piezo1.config import RESOURCE_DIR
from piezo1.structure.geometry import fit_sphere, measure_dome, radial_profile
from piezo1.structure.superpose import detect_c3_axis
from conftest import protomer_blocks


def test_fit_sphere_on_a_synthetic_sphere():
    rng = np.random.default_rng(3)
    centre = np.array([3.0, -7.0, 11.0])
    radius = 42.0
    v = rng.normal(size=(500, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    pts = centre + v * radius
    fit = fit_sphere(pts)
    assert fit.radius == pytest.approx(radius, rel=1e-6)
    assert np.allclose(fit.center, centre, atol=1e-6)
    assert fit.rmse < 1e-6


def test_fit_sphere_on_a_cap_only():
    """A spherical cap, not a whole sphere — the actual PIEZO1 situation."""
    rng = np.random.default_rng(4)
    radius = 100.0
    theta = rng.uniform(0, 0.6, 800)          # a cap, up to ~34 degrees
    phi = rng.uniform(0, 2 * np.pi, 800)
    pts = radius * np.stack([np.sin(theta) * np.cos(phi),
                             np.sin(theta) * np.sin(phi),
                             np.cos(theta)], axis=1)
    fit = fit_sphere(pts, iterations=4)
    assert fit.radius == pytest.approx(radius, rel=1e-3)


def test_radial_profile_of_a_known_cap():
    r = np.linspace(0, 50, 400)
    theta = np.linspace(0, 2 * np.pi, 40, endpoint=False)
    R, T = np.meshgrid(r, theta, indexing="ij")
    Z = -R ** 2 / 200.0
    pts = np.stack([(R * np.cos(T)).ravel(), (R * np.sin(T)).ravel(),
                    Z.ravel()], axis=1)
    axis = detect_c3_axis([pts, pts, pts])
    prof = radial_profile(pts, axis, n_bins=20)
    assert prof.count.sum() > 0
    assert np.isfinite(prof.z[prof.valid]).all()


def _tm_surface(st, species):
    tms = json.loads((RESOURCE_DIR / f"uniprot_{species}.json").read_text())["transmembrane"]
    pts = []
    for ch in st.chains:
        m = st.mask_ca() & (st.chain == ch)
        if m.sum() < 300:
            continue
        xyz, seq = st.xyz[m], st.res_seq[m]
        for tm in tms:
            mid = 0.5 * (tm["start"] + tm["end"])
            half = max(2.0, (tm["end"] - tm["start"]) / 6.0)
            sel = (seq >= mid - half) & (seq <= mid + half)
            if sel.sum() >= 3:
                pts.append(xyz[sel].mean(axis=0))
    return np.array(pts)


def test_curved_dome_matches_published_curvature(curved_structure):
    """Regression against Haselwandter & MacKinnon 2018: R_c = 10.2 nm.

    Our measurement of PDB 7WLT should land near that, and must be clearly
    distinguishable from the flattened state.
    """
    blocks, _ = protomer_blocks(curved_structure)
    dome = measure_dome(blocks, _tm_surface(curved_structure, "mouse"))
    r_nm = dome.radius_of_curvature / 10.0
    assert 8.0 < r_nm < 13.0, f"radius of curvature {r_nm:.1f} nm is off"
    assert dome.dome_depth / 10.0 > 3.0
    assert dome.notes["c3_angle_deg"] == pytest.approx(120.0, abs=0.05)


def test_flat_state_is_measurably_flatter(curved_structure, flat_structure):
    curved = measure_dome(protomer_blocks(curved_structure)[0],
                          _tm_surface(curved_structure, "mouse"))
    flat = measure_dome(protomer_blocks(flat_structure)[0],
                        _tm_surface(flat_structure, "mouse"))
    assert flat.radius_of_curvature > curved.radius_of_curvature * 1.4
    assert flat.dome_depth < curved.dome_depth
