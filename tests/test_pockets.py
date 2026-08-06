"""Alpha-sphere pocket detection and ligand contact mapping."""

import numpy as np
import pytest

from piezo1.analysis.pockets import (alpha_spheres, find_pockets,
                                     ligand_contact_residues)
from piezo1.core.annotations import load_annotations


@pytest.fixture(scope="module")
def enclosed(human_structure):
    """Default parameters: enclosed cavities only."""
    return find_pockets(human_structure, min_neighbours=30, r_max=5.5,
                        cluster_distance=2.0)


@pytest.fixture(scope="module")
def grooves(human_structure):
    """Burial filter off: surface grooves included as well."""
    return find_pockets(human_structure, min_neighbours=0, r_max=6.5,
                        cluster_distance=3.0)


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

def test_circumsphere_of_a_known_tetrahedron():
    """A regular tetrahedron on the unit sphere must give back radius 1.

    The circumsphere solve is the whole method; if it drifts, every pocket is
    wrong in a way nothing downstream would catch.
    """
    verts = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
                     dtype=float) / np.sqrt(3)
    # Surround it so Delaunay produces that tetrahedron among others.
    rng = np.random.default_rng(0)
    shell = rng.normal(size=(60, 3))
    shell /= np.linalg.norm(shell, axis=1, keepdims=True)
    coords = np.vstack([verts, shell * 8.0])
    spheres = alpha_spheres(coords, r_min=0.5, r_max=2.0, min_neighbours=0)
    assert len(spheres) >= 1
    assert np.isclose(spheres.radii, 1.0, atol=1e-6).any()


def test_radius_filter_bounds_are_respected():
    rng = np.random.default_rng(1)
    coords = rng.normal(scale=12.0, size=(400, 3))
    s = alpha_spheres(coords, r_min=3.0, r_max=5.5, min_neighbours=0)
    assert (s.radii >= 3.0).all() and (s.radii <= 5.5).all()
    assert s.n_total > len(s)


def test_too_few_atoms_raises():
    with pytest.raises(ValueError, match="at least five"):
        alpha_spheres(np.zeros((4, 3)))


# --------------------------------------------------------------------------
# Percolation
# --------------------------------------------------------------------------

def test_burial_filter_prevents_percolation(enclosed, grooves):
    """Without it, single-linkage merges the whole exterior into one "pocket".

    On PIEZO1 that produced a top pocket of 408 000 A^3 with 601 lining
    residues, which is the protein's outside surface, not a cavity. The filter
    is what makes the output mean anything.
    """
    with_filter, without = enclosed, grooves
    assert max(len(p.residues) for p in with_filter) < 150
    assert max(len(p.residues) for p in without) > 300
    assert max(p.volume for p in with_filter) < max(p.volume for p in without) / 10


def test_pocket_volumes_are_physically_plausible(enclosed):
    """A druggable cavity is hundreds to a few thousand cubic Angstrom."""
    pockets = enclosed
    assert pockets
    volumes = np.array([p.volume for p in pockets])
    assert volumes.max() < 20000
    assert np.median(volumes) > 100


def test_volume_is_a_union_not_a_sum(enclosed):
    """Alpha spheres overlap heavily; summing 4/3 pi r^3 overcounts badly."""
    p = enclosed[0]
    naive = float((4.0 / 3.0 * np.pi * p.radii ** 3).sum())
    assert p.volume < naive / 2


def test_pockets_are_ranked_and_indexed(enclosed):
    pockets = enclosed
    assert [p.index for p in pockets] == list(range(1, len(pockets) + 1))
    scores = [p.volume * (0.5 + p.buriedness) for p in pockets]
    assert scores == sorted(scores, reverse=True)


def test_buriedness_is_a_fraction(enclosed):
    for p in enclosed:
        assert 0.0 <= p.buriedness <= 1.0


# --------------------------------------------------------------------------
# Recovering annotated sites
# --------------------------------------------------------------------------

def test_gate_and_anchor_sites_are_recovered(enclosed):
    """Two annotated sites should fall out of the geometry unprompted."""
    ann = load_annotations("human")
    pockets = enclosed
    for group_id, expected in (("hydrophobic_gate", 2), ("anchor_brake", 2)):
        target = set(ann.group(group_id).residues)
        best = max(len(p.contains_residues(target)) for p in pockets)
        assert best >= expected, f"{group_id}: only {best} of {len(target)}"


def test_yoda1_site_is_a_groove_not_an_enclosed_cavity(enclosed, grooves):
    """A negative worth pinning, because it is informative rather than a bug.

    The Yoda1 site is mapped by mutagenesis and docking, never by a
    co-structure, and this project's own annotation labels its evidence as
    "predicted". Searching for enclosed cavities recovers at most one of its
    three residues; allowing surface grooves recovers two. The site is
    interfacial, which is consistent with Yoda1 acting as a wedge from the
    lipid phase — and with a lipid occupying part of it in PDB 7WLT.
    """
    ann = load_annotations("human")
    target = set(ann.group("yoda1_pocket").residues)
    assert ann.group("yoda1_pocket").evidence == "predicted"

    best_enclosed = max(len(p.contains_residues(target)) for p in enclosed)
    best_groove = max(len(p.contains_residues(target)) for p in grooves)
    assert best_enclosed <= 1
    assert best_groove >= 2


# --------------------------------------------------------------------------
# Ligands
# --------------------------------------------------------------------------

def test_ligands_are_excluded_from_the_geometry_by_default(human_structure, enclosed):
    """A bound lipid fills the very cavity being searched for."""
    with_lig = find_pockets(human_structure, protein_only=False)
    without = enclosed
    assert with_lig and without
    # Excluding the lipid can only open space up, never close it.
    assert max(p.volume for p in without) >= max(p.volume for p in with_lig) * 0.9


def test_ligand_contacts_are_mapped(human_structure):
    contacts = ligand_contact_residues(human_structure)
    assert "L9Q" in contacts, "8YEZ contains the L9Q lipid"
    info = contacts["L9Q"]
    assert info["n_copies"] >= 1
    assert len(info["residues"]) > 5
    assert all(isinstance(r, int) for r in info["residues"])


def test_ligand_contacts_report_the_cutoff(human_structure):
    tight = ligand_contact_residues(human_structure, cutoff=3.5)
    loose = ligand_contact_residues(human_structure, cutoff=6.0)
    assert tight["L9Q"]["cutoff"] == 3.5
    assert len(loose["L9Q"]["residues"]) > len(tight["L9Q"]["residues"])
