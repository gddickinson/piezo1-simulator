"""Surface electrostatics, calibrated before any PIEZO1 number is believed.

The calibration here is not decoration. The first draft of
:mod:`piezo1.physics.electrostatics` had ``e^2 / 4 pi eps0`` written in
Joule-Angstrom and then converted from metres a second time, so the Bjerrum
length came out at 7.1e10 A, the Debye length underflowed to zero, every
potential was exactly zero, and the truncation-error check reported a flawless
0.000% because it was comparing zero with zero. Nothing raised. The
single-point-charge case below is what said no.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.physics.electrostatic_patches import (compare_conventions,
                                                  patch_interaction)
from piezo1.physics.electrostatics import (ChargeSet, bjerrum_length,
                                           debye_length, formal_charges,
                                           potential_at, surface_points,
                                           surface_potential)


# --------------------------------------------------------------------------
# Known answers
# --------------------------------------------------------------------------

def test_bjerrum_and_debye_lengths_match_the_published_values():
    """7.14 A and 7.86 A for water at 298 K and 150 mM 1:1 salt."""
    assert bjerrum_length() == pytest.approx(7.14, abs=0.01)
    assert debye_length() == pytest.approx(7.86, abs=0.02)
    # The textbook scaling: lambda_D = 3.04 / sqrt(I) Angstrom.
    for ionic in (0.01, 0.05, 0.15, 0.5):
        assert debye_length(ionic) == pytest.approx(3.04 / np.sqrt(ionic),
                                                    rel=0.01)
    assert debye_length(0.0) == float("inf")


def test_a_single_point_charge_reproduces_the_closed_form():
    """The case the units bug could not survive."""
    charge = ChargeSet(xyz=np.zeros((1, 3)), charge=np.array([1.0]))
    l_b, l_d = bjerrum_length(), debye_length()
    for r in (2.0, 5.0, 10.0, 25.0):
        got = potential_at(np.array([[r, 0.0, 0.0]]), charge)[0]
        assert got == pytest.approx(l_b * np.exp(-r / l_d) / r, rel=1e-12)


def test_a_shell_of_charges_reproduces_its_centre_potential():
    """N charges at radius a give N z l_B exp(-a/lambda)/a at the centre."""
    rng = np.random.default_rng(0)
    direction = rng.normal(size=(200, 3))
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    a = 12.0
    shell = ChargeSet(xyz=a * direction, charge=np.full(200, -1.0))
    got = potential_at(np.zeros((1, 3)), shell)[0]
    want = 200 * bjerrum_length() * (-1.0) * np.exp(-a / debye_length()) / a
    assert got == pytest.approx(want, rel=1e-9)


def test_removing_the_salt_gives_bare_coulomb():
    charge = ChargeSet(xyz=np.zeros((1, 3)), charge=np.array([1.0]))
    got = potential_at(np.array([[10.0, 0.0, 0.0]]), charge,
                       debye=float("inf"))[0]
    assert got == pytest.approx(bjerrum_length() / 10.0, rel=1e-12)


@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_the_sign_of_the_potential_follows_the_sign_of_the_charge(sign):
    """A check that would fail on a sign error anywhere in the sum."""
    rng = np.random.default_rng(1)
    charges = ChargeSet(xyz=rng.normal(scale=10.0, size=(50, 3)),
                        charge=np.full(50, sign))
    phi = potential_at(rng.normal(scale=25.0, size=(300, 3)), charges)
    assert np.all(np.sign(phi) == sign)


def test_the_truncation_cutoff_costs_a_measured_amount(structure_6b3r):
    """Against the exact all-pairs sum on the real surface, not on a bound."""
    charges = formal_charges(structure_6b3r)
    points, _ = surface_points(structure_6b3r)
    rng = np.random.default_rng(0)
    sample = points[rng.choice(len(points), 2000, replace=False)]
    exact = potential_at(sample, charges, cutoff=None)

    worse = np.abs(potential_at(sample, charges, cutoff=20.0) - exact).max()
    default = np.abs(potential_at(sample, charges, cutoff=30.0) - exact).max()
    assert default < worse, "a longer cutoff must be more accurate"
    assert default < 0.15, (
        f"30 A truncation costs {default:.3f} kT/e, more than documented")


def test_surface_points_reproduce_the_shrake_rupley_areas(structure_6b3r):
    """The generator is the same construction as ``sasa``, so it must agree.

    A surface generator that quietly kept buried points would still produce a
    plausible coloured picture; this is what stops that.
    """
    from piezo1.analysis.measure import sasa

    mask = structure_6b3r.mask_protein() & (~structure_6b3r.hetero)
    n_points = 64
    points, atom_index = surface_points(structure_6b3r, n_points=n_points,
                                        mask=mask)
    selected = np.flatnonzero(mask)
    radii = structure_6b3r.vdw_radii()[selected] + 1.4
    counts = np.bincount(atom_index,
                         minlength=structure_6b3r.n_atoms)[selected]
    implied = float((4 * np.pi * radii ** 2 * counts / n_points).sum())

    reference = sasa(structure_6b3r, mask=mask, n_points=n_points)
    assert implied == pytest.approx(float(reference.atom.sum()), rel=1e-6)


# --------------------------------------------------------------------------
# PIEZO1
# --------------------------------------------------------------------------

def test_piezo1_is_net_positive_as_the_inside_rule_expects(structure_6b3r):
    """Guo & MacKinnon invoke von Heijne's positive-inside rule for Figure 4c."""
    charges = formal_charges(structure_6b3r)
    assert len(charges) > 700
    assert charges.total > 0, "the trimer's formal charge should be positive"


def test_the_surface_never_reaches_the_panels_saturation(structure_6b3r):
    """The documented consequence of leaving out the dielectric boundary.

    Figure 4c saturates at +-5 k_BT/e and visibly does saturate. Ours does not
    get near it, and that is the *expected* direction of the approximation —
    a uniform solvent dielectric does not focus field lines out of a
    low-dielectric interior. If this ever starts saturating, the model has
    changed and the caveat has to change with it.
    """
    result = surface_potential(structure_6b3r)
    assert result.scale == pytest.approx(5.0)
    assert result.fraction_saturated() < 0.01
    assert np.abs(result.potential).max() < 5.0
    assert "APBS" in result.meta["not_apbs"]


def test_the_two_charged_patches_attract_and_do_so_across_protomers(
        structure_6b3r):
    """Figure 4-figure supplement 1's actual claim, which is not a colour."""
    from piezo1.analysis.guo2017_panels import (CED_ACIDIC_PATCH,
                                                LOOP_BASIC_PATCH)

    result = patch_interaction(structure_6b3r, CED_ACIDIC_PATCH,
                               LOOP_BASIC_PATCH)
    assert result["attractive"], "acidic against basic must be attractive"
    assert result["domain_swapped"], (
        "the paper's claim is specifically that the contacts are "
        "domain-swapped, so a same-chain term must not dominate")
    assert abs(result["same_chain_kT"]) < 0.1 * abs(result["cross_chain_kT"])


def test_the_interaction_instrument_can_say_no(structure_6b3r):
    """Two controls whose answers are known before the measurement.

    Without these, "the patches attract" is a number the function was always
    going to produce for any pair of residues near each other.
    """
    from piezo1.analysis.guo2017_panels import CED_ACIDIC_PATCH

    # A patch against itself is like-charged and must repel.
    same = patch_interaction(structure_6b3r, CED_ACIDIC_PATCH,
                             CED_ACIDIC_PATCH)
    assert not same["attractive"], "like charges must repel"
    assert same["energy_kT"] > 0

    # Something far away must contribute essentially nothing.
    charges = formal_charges(structure_6b3r)
    patch = np.array([charges.xyz[i] for i, label in enumerate(charges.label)
                      if label.split("/")[1][3:] in {"2257", "2258", "2264"}])
    far = sorted({int("".join(c for c in label.split("/")[1] if c.isdigit()))
                  for i, label in enumerate(charges.label)
                  if charges.charge[i] > 0
                  and np.linalg.norm(charges.xyz[i] - patch.mean(axis=0)) > 80})[:3]
    distant = patch_interaction(structure_6b3r, CED_ACIDIC_PATCH, far)
    assert abs(distant["energy_kT"]) < 0.01, (
        "screened Coulomb at 60+ A must be negligible; if it is not, the "
        "screening is not being applied")


def test_a_wrong_numbering_system_is_refused_not_answered(structure_6b3r):
    """Residues that are not ionisable in this entry must raise."""
    with pytest.raises(ValueError, match="numbering"):
        patch_interaction(structure_6b3r, (1, 2, 3), (4, 5, 6))


def test_the_open_conventions_are_reported_rather_than_chosen(structure_6b3r):
    """Salt and histidine move the surface; the table shows by how much."""
    mask = structure_6b3r.mask_residues([2257, 2258, 2264])
    rows = {r["convention"]: r for r in
            compare_conventions(structure_6b3r, mask=mask)}
    assert len(rows) == 6
    unscreened = rows["no salt (unscreened)"]
    default = rows["Figure 4c: 150 mM, neutral His"]
    assert unscreened["debye_length_A"] == float("inf")
    assert abs(unscreened["mean_kT_per_e"]) > abs(default["mean_kT_per_e"]), (
        "removing the screening must increase the magnitude")
    assert (rows["His fully protonated"]["n_charges"]
            > default["n_charges"]), "protonating histidine must add charges"
