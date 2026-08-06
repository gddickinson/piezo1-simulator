"""The nonlinear footprint propagated into the gating energetics.

Round 18 built the elastica solver and showed the linearised footprint is 3.5×
too large at PIEZO1's 63° contact slope. Only `DomeModel.footprint_nonlinear`
consumed it; the two-state model's ΔA still came from linear numbers.

The roadmap asked for one thing explicitly: report the change **including if
the linear version happened to agree better with the measurement**, because a
wrong model can fit a right number. It did not — but the test says so rather
than leaving it implied.
"""

import numpy as np
import pytest

from piezo1.parameters import PARAMETERS
from piezo1.physics.dome import DomeGeometrySummary, DomeModel
from piezo1.physics.membrane import kt_per_nm2_to_mnm

#: Measured from the deposited structures in this round. Curved 7WLT and
#: flattened 7WLU, the pair the whole gating analysis rests on.
CLOSED = DomeGeometrySummary(radius_of_curvature=9.72, dome_area=493.2,
                             projected_area=237.3)
OPEN = DomeGeometrySummary(radius_of_curvature=18.38, dome_area=816.8,
                           projected_area=438.1)


@pytest.fixture(scope="module")
def model():
    return DomeModel(geometry=CLOSED)


# --------------------------------------------------------------------------
# The area change
# --------------------------------------------------------------------------

def test_delta_area_is_a_change_not_an_absolute_area(model):
    """The two-state model is driven by the *change* in projected area.

    Quoting a footprint's absolute excess area as ΔA — which is what a naive
    reading of Round 3 invites — uses the wrong quantity entirely.
    """
    result = model.gating_area_change(OPEN)
    assert result["dome_term_nm2"] == pytest.approx(
        OPEN.projected_area - CLOSED.projected_area, rel=1e-9)
    assert result["delta_area_nm2"] == pytest.approx(
        result["dome_term_nm2"] + result["footprint_term_nm2"], rel=1e-9)


def test_the_footprint_releases_area_as_it_flattens(model):
    """Both endpoints store excess area, and the closed one stores more."""
    result = model.gating_area_change(OPEN)
    assert result["closed_footprint_nm2"] > result["open_footprint_nm2"] > 0
    assert result["footprint_term_nm2"] > 0


def test_the_nonlinear_correction_is_large(model):
    """The headline. Because the closed state sits at 63° where the linear
    theory is badly wrong and the open state at 40° where it is much less so,
    the *difference* is overstated far more than either endpoint."""
    result = model.gating_area_change(OPEN)
    linear = result["footprint_term_linear_nm2"]
    nonlinear = result["footprint_term_nonlinear_nm2"]
    assert linear == pytest.approx(463.0, abs=10.0)
    assert nonlinear == pytest.approx(71.0, abs=5.0)
    assert linear / nonlinear > 5.0


def test_total_delta_area_shrinks_by_more_than_half(model):
    linear = model.gating_area_change(OPEN, nonlinear=False)
    nonlinear = model.gating_area_change(OPEN, nonlinear=True)
    assert linear["delta_area_nm2"] == pytest.approx(664.0, abs=15.0)
    assert nonlinear["delta_area_nm2"] == pytest.approx(272.0, abs=10.0)


# --------------------------------------------------------------------------
# What it does to T50 — and what it does not
# --------------------------------------------------------------------------

def test_t50_moves_toward_the_measurement_not_away(model):
    """The question the roadmap asked. A wrong model can fit a right number,
    so this could have gone the other way."""
    measured = PARAMETERS.value("kinetics.t50_measured")     # 2.7 mN/m
    linear = model.gating_area_change(OPEN, nonlinear=False)["t50_mnm"]
    nonlinear = model.gating_area_change(OPEN, nonlinear=True)["t50_mnm"]
    assert abs(nonlinear - measured) < abs(linear - measured), (
        f"linear {linear:.3f} was closer to {measured} than nonlinear "
        f"{nonlinear:.3f} — the correction made the agreement worse")
    assert nonlinear > linear


def test_the_correction_does_not_close_the_gap(model):
    """The honest headline, and the more useful half of the result.

    Improving the membrane physics by a factor of six moves T₅₀ by a factor of
    2.4 and leaves it ~18× below measurement. So the structural-versus-
    functional discrepancy is **not** a membrane-modelling error — it is about
    which quantity each number measures.
    """
    measured = PARAMETERS.value("kinetics.t50_measured")
    nonlinear = model.gating_area_change(OPEN, nonlinear=True)["t50_mnm"]
    assert nonlinear < measured / 10.0, (
        "if this ever passes, the gap has closed and SCIENCE.md needs "
        "rewriting")


def test_structural_delta_area_still_dwarfs_the_functional_one(model):
    functional = PARAMETERS.value("dome.delta_area")          # 8 nm^2
    structural = model.gating_area_change(OPEN)["delta_area_nm2"]
    assert structural / functional > 20.0


# --------------------------------------------------------------------------
# The comparison table
# --------------------------------------------------------------------------

def test_every_route_is_reported_with_its_t50(model):
    rows = model.compare_gating_area_routes(OPEN)
    assert len(rows) == 4
    names = [r["route"] for r in rows]
    assert any("functional" in n for n in names)
    assert any("linear footprint" in n for n in names)
    assert any("nonlinear footprint" in n for n in names)
    for row in rows:
        assert row["delta_area_nm2"] > 0
        assert row["t50_mnm"] > 0
        assert row["note"]


def test_t50_follows_from_delta_area_consistently(model):
    """Every row must satisfy T50 = dG0/dA, so the table cannot drift from the
    model it is describing."""
    for row in model.compare_gating_area_routes(OPEN):
        expected = kt_per_nm2_to_mnm(model.delta_g0 / row["delta_area_nm2"])
        assert row["t50_mnm"] == pytest.approx(expected, rel=1e-9)


def test_the_functional_route_still_reproduces_cox(model):
    """The regression that must not move: 9.7 k_BT over 8 nm² gives 4.99 mN/m
    against Cox's measured 5.1 ± 0.2."""
    rows = model.compare_gating_area_routes(OPEN)
    functional = next(r for r in rows if "functional" in r["route"])
    assert functional["t50_mnm"] == pytest.approx(4.99, abs=0.02)


def test_geometry_is_required(model):
    bare = DomeModel()
    with pytest.raises(ValueError, match="no geometry"):
        bare.gating_area_change(OPEN)


def test_flatter_open_states_release_less(model):
    """A sanity check on the direction of the physics: the flatter the open
    state, the more area has been released getting there."""
    flatter = DomeGeometrySummary(radius_of_curvature=21.59, dome_area=825.4,
                                  projected_area=490.5)      # 11ZC
    to_flattened = model.gating_area_change(OPEN)["delta_area_nm2"]
    to_flat = model.gating_area_change(flatter)["delta_area_nm2"]
    assert to_flat > to_flattened
