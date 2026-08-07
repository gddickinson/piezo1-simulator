"""The score is named for what it measures, not for what it failed to predict.

Round 39 recorded that this number has a legitimate use — finding mechanically
coupled positions — and an illegitimate one. Round 58 makes that structural.
Until now the output was a class called ``VariantPrediction`` carrying a
``direction`` property whose docstring read "stiffening (LoF-like) or softening
(GoF-like)", which is the illegitimate reading written into the API.

Five pre-registered tests failed to support that reading, and Rounds 47 and 54
showed it cannot be supported by data this project could obtain. So the name
goes and the coupling analysis stays.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest

from piezo1.analysis import variant_impact

ROOT = Path(__file__).resolve().parents[1]


# ------------------------------------------- the misleading names are gone

def test_the_output_is_not_called_a_prediction():
    assert hasattr(variant_impact, "CouplingScore")
    assert not hasattr(variant_impact, "VariantPrediction"), (
        "the old name reads as a direction prediction and must not survive "
        "as an alias — an alias is how the reading comes back")


def test_no_direction_property_exists():
    """`direction` returning GoF/LoF language was the specific defect."""
    score = variant_impact.CouplingScore(residue=1, wt_aa="R", mut_aa="H",
                                         gating_cost_change=0.5)
    assert not hasattr(score, "direction")
    assert score.sign == "stiffening"


def test_the_sign_property_refuses_gain_loss_language():
    doc = inspect.getdoc(type(
        variant_impact.CouplingScore(residue=1, wt_aa="R", mut_aa="H")).sign)
    lowered = doc.lower()
    assert "not" in lowered
    assert "r2456" in lowered, (
        "the docstring should name the position that demonstrates the limit")


def test_the_field_is_named_for_the_quantity_not_for_a_free_energy():
    fields = set(variant_impact.CouplingScore.__dataclass_fields__)
    assert "gating_cost_change" in fields
    assert "ddg_gating" not in fields, (
        "'ddG' implies a thermodynamic free energy of folding, which this is not")


def test_nothing_in_the_package_still_uses_the_old_names():
    # Two files may say the old names, for opposite reasons:
    #  - hazards.py records what the name USED to be, which is the point;
    #  - run_validation.py keeps "ddg_normalised" as a STORED JSON key,
    #    because the Round 7 record is frozen.
    allowed = {"piezo1/ui/hazards.py", "scripts/run_validation.py"}
    offenders = []
    for path in list((ROOT / "piezo1").rglob("*.py")) + \
                list((ROOT / "scripts").rglob("*.py")):
        rel = str(path.relative_to(ROOT))
        if rel in allowed:
            continue
        text = path.read_text()
        for name in ("VariantPrediction", "ddg_gating", "ddg_normalised"):
            if name in text:
                offenders.append(f"{rel}: {name}")
    assert not offenders, offenders

    # And in run_validation.py the old name may appear ONLY as a stored key,
    # never as an attribute read.
    runner = (ROOT / "scripts" / "run_validation.py").read_text()
    assert ".ddg_normalised" not in runner
    assert ".ddg_gating" not in runner


def test_the_class_docstring_states_what_it_is_not():
    doc = inspect.getdoc(variant_impact.CouplingScore)
    assert "not a prediction" in doc.lower()
    assert "pre-registered" in doc.lower()


# ------------------------------------------ the frozen record still round-trips

def test_the_round7_runner_still_writes_its_recorded_keys():
    """Renaming a stored key would break the reproduction the script exists for.

    The Round 7 result is frozen. Its JSON uses "ddg" and "direction", and the
    runner must keep writing those names even though the attributes changed.
    """
    text = (ROOT / "scripts" / "run_validation.py").read_text()
    assert '"ddg": r["pred"].gating_cost_change' in text
    assert '"direction": r["pred"].sign' in text
    assert "frozen" in text.lower()


# ------------------------------------------------- the coupling map stays

def test_the_coupling_analysis_is_kept():
    """The legitimate use survives the retirement of the illegitimate one."""
    from piezo1.analysis.features import FEATURE_NOTES

    for column in ("prs_gate_response", "prs_coupling", "dcc_to_gate",
                   "betweenness"):
        assert column in FEATURE_NOTES, f"{column} is the coupling map"


def test_the_score_still_measures_coupling():
    """It must still do the thing it is now named for.

    A retirement that broke the calculation would be a deletion, not a rename.
    """
    rng = np.random.default_rng(0)
    coords = np.repeat(rng.normal(scale=9.0, size=(40, 3)), 3, axis=0)
    residues = np.tile(np.arange(1, 41), 3)
    model = variant_impact.VariantImpactModel(
        coords=coords, residues=residues,
        gating_vector=rng.normal(scale=0.1, size=coords.shape))
    score = model.predict(10, "A", "W")
    assert isinstance(score, variant_impact.CouplingScore)
    assert np.isfinite(score.gating_cost_change)
    assert score.sign in {"stiffening", "softening", "neutral"}


def test_no_user_facing_path_computes_the_score():
    """Measured in Round 58: the GUI and CLI never show this number.

    That narrows the hazard rather than removing it — a notebook user still
    reaches it — and the register says so rather than claiming an exposure the
    software does not have.
    """
    for tree in ("ui", ):
        for path in (ROOT / "piezo1" / tree).rglob("*.py"):
            text = path.read_text()
            assert "VariantImpactModel(" not in text, (
                f"{path.name} computes a coupling score; the hazard register "
                f"and this test both assume no GUI path does")

    cli = (ROOT / "piezo1" / "cli.py").read_text()
    assert "VariantImpactModel(" not in cli


def test_the_hazard_register_describes_the_real_exposure():
    from piezo1.ui.hazards import HAZARDS

    hazard = next(h for h in HAZARDS if h.key == "prediction_read_as_validated")
    assert "no GUI or" in hazard.scenario
    assert "CouplingScore" in hazard.guard
    assert "VariantPrediction" in hazard.guard, (
        "the register should record what the name used to be")
