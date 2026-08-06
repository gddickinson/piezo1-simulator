"""The parameter registry, its provenance gate, and the audit that enforces it.

The rule: **any number a calculation depends on is a registered parameter with a
unit, bounds and a citation.** A constant written into a function default is
invisible — it cannot be listed, shown to a user, or traced to a paper — and
this project has already had to correct several numbers that were invisible in
exactly that way.

These tests are the enforcement. The audit in particular is what stops the rule
decaying into an aspiration the moment somebody adds a new constant.
"""

import json

import pytest

from piezo1.config import RESOURCE_DIR
from piezo1.parameter_audit import EXEMPT, EXEMPT_NAMES, MAPPED, audit
from piezo1.parameters import PARAMETERS, Parameter, resolve


@pytest.fixture(autouse=True)
def pristine():
    """Every test starts and ends with the registry at its defaults."""
    PARAMETERS.reset()
    yield
    PARAMETERS.reset()


# --------------------------------------------------------------------------
# The registry itself
# --------------------------------------------------------------------------

def test_registry_loaded():
    assert len(PARAMETERS) > 50
    assert "membrane.kappa" in PARAMETERS
    assert PARAMETERS.value("membrane.kappa") == 20.0


def test_every_parameter_is_complete():
    for parameter in PARAMETERS.parameters.values():
        assert parameter.name and parameter.description
        assert parameter.kind in ("physical", "empirical", "method", "convention")
        assert parameter.category
        assert isinstance(parameter.default, (int, float))


def test_every_citation_resolves_or_declares_itself():
    """The provenance gate. A parameter either cites a real paper, or says
    explicitly that it is a method choice and why."""
    refs = {e["key"] for e in json.loads(
        (RESOURCE_DIR / "references.json").read_text())["references"]}
    sentinels = set(PARAMETERS.sentinels)
    for parameter in PARAMETERS.parameters.values():
        if parameter.citation in sentinels:
            assert parameter.source_note, (
                f"{parameter.key} claims '{parameter.citation}' but does not "
                f"say why")
        else:
            assert parameter.citation in refs, (
                f"{parameter.key} cites {parameter.citation!r}, which is not "
                f"in references.json")


def test_defaults_lie_within_their_own_bounds():
    for parameter in PARAMETERS.parameters.values():
        if parameter.minimum is not None:
            assert parameter.default >= parameter.minimum, parameter.key
        if parameter.maximum is not None:
            assert parameter.default <= parameter.maximum, parameter.key


def test_physical_and_empirical_parameters_cite_literature():
    """A measured quantity may not be a 'method choice'."""
    for parameter in PARAMETERS.parameters.values():
        if parameter.kind in ("physical", "empirical"):
            assert parameter.cited or parameter.citation == "measured_here", (
                f"{parameter.key} is {parameter.kind} but cites "
                f"{parameter.citation!r}")


# --------------------------------------------------------------------------
# Overriding
# --------------------------------------------------------------------------

def test_override_changes_the_value_and_is_tracked():
    assert PARAMETERS.is_default("membrane.kappa")
    PARAMETERS.set_value("membrane.kappa", 25.0)
    assert PARAMETERS.value("membrane.kappa") == 25.0
    assert not PARAMETERS.is_default("membrane.kappa")
    assert PARAMETERS.overrides() == {"membrane.kappa": 25.0}
    assert PARAMETERS.modified


def test_setting_back_to_the_default_clears_the_override():
    """Otherwise a user who fixes their own typo would leave the application
    permanently claiming its numbers are non-standard."""
    PARAMETERS.set_value("membrane.kappa", 25.0)
    PARAMETERS.set_value("membrane.kappa", 20.0)
    assert not PARAMETERS.modified


def test_out_of_range_values_are_clamped_not_raised():
    """These are edited in a UI with free-text entry; a typo must not take the
    application down, but the caller has to be told what was applied."""
    applied = PARAMETERS.set_value("membrane.kappa", 1e6)
    assert applied == PARAMETERS.get("membrane.kappa").maximum
    assert PARAMETERS.value("membrane.kappa") == applied


def test_unknown_parameter_raises_with_instructions():
    with pytest.raises(KeyError, match="build_parameters"):
        PARAMETERS.value("membrane.not_a_parameter")


def test_overrides_take_effect_on_the_next_call():
    """Modules resolve at call time, not import time, so a change does not
    require a reimport to bite."""
    from piezo1.physics.dome import DomeModel
    from piezo1.physics.membrane import MembraneParameters

    assert MembraneParameters().kappa == 20.0
    PARAMETERS.set_value("membrane.kappa", 30.0)
    assert MembraneParameters().kappa == 30.0

    before = DomeModel().half_activation_mnm
    PARAMETERS.set_value("dome.delta_area", 16.0)
    assert DomeModel().half_activation_mnm == pytest.approx(before / 2, rel=1e-6)


def test_resolve_prefers_an_explicit_argument():
    assert resolve(None, "anm.cutoff") == PARAMETERS.value("anm.cutoff")
    assert resolve(12.0, "anm.cutoff") == 12.0


def test_provenance_rows_carry_everything_a_reader_needs():
    rows = {r["key"]: r for r in PARAMETERS.provenance_rows()}
    row = rows["hydration.closed_cutoff"]
    assert row["citation"] == "rao2019heuristic"
    assert row["unit"] == ""
    assert row["overridden"] is False
    assert row["default"] == 0.55


# --------------------------------------------------------------------------
# Interaction with recorded results
# --------------------------------------------------------------------------

def test_claims_refuse_to_run_against_modified_parameters():
    """The important one.

    Every documented number was produced at the defaults. Recomputing with a
    changed value would report drift the user caused, and the obvious reading
    of that report is that the code is broken.
    """
    from piezo1.analysis.claims import claims_by_cost, verify_claims

    PARAMETERS.set_value("membrane.kappa", 25.0)
    with pytest.raises(RuntimeError, match="modified parameters"):
        verify_claims(claims_by_cost("fast"), verbose=False)

    results = verify_claims(claims_by_cost("fast"), verbose=False,
                            allow_overrides=True)
    assert all(not r.comparable for r in results)
    assert all(r.parameters == {"membrane.kappa": 25.0} for r in results)


def test_reports_record_and_flag_non_default_parameters():
    from piezo1.analysis.report import AnalysisReport, collect_provenance

    clean = AnalysisReport(provenance=collect_provenance(), title="t")
    assert clean.provenance.parameter_overrides == {}
    assert "all at documented defaults" in clean.to_markdown()

    PARAMETERS.set_value("dome.delta_area", 12.0)
    dirty = AnalysisReport(provenance=collect_provenance(), title="t")
    assert dirty.provenance.parameter_overrides == {"dome.delta_area": 12.0}
    assert dirty.provenance.warnings
    markdown = dirty.to_markdown()
    assert "Non-default parameters" in markdown
    # The banner must be near the top, not buried in the provenance footer.
    assert markdown.index("Non-default parameters") < markdown.index("Provenance")


# --------------------------------------------------------------------------
# The audit — what makes this a rule rather than an aspiration
# --------------------------------------------------------------------------

def test_no_unregistered_numbers_in_the_scientific_packages():
    """Every number in physics/, structure/ and analysis/ is either registered
    or exempt with a stated reason.

    If this fails, the fix is one of two things and never a third: register the
    number in `scripts/build_parameters.py` with a unit and a citation, or add
    it to EXEMPT in `piezo1/parameter_audit.py` saying why it is not a
    scientific parameter.
    """
    findings = audit()
    assert not findings, "\n".join(str(f) for f in findings)


def test_every_exemption_states_a_reason():
    for reason in list(EXEMPT.values()) + list(EXEMPT_NAMES.values()):
        assert isinstance(reason, str) and len(reason) > 10


def test_every_mapping_points_at_a_real_parameter():
    for (module, name), key in MAPPED.items():
        assert key in PARAMETERS, f"{module}:{name} maps to unknown {key}"


def test_the_audit_would_catch_a_new_unregistered_number(tmp_path):
    """A detector that has never detected is decoration."""
    package = tmp_path / "piezo1" / "physics"
    package.mkdir(parents=True)
    (package / "invented.py").write_text(
        "def compute(binding_energy: float = 7.25):\n    return binding_energy\n")
    findings = audit(tmp_path)
    assert any(f.name == "binding_energy" and f.value == 7.25 for f in findings)
