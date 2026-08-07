"""Can this predictor ever be validated with the data that could exist?

Round 26 made the mechanical predictor sensitive to *which* substitution
occurred — within-position variance 4.9% → 52.5% — and Round 36 tested it under
pre-registration and got δ = −0.249, p = 0.405. The question this module answers
is the one that follows, and it is a **design** question rather than another
test: given the variants this project has, and the most it could ever have, is
an effect that size detectable at all?

**No new comparison is run here.** The effect sizes come from the recorded
validation results; running a fresh test would need its own pre-registration.
This asks only what those recorded numbers imply about what is achievable.

**The answer is no, and by a clear margin.**

===========================  ====  ===  ==============================
scenario                        n  MDE  power at the observed |δ|=0.249
===========================  ====  ===  ==============================
today                          34  .47  0.32
every directional variant      46  .41  0.45
optimistic ceiling             59  .36  **0.51**
required for 80% power        134  .25  0.80
===========================  ====  ===  ==============================

The ceiling assumes **every** one of Round 45's 35 harvest candidates is
hand-curated with a measured direction and survives the modelling gate at the
same 74% rate Round 36 saw. Even then the minimum detectable effect is 0.356,
larger than the 0.249 actually observed, and power is a coin flip.

So the honest position is not "we need more data" but "the data that could
exist is not enough for this effect size". Those are different statements, and
only the second is actionable — it says a fifth test on this variant set should
not be run, whatever predictor is put into it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["Scenario", "FeasibilityReport", "assess", "recorded_round",
           "observed_effect", "modelling_survival", "PairedRequirement",
           "paired_positions_required", "paired_feasibility",
           "shared_positions_available"]

#: Which recorded validation this round reasons about. The substitution-aware
#: predictor is Round 36; Round 26 is what made it substitution-aware.
ROUND = 36


def recorded_round(number: int = ROUND):
    """The recorded result, read rather than restated.

    Every number below is derived from :data:`prediction_record.VALIDATION_RECORD`
    instead of being written out again here. That is deliberate: a copied effect
    size can drift away from the document it came from, and this module's whole
    argument is about that effect size.
    """
    from .prediction_record import VALIDATION_RECORD

    for record in VALIDATION_RECORD:
        if record.round == number:
            return record
    raise KeyError(f"no recorded validation for round {number}")


def observed_effect(number: int = ROUND) -> float:
    """The effect the predictor actually produced. Not re-measured here."""
    return float(recorded_round(number).cliffs_delta)


def modelling_survival(number: int = ROUND) -> float:
    """Fraction of directional variants that survived the inclusion criteria.

    Mostly the requirement that the residue be modelled in the reference pair.
    Used to project a ceiling rather than assuming every curated variant would
    be testable.
    """
    from .variant_sets import build_analysis_set

    record = recorded_round(number)
    tested = record.n_gof + record.n_lof
    directional = len(build_analysis_set(
        levels=("measured", "disease_mechanism")).missense())
    return tested / directional


def observed_split(number: int = ROUND) -> float:
    """The gain-of-function share of the tested set.

    Every projection keeps this imbalance rather than assuming a balanced set,
    because gain-of-function variants are the easier ones to find and a skewed
    set has *less* power, not more.
    """
    record = recorded_round(number)
    return record.n_gof / (record.n_gof + record.n_lof)


@dataclass
class Scenario:
    """One hypothetical dataset size, and what it could detect."""

    label: str
    n: int
    n_a: int
    n_b: int
    minimum_detectable: float = float("nan")
    power_at_observed: float = float("nan")
    reachable: bool = True
    note: str = ""

    observed_effect: float = float("nan")

    @property
    def detects_observed(self) -> bool:
        return abs(self.observed_effect) >= self.minimum_detectable

    def summary(self) -> str:
        return (f"{self.label}: n={self.n} ({self.n_a}/{self.n_b}), "
                f"MDE {self.minimum_detectable:.3f}, power at the observed "
                f"effect {self.power_at_observed:.2f}"
                + ("" if self.reachable else " — NOT reachable"))


@dataclass
class FeasibilityReport:
    """Whether the central claim is testable with data that could exist."""

    scenarios: list = field(default_factory=list)
    observed_effect: float = float("nan")
    ceiling_n: int | None = None
    required_n: int | None = None
    meta: dict = field(default_factory=dict)

    def get(self, label: str) -> Scenario | None:
        return next((s for s in self.scenarios if s.label == label), None)

    @property
    def harvest_available(self) -> bool:
        """Whether the literature corpus was present when this was computed.

        False means the ceiling is an underestimate and must not be quoted:
        the harvest contributes zero candidates rather than its real 35.
        """
        return bool(self.meta.get("harvest_available", True))

    @property
    def achievable(self) -> bool:
        """Could any reachable dataset detect the observed effect?"""
        return any(s.detects_observed for s in self.scenarios if s.reachable)

    def summary(self) -> str:
        ceiling = self.get("optimistic ceiling")
        if not self.harvest_available:
            return ("the open-access corpus is not downloaded, so the ceiling "
                    "cannot be computed — run `python -m piezo1.io.fetch`")
        return (f"observed |δ| = {abs(self.observed_effect):.3f}; "
                f"{self.required_n} variants would be needed for 80% power; "
                f"the optimistic ceiling is {self.ceiling_n}, where power is "
                f"{ceiling.power_at_observed:.2f} and the minimum detectable "
                f"effect is {ceiling.minimum_detectable:.3f}. "
                f"Achievable: {self.achievable}.")


def assess(n_simulations: int = 1500, split: float | None = None,
           number: int = ROUND) -> FeasibilityReport:
    """Project what each achievable dataset size could detect.

    ``split`` defaults to the gain-of-function share the recorded round
    actually had; every enlargement is assumed to keep that imbalance.
    """
    from .design import minimum_detectable_effect, power_curve, sample_size_for
    from .harvest import harvest
    from .variant_sets import build_analysis_set

    effect = observed_effect(number)
    survival = modelling_survival(number)
    if split is None:
        split = observed_split(number)
    record = recorded_round(number)
    tested = record.n_gof + record.n_lof

    directional = len(build_analysis_set(
        levels=("measured", "disease_mechanism")).missense())
    report = harvest()
    fresh = len([c for c in report.passing()
                 if c.human_label and not c.already_curated])

    # The harvest needs the downloaded open-access corpus. Without it the
    # scan finds nothing and the ceiling would silently fall from 59 to 34 —
    # a documented number quietly changing with the cache state, which Round 60
    # found on an empty clone. Report the absence instead of the smaller number.
    harvest_available = report.n_papers > 0
    if not harvest_available:
        fresh = 0
    ceiling = int((directional + fresh) * survival)
    required = 2 * sample_size_for(abs(effect), n_simulations=500,
                                   n_permutations=299)

    plan = [
        ("today", tested, True,
         f"the set Round {number} actually tested"),
        ("every directional variant", directional, True,
         "if every curated directional variant were modelled"),
        ("optimistic ceiling", ceiling, True,
         f"if all {fresh} harvest candidates were hand-curated with a measured "
         f"direction and survived the modelling gate"),
        ("required for 80% power", required, False,
         "no route to this number exists"),
    ]

    scenarios = []
    for label, n, reachable, note in plan:
        n_a = int(round(n * split))
        n_b = n - n_a
        scenarios.append(Scenario(
            label=label, n=n, n_a=n_a, n_b=n_b,
            minimum_detectable=float(minimum_detectable_effect(n_a, n_b)),
            power_at_observed=float(power_curve(
                n_a, n_b, deltas=[effect],
                n_simulations=n_simulations, seed=1).power[0]),
            observed_effect=effect, reachable=reachable, note=note))

    return FeasibilityReport(
        scenarios=scenarios, ceiling_n=ceiling, required_n=required,
        observed_effect=effect,
        meta={"round": number, "tested": tested, "split": split,
              "harvest_available": harvest_available,
              "n_papers_scanned": report.n_papers,
              "directional_variants": directional,
              "fresh_harvest_candidates": fresh,
              "modelling_survival": survival,
              "note": "no new comparison is run here; the effect size is the "
                      "one recorded in " + record.document})


def _main() -> None:
    report = assess()
    print(f"Round {report.meta['round']} recorded delta = "
          f"{report.observed_effect:+.3f}\n")
    for scenario in report.scenarios:
        print("  " + scenario.summary())
        print(f"      {scenario.note}")
    print("\n" + report.summary())


if __name__ == "__main__":
    _main()


# --------------------------------------------------------------------------
# The other design: comparing two variants at the SAME position
# --------------------------------------------------------------------------
#
# Round 47 costed the across-position route and closed it. The within-position
# route was the one it left open, because pairing removes the between-position
# variance that consumed 99.8% of Round 7's predictor. Round 61 costs it.
#
# The natural statistic is a **sign test**: at each position carrying variants
# of both directions, does the predictor rank the gain-of-function one above
# the loss-of-function one? Under the null that is a coin flip, so no
# distributional assumption is needed — which matters, because there is no
# sample to estimate a distribution from.
#
# Cliff's delta for a paired ordering is 2p - 1, so the requirement can be
# stated against the same effect scale the other rounds use.

def shared_positions_available() -> int:
    """Positions carrying variants of both directions — measured, not stated.

    Derived from :func:`data_routes.discriminating_positions` rather than
    written down, so the requirement and the supply cannot disagree. Today it
    is one: R2456.
    """
    from .data_routes import discriminating_positions

    return len(discriminating_positions())


@dataclass
class PairedRequirement:
    """How many shared positions a within-position test would need."""

    delta: float
    positions: int | None       # None when the effect is undetectable at any n
    available: int = field(default_factory=shared_positions_available)

    @property
    def reachable(self) -> bool:
        return self.positions is not None and self.positions <= self.available

    def summary(self) -> str:
        if self.positions is None:
            return f"paired δ {self.delta:+.2f}: not detectable at any sample size"
        return (f"paired δ {self.delta:+.2f}: needs {self.positions} shared "
                f"positions, {self.available} available")


def paired_positions_required(delta: float, alpha: float | None = None,
                              power: float | None = None, n_max: int = 400,
                              n_simulations: int = 4000,
                              seed: int = 0) -> PairedRequirement:
    """Shared positions needed for a sign test to detect a paired effect.

    ``delta`` is Cliff's delta on the paired ordering, so the probability of
    ranking a position's pair correctly is ``0.5 + delta / 2``. At ``delta = 0``
    the predictor is a coin flip and no sample size suffices; the function
    returns ``None`` rather than ``n_max``, because returning the search bound
    would look like an answer.
    """
    import numpy as np
    from scipy import stats

    from ..parameters import PARAMETERS

    if alpha is None:
        alpha = PARAMETERS.value("stats.alpha")
    if power is None:
        power = PARAMETERS.value("stats.target_power")

    p = 0.5 + float(delta) / 2.0
    if not 0.0 < p <= 1.0:
        raise ValueError(f"delta {delta} puts p outside (0, 1]")

    rng = np.random.default_rng(seed)
    for n in range(4, n_max + 1):
        hits = rng.binomial(n, p, size=n_simulations)
        critical = stats.binom.ppf(1.0 - alpha, n, 0.5)
        if (hits > critical).mean() >= power:
            return PairedRequirement(delta=float(delta), positions=n)
    return PairedRequirement(delta=float(delta), positions=None)


def paired_feasibility(deltas=(0.249, 0.35, 0.5, 0.7, 0.8), **kw) -> list:
    """The requirement across a range of paired effect sizes.

    A range rather than a point, because the paired effect **cannot be
    measured without running the comparison** — which the pre-registration
    protocol forbids until a design is registered. So this reports what each
    hypothetical effect would need, and the reader compares that with what
    exists.
    """
    return [paired_positions_required(d, **kw) for d in deltas]
