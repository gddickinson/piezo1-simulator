"""Which interval to publish, when several kinds disagree about the width.

:mod:`uncertainty` computes three kinds of spread and refuses to conflate them:
a **bootstrap** confidence interval, a **sensitivity** range over a method
choice, and a **parameter range** propagated from a registered input.
:mod:`model_error` adds a fourth, the spread over the *form* of the model.

Round 38 then measured something awkward: for the dome radius the model spread
is **six times** the bootstrap interval. So quoting ±0.9 nm — which is what a
bootstrap on the surface points gives — states the precision of a sphere fit
while the open question is whether a sphere is the right shape at all. The
number was not wrong; the interval on it was answering a different question.

**The rule this module applies.**

1. Compute every term that applies to the quantity.
2. Publish the **widest** one, and name its kind in the same breath. A number
   quoted with a confidence interval when a larger model spread exists is
   overconfident even though every individual figure is correct.
3. Never call a sensitivity range or a model spread a confidence interval. A
   network cutoff has no sampling distribution.
4. Where the widest term is a model spread, say that it is a **lower bound**:
   two models agreeing does not bound the error from above, because both may
   be wrong in the same direction.

The point estimate never moves. This changes only what is claimed *about* it,
which is why ``verify_claims`` is unaffected — the claim tolerances detect code
drift and are a separate question from what the science section publishes.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["Term", "PublishedInterval", "HEADLINE", "publish", "KINDS"]

#: The four kinds of spread, and whether each may be called a confidence
#: interval. Only a bootstrap resamples data, so only a bootstrap may.
KINDS = {
    "bootstrap": True,
    "sensitivity": False,
    "parameter": False,
    "model": False,
}


@dataclass(frozen=True)
class Term:
    """One kind of spread on one quantity."""

    kind: str
    low: float
    high: float
    source: str              # what was varied
    note: str = ""

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"unknown kind {self.kind!r}")
        if self.high < self.low:
            raise ValueError(f"{self.kind}: high {self.high} below low {self.low}")

    @property
    def width(self) -> float:
        return self.high - self.low

    @property
    def is_confidence_interval(self) -> bool:
        return KINDS[self.kind]

    def describe(self) -> str:
        label = {"bootstrap": "95% CI",
                 "sensitivity": "range over a method choice — NOT a CI",
                 "parameter": "propagated from a registered input — NOT a CI",
                 "model": "spread over the model form — NOT a CI, and a LOWER "
                          "BOUND"}[self.kind]
        return f"[{self.low:.4g}, {self.high:.4g}] ({label}; {self.source})"


@dataclass
class PublishedInterval:
    """A headline number and the interval that should be quoted with it."""

    quantity: str
    estimate: float
    unit: str
    terms: list = field(default_factory=list)
    published_value: str = ""

    @property
    def dominant(self) -> Term:
        """The widest term. This is what gets published."""
        return max(self.terms, key=lambda t: t.width)

    @property
    def low(self) -> float:
        return self.dominant.low

    @property
    def high(self) -> float:
        return self.dominant.high

    @property
    def overconfident_by(self) -> float:
        """How many times too narrow the narrowest term would have been.

        The number Round 38 reported for the dome. A value near 1 means the
        choice of interval does not matter; a large value means quoting the
        wrong kind would have materially overstated the precision.
        """
        widths = [t.width for t in self.terms if t.width > 0]
        return max(widths) / min(widths) if len(widths) > 1 else 1.0

    def statement(self) -> str:
        term = self.dominant
        text = (f"{self.quantity}: {self.estimate:.4g} {self.unit}, "
                f"{term.describe()}")
        if self.published_value:
            text += f"; published {self.published_value}"
        if term.note:
            text += f". {term.note}"
        return text


#: The four headline numbers Round 52 was asked about, with every term measured.
#: Regenerate with ``scripts/report_uncertainty.py`` and ``scripts/model_error.py``.
HEADLINE: tuple = (
    PublishedInterval(
        quantity="Dome radius of curvature (7WLT)", estimate=9.72, unit="nm",
        published_value="10.2 nm (Haselwandter & MacKinnon 2018)",
        terms=[
            Term("bootstrap", 8.80, 10.30,
                 "resampling the 66 transmembrane surface points"),
            Term("parameter", 9.43, 9.73,
                 "outlier trim fraction from 0 to 0.25 (geometry.sphere_trim)",
                 note="Measured while writing this: the registered trim moves "
                      "the radius by 0.30 nm. It matters because the model "
                      "comparison below is anchored on the UNTRIMMED fit "
                      "(9.45 nm) while the published number is trimmed "
                      "(9.72 nm). The gap is 0.27 nm against a 5.54 nm model "
                      "spread, so the conclusion is unaffected — but the two "
                      "were not like-for-like and that had not been noticed."),
            Term("model", 9.45, 14.99,
                 "sphere versus oblate spheroid (apex curvature), both untrimmed",
                 note="The spheroid fits BETTER (rmse 5.24 vs 6.18 A) and has "
                      "flattening +0.431, so the surface is not spherical. The "
                      "sphere stays the comparator because the published value "
                      "is itself a sphere fit — but the shape assumption, not "
                      "the point scatter, is what limits this number."),
        ]),

    PublishedInterval(
        quantity="Lowest A-mode gating overlap", estimate=0.705, unit="",
        terms=[
            Term("sensitivity", 0.554, 0.723,
                 "elastic-network cutoff from 10 to 20 A",
                 note="The qualitative result — a single symmetric mode "
                      "captures most of the transition — survives every "
                      "cutoff. The third digit does not."),
            Term("model", 0.890, 0.937,
                 "uniform / inverse-square / inverse-sixth spring models"),
        ]),

    PublishedInterval(
        quantity="Half-activation tension T50", estimate=2.711, unit="mN/m",
        published_value="2.7 +- 0.1 mN/m (Lewis & Grandl 2015)",
        terms=[
            Term("sensitivity", 2.711, 2.727,
                 "matrix exponential versus adaptive ODE integration",
                 note="The two solvers agree to 0.6%, so the numerical route "
                      "is not what limits this."),
            Term("parameter", 2.584, 2.838,
                 "Young 2023 rate constants at +-20%",
                 note="The published measurement lies inside this range, so "
                      "the agreement survives the uncertainty on the inputs "
                      "rather than depending on their exact values."),
        ]),

    PublishedInterval(
        quantity="Nonlinear footprint energy", estimate=25.27, unit="k_BT",
        terms=[
            Term("parameter", 25.27, 26.94,
                 "membrane bending modulus kappa from 20 to 25 k_BT",
                 note="kappa is a measured quantity with a literature range, "
                      "so this is propagated input uncertainty rather than a "
                      "modelling choice."),
        ]),
)


def publish(quantity: str) -> PublishedInterval | None:
    """The entry whose quantity name contains ``quantity``, case-insensitively."""
    key = quantity.lower()
    return next((h for h in HEADLINE if key in h.quantity.lower()), None)
