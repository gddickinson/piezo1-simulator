"""Tension-dependent Markov gating kinetics.

This module closes the link between structure and electrophysiology: given a
membrane tension it predicts open probability, macroscopic current and
single-channel behaviour, so a structural change can be turned into something a
patch-clamp rig would actually measure.

The implementation is **Model C of Young, Sindoni, Lewis, Zauscher & Grandl
2023 (PNAS)** — a four-state scheme C ⇌ O ⇌ I₁ ⇌ I₂ parameterised directly in
membrane tension rather than pipette pressure. It was chosen over the four
other published PIEZO1 models because it is the only one whose rates *and*
their tension dependence are all given numerically; the others have rate
constants that exist only inside figure panels.

Topology and rates (σ in mN/m, all rates s⁻¹)::

        k₋₃ ← C → k₁(σ)              k₁(σ) = 5.1  · exp( σ/b)     C  → O
              ↑ ↓                    k₋₁   = 5    · exp(σ₅₀/b)    O  → C
        k₃(σ) │ │ k₋₁                k₂    = 8.0                  O  → I₁
              │ ↓                    k₋₂   = 0.4                  I₁ → O
        I₁ ⇌ O ⇌ I₂                  k₃(σ) = 34.6 · exp(−σ/b)     I₁ → C
                                     k₄    = 4.0                  O  → I₂
                                     k₋₄   = 0.6                  I₂ → O

with σ₅₀ = 1.4 mN/m and slope factor b = 0.8 mN/m. The remaining rate C → I₁ is
**not free**: it is fixed by microscopic reversibility around the C-O-I₁ cycle,
k₋₃ = k₁·k₂·k₃ / (k₋₁·k₋₂). Setting it independently would let the model
generate energy from nothing.

Two caveats carried from the source literature and worth respecting:

* The topology of the **second inactivated state is genuinely ambiguous**. An
  equally good fit (2.0% vs 1.9% residual) places I₂ downstream of I₁ rather
  than of O. ``GatingModel.i2_downstream_of_i1`` switches between them.
* The gain-of-function **activation latency (344 ± 133 ms) cannot be
  represented by any Markov model**, per the authors. It is deliberately absent
  here rather than faked.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.linalg import expm, null_space
from ..parameters import PARAMETERS as _P

__all__ = ["GatingModel", "GatingResult", "STATE_NAMES", "MUTANT_PRESETS"]

STATE_NAMES = ("C", "O", "I1", "I2")
C, O, I1, I2 = 0, 1, 2, 3


@dataclass
class GatingResult:
    """Time course of a simulated stimulus."""

    time: np.ndarray                 # seconds
    occupancy: np.ndarray            # (n_time, 4)
    tension: np.ndarray              # mN/m at each time point
    current: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    @property
    def open_probability(self) -> np.ndarray:
        return self.occupancy[:, O]

    def peak_open(self) -> float:
        return float(self.occupancy[:, O].max())

    def peak_time(self) -> float:
        return float(self.time[int(np.argmax(self.occupancy[:, O]))])

    def inactivation_tau(self) -> float | None:
        """Single-exponential time constant of decay from the peak.

        Fitted by linear regression on log(P_open) between the peak and the
        point where it has fallen to 1/e of the way to its final value, which
        is the same operational definition used in the electrophysiology
        literature.
        """
        p = self.occupancy[:, O]
        i0 = int(np.argmax(p))
        if i0 >= len(p) - 3:
            return None
        tail = p[i0:]
        floor = tail[-1]
        amp = tail[0] - floor
        if amp <= 1e-9:
            return None
        target = floor + amp / np.e
        below = np.flatnonzero(tail <= target)
        if len(below) == 0:
            return None
        i1 = i0 + int(below[0])
        seg = p[i0:i1 + 1] - floor
        t = self.time[i0:i1 + 1]
        good = seg > 1e-12
        if good.sum() < 3:
            return float(t[-1] - t[0])
        slope = np.polyfit(t[good], np.log(seg[good]), 1)[0]
        return float(-1.0 / slope) if slope < 0 else None


@dataclass
class GatingModel:
    """Four-state tension-dependent PIEZO1 gating model (Young et al. 2023)."""

    sigma_50: float = field(
        default_factory=lambda: _P.value("kinetics.sigma_50"))        # mN/m
    b: float = field(
        default_factory=lambda: _P.value("kinetics.b"))               # mN/m, slope factor
    k1_0: float = field(
        default_factory=lambda: _P.value("kinetics.k1_0"))            # s^-1, C -> O prefactor
    km1_0: float = field(
        default_factory=lambda: _P.value("kinetics.km1_0"))           # s^-1, O -> C prefactor
    k2: float = field(
        default_factory=lambda: _P.value("kinetics.k2"))              # s^-1, O -> I1
    km2: float = field(
        default_factory=lambda: _P.value("kinetics.km2"))             # s^-1, I1 -> O
    k3_0: float = field(
        default_factory=lambda: _P.value("kinetics.k3_0"))           # s^-1, I1 -> C prefactor
    k4: float = field(
        default_factory=lambda: _P.value("kinetics.k4"))              # s^-1, O -> I2
    km4: float = field(
        default_factory=lambda: _P.value("kinetics.km4"))             # s^-1, I2 -> O
    #: Alternative topology in which I2 sits downstream of I1 rather than O.
    i2_downstream_of_i1: bool = False
    #: Single-channel conductance (pS) and holding potential (mV).
    conductance_pS: float = field(
        default_factory=lambda: _P.value("kinetics.conductance_pS"))
    holding_mV: float = -80.0
    label: str = "wild type"
    provenance: str = ("Young, Sindoni, Lewis, Zauscher & Grandl 2023, PNAS — "
                       "Model C, four-state, parameterised in tension")

    # ------------------------------------------------------------- rates

    def rates(self, tension: float) -> dict[str, float]:
        """All eight rate constants at a given tension, in s⁻¹."""
        k1 = self.k1_0 * np.exp(tension / self.b)
        km1 = self.km1_0 * np.exp(self.sigma_50 / self.b)
        k3 = self.k3_0 * np.exp(-tension / self.b)
        # Microscopic reversibility around C -> O -> I1 -> C fixes C -> I1.
        km3 = (k1 * self.k2 * k3) / (km1 * self.km2)
        return {"k1": float(k1), "km1": float(km1), "k2": self.k2,
                "km2": self.km2, "k3": float(k3), "km3": float(km3),
                "k4": self.k4, "km4": self.km4}

    def rate_matrix(self, tension: float) -> np.ndarray:
        """Generator matrix Q with ``Q[i, j]`` the i→j rate and zero row sums."""
        r = self.rates(tension)
        q = np.zeros((4, 4))
        q[C, O] = r["k1"]
        q[O, C] = r["km1"]
        q[O, I1] = r["k2"]
        q[I1, O] = r["km2"]
        q[I1, C] = r["k3"]
        q[C, I1] = r["km3"]
        if self.i2_downstream_of_i1:
            q[I1, I2] = r["k4"]
            q[I2, I1] = r["km4"]
        else:
            q[O, I2] = r["k4"]
            q[I2, O] = r["km4"]
        np.fill_diagonal(q, 0.0)
        np.fill_diagonal(q, -q.sum(axis=1))
        return q

    def detailed_balance_residual(self, tension: float) -> float:
        """Cycle-product mismatch around C-O-I₁. Should be ~0 by construction."""
        r = self.rates(tension)
        fwd = r["k1"] * r["k2"] * r["k3"]
        rev = r["km3"] * r["km2"] * r["km1"]
        return float(abs(fwd - rev) / max(fwd, 1e-30))

    # -------------------------------------------------------- equilibrium

    def steady_state(self, tension: float) -> np.ndarray:
        """Equilibrium occupancy of the four states."""
        q = self.rate_matrix(tension)
        ns = null_space(q.T)
        if ns.shape[1] == 0:                       # numerically degenerate
            vals, vecs = np.linalg.eig(q.T)
            ns = np.real(vecs[:, [int(np.argmin(np.abs(vals)))]])
        p = np.abs(ns[:, 0])
        return p / p.sum()

    def open_probability(self, tension: float | np.ndarray) -> np.ndarray:
        """Steady-state open probability at one or many tensions."""
        t = np.atleast_1d(np.asarray(tension, dtype=float))
        return np.array([self.steady_state(x)[O] for x in t])

    def peak_open_probability(self, tension: float | np.ndarray,
                              resting_tension: float = 0.0,
                              duration: float = 0.5) -> np.ndarray:
        """Peak open probability after a step from rest, as measured in a patch.

        Experiments report the *peak* of the transient, not the steady state,
        because PIEZO1 inactivates. Comparing a model's steady state with a
        measured peak is a common and serious error.
        """
        t = np.atleast_1d(np.asarray(tension, dtype=float))
        out = np.empty(len(t))
        p0 = self.steady_state(resting_tension)
        for i, x in enumerate(t):
            res = self.step(x, duration=duration, resting=resting_tension,
                            p0=p0, n_points=400)
            out[i] = res.peak_open()
        return out

    def half_activation(self, peak: bool = True, lo: float = 0.0,
                        hi: float = 20.0, n: int = 200) -> float:
        """Tension at which open probability reaches half its maximum."""
        grid = np.linspace(lo, hi, n)
        po = self.peak_open_probability(grid) if peak else self.open_probability(grid)
        target = 0.5 * po.max()
        idx = np.flatnonzero(po >= target)
        if len(idx) == 0:
            return float("nan")
        i = int(idx[0])
        if i == 0:
            return float(grid[0])
        # Linear interpolation between the bracketing grid points.
        x0, x1 = grid[i - 1], grid[i]
        y0, y1 = po[i - 1], po[i]
        return float(x0 + (target - y0) * (x1 - x0) / max(y1 - y0, 1e-12))

    # ------------------------------------------------------------ dynamics

    def step(self, tension: float, duration: float = 0.5,
             resting: float = 0.0, p0: np.ndarray | None = None,
             n_points: int = 1000, n_channels: int = 1) -> GatingResult:
        """Response to a step of tension held for ``duration`` seconds."""
        q = self.rate_matrix(tension)
        t = np.linspace(0.0, duration, n_points)
        if p0 is None:
            p0 = self.steady_state(resting)
        occ = np.empty((n_points, 4))
        # Propagate with the matrix exponential of a constant generator; exact
        # for a square stimulus, and far more stable than integrating the ODE.
        dt = t[1] - t[0]
        step_op = expm(q.T * dt)
        p = p0.copy()
        for i in range(n_points):
            occ[i] = p
            p = step_op @ p
        res = GatingResult(time=t, occupancy=occ,
                           tension=np.full(n_points, tension),
                           meta={"protocol": "step", "tension": tension,
                                 "resting": resting, "model": self.label})
        res.current = self.current(occ[:, O], n_channels)
        return res

    def ramp(self, tension_max: float, duration: float = 1.0,
             n_points: int = 1000, n_channels: int = 1) -> GatingResult:
        """Response to a linear tension ramp from zero."""
        t = np.linspace(0.0, duration, n_points)
        sigma = np.linspace(0.0, tension_max, n_points)
        occ = np.empty((n_points, 4))
        p = self.steady_state(0.0)
        dt = t[1] - t[0]
        for i in range(n_points):
            occ[i] = p
            p = expm(self.rate_matrix(sigma[i]).T * dt) @ p
        res = GatingResult(time=t, occupancy=occ, tension=sigma,
                           meta={"protocol": "ramp", "model": self.label})
        res.current = self.current(occ[:, O], n_channels)
        return res

    def current(self, open_probability: np.ndarray, n_channels: int = 1) -> np.ndarray:
        """Macroscopic current in pA from open probability.

        i = γ·V, with γ in pS and V in mV, giving fA; scaled to pA.
        """
        unitary_pA = self.conductance_pS * self.holding_mV * 1e-3
        return np.asarray(open_probability) * n_channels * unitary_pA

    # ----------------------------------------------------- single channel

    def simulate_single_channel(self, tension: float, duration: float = 1.0,
                                seed: int | None = None,
                                start_state: int | None = None
                                ) -> tuple[np.ndarray, np.ndarray]:
        """Gillespie simulation of one channel. Returns ``(times, states)``.

        ``times[i]`` is when the channel *entered* ``states[i]``.
        """
        rng = np.random.default_rng(seed)
        q = self.rate_matrix(tension)
        state = (int(rng.choice(4, p=self.steady_state(tension)))
                 if start_state is None else int(start_state))
        t = 0.0
        times, states = [0.0], [state]
        while t < duration:
            out = q[state].copy()
            out[state] = 0.0
            total = out.sum()
            if total <= 0:
                break
            t += rng.exponential(1.0 / total)
            if t >= duration:
                break
            state = int(rng.choice(4, p=out / total))
            times.append(t)
            states.append(state)
        return np.array(times), np.array(states, dtype=int)

    def mean_open_time(self, tension: float) -> float:
        """Expected dwell time in the open state, in seconds."""
        q = self.rate_matrix(tension)
        return float(-1.0 / q[O, O])

    # ------------------------------------------------------------- mutants

    def with_modification(self, label: str, **scale) -> "GatingModel":
        """Return a copy with named rate constants scaled.

        For example ``model.with_modification("R2456H", k2=0.39)`` slows entry
        into inactivation by the measured factor. Scaling a rate is a
        phenomenological stand-in for a structural change, and is labelled as
        such wherever it is displayed.
        """
        fields = {}
        mapping = {"k1": "k1_0", "km1": "km1_0", "k2": "k2", "km2": "km2",
                   "k3": "k3_0", "k4": "k4", "km4": "km4"}
        for name, factor in scale.items():
            if name not in mapping:
                raise ValueError(f"unknown rate {name!r}; "
                                 f"choose from {sorted(mapping)}")
            attr = mapping[name]
            fields[attr] = getattr(self, attr) * float(factor)
        return replace(self, label=label,
                       provenance=self.provenance + f" — modified: {scale}",
                       **fields)

    def calibrate_k2_for_tau(self, target_tau_s: float, tension: float = 5.0,
                             duration: float = 3.0, lo: float = 0.02,
                             hi: float = 8.0, tol: float = 1e-3) -> float:
        """Find the O→I₁ scale factor giving a target inactivation τ.

        Scaling ``k2`` by the ratio of measured time constants is the obvious
        thing to do and it is wrong: τ is a property of the whole four-state
        system, not of one rate, so a 2.6-fold slowing of the measured τ needs
        a rather different factor. This solves for it by bisection instead,
        which makes the mutant presets quantitatively faithful to the
        electrophysiology rather than merely qualitatively so.
        """
        def tau_of(scale: float) -> float:
            m = self.with_modification("probe", k2=scale)
            t = m.step(tension, duration=duration, n_points=3000).inactivation_tau()
            return float("inf") if t is None else t

        # tau decreases as k2 increases, so bracket accordingly.
        f_lo, f_hi = tau_of(lo), tau_of(hi)
        if not (min(f_lo, f_hi) <= target_tau_s <= max(f_lo, f_hi)):
            raise ValueError(
                f"target tau {target_tau_s * 1e3:.1f} ms is outside the "
                f"reachable range {min(f_lo, f_hi) * 1e3:.1f}-"
                f"{max(f_lo, f_hi) * 1e3:.1f} ms at {tension} mN/m")
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            t_mid = tau_of(mid)
            if abs(t_mid - target_tau_s) < tol * target_tau_s:
                return mid
            if (t_mid > target_tau_s) == (f_lo > target_tau_s):
                lo, f_lo = mid, t_mid
            else:
                hi, f_hi = mid, t_mid
        return 0.5 * (lo + hi)

    def mutant(self, name: str, tension: float = 5.0) -> "GatingModel":
        """Build a mutant model from :data:`MUTANT_PRESETS`, calibrated to τ.

        Where the preset gives a measured time constant, the O→I₁ rate is
        solved for rather than guessed. These remain *phenomenological* models
        of the mutant — a rate scaled to match an observation, not a structural
        prediction — and everything that displays them says so.
        """
        if name not in MUTANT_PRESETS:
            raise KeyError(f"unknown mutant {name!r}; "
                           f"choose from {sorted(MUTANT_PRESETS)}")
        preset = MUTANT_PRESETS[name]
        wt_tau = self.step(tension, duration=3.0,
                           n_points=4000).inactivation_tau()
        if wt_tau is None:
            raise RuntimeError("wild-type tau could not be measured")
        target = wt_tau * preset["tau_ratio"]
        scale = self.calibrate_k2_for_tau(target, tension=tension)
        model = self.with_modification(name, k2=scale)
        model.provenance = (
            f"{name}: phenomenological — O->I1 rate scaled by {scale:.3f} so "
            f"that tau is {preset['tau_ratio']:.2f}x the model's own wild-type "
            f"value, matching the measured fold change. {preset['note']} "
            f"[{preset['source']}]")
        return model


#: Phenomenological mutant and treatment parameterisations.
#:
#: ``tau_ratio`` is the measured **fold change** in inactivation time constant
#: relative to wild type, and the O→I₁ rate is solved for so the model
#: reproduces that ratio against *its own* wild-type value.
#:
#: Absolute time constants are deliberately not used as targets. The Young 2023
#: parameterisation gives a wild-type tau of 35-80 ms depending on tension,
#: whereas Bae's whole-cell measurement at -80 mV gives 8.6 ms — different
#: preparations and different stimulus protocols. Calibrating a mutant to an
#: absolute value from a different experiment made R2456H come out *faster*
#: than the model's wild type, which is the exact opposite of the biology.
#: Fold changes are the quantity that transfers between preparations.
#:
#: None of these is a structural prediction — each is a rate tuned to match an
#: observation, and every display of them says so.
WT_TAU_MS = _P.value("kinetics.wt_tau_ms")          # Bae et al. PNAS 2013, whole cell at -80 mV
MUTANT_PRESETS = {
    "R2456H": {
        "tau_ratio": 22.2 / 8.6, "classification": "GoF",
        "note": "DHS1 gain of function. Inactivation tau 22.2 +/- 2.1 ms vs "
                "8.6 +/- 0.4 ms wild type (2.6x slower).",
        "source": "Bae et al. PNAS 2013 (PMID 23487776); "
                  "Zarychanski et al. Blood 2012 (PMID 22529292)",
        "citation_key": "bae2013",
    },
    "M2225R": {
        "tau_ratio": 2.2, "classification": "GoF",
        "note": "DHS1 gain of function, slowed inactivation. The published tau "
                "is figure-embedded, so this target is approximate.",
        "source": "Albuisson et al. Nat Commun 2013 (PMID 23695678)",
        "citation_key": "albuisson2013",
    },
    "A2020T": {
        "tau_ratio": 2.0, "classification": "GoF",
        "note": "DHS1 gain of function, slowed inactivation. Approximate — the "
                "published value is figure-embedded.",
        "source": "Andolfo et al. Blood 2013 (PMID 24004669)",
        "citation_key": "andolfo2013",
    },
    "EPA-treated": {
        "tau_ratio": 20.0 / 34.0, "classification": "treatment",
        "note": "Eicosapentaenoic acid speeds inactivation (34 -> ~20 ms at "
                "100 uM) and normalises the R2456H phenotype to wild type.",
        "source": "Romero et al. Nat Commun 2019 (PMID 30867417)",
        "citation_key": "romero2019",
    },
}
