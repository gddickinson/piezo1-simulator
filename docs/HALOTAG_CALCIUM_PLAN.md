# Plan: HaloTag constructs, calcium permeation, and the two together

Written 2026-08-06 after reviewing `/Users/george/claude_test/halotag_binding_sim`.
This is a design document, not a record of work done. Nothing below has been
implemented.

---

## 1. What the HaloTag project already provides

`halotag_binding_sim` is a well-formed kinetic model of covalent HaloTag-ligand
labelling, with the same provenance discipline this project uses. Its structure:

| Module | What it gives us |
|---|---|
| `references.py` | Citation registry. Source experiment is **Bertaccini et al. 2025**, *Nat Commun* 16:5556 — PIEZO1-HaloTag in hiPSC-derived cells, using JF646-HTL and **JF646-BAPTA-HTL** |
| `parameters.py` | `k_on = 2.7 × 10⁶ M⁻¹s⁻¹`; cell conditions differing in `k_perm`, `partition`, `active_fraction` |
| `permeation.py` | Cumulative exposure `E(t) = partition·[L]·(t − (1−e^{−k_perm t})/k_perm)` |
| `kinetics.py` | Per-site `p(t) = active_fraction·(1 − e^{−k_on E(t)})` |
| `trimer.py` | `Binomial(3, p)`, fully-labelled `= p³`, predicted 1:2:3-dye brightness mixture |
| `stochastic.py` | Monte-Carlo population of individual channels |

**The model is complete and self-consistent, and it is entirely non-spatial.**
It knows a PIEZO1 trimer has three sites; it does not know where they are. That
is exactly the gap this project can fill, and the reason combining them is worth
doing rather than merely convenient.

The clean way to combine: **import the kinetics, do not reimplement it.** The
existing package becomes a dependency (or a vendored module with its provenance
intact), and this project supplies geometry, structure and the ion physics.

---

## 2. The four things requested, in order of tractability

### 2a. Put the HaloTags on the structure — *tractable, mostly bookkeeping*

PIEZO1's C-terminus is residue 2521 (human). Measured on 8YEZ, it sits at
**z = −60.4 Å on the conduction axis, 23.4 Å off-axis — only 8 Å below the CTD
constriction** at P2510/E2511. Three of them, one per protomer, related by C3.

HaloTag is 33 kDa (~300 aa). Structures are available and good: **6U32** (1.8 Å,
HaloTag with a TMR-HaloTag ligand bound — so the covalent conjugate is
resolved) and **6U2M** (2.0 Å, HaloCaMP, a HaloTag-based calcium indicator).

What has to be honest here: **the tag's position and orientation relative to
PIEZO1 are not experimentally determined.** There is no structure of the fusion.
Whatever is rendered is a model, and the linker is flexible, so the tag samples
an ensemble rather than sitting somewhere. The project already has the
convention for this — `structure/hybrid.py` records a seam and renders it
visibly — and the same must apply: the tag is drawn in a distinct style, its
placement is labelled a model, and an accessible-volume envelope is shown rather
than one arbitrary pose.

Deliverable: a fusion-model builder, a rendering style that does not pretend the
placement is measured, and an accessible-volume cloud for the linker.

### 2b. Dynamic binding on the actual trimer — *tractable*

The kinetics already exist. What this project adds is showing them on the
structure: three sites, each stochastically labelled over time, the trimer
brightness stepping 0 → 1 → 2 → 3 dyes.

The scientific point this makes visible is the one the source paper raises.
`p = 0.9` per site gives only `0.9³ = 0.73` of channels fully labelled, so a
population of puncta shows a discrete brightness mixture. Rendering that on a
real trimer, alongside the concentration-time sweep, connects a labelling
protocol to what a microscope actually sees.

Deliverable: a labelling controller driving per-site occupancy from
`halotag_sim.kinetics`, an animation, and a brightness-histogram readout
comparable with the paper's all-points amplitude histograms.

### 2c. Calcium permeation, open versus closed, across variants — *tractable if scoped honestly*

This is where the temptation to overreach lives. All-atom MD of permeation is
explicitly a non-goal in `CLAUDE.md`, and rightly: it would not be interactive
and it would not be validated.

The defensible approach is **1-D Poisson–Nernst–Planck electrodiffusion along
the conduction axis**, using the pore radius profile this project already
measures. That is a standard, citable treatment for channel conductance, it runs
in milliseconds, and — critically — it has a **validation target already in the
ground-truth table**: unitary conductance **25–30 pS** (Shi 2020; Vaisey &
MacKinnon 2026).

It also connects directly to Round 19. A pore that the wetting heuristic calls
dewetted should carry no current, and the mechanism is not steric — it is that
liquid water, and therefore the hydration shell an ion needs, is absent. So the
model should be: **conductance from PNP through the radius profile, gated by the
wetting verdict.** Closed 8YEZ is blocked twice over (0.095 nm bottleneck *and*
Σd = 0.82); flat 11ZC is open on both counts. Reporting which mechanism blocks
which structure is more informative than a single number.

Across variants: honest scope is that only variants with deposited structures
can be done this way — 8ZU3, 8ZU8, 8YFC, 8YFG carry disease variants. For the
other 60-odd curated variants there is no structure, and the project should say
so rather than extrapolating. This is the same data limitation Round 22 ran into
from the other direction.

Validation targets, all already cited in the project:
- unitary conductance 25–30 pS
- P_Ca/P_Na ≈ 3–6 (PIEZO1 is a non-selective cation channel with modest Ca²⁺
  preference)
- closed structures must give ~0 pS

Deliverable: a `physics/permeation.py` implementing 1-D PNP over
`PoreProfile`, gated by `WettingPrediction`, validated against 25–30 pS, plus a
particle animation whose flux is *set by* the computed current rather than
chosen for looks.

### 2d. Calcium meeting the HaloTag — *tractable, and the most interesting*

JF646-BAPTA-HTL is a Ca²⁺-sensitive ligand. So the question is: what [Ca²⁺] does
a tag at the C-terminus actually see when the channel opens?

This is the classic **calcium nanodomain** calculation (Neher 1998; Stern 1992):
linearised buffered diffusion from a point source,

```
[Ca²⁺](r) = i_Ca / (4π F D_Ca r) · exp(−r/λ),    λ = √(D_Ca / (k_on^B [B]))
```

Two things make this a good fit for this project rather than a bolt-on. First,
it is **the same mathematical form as the membrane footprint** already
implemented in `physics/membrane.py` — a screened Green's function with a decay
length set by a ratio of material constants. Second, it needs exactly the two
numbers the other pieces produce: `i_Ca` from §2c, and `r` from §2a.

**A prediction falls out immediately, and it is testable.** With γ = 30 pS at
−80 mV, a 15% Ca²⁺ fraction gives i_Ca ≈ 0.36 pA. At the tag's distance from the
pore mouth — 4–6 nm, given the C-terminus sits 8 Å below the CTD constriction
and HaloTag's centre is another 25–45 Å out — the nanodomain concentration is of
order **200 µM**. BAPTA's Kd for Ca²⁺ is **~0.2 µM**, a thousandfold lower.

So the model predicts that **a JF646-BAPTA HaloTag on the PIEZO1 C-terminus is
saturated whenever its own channel is open.** It reports opening as a binary
event, not as graded calcium. If true, that changes how the published puncta
traces should be read: brightness steps would reflect *how many tags are
labelled and how often the channel opens*, not local calcium amplitude.

That is a real, falsifiable, experimentally relevant claim, it comes from
combining the two projects, and neither could produce it alone. It is also the
kind of claim this project must state with its uncertainty attached — the tag
distance is modelled, the Ca²⁺ fraction is uncertain, and buffering in an hiPSC
cytosol is not well constrained. A sensitivity sweep over all three is part of
the deliverable, not an afterthought.

---

## 3. Where this could go wrong, recorded in advance

1. **The tag position is a model, not a measurement.** Every number that depends
   on `r` inherits that. The sensitivity sweep must span 2–20 nm, not report a
   point estimate.
2. **PNP is a mean-field theory in a pore two ions wide.** It gets conductance
   about right and it does not resolve single-ion energetics. If it reproduces
   25–30 pS that is encouraging; it does not make it a permeation mechanism.
3. **Variant coverage will be poor, again.** Four structures against sixty-eight
   variants. Round 22 already established that the binding constraint on this
   project is structural and phenotypic data, not modelling.
4. **The wetting gate is a heuristic**, with AUROC 0.91 on a training set of
   ~200 channels. Using it as an on/off switch for current is a stronger claim
   than it was validated for, and the coupling should be reported as a factor
   rather than hidden inside a conductance.
5. **Do not let the animation outrun the physics.** A calcium particle animation
   is easy to make look impressive and easy to make meaningless. The flux must
   be set by the computed current, and the HUD must say what the frame rate
   corresponds to in real time — the same discipline the morph clock already
   applies.

---

## 4. Proposed rounds

Sequenced so each has an independent validation target.

- **Round 31 — HaloTag fusion geometry.** Fetch 6U32; build the fusion model at
  each C-terminus; accessible-volume envelope for the linker; render with a
  visible seam. *Validate:* C3 symmetry of the three tags preserved; tag centre
  4–6 nm from the pore exit; no steric clash with the CTD.
- **Round 32 — Labelling on the structure.** Import `halotag_sim` kinetics;
  per-site stochastic occupancy; brightness animation and histogram.
  *Validate:* reproduce the source project's `p³` curve and the 1:2:3 mixture
  exactly — a divergence means the import went wrong.
- **Round 33 — Calcium permeation.** 1-D PNP over `PoreProfile`, gated by the
  wetting verdict. *Validate:* **25–30 pS** for an open structure, ~0 for
  closed, and report which mechanism blocks which.
- **Round 34 — Variant permeation.** Apply to the four variant structures.
  *Validate:* direction of change against measured phenotype where published,
  and state coverage honestly.
- **Round 35 — The calcium nanodomain at the tag.** Buffered-diffusion Green's
  function; sensitivity over distance, Ca²⁺ fraction and buffering.
  *Validate:* the saturation prediction, with intervals, and a clear statement
  of what would falsify it.
