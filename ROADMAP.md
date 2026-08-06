# ROADMAP

Planned work, organised into ~20-minute rounds. Each round: implement, test,
fix, update docs, commit. Items are marked `[ ]` planned, `[~]` in progress,
`[x]` done. New improvements are appended after every fifth round.

**Status:** Round 3 of 5 (block A)

---

## The destination

Everything here builds toward one capability, which is the reason the project
exists:

> **Predict, from structure alone, whether a PIEZO1 variant is gain- or
> loss-of-function — and validate that prediction blind against the 68 curated
> variants whose phenotypes are already known.**

If elastic-network perturbation at a mutated residue systematically shifts the
predicted gating energetics in the direction of the measured phenotype, that is
a real, testable, publishable result. Every item below is either a component of
that pipeline or a tool for interrogating it.

The chain is: **structure → elastic network → gating coordinate → dome area
change → tension-dependent free energy → open probability → comparison with
measured P50 and inactivation kinetics.** Each round closes one link.

---

## Block A — complete the physics chain (rounds 1–5)

### Round 1 — Pore geometry  ✅
- [x] `structure/pore.py`: pore-radius profile along the conduction axis.
      Apollonius-clearance maximisation per slice, with the mandatory leash
      constraint (unconstrained, the probe sphere escapes to R ≈ 6188 Å).
- [x] Bottleneck detection and per-residue lining assignment.
- [x] Validate: the hydrophobic gate (human I2447/V2450/F2454) must fall at or
      near the measured constriction in a closed structure.
- [x] Tests; docs; commit.
- **Result:** closed human 8YEZ bottleneck **0.76 Å** (non-conductive); flat
  11ZC bottleneck **3.25 Å** (conductive). The profiler independently
  rediscovered all three curated constrictions — the V2450 hydrophobic gate at
  3.0 Å and the CTD constrictions at M2467 (1.2 Å) and P2510 (1.4 Å) — from
  coordinates alone, with no knowledge of the annotation.
- [ ] *Deferred to a later round:* expose the pore profile in the GUI.

### Round 2 — Gating kinetics  ✅
- [x] `physics/kinetics.py`: the Young et al. 2023 PNAS four-state model
      (C ⇌ O ⇌ I₁ ⇌ I₂), parameterised directly in tension.
      σ₅₀ = 1.4 mN/m, b = 0.8 mN/m, k₁ = 5.1·exp(σ/b), k₃ = 34.6·exp(−σ/b),
      k₋₃ by microscopic reversibility.
- [x] Deterministic solution (matrix exponential) and stochastic single-channel
      simulation (Gillespie).
- [x] Simulated macroscopic current traces and pressure–response curves.
- [x] Validate: recover the published P50 and inactivation time course.
- [x] Tests; docs; commit.
- **Result:** emergent half-activation **2.71 mN/m** against a measured
  cell-attached T50 of **2.7 ± 0.1 mN/m** (Lewis & Grandl 2015). Detailed
  balance exact to 1e-16. Mutants calibrated by *fold change*, reproducing
  R2456H at 2.58x wild-type τ.
- **Bug caught:** calibrating mutants to an *absolute* τ measured in a different
  preparation (Bae's 8.6 ms whole-cell) against a model whose own wild-type τ is
  35–80 ms made R2456H come out **faster** than wild type — the opposite of the
  biology. Fold changes transfer between preparations; absolute time constants
  do not. Now pinned by a test.

### Round 3 — Membrane mechanics
- [ ] `physics/membrane.py`: Monge-gauge Helfrich solver for the membrane
      footprint. 1-D radial solve validated against the exact modified-Bessel
      K₀ solution, then revolved.
- [ ] `physics/dome.py`: ΔE = −T·ΔA energetics, tension–area coupling, and the
      footprint contribution — which dominates tension sensitivity
      (λ = 14 nm, κ = 20 k_BT ⟹ γ = 0.42 mN/m).
- [ ] Validate against λ = √(κ/γ) and the published dome free energies.
- [ ] Tests; docs; commit.

### Round 4 — Experimental conformational space
- [ ] `analysis/ensemble.py`: superpose all structures onto a common residue
      basis and run PCA. The principal components are the conformational space
      the experiments actually sampled.
- [ ] Compare ensemble PCs with ANM modes — a strong, independent validation of
      the elastic network, and a map of which motions are real.
- [ ] Handle the protomer-correspondence and coverage traps across 21 entries.
- [ ] Tests; docs; commit.

### Round 5 — Allostery and force transmission
- [ ] `analysis/allostery.py`: Perturbation Response Scanning — apply unit
      forces at every residue, measure the response at the gate. Identifies
      sensor and effector residues from first principles.
- [ ] Dynamic cross-correlation from the ANM covariance; shortest-path
      allosteric pathway from blade tip to pore gate.
- [ ] Validate: the beam and anchor should appear on the dominant pathway, as
      the lever-like transduction model predicts.
- [ ] Tests; docs; commit.

---

## Block B — the variant pipeline and research tooling (rounds 6–10)

*Added after Round 5 review.*

### Round 6 — Variant impact prediction
- [ ] `analysis/variant_impact.py`: perturb the elastic network at a mutated
      residue (contact-weighted spring modification), recompute the low-frequency
      A-mode spectrum, and report the shift in the gating coordinate.
- [ ] Predicted ΔΔG of gating per variant.
- [ ] Tests; docs; commit.

### Round 7 — Blind validation
- [ ] Run the pipeline over all 68 curated variants.
- [ ] Ask whether GoF and LoF separate in the predicted direction. Report the
      effect size and a p-value honestly, including if the answer is no.
- [ ] `docs/VALIDATION.md` with the full result.
- [ ] Tests; docs; commit.

### Round 8 — Pockets and ligands
- [ ] `analysis/pockets.py`: Delaunay alpha-sphere pocket detection.
- [ ] Validate by recovering the Yoda1 pocket de novo and checking it against
      the mutagenesis-mapped residues (human A1718/A2075/A2078).
- [ ] Map resolved lipid densities (L9Q, PLX, P5S, PEE, D12) to contact residues.
- [ ] Tests; docs; commit.

### Round 9 — Conservation and constraint
- [ ] `analysis/conservation.py`: fetch orthologs, align, per-residue
      conservation, overlay on structure.
- [ ] Cross with variant density to find constrained regions with no reported
      variants — candidate untested functional sites.
- [ ] Tests; docs; commit.

### Round 10 — Research workflow
- [ ] Session save/load; analysis report export with full provenance.
- [ ] Headless CLI (`python -m piezo1.cli`) for batch analysis.
- [ ] Documented notebook API.
- [ ] Tests; docs; commit.

---

## Standing per-round checklist

1. Run the full test suite; fix anything red before adding features.
2. Run the scripted GUI smoke test — mechanical refactors of Qt code have
   broken it twice, silently.
3. Check no file exceeds 500 lines.
4. Update `INTERFACE.md` status, `SESSION_LOG.md` reasoning, `docs/SCIENCE.md`
   if a parameter or source changed, and this file's checkboxes.
5. Commit with a message explaining *why*, not just what.

---

## Deliberately not doing

- All-atom MD in the interactive loop. OpenMM is available for offline
  refinement, but interactive dynamics stay coarse-grained by design.
- Cooperativity between channels: measured P50 and open probability are
  invariant from 1 to 100 channels/µm² (Lewis & Grandl 2021), so modelling it
  would be inventing physics.
- Modelling the gain-of-function activation latency (344 ± 133 ms) inside the
  Markov scheme — the authors state explicitly that no Markov model reproduces
  it.
- Full C3 block-diagonalisation of the Hessian: benchmarked at only 1.76× on
  top of the sparse solve, not worth the complexity. Symmetry *labelling* gives
  the scientific payoff already.
