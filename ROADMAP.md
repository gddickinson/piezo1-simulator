# ROADMAP

Planned work, organised into ~20-minute rounds. Each round: implement, test,
fix, update docs, commit. Items are marked `[ ]` planned, `[~]` in progress,
`[x]` done. New improvements are appended after every fifth round.

**Status:** Round 5 of 5 (block A)

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

### Round 3 — Membrane mechanics  ✅
- [x] `physics/membrane.py`: Monge-gauge Helfrich solver for the membrane
      footprint. 1-D radial solve validated against the exact modified-Bessel
      K₀ solution, then revolved.
- [x] `physics/dome.py`: ΔE = −T·ΔA energetics, tension–area coupling, and the
      footprint contribution — which dominates tension sensitivity
      (λ = 14 nm, κ = 20 k_BT ⟹ γ = 0.42 mN/m).
- [x] Validate against λ = √(κ/γ) and the published dome free energies.
- [x] Tests; docs; commit.
- **Result:** λ = 14 nm with κ = 20 k_BT gives γ = **0.420 mN/m** (published
  0.42), and the unit conversion comes out at **4.116 mN/m per k_BT/nm²**
  (published 4.114). The solver is second-order convergent against the exact
  K₀ solution and recovers λ = **13.998 nm** from its own profile.
  Cox et al.'s ΔG₀ = 9.7 k_BT and ΔA = 8 nm² reproduce **T₅₀ = 4.99 mN/m**
  against their measured 5.1 ± 0.2 mN/m.
- **The footprint does dominate**, as Haselwandter & MacKinnon argued: around
  the measured 7WLT dome it stores **622 nm² of excess area against the dome's
  256 nm²** — 243% on top.
- **Two things had to be fixed, and one flagged.** Building the biharmonic
  operator as `L @ L` squares the Laplacian's condition number: it converged to
  a 47 nm decay length where the answer is 14 nm, with a 59% energy error that
  did not improve on refinement. Rewritten as a coupled second-order system.
  And the closed-form energy had the Bessel ratio inverted (K₁/K₀ rather than
  K₀/K₁), which is 2.5x too large at PIEZO1's r₀/λ — caught only by integrating
  the functional over the exact profile.
  **Flagged:** PIEZO1's measured contact slope is ~2.0 (63°), far outside the
  small-slope regime the linearised Monge gauge assumes. The solution reports
  this rather than presenting the numbers as quantitative.

### Round 4 — Experimental conformational space  ✅
- [x] `analysis/ensemble.py`: superpose all structures onto a common residue
      basis and run PCA. The principal components are the conformational space
      the experiments actually sampled.
- [x] Compare ensemble PCs with ANM modes — a strong, independent validation of
      the elastic network, and a map of which motions are real.
- [x] Handle the protomer-correspondence and coverage traps across 21 entries.
- [x] Tests; docs; commit.
- **Result — the strongest validation so far.** PCA over 10 mouse PIEZO1
  structures on a shared 1091-residue basis gives **PC1 = 90.0% of variance**,
  and PC1 orders every structure by gating state with no labels supplied:
  seven curved entries all negative, the 8IXO intermediate at +334, flattened
  7WLU at +678, flat 11ZC at +1045. The dominant axis of experimental
  variation *is* the gating transition.
- **PC1 overlaps 0.804 with ANM mode 6**, cumulative 0.960 over 30 modes, and
  all of the top three PCs match **A**-symmetric modes despite E modes
  outnumbering A two to one. RWSIP 0.555 against a random-vector control of
  0.001. The symmetry selection rule now shows up in the deposited structures,
  not just in one pairwise transition.
- **Four traps handled, each of which silently returns a number rather than
  erroring:** cross-species numbering (mouse converted to human, offsets not
  constant); coverage (all 20 entries share only 325 residues, so poor-coverage
  entries are dropped and reported); protomer correspondence (four entries are
  labelled in reversed rotational order); and paralogues — **6KG7 is PIEZO2**
  and is excluded by default.
- **6LQI excluded as well**, and this mattered: it is the Piezo1.1 splice
  isoform missing residues 1382–1405, so its difference from the rest is a
  *sequence* difference, not a conformational one. Included, it dominates a
  whole component on its own and splits the gating coordinate across PC1 (58%)
  and PC2 (36%); excluded, PC1 is a single clean 90%.

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

## Block C — measurement tools  *(added on request, 2026-08-05)*

### Round 11 — Geometric measurement toolkit  ✅
- [x] `analysis/measure.py`: distance, angle, dihedral between picked atoms or
      residue centroids; radius of gyration; end-to-end and domain–domain
      distances; helix tilt and crossing angles; per-residue RMSF from a mode
      set or trajectory; principal-axis and inertia tensor.
- [x] Solvent-accessible surface area (Shrake–Rupley) per atom, residue and
      domain; buried surface area at interfaces.
- [x] Pore hydrophobicity profile alongside the radius profile — CHAP's key
      insight is that radius alone does not predict conduction; a wide but
      hydrophobic constriction can still dewet and block.
- [x] Tests; docs; commit.
- **Result:** independently reproduced two structural facts. TM38, the
  pore-lining helix, is the least tilted of the pore-proximal helices at
  **6.9°** from the three-fold axis. Note the claim had to be narrowed: blade
  helices out at 50–60 Å radius are *also* near-vertical (TM30 is 3.1°), so
  "least tilted overall" is false and "least tilted of the pore module" is what
  the data supports.

### Round 12 — Interaction detection  ✅
- [x] `analysis/interactions.py`: hydrogen bonds, salt bridges, hydrophobic
      contacts, π-stacking, cation–π and disulfides, with explicit published
      geometric criteria rather than a single distance cutoff.
- [x] Protein–ligand interaction profiling for the resolved lipids
      (L9Q, PLX, P5S, PEE, D12) and any docked compound.
- [x] Interaction *changes* between two states — which contacts break when the
      dome flattens is exactly the mechanotransduction question.
- [x] Tests; docs; commit.
- **Result:** the UniProt-annotated C2411–C2415 disulfide is recovered in
  all three protomers at 2.04 Å, from coordinates alone. And a finding that
  matters for the destination: **R2456 forms an inter-protomer salt bridge with
  E2117** — the selectivity glutamate — at 3.66–3.91 Å in all three protomers.
  R2456H is the archetypal gain-of-function variant, so this is a concrete
  structural route from mutation to phenotype.
- **Two criteria had to be corrected.** PLIP's 4.1 Å hydrogen-bond cutoff is
  only valid *with* hydrogens and an angle check; applied to heavy atoms alone
  it gave 8005 "bonds" per trimer. Tightened to the conventional 3.5 Å. And
  N···N pairs were being admitted, which is donor–donor and impossible; now
  excluded except for histidine, whose ring nitrogens may be unprotonated.

### Round 13 — Measurement in the GUI
- [ ] Click-to-measure tool: pick two atoms for a distance, three for an angle,
      four for a dihedral; persistent labelled measurements in the 3D view.
- [ ] Measurement panel listing everything measured, with copy-to-clipboard
      and CSV export.
- [ ] Tests; docs; commit.

---

## Block D — animation  *(added on request, 2026-08-05)*

### Round 14 — Animation engine  ✅
- [x] `render/animation.py`: a timeline that interpolates camera, coordinates,
      colours and visibility; frame-accurate offscreen capture.
- [x] GIF and MP4 export (imageio/ffmpeg), with a legibility-first default:
      slow ease-in-out, a held first and last frame, and an on-frame caption
      stating what is being shown.
- [x] Tests; docs; commit.

### Round 15 — The animation library  ✅
- [x] Normal-mode animations, per mode, with porcupine displacement arrows.
- [x] The gating morph, with the dome profile and pore radius plotted live
      alongside the structure.
- [x] Tension-driven gating: tension ramp → dome flattening → pore opening →
      simulated current trace, all in one synchronised animation.
- [x] Ligand and lipid interaction animations: the Yoda1 pocket, the PIP2
      lysine cluster, and the pore lipid that contacts R2456.
- [x] Variant comparison animations: wild type against a gain-of-function
      mutant, side by side on the same clock.
- [x] Tests; docs; commit.
- **Result:** seven animations render offscreen. MP4 is ~10x smaller than GIF
  for the same content (3.8 MB vs 34 MB for the gating morph at full size), so
  it is the recommended format. `docs/anim/` is git-ignored and regenerable,
  consistent with the rule that nothing generated is committed; one small GIF
  is kept for the README.

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
