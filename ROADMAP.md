# ROADMAP

Planned work, organised into ~20-minute rounds. Each round: implement, test,
fix, update docs, commit. Items are marked `[ ]` planned, `[~]` in progress,
`[x]` done. New improvements are appended after every fifth round.

**Status:** Block A complete (rounds 1–5). Blocks B–E outstanding.

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

### Round 5 — Allostery and force transmission  ✅
- [x] `analysis/allostery.py`: Perturbation Response Scanning — apply unit
      forces at every residue, measure the response at the gate. Identifies
      sensor and effector residues from first principles.
- [x] Dynamic cross-correlation from the ANM covariance; shortest-path
      allosteric pathway from blade tip to pore gate.
- [x] Validate: the beam and anchor should appear on the dominant pathway, as
      the lever-like transduction model predicts.
- [x] Tests; docs; commit.
- **Result:** the anchor domain sits **on the optimal route** — forcing the
  blade→gate path through it costs a detour penalty of −0.000, i.e. it is
  already there — and it ranks second by path betweenness (5.19) behind only
  the CTD (7.67). The cap is not a transmission route (+0.055 penalty).
- **The beam result is more nuanced than the prediction.** It does *not* appear
  on the single shortest path, but forcing the route through it costs only
  **+0.010**, so it is a near-degenerate parallel channel rather than an
  excluded one. Its betweenness is real but low (1.30). Reported as measured
  rather than rounded towards the lever model.
- **An error worth recording:** asking "does the path go through X" by computing
  source→X and X→target separately and adding is wrong. Each leg picks its own
  best endpoints, which on a C3 trimer can be in *different protomers*, so the
  legs never join. Done that way the detour came out **cheaper** than the
  unconstrained shortest path, which is impossible. `detour_cost()` now shares
  a single via-point, and an invariant test asserts a constrained path can
  never beat the free one.

---

## Block B — the variant pipeline and research tooling (rounds 6–10)

*Added after Round 5 review.*

### Round 6 — Variant impact prediction  ✅
- [x] `analysis/variant_impact.py`: perturb the elastic network at a mutated
      residue (contact-weighted spring modification), recompute the low-frequency
      A-mode spectrum, and report the shift in the gating coordinate.
- [x] Predicted ΔΔG of gating per variant.
- [x] Tests; docs; commit.
- [x] **`docs/PREREGISTRATION.md` written now, before Round 7** — pulling Round
      20's protocol forward, because a test is only blind if the rule was fixed
      first.
- **Method:** ΔΔG_gating = ½·dᵀ(H_mut − H_wt)·d, the change in elastic cost of
  the *observed* gating motion. Since H_mut − H_wt is non-zero only at the
  mutated residue's contacts, this costs O(contacts) rather than a
  re-diagonalisation, and it is **exact to 7e-16** against an explicitly
  rebuilt Hessian.
- **Coverage:** 48 of 68 variants scored; 17 residues fall outside the resolved
  range and 3 are not single substitutions. Reported, not silently dropped.
- **No phenotype comparison was performed in this round**, deliberately. The
  ΔΔG distribution is balanced a priori (23 stiffening, 21 softening), which is
  what an untuned predictor should look like.

### Round 7 — Blind validation  ✅  **NULL RESULT**
- [x] Run the pipeline over all 68 curated variants.
- [x] Ask whether GoF and LoF separate in the predicted direction. Report the
      effect size and a p-value honestly, including if the answer is no.
- [x] `docs/VALIDATION.md` with the full result.
- [x] Tests; docs; commit.
- **Result: H0 not rejected.** 25 variants survived the pre-registered
  inclusion criteria (16 GoF, 9 LoF). One-sided permutation test **p = 0.234**;
  Cliff's delta **−0.083**, CI [−0.528, +0.403] (negligible, spans zero);
  **AUROC 0.542**. The mean difference points the predicted way but the effect
  is negligible and the ranking is barely above chance.
- **Post-hoc diagnostic explains why:** only **0.2% of ΔΔG variance is
  within-position**. The score is dominated by *where a residue sits*, not
  *which substitution occurred* — ΔΔG scales with local gating strain and
  contact count, both properties of the position, while the substitution enters
  through a single scalar. Four variants at R2456 spanning GoF and LoF all get
  "softening", the largest belonging to the **LoF** one.
- **This does not invalidate the physics chain**, every link of which was
  validated against an independent published number. It means a single scalar
  from a volume-scaled elastic network is not sufficient to call a phenotype at
  this sample size.
- [x] `docs/VALIDATION.md` written in the pre-registered order, with the null
      result stated first.

### Round 8 — Pockets and ligands  ✅
- [x] `analysis/pockets.py`: Delaunay alpha-sphere pocket detection.
- [x] Validate by recovering the Yoda1 pocket de novo and checking it against
      the mutagenesis-mapped residues (human A1718/A2075/A2078).
- [x] Map resolved lipid densities (L9Q, PLX, P5S, PEE, D12) to contact residues.
- [x] Tests; docs; commit.
- **Recovered de novo:** the transmembrane hydrophobic gate (2/3 residues) and
  the anchor-domain apex brake (2/2), from geometry alone.
- **The Yoda1 site is an interfacial groove, not an enclosed cavity.** Searching
  for enclosed cavities finds at most 1 of its 3 residues; allowing surface
  grooves finds 2. That is consistent with Yoda1 acting as a wedge from the
  lipid phase, with a PLX lipid occupying part of the site in 7WLT, and with
  this project's own annotation labelling the site's evidence as *predicted* —
  it has never been seen in a co-structure. Reported as a nuanced negative
  rather than forced into a positive.
- **A percolation trap, fixed.** With a radius filter alone, single-linkage
  merged PIEZO1's whole exterior into one "pocket" of **408 000 Å³ with 601
  lining residues** — the protein's outside surface. A per-sphere burial filter
  (≥30 atoms within 8 Å) plus a tighter r_max brings the largest pocket to
  6 691 Å³ / 63 residues. Parameters were chosen on pocket-size plausibility,
  **before** checking any site recovery, to avoid tuning to the answer.

### Round 9 — Conservation and constraint  ✅
- [x] `analysis/conservation.py`: fetch orthologs, align, per-residue
      conservation, overlay on structure.
- [x] Cross with variant density to find constrained regions with no reported
      variants — candidate untested functional sites.
- [x] Tests; docs; commit.
- **62 vertebrate orthologs**, one per species, reference-anchored to human
  numbering. Mean conservation 0.770 over well-covered positions; 594 invariant.
- **An independent confirmation of Round 5.** Ranked by mean conservation, the
  **anchor domain is the most constrained of all (0.987)**, ahead of the inner
  helix (0.980) and CTD (0.960), while the distal blade THU1 is the least
  (0.719). Round 5 identified the anchor as the force-transmission hub from
  mechanics alone; evolution agrees, by a completely separate line of evidence.
- Annotated sites score as they should: anchor brake **1.000** (invariant
  across all 62 species), selectivity glutamates 0.986, PIP2 cluster 0.986,
  hydrophobic gate 0.934. The Yoda1 pocket is the *least* conserved at 0.859
  (A2075 only 0.63) — consistent with a synthetic agonist acting at a site not
  under strong selection.
- **Conservation alone is too blunt to be a hypothesis:** 426 positions are
  invariant, never mutated in the literature, and structurally resolved — a
  quarter of the protein. `rank_candidates()` therefore crosses it with the
  Round 5 mechanical coupling, which neither could do alone.
- **Nominated for testing:** residues **2021 and 2034** are invariant across 62
  species, carry no reported variant, and lie *on the blade→gate allosteric
  path* computed in Round 5. Of the top 40 distal candidates, 20 are in the
  anchor. Written to `data/derived/conservation_candidates.json`.

### Round 10 — Research workflow  ✅
- [x] Session save/load; analysis report export with full provenance.
- [x] Headless CLI (`python -m piezo1.cli`) for batch analysis.
- [x] Documented notebook API.
- [x] Tests; docs; commit.
- **`python -m piezo1.cli`** with `list`, `dome`, `pore`, `modes`, `pockets`,
  `interactions`, `variants`, `conservation`, `report` and `batch`. The CLI and
  the report share one analysis registry, so they cannot diverge — a test
  asserts every analysis is reachable from both.
- **Batch across all 20 structures in one command.** The curved entries cluster
  at R_c = 9.3–12.5 nm against a published 10.2; 8IXO (intermediate) sits at
  16.5; 11ZC (flat) at 21.6 and is the only entry called conductive. The run
  independently flagged **3JAC as an outlier** (R_c 5.3 nm, spuriously
  conductive) — the same entry the ensemble analysis excluded for poly-UNK
  numbering, found again by a different route.
- Provenance on every report: software version, input file, parameters, library
  versions, and any warnings. Sessions store *what you were looking at*, never
  coordinates or results, so a saved file cannot drift out of step with the code.
- `docs/NOTEBOOK.md` documents the headless API, including a "things that will
  bite you" table.

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

### Round 13 — Measurement in the GUI  ✅
- [x] Click-to-measure tool: pick two atoms for a distance, three for an angle,
      four for a dihedral; persistent labelled measurements in the 3D view.
- [x] Measurement panel listing everything measured, with copy-to-clipboard
      and CSV export.
- [x] Tests; docs; commit.
- **Validated against an independent measurement:** picking the two sulfurs of
  C2411/C2415 in the live GUI returns **2.04 Å**, exactly what
  `analysis.interactions` finds for that disulfide by a completely different
  code path. That check is now part of the scripted GUI smoke test.
- **A Qt trap worth recording.** Drawing labels with `QPainter` inside
  `paintGL` alongside the moderngl render does not work — the two fight over GL
  state and the text simply never appears, with no error. Labels are now drawn
  by a transparent child widget with its own paint event, which cannot
  interfere with the scene at all.
- Measurements are cleared when a new structure loads: atom indices are
  per-structure, so a retained measurement would silently point at different
  atoms.

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

## Block E — sharpening the prediction  *(added after the Block A review, 2026-08-06)*

**Where Block A leaves us.** The physics chain is closed: structure → elastic
network → gating coordinate → dome and footprint energetics → open probability,
with every link validated against a published number. Two results are strong
enough to build on — PC1 of the experimental ensemble *is* the gating
coordinate and matches an A-symmetric mode at 0.804, and the anchor domain is
the dominant force-transmission hub. What is still missing for the destination
is a per-variant number and an honest test of it.

### Round 16 — Allostery-derived variant features
- [ ] Per-residue features from the machinery now in place: PRS response at the
      gate, path betweenness, DCC coupling to the gate, mode-6 displacement
      amplitude, burial (SASA), and domain membership.
- [ ] These are the structure-derived predictors that feed Round 6, and they
      exist now only because Round 5 built them.
- [ ] Tests; docs; commit.

### Round 17 — External predictors, licence-clean
- [ ] Integrate the **ProtVar API** (`https://www.ebi.ac.uk/ProtVar/api`,
      **CC BY 4.0**), which serves AlphaMissense, EVE, ESM-1b, conservation and
      **precomputed FoldX ΔΔG** for a UniProt accession and position.
- [ ] This sidesteps the licensing traps entirely: FoldX itself is not
      redistributable, SIFT4G is GPL-3.0 copyleft, and every tool on
      biosig.lab.uq.edu.au carries no licence at all. VarSite and VarMap are
      both retired. Cache responses locally; degrade gracefully offline.
- [ ] Tests; docs; commit.

### Round 18 — Nonlinear membrane mechanics
- [ ] Round 3 flagged that PIEZO1's contact slope of ~2.0 (63°) is far outside
      the small-slope regime the Monge gauge assumes. Implement an axisymmetric
      **Euler–elastica / full-curvature** solve so the footprint energies
      become quantitative rather than indicative.
- [ ] Validate the nonlinear solver against the linear one in the small-slope
      limit, where they must agree.
- [ ] Tests; docs; commit.

### Round 19 — Pore hydration and wetting
- [ ] Radius alone does not predict conduction: a wide but hydrophobic neck can
      dewet and block. Combine the pore radius profile with the hydrophobicity
      profile into a CHAP-style conduction prediction.
- [ ] Check it calls closed 8YEZ non-conductive and flat 11ZC conductive for
      the *right reason*, not just by bottleneck radius.
- [ ] Tests; docs; commit.

### Round 20 — Statistical rigour for the blind test
- [ ] Before Round 7 reports anything: permutation test for the GoF/LoF
      separation, leave-one-out cross-validation, effect size with a confidence
      interval, and a pre-registered decision rule.
- [ ] **Write the negative-result protocol first**, so that "the predictor does
      not separate them" is a publishable outcome recorded in
      `docs/VALIDATION.md` rather than a prompt to keep tuning until it does.
- [ ] Tests; docs; commit.

---

## Block F — closing the loop  *(added after the Block B review, 2026-08-06)*

**Where Block B leaves us.** The blind test came back null, and its diagnostic
was precise: the mechanical predictor reports *where a residue sits*, not
*which substitution occurred*. Meanwhile Round 9 showed that conservation
crossed with mechanics is sharp where either alone is blunt, and Round 5 and
Round 9 independently converged on the anchor. The engine is now well ahead of
the interface: the GUI can show a structure, a dome measurement and normal
modes, but cannot reach the pore profile, pockets, conservation, allostery or
any of the reporting.

### Round 21 — Let the GUI reach the engine
- [ ] Expose the pore profile (with the hydrophobicity trace alongside it),
      pockets, conservation colouring and the allostery maps in the interface.
- [ ] Colour-by-conservation and colour-by-PRS as first-class colour schemes.
- [ ] Session save/load wired to the window; report export from a menu.
- [ ] Tests; docs; commit.

### Round 22 — A second, stated hypothesis for the variant test
- [ ] Write a *new* pre-registration before touching the labels again. The
      Round 7 result stands as recorded and is not to be revised.
- [ ] Combine substitution-aware evidence (ProtVar: AlphaMissense, EVE, ESM-1b,
      precomputed FoldX ΔΔG, all CC BY 4.0) with the mechanical and
      conservation features this project supplies and they do not.
- [ ] Report against the Round 7 baseline honestly, including if the addition
      does not help.
- [ ] Tests; docs; commit.

### Round 23 — Packaging and one-command reproduction
- [ ] `pyproject.toml`, a pinned environment lock, and a single
      `make reproduce` that fetches data, rebuilds every resource, runs the
      suite and regenerates every figure and number in the docs.
- [ ] Tests; docs; commit.

### Round 24 — Performance
- [ ] Profile the slow paths. Pocket detection is ~10 s, PRS builds an N×N
      matrix, and the ensemble PCA reloads every structure. None is fatal, all
      are avoidable.
- [ ] Tests; docs; commit.

### Round 25 — The teaching layer
- [ ] Project aim A1 is that this be a *learning* instrument, and it has had the
      least attention of any aim. Add a guided tour that walks the mechanism —
      dome, blades, lever, gate — with each step tied to the live measurement
      it corresponds to.
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
