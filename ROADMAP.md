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
- ~~**The footprint does dominate**, as Haselwandter & MacKinnon argued: around
  the measured 7WLT dome it stores **622 nm² of excess area against the dome's
  256 nm²** — 243% on top.~~
  **⚠ OVERTURNED BY ROUND 18.** The 622 nm² is linear-theory output at a 63°
  contact slope, where the theory does not apply; the nonlinear value is
  **179 nm²**, i.e. **0.70× the dome, not 2.4×**. The comparison was also not
  like for like — the dome's 256 nm² is an exact area, the footprint's was
  linearised. Left visible rather than edited away, because the caveat *was*
  recorded at the time and the lesson is that a caveat is not a substitute for
  doing the calculation.
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

### Round 16 — Allostery-derived variant features  ✅
- [x] Per-residue features from the machinery now in place: PRS response at the
      gate, path betweenness, DCC coupling to the gate, mode-6 displacement
      amplitude, burial (SASA), and domain membership.
- [x] These are the structure-derived predictors that feed Round 6, and they
      exist now only because Round 5 built them.
- [x] Tests; docs; commit.
- **`analysis/features.py`**: 1279 residues × 11 features for 8YEZ in 9 s —
  PRS gate response and coupling, path betweenness, DCC to the gate, amplitude
  along the lowest **A**-symmetric mode, mean-square fluctuation, relative
  SASA, conservation, distance to gate and axis, contact count, domain.
- **Structural validation only. No phenotype comparison was run**, deliberately;
  the Round 7 null result stands and re-testing needs the new pre-registration
  in Round 22. Checks used instead: gate response falls off with distance to
  the gate (r = −0.55), relative SASA lies in [0, 1], the gating mode used is
  confirmed symmetric, and the conservation join reproduces Round 9's domain
  ranking exactly (anchor highest at 0.987) — which it would not if the residue
  join were off by even one.
- **Two redundant features found and removed**, by a test that forbids any pair
  of columns correlating above 0.99. The PRS response matrix is *symmetric*
  (‖C_ij‖_F = ‖C_ji‖_F, the Frobenius norm being transpose-invariant), so
  "effectiveness" and "sensitivity" are literally the same numbers — effector
  and sensor coincide under this formulation. Row normalisation appears to
  separate them but makes effectiveness near-constant (spread/mean 0.0018) and
  leaves the normalised sensitivity correlating with raw coupling at **0.998**.
  Shipping both would have looked like two independent lines of evidence and
  been one. The table now carries a single honestly named `prs_coupling`.

### Round 17 — External predictors, licence-clean ✅
- [x] Integrated the **ProtVar API** (`https://www.ebi.ac.uk/ProtVar/api`,
      **CC BY 4.0** — confirmed from the service's own OpenAPI `info.license`,
      not assumed) in `piezo1/analysis/external.py`. Serves AlphaMissense, EVE,
      ESM-1b, conservation and **precomputed FoldX ΔΔG** per accession/position.
- [x] Sidesteps the licensing traps: FoldX itself is not redistributable,
      SIFT4G is GPL-3.0 copyleft, everything on biosig.lab.uq.edu.au carries no
      licence at all, VarSite and VarMap are retired. Responses cache to
      `data/cache/protvar/`; an offline client returns `None` rather than
      raising, so a missing score weakens an analysis instead of aborting it.
- [x] **Live annotation of the curated variants: 64/65 single substitutions in
      77 s.** Coverage — conservation 64, AlphaMissense 51, EVE 51, ESM-1b 51,
      FoldX ΔΔG 50. The 13 without missense scores are nonsense/frameshift
      variants, where a missense predictor correctly has nothing to say.
- [x] **Cross-validation, unplanned and worth recording: 0 wild-type mismatches
      out of 64.** ProtVar reports the wild-type residue it holds for each
      position, so this is an *external* confirmation that every variant in our
      table is numbered correctly against Q92508 — an independent check on the
      numbering work from Rounds 1–2 that we had no other way to make.
- [x] **A correctness trap found and closed.** `/score/{acc}/{pos}` returns
      nineteen entries per predictor — one per possible substitution — and the
      payload contains **no field saying which is which**. Reading them in
      array order, or assuming alphabetical, would have silently attributed the
      wrong pathogenicity to every variant. The undocumented `mt=` query
      parameter resolves it; a position-only query now keeps *only*
      conservation, which is genuinely position-level. `/prediction/foldx/`
      needs no such care — it labels each entry with `mutatedType`.
      (Also found: the documented `/prediction/interaction/` endpoint 404s.)
- [x] **Recorded limitation.** AlphaMissense, EVE and ESM-1b each emit a single
      *pathogenicity* axis, which by construction cannot express **direction**.
      Observed concretely: all four R2456 substitutions score PATHOGENIC, yet
      R2456H/K/P are gain-of-function and R2456C is loss-of-function. These
      predictors are therefore complementary to the mechanical features, not a
      replacement — which is exactly the hypothesis Round 22 must pre-register.
- [x] 11 tests (`tests/test_external.py`), all running **offline from the
      cache**, so the suite needs no network and does not hammer a public
      service. Suite 243 → **254 passed**; GUI smoke test clean.
- [x] **No phenotype comparison was run.** Round 7's null result stands as
      recorded; re-testing against the 68 labels requires Round 22's new
      pre-registration to be written *first*.

### Round 18 — Nonlinear membrane mechanics ✅
- [x] Implemented `piezo1/physics/elastica.py`: an axisymmetric **Euler–elastica**
      solve in arc-length parametrisation, with exact principal curvatures
      `c₁ = ψ̇`, `c₂ = sin ψ / r` and no small-slope expansion anywhere. The
      Euler–Lagrange equations reduce to a first-order system in
      `(r, z, ψ, M, η)` solved as a BVP, with slope continuation as a safety
      net because a BVP that fails from a bad guess returns a *plausible wrong
      shape* rather than an error.
- [x] **Free accuracy diagnostic.** The Lagrangian has no explicit `s`
      dependence, so its Hamiltonian is conserved and equals the axial force —
      zero for an inclusion nobody pulls on. Imposed as a boundary condition,
      its drift along the solution then measures integration error directly:
      **|H|max ≈ 7e-11**.
- [x] **Small-slope validation, the roadmap's own criterion: they agree.**
      Relative discrepancy over slope² converges to a constant **0.746**, so
      the error is exactly the O(|∇h|²) the Monge expansion discards rather
      than merely something that shrinks. At slope 0.02 the energies differ by
      0.02%; at 0.05, by 0.2%.
- [x] **The measured result, and it disagrees with our own Round 3 number.**
      At the 7WLT geometry (inclusion radius 8.69 nm, contact slope 1.99 = 63°):
      footprint energy **92.2 → 25.3 k_BT** and excess area **622 → 179 nm²**.
      The linearised values are **3.65× and 3.48× too large**.
- [x] **A conclusion reverses.** Round 3 compared the dome's *exact* 256 nm²
      against the footprint's *linearised* 622 nm² and reported the footprint
      holding "2.4× as much deformable area as the dome". Measured
      consistently the footprint holds **0.70× the dome** — less, not more.
      `docs/SCIENCE.md` now states this and separates it from Haselwandter &
      MacKinnon's actual claim, which is about tension *sensitivity* (area
      released between states) and which absolute stored area never tested.
- [x] Robustness: invariant to domain truncation 8λ→40λ (6 s.f.), grid,
      tolerance and continuation path; correction factor 3.46–3.67× across
      κ = 20–25 k_BT and γ = 0.42–3.0 mN/m, so not a parameter artefact.
- [x] `DomeModel` gained `contact_slope()`, `footprint_nonlinear()`,
      `footprint_area_nonlinear()` and `compare_footprint_theories()`;
      `footprint_area()` now documents that it is not quantitative here.
- [x] 18 tests (`tests/test_elastica.py`). The load-bearing ones do **not**
      reuse the hand-derived shape equations: agreement with linear theory,
      axial-force conservation, re-evaluation of the exact functional in the
      Monge gauge, and a perturbation test that no nearby admissible shape has
      lower energy. Suite 254 → **272 passed**; GUI smoke test clean.

### Round 19 — Pore hydration and wetting ✅
- [x] `piezo1/analysis/hydration.py` implements the **Rao et al. 2019** (PNAS
      116:13989, PMID 31235590) hydrophobic-gating heuristic: kernel-smoothed
      pore hydrophobicity on the normalised Wimley–White scale, joined to the
      radius profile, looked up against their MD-derived water free-energy
      landscape; residues above **1 RT = 2.6 kJ/mol** are flagged and the score
      is the sum of shortest distances to that contour, **Σd > 0.55 ⟹ closed**.
- [x] **We use the published landscape, not a redrawing of it.** The 100×100
      grid ships in the CHAP repository under the **MIT licence**, so it is
      downloaded (`python -m piezo1.io.fetch`) and used directly. Reading a
      boundary off a figure would have been exactly the silent correctness bug
      Round 17 was about. Reported AUROC for this heuristic is **0.91 against
      0.59 for minimum radius alone** — the reason the round exists.
- [x] **Independent check that we read the grid correctly:** our extracted
      1 RT contour recovers the paper's stated critical radii — **0.10 nm** at
      the hydrophilic end rising to **0.43 nm** at the hydrophobic end, against
      their "hydrophilic pores wet below 0.2 nm, hydrophobic ones can hold a
      barrier out to ~0.4 nm".
- [x] **The result asked for: 8YEZ Σd = 0.82 → non-conductive; 11ZC Σd = 0.00
      → conductive.** Both correct.
- [x] **And for the right reason, demonstrated by control rather than
      asserted.** Holding every radius fixed and replacing the hydrophobicity
      scale with a uniform hydrophilic value collapses 8YEZ from 0.82 to
      **0.00 (conductive)**. The verdict is chemistry, not a radius threshold
      in disguise. Concretely: 8YEZ's F2451/V2454 sit at **0.325 nm** and are
      called dewetted, while 11ZC's bottleneck at **0.330 nm** is called wet —
      the same radius, opposite verdict.
- [x] **The heuristic rediscovers the curated gate.** Seeing only coordinates
      and a hydrophobicity scale, it flags F2451 and V2454 (curated hydrophobic
      gate / pore-lining) and R2467, F2468 (curated cytoplasmic constrictions).
- [x] **A limitation found by testing beyond the two structures asked for, and
      reported rather than hidden.** 7WLU and 8IXO have 0.098 nm bottlenecks —
      far too narrow for water — but hydrophilic linings, so their Σd is small
      and the Rao score *alone* calls them open. The heuristic answers "would
      water dewet here?", not "does water fit here?". `WettingPrediction`
      therefore exposes `hydrophobic_gate` and `sterically_occluded`
      separately, and `conductive` requires neither. With both, all five states
      come out right: 8YEZ and 7WLT non-conductive on both counts, 7WLU and
      8IXO non-conductive on sterics, 11ZC conductive.
- [x] Reachable from the CLI as `piezo1 hydration <PDB>` and through the
      `ANALYSES` registry; 16 tests (`tests/test_hydration.py`), which skip
      rather than fail when the grid is not downloaded. Suite 272 →
      **288 passed**; GUI smoke test clean.

### Round 20 — Statistical rigour for the blind test ✅
- [x] Permutation test, effect size with a bootstrap CI and a pre-registered
      decision rule were delivered ahead of schedule in Rounds 6–7
      (`analysis/validation.py`, `docs/PREREGISTRATION.md`). This round adds
      what was genuinely missing: **power, multiplicity and cross-validation**,
      in a new `piezo1/analysis/design.py`. The split is deliberate —
      `validation.py` answers "did it work?", `design.py` answers "could it
      have worked, and did we look too many times?".
- [x] **Leave-one-out cross-validation**, with every label-consuming step
      inside the fold. On the Round 7 predictors: AUROC **0.535 out-of-sample
      against 0.542 in-sample**, optimism +0.007 — so there was no hidden
      overfitting inflating the original number.
- [x] **Benjamini–Hochberg FDR** with a named primary endpoint. BH rather than
      Bonferroni because the sequence predictors read the same evolutionary
      signal and are strongly correlated. Worked illustration: of six candidate
      predictors, **three clear p < 0.05 and none survives correction**.
- [x] **Power analysis, and it produced the round's real finding.** Simulating
      the pre-registered test at Round 7's actual group sizes (16 vs 9), under
      both a normal model and resampling from the observed heavy-tailed ΔΔG
      pool — which agree — **80% power is reached only at |δ| ≥ 0.55**, past
      'large'. Power at the observed effect was **0.13**; at a medium effect,
      0.35; at a large effect, 0.60.
- [x] **This qualifies the recorded null without revising it.** Round 7
      excludes a *large* mechanical effect and is close to uninformative about
      a small or medium one. Added as `docs/VALIDATION.md` §6b, explicitly
      marked as not amending §§1–3. The §6 diagnostic (99.8% of ΔΔG variance is
      between-position) remains the mechanistic explanation and is independent
      of this; both are true, and the power limit was the one not stated at the
      time.
- [x] **Sample sizes needed** at 80% power, equal groups: **42** variants for a
      large effect, **98** for medium, ≥600 for small. Only 25 of the 68
      curated variants survive the inclusion criteria and relaxing them cannot
      reach ~45 — so **a confirmatory test of anything below a large effect is
      not available from this variant set.** Round 22 is bound by that.
- [x] **`docs/NEGATIVE_RESULT_PROTOCOL.md` written first**, before the Round 22
      hypothesis it governs, so the rule cannot be tuned to a result already
      seen. Fixes what must exist before a test is run, that a recorded result
      is never revised (only superseded by a new entry that points back, as
      Round 18 did to Round 3), and that re-running a hypothesis with a changed
      predictor is a **new** test needing a **new** pre-registration.
- [x] A sign error caught by its own diagnostic: `power_curve` injected the
      effect into the wrong group, reporting the power to detect the *opposite*
      direction. The achieved-vs-target Cliff's delta check exists precisely
      because a power curve looks perfectly sensible when inverted.
- [x] 23 tests (`tests/test_design.py`), including that the fast subset-sum
      permutation path agrees with the real test, and that the rejection rate
      under a true null is α rather than more. Suite 288 → **311 passed**; GUI
      smoke test clean.

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

### Round 21 — Let the GUI reach the engine ✅
- [x] New **Analysis dock** (`ui/panels/analysis_panel.py`) with three tabs —
      Pore, Pockets, Residue maps — driven by `ui/analysis_controller.py`.
      Every analysis runs on a worker thread: a pocket search or PRS scan takes
      seconds, and a window that stops repainting during a calculation is
      indistinguishable from one that has crashed.
- [x] **Pore profile with the hydrophobicity trace alongside it**, on
      independent left/right axes, with the bottleneck and the 1.5 Å water
      radius marked. Clicking the plot selects the residues lining the pore at
      that height. Written as a QPainter widget (`ui/profile_plot.py`) rather
      than adding matplotlib or pyqtgraph — one plot type, needs to repaint at
      interactive rates inside a dock, and must match the dark theme.
- [x] **Verified the GUI is not a second implementation.** A test asserts the
      pore worker reproduces `pore_profile` and `predict_wetting` exactly, and
      the smoke test now reads the panel's own label: bottleneck **0.95 Å**,
      Rao score **0.82 → non-conductive (sterically occluded + hydrophobic
      gate)** — identical to `piezo1 hydration 8YEZ`.
- [x] Pockets listed with volume, buriedness and lining, selecting one
      highlights and frames it. Verified top pocket 6594 Å³, buriedness 1.00.
- [x] **Colour-by-conservation and colour-by-PRS**, through the existing
      `ColorBy.VALUE` path rather than a new one. They appear in the Analysis
      dock's scalar menu once computed rather than in the main colour dropdown,
      because unlike domain or chain they do not exist until an analysis has
      run — an inert entry in the main list would be a dead control.
      Conservation returns 2484 residues (mean 0.770, 901 above 0.95) after
      dropping positions with ortholog coverage below 0.7, since those measure
      the alignment rather than selection pressure. PRS folds the three
      protomers, giving 1279 residues, coupling 0.473–20.2.
- [x] **Session save/load and report export** on the File menu
      (`ui/session_controller.py`, kept out of `main_window.py` which was
      already near the 500-line limit). Sessions record *what was being looked
      at* and never results — a test asserts no eigenvalues or coordinates
      reach the file. Round-tripped in the smoke test, not just unit-tested.
- [x] Two traps handled: restoring four appearance settings blocks signals and
      re-emits once, so a 120k-atom trimer rebuilds once rather than four times
      and does not flicker; and residues with no computed value take the map's
      **floor rather than zero**, since injecting a zero into a conservation
      map running 0.6–1.0 would rescale the whole legend around positions that
      were never measured.
- [x] 19 tests (`tests/test_ui_analysis.py`), running the real widgets on the
      **offscreen** Qt platform rather than mocking them — mocks would not have
      caught either of the two silent GUI breakages this project has had.
      `scripts/screenshot_app.py --analysis` exercises the new panel end to
      end. Suite 311 → **330 passed**.

### Round 22 — A second, stated hypothesis for the variant test ✅
- [x] `docs/PREREGISTRATION_ROUND22.md` written **and committed in its own
      commit before the test was run**. Round 7's null stands unrevised.
- [x] **The design finding, obtained blind and worse than Round 7.** Of 39
      variants with a directional label only 26 are single-residue
      substitutions, and **11 of the 13 dropped are loss-of-function** nonsense,
      frameshift or deletion variants. That is biology, not curation: loss of
      function is commonly achieved by truncation and no missense predictor can
      score a stop codon. Usable design **20 GoF vs 6 LoF**, so 80% power needs
      **|δ| ≥ 0.61** against Round 7's 0.55; power at a large effect is 0.52.
- [x] **Declared EXPLORATORY** under `NEGATIVE_RESULT_PROTOCOL.md` §3. The
      confirmatory alternative was declined because at this n a "confirmatory
      null" would exclude only effects beyond large.
- [x] Primary endpoint named in advance: **FoldX ΔΔG**, hypothesis *LoF more
      destabilising than GoF* — loss can be achieved by breaking the protein,
      gain cannot, since a channel that opens too readily must still fold,
      traffic and gate.
- [x] **Result: nothing separates the classes, and the primary runs opposite to
      the hypothesis.** Cliff's δ = **−0.211**, CI **[−0.684, +0.298]**,
      AUROC **0.395**; mean ΔΔG **0.767 (LoF) vs 1.309 (GoF)**. Three
      statistics agree on direction, so it is not a sign error. The interval
      spans zero, so this is uninformative rather than a reversal.
- [x] **The pre-registered objection explains it, and was recorded first.**
      Excluding the truncating LoF variants removes exactly the "break the
      protein" mechanism, leaving a LoF subset selected for *not* being
      truncating — a destabilisation predictor asked about the one subset of
      loss-of-function variants that does not act by destabilisation.
- [x] Secondary family (BH): **nothing significant, nothing close** — smallest
      q = 0.448. AlphaMissense +0.183, EVE +0.133, ESM-1b −0.350, conservation
      +0.200, mechanical +0.125, every interval spanning zero. The
      pre-registered expectation of *no separation* from the three pathogenicity
      predictors **held**.
- [x] Combined score, equal weights, leave-one-out **AUROC 0.535**, optimism
      **0.000**. Combining recovers nothing the parts lacked.
- [x] `docs/VALIDATION_ROUND22.md`; 10 tests pinning the recorded numbers.
      Suite 342 → **352 passed**.
- **Standing position:** two nulls from two different predictors, both
  pre-registered, neither revised. The central claim remains **untested at
  adequate power**, and the binding constraint is data (Round 27), not method.

### Round 23 — Packaging and one-command reproduction ✅
- [x] `pyproject.toml` with metadata, classifiers, console entry points
      (`piezo1`, `piezo1-gui`) and tool config. **GUI dependencies are
      optional**: everything below `render` runs headless, which is what makes
      the science testable without a display, so PyQt is an extra rather than a
      requirement. Curated resources ship as package data; `ref/` and `data/`
      downloads never do. Verified with an editable install — the `piezo1`
      console script resolves and runs.
- [x] Pinned locks: `environment.lock.yml` (230 conda entries) and
      `requirements.lock.txt` (25 pip entries, all `==` or VCS-pinned).
- [x] `Makefile` with `env`, `lock`, `fetch`, `resources`, `test`, `lint`,
      `sizes`, `gui`, `figures`, `validate`, `verify`, `quick`, `reproduce`.
      Each wraps `conda run`, so `make test` works from a bare shell.
      `make sizes` enforces the project's 500-line limit as a build step
      rather than a habit.
- [x] `scripts/reproduce.py` runs every step in dependency order — fetch,
      resources, tests, both pre-registered validations, figures, GUI
      screenshots — and finishes by verifying the documentation.
- [x] **The part that earns the round: `piezo1/analysis/claims.py`.** The
      project asserts a lot of specific numbers in prose, and prose does not
      fail a test suite; a solver rewrite can leave `docs/SCIENCE.md`
      confidently stating a value the code stopped producing. Seventeen claims
      now name the documented number, its tolerance, the document it appears
      in, and a callable that recomputes it from scratch.
- [x] **All 17 reproduce**, in about 10 s: dome 9.7245 nm, gating overlap
      0.7048, T₅₀ 2.7109 mN/m, γ 0.4200 mN/m, λ 13.9897 nm, elastica 25.2702
      k_BT / 178.93 nm² / 3.649× overestimate, pore 0.9518 and 3.2973 Å,
      wetting 0.8227 and 0.0000, CDS identity 1.0000, and the four frozen
      validation results.
- [x] **Frozen claims cannot be fixed by editing prose.** The four recorded
      validation numbers are marked, and drift in one prints an explicit
      instruction not to edit the document to match but to work out why the
      computation changed.
- [x] The drift detector is itself tested — a deliberately wrong claim must be
      reported, and a claim that cannot run for want of downloaded data must be
      reported as *skipped* rather than as drift, so a fresh clone does not
      look broken. 14 tests; suite 376 → **390 passed**.

### Round 24 — Performance ✅
- [x] **Profiled first, and the roadmap's premises were partly wrong.** PRS is
      0.52 s and was never a problem; pocket detection was 4.2 s, not ~10 s.
      The genuine costs were **SASA at 7.5 s**, which the roadmap did not
      mention at all, and the ensemble at 12.4 s — of which **99% was mmCIF
      parsing**, not the PCA.
- [x] **SASA 7.54 → 1.27 s (5.9×), bit-identical.** The hot loop built a
      (256 × neighbours × 3) array per atom and took a square root it did not
      need. Expanding `|t−x|² = |v|² + r² + 2r(p·v)` turns that into one BLAS
      product, and `d ≥ r` and `d² ≥ r²` decide the same way for non-negative
      values. `np.array_equal` on all 31,599 atoms.
- [x] **Ensemble 12.35 → 2.05 s (6.0×).** The mmCIF tokenizer walked characters
      in Python; 99.5% of lines contain no quote or comment, where
      `str.split()` is exactly equivalent and runs in C. Verified on **245,528
      lines of deposited structure with zero mismatches** before relying on it.
      The careful path — where the whitespace bug once shifted every column by
      one — is untouched and merely bypassed. Also replaced an O(n) `pop(0)`
      running several million times per load.
- [x] **Conservation 3.67 s → 0.003 s** via a disk cache keyed on a **content
      hash** of the reference and ortholog sequences, so re-fetching or
      changing the reference invalidates it automatically. A cache that can go
      stale is worse than none — it would report last week's conservation
      against this week's alignment.
- [x] Pockets 4.21 → 3.60 s (squared distances, and points already inside are
      not carried into the next block — `inside` is monotone, so this cannot
      change the answer). Structure load 0.37 → 0.19 s. Feature table
      7.99 → 4.07 s.
- [x] **Overall 33.8 → 12.5 s (2.7×) with every number unchanged**: SASA total
      197490.5582 Å², top pocket 6593.6 Å³, PC1 0.9000, bottleneck 0.9518 Å.
      The test suite itself dropped 118 → 97 s.
- [x] 11 tests (`tests/test_performance.py`) asserting **identity, not
      closeness** — the fast tokenizer against the careful one over a whole
      deposited file, SASA against the direct formulation, the Monte-Carlo
      volume against the un-optimised loop, and the conservation cache against
      an uncached rebuild. Timing assertions are loose ceilings only; pinning a
      runtime would fail on a slower machine for no scientific reason.
      Suite 408 → **419 passed**; all 17 documented numbers still reproduce.

### Round 25 — The teaching layer ✅
- [x] An **11-step guided tour** (`piezo1/tour.py`, Qt-free; `ui/tour_panel.py`
      and `ui/tour_controller.py` for the GUI half) walking the mechanism:
      trimer → blades → dome → footprint → lever → gate → open state → normal
      modes → gating energetics → a variant → what the project cannot do.
      Reachable from Help → Guided tour (F2) or its own dock.
- [x] **Every number a step states is computed when the step runs.** None is
      written into the prose. A tour that narrated "the dome radius is 9.7 nm"
      would be a fourth place for that number to live and go stale, beside the
      code, the documentation and the claims registry. Published comparisons
      come from the **parameter registry**, so a step inherits the provenance
      rule rather than side-stepping it — there is a test that changing
      `dome.published_radius_closed` changes the tour text.
- [x] The controller calls the **same** controllers the panels use rather than
      computing anything itself. A teaching tool that quietly disagreed with
      the application it is teaching would be worse than none.
- [x] **The tour ends on the two null results**, with a test asserting the
      p-value, AUROC and effect size are present. A learning instrument that
      only shows its successes teaches the wrong lesson.
- [x] 12 tests (`tests/test_tour.py`), including that every step degrades
      gracefully before anything is computed and none raises on junk input —
      a tour must never crash the application. The GUI smoke test now walks all
      11 steps. Suite 419 → **431 passed**.

---

## Block G — after the Round 20 review  *(added 2026-08-06)*

**Where twenty rounds leave us.** The physics chain is closed and every link
validated against an independent published number. The engine now also knows
its own limits, which took three rounds to establish and is the more useful
half of the progress:

- Round 18 overturned a Round 3 headline — the linearised footprint is 3.5×
  too large at PIEZO1's 63° contact slope, and the footprint holds *less*
  excess area than the dome, not 2.4× more.
- Round 19 showed the pore verdict must separate steric occlusion from
  hydrophobic gating, because two deposited states are shut for one reason and
  not the other.
- Round 20 showed the Round 7 null excludes only a *large* effect, and that no
  confirmatory test below that is available from 25 variants.

The through-line: every round that looked hard at a previously reported number
found the number was reported with more confidence than it had earned. That is
the failure mode to keep hunting.

### Round 26 — Substitution-aware mechanics ✅
- [x] **The cause was algebraic, not statistical.** The old model scaled every
      contact of a residue by *one number*, so ΔΔG = (s−1)·Q(position) — a
      rank-one product in which the substitution enters only as a multiplicative
      scalar. Four substitutions at R2456 could then differ solely by a factor
      and had to rank every position identically. No refinement of `s` could
      have fixed that; the separability itself had to go.
- [x] `piezo1/analysis/substitution.py` scales **each contact individually**, by
      properties of the new residue *and the partner it touches*: packing,
      charge (felt only at charged partners), hydrogen-bond complementarity,
      proline stiffening restricted to sequence-local contacts, and glycine
      softening. Different substitutions now perturb different *subsets* of
      contacts. Measured directly: the per-contact patterns for R2456H/K/P/C
      correlate 0.62–0.98 where they were 1.00 by construction.
- [x] **The pre-registered criterion is met.** On the multiply-substituted
      positions — six of them, including all four R2456 substitutions, as the
      criterion specified — within-position variance goes from **4.9% to
      52.5%** against a threshold of 20%.
- [x] **Reported honestly alongside it:** across *all* 35 substituted positions
      the figure is 2.4%, up from 0.8%. That is not a contradiction — 29 of
      those positions carry a single substitution and contribute exactly zero
      within-variance by construction, so including them drives the statistic
      down for reasons that have nothing to do with the model.
- [x] **Round 7's frozen result is untouched.** Its script passes no sequence,
      so the model keeps the uniform scale; there is a test asserting the two
      agree to 1e-12, and all 18 documented numbers still reproduce.
- [x] **No phenotype comparison was run.** This is method development against a
      variance criterion, not a hypothesis test. Whether the new distinctions
      point in the *right* direction is untested and needs a new
      pre-registration under `docs/NEGATIVE_RESULT_PROTOCOL.md` §7.
- [x] Eight new registered parameters for the weights, each with bounds and a
      stated basis. 17 tests (`tests/test_substitution.py`), including that a
      spring may weaken but never invert — a negative spring makes the Hessian
      indefinite and the quadratic form stops being an energy. Suite 445 →
      **462 passed**.

### Round 27 — Expand the phenotyped variant set ✅
- [x] **ClinVar gives pathogenicity, not direction**, and this project needs
      direction. What makes direction recoverable for PIEZO1 is that its two
      diseases have opposite mechanisms: dehydrated hereditary stomatocytosis
      is dominant gain-of-function, generalised lymphatic dysplasia is
      recessive loss-of-function. So a condition can imply a direction — and
      that is **weaker evidence than measuring the current**, recorded per
      variant rather than pooled.
- [x] 354 pathogenic/likely-pathogenic records fetched; **232 pass the
      wild-type gate against Q92508**, with 3 rejected for disagreeing (P481,
      I505, C463) and 117 for an unparseable protein change.
- [x] **The inter-curator ambiguity is real and substantial: 11 of 63 directed
      records are reported under *both* diseases.** ClinVar submitters routinely
      attach the whole gene's disease list to a variant. Those carry no
      direction and are excluded, not resolved by preferring one.
- [x] **An independent check on the inference.** Nine of the ClinVar variants
      are already curated from electrophysiology; the condition-based direction
      agrees with the measured one **8 times out of 9**. The single
      disagreement — V598M, curated GoF, inferred LoF — is *reported, not
      resolved*: our own record reads "increased opening (one report); no
      change in another", so the literature is genuinely mixed.
- [x] **Achieved n: the directional missense set goes from 26 (20 GoF, 6 LoF)
      to 46 (27 GoF, 19 LoF).** The loss-of-function class — the thing that
      made Round 22 underpowered — more than triples, 6 → 19.
- [x] **Minimum detectable effect: 0.61 (Round 22) → 0.41.** Now inside
      "large" rather than beyond it, and power at a large effect rises from
      0.50 to 0.83.
- [x] **But a medium effect is still out of reach, and the roadmap asked to say
      so: power at δ = 0.28 is 0.49, not 0.80.** Reaching it at this 27:19
      ratio would need **104** variants. The constraint has loosened, not
      lifted.
- [x] `piezo1/analysis/variant_sets.py` assembles a set at a **stated evidence
      level**, defaulting to the conservative `measured` one so a caller who
      does not think about evidence strength gets the smaller answer rather
      than the larger. The original `variants.json` is **untouched** — Round 7
      and Round 22 reference it, and growing it underneath a frozen result
      would invalidate it with nothing appearing to change.
- [x] 15 tests (`tests/test_variant_sets.py`), including that the wild-type
      gate still rejects and that the evidence levels cannot be pooled by
      accident. Suite 462 → **477 passed**.

### Round 28 — Nonlinear footprint in the gating energetics ✅
- [x] **ΔA is a change, not an absolute area** — a distinction Round 3's
      framing invited getting wrong. `DomeModel.gating_area_change` measures
      the closed→open difference: the dome's projected area grows, and the
      footprint releases the excess area it was storing.
- [x] Both endpoints measured from deposited coordinates: **7WLT** R_c 9.72 nm,
      contact slope 1.992 (63.3°); **7WLU** R_c 18.38 nm, slope 0.839 (40.0°).
- [x] **The correction is larger on the difference than on either endpoint.**
      The closed state sits at 63°, where the linear theory is badly wrong; the
      open state at 40°, where it is much less so. So the footprint's
      contribution to ΔA falls from **463 nm² to 71 nm² — a factor of 6.5**,
      against the 3.5× Round 18 found for the closed state alone. Total ΔA
      664 → 272 nm².
- [x] **The roadmap's question, answered explicitly: the linear version did
      *not* agree better.** T₅₀ moves from 0.060 to 0.147 mN/m — toward the
      measured 2.7 ± 0.1, not away. There is a test asserting that direction,
      because a wrong model fitting a right number was a real possibility.
- [x] **But the correction does not close the gap, and that is the more useful
      half of the result.** Improving the membrane physics 6.5× moves T₅₀ by
      2.4× and leaves it **~18× below measurement**. The structural ΔA is still
      34× the functional 8 nm². So the structural-versus-functional
      discrepancy is **not a membrane-modelling error** — it is about which
      quantity each number measures, and no refinement of the footprint will
      fix it.
- [x] `compare_gating_area_routes` reports all four routes side by side with
      the T₅₀ each implies, and a test requires every row to satisfy
      T₅₀ = ΔG₀/ΔA so the table cannot drift from the model. The functional
      route still reproduces Cox at 4.99 against 5.1 ± 0.2.
- [x] 12 tests (`tests/test_gating_area.py`); the 71 nm² is pinned as a
      documented claim. Suite 477 → **489 passed**, 19 numbers reproduce.

### Round 29 — Uncertainty on every reported quantity ✅
- [x] `piezo1/analysis/uncertainty.py` with **three kinds of spread kept
      apart**, because conflating them would be a second kind of
      overconfidence rather than a cure for the first: a **bootstrap**
      confidence interval (resampling data), a **sensitivity range** (varying a
      method choice — a network cutoff has no sampling distribution), and a
      **parameter range** (propagating a registered input over its published
      values). The class names and the printed summaries all say which.
- [x] **The dome radius reframes a headline claim.** 9.73 nm with a 95% CI of
      **[8.83, 10.34]** over 66 surface points — and the published 10.2 nm sits
      *inside* it. The project has been reporting this as "close but not
      exact"; the honest statement is that the two are **statistically
      indistinguishable**, which is a stronger claim of consistency and a
      weaker claim of precision.
- [x] **The gating overlap is not robust to three digits.** Across network
      cutoffs from 10 to 20 Å it ranges **0.554–0.723**, a 24% spread, and
      non-monotonically. The qualitative result — a substantial overlap carried
      entirely by A-symmetric modes — survives every cutoff; the specific 0.705
      does not, and was being quoted as though it did.
- [x] **Ensemble PC1 is 0.900 [0.796, 0.972]** over ten structures. Leaving out
      11ZC alone moves it to 0.832, the largest single influence, but no entry
      dominates.
- [x] Footprint energy over the published κ range (20–25 k_BT): 25.27–26.94
      k_BT, a 6.6% propagated spread.
- [x] The confidence level is **derived from `stats.alpha`** rather than
      repeated, so the two cannot drift; setting α to 0.01 widens every
      interval. The parameter sweep restores the registry even when the
      statistic raises, since leaving it modified would make every later number
      in the session incomparable with the documentation.
- [x] **Stated on every result: none of this captures model error.**
      Bootstrapping a sphere fit says how well the sphere is determined, not
      whether a sphere was the right question.
- [x] `scripts/report_uncertainty.py`; 18 tests (`tests/test_uncertainty.py`),
      including that a statistic failing on most resamples raises rather than
      silently narrowing the interval. Suite 489 → **507 passed**.

### Round 30 — Adversarial review of the whole chain ✅
- [x] `piezo1/analysis/crosscheck.py` re-derives each headline by a route
      sharing no machinery with the original, plus
      `scripts/crosscheck_chain.py` to run them all.
- [x] **Dome curvature, two alternatives, and the first one disagreed.** A
      parabola through the radial height profile gave **8.12 nm against the
      sphere fit's 9.72 — 16.6% apart.** Diagnosed on synthetic caps of *known*
      radius: the parabola is a **shallow-cap approximation**, accurate to 0.6%
      at an 8.6° contact angle and **25.8% low at 63.4°**, which is exactly
      where PIEZO1 sits. This is the Round 18 small-slope lesson reappearing in
      a geometry method. The sphere fit is exact at every slope on the same
      synthetics.
- [x] So a genuinely valid alternative was built — the **exact cap relation
      R = −(h²+r²)/2h** inverted per point, no fitting, no expansion. It gives
      **10.17 nm against 9.72, 4.5% apart**, both inside Round 29's bootstrap
      interval [8.83, 10.34] and both consistent with the published 10.2.
- [x] **Mode overlap without superposing anything.** Distances are invariant to
      rotation and translation, so comparing the transition and each mode in
      *pairwise distance changes* uses no Kabsch fit and no protomer matching.
      **0.641 against 0.705 — 9.0%**, so the overlap is not an artefact of the
      superposition.
- [x] **T₅₀: a disagreement that turned out to be my error, and worth
      recording.** The first alternative used the analytic steady state and
      returned 0. Investigation showed why: at equilibrium the channel sits
      **~96% inactivated at every tension**, so steady-state open occupancy
      runs only 0.030–0.036 and has no half-maximum. It was computing a
      *different quantity*. T₅₀ is necessarily a property of the peak
      transient — which is also what a patch-clamp measures. Replaced with
      adaptive Runge–Kutta integration of the same master equation: **2.727
      against 2.711, 0.6%**.
- [x] The parameter audit caught the new module twice, and the second catch
      prompted deriving the open-state index from `STATE_NAMES` rather than
      writing `1` — a hardcoded index would silently read the wrong occupancy
      if the state order ever changed.
- [x] 11 tests (`tests/test_crosscheck.py`), including that the sphere fit is
      exact on synthetic caps at every slope, that the parabola's error grows
      monotonically with slope and always underestimates, that the
      distance-space overlap is invariant to a rigid motion, and that the
      steady-state route must **not** reproduce the peak-based T₅₀. Suite
      507 → **518 passed**.

---

## Block H — HaloTag constructs and calcium  *(added 2026-08-06, on request)*

Plan and feasibility review in **`docs/HALOTAG_CALCIUM_PLAN.md`**. The companion
project `halotag_binding_sim` supplies a complete, well-provenanced kinetic
model of covalent HaloTag labelling that is **entirely non-spatial** — it knows
a PIEZO1 trimer has three sites, not where they are. This project supplies the
geometry, the structure and the ion physics. Import that kinetics; do not
reimplement it.

Feasibility already checked: HaloTag structures **6U32** (1.8 Å, with ligand
bound) and **6U2M** (HaloCaMP calcium indicator) exist; PIEZO1's C-terminus sits
at **z = −60.4 Å, 8 Å below the CTD constriction**, putting a tag centre 4–6 nm
from the pore exit; and the nanodomain there is of order **200 µM** against
BAPTA's **~0.2 µM** Kd.

### Round 31 — HaloTag fusion geometry
- [x] Fetched 6U32 (1.8 Å, TMR-HaloTag ligand bound). `structure/fusion.py`
      builds the fusion at all three C-termini as an **accessible-volume
      envelope**, not a pose. Measured tag inputs: Rg **17.6 Å**, N-terminus
      **19.9 Å** from the centre, ligand 21.8 Å from that N-terminus. On 8YEZ
      the envelope holds 30,698 positions / **246 nm³**, with **65%** of the
      tether's reach blocked by the channel. `FusionModel.seams()` returns the
      anchor→tag segments for a renderer to draw differently; reachable as
      `python -m piezo1.cli fusion 8YEZ`. Drawing it in the GUI is still open —
      `hybrid.py`, the stated model for the seam, is itself unimplemented.
- [x] *Validate:* **two of three pass.**
      **C3 symmetry — PASS.** Deviation `0.0000 Å` on all 20 entries (exact by
      construction: one envelope is solved and rotated).
      **No steric clash — PASS on 18/20.** Clearance 21.5 Å against the tag's
      17.6 Å radius on 8YEZ. 3JAC (17.6 Å) and 11ZC (15.7 Å) are marginal.
      **Tag centre 4–6 nm from the pore exit — MISSES.** Measured **3.95 nm**
      on 8YEZ and **3.27–4.21 nm (mean 3.81)** across all 20 structures. Not an
      artefact of the unverified linker: sweeping it 1→30 residues moves the
      answer only 3.0–4.0 nm, and *downward*, because a longer tether wraps
      further round the channel. The 4–6 nm estimate came from adding the tag's
      ~2 nm anchor-to-centre offset to the anchor's 2.6 nm from the pore exit,
      which assumes the tag points straight away from the channel; averaged over
      accessible directions the mean is pulled in. The band is not unreachable —
      the envelope spans 1.7–7.9 nm and **51%** of it lies inside 4–6 nm — so
      the window describes an achievable position, not the ensemble mean.
      Round 35's nanodomain estimate should use the envelope, not the centroid.
- [x] Two sign faults found and fixed, each of which gave a confident wrong
      answer: the pore exit taken as the lowest protein atom is a distal blade
      tip, and `SymmetryAxis.direction` has no fixed sign (it returns −z for
      7WLT and 8YFG). Together they put those structures' tags 15–16 nm from
      the "pore exit" against 3.9 nm for the same construct on 8YEZ.

### Round 32 — Labelling on the structure
- [x] `analysis/labelling.py` imports the kinetics from `halotag_binding_sim`:
      exposure `E(t)`, per-site `p(t) = a(1 − e^{−k_on·E})`, and Binomial(3, p)
      over the trimer. `label_sites()` drives per-site occupancy on the **real**
      tag positions from Round 31's fusion model. The 1:2:3-dye histogram is in
      `docs/img/labelling.png` and `python -m piezo1.cli labelling 8YEZ`.
      Three references added through the title-verification gate
      (`los2008halotag`, `grimm2015jf`, `bertaccini2025piezo1`); 8 parameters
      registered. **The brightness *animation* is not done** — the histogram and
      per-site occupancy are, but nothing is rendered on the trimer over time.
- [x] *Validate:* **exact, to the last bit.** `compare_with_source()` re-runs
      the original functions over 241 time points: max |Δ| = **0.0** for
      `p_site`, for `p³`, for the occupancy distribution and for the sampled
      histogram, and the Monte-Carlo dye counts are **identical channel for
      channel** (`dye_counts_identical: true`). The sampler reproduces the
      source's two `rng.random((n,3))` draws in the same order, because a
      different order gives a statistically identical population and a
      numerically different one — which would have hidden a divergence behind
      sampling noise.
- [x] **What the model says, which is not what the round expected.** At the
      standard protocol (200 nM JF646, 30 min, live cell) labelling is complete
      in **54 s** to 99%: per-site p = 1.0000 and the mixture is **100%
      three-dye**. So at any realistic concentration there is *no* kinetic dye
      mixture — a mixture needs sub-nanomolar ligand or ≤1 min incubation. A
      population of chemically unreactive tags does produce one at every time,
      because the ceiling is `active_fraction³`: at 90% reactive the mixture is
      72.9% three-dye / 24.3% two-dye and no incubation removes it. **Two
      different things get called "sub-saturation labelling", and under a
      saturating protocol only the unreactive-tag route is available** — so an
      observed 1:2:3 brightness mixture argues for unreactive tags, not for a
      short incubation. `labelling.k_perm_live` and `labelling.active_fraction`
      are both registered `unverified`; they are the assumptions this rests on.

### Round 33 — Calcium permeation
- [x] `physics/permeation.py`: steady-state drift-diffusion per species over the
      measured `PoreProfile`, Scharfetter–Gummel discretised, with Hall access
      resistance in series at both mouths. Gated by the Round 19 wetting
      verdict. Reachable as `python -m piezo1.cli permeation 11ZC`.
- [x] **The Poisson half does not converge, and the reason is physical.**
      Direct Gummel diverged (−0.37 V → −171 V → −2×10¹⁶ V); a proper Newton
      step with the screening derivative still would not converge. In 150 mM the
      **Debye length is 5.7–8.1 Å against an open bottleneck radius of 3.3 Å**,
      so the double layers overlap completely and the pore has no electroneutral
      core to relax onto. The potential is therefore solved in the
      **electroneutral limit** (current continuity), which converges and agrees
      with the independent closed-form `series_conductance` to **1.5%**
      (41.0 vs 40.4 pS). `debye_length()` is reported on every result.
- [x] *Validate:* **41.0 pS** for the open 11ZC against a published **25–30 pS**
      — high by about half. Closed structures give exactly 0.
      **Which mechanism blocks which**, as asked: 8YEZ is shut by **2**
      mechanisms (sterically, 0.95 Å bottleneck; *and* hydrophobic gate, wetting
      0.82), 7WLU by **1** (sterically only, 0.98 Å, wetting 0.11).
      `blocking_mechanisms()` returns every reason rather than the first,
      because returning early collapsed exactly this distinction.
- [x] **The agreement cannot be claimed as a prediction.** Sweeping the two
      *unmeasured* confinement parameters over plausible ranges — in-pore
      diffusivity 0.25–1.0 of bulk, ion radius 1.0–2.0 Å — moves the answer
      across **16–94 pS**, a 5.8× span that straddles the measurement. The model
      can be made to agree, but only by choosing values nobody has measured, so
      that would be tuning. Both are registered `unverified`.
      Calcium at 2 mM carries **<5%** of the current, consistent with PIEZO1's
      weak selectivity.
- [ ] Particle animation whose flux is **set by** the computed current, with the
      HUD stating what the frame rate is in real time. The morph clock's
      discipline applies. **Not done** — the physics is in place and the current
      is available to drive it, but nothing is animated yet.

### Round 34 — Variant permeation
- [x] `analysis/variant_structures.py` runs the same pore / wetting / permeation
      pipeline over every deposited human entry, so a difference between them
      cannot come from a difference in treatment. Reachable as
      **Analysis → Variant structures…**
- [x] *Validate:* **the intended comparison is not available, and that is the
      result.** Three measured facts, each pinned by a test:
      1. **Every deposited human PIEZO1 structure is closed.** Bottlenecks
         0.67–0.93 Å against a 1.38 Å cation. All conductances are exactly 0,
         so a *difference* in conductance cannot be measured.
      2. **Three of the four variant entries do not contain their variant.**
         A1988 is unmodelled in both entries named for A1988V (8ZU8, 8YFC) and
         E756 is unmodelled in 9VMX. Only **8YFG (R2456H)** shows its mutation —
         HIS there, ARG in every other entry.
      3. **8ZU3, 8YFC and 9VMX share one model.** Byte-identical protein
         coordinates (31,839 atoms, same hash, 0.000 Å RMSD) across three
         separate depositions with different titles and different file
         checksums — verified not to be a download artefact.
- [x] **Coverage, stated plainly.** 4 deposited variant entries → **1** resolves
      its own mutation → **1** informative, against **68** curated variants (39
      with a direction). And all four are **gain-of-function**: there is *no*
      deposited loss-of-function structure, so this route cannot discriminate
      direction even in principle. The same data limit as Round 22, from the
      structural side.

### Out of band — everything computable reachable from the GUI
- [x] Audited on request. `permeation` and `interactions` were in the shared
      `ANALYSES` registry and wired into the CLI but absent from every menu —
      invisible to a GUI user. `ui/result_dialog.py` and
      `ui/tabular_analyses.py` add **Analysis → Ion permeation / Interactions /
      Variant structures / HaloTag labelling / HaloTag geometry**, each with a
      tooltip and a caveat shown above the numbers.
- [x] Two help topics added — *HaloTag and ion current*, *Framing and multiple
      structures* — since the guide had stopped at Round 30 and therefore
      misdescribed the application.
- [x] Three tests keep it from decaying: every registry analysis must have a GUI
      entry point, every menu action must carry a tooltip, and the guide must
      still mention what the recent rounds added **and** what they cannot do.

### Round 35 — The calcium nanodomain at the tag
- [x] `physics/nanodomain.py` implements the screened Green's function, taking
      exactly the two numbers Rounds 31 and 33 produce: unitary current
      **2.46 pA** (11ZC) and tag distance **3.95 nm** (8YEZ). Screening length
      **λ = 148 nm**, so at the tag the exponential is ≈1 and the answer is set
      by geometry, not by buffering. Reachable as **Analysis → Calcium
      nanodomain…** and `python -m piezo1.cli nanodomain 8YEZ`.
- [x] **The prediction holds.** **113.8 µM at 3.95 nm** — half the ~200 µM the
      roadmap expected, because Round 31 moved the tag closer than the 4–6 nm
      assumed, but the same order and far above the sensor's 0.2 µM Kd either
      way. Occupancy **99.82%**: the sensor is **saturated whenever its own
      channel opens**. So published puncta brightness reports *how many tags are
      labelled and how often the channel opens*, not calcium amplitude —
      which, joined to Round 32's result that a saturating protocol labels all
      three tags, means brightness heterogeneity points at **unreactive tags and
      open probability**, not at sub-saturating dye or graded calcium.
- [x] *Validate:* an **80-combination** sweep over tag distance (2–20 nm),
      calcium fraction (0.5–20%) and buffering (10 µM–10 mM). **78 of 80 stay
      saturated.** The two exceptions need 20 nm *and* 0.5% calcium *and*
      ≥1 mM buffer simultaneously. Across Round 31's full envelope
      (1.74–7.89 nm) the range is 55.6–263 µM, occupancy 99.64–99.92%; across
      Round 33's conductance range (16–94 pS), 99.55–99.92%.
- [x] **What would falsify it**, as numbers: the tag would have to sit at
      **373 nm** (≈100× further), or calcium carry **4.4×10⁻⁵** of the current
      (≈1000× less), or free buffer reach **0.14 M** (≈1400× physiological).
- [x] Found on the way: `resting_occupancy` — 100 nM resting calcium against a
      0.2 µM Kd already holds the sensor **33%** occupied, so its dynamic range
      is 33–100%, not 0–100%. Any occupancy target below that floor is
      unreachable at any distance, and the solver now says so.
- [x] **A silent frame bug, caught and pinned.** The report entry detected the
      C3 axis on the *unframed* structure and applied it to the *framed* one, so
      the pore was measured along a line that misses the pore. It reported the
      **closed** 8YEZ as carrying **32 pA** and making a 1.5 mM nanodomain —
      every number finite and plausible. Fixed, with a test asserting a closed
      structure never reports its own current.

---

## Block I — after the Round 30 review  *(added 2026-08-06)*

*Numbered from 36; Block H already claims 31–35 for the HaloTag work,
and the two blocks are independent — either can be worked first.*

**Where thirty rounds leave us.** The physics chain is closed, every link
validated against an independent published number, and — since Round 18 — every
link also checked by a route that does not reuse its own derivation. What the
last ten rounds mostly did was find out how much the project did *not* know:

- Round 18 overturned a Round 3 headline by 3.5×.
- Round 20 showed the Round 7 null could only ever have excluded a large effect.
- Round 22 added a second null, with its explanation recorded in advance.
- Round 26 found the mechanical ΔΔG's blindness was **algebraic** — a rank-one
  product — and repaired it (within-position variance 4.9% → 52.5%).
- Round 27 tripled the loss-of-function class but left a medium effect out of
  reach: 46 variants against the 104 needed.
- Round 29 found the dome radius and the published value are statistically
  **indistinguishable**, and that the gating overlap's third digit was never
  meaningful.
- Round 30 found a cross-check that failed because the *checking* route was the
  invalid one.

The through-line has not changed since Block G: every round that looked hard at
a recorded number found it stated with more confidence than it had earned. The
difference now is that the machinery to catch that — claims, parameters,
intervals, cross-checks — exists and runs on every commit.

**The destination remains untested at adequate power.** That is the honest
summary of thirty rounds, and the next block should be judged on whether it
changes that.

### Round 36 — Pre-register the third variant test
- [ ] Rounds 26 and 27 changed both sides of the question: the predictor can now
      distinguish substitutions, and the set has 46 directional missense
      variants at |δ| ≥ 0.41 detectable. That is a genuinely new situation and
      warrants a new pre-registration — **written and committed first**, under
      `docs/NEGATIVE_RESULT_PROTOCOL.md`.
- [ ] Declare confirmatory-for-a-large-effect or exploratory **before** looking,
      and state how the two evidence levels (measured vs disease-mechanism) are
      handled — pooling them silently would let 20 inferred labels outvote 26
      measured ones.

### Round 37 — Cross-check the remaining chain
- [ ] Round 30 covered three links. Do the rest: pore radius without Apollonius
      maximisation, SASA without Shrake–Rupley, conservation without pairwise
      alignment, PCA without SVD.
- [ ] Any disagreement is either a bug or a diagnosed approximation, and both
      are worth the round.

### Round 38 — Systematic error, not just statistical
- [ ] Round 29 attached intervals but stated plainly that none captures model
      error. Estimate it where possible: fit the dome with a spheroid as well as
      a sphere, run the ANM with two spring models, measure the pore with two
      probe conventions.
- [ ] Report the model spread beside the sampling interval, and say which
      dominates.

### Round 39 — The GUI reaches the variant pipeline
- [ ] The variant work is the project's central claim and is still CLI-only.
      Surface the evidence levels, the power statement and both null results
      where a user meets a variant.

### Round 40 — Reproduce a published figure end to end
- [ ] Pick one figure from Haselwandter & MacKinnon 2018 or Young 2023 and
      regenerate it from this codebase. Agreement is a strong integration test;
      disagreement is a finding.

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

### Out of band — seeing the HaloTag fusion
- [x] Rounds 31–32 computed the fusion and the labelling but drew nothing; the
      tags were reachable only through the CLI. `ui/fusion_controller.py` now
      draws them under **View → HaloTag fusion**: the three tag bodies, the
      linker seams in a colour the channel never uses, the accessible-volume
      cloud, and the dyes the labelling model says are bound. Everything is
      drawn so as to read as a model — sphere-of-gyration bodies, straight
      seams, and the envelope shown precisely so one sphere is not mistaken for
      a determined position. Figure: `docs/img/halotag_fusion.png`.
- [x] Corrected: the Round 31/32 display options were reported as living under
      *View*; they had in fact been added to *Options*. Moved to *View*, which
      is where they belong and where they were documented.


---

## Block J — review after Rounds 31–35

Written after five rounds of tag, labelling, permeation and variant work. What
those rounds actually established, and what it implies for the destination.

**The structural route to direction is empirically blocked, not just unproven.**
Round 22 found there were not enough phenotyped variants; Round 34 found there
are not enough *structures* — one informative variant structure, all
gain-of-function, every entry closed. Both ends of the comparison are
data-limited, and no amount of method work moves either.

**Three headline numbers now carry an explicit "this is tuning, not prediction"
label**: the fusion distance (robust to its assumption, so trustworthy), the
labelling curve (imported and exact), and the conductance (spans 16–94 pS across
unmeasured parameters, so not an independent prediction). Keeping those three
distinguishable from one another is worth more than improving any of them.

### Round 36 — Where the model *can* be tested
- [ ] The one comparison Round 34 leaves open: **8YFG (R2456H) against 8YEZ
      (wild type)**, both human, both closed, mutation resolved in one. Report
      bottleneck, wetting score and blocking mechanisms as a *paired* structural
      comparison, and state plainly that n = 1 supports no inference.
- [ ] *Validate:* whether the R2456H structure differs from wild type by more
      than the wild-type entries differ among themselves — the only control that
      makes a single pair interpretable.

### Round 37 — A predictor that could survive its own data limit
- [ ] Round 26 raised within-position variance from 4.9% to 52.5%. Ask the
      question that follows: given 39 directioned variants, what effect size is
      now detectable, and does the substitution-aware predictor reach it?
- [ ] *Validate:* against `design.minimum_detectable_effect`, and **do not run
      the comparison** unless the pre-registration protocol is followed first.

### Round 38 — The LoF gap, addressed rather than lamented
- [ ] Loss-of-function variants are absent from the structures but present in
      the curated set. Test whether they are also structurally distinguishable
      *in the wild type* — do LoF positions differ from GoF positions in burial,
      conservation, or coupling to the gate?
- [ ] *Validate:* this is a position-level question, so it is vulnerable to the
      exact confound Round 7 died of. Pre-register, and report the
      between-position variance share alongside any result.

### Round 39 — Provenance of the whole chain
- [ ] One command that walks a claim from the figure back to the file, the
      parameter set and the commit. `verify_claims` checks numbers; this checks
      that the *path* to each number is reconstructible.

### Round 40 — What a user should not be able to do
- [ ] Audit the UI for ways to produce a confident wrong number: analyses run on
      a cross-species overlay, a modified registry left unmarked, a companion
      structure mistaken for the primary. Round 33's menu audit found real gaps;
      this is the same exercise pointed at correctness rather than reachability.


---

## Block K — external resources not yet used  *(added on request)*

An audit of what this project pulls from the internet against what exists. It
fetches RCSB, AlphaFold DB, UniProt, PubChem, Europe PMC, ClinVar, Ensembl
orthologs and the CHAP grid. What follows is what it does **not**, ordered by
what would most move the destination — predicting direction from structure.

**The binding constraint is still data, and two of these attack it directly.**

### Round 41 — gnomAD, for the missing loss-of-function direction
- [ ] Round 34 established there is no deposited loss-of-function *structure*,
      and Round 22 that there are too few phenotyped *variants*. gnomAD v4 has
      ~800k exomes: pLoF variants, constraint scores (LOEUF, pLI, mis_z) and
      per-residue missense depletion. It cannot give a *measured* direction, but
      **regional constraint** is a phenotype-adjacent signal that exists for
      every residue rather than for 39.
- [ ] *Validate:* does missense depletion in gnomAD separate the curated GoF
      from LoF positions? Pre-register first — this is the same position-level
      confound Round 7 died of.

### Round 42 — MD trajectories other people have already run
- [ ] Two sources make this cheap. **MemProtMD** (Oxford) hosts coarse-grained
      and atomistic simulations of membrane proteins including PIEZO1, with
      lipid-contact occupancies already computed. **GPCRmd/MDDB** and Zenodo
      carry deposited PIEZO trajectories from published papers.
- [ ] *Validate:* compare their lipid-contact occupancies against this
      project's curated PIP2 cluster and the pockets found geometrically. An
      independent method agreeing is worth more than a better version of ours.

### Round 43 — the ligands that have no structure
- [ ] `ligands.json` is still 📋 and the registry contains **no** Yoda1-, Jedi-,
      Dooku1- or GsMTx4-bound entry, because none has been deposited. What does
      exist: **ChEMBL / PubChem BioAssay** dose-response data, **BindingDB**
      affinities, and the published mutagenesis that defines the Yoda1 pocket.
- [ ] Build `ligands.json` from those with provenance, and state plainly that
      every binding site in it is **inferred from mutagenesis and geometry, not
      from a bound structure**.

### Round 44 — predicted structures beyond AlphaFold2
- [ ] The project already fetches AlphaFold DB. **AlphaFold3 / Boltz-2 / Chai-1**
      predict complexes and ligands, and the **AlphaFold Protein Structure
      Database** now carries PAE matrices this project does not read. PAE is the
      honest way to say which inter-domain distances the prediction actually
      constrains — directly relevant to the unresolved distal blade.
- [ ] *Validate:* does PAE-weighted confidence agree with where the hybrid model
      seam had to be placed?

### Round 45 — the electrophysiology that is already published as data
- [ ] **IonChannelGenealogy**, **Channelpedia** and the supplementary tables of
      the PIEZO1 mutagenesis literature carry T50, τ_inact and conductance per
      mutant. Round 22 needed 104 variants for a medium effect and had 46.
      Harvesting published supplementary tables is the only route to that
      number that does not require new experiments.
- [ ] *Validate:* every harvested value must pass the same wild-type-residue
      gate the curated set already uses, and carry its PMID and its recording
      conditions — a T50 from a different preparation is not comparable.
