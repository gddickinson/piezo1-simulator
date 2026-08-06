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

### Round 26 — Substitution-aware mechanics
- [ ] The Round 7 diagnostic said the mechanical ΔΔG reports *position*, not
      *substitution* (99.8% between-position variance). Fix the cause rather
      than adding predictors around it: perturb the network by charge change,
      hydrogen-bonding capacity and proline backbone disruption, not volume
      alone, so that four substitutions at R2456 can differ.
- [ ] Success criterion fixed now: **within-position variance must exceed 20%**
      of total, measured on the four R2456 substitutions and the other
      multiply-substituted positions. If it does not, the approach has failed
      and that is the reported result.

### Round 27 — Expand the phenotyped variant set
- [ ] Round 20 is unambiguous: 42 variants for a large effect, 98 for a medium
      one, against 25 available. The binding constraint on this project's
      central claim is **data, not method**.
- [ ] Curate from ClinVar, the ProtVar cross-check and the primary
      electrophysiology literature, with the same wild-type verification gate
      the existing 68 passed. Record inter-curator ambiguity explicitly rather
      than resolving it silently.
- [ ] Report the achieved n and recompute the minimum detectable effect. If it
      still cannot reach a medium effect, say so.

### Round 28 — Nonlinear footprint in the gating energetics
- [ ] Round 18 built the elastica solver but only `DomeModel` consumes it. The
      two-state model's ΔA and the footprint contribution to T₅₀ still use
      linear numbers that are known to be 3.5× too large.
- [ ] Propagate the nonlinear areas through `dome.py` and re-derive T₅₀.
      Compare with the measured 2.7 ± 0.1 and 5.1 ± 0.2 mN/m and report the
      change, **including if the linear version happened to agree better** —
      a wrong model can fit a right number.

### Round 29 — Uncertainty on every reported quantity
- [ ] Dome curvature, pore radius, mode overlaps and ΔΔG are all reported as
      point estimates. Add intervals: bootstrap over atoms for the sphere fit,
      over structures for the ensemble PCs, over network cutoff for the ANM.
- [ ] A number without an interval invites exactly the overconfidence Rounds
      18–20 kept finding.

### Round 30 — Adversarial review of the whole chain
- [ ] Re-derive each headline result by a deliberately different route and
      record where the two disagree: dome curvature without sphere fitting,
      mode overlap without superposition, T₅₀ without the Markov scheme.
- [ ] The Round 18 lesson generalised: the useful check is the one that does
      not reuse the derivation being checked.

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
