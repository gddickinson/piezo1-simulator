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
- [x] *Deferred to a later round:* expose the pore profile in the GUI.
      **Done** — Analysis dock → Pore tab, drawn by `ui/profile_plot.py` with
      radius against hydrophobicity and click-to-locate. The checkbox was
      simply never updated when the panel landed.

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
- [x] Particle animation whose flux is **set by** the computed current, with the
      HUD stating what the frame rate is in real time. The morph clock's
      discipline applies.
      **`render/flux.py` (Qt-free) + `ui/ion_flux_controller.py`, under View →
      Ion flux animation.** The measured fact that makes the label
      non-optional: a channel passes **1.5 × 10⁷ ions/s**, so a watchable
      stream runs about **10⁶× slower than real time** and the HUD says so.
      The morph clock refuses a seconds axis because an interpolation is not a
      trajectory; here the current genuinely *is* a rate, so the honest move is
      the opposite — name the factor rather than avoid it.
- [x] *A shut pore animates nothing.* Gated by the wetting verdict, so the GUI
      cannot reach a different conclusion from the headless pipeline. Measured:
      **8YEZ shows no ions** and states why (score 0.82 > 0.55 cutoff,
      bottleneck 0.095 nm, occluded *and* hydrophobic); **11ZC gives 2.44 pA**.
      Drawing a trickle through a closed gate would contradict the project's own
      structural result while looking like a demonstration of it.
- [x] *The rate declares the disagreement it inherits.* The solver gives 41 pS
      against a published 25–30 — a recorded result, reported not tuned — so the
      animation runs ~1.5× fast and the HUD states that too. An animation
      calibrated to a number 1.5× too large, shown without saying so, is exactly
      the confident-wrong-picture failure Round 50 audited for.
- [x] *Two faults in my own work, both instructive.* The controller first read
      `result.current_pA`, which does not exist — and `PermeationResult.current`
      is in **amperes**, so had the field existed the rate would have been wrong
      by 10¹². The test that was supposed to cover this only inspected the
      controller's *source* for `predict_wetting`; it passed throughout, then
      broke when the logic moved somewhere testable. The physics now lives in
      `timebase_for_structure`, exercised on real 8YEZ and 11ZC coordinates.

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

### Round 36 — Pre-register the third variant test — **and the third null**
- [x] `docs/PREREGISTRATION_ROUND36.md` written and committed **first, in its
      own commit** (`af37a82`), before any comparison ran. One primary endpoint,
      a six-test secondary family under Benjamini–Hochberg including a negative
      control, the decision rule, the inclusion criteria and five caveats — all
      fixed in advance. V598M excluded in writing beforehand as the one variant
      whose curated and inferred directions disagree.
- [x] Declared **confirmatory for a large effect, exploratory below it** (84%
      power at δ = −0.43, 50% at −0.28). The two evidence levels are **not
      pooled silently**: the primary uses the combined set, the document says
      why in advance, and measured-only is a pre-declared secondary so a reader
      can see whether the weaker labels carried it.
- [x] **Executed. THIRD NULL.** Primary Cliff's δ = **−0.249**, CI
      **[−0.628, +0.151]**, p = **0.405**, AUROC 0.625, on 19 GoF vs 15 LoF.
      Direction is as hypothesised; significance is not. **Fail to reject H₀**,
      and per the protocol the predictor is *not* adjusted and re-run.
      Nothing in the secondary family survives correction (min q = 0.591).
      Recorded in `docs/VALIDATION_ROUND36.md` in the prescribed order.
- [x] **The predictor did improve, visibly.** Substitution-aware δ = −0.249
      against the volume-only control's −0.025 on the same 34 variants — tenfold
      larger, consistent with Round 26. Across the three tests the effect has
      grown monotonically (−0.083 → −0.211 → −0.249). That is **suggestive and
      not evidence**: at δ = −0.25 roughly **130 variants** would be needed for
      80% power, against 34 available.
- [x] Two defects found and fixed: `design.sample_size_for` returned `max_n` for
      any *positive* effect size (sign convention against a one-sided
      alternative), so "how many variants for a large effect?" got a confident
      wrong answer that looked like a finding; and the AUROC call passed two
      groups to a function taking scores plus a boolean mask, yielding `nan`.
- [x] **FoldX could not be run at all** — 0 of 34 variants carry a value in the
      offline cache. Recorded as untestable rather than dropped, since "could
      not be run" and "was not significant" are different statements.

### Round 37 — Cross-check the remaining chain
- [x] `analysis/crosscheck_methods.py` re-derives all four by routes sharing no
      machinery with the pipeline. Run with `python scripts/crosscheck_methods.py`.
      Each alternative is first tested against a case with an analytic answer,
      because two routes agreeing only means something if neither is guessing.
- [x] **PCA — exact.** Power iteration on XᵀX (no LAPACK eigensolver, never
      forming the 7389×7389 covariance) reproduces the SVD eigenvalue to
      **0.0%** with **|cos| = 1.000000** between the two PC1 directions.
- [x] **SASA — 0.1%.** Shrake–Rupley's 256-point golden spiral gives
      11707.9 Å² on the 4RAX cap; Monte-Carlo with 4000 independent random
      directions gives 11692.0 Å².
- [x] **Pore radius — 5.2%, and the sign is the informative part.** The
      pipeline's polar-grid-plus-pattern-search gives 0.9300 Å at the 8YEZ
      bottleneck; 20k uniform random probe centres give **0.9783 Å**. A brute
      force can only match or beat a local optimiser, so the larger value means
      the pattern search **stops slightly short** — under-convergence, not a
      wrong answer. Bounded by a test so a real regression widens the gap.
- [x] **Conservation — correlation 0.817, and the residual is a bias in the
      *alternative*.** k-mer anchoring with no dynamic programming and no gap
      penalties agrees at invariant positions (**0.993** where the pipeline says
      1.00) but reads **0.653** where the pipeline says below 0.50, and its
      floor is **0.36 rather than 0**. The cause is selection: anchoring by
      *maximum exact matches* preferentially lines up residues that agree, so it
      inflates conservation exactly where a position is variable. Needleman–
      Wunsch uses gap penalties and a substitution matrix rather than raw match
      counts and is not subject to it. **The k-mer route is the weaker
      instrument** — the same verdict Round 30 reached about the parabola.
- [x] My first explanation of that residual (that it would concentrate near
      indels) was **wrong**: the eight worst positions all have coverage 1.00.
      Diagnosed properly and the module docstring corrected, rather than leaving
      a plausible story in place.

### Round 38 — Systematic error, not just statistical
- [x] `analysis/model_error.py` adds a fourth kind of spread, deliberately a
      separate type from `Sensitivity`: changing a spring exponent is a knob,
      changing a sphere into a spheroid is a different claim about the object.
      Run with `python scripts/model_error.py`.
- [x] **Dome — model error dominates 6×, and this is the round's headline.**
      A sphere gives **9.454 nm** (geometric rmse 6.180 Å); an oblate spheroid
      fits *better* — 5.243 Å, as it must with one more parameter — with
      flattening **+0.431** and apex curvature **14.991 nm**. Model spread
      **5.54 nm (58.6%)** against a bootstrap interval of **0.92 nm**:
      **6.0×**. The published confidence interval measures how well a sphere is
      determined, not whether a sphere is the right shape.
- [x] **Elastic network — 5.2%.** Cumulative gating overlap is 0.890 (uniform),
      0.912 (inverse_square, what the project reports) and 0.937
      (inverse_sixth). All three find the transition; the network is far less
      model-sensitive than the dome geometry.
- [x] **Pore — a null with a mechanism.** Apollonius and a uniform probe agree
      **exactly** at 1.70 Å because 7WLT's bottleneck lining is carbon; off that
      radius the gap is precisely the offset (0.300 Å at both 1.40 and 2.00 Å).
      So this is not a fixed systematic error but a restatement of the probe
      radius, and the per-atom refinement buys nothing at a carbon-lined
      constriction. Proven live rather than silently returning one number.
- [x] **The spheroid fitter was wrong first, and would have been reported as
      science.** Version one alternated centre and semi-axes with a hand-rolled
      gradient step: on a *known* spheroid it returned a = 163 for a true 100
      and c = 98 for a true 60, both inflated by the same 1.63× — the signature
      of a drifting centre. Replaced with the exact linear solution (null vector
      of the implicit quadric design matrix), which recovers known shapes to
      **0.01 Å**. Caught only by testing the fitter on knowns before using it.
- [x] Every result carries "this is a **lower bound**" — two models disagreeing
      bounds model error from below; two agreeing does not bound it from above,
      since both may be wrong the same way.

### Round 39 — The GUI reaches the variant pipeline
- [x] `analysis/prediction_record.py` holds the central claim's record as
      **data, Qt-free**, so the GUI, the CLI and the tests read the same numbers
      and cannot drift. Three tests, three nulls: Round 7 (δ −0.083, p 0.234),
      Round 22 (δ −0.211), Round 36 (δ −0.249, p 0.405).
- [x] **Evidence level now appears beside every variant** in the Annotation
      panel, colour-coded: `measured` (electrophysiology) in green,
      `disease_mechanism` (inferred from the disease) in amber, with the
      sentence explaining the difference. A classification alone reads as a
      fact, and for 20 of the 46 directional variants it is inferred.
- [x] **The one conflicting variant is flagged where it is met.** V598M shows
      "sources disagree: curated says GoF, the disease mechanism implies LoF —
      this project reports the disagreement rather than resolving it". Round 36
      excluded it in writing beforehand; a user should see why.
- [x] **Analysis → Variant prediction record…** shows all three tests, the power
      statement, and five standing caveats — including the one that says what
      the score may still legitimately be used for (finding mechanically coupled
      positions) rather than only what it cannot do.
- [x] `verify_record()` re-reads the stored Round 36 run and fails if the frozen
      numbers drift; measured agreement is **1.2×10⁻⁴**. The same discipline
      `analysis.claims` applies to the documented numbers.
- [x] The GUI-reachability guard written in Round 34 caught this round's new
      analysis the moment it entered the registry without a menu entry — the
      second time it has paid for itself.

### Round 40 — Reproduce a published figure end to end
- [x] Chose **Young et al. 2023's four-state tension response**, because their
      full rate set is published and registered, so the output can be checked
      against **two other papers** rather than against the numbers the model was
      built from. `scripts/reproduce_young2023.py` →
      `docs/img/young2023_response.png`.
- [x] **Half-activation AGREES, and strongly.** Young's rates through this
      project's solver give **T₅₀ = 2.711 mN/m** against Lewis et al. 2015's
      measured **2.7 ± 0.1** — a **0.4%** difference, from three independent
      sources (their rates, our solver, a third group's measurement).
- [x] **Inactivation DISAGREES by 8.5×.** τ at 5 mN/m is **73.3 ms** against Bae
      et al. 2013's **8.6 ± 0.4 ms**. The decay is cleanly mono-exponential (a
      bi-exponential fit adds nothing), so this is not a fitting artefact. k2,
      the O→I₁ rate, carries the timescale: at 8 s⁻¹ it sets ~125 ms before the
      rest of the system pulls it to 73. Reaching 8.6 ms needs a **12.8×**
      increase in k2, to ~103 s⁻¹.
- [x] **This justifies a policy the project already had.** `kinetics.wt_tau_ms`
      already carried the note that mutants are calibrated by *fold change*
      against the wild-type τ, "never by absolute τ across preparations". That
      was written before this measurement; the measurement is what makes it more
      than caution.
- [x] `calibrate_k2_for_tau` **refuses** an out-of-reach target and states the
      reachable range (13.3–211.2 ms at the default bounds) rather than
      returning its search bound — checked, because clipping would have looked
      like an answer.
- [x] An API I misread, pinned so nobody repeats it: `with_modification` takes
      **fold changes, not absolute rates**. Passing `k2=8.0` to a model whose k2
      is already 8 gives **64**. Reading it as a setter produced τ = 13 ms where
      the model really gives 73 — a plausible number, wrong by roughly the very
      factor this round measures.

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

Written after five rounds of tag, labelling, permeation and variant work.
*Numbered from 46: Block I already claimed 36–40 and Block K 41–45. The
collision was mine, and the blocks are independent — any can be worked
first.*
 What
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

### Round 46 — Where the model *can* be tested
- [x] `analysis/paired_variant.py` measures 8YFG (R2456H) against the wild-type
      entries, in the canonical frame, through the same pore and wetting
      pipeline. Reachable as **Analysis → R2456H vs wild type…** and
      `python -m piezo1.cli paired_variant`.
- [x] **The control was widened, not narrowed.** The comparison is against
      *three* independent wild-type entries (8YEZ, 8ZU3, 8ZU8), not one.
      8YFC and 9VMX are excluded **by coordinate fingerprint** — they are
      byte-identical to 8ZU3, and including them would have added two
      zero-difference pairs, shrinking the wild-type spread and making the
      variant look more distinct than it is.
- [x] *Validate:* **R2456H is not structurally distinguishable.**

      | measure | variant | wild-type range | WT spread | largest variant difference |
      |---|---|---|---|---|
      | bottleneck | **0.808 Å** | 0.673–0.930 | 0.257 | **0.135** |
      | wetting score | **0.904** | 0.457–0.986 | 0.529 | **0.446** |

      The variant falls *inside* the wild-type range on both, and its largest
      difference from any wild-type entry is **smaller** than the largest
      difference between two wild-type entries. Tested generously — *any*
      measure exceeding the spread would have counted — and it still does not.
- [x] **Why that is unsurprising once stated.** Every deposited human entry is
      closed, and R2456H's phenotype is *slowed inactivation*. A closed
      structure has no obligation to show a gating defect. That is part of the
      result rather than an excuse for it.
- [x] The control is itself tested both ways: a synthetic variant far outside
      the wild-type set is detected, one just inside is not, and a wild-type set
      of *identical* structures makes even a 0.001 difference "distinguishable"
      — which is exactly what including the duplicates would have moved towards.
- [x] **n = 1**, stated in the summary, the caveat and the note. This says what
      the deposited structures show, not what R2456H does.

### Round 47 — A predictor that could survive its own data limit
- [x] Round 26 raised within-position variance from 4.9% to 52.5%. Ask the
      question that follows: given 39 directioned variants, what effect size is
      now detectable, and does the substitution-aware predictor reach it?
      **No — and no reachable dataset would.** Round 26 did help: the effect
      grew −0.083 → −0.249 and the requirement fell from >800 variants to 134.
      But the optimistic ceiling is **59** (46 directional + Round 45's 35
      harvest candidates, × the 74% modelling-gate survival), where the minimum
      detectable effect is **0.356** against the observed **0.249** and power is
      **0.51** — a coin flip. `analysis/feasibility.py`, docs/FEASIBILITY_ROUND47.md.
      *(The roadmap's "39" is stale; the count is 46 directional, 34 tested.
      The conclusion holds for any of the three.)*
- [x] *Validate:* against `design.minimum_detectable_effect`, and **do not run
      the comparison** unless the pre-registration protocol is followed first.
      **No comparison was run.** Every effect size is read from
      `prediction_record.VALIDATION_RECORD` rather than recomputed, and a test
      asserts the module imports no statistic that could produce a fresh one.
      Both numbers are now guarded in the claims registry.

### Round 48 — The LoF gap, addressed rather than lamented
- [x] Loss-of-function variants are absent from the structures but present in
      the curated set. Test whether they are also structurally distinguishable
      *in the wild type* — do LoF positions differ from GoF positions in burial,
      conservation, or coupling to the gate? **No — the fifth null, and the
      flattest.** Primary (relative SASA, 14 LoF vs 16 GoF): Cliff's δ **+0.036**,
      CI [−0.384, +0.473], p = 0.509, AUROC 0.482 — all three decision clauses
      fail. Nothing in the six-endpoint BH family survives; smallest q = **0.930**.
      Gate distance separates *exactly* nothing (δ = +0.000).
      docs/VALIDATION_ROUND48.md.
- [x] *Validate:* this is a position-level question, so it is vulnerable to the
      exact confound Round 7 died of. Pre-register, and report the
      between-position variance share alongside any result.
      **Pre-registered in its own commit (`7ffb008`) before the run, with the
      ceiling stated up front rather than as a closing caveat.** The variance
      share is now *measured*, not asserted: between-position **1.000000**,
      within-position **0.000000**. R2456's four variants — R2456H/K/P (GoF) and
      R2456C (LoF) — all take the identical value 0.127326, because the feature
      never sees the substitution. Against 4.9% (Round 7) and 52.5% (Round 26),
      a wild-type positional feature has **0%**, so a positive result could
      never have become a variant-direction predictor.
      The pre-registered **negative control** (distance from the C3 axis,
      δ = +0.268) has a **larger effect than every mechanistic endpoint** — so
      the spread across endpoints is noise, exactly as in Round 41.

### Round 49 — Provenance of the whole chain
- [x] One command that walks a claim from the figure back to the file, the
      parameter set and the commit. `verify_claims` checks numbers; this checks
      that the *path* to each number is reconstructible.
      **`make provenance` / `python -m piezo1.analysis.provenance_chain`.** Five
      links per claim — document, code, parameters, data, commit — with the
      parameter and data links **measured while the claim runs** rather than
      declared. **It found two real defects on its first two runs.**
- [x] *Defect 1, in the checker itself:* it reported Round 22's δ as missing
      from `VALIDATION_ROUND22.md`. The document writes `−0.211` with U+2212
      MINUS SIGN, and the pattern only matched the ASCII hyphen. Fixed, and
      pinned — along with a guard that an en-dash range (`2.7–4.7 mN/m`) is
      *not* read as a negative.
- [x] *Defect 2, the real one:* **26 of 101 registered parameters were read by
      no code at all.** They appear in the parameters dialog with a unit and a
      citation, an override on them is recorded, reports carry the non-default
      banner, and `verify_claims` refuses to run — while the number does not
      move. Proved on `pore.step`: overriding 1.0 → 0.25 left the 8YEZ
      bottleneck bit-identical at 0.951756 Å. This is strictly worse than an
      unregistered literal, which is at least honestly invisible, and the
      parameter audit cannot see it — the audit checks a literal is *declared*
      to correspond to a parameter, and **a declaration is not a wire**.
- [x] *Fixed here:* all five `pore.*` parameters, end to end — the callee
      resolves at call time and the six callers that passed literals no longer
      do. Overriding `pore.step` now moves the bottleneck (0.9518 → 0.7649 Å)
      and the default is unchanged to every digit. `analysis_pore` had been
      using `step=1.5` while the registry advertised 1.0; it now uses the
      registered value. Unwired count 26 → 21, ratcheted by a test.

### Round 49b — The other twenty-one dead parameters
- [x] Wire the remaining 21, in the same way and with the same proof: override
      it, show the number moves, show the default is unchanged.
      **All 21 wired across 11 modules and 28 call sites. Unwired count 21 → 0:
      every one of the 101 registered parameters is now read by code.** Every
      documented number is unchanged — 759 tests pass and all 18 claims verify
      with zero drift.
- [x] *Validate:* the ratchet is now `test_no_registered_parameter_is_read_by
      _nothing`, asserting **zero**, plus `test_parameter_effect` which
      *empirically* proves 11 representative parameters move a number and
      restore exactly. None were deleted: each had a real call site.
- [x] *A defect the wiring itself introduced, caught by the empirical probe.*
      `value()` returns a float, so `n_permutations` arrived as `10000.0` and
      numpy rejected it outright. Eleven count-valued parameters now cast at
      the point of resolution. **Static wiring would have called these done** —
      the probe is what found it, which is the argument for the probe.
- [x] *A pre-existing bug the probing uncovered.* `ConservationProfile.
      top_conserved` sorted residues failing the coverage filter to the bottom
      with `−inf` but **still returned them** when fewer than `n` qualified —
      carrying their real conservation value with nothing to mark them, and
      reachable from the CLI as `conservation --top`. With 1 residue passing at
      coverage ≥ 0.99 it returned 10, of which 9 failed. Fixed to return *at
      most* n, and pinned.

### Round 50 — What a user should not be able to do
- [x] Audit the UI for ways to produce a confident wrong number: analyses run on
      a cross-species overlay, a modified registry left unmarked, a companion
      structure mistaken for the primary. Round 33's menu audit found real gaps;
      this is the same exercise pointed at correctness rather than reachability.
      **`piezo1/ui/hazards.py` — ten hazards as data, each with the scenario,
      what would be wrongly concluded, and the guard.** Every one is driven by a
      **positive control** in `test_hazards.py`: the dangerous situation is
      constructed and the guard asserted to fire, because a guard nobody has
      watched fail is not evidence.
- [x] *Of the three named suspects, one was already guarded and two were not.*
      The cross-species overlay refusal was real and fires. The other two were
      open: **a result window named neither the structure the numbers came from
      nor the parameter set they were computed under.** With companions
      displayed there was nothing on the window saying which structure it was —
      the exact "companion mistaken for the primary" failure — and an
      overridden registry produced numbers that looked documented. Both are now
      stamped on every result window, recorded **at compute time** so a
      non-modal window cannot silently agree with a registry that moved after it.
- [x] *One further gap:* `CAVEATS["interactions"]` was the empty string, so the
      interaction inventory was the one tabular analysis displayed with no
      warning at all. It now states that contacts are those of *this* structure
      in *this* state, and that criteria are heavy-atom based because deposited
      entries carry no hydrogens. A test forbids any analysis being shown
      without a caveat.
- [x] *Reaches the user, not just the suite:* the register generates a tenth
      help topic, "Wrong numbers, and what stops them", rather than being
      re-written as HTML that could drift from it.


---

## Block K — external resources not yet used  *(added on request)*

An audit of what this project pulls from the internet against what exists. It
fetches RCSB, AlphaFold DB, UniProt, PubChem, Europe PMC, ClinVar, Ensembl
orthologs and the CHAP grid. What follows is what it does **not**, ordered by
what would most move the destination — predicting direction from structure.

**The binding constraint is still data, and two of these attack it directly.**

### Round 41 — gnomAD, for the missing loss-of-function direction
- [x] `analysis/gnomad.py`: cached, offline-tolerant GraphQL client following
      `analysis/external.py`. 6,708 missense variants placed on 2,521 residues,
      none unplaced. Pre-registered in `docs/PREREGISTRATION_ROUND41.md`,
      committed in its own commit (`eb9ae1f`) before any comparison.
- [x] **The gene-level answer arrived first and it is discouraging** — recorded
      in the pre-registration rather than after: **LOEUF 1.10**, **pLI 1.5×10⁻¹⁰⁰**,
      `oe_mis` **1.45**, `mis_z` **−11.3**. PIEZO1 is not merely unconstrained
      but missense-**enriched**. Looking for regional depletion inside such a
      gene was always a long shot, and saying so beforehand is what makes the
      null interpretable.
- [x] *Validate:* **FOURTH NULL.** Primary Cliff's δ **−0.269**, CI
      **[−0.595, +0.074]**, **p = 0.0477**, AUROC 0.634 on 18 LoF vs 24 GoF
      positions. Medians differ in the hypothesised direction (LoF 2.275 vs GoF
      2.520 missense/residue) but the interval contains zero. Nothing in the
      five-test secondary family survives BH (min q = 0.140).
- [x] **The three-part decision rule earned its place.** p = 0.0477 is *below*
      threshold. A rule written as "p < 0.05" alone would have produced this
      project's first positive result, on a p-value 0.0023 under the line, with
      an interval comfortably containing no effect. The conjunction — p, sign,
      **and** an interval excluding zero — was fixed in advance for exactly this.
- [x] **The pre-registered negative control is the informative row.** The raw
      per-residue count, declared beforehand as something that "should show
      nothing", gives δ **−0.231, p 0.078** — indistinguishable from the ±25
      predictor. The windowing that makes this a *regional* estimate contributes
      nothing, so the predictor cannot be told apart from its own control. That
      is a stronger reason to disbelieve the primary than the interval alone.
- [x] Position **2456 excluded in writing beforehand** as the only residue
      carrying both directions (R2456H/K/P gain, R2456C loss) — simultaneously
      the reason for the exclusion and the cleanest demonstration that a
      position-level predictor has a ceiling.
- [x] Named the biggest methodological weakness for a successor: observed counts
      conflate "few variants seen" with "few expected". gnomAD's own regional
      missense constraint model fits expected counts and would separate them.

### Round 42 — MD trajectories other people have already run
- [x] **The premise was wrong, and that is the result.** "Two sources make this
      cheap" does not survive contact with the data. `analysis/external_md.py`
      implements the *check* rather than the comparison, so the conclusion is
      reproducible and will change by itself when the situation does.
- [x] **MemProtMD holds 1 of 21** catalogued PIEZO entries — only **3JAC**, the
      2015 structure. Absent are 7WLT, 7WLU, 6B3R and 8YEZ, i.e. every structure
      this project actually uses. Measured with a **working control**: 2RH1 and
      1M0L return 200 on the same probe, so the absence is about PIEZO and not
      about the request.
- [x] **The one entry cannot answer the question.** 3JAC resolves 918 of 2,547
      residues (36%), and of the 15 curated lipid-associated residues it
      resolves **4** — the PIP2 cluster in full and **none** of the three blade
      basic clusters. A simulation of a model that omits the lipid-binding
      residues cannot report their lipid contacts.
- [x] The other two named sources do not help either: Zenodo's PIEZO1 records
      are microscopy TIFFs and PDFs rather than trajectories, and GPCRmd is
      GPCR-specific — PIEZO1 is not a GPCR.
- [x] Also recorded: MemProtMD's analysis is **browsable but not fetchable**.
      The site is a single-page app and no public API endpoint could be found;
      ten candidate paths returned 404 against working page URLs.
- [x] A test asserts that an **offline run returns "not checked" rather than an
      empty coverage** — otherwise a network failure would silently manufacture
      exactly this round's conclusion.

### Round 43 — the ligands that have no structure
- [x] `resources/ligands.json` built by `scripts/build_ligands.py` behind a
      provenance gate, with the table split into `scripts/ligand_table.py` so
      the data can be read without the machinery. Six modulators: Yoda1, Yoda2,
      Jedi1, Jedi2, Dooku1, GsMTx4.
- [x] Chemistry fetched from PubChem and **verified**: the build compares the
      returned InChIKey against the recorded one, so a wrong CID cannot pass
      silently. Yoda1 C13H8Cl2N4S2 (CID 2746822), Dooku1 C13H9Cl2N3OS,
      Jedi1 C12H10O3, Jedi2 C10H8O3S, Yoda2 C16H9Cl2KN2O2S2; GsMTx4 is a
      peptide, recorded by UniProt Q7YT39.
- [x] **Every binding site is labelled as inferred, and the build verifies the
      claim rather than asserting it.** `deposited_modulators()` scans the
      heteroatoms of all 21 downloaded structures: nothing outside lipid,
      detergent, glycan and ion codes is present, so no bound modulator exists
      in any of them. If one is ever deposited, the build **fails** and the
      resource is correctly marked out of date.
- [x] Site evidence is a graded field — `mutagenesis`, `docking_md`,
      `geometric`, `none` — and `bound_structure` is **rejected by the build**.
      Only **one** of six carries a residue-level site: Yoda1's
      1718/2075/2078, at `docking_md`, from MD rather than from contact.
      Each of the other five records *why* it has no site, so silence cannot
      read as "not looked at".
- [x] *Validated* against the project's own ground-truth table: Yoda1
      **EC50 26.6 µM** (Syeda 2015) and GsMTx4 **Kd 155 nM** (Bae 2011), both
      matching, with every citation required to resolve in `references.json`.
- [x] Reachable as **Analysis → Modulators…** and
      `python -m piezo1.cli ligands`, with the caveat shown above the numbers.

### Round 44 — predicted structures beyond AlphaFold2
- [x] `analysis/prediction_confidence.py` reads the **PAE matrix** the project
      has been downloading and ignoring. pLDDT says how well a residue's *local*
      environment is predicted; PAE says how well residue *i* is placed when the
      model is aligned on *j* — the only one that answers a hybrid model's
      question. AlphaFold DB v6, 2521×2521, max 31.75 Å.
- [x] *Validate:* **half yes, half no, and the "no" is the more useful half.**
      **pLDDT agrees with the seam**: the unresolved distal blade (1–569)
      averages **64.5** against the core's **74.2**, with **52.2%** of blade
      residues below 70 versus **27.0%** of the core.
      **PAE does not single out the seam.** The raw block comparison looks
      decisive — 27.3 Å across versus 16.1/20.7 Å within — but that conflates
      "across the seam" with "far apart in sequence". Controlled for separation
      the penalty is at most **+4.3 Å on a 31.75 Å scale**, and at 50–150
      separation it **reverses**: cross-seam pairs score **13.25** against
      **15.82** within.
- [x] **The stronger finding.** PAE is **85% saturated** beyond 800 residues of
      separation — and **80% saturated within the cryo-EM-resolved core alone**,
      a region experiment places confidently. AlphaFold does not determine
      PIEZO1's long-range architecture *anywhere*. So for `hybrid.py` the seam
      is **not** the weak point: the global arrangement is unconstrained
      wherever the cut is made, which argues for placing the distal blade by the
      experimental C3 symmetry and dome geometry rather than trusting the
      prediction's relative placement at all.
- [x] The separation control is itself tested against a **planted** seam
      penalty (recovers 6.0 as 6.5) and against a matrix with none (residual
      under 1.0), because "no penalty found" and "the control does not work"
      would otherwise be indistinguishable — the same discipline as Round 42's
      probe control.
- [x] Not pursued: AlphaFold3 / Boltz-2 / Chai-1. The PAE result says the
      limitation is the prediction's global architecture rather than its
      vintage, so a newer predictor is not obviously the fix and would need its
      own confidence readout to be worth trusting.

### Round 45 — the electrophysiology that is already published as data
- [x] `analysis/harvest.py` scans the 38 open-access JATS full texts the project
      already downloads, extracts candidate substitutions with the sentence they
      appeared in, gates them on the wild-type residue, and resolves which
      numbering system each is in.
- [x] *Validate:* **the gate is not a formality — 23% of raw hits fail it.**
      The funnel, measured:

      | stage | n | lost |
      |---|---|---|
      | raw regex hits, 15 papers | 86 | |
      | pass the wild-type gate | 66 | −20 (cDNA changes written like protein ones) |
      | mappable to human numbering | 66 | −0 |
      | not already in the curated 68 | 35 | **−31 already held** |
      | carry an extractable measurement | **2** | −33 prose only |

- [x] **The premise does not hold.** Round 36 needed ~130 directional variants
      and had 34. This harvest adds **2**, and **neither carries a direction**.
      The bottleneck is not the gate: the numbers live in *prose* and in
      non-open-access supplements, not in the machine-readable tables the round
      assumed. Across all 38 papers the tables contained **four** substitution
      strings, two of them cDNA.
- [x] **The existing curation is better than assumed** — 31 of 66 gated
      candidates are already in the curated 68. That bounds what any harvest of
      this corpus could ever have added.
- [x] **40 of 66 are mouse-numbered**, against 18 human. Conversion goes through
      the alignment in `core.sequence`, and a test asserts the offsets are *not*
      constant, so the harvest genuinely exercises the numbering map rather than
      accidentally passing on a fixed shift.
- [x] **No direction is ever assigned.** `Candidate` has no direction field, and
      a test enforces its absence. Reading "slowed inactivation" out of prose
      and calling it gain-of-function would put unreviewed labels into the set
      the blind tests depend on — the one thing that must stay hand-checked.


---

## Block L — review after Rounds 35–39

**What these five rounds were actually about: knowing what the numbers are worth.**
None of them added a new capability to the model. Rounds 35–39 measured the
*trustworthiness* of what was already there, and in four of the five the
instrument being checked turned out to be the thing that was wrong.

- Round 35's nanodomain prediction **held**, robustly — the one positive result.
- Round 36 ran the third pre-registered test and returned the **third null**,
  and found `sample_size_for` returning a confident wrong answer.
- Round 37 cross-checked four methods; the **k-mer conservation route** was the
  biased one, not the pipeline.
- Round 38 estimated model error and found the **spheroid fitter** wrong before
  it found anything about the dome — and then that the dome interval is **6×
  too narrow**.
- Round 39 surfaced all of it where a user meets a variant.

**The recurring failure mode is worth naming.** Three times in five rounds, an
*alternative* built to check the pipeline was itself defective, and each time it
would have produced a plausible number rather than an error: a 15-16 nm tag
distance, an 89% dome model error, a 32 pA current through a closed channel.
What caught all three was calibrating the alternative against a known answer
before pointing it at the unknown. That should be a standing rule, not a habit.

### Round 51 — Calibrate every alternative before it is believed
- [x] Audit the cross-check and model-error modules for any route not tested
      against an analytically known case first. Add the missing calibrations.
      **Audited 8 modules and 42 public checking callables. Four had no
      known-answer case** and are now calibrated:
      `model_error.spring_model_error` (plant a mode from one spring model —
      it recovers 0.990 against 0.107 and 0.459 for the other two);
      `design.minimum_detectable_effect` (the δ it returns must actually
      deliver its stated power — measured 0.798 and 0.803 against a 0.80
      target, closing a loop that had only ever been checked against this
      project's own recorded results); `uncertainty.sensitivity` (the range
      must be exactly the min and max of known outputs — it had only been
      tested on its *wording*); `validation.permutation_test` (against
      exhaustive enumeration of all C(8,4) = 70 partitions — exact 0.01429 vs
      sampled 0.01515, and the (r+1)/(n+1) convention asserted to err
      conservative rather than merely tolerated).
- [x] Write the rule into `CLAUDE.md`: a checking instrument is a measuring
      instrument, and an uncalibrated one is worse than none because its
      disagreement looks like a finding.
      **Written, with the four instances that motivated it, and made
      enforceable:** `tests/test_calibration.py` holds the register, and
      `test_every_checking_instrument_has_a_calibration` fails if an
      instrument is added without one while `test_named_calibrating_tests_exist`
      fails if the named test does not exist — which caught two test names I
      had guessed wrongly on its first run.
- [x] *A reproducibility finding, checked before being reported as harmless.*
      Writing the spring calibration exposed that `ANM.calc_modes` is **not
      run-to-run deterministic on near-degenerate systems** — identical inputs
      gave overlaps from 0.954 to 0.997, because ARPACK starts from an unseeded
      random vector. The obvious worry is `anm.gating_overlap`, a documented
      claim. Measured over four runs it is **bit-identical** (0.70482207), since
      the real structures have well-separated low modes and only the artificial
      random geometry is degenerate. So no documented number is affected and
      the fix belonged in the test, which now asserts the *separation* between
      spring models rather than an absolute threshold.
- [x] *The audit's own instrument needed calibrating, which is the round in
      miniature.* My first keyword scan reported 12 instruments as
      uncalibrated. Reading the tests showed it was missing calibrations named
      in the test *name* rather than the body — `test_auroc_known_cases`,
      `test_cliffs_delta_extremes`. Corrected, the real count was four.
      An uncalibrated audit would have produced eight false findings.

### Round 52 — Widen the intervals that Round 38 showed are too narrow
- [x] The dome radius is quoted with a sampling interval 6× smaller than its
      model spread. Decide what to publish: a wider interval, a stated model
      caveat, or both — and apply the same question to the footprint, T₅₀ and
      the gating overlap.
      **Decided: both, by a stated rule — publish the *widest* term and name
      its kind.** `analysis/published_interval.py` holds the rule and all four
      numbers with every term measured:

      | Quantity | Value | Published | Kind | Narrowest would have been |
      |---|---|---|---|---|
      | Dome radius | 9.72 nm | **[9.45, 14.99] nm** | model form (lower bound) | **18.5× too tight** |
      | Gating overlap | 0.705 | **[0.554, 0.723]** | cutoff sensitivity | 3.6× |
      | T₅₀ | 2.711 mN/m | **[2.584, 2.838]** | rates at ±20% | **15.9×** |
      | Footprint energy | 25.27 k_BT | **[25.27, 26.94]** | κ = 20–25 k_BT | 1.0× |

      **T₅₀ was the surprise**: the two solvers agree to 0.6%, so the number
      looked well determined, but its input rates make it 16× less certain than
      the numerical route suggests. The measured 2.7 ± 0.1 still lies inside, so
      the agreement with Lewis 2015 survives the input uncertainty rather than
      depending on the exact published rates — a stronger statement than before.
- [x] *Validate:* `verify_claims` must still pass, which means the documented
      numbers and their stated uncertainties have to move together.
      **18 claims, 0 drift.** The point estimates deliberately did not move —
      only what is claimed *about* them — so claim tolerances (which detect
      code drift) and published uncertainty (a scientific statement) stay
      separate questions.
- [x] *A mismatch found while doing this, kept rather than quietly repaired.*
      The dome's model comparison is anchored on the **untrimmed** sphere fit
      (9.45 nm) while the published number is **trimmed** (9.72 nm) — the two
      were never like-for-like. Measured: the registered `geometry.sphere_trim`
      moves the radius by only 0.30 nm across 0–0.25, against a 5.54 nm model
      spread, so the conclusion is unaffected. Recorded as its own term.

### Round 53 — The tour should end on the record, not the mechanism
- [x] The guided tour ends on two null results. There are now three, plus a
      quantified model error and a data limit with a number attached. Rewrite
      the closing steps so a student leaves knowing what the project does *not*
      know, with figures.
      **The roadmap item was itself stale: there are now *five* nulls, not
      three, and the tour's closing step still said "tested twice".** Rewritten
      as three closing steps (11 steps → 13):
      **11 · The central claim, and five attempts on it** — every
      pre-registered test with its effect size, and the sharpest result: a
      wild-type positional feature has *exactly zero* within-position variance,
      so R2456H/K/P/C receive the identical value and no such feature could
      ever assign a direction. Figure: forest plot, every interval crossing zero.
      **12 · Why more data would not settle it** — 134 variants needed against
      a ceiling of 59, and the distinction between "we need more data" and "the
      data that could exist is not enough". Figure: reachable versus required.
      **13 · What this application cannot do** — Round 52's model error where a
      student meets the dome number (9.72 nm, but [9.45, 14.99] across shapes,
      ~18× the bootstrap interval), plus the clinical disclaimer.
- [x] *Figures:* `scripts/make_record_figure.py` writes both from **recorded**
      results, recomputing nothing frozen. `TourStep` gained `image` /
      `image_caption`, and a missing PNG degrades to prose rather than raising —
      the figures are git-ignored regenerable outputs, so a fresh clone must
      still be able to take the tour.
- [x] *No number is stated in prose.* The closing measurements read
      `ALL_PREREGISTERED` and the claims registry, so a sixth test updates the
      tour instead of leaving it stale — which is precisely how it came to say
      "tested twice" after five tests had run. A test asserts `_data_limit`'s
      source contains neither 134 nor 59.
      `prediction_record` gained `OTHER_PREREGISTERED` (Rounds 41 and 48) and
      `ALL_PREREGISTERED`; `VALIDATION_RECORD` stays scoped to the ΔΔG score
      the GUI shows beside a variant, so its caveat does not start describing a
      number the user is not looking at.

### Round 54 — Make the data limit actionable
- [x] Every route to the central claim now ends in "not enough phenotyped
      variants". Block K §41 (gnomAD constraint) and §45 (published
      supplementary tables) are the two that need no new experiments. Cost them
      honestly and do the cheaper one.
      **Both named routes were already spent** — Round 41 ran gnomAD constraint
      (null, with the negative control indistinguishable from the predictor) and
      Round 45 harvested the supplementary literature (35 candidates, **none**
      directional, only 2 with any measurement). So the costing is of what is
      left, in `analysis/data_routes.py`:

      | Route | Yield | Cost | Status |
      |---|---|---|---|
      | Population constraint | +0 | already spent | done |
      | Literature harvest | +0 | already spent | done |
      | Within-position pairs as they stand | **+1** | none | open |
      | Curate the variants one label from a pair | **+3** | 3 literature reads | open |
      | Admit the engineered variants | +15 | a scientific decision | blocked |

- [x] *"Not enough variants" is now a named list rather than a lament.* Exactly
      **three** variants would each unlock a new within-position pair:
      **M870V** (position has M870I, LoF), **R1358C** (has R1358P, GoF) and
      **A2020V** (has A2020T, GoF). That is the cheapest route with any yield
      and the only one that is a finite list rather than a search — but the
      yield is an **upper bound**, since two of the three are curated as VUS
      *precisely because* the evidence to direct them was not found. Even at
      full yield this reaches four positions, which is not a design.
- [x] *The important result is a correction to the Round 50 review, not a new
      route.* That review counted 40 positions carrying more than one variant
      and called it "a real design". Filtered to what the design needs — two or
      more **missense** variants, each **directional**, with no **source
      conflict** — the count is **one** (R2456), which is what Round 48 said.
      The 40 included nonsense variants (Q1009\*), insertions (E2496ELE),
      positions whose second variant has no direction, and V598M, which is
      curated GoF and ClinVar LoF. That conflict was already being reported by
      `variant_sets.disagreements()`, so the existing machinery was right and
      the review was not.

### Round 55 — Retire what does not earn its place
- [x] The codebase has grown to ~100 modules and 651 tests. Find the analyses
      nothing depends on, that no round cites, and that no test pins to a
      result — and delete them. A smaller project that is entirely load-bearing
      is easier to trust than a large one that is mostly scaffolding.
      **Measured across 109 modules and 568 top-level definitions: the project
      is already load-bearing at module level — every module is imported,
      tested or a documented entry point.** At definition level, **five** things
      earned nothing and were deleted: `_poisson_newton_step` (42 lines, the
      Gummel-loop route abandoned for the electroneutral limit),
      `PredictionContext` (a dataclass never constructed), `permutation_p`,
      `distinct_colors` and `_published` (dead since Round 25). Suite 822 green
      after removal.
- [x] *The hard part was the detector, not the deleting, and two versions of it
      would have done real damage.* A grep over `__all__` reported **102**
      unused public names — including `format_result`, used inside its own
      module, and every return-type dataclass, which is constructed but never
      named elsewhere. An AST version that collapsed same-file references into
      a set reported **129** dead functions including `fetch_pdb`, `cmd_list`
      and `_optimise_slice`; acting on it would have deleted the CLI. Round 51's
      rule is why neither was believed.
- [x] *`piezo1/dead_code.py` is the standing guard*, registered in
      `test_calibration.py` and calibrated on every run: known-used names must
      not be flagged, and a planted unreferenced name must be. It counts
      *occurrences* rather than files (the bug that hid same-file calls) and
      counts bare words in string literals, because this project dispatches by
      string through `ANALYSES` and the CLI.
- [x] *One deletion needed a check first.* `_poisson_newton_step` carried a real
      numerical finding in its docstring — the Gummel loop going −0.37 V →
      −171 V → −2×10¹⁶ V. It was removed only after confirming that divergence
      is recorded in the module, in `test_permeation.py` and in `SCIENCE.md`.
      Deleting the last copy of a finding would be worse than keeping dead code,
      and a test now pins that all three records survive.


---

## Block M — review after Rounds 41–45

**Five rounds, five answers about data rather than method.** Round 41 (gnomAD:
an unconstrained gene, fourth null), 42 (deposited MD: 1 of 21 structures), 43
(modulators: no bound structure exists for any of them), 44 (AlphaFold: no
long-range constraint anywhere), 45 (literature harvest: 2 usable candidates
against ~96 needed).

**The destination is now measurably out of reach with the available data, and
every route has been tried and costed.** That is a result. What the project
should not do is a fifth predictor on 34 variants.

**Two things these rounds did well and should be kept.** Every negative result
was measured with a *control* — a probe that finds known-present entries, a seam
penalty planted in a synthetic matrix, a negative-control endpoint pre-registered
before the test. And each was implemented as a **check that will change by
itself**: if MemProtMD ingests a PIEZO structure, if a Yoda1-bound entry is
deposited, if the corpus gains a supplementary table, the tests fail and the
conclusion is revisited without anyone remembering to.

### Round 56 — Say the conclusion once, at the top
- [ ] The four nulls, the data limits and the model-error result are spread
      across `SCIENCE.md`, four `VALIDATION_*.md` files and the roadmap. Write a
      single page — `docs/CONCLUSION.md` — that states what this project set out
      to do, what it established, and what it could not, with the numbers.
- [ ] Link it from `README.md` and the in-application help, so it is the first
      thing a reader meets rather than the last thing they assemble.

### Round 57 — Hand-curate the 35 fresh candidates
- [ ] Round 45 found 35 substitutions not in the curated set, each with its
      sentence and source. Curating them by hand is a bounded task with a known
      denominator, and it is the only remaining route that adds *measured*
      directions.
- [ ] *Validate:* how many of the 35 have a direction recoverable by a human
      from the sentence alone? That number decides whether a fifth test is ever
      possible.

### Round 58 — Retire the predictor, keep the coupling map
- [ ] Round 39 recorded that the score has a legitimate use — finding
      mechanically coupled positions — and an illegitimate one. Make that
      structural: rename the variant-impact output so it cannot be read as a
      direction prediction, and keep the coupling analysis.

### Round 59 — The tour and the README should end where the science does
- [ ] Both still present the project as pursuing the central claim. It has been
      tested four times. Rewrite the closing steps and the README summary to
      match Round 56's page.

### Round 60 — A reproducibility run from an empty clone
- [ ] `make reproduce` has never been run from a genuinely empty state in CI.
      Do it: fresh clone, `create_env.sh`, `python -m piezo1.io.fetch`, full
      suite, every figure, `verify_claims`. Anything that only works because of
      a stale cache is a reproducibility bug and this is how it surfaces.

---

## Review after Rounds 46–50

**What the five rounds established, in one line each.** Round 46: the one
variant structure is indistinguishable from wild type, with a control proving
the comparison could have detected a difference. Round 47: no reachable dataset
can resolve the effect the mechanical predictor produces — 134 variants needed,
59 the optimistic ceiling. Round 48: a fifth null, and wild-type positional
features have *exactly* 0% within-position variance. Round 49/49b: 26 of 101
registered parameters were inert while advertised, now all wired and each
proved to move a number. Round 50: two of three named UI hazards were open.

**The pattern worth naming.** Four of these five rounds found the defect in the
*instrument*, not the science: a spheroid fitter that would have reported 89%
model error, a checker that could not read its own documents' minus sign, a
registry whose parameters did nothing, a probe whose "no effect" came from
badly chosen coordinates. The standing habit — calibrate an alternative against
a known answer before believing its disagreement — has now paid for itself
often enough to be the project's most reliable rule.

**Where the destination stands.** Five pre-registered tests, five nulls, five
predictor families. Round 47 showed this is not a "more data" problem for the
across-position design: the data that *could* exist is not enough. So the
across-position route is closed, and saying so is a result rather than a
failure.

**~~But the within-position route is far more open than Round 48 suggested~~ —
CORRECTED BY ROUND 54, and the correction is the more useful result.** This
review counted 40 positions carrying more than one variant and called that "a
real design". Round 54 applied the filters such a design actually needs — two or
more **missense** variants at one position, each carrying a **direction**, with
no **source conflict** — and the count is **one**, which is what Round 48 said.

| Filter | Positions surviving |
|---|---|
| carrying > 1 variant of any kind (this review's figure) | 40 |
| … both missense (drops nonsense, insertions, deletions) | fewer |
| … both carrying a direction | 1 |
| … no curated-versus-ClinVar conflict (drops V598M) | **1** — R2456 |

The forty included nonsense variants (Q1009\*), insertions (E2496ELE), positions
whose second variant has no direction at all, and one where the two sources
disagree. A within-position comparison still removes the between-position
variance that consumed 99.8% of Round 7's predictor — the reasoning was right —
but the data to run it does not exist, and this review's optimism came from a
number that does not survive its own filters.

There is also an unused asset: **15 engineered variants, every one with a
measured functional effect**, two of them at positions shared with natural
variants (1335, 2117). They are excluded from every analysis set because
`engineered` is not `GoF`/`LoF` — which is correct for the disease question and
wasteful for the mechanism question.

---

## Block N — the within-position route (Rounds 61–65)

### Round 61 — How many shared positions would be enough
- [ ] The Round 47 question, asked of the design that is actually open: for a
      within-position comparison, what effect size is detectable at 40 shared
      positions, and how does that compare with the across-position 134?
- [ ] *Validate:* against `design.minimum_detectable_effect` and a paired
      statistic, and **run no comparison** — this is feasibility, as Round 47
      was. If the answer is that 40 positions suffice, that is the first time
      this project has had a testable route to its central claim.

### Round 62 — Direction at the shared positions
- [ ] The 40 positions are only usable if both variants at a position carry a
      direction. Measure how many do, at each evidence level, and what it would
      take to resolve the rest. Round 45 costed the literature harvest; this
      costs the much smaller, targeted question.
- [ ] *Validate:* report the count honestly, including if it collapses the
      Round 61 design. A route that looks open until the directions are counted
      is exactly the kind of optimism this project has learned to check early.

### Round 63 — The engineered variants, used or explicitly refused
- [ ] Fifteen engineered variants carry measured functional effects and are
      excluded from every analysis set. Decide, in writing, whether a
      conductance or selectivity change can stand in for gain/loss of
      mechanosensitive function — and if it cannot, record *why* rather than
      leaving them silently unused.
- [ ] *Validate:* if they are admitted, they enter as their own evidence level,
      never pooled with `measured`; `variant_sets` already enforces that.

### Round 64 — Pre-register the within-position test
- [ ] Only if Rounds 61–63 leave a design with adequate power. Written under
      `NEGATIVE_RESULT_PROTOCOL.md` §2, committed alone before anything runs,
      with the power statement up front as Round 48 did.
- [ ] *Validate:* the decision rule keeps the three clauses — Round 41 proved
      the interval clause earns its place, and Round 48 that stating the
      ceiling first stops a positive being over-read.

### Round 65 — Finish the modules still marked planned
- [ ] `structure/hybrid.py`, `physics/modes.py`, `analysis/contacts.py`,
      `analysis/variants.py` and `analysis/docking.py` have been 📋 for the
      whole project. Either implement them or delete the rows — an INTERFACE
      that promises five modules it does not have is the documentation
      equivalent of a registered parameter nothing reads.
- [ ] *Validate:* `make provenance` reports 5/18 chains complete; most breaks
      are legitimate, but the count should be explained rather than left as a
      number that looks like failure.

---

## Review after Rounds 51–55

**What the five rounds did.** 51: calibrated every checking instrument and put
the rule in `CLAUDE.md`. 52: decided what interval to publish beside each
headline number — the widest term, named for its kind. 53: rewrote the tour's
ending on the record, with figures. 54: costed every remaining route to more
phenotyped variants. 55: audited for scaffolding and found the project already
load-bearing.

**The pattern, and it has hardened into the project's main risk.** In *four* of
these five rounds the defect was in the checking apparatus, not the science:
an audit that missed calibrations named in test names (51), a model comparison
anchored on a differently-trimmed fit (52), a review that counted 40 positions
without asking what each row was (54), and two dead-code detectors that would
have deleted the CLI (55). Every one produced a plausible number rather than an
error. The standing rule — calibrate against a known answer, and suspect the
checker before the pipeline — has now caught six separate incidents and is the
single most valuable convention here.

**Two corrections to my own earlier writing**, both in this block of rounds:
the Round 50 review's "40 positions is a real design" (actually one), and the
Round 53 discovery that the tour still said "tested twice" after five tests.
Both were optimistic readings that survived because nobody re-derived them.
That is an argument for the claims registry and the ratcheting tests, not for
more prose.

**Where the science stands.** Five pre-registered tests, five nulls, five
predictor families. Round 47 closed the across-position route (134 variants
needed, 59 reachable). Round 54 closed the within-position route (one usable
position, at most four if three named variants could be directed). Both are now
measured costs rather than impressions. **The central claim is not merely
unproven; it has been shown to be unprovable with data that could exist.** That
is a result, and Block O should treat it as one rather than looking for a sixth
predictor.

---

## Block O — say what was established, and make it reusable (Rounds 66–70)

### Round 66 — The conclusion document
- [ ] `docs/CONCLUSION.md`: what the project set out to do, what it established,
      what it could not, with the numbers and the figures Round 53 built. Block
      M §56 proposed this before Rounds 47 and 54 existed; it can now state the
      *unprovability* result rather than a list of nulls.
- [ ] *Validate:* every number in it must come from the claims registry or the
      validation record, and a test must assert that — this project has twice
      shipped prose that went stale.

### Round 67 — What the negative result is worth to someone else
- [ ] The reusable output is not the predictor; it is the machinery that showed
      the predictor could not be validated: pre-registration discipline,
      `feasibility`, `data_routes`, `published_interval`, `calibration`. Write
      the short methods note that would let another structural-biology project
      apply the same test to its own central claim.
- [ ] *Validate:* it must be honest that this pipeline's main output was five
      nulls, and say why that is the point rather than a disclaimer.

### Round 68 — The engineered variants, decided
- [ ] Round 54 left 15 measured functional effects marked `blocked` on a
      scientific question: may a conductance or selectivity change stand for
      gain or loss of mechanosensitive function? Answer it in writing, with the
      literature, and either admit them at their own evidence level or record
      why not.
- [ ] *Validate:* `variant_sets` already refuses to pool evidence levels; if
      they are admitted, that refusal must still hold.

### Round 69 — Retire the planned modules
- [ ] `structure/hybrid.py`, `physics/modes.py`, `analysis/contacts.py`,
      `analysis/variants.py`, `analysis/docking.py` have been 📋 for the whole
      project. Round 55 showed everything that exists is load-bearing; these are
      the opposite problem. Implement or delete the rows.
- [ ] *Validate:* `dead_code.audit()` must stay at zero, and INTERFACE must not
      promise a module that does not exist.

### Round 70 — A fresh clone, start to finish
- [ ] `make reproduce` on a clean checkout, timed, with every step's failure
      mode recorded. The project claims reproducibility as aim A5 and has never
      measured how long it takes or what breaks first without a warm cache.
- [ ] *Validate:* report the wall-clock and the first thing that fails, rather
      than fixing forward until it passes.

