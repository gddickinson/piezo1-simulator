# SESSION LOG

Running record of what was done and — more importantly — *why*. Newest first.

---

## Session 1 — 2026-08-05

### Goal
Establish the project: research PIEZO1 thoroughly, acquire all structural and
sequence data, choose the technology stack, and build the foundation layers.

### Environment decisions

**Conda environment `piezo1` (Python 3.11.15).** Created by
`scripts/create_env.sh`. Scientific core from conda-forge (numpy 2.4.6,
scipy 1.17.1, numba 0.66, MDAnalysis 2.10, mdtraj 1.11.1, OpenMM 8.5.2,
pdbfixer, biotite 1.4, Biopython 1.87, scikit-image 0.26, networkx 3.6);
GUI/GL and structural-bioinformatics layer from pip (PyQt6, moderngl 5.12,
PyOpenGL, pyqtgraph 0.14, ProDy 2.6.1, pydssp, freesasa). All 21 imports
verified working on macOS ARM.

**Renderer: moderngl + QOpenGLWidget on PyQt6, OpenGL 4.1 core.**
This was the highest-risk decision, so it was de-risked first with a probe
before any renderer code was written. Result on this machine (Apple M1 Max):
context reports `4.1 Metal - 90.5`, and a fragment shader writing
`gl_FragDepth` compiles. That last point is what matters — it means ray-cast
**impostor** rendering is available, so spheres and cylinders can be drawn as
screen-space quads with per-pixel-exact geometry instead of tessellated meshes.
That is how PyMOL, VMD and ChimeraX get their speed, and it is the difference
between 120k atoms being interactive and being a slideshow.
Rejected: VTK (heavy, harder to style, awkward instancing), pyqtgraph.opengl
(too limited for custom shaders), Qt3D (immature Python bindings), embedding a
web renderer (loses direct access to our numpy arrays).

### Data acquired

- **UniProt**: human Q92508 (2521 aa, verified) and mouse E2JF22 (2547 aa),
  FASTA + full JSON. 38 transmembrane segments annotated per protomer, plus
  topology, PTMs, disulfide C2411–C2415, coiled coil 1339–1368, and 26 natural
  variants with disease annotations.
- **RCSB**: a targeted query (by UniProt accession *and* by entity description,
  rather than free text, which returned 236 mostly irrelevant hits) found
  **28 PIEZO entries**. 21 downloaded as mmCIF. The find that matters most:
  **human PIEZO1 structures now exist** — 8YEZ (3.3 Å apo), 8ZU3 (3.1 Å
  PIEZO1–MDFIC), and, remarkably, **cryo-EM structures of three disease
  variants**: 8YFG (R2456H), 8ZU8/8YFC (A1988V) and 9VMX (E756del). Earlier
  work in this field was almost entirely mouse.
- **Curved/flattened pairs** suitable as morphing endpoints: 7WLT/7WLU (mouse in
  bilayer, 2022) and 11YE/11ZC (mouse in plasma-membrane vesicles, 2026).
- **AlphaFold DB**: `AF-Q92508-F1-model_v6`. Note the v4 URLs 404 — the current
  version is **v6** (created 2025-08-01) and must be discovered via the API
  endpoint, not guessed. This model covers residues 1–2521, which matters
  because no experimental structure resolves the distal blade below ~570.

### Code written

`config.py`, `io/cif_reader.py`, `core/structure.py`,
`structure/superpose.py`, `structure/geometry.py`,
`scripts/build_uniprot_annotations.py`.

**Why a custom mmCIF reader rather than Biopython.** Biopython allocates a
Python object per atom. For a 34k-atom PIEZO1 trimer — and we will load several
at once, plus morph trajectories — that is the wrong shape entirely. The custom
reader walks the file once into contiguous numpy arrays and takes ~0.6 s.
Selections then become boolean masks and coordinate maths becomes vectorised,
which is what both the renderer and the elastic-network model want anyway.

*Bug worth remembering:* the first version of the tokenizer treated only space
and tab as whitespace, so the trailing newline of each `_atom_site` row became a
22nd token and silently shifted every subsequent row by one column. It surfaced
as `invalid literal for int(): 'ATOM'`. Whitespace handling in hand-written
parsers deserves a test — one is planned in `tests/test_cif_reader.py`.

### Scientific validation achieved

The dome-geometry pipeline was checked against the literature rather than
assumed correct. Taking the mid-point of each of the 38 transmembrane helices in
each protomer as a sample of the mid-membrane surface, recovering the three-fold
axis, and fitting a sphere gives:

| Structure | State | Radius of curvature | Dome depth | Excess area |
|---|---|---|---|---|
| 7WLT | curved, bilayer | **9.7 nm** | 4.9 nm | 256 nm² |
| 11YE | curved, PM vesicle | **10.4 nm** | 4.6 nm | 293 nm² |
| 8YEZ | human apo | 12.0 nm | 5.8 nm | 279 nm² |
| 8ZU3 | human + MDFIC | 12.5 nm | 5.4 nm | 270 nm² |
| 7WLU | flattened | 18.4 nm | 2.5 nm | 379 nm² |
| 11ZC | flat, PM vesicle | 21.6 nm | 3.5 nm | 335 nm² |

Published value for the closed state is **10.2 nm** (Haselwandter & MacKinnon
2018, eLife), and ~11.8 nm outside-in (Vaisey & MacKinnon 2026). Our 9.7 and
10.4 nm land squarely there, and the curved→flat contrast is unambiguous. C3
axis recovery is exact (120.00°, 0.00 Å RMSD). **This is now the standing
regression test for the geometry pipeline.**

*Known caveat:* the curved and flattened entries resolve different residue
ranges (7WLT 784–2547, 7WLU 576–2547), so their footprint radii are not directly
comparable. A fair comparison must restrict both to the commonly resolved
residues — to be implemented alongside the morphing code.

### Research

Six parallel literature-research agents were dispatched, writing dossiers to
`ref/research/`. The first wave hit an API session limit and was relaunched
against REST endpoints (Europe PMC, UniProt, RCSB, PubChem) instead of web
search, which had also been exhausted. The lipid-modulation dossier completed
and is a strong result: it supplies the membrane parameters (κ = 20–25 k_BT,
footprint decay λ = 14 nm, ΔA values), ligand potencies, and the important
2026 finding from Vaisey & MacKinnon that mechanical force alone is *not*
sufficient to gate PIEZO1 — a specific lipid cofactor is also required.

### Next

Sequence numbering map (human↔mouse), the annotation layer, then the elastic
network model, then the renderer.

---

## Variants & disease dossier — `ref/research/04_variants_disease.md`

Written 2026-08-05. 68-entry curated JSON variant table (22 GoF, 17 LoF, 8 VUS,
6 blood-group, 15 engineered) plus narrative. Every wild-type residue in the
table was validated programmatically against `ref/sequences/Q92508_human_PIEZO1.fasta`
— zero mismatches.

**The human↔mouse numbering map now exists** (Biopython global alignment,
BLOSUM62, 82.47% identity) and is tabulated in §8.2. The headline result is that
**the offset is not constant**: −6 at E756, −5 through the beam, **+16** across
THU9/anchor/outer-helix, and **+26** from the cap onward. Validated against the
known anchor mouse R2482 ≡ human R2456 (Ma et al. 2018). Traps to remember:
human E2496 ≠ mouse E2496 (mouse E2496 = human **E2470**); human S2446 = mouse
S2472; human E2117 = mouse E2133; human S1335 = mouse S1330.

**Two literature errors found and documented:** the Open Biology 2025 review
gives T2127's mouse equivalent as 2142 (correct: **2143**; the local alignment is
gapless), and Albuisson 2013 gives R1358P as c.4072G>C (ClinVar: **c.4073G>C**).
Also note a WebFetch summarisation hazard — the summariser fabricated a mouse
numbering column for that review; the raw XML had to be parsed directly to catch
it. Prefer parsing tables from XML over summarised fetches for numeric data.

**Blocking finding for the renderer/annotation layer:** the mutated residue is
often *absent* from its own structure. Verified against the local CIFs — all six
human PIEZO1 entries model chain A from residue **570** only, and **E756 is not
modelled in 9VMX**, nor **A1988 in 8ZU8 or 8YFC**. Only R2456 (8YFG) is present.
Residues 2060–2521 are continuously modelled, so the whole pore module renders
cleanly; G253, D669, G718, E756, C1064, K1877 and A1988 need AlphaFold and must
be flagged model-only in the UI.

Population genetics: gnomAD v4 gives pLI ≈ 0 and LOEUF 1.097 — PIEZO1 is
unconstrained for heterozygous LoF, exactly as expected for recessive LMPHM6.
The missense Z of −11.3 is *not* usable: synonymous Z is −14.1, so the mutational
model fails at this locus. ClinVar (2858 records) yields 133 sequence-level P/LP;
truncating alleles are spread uniformly while missense P/LP density is 6–8×
enriched in the pore module (inner helix 5.88 per 100 aa vs blade 0.69).

---

## Session 1, part 2 — renderer, GUI, validation

### The result that matters

An elastic network model built from the **closed** structure alone reproduces
the experimentally observed gating transition. Comparing curved 7WLT with
flattened 7WLU over the 1274 residues common to all six protomers (19.7 Å
trimer RMSD):

| | |
|---|---|
| Best single mode overlap | **0.705** (mode 3, symmetry A, collectivity 0.610) |
| Cumulative overlap over 40 modes | **0.964** |
| Best E-mode overlap | **0.0011** |
| Fraction of overlap² in A modes | **100.00%** |

The symmetry result is the part worth dwelling on. PIEZO1 is a C3 trimer, so
every mode carries an irreducible-representation label. Isotropic membrane
tension is itself C3-symmetric, so only A modes can couple to it at first
order. Every E mode scores essentially zero overlap with the real transition —
the analysis recovers the selection rule without being told about it. That is a
strong internal consistency check, and it means the app can tell a user which
modes are candidate gating coordinates on principle rather than by eye.

### Traps found by checking rather than trusting

**Protomer labels lie.** 7WLT and 7WLU label their three chains in opposite
rotational order around the symmetry axis. Superposing by chain label gave
71.2 Å RMSD instead of 19.7 Å, and the difference vector built from it was
meaningless — the first overlap calculation returned a misleading 0.213 with no
obvious sign of error. `match_protomers()` now always determines correspondence
by superposition, and a test pins it.

**Secondary structure came out 100% coil.** The distance criteria were fine;
the C-alpha pseudo-torsion had an inverted sign, giving −51° for the PIEZO1
inner helix where the IUPAC convention gives +51°, so the helix test never
fired. Fixed by correcting the cross-product order. PIEZO1 now assigns as 77%
helix / 10% strand / 13% coil, with OH, IH and the beam coiled coil all 100%
helix and the cap correctly β-rich.

**Camera framing by bounding sphere** left the molecule filling ~55% of the
viewport. PIEZO1 is a wide flat propeller, so the bounding-sphere radius badly
overestimates what needs to fit. `frame()` now projects into the current camera
orientation and solves for the exact containing distance.

**Disconnected networks.** Writing the ANM tests surfaced that a contact network
in several pieces contributes six rigid-body modes *per piece*. A model with a
detached fragment would have returned rigid-body motions as its lowest
functional modes — silent and thoroughly misleading. `ANM.n_components()` now
detects this and `calc_modes` discards `6 × n_components` by default.

### Research integration

All six literature agents completed on the second attempt (the first wave hit an
API session limit; relaunching them against REST endpoints rather than web
search worked). The most valuable outcome was **independent convergence**: an
agent's alignment of human and mouse PIEZO1 reproduced this project's offset map
block-for-block, and its Yoda1 pocket conversion (human A1718/A2075/A2078) and
selectivity residues (E2117/E2470) matched the values computed here exactly.

Two corrections were applied as a result. The "clasp" and "latch" domains were
**removed** — neither survived verification against primary sources, and
inventing boundaries would colour residues with false confidence. They were
replaced with two verifiable elements: the UniProt-annotated beam coiled coil
and the Piezo1.1 spliced segment. The anchor description was corrected:
P2113/F2114 is primarily the anchor apex brake on the inner helix, with the
cholesterol context secondary.

### Coverage honesty

Building the variant resource revealed something that shapes the whole UI: all
six human PIEZO1 structures model from residue 570 only, and **14 of 68 curated
variants — including the E756del malaria-associated allele — are resolved in no
human structure at all.** Only R2456 appears in its own structure (8YFG). The
viewer now greys those out, states the count, and warns on selection, rather
than highlighting nothing and letting the user assume it worked.

### Renderer

moderngl + QOpenGLWidget, OpenGL 4.1 core, ray-cast impostors. 31 599 atoms and
275k ribbon vertices render in 14–20 ms on an M1 Max. One Qt-specific trap worth
recording: `QOpenGLWidget` does **not** render to framebuffer 0, so moderngl
must be pointed at `defaultFramebufferObject()` every frame or nothing appears.

### Verification

43 tests, ~10 s. They pin the tokenizer's whitespace handling, the reversed-
handedness detection, dome curvature against the published 10.2 nm, the ANM
symmetry characters, the gating-overlap result itself, ten cross-species residue
equivalences (each also checked for matching amino-acid identity), and an
assertion that the numbering offset is *not* constant.

A scripted GUI smoke test (`scripts/screenshot_app.py`) drives the real
application and checks its outputs, so "the app still starts and computes" is a
test rather than a hope.

### Next

Helfrich membrane-footprint solver, tension-dependent Markov gating kinetics
(the Young et al. 2023 PNAS four-state model is fully parameterised in the
research dossier and is the one to implement), conformational morphing between
curved and flat endpoints, and pore-radius profiling.

---

## Session 1, part 3 — morphing, and the repo

`piezo1/structure/morph.py` interpolates between two experimental endpoints.
Three methods, each reporting its own geometric error so the cost of the
interpolation is visible rather than hidden:

| method | worst C-alpha bond error | note |
|---|---|---|
| linear | 2.94 A | the chord artefact |
| restrained | **0.00 A** | distances restored to interpolated targets |
| modal | 1.60 A | 30 modes capture 95.2% of the change |

The chord artefact is worth naming because it is not obvious: under
straight-line interpolation atoms cut chords through space, so C-alpha-C-alpha
distances contract wherever the local motion is rotational. PIEZO1's blades
swing through large arcs, so mid-path frames are measurably wrong. For
comparison, ProDy's `calcAdaptiveANM` was benchmarked on this same 3822-site
trimer at 18 minutes and still left bonds stretched to 5.39 A.

Dome geometry tracks along the restrained path — radius of curvature
9.2 -> 13.1 nm, dome depth 4.6 -> 2.7 nm, excess area 278 -> 316 nm2 — so the
physics follows the morph rather than being imposed on it.

**Refactor.** `main_window.py` reached 657 lines, past the project's 500-line
rule, so dome/mode handling moved to `ui/physics_controller.py` and morph
handling to `ui/morph_controller.py`; the window is now 364 lines. The split
introduced two bugs, both caught by re-running the scripted GUI test rather
than by assuming the move was safe: a guard still testing `hasattr` on an
attribute now initialised to `None`, and a method left with its pre-split name.
Mechanical refactors of GUI code need the smoke test run afterwards, every time.

**Final research agent** returned benchmarks that independently validated the
architecture: sparse ANM at N=7500 takes 4.4 s versus 311 s dense (71x); full
C3 block-diagonalisation would give only 1.76x on top of that, confirming the
decision to do symmetry *labelling* rather than symmetry-adapted solving;
ProDy's `imANM` is an O(N^2) Python loop and unusable at this scale;
MDAnalysis's `hole2` is an empty stub in 2.10 and HOLE has no arm64 build, so a
native pore profiler is required.

**Repository** published to github.com/gddickinson/piezo1-simulator (private),
eight commits, 86 tracked files. `ref/` and `data/` are git-ignored and fully
regenerable with `python -m piezo1.io.fetch`.

### State at end of session

Working: data layer, structure model, cross-species numbering, dome geometry,
elastic network models with C3 symmetry labelling, conformational morphing, the
renderer, and the GUI. 51 tests, ~15 s.

Next, in order of value: the tension-dependent Markov gating model (Young et al.
2023 PNAS four-state, fully parameterised in `ref/research/03a_kinetic_models.md`
— sigma_50 = 1.4, b = 0.8 mN/m, all rates given); the Helfrich membrane
footprint solver (1-D radial validates to 3.8e-3 against the exact K_0 in
0.9 ms, then revolve; lambda = 14 nm with kappa = 20 kT implies gamma =
0.42 mN/m); a native pore-radius profiler (a leash is mandatory — unconstrained
the probe sphere escapes to R = 6188 A); and the hybrid full-length model,
noting that AlphaFold **cannot** place residues 1-570 relative to the core
(PAE 25-29 A against a 31.75 A maximum), so PIEZO2 6KG7 is the better guide.

---

## Session 2 — 2026-08-05/06 — rounds, references, measurement, animation

Set up a ROADMAP.md organised into rounds, and a recurring 20-minute loop that
works the next unchecked one.

### Rounds completed

**Round 1 — pore geometry.** Closed human 8YEZ has a 0.76 Å bottleneck and is
non-conductive; flat 11ZC has 3.25 Å and is conductive. The profiler
rediscovered all three curated constrictions from coordinates alone: the V2450
hydrophobic gate (3.0 Å) and the CTD constrictions at M2467 (1.2 Å) and P2510
(1.4 Å). *The leash is a correctness requirement, not a convenience* — the
clearance function has no interior maximum, and an untethered probe escapes to
R ≈ 6188 Å, which is a true maximum and a useless answer.

**Round 2 — gating kinetics.** Young et al. 2023 four-state tension model.
Emergent half-activation **2.71 mN/m** against a measured cell-attached T50 of
**2.7 ± 0.1 mN/m**. Detailed balance exact to 1e-16, because C→I₁ is fixed by
microscopic reversibility rather than fitted.

**Rounds 11–12 — measurement and interactions.** The C2411–C2415 disulfide is
recovered in all three protomers at 2.04 Å. TM38 is the least tilted
pore-proximal helix at 6.9°. And **R2456 salt-bridges to E2117 of the
neighbouring protomer** in all three copies — the archetypal gain-of-function
residue pairing with the selectivity glutamate, a concrete structural route
from mutation to phenotype.

**Rounds 14–15 — animation.** Seven animations, offscreen-rendered.

### Errors caught, and what they teach

**Mutant calibration was inverted.** Presets originally solved for the rate that
reproduced an *absolute* measured inactivation τ. But Young's parameterisation
gives a wild-type τ of 35–80 ms while Bae's whole-cell measurement gives 8.6 ms
— different preparations. Calibrating R2456H to 22.2 ms therefore made it
*faster* than the model's own wild type, exactly inverting the phenotype of the
best-known gain-of-function variant. Fold changes transfer between
preparations; absolute time constants do not.

**Six citations were confidently wrong.** PMIDs entered from memory resolved
cleanly to unrelated papers — a Piezo1 structure citation came back as a
bone-marrow transplantation study, another as stem-cell reprogramming. Europe
PMC returns whatever the identifier points at. Every seed entry now carries an
`expect` keyword checked against the resolved title, and anything failing is
reported rather than written out. 51/51 now verify.

**Interaction criteria were too loose.** PLIP's 4.1 Å hydrogen-bond cutoff is
only valid *with* hydrogens and an angle test; on heavy atoms alone it produced
8005 "bonds" per trimer including donor–donor N···N pairs. Tightened to 3.5 Å
with N···N excluded except for histidine.

**A test claim was stronger than the data.** "TM38 is the least tilted helix"
is false — blade helices 50–60 Å out are also near-vertical. Narrowed to "least
tilted of the pore module", which is what the coordinates support.

The pattern across all four: the failure mode is never a crash, it is a
confident wrong number. Every one was caught by checking the output against
something independent rather than by the code raising an error.

### Bibliography

51 references resolved from Europe PMC into a committed JSON plus a generated
`docs/REFERENCES.md`, each recording what the project uses it *for*. 29
open-access full texts downloaded to `ref/papers/` (git-ignored — other
people's copyright, and the bibliography suffices to retrieve them).

### Operational note

The 20-minute cron cadence is shorter than a round actually takes, so
invocations stack up. Harmless — each fire simply picks up the next unchecked
item — but worth knowing when reading the history.

---

## Round 3 — membrane mechanics (2026-08-06)

Implemented the Helfrich footprint solver and the dome energetics. Four
published numbers reproduced: γ = 0.420 mN/m from λ = 14 nm and κ = 20 k_BT;
4.116 mN/m per k_BT/nm²; λ = 13.998 nm recovered from the solver's own output;
and T₅₀ = 4.99 mN/m from Cox's ΔG₀ and ΔA against their measured 5.1 ± 0.2.

The footprint result is the one that matters scientifically: around the
measured 7WLT dome it stores **622 nm² of excess area against the dome's own
256 nm²**. Haselwandter & MacKinnon argued the footprint dominates tension
sensitivity; this puts a number on it.

### Two errors that produced plausible numbers

**`L @ L` is not how you build a biharmonic operator.** Squaring the discrete
Laplacian squares its condition number. The solver converged — to a profile
with a 47 nm decay length where the exact answer is 14 nm, and a 59% energy
error that stayed at 59% under grid refinement. That last detail is the tell:
an error that does not shrink when the grid does is not a discretisation error,
it is convergence to the wrong problem. Rewritten as a coupled second-order
system, now second-order convergent.

**The closed-form energy had the Bessel ratio upside down.** K₁/K₀ instead of
K₀/K₁, which is 2.5× too large at PIEZO1's r₀/λ. Both the formula and the
solver looked reasonable in isolation; the disagreement only surfaced when the
functional was integrated over the *exact* analytic profile — a third,
independent route to the same number. Two implementations agreeing is weak
evidence when they share an author; three disagreeing is what localises the
fault.

### One thing flagged rather than fixed

PIEZO1's dome meets the bilayer at a contact slope near 2.0 — about 63°. The
Monge gauge assumes |∇h| ≪ 1 and drops terms of order |∇h|², so at that slope
the neglected terms exceed the retained ones. The code solves it anyway,
because the trend and scale are still informative, but
`FootprintSolution.validity_note()` states plainly that the numbers are not
quantitative and that a nonlinear Helfrich or Euler–elastica treatment is what
the problem actually needs. Reporting a number with a caveat is better than
either silently reporting it or refusing to compute it.

---

## Round 4 — experimental conformational space (2026-08-06)

PCA over the deposited structures, compared with the elastic-network modes.
This is the strongest validation the project has produced.

**PC1 = 90.0% of variance, and it is the gating coordinate.** The PCA sees only
coordinates — no state labels — yet PC1 orders every structure correctly:
seven curved entries negative, the 8IXO intermediate at +334, flattened 7WLU at
+678, flat 11ZC at +1045. It overlaps **0.804** with ANM mode 6, cumulative
0.960 over 30 modes, RWSIP 0.555 against a random control of 0.001.

And the top three principal components all match **A**-symmetric modes, even
though E modes outnumber A two to one in the mode set. The symmetry selection
rule — only C3-symmetric modes can couple to isotropic tension — now shows up
in the deposited structural record, not merely in one pairwise transition.

### Four traps, all of which return a number rather than an error

1. **Species.** Human entries are numbered by Q92508 and mouse by E2JF22, with
   a non-constant offset. Everything converts to human numbering first.
2. **Coverage.** All 20 usable PIEZO1 entries share only 325 residues, because
   a couple of poorly-ordered structures drag the intersection down. Entries
   are dropped worst-first with the cost recorded.
3. **Protomer correspondence.** Four entries label their protomers in reversed
   rotational order. A test now asserts that at least one is detected — if that
   ever returns none, matching has silently broken.
4. **Paralogues.** 6KG7 is *PIEZO2*. Putting a 40%-identity paralogue into an
   ensemble meant to describe one protein's motion would be a category error,
   so it is excluded by default.

### The exclusion that changed the answer

6LQI is the Piezo1.1 splice isoform, missing residues 1382–1405. What sets it
apart from the rest is a **sequence** difference, not a conformational one, but
PCA cannot tell those apart — it sees only coordinate variance. Included, it
dominates an entire component by itself and splits the gating coordinate across
PC1 (58%) and PC2 (36%). Excluded, PC1 is a single clean 90%.

The general lesson: an ensemble method will happily report the largest axis of
variation without caring whether that variation is the biology you were asking
about. Deciding what belongs in the ensemble is part of the analysis, not a
preliminary to it, and every exclusion here carries its reason in the code.

---

## Round 5 — allostery and force transmission (2026-08-06)

Completes Block A. Perturbation response scanning, dynamic cross-correlation
and correlation-weighted allosteric pathways, all from the elastic network's
covariance — which is never formed in full, since for a PIEZO1 trimer it would
be an 11466² matrix of about a gigabyte.

**The anchor is the transmission hub.** Forcing the blade→gate path through it
costs a detour penalty of −0.000: it is already on the optimal route. It ranks
second by betweenness (5.19) behind only the CTD (7.67), and the cap is clearly
not a transmission route (+0.055).

**The beam result is softer than the prediction, and is reported that way.**
The lever model says the beam carries blade motion to the pore. It does not
appear on the single shortest path; but forcing the route through it costs only
+0.010, so it is a near-degenerate parallel channel rather than an excluded
one, with low but real betweenness (1.30). The honest statement is "viable
parallel route, not the dominant one" — not "confirms the lever model" and not
"refutes it".

Worth noting the beam is fully resolved (all 66 residues) in 8YEZ, so this is a
genuine negative rather than an artefact of missing density. That was checked
before drawing the conclusion.

### The error: a detour cheaper than the direct path

Asking whether the signal passes through a region by computing source→region
and region→target separately and adding their costs is wrong. Each leg
independently picks its best endpoints, and on a C3 trimer those can be in
*different protomers* — so the two legs never join into a path at all. Done
that way, routing "via the beam" came out at 0.101 against a direct path of
0.223: cheaper than the shortest path, which is impossible by definition.

The tell was the impossibility itself. A constrained optimum can never beat an
unconstrained one over the same feasible set, so a negative penalty is not a
surprising result, it is proof of a bug. `detour_cost()` now minimises
`d(source→v) + d(v→target)` over a shared via-point `v`, and an invariant test
asserts the constrained path is never cheaper.

Single shortest paths also turned out to be fragile — one marginally better
edge reroutes the whole thing — so `path_betweenness()` aggregates over many
source/target pairs, which is the standard dynamical-network-analysis answer
and far more stable.

### Block A review

The physics chain is closed and every link is validated against a published
number. The two results strong enough to build on: PC1 of the experimental
ensemble *is* the gating coordinate and matches an A-symmetric mode at 0.804,
and the anchor is the dominant force-transmission hub. Block E has been added
to the roadmap: allostery-derived per-variant features, licence-clean external
predictors via the ProtVar API, nonlinear membrane mechanics to fix the flagged
small-slope violation, pore wetting prediction, and — written *before* the
blind test rather than after — a statistical protocol under which "the
predictor does not separate GoF from LoF" is a recordable outcome rather than a
prompt to keep tuning.

---

## Rounds 6-7 — the blind test, and a null result (2026-08-06)

Round 6 built the predictor: ΔΔG_gating = ½dᵀ(H_mut − H_wt)d, the change in
elastic cost of the *observed* gating motion. Exact to 7e-16 against an
explicitly rebuilt Hessian, and cheap, because H_mut − H_wt is non-zero only at
the mutated residue's contacts.

**No phenotype comparison was made in Round 6, deliberately**, and
`docs/PREREGISTRATION.md` was written before Round 7 ran — pulling Round 20's
protocol forward, because a blind test is only blind if the rule was fixed
first.

### The result

**H0 not rejected.** 25 variants (16 GoF, 9 LoF) after pre-registered
inclusion. Permutation p = 0.234, Cliff's delta −0.083 with CI spanning zero,
AUROC 0.542. The mean difference points the predicted way; the effect is
negligible.

Writing the pre-registration turned out to matter more than expected. Two
numbers were close enough to tempt: the secondary analysis on the normalised
score came in at **p = 0.054**, and the primary difference had the right sign.
Without a decision rule fixed in advance, "the normalised version is basically
significant" is exactly the sentence one writes. The rule said p ≥ 0.05 is a
null result, and it is reported as one.

### Why it fails, and why that is the useful part

Post-hoc, labelled as such: **99.8% of the ΔΔG variance is between-position,
0.2% within-position.** The predictor reports *where a residue is*, not *what
happened to it*. That falls straight out of the construction — ΔΔG scales with
local gating strain and contact count, both properties of the position, while
the substitution enters through a single scalar spring multiplier, a far weaker
lever.

The R2456 series makes it concrete: four substitutions at one position, GoF and
LoF among them, all predicted to soften, with the *largest* softening belonging
to the loss-of-function variant R2456C. Positions do not have phenotypes;
variants do.

### What this does not mean

It does not touch the physics chain. Dome curvature 9.7 nm against a published
10.2; T₅₀ 4.99 mN/m against 5.1 ± 0.2; footprint decay 13.998 nm against 14.0;
PC1 of the experimental ensemble matching an A-symmetric mode at 0.804. Those
stand. What fails is a single scalar as a phenotype call.

The right reading is that the elastic network models the *machine* well and the
*substitution* badly. Round 17's sequence-based predictors (AlphaMissense, EVE,
ESM-1b via ProtVar, all substitution-aware by construction) are the natural
complement — they lack the mechanism, which is exactly what this supplies.
`docs/VALIDATION.md` §8 records what a fair next test would look like, written
now so that whatever is tried next is a stated new hypothesis rather than a
retrofit.

---

## Round 8 — pockets and ligands (2026-08-06)

Delaunay alpha spheres, the fpocket construction, reimplemented in numpy. Every
tetrahedron of the Delaunay triangulation has an empty circumsphere; its radius
says whether it sits in packed interior, in bulk solvent, or in a cavity of
small-molecule size.

**Two annotated sites recovered from geometry alone:** the transmembrane
hydrophobic gate (2/3 residues) and the anchor-domain apex brake (2/2). Neither
the detector nor the annotation knows about the other.

### The percolation trap

The first run produced a top "pocket" of **408 000 Å³ with 601 lining
residues** — which is the outside of the protein. On a large open structure a
radius filter alone is not enough: PIEZO1 is a curved propeller with enormous
grooves between its blades, and single-linkage clustering happily percolates
the whole exterior into one object.

Requiring each alpha sphere to have at least 30 atoms within 8 Å discards the
surface spheres and stops the merge; the largest pocket becomes 6 691 Å³ with
63 residues. **The parameters were chosen on pocket-size plausibility — a
druggable cavity is hundreds to a couple of thousand cubic Ångström — and
fixed before any site recovery was checked.** Tuning them until Yoda1 appeared
would have made the recovery meaningless.

### The Yoda1 result, and why the negative is the interesting part

The detector does *not* recover the Yoda1 site as an enclosed cavity: at most
one of its three residues, in either the human apo structure or the mouse
lipid-bound one. Allowing surface grooves recovers two.

Rather than report that as a failure, it is worth asking what it means, and
three independent facts line up. Yoda1 is proposed to act as a **molecular
wedge from the lipid phase**, which is a description of an interfacial site
rather than a pocket. A **PLX lipid occupies part of the site** in 7WLT — the
contact mapping here confirms it touches A2091. And the site has never been
seen in a co-structure: every PIEZO entry in the PDB contains only lipids, so
the mapping rests on mutagenesis and docking, which is exactly why this
project's annotation labels its evidence as *predicted*.

So the honest statement is that the Yoda1 site is interfacial, not enclosed,
and a cavity detector is the wrong instrument for it. That is a more useful
sentence than either "not recovered" or a tuned-up "recovered".

### Also fixed

`np.linalg.solve` in numpy 2 treats a 2-D right-hand side as a single matrix
rather than a batch of vectors, so the batched circumsphere solve needed an
explicit trailing axis. And the pocket tests originally took 137 s because each
one recomputed the pockets; module-scoped fixtures and a smaller Monte-Carlo
sample bring that to 16 s.
