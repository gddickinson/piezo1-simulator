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

---

## Round 9 — conservation and constraint (2026-08-06)

62 vertebrate orthologs from UniProt, deduplicated to one per species,
reference-anchored to human numbering.

**The result worth reporting is a convergence.** Ranked by mean conservation,
the **anchor domain is the most constrained region of PIEZO1 (0.987)**. Round 5
had already identified the anchor as the force-transmission hub — the only
region with a zero detour penalty on the blade-to-gate path. Those are entirely
independent lines of evidence: one is elastic mechanics on a single structure,
the other is 62 genomes' worth of selection. They agree.

The gradient across the whole protein is consistent with that picture: pore
module and anchor at 0.95–0.99, the cap and beam in the middle, and the distal
blade THU1 at 0.719 — least constrained, and also the region no experimental
structure resolves.

### Conservation alone is not a hypothesis

426 positions are invariant, carry no reported variant, and are structurally
resolved. That is about a quarter of the modelled protein, and as a "candidate
functional site" list it is useless — PIEZO1 is simply very conserved.

What makes it a hypothesis is the intersection with mechanics, which is
something this project can do and a sequence-only method cannot. Crossing
conservation with the Round 5 perturbation-response and path-betweenness
profiles narrows it sharply, and the top distal candidates are dominated by the
anchor (20 of 40). Two residues, **2021 and 2034**, are invariant across all 62
species, have never been reported as variants, and lie *on the blade-to-gate
allosteric path* — they appear literally in that path's residue list. Those are
well-motivated mutagenesis targets.

Also worth noting as a negative control that behaved: the **Yoda1 pocket is the
least conserved annotated site** (0.859, with A2075 at 0.63). A synthetic
agonist's binding site has no reason to be under selection, and this is the
third independent observation pointing at that site being unusual — the pocket
detector found it interfacial rather than enclosed, the PDB contains no
agonist co-structure, and now evolution is indifferent to it.

### The bug

The reference sequence was skipped by testing `member.sequence == ref`, which
discards *every* sequence identical to the reference rather than the reference
entry itself. Two closely related species with identical sequences would both
vanish, and on a small set that empties the alignment completely — the unit
test with three identical sequences returned conservation 0 everywhere. Fixed
to skip at most one exact match. Real-data numbers are unchanged, which is the
point: the bug was only reachable with near-identical inputs, which is exactly
what a synthetic test provides and real orthologs do not.

---

## Round 10 — research workflow (2026-08-06) — Block B complete

A headless CLI, provenance-stamped reports, session save/load and a documented
notebook API.

The design decision worth recording is that **the CLI and the report share one
analysis registry**. Adding an analysis to `ANALYSES` makes it available in
both at once, and a test asserts every registered analysis is reachable from
the command line. Two parallel dispatch tables would have drifted apart within
a round or two.

**Sessions store what you were looking at, never the data.** Structure, style,
camera, selection, which analyses had been run and with what parameters — but
no coordinates and no results. A session carrying its own copy of the numbers
would let a saved file drift silently out of step with the code that produced
them, which is precisely backwards for a reproducibility feature.

### What the batch run turned up

`python -m piezo1.cli batch` over all 20 structures reproduces the whole gating
series in one command: curved entries clustered at R_c 9.3–12.5 nm against a
published 10.2, the 8IXO intermediate at 16.5, and flat 11ZC at 21.6 as the
only entry called conductive.

It also **independently flagged 3JAC** — R_c 5.3 nm and spuriously conductive.
That is the same entry the Round 4 ensemble excluded for having poly-UNK
regions with arbitrary residue numbering, found again by a completely different
route. A batch mode earns its place partly by surfacing this kind of thing
without being asked.

### The argparse trap

`--json` was a top-level flag, so `cli dome 8YEZ --json` — which is how anyone
would actually type it — failed with "unrecognized arguments". Adding it to
each subparser via a shared parent fixed that but broke the other order:
a subparser writes its own default over whatever the parent already parsed, so
`cli --json dome 8YEZ` silently came back as `json=False`. `default=SUPPRESS`
plus a default supplied in `main()` makes both orders work, and a test pins
both.

### Block B review

The blind test came back null and its diagnostic was precise: the mechanical
predictor reports *where* a residue sits, not *which* substitution occurred.
Round 9 then showed conservation crossed with mechanics is sharp where either
alone is blunt, and Rounds 5 and 9 converged independently on the anchor.

The clearest gap now is that **the engine has raced ahead of the interface**.
The GUI can show a structure, a dome measurement and normal modes; it cannot
reach the pore profile, pockets, conservation, allostery or any reporting.
Block F leads with that, then a *new* pre-registration for a second variant
hypothesis (the Round 7 result stands as recorded and is not to be revised),
packaging, performance, and the teaching layer — which project aim A1 asks for
and which has had the least attention of anything.

---

## Round 16 — the feature table (2026-08-06)

`analysis/features.py` assembles everything the project computes into one row
per residue: PRS gate response and coupling, path betweenness, correlation to
the gate, amplitude along the lowest symmetric mode, fluctuation, relative
SASA, conservation, geometry and domain. 1279 residues × 11 features for 8YEZ
in 9 s.

**No phenotype comparison was run, deliberately.** The Round 7 blind test
returned a null result which stands as recorded; re-testing these features
against the variant labels needs the new pre-registration scheduled for Round
22. Assembling predictors and evaluating them in the same breath is exactly how
a blind test stops being blind, so the validation here is entirely structural:
gate response falls off with distance to the gate (r = −0.55), relative SASA
lies in [0, 1], the mode used as the gating coordinate is confirmed to be
A-symmetric, and the conservation join reproduces Round 9's domain ranking
exactly — which it would not if the residue join were off by a single position.

### The redundancy the tests caught

I wrote a test asserting that no two feature columns correlate above 0.99, on
the general principle that a near-perfect correlation means a duplicated
column. It fired immediately, and what it exposed is worth stating.

**The PRS response matrix is symmetric.** Its entries are Frobenius norms of
covariance blocks, and the Frobenius norm is invariant under transpose, so
‖C_ij‖ = ‖C_ji‖. Row and column means are therefore the *same numbers*:
"effectiveness" and "sensitivity" — effector and sensor — are not two
populations to compare, they are one quantity written twice.

Row normalisation appears to break the symmetry, and I initially "fixed" the
problem that way. It does not work either. Normalisation forces every row mean
to 1, so effectiveness becomes near-constant (spread/mean 0.0018, a range of
0.9655–0.9999), and the normalised *column* mean still correlates with the raw
row mean at **0.998**. Two columns, one quantity, wearing different scales.

Shipping both would have looked like two independent lines of mechanical
evidence and been one — the sort of thing that quietly inflates a combined
score. The table now carries a single, honestly named `prs_coupling`, and the
symmetry is documented on `PRSResult.is_symmetric` so the next person does not
re-derive it.

Highest remaining inter-feature correlation is 0.953 (`prs_coupling` against
`gating_amplitude`), which is high but genuinely two things: how strongly a
residue couples to the whole protein, and how far it moves along the specific
gating coordinate.

## Round 17 — external predictors, and a payload that lies by omission (2026-08-06)

The mechanical ΔΔG cannot see the substitution (Round 7). The obvious response
is to bring in predictors that can — AlphaMissense, EVE, ESM-1b, FoldX. The
obstacle was never the science, it was the licences.

### Why an API rather than local tools

FoldX is not redistributable at all: academic use needs a signed agreement, and
the community Python wrapper carries no licence file, which means all rights
reserved rather than permissive. SIFT4G is GPL-3.0, so calling it would drag
this whole PyQt application under copyleft. Everything on biosig.lab.uq.edu.au
— mCSM, DynaMut2 and relatives — publishes no licence text whatsoever, so
nothing is granted and using it would be a bet, not a permission. VarSite and
VarMap are both retired.

ProtVar (EMBL-EBI) serves all four predictors plus per-position conservation
from one endpoint, and I confirmed **CC BY 4.0** by reading `info.license` in
the service's own OpenAPI document rather than assuming it from the EBI's usual
terms. One source, one licence, no local models, attribution recorded on every
cached response and in `docs/REFERENCES.md`.

### The trap

`/score/{acc}/{pos}` returns nineteen entries per predictor — one per possible
substitution — and **no field anywhere in the payload says which entry is which
mutation**. Nothing errors. You get a well-formed JSON array of plausible
pathogenicity scores, and if you read them in array order, or assume they are
alphabetical by residue, you attribute the wrong score to every variant in the
study and never find out.

The fix is an undocumented `mt=` query parameter, which returns the single
score for the substitution asked for. A position-only query now keeps *only*
the conservation value, which is genuinely position-level; the missense scores
are dropped rather than guessed. `/prediction/foldx/` needed no such care — it
labels every entry with `mutatedType`, which is how the endpoint should behave
and how I noticed the other one did not. (The documented
`/prediction/interaction/` endpoint 404s; not needed here.)

This is the second time in this project that a silent field-alignment error has
been the real danger — the first was the mmCIF tokenizer shifting every column
by one because a trailing newline was not treated as whitespace. Both produce
confident numbers. Neither raises.

### What came back

64 of 65 single substitutions annotated in 77 s: conservation 64,
AlphaMissense 51, EVE 51, ESM-1b 51, FoldX ΔΔG 50. The 13 without missense
scores are nonsense and frameshift variants, where a missense predictor
correctly has nothing to say — absence there is right behaviour, and there is a
test asserting it stays that way.

### An unplanned validation worth more than the scores

ProtVar reports the wild-type residue it holds at each position, so annotating
the variants cross-checked our numbering against Q92508 from outside the
project. **Zero mismatches out of 64.** Rounds 1–2 established that the human
and mouse numbering differ by a non-constant offset across twelve blocks, and
that work has until now been checked only against itself. This is the first
independent confirmation that every curated variant sits on the residue we say
it does, and it is now a test.

### What these predictors still cannot do

All three missense predictors emit a single *pathogenicity* axis, benign to
damaging. That axis has no dimension in which to express **direction**. The
demonstration is R2456: all four substitutions score PATHOGENIC, but R2456H,
R2456K and R2456P are gain-of-function and R2456C is loss-of-function. A
predictor that calls all four damaging is not wrong — it is answering a
different question from the one this project is aiming at.

So the position after this round is that we hold two families of features that
fail in *opposite* directions: mechanical ones that see the position but not
the substitution, and sequence ones that see the substitution but not the
direction. Whether combining them recovers direction is a real hypothesis and
the only reason to run Round 22.

**No phenotype comparison was run.** Round 7's null result stands as recorded.
Touching the 68 labels again requires the new pre-registration to be written
first — that is the whole point of pre-registering, and the temptation to
"just peek" at 51 fresh AlphaMissense scores is exactly what the discipline is
for. Suite 243 → 254 passing, all of the new tests offline from cache.

## Round 18 — the nonlinear footprint, and a Round 3 result overturned (2026-08-06)

Round 3 built the linearised Helfrich footprint solver, got it right against the
exact Bessel solution, and then *flagged* that PIEZO1's contact slope of ~2.0
(63°) is far outside the regime where that theory holds. Having flagged it, it
quoted the numbers anyway: 622 nm² of footprint excess area against the dome's
256 nm², "about 2.4× as much deformable area as the dome". This round did the
calculation the caveat was standing in for, and the caveat was not strong
enough. The linear number is **3.5× too large**, and the conclusion reverses.

### The formulation

The Monge gauge writes the surface as a height field h(r) and expands in |∇h|.
There is no fixing that at 63°; the terms it drops are larger than the ones it
keeps. So parametrise the meridian by arc length instead, with ψ(s) the tangent
angle, which makes the principal curvatures exact with no expansion at all:
c₁ = ψ̇ and c₂ = sin ψ / r. Minimising the Helfrich energy subject to ṙ = cos ψ
with a multiplier η gives a first-order system in (r, z, ψ, M, η), solved as a
boundary-value problem.

The nice part is free. The Lagrangian has no explicit s dependence, so its
Hamiltonian is conserved, and that Hamiltonian *is* the axial force transmitted
through the membrane — zero for an inclusion nobody is pulling on. Imposing
H = 0 as a boundary condition and then measuring how far H drifts along the
solved profile gives an error estimate that costs nothing and depends on
nothing I derived. It runs at 7e-11.

### Checking a derivation I did by hand

A BVP solver will converge happily onto the wrong equations. So the checks that
count are the ones that do not reuse the derivation:

1. **Small-slope agreement.** The roadmap's own criterion. The relative
   discrepancy divided by slope² converges to a constant, 0.746 — not merely
   shrinking, but shrinking at exactly the order the Monge expansion discards.
2. **The exact functional in another gauge.** Re-evaluate the energy of the
   solved shape in Monge form with the unexpanded expressions
   (dA = √(1+h'²)·2πr·dr, c₁+c₂ = h''/(1+h'²)^{3/2} + h'/(r√(1+h'²))). Agrees
   to 1.3e-3, which is finite-difference noise.
3. **No nearby shape is cheaper.** Perturb the profile with bumps that respect
   both boundary conditions; the energy may only rise.

Check 2 is where I nearly went wrong in an instructive way. My first version of
it minimised the exact functional from scratch, starting from a linear guess,
and it *disagreed with the solver by 142%* at slope 2. That looks like a
refutation. It was not: the minimiser returned a **higher** energy than the BVP,
and a minimiser that has not converged always errs high. Starting it from the
BVP solution, it could not improve on it in 151 iterations. The independent
check was the weaker instrument, and the tell was the sign of the disagreement.

### The result

At the measured 7WLT geometry (inclusion radius 8.69 nm, contact slope 1.99):

| | linear | nonlinear |
|---|---|---|
| footprint energy | 92.2 k_BT | **25.3 k_BT** |
| footprint excess area | 622 nm² | **179 nm²** |

Invariant to domain truncation from 8λ to 40λ at six significant figures, to
grid, to solver tolerance and to whether the slope is walked up by continuation
or hit directly. The correction factor is 3.46–3.67× across κ = 20–25 k_BT and
γ = 0.42–3.0 mN/m, so it is not an artefact of one parameter choice.

### What was actually wrong with Round 3

Two things, and the second is the one worth remembering.

The first is the obvious one: linear theory used outside its range.

The second is that the comparison was never like for like. The dome's 256 nm² is
an **exact** area difference, measured from the fitted spherical cap. The
footprint's 622 nm² was **linearised**. Putting them in a ratio compares a
quantity to an approximation of a different quantity. Measured consistently the
footprint holds 179 nm², which is **0.70× the dome — less than the dome, not
2.4× more**.

I have left the Round 3 entry in ROADMAP.md struck through rather than edited
away, with a pointer to this round. The caveat was recorded honestly at the
time; the lesson is that recording a caveat is not a substitute for doing the
calculation, and deleting the evidence of that would remove the lesson.

### What this does *not* say

It does not refute Haselwandter & MacKinnon. Their argument concerns the
footprint's contribution to tension *sensitivity* — the area released between
closed and open states — and absolute stored area was never a test of it. Round
3's error was rhetorical as much as numerical: it presented our absolute-area
ratio as "the quantitative form of" their claim, which it is not. The defensible
statement is narrower: at PIEZO1's contact slope the linearised footprint is not
quantitatively usable, and corrected, dome and footprint hold comparable excess
areas with the dome slightly larger.

Suite 254 → 272 passing; GUI smoke test clean.

## Round 19 — hydrophobic gating, and what a heuristic is actually for (2026-08-06)

Round 1 built the pore profiler and it worked: 8YEZ's bottleneck is 0.095 nm,
11ZC's is 0.330 nm, and the profiler rediscovered the curated gate and CTD
constrictions from coordinates alone. But radius is a weak predictor of
conduction — Rao et al. put its AUROC at **0.59**, barely better than a coin —
because a pore can be wide enough for a hydrated ion and still block when a
hydrophobic neck expels liquid water. Their heuristic, combining radius with
local hydrophobicity, reaches **0.91**. This round implements it.

### Getting the boundary rather than redrawing it

The paper gives the construction, the 1 RT = 2.6 kJ/mol contour and the
Σd > 0.55 cutoff, but not the classification line itself — that lives in a
figure. Digitising a figure by eye would have been precisely the kind of silent
correctness bug Round 17 was about.

It turned out not to be necessary. CHAP is **MIT licensed** and its repository
ships `heuristic_grid.json`: the actual 100×100 water free-energy landscape over
(hydrophobicity, radius), built from ~600 MD simulations. It also ships the
exact normalised Wimley–White scale the landscape is indexed by, and the default
kernel bandwidth (0.35 nm). So the published artefact is used directly, as a
download rather than a commit, and analyses degrade to "unavailable" without it.

An independent check that we read it correctly: our extracted 1 RT contour gives
a critical radius rising from **0.10 nm** at the hydrophilic end to **0.43 nm**
at the hydrophobic end, against the paper's prose "hydrophilic pores wet below
0.2 nm; hydrophobic ones can hold a barrier out to ~0.4 nm". Nothing in our code
was fitted to those numbers.

### The bug that returned confident nonsense

My first hydrophobicity profile averaged residues in a sphere around each probe
centre. It ran, produced smooth plausible values, and gave 8YEZ a score of 0.45
— just under the cutoff, so the closed structure was called **conductive**.

The tell was the range. CHAP smooths **along the pore coordinate** over
pore-facing residues, and the published grid spans −0.45 to +0.30. My 3-D
neighbourhood was 1.85 nm wide, which pulls in the entire shell of residues
surrounding the lumen, and the profile collapsed into a band from −0.12 to
+0.02. Every energy was then read out of the landscape at a coordinate the
landscape was never built on. Nothing errored; the numbers merely meant nothing.

Rewritten as a proper Nadaraya–Watson average along the axis, using side-chain
centroids because it is the side chain that faces the lumen, the range opens to
−0.635…+0.229 and 8YEZ scores 0.82 — non-conductive. There is now a test
asserting the profile uses most of the grid's range, because "plausible but
compressed" is what this failure looks like from the outside.

### Right answer, right reason

The roadmap asked for the right answer *for the right reason*, which is a
demand for a control rather than an assertion. So: hold every radius fixed and
replace the hydrophobicity scale with a uniform hydrophilic value. If the closed
verdict were a radius threshold wearing a disguise, the score would not move.
It goes from **0.82 to 0.00** — conductive. The verdict is chemistry.

The sharpest single fact: 8YEZ's F2451 and V2454 sit at **0.325 nm** and are
called dewetted; 11ZC's *bottleneck* sits at **0.330 nm** and is called wet.
Same radius, opposite verdict. And the flagged set — F2451, V2454, R2467,
F2468 — is the curated hydrophobic gate and cytoplasmic constrictions, which the
heuristic never sees.

### The limitation I would have missed by testing only what was asked

The round specified two structures. Running five exposed something the two
would not have: **7WLU and 8IXO have 0.098 nm bottlenecks — less than a water
molecule's 0.15 nm radius — yet score 0.11 and 0.30, i.e. open.**

That is not a bug in the heuristic; it is what the heuristic is *for*. It
answers "would water dewet here?", not "does water fit here?". Rao et al. built
it to find hydrophobic gates, which are by definition blockages *without* steric
occlusion. A pore too narrow for water never poses the wetting question.

Merging the two into one verdict would have hidden this, and would have made the
heuristic look like a general conduction predictor it does not claim to be. So
`WettingPrediction` exposes `hydrophobic_gate` and `sterically_occluded`
separately and `conductive` requires neither. With both, all five states come
out right: 8YEZ and 7WLT shut on both counts, 7WLU and 8IXO shut on sterics
alone, 11ZC open.

Suite 272 → 288 passing; GUI smoke test clean.

## Round 20 — what the null was entitled to claim (2026-08-06)

Most of this round's list was already done. Rounds 6–7 delivered the permutation
test, the bootstrap effect size and the pre-registered decision rule; Round 7
recorded the null. What was missing were the two questions asked *around* a
result rather than by it — could the design have detected the effect, and did we
look in too many places — plus cross-validation. Those went into a new
`analysis/design.py`, kept separate from `validation.py` because "did it work?"
and "could it have worked?" are different questions and conflating them is how
underpowered nulls get overclaimed.

### The finding

Simulating the pre-registered test at Round 7's actual group sizes, 16
gain-of-function against 9 loss-of-function, one-sided at α = 0.05:

| effect | Cliff's δ | power |
|---|---|---|
| **observed in Round 7** | −0.083 | **0.13** |
| small | −0.11 | 0.16 |
| medium | −0.28 | 0.35 |
| large | −0.43 | 0.60 |

**80% power is reached only at |δ| ≥ 0.55**, past 'large' on the usual
thresholds. I ran it twice — once from a normal model, once resampling the
observed heavy-tailed ΔΔG values, because a difference-in-means test loses power
on heavy tails and the normal model would flatter the design. They agree: 0.55
and 0.56.

So the Round 7 null **excludes a large mechanical effect and is close to
uninformative about a small or medium one.** That is a real qualification and it
was not stated at the time. I have added it as `VALIDATION.md` §6b, explicitly
marked as not amending §§1–3: the result stands, its scope is now bounded. The
existing §6 diagnostic — 99.8% of ΔΔG variance is between-position — remains the
mechanistic explanation and is independent of this. Both are true. The power
limit is simply the one nobody wrote down.

### The constraint that binds the rest of the project

At 80% power with equal groups: **42 variants for a large effect, 98 for a
medium one, 600+ for a small one.** Twenty-five survive Round 7's inclusion
criteria, and relaxing every criterion cannot reach forty-five.

So a *confirmatory* test of anything below a large effect is not available from
this variant set, however good the predictor gets. Round 22 is now explicitly
bound to declare in advance whether it is confirmatory-for-a-large-effect or
exploratory, and I have added Round 27 — expand the phenotyped set — because the
binding constraint on this project's central claim turns out to be **data, not
method**. That is worth knowing at round 20 rather than round 40.

### Why the protocol had to be written now

`docs/NEGATIVE_RESULT_PROTOCOL.md` is deliberately written *between* Round 7's
null and Round 22's hypothesis. Written after Round 22 it would be a rule fitted
to a result already seen, which is the thing it exists to prevent.

The failure mode it guards against is not dishonesty, it is drift: a null comes
back, and there is always one more defensible adjustment — drop the noisy
variants, use the normalised score, go two-sided. Each is reasonable alone.
Together they are an unrecorded search, and the p-value at the end means nothing
like what it appears to. With 68 curated variants and one primary claim, that
resource is finite and non-renewable.

The multiplicity worked example makes the point better than the argument does.
Six candidate predictors with plausible p-values: three clear 0.05, and after
Benjamini–Hochberg **none survives**. Reporting AlphaMissense at p = 0.012
without its family would be a false discovery manufactured by looking six times.

### A sign error caught by its own diagnostic

`power_curve` injected the effect into the wrong group. `shift_for_delta` is
defined with a ~ N(0,1) and b ~ N(shift,1), so the displacement belongs on group
b; I added it to a. The resulting power curve rose monotonically with effect
size, saturated at 1.0, and was entirely plausible — it was simply the power to
detect the opposite direction.

It was caught because `power_curve` measures the Cliff's delta it actually
achieved rather than assuming it got what it asked for. That diagnostic existed
only because Round 18 had just taught the same lesson: the useful check is the
one that does not reuse the thing being checked. It is now a test.

### Cross-validation

Leave-one-out on the Round 7 predictors gives AUROC 0.535 out-of-sample against
0.542 in-sample, optimism +0.007. Small because the default combination is
unsupervised, and reassuring: there was no hidden overfitting inflating the
original number. The machinery matters more for Round 22, where a *fitted*
combination of six predictors on 25 variants would otherwise measure how well 25
points can be fitted.

Suite 288 → 311 passing; GUI smoke test clean.

## Round 21 — the engine reaches the interface (2026-08-06)

Twenty rounds built an engine the GUI could not reach. It could show a
structure, measure a dome and animate normal modes; the pore profiler, the
pockets, conservation, allostery and every report were CLI-only. This round
closes that gap.

### The thing worth guarding against

A GUI that recomputes what the engine already computes becomes a second
implementation, and second implementations diverge. So the workers call the
same functions the CLI calls, and a test asserts the pore worker reproduces
`pore_profile` and `predict_wetting` *exactly* rather than producing something
plausible. The smoke test then reads the panel's own label back: bottleneck
0.95 Å, Rao score 0.82, "non-conductive (sterically occluded + hydrophobic
gate)" — identical to `piezo1 hydration 8YEZ`.

### Two axes, because that is the whole point

Round 19's finding is that radius alone predicts conduction at AUROC 0.59 and
radius-with-hydrophobicity at 0.91. A plot that puts both on one y-scale makes
the hydrophobicity trace a flat line at the bottom of a 0–10 Å axis, which
destroys exactly the comparison the plot exists to show. So the widget has
independent left and right axes, with the left anchored at zero (a radius of
zero is meaningful) and the right not (hydrophobicity is signed, and anchoring
it at zero squashes the range carrying the signal).

I wrote it as a QPainter widget rather than adding matplotlib or pyqtgraph.
One plot type, has to repaint inside a dock at interactive rates, has to match
the dark theme, and both libraries drag in their own event-loop integration —
more to configure than to write.

### Where the colour schemes went, and why not where the roadmap said

The roadmap asked for conservation and PRS as "first-class colour schemes",
which I read as the main colour dropdown next to Domain and Chain. I put them
in the Analysis dock instead, appearing once computed.

The reason is that they are not the same kind of thing. Domain and chain are
properties of a loaded structure and are available the instant it loads.
Conservation needs 61 orthologs fetched and aligned; PRS needs normal modes
computed first. An entry in the main dropdown that silently does nothing until
some other panel has been used is a dead control, and a user who selects it and
sees no change has been told the feature is broken.

They do go through the existing `ColorBy.VALUE` path rather than a new one —
the same mechanism mode-displacement colouring already used.

### Two traps

**Unmeasured residues take the floor, not zero.** Conservation runs about
0.6–1.0. Giving an unmeasured residue a zero stretches the colour scale across
a range containing no data and washes out every real difference — while looking
like a perfectly valid map. The value map now fills with the minimum of what it
was actually handed.

**Restoring a session rebuilt the scene four times.** Style, colouring, ligands
and atom size each emit a signal. Setting all four from a saved session fires
four rebuilds of a 120k-atom trimer, which is visible as flicker through three
intermediate views nobody asked for. `set_state` blocks signals, sets the
widgets, then emits once.

Also worth recording: conservation drops positions with ortholog coverage below
0.7. Those values measure how well the alignment covered the position, not
selection pressure on it, and left in they paint a confident band across the
unresolved distal blade. 2484 residues survive, mean 0.770, 901 above 0.95.

### Testing a GUI

Qt refactors have silently broken this application twice, so the tests run real
widgets on the offscreen platform rather than mocks. Mocks would not have caught
either breakage — both were real widgets behaving differently from how the code
assumed, which is precisely what a mock encodes rather than tests.
`scripts/screenshot_app.py --analysis` now drives the new panel end to end and
round-trips a session through the controller.

Suite 311 → 330 passing.

## Interlude — the interface, driven by using it (2026-08-06)

Round 21 put the engine behind panels; actually running the application turned
up four things that no amount of testing from inside had shown. All four came
from the user opening it on their own machine, which is the only way some of
these surface.

### The window did not fit the screen

`resize(1680, 1000)` was hard-coded. On a laptop or a 1080p display that puts
the title bar above the top of the screen, and on some window managers the
window then cannot be dragged back or resized. It now opens at the smaller of
the preferred size and 95% of the *available* geometry — which already excludes
the menu bar and dock — and centres itself.

That was necessary but not sufficient. Qt will not shrink a dock below its
content's minimum, so the tallest panel still set a floor on the whole window.
Every panel is now wrapped in a `QScrollArea`, and the window's minimum size
hint dropped to 535×258. There is also `--geometry WxH` for when screen
detection misjudges a multi-monitor setup.

### The model jumped whenever anything was selected

Selecting a domain, site or variant called `_focus_residues`, which moved the
camera pivot to the selection's centroid. On a 2500-residue trimer that means
the whole structure lurches every time you click a list entry.

The honest fix was not to remove it — on a structure this size a highlight can
easily be behind the protein, and centring is genuinely useful — but to make it
a choice. Options → *When something is selected* now offers **keep the view
still** (the new default), **centre on the selection**, or **centre and zoom**.
The zoom variant deliberately preserves the current *orientation*: reframing
rotation as well would throw away the view the user had set up, which is the
same complaint arriving by another route.

### Panels could not be rearranged

Docks were created with Qt's defaults, which confine each one to the two areas
it was born in and give no way to float or close it. All five are now movable to
any edge, floatable into their own windows, closable, and tabbable.

Making them closable creates the obvious problem: a user who closes two panels
and floats two more has no route back. So the shipped arrangement is captured
once at startup — *before* any remembered layout is applied, so Reset is always
the application's own layout rather than whatever was last left — and
View → Reset layout (Ctrl+R) restores it. The trap here is `objectName`:
`saveState`/`restoreState` key on it, and without one Qt silently declines to
restore that dock, which looks like a corrupt settings file rather than a
missing property.

Layouts persist between runs, with the geometry clamped to the current screen on
restore — otherwise a layout saved on a large external monitor reopens mostly
off-screen, which is the first bug again by a different path.

### Menus, help and tooltips

Options and Help menus now exist. The help window is non-modal, because the
point of a feature guide is to be read *while* driving the application.

The guide has seven topics, and the one that matters most is **"What this
application will not do"** — the pre-registered blind test that returned
p = 0.234 and AUROC 0.542, the diagnostic that 99.8% of the predictor's variance
is between-position, the power limit that lets the null exclude only a large
effect, and the footprint number that turned out 3.5× too large. A test asserts
those figures are present, because a help file is exactly where an inconvenient
result would quietly stop being mentioned.

Tooltips now carry provenance rather than restating labels: the dome button says
what number it should produce and against which publication, the hydrophobicity
checkbox gives 0.59 against 0.91, the mode selector explains why only A modes
can couple to isotropic tension.

Suite 330 → 342 passing, and the GUI smoke test now checks the two behaviours
that were reported — that focus-off leaves the camera still and focus-centre
moves it, and that reset restores hidden and floated panels.

## Round 22 — the second null, and why it was foreseen (2026-08-06)

Round 7 tested a mechanical predictor and returned a null whose diagnostic was
precise: 99.8% of its variance was between-position, so it reported *where a
residue sits* rather than *which substitution occurred*. Round 17 brought in
four substitution-aware predictors. This round asked whether they do better.

The pre-registration was written and **committed in its own commit before the
test ran**. That is not ceremony. A rule written after seeing a p-value is a
rule fitted to it, and with 68 curated variants and one primary claim, the
resource being spent is finite.

### The design was worse than Round 7, and I knew that before testing

Of 39 variants carrying a directional label, only 26 are single-residue
substitutions. **Eleven of the thirteen dropped are loss-of-function** — stop
codons, frameshifts, deletions.

That is not a curation artefact. It is how loss of function actually happens:
you break the protein. And no missense predictor can score a stop codon.

So the usable design is 20 GoF against 6 LoF, and 80% power needs |δ| ≥ 0.61
where Round 7 needed 0.55. Power at a conventionally large effect is 0.52 — a
coin toss. Under the Round 20 protocol that forces a declaration, and I declared
the round **exploratory**: at this n a "confirmatory null" would exclude only
effects beyond large, and calling that confirmation of anything would overstate
it.

### The hypothesis, and the objection recorded alongside it

Primary endpoint: FoldX ΔΔG. The mechanism is that loss of function can be
achieved by breaking the protein, while gain of function cannot — a channel that
opens too readily must still fold, traffic and gate, so gain-of-function
variants are constrained to be structurally tolerable. Prediction: LoF more
destabilising.

I also recorded, in the same document, the objection: excluding the truncating
LoF variants removes precisely the "break the protein" mechanism, leaving a LoF
subset selected for **not** being truncating. Writing that down first is what
stops it later being a post-hoc rescue.

### The result

**Cliff's δ = −0.211, CI [−0.684, +0.298], AUROC 0.395.** Mean ΔΔG 0.767 for
LoF against 1.309 for GoF — the point estimate runs *opposite* to the
hypothesis. Mean difference, Cliff's delta and AUROC all agree on the direction,
so it is not a sign error. The interval spans zero, so it is uninformative
rather than a reversal.

The secondary family: nothing significant, nothing close, smallest adjusted
q = 0.448. The pre-registered expectation that AlphaMissense, EVE and ESM-1b
would show *no* separation held — a single pathogenicity axis has no room for
direction when both classes are pathogenic. Combining all five by equal weights
gives leave-one-out AUROC 0.535 with zero optimism: nothing recovered.

### What the foreseen objection now looks like

It looks like the explanation. A destabilisation predictor was asked about the
one subset of loss-of-function variants that does not act by destabilisation,
because the subset that does was excluded for being unscoreable. That it was
written down before the run is the only reason it can be offered at all.

### Where this leaves the project

Two nulls, two different predictors, both pre-registered, neither revised. The
central claim is not refuted — it is **untested at adequate power**, and after
two attempts the constraint is clearly data rather than method. Fifty-three
variants would be needed for a large effect at this class ratio, 148 for a
medium one, against 26 available. Round 27 is now the round that matters.

There is also a sharper question hiding in the counts. Eleven of seventeen
loss-of-function variants are truncating, and predicting *those* from structure
is close to trivial. The hard and interesting question — distinguishing "opens
too easily" from "does not open" among missense variants — has a sample size of
six. That is the honest scope of what this project can currently ask.

Suite 342 → 352 passing.

## Interlude II — measurement furniture, sequences and overlays (2026-08-06)

Four more requests from actually using the application. Each turned out to have
a correctness question hiding inside what sounded like a display feature.

### Scale bars and a clock

A structure screenshot without a scale bar states nothing quantitative; a reader
cannot recover the scale of a perspective projection by eye. The bar rounds to a
tidy length — 1, 2, 5, 10, 20 … Å — because a bar labelled 47 Å is unreadable.
It is exact in the plane through the camera pivot, which is where the molecule
is, and the docstring says so rather than implying a perspective bar is exact
everywhere.

The animation clock raised the more interesting question. For a **mode sweep**
it reports elapsed seconds and the fraction of a display cycle, explicitly
labelled *not a physical period*: the frequency of an elastic-network mode is
not time. For a **morph** it refuses to report seconds at all and gives percent
along the path, because a morph is an interpolation between two endpoints and a
seconds axis would imply kinetics the model does not contain.

### Presentation mode

Full screen hides the panels and menu. The detail worth recording is that
leaving it restores each panel's **previous** visibility rather than showing
them all — someone who had closed a panel before presenting should not find it
reopened afterwards. And Escape leaves, because with the menu bar hidden there
is otherwise no visible way out.

### The sequence window

Three kinds of sequence are involved and they are **not** interchangeable:
UniProt (the reference numbering variants are quoted in), structure (only
resolved residues, with gaps, starting at 570 for 8YEZ), and translated CDS.
`NamedSequence` carries its own positions for exactly this reason — a viewer
that showed "the sequence" without saying which one is a numbering bug waiting
to happen, and this project has already been bitten by three of those.

For DNA I fetched the real Ensembl transcript rather than back-translating. A
back-translation would look like a gene, invent every codon choice, and make
silent variants unrepresentable. The check that the right transcript was
fetched: **the human CDS translates to Q92508 at 100% identity over 2521
residues.** Mouse gives 99.9% — three differences at 147, 229 and 1572, being
the reference genome against the TrEMBL submission. That is recorded in a test
rather than smoothed over.

I also guessed the mouse transcript ID from memory and it 400'd, which is the
same failure mode as the six wrong PMIDs in Round 8. Looked it up properly.

Comparison offers global alignment or pairing by residue number, and the second
is not a lesser option: when two sequences already share a numbering, an aligner
can slide a run of residues to buy score and manufacture differences that are
alignment artefacts. Pairing by number is refused across species, where it would
be meaningless.

### The structure overlay

The comparison this project needs constantly — curved against flattened, wild
type against variant — is two structures in one frame.

The correctness question is protomer correspondence. Round 4 found four
deposited entries labelled in the reverse rotational order. Overlaying 7WLU on
7WLT, correspondence search rematches to (2, 1, 0) and gives 12.3 Å; taken at
chain-label face value the two sit **90.7 Å** apart. A viewer that trusted
labels would display an enormous conformational change that does not exist —
and it is exactly the conclusion someone would want to draw. Both numbers are
reported.

A second guard: cross-species overlays are refused outright. Residue numbers are
the join key, and human and mouse numbering do not correspond, so the fit would
be confidently wrong rather than obviously wrong.

Building the protomer blocks per structure also failed, because two entries
rarely resolve the same residues and `match_protomers` cannot compare blocks of
different length. The basis has to be shared across both.

Suite 352 → 376 passing.

## Round 23 — one command, and a guard against the docs going stale (2026-08-06)

Aim A5 says a fresh clone plus an environment plus a fetch should reproduce the
whole working state. The packaging half of that is routine: `pyproject.toml`,
two pinned lock files, a `Makefile` whose targets each wrap `conda run` so they
work from a bare shell. One decision worth recording — **the GUI dependencies
are an optional extra, not a requirement**. Everything below `render` runs
headless, and that is precisely what makes the science testable without a
display; making PyQt mandatory would quietly give that up.

### The half that actually matters

This project states a great many specific numbers in prose. A dome radius of
9.7 nm. A mode overlap of 0.705. A half-activation tension of 2.71 mN/m. A
footprint of 179 nm² that replaced one of 622. Two null results with their
p-values and effect sizes.

**Prose does not fail a test suite.** A solver rewrite, a re-fetched structure
or a changed default can leave `docs/SCIENCE.md` confidently asserting a value
the code no longer produces, and nothing in the repository would notice. Given
how many recorded numbers have already needed correcting here — the footprint
was wrong by 3.5×, the Bessel ratio by 2.5×, the decay length by a factor of
three — the risk is not hypothetical.

So `piezo1/analysis/claims.py` holds seventeen claims. Each names the documented
value, its tolerance, the document it appears in, and a callable that recomputes
it from scratch. `make verify` runs them in about ten seconds. All seventeen
currently reproduce.

A claim is not a test, and the distinction is the point. Tests assert the code
behaves; claims assert **the documentation is still true of the code**. They
fail for different reasons and both are worth having.

### Two design decisions inside that

**Frozen claims.** The four recorded validation numbers — Round 7's p-value and
AUROC, Round 22's effect size and AUROC — are marked `frozen`. If one drifts,
the report prints an explicit instruction *not* to edit the document to match,
but to work out why the computation changed. The obvious way to resolve a
failing claim is to update the prose, and for a recorded null result that would
be exactly the wrong move.

**Skipped is not drifted.** A claim that cannot run because a structure has not
been downloaded reports as skipped. Conflating that with drift would make every
fresh clone look like it had a broken scientific record, and the person seeing
it would learn to ignore the report.

### Testing the detector

A drift detector that has never detected drift is decoration. There is a test
that feeds it a deliberately wrong claim and requires it to complain, and
another that feeds it a claim raising `FileNotFoundError` and requires *skipped*
rather than *drifted*. Both would have passed vacuously if the registry simply
returned success.

Suite 376 → 390 passing.

## Parameter registry — every number gets a citation (2026-08-06)

An audit of `physics/`, `structure/` and `analysis/` found **203 numeric
literals**. Some were physical constants from named papers; some were
convergence tolerances; and from the outside there was no way to tell which was
which. That is the problem: a constant written into a function default is
invisible. You cannot list them, show them to a user, or trace one to a paper.

This project has already had to correct several numbers that were invisible in
exactly that way — the footprint area wrong by 3.5×, the Bessel ratio by 2.5×,
the biharmonic decay length by a factor of three.

### What was built

`piezo1/resources/parameters.json` holds **61 parameters** across 12
categories, each with a unit, bounds, a kind (physical / empirical / method /
convention), a description and a citation. It is authored by
`scripts/build_parameters.py`, in the same shape as `build_variants.py`:
authored content, validated on the way out, committed as a resource. Keeping it
as data means the whole parameter set can be read and diffed without opening a
module.

**The gate is the point.** A citation must resolve to a key in
`references.json`, or be one of five sentinels — `method_choice`,
`measured_here`, `derived`, `convention`, `unverified` — each of which *obliges
the entry to say why* in `source_note`. The build refuses to write otherwise. It
caught eight entries where I had left the note blank, and I filled them in
rather than weakening the gate. Two references had to be added first (fpocket
and Shrake–Rupley); the reference builder's own `expect` gate then rejected my
first attempt at Shrake–Rupley because I picked a word that is not in the title.

31 of the 61 cite a published paper. The other 30 are method choices, and saying
so explicitly is more honest than dressing them up.

### Overrides are tracked, not silent

Modules resolve through the registry **at call time**, so a change takes effect
on the next call. That makes the values genuinely editable — and immediately
raises the real problem: a number computed with a changed parameter is not the
number in `docs/SCIENCE.md`, and nothing would have said so.

So three things enforce it:

- `verify_claims` **refuses to run** against a modified registry. Every
  documented number was produced at the defaults; recomputing with a changed
  value would report drift the *user* caused, and the obvious reading of that
  report is that the code is broken.
- Reports carry a warning banner **at the top**, not in the provenance footer,
  and `Provenance` gained a `parameter_overrides` field.
- The GUI keeps a persistent amber indicator in the status bar.

### The audit is what makes it a rule

A convention nobody checks decays at the first hurry. `piezo1/parameter_audit.py`
scans the three scientific packages and fails on any numeric literal that is
neither registered nor listed in `EXEMPT` **with a stated reason**. It started
at 18 unaccounted; 11 were registered-but-unmigrated and 7 were genuinely
unregistered — those became new parameters rather than new exemptions.

The exemption categories matter as much as the registrations. A convergence
tolerance, an iteration cap, a random seed and a zero-initialised dataclass
field are implementation details, and pretending otherwise would bury the real
parameters in noise. What the audit insists on is that the exemption is a
decision someone made and wrote down.

There is a test that invents a module containing `binding_energy = 7.25` and
requires the audit to catch it — a detector that has never detected is
decoration.

`make audit` and `make params` are build steps, and both run inside
`make reproduce`. The rule is written into `CLAUDE.md`. Suite 390 → 408 passing,
and all 17 documented numbers still reproduce, which is the check that the
migration changed nothing scientific.

## Round 24 — performance, and profiling before optimising (2026-08-06)

The roadmap named three slow paths. Profiling found that two of them were not
the problem and the real one was not on the list.

- **PRS**: 0.52 s. Never a problem.
- **Pocket detection**: 4.2 s, not the ~10 s the roadmap assumed.
- **Ensemble PCA**: 12.4 s — but **99% of that was mmCIF parsing**, not the
  PCA at all.
- **SASA**: 7.5 s, and the roadmap did not mention it.

Recording that is worth more than the speedups. A round spent optimising PRS
because a note said it was slow would have been a round wasted.

### The rule for this round

An optimisation that alters a number is a bug, not a speedup. Every change here
is a reformulation with the same value, and the tests assert **identity**, not
closeness.

### What was actually slow

**SASA, 7.54 → 1.27 s.** The inner loop built a (256 × neighbours × 3) array per
atom and took a square root. Two observations remove both: `d ≥ r` and `d² ≥ r²`
decide the same way for non-negative values, so the root was never needed; and
expanding `|t − x_j|² = |v_j|² + r_i² + 2r_i(p_k · v_j)` turns the 3-D broadcast
into a single BLAS product. Bit-identical across all 31,599 atoms —
`np.array_equal`, not `allclose`.

**The mmCIF tokenizer, and through it the ensemble, 12.35 → 2.05 s.** The
tokenizer walked characters in a Python loop, 321,913 times per ensemble load.
But 99.5% of mmCIF lines contain no quote and no comment, and for those
`str.split()` is exactly equivalent — its default whitespace set is the same one
the function uses, and it discards empty fields the same way.

I checked that claim rather than asserting it: 245,528 lines of deposited
structure, **zero mismatches**. This is the function whose whitespace handling
once shifted every column by one and produced `invalid literal for int(): 'ATOM'`,
so the careful path is left exactly as it was and merely bypassed.

**Conservation, 3.67 s → 0.003 s.** Aligning 61 orthologs dominated anything
needing conservation, and the result depends only on the sequences. The cache is
keyed on a **content hash** of the reference and the ortholog sequences rather
than on a filename or timestamp, so re-fetching or changing the reference
invalidates it automatically. A cache that can go stale is worse than no cache:
it would report last week's conservation against this week's alignment, and
nothing would look wrong.

### The check that matters

Overall 33.8 → 12.5 s, and every number unchanged: SASA total 197490.5582 Å²,
top pocket 6593.6 Å³, ensemble PC1 0.9000, pore bottleneck 0.9518 Å. All 17
documented numbers still reproduce, which is the guard Round 23 built for
exactly this kind of change. The suite itself went 118 → 97 s.

Timing assertions in the tests are loose ceilings only. Pinning a runtime would
fail on a slower machine for no scientific reason; what is worth pinning is that
the fast path and the careful path agree.

## Round 25 — the teaching layer, with nothing narrated (2026-08-06)

Aim A1 says this should be a learning instrument, and of the six aims it had
received the least attention. The tour walks the mechanism in eleven steps:
trimer, blades, dome, footprint, lever, gate, open state, normal modes, gating
energetics, a variant, and what the project cannot do.

### The rule that shapes it

**Every number a step states is computed when the step runs.**

A tour is prose, and prose is where numbers go to rot. Writing "the dome radius
is 9.7 nm" into a tour step would create a fourth place for that value to live —
beside the code, `docs/SCIENCE.md` and the claims registry — and the fourth copy
is the one nobody remembers to update. So each step carries a callable that
reads whatever the application has actually computed, and where a step quotes a
published comparison it reads it from the **parameter registry**. There is a
test that changing `dome.published_radius_closed` changes the tour text, which
is the proof it is not a literal.

The controller calls the same controllers the panels use. A teaching tool that
quietly disagreed with the application it was teaching would be worse than none.

### Where it ends

On the failures. The last step states both null results — p = 0.234, AUROC
0.542, and the second test's −0.211 with an interval spanning zero — and the
power limit that means the first excludes a large effect and little else. There
is a test asserting those figures are present, because a tour is exactly the
sort of document where an inconvenient result quietly stops being mentioned.

A learning instrument that only shows its successes teaches the wrong lesson,
and the lesson worth teaching here is that the interesting question about R2456
is why four substitutions at one residue do not all do the same thing.

### Robustness

Two tests exist because a tour must never take the application down: every step
must degrade to a readable message before anything has been computed, and none
may raise when handed junk. The GUI smoke test walks all eleven steps.

Suite 419 → 431 passing.

## Structure composition — what is actually in these files (2026-08-06)

Asked to check that loading is consistent across structures and that bound
elements are handled correctly. All 23 entries load, but auditing what is *in*
them turned up something worth fixing.

### Seven entries contain a second protein

Six carry three copies of **MDFIC** — a 21-residue cysteine-rich peptide
(`CCESSDCLEICMECCGICFPS`) that is a genuine auxiliary subunit inserting into the
pore (Zhou et al. 2023, already in our bibliography). 6B3R carries three
16-residue poly-UNK chains.

These are *protein*, so a protein mask includes them. Worse, MDFIC is numbered
**226–247**, which sits inside PIEZO1's own numbering. A selection keyed on
residue number alone would pool the two.

Nothing has actually been wrong: PIEZO1 resolves from residue 570 upward in
exactly the entries that carry MDFIC, so the ranges happen not to overlap. But
that is luck, not design, and the `> 300 C-alpha` rule that kept them out of
protomer blocks was a coincidence of scale rather than a statement about what
they are. There is now a test asserting no auxiliary residue number reaches a
protomer basis, so if a future entry resolves further into the N-terminus it
fails loudly instead of quietly averaging MDFIC into a mode.

`core/entities.py` classifies every atom. The principal-chain rule is
**relative to the largest chain**, not an absolute threshold — 4RAX is a lone
227-residue domain and is the entire structure, while a 21-residue peptide
beside a 1,280-residue protomer is not a protomer. Across all 23 entries the
classifier returns 1 or 3 protomers, never anything else.

The heterogens were looked up in the PDB chemical component dictionary rather
than inferred: PLX, PEE, P5S and L9Q are phospholipids, D12 is dodecane, NAG is
a glycan on PIEZO2.

### Display is now a choice, and only display

The Model panel lists the categories actually present with their atom counts,
each independently switchable. Only present categories get a control — a
permanent list mostly greyed out would say nothing about what you are looking
at, and what is in *this* file is the thing worth surfacing.

There is a test that hiding a category cannot change an analysis. Display and
computation are separate questions, and the analyses always use the channel
protomers whatever is drawn.

### Three structures that legitimately cannot answer

Extending the GUI smoke test beyond 8YEZ exposed that it assumed human
numbering and full side chains. Fixed, and the three failures turned out to be
correct behaviour needing to be *reported* rather than *fixed*:

- **11ZC** is a backbone-only model — C, CA, N, O and nothing else — so there
  are no sulfur atoms to measure.
- **6KG7** is PIEZO2. PIEZO1 numbering does not transfer to a paralogue; the
  second position is isoleucine.
- **4RAX** is one isolated domain, so there is no trimer to fit a dome to.

Each is now skipped with its reason. Conflating "cannot answer" with "answered
wrongly" is the same mistake the claims registry avoids, and it makes correct
behaviour look broken.

### An observation about the disulfide

Converting the curated C2411–C2415 disulfide into mouse numbering (2437/2441,
via the sanctioned path rather than a constant) shows both residues are cysteine
in every mouse entry — a 26-residue offset landing on cysteine twice is not
chance, so the numbering map is doing its job.

But the **bond is only modelled as formed in some of them**. Human 8YEZ and
8YFC give 2.04 Å, and so does mouse 4RAX — the 1.45 Å X-ray structure of the
isolated CED domain. The full-length mouse cryo-EM models put the same sulfurs
**5.2 Å apart** (7WLT, 7WLU) or 6.7 Å (6B3R). That is a resolution and
modelling difference between depositions rather than a biological one, and the
smoke test now records it instead of failing on it.

Suite 431 → 445 passing.

## Round 26 — the substitution finally enters the mechanics (2026-08-06)

Round 7's blind test failed and its diagnostic was precise: 99.8% of the
mechanical ΔΔG's variance was between position, not between substitutions. Two
rounds since have worked around that by adding other predictors. This one
attacks the cause.

### The cause was algebra, not statistics

The old model scaled **every contact of the mutated residue by one number**. So

    ΔΔG = (s − 1) · Q(position)

which is a rank-one product: the substitution enters only as a multiplicative
scalar, and the positional factor Q — contact count times local strain — varies
enormously between positions and not at all within one. Four substitutions at
R2456 could therefore differ only by a factor and had to rank every position
identically.

Seeing that written out matters, because it says no amount of refining `s` could
ever have helped. A better volume term, a substitution matrix, a fitted
sensitivity — all of them are still a scalar. The separability itself had to go.

### The repair

Scale each contact individually, by properties of the new residue *and of the
partner it touches*:

- **packing** — as before, but weighted by how close the contact is;
- **charge** — felt only at contacts with charged partners, so losing R2456's
  salt bridge softens the contact to an aspartate and not the one to a leucine;
- **hydrogen bonding** — donor against acceptor, both ways;
- **proline** — stiffens sequence-local contacts specifically, because that is
  where a backbone restraint is felt;
- **glycine** — softens whatever the side chain was mediating, since there is no
  longer a side chain.

Different substitutions now perturb *different subsets* of contacts. Measured
directly rather than assumed: the per-contact patterns for R2456H/K/P/C
correlate 0.62–0.98, where under the old model they were 1.00 by construction.

An elastic network is a mechanical model, and its springs are an effective
stiffness standing in for packing, hydrogen bonds and ion pairs together.
Letting that stiffness depend on charge is a statement about what the springs
represent, not a claim to have added electrostatics, and the module says so
rather than dressing it up.

### The criterion, and reporting both numbers

The success criterion was fixed in the roadmap before the work: within-position
variance above 20%, measured on the multiply-substituted positions. **It is met:
4.9% → 52.5%.**

Across *all* 35 substituted positions the figure is 2.4%, up from 0.8%. Both are
reported because the difference is instructive rather than embarrassing: 29 of
those positions carry a single substitution and contribute exactly zero
within-variance by construction, so including them pushes the statistic down for
arithmetic reasons that say nothing about the model. The criterion named the
multiply-substituted subset for exactly this reason.

### What this does not license

**No phenotype comparison was run.** The criterion was about variance
decomposition; whether the new distinctions point in the *right* direction is a
different question, it has not been asked, and asking it needs a new
pre-registration under `NEGATIVE_RESULT_PROTOCOL.md` §7. The temptation to
"just check" now that the model can tell substitutions apart is precisely what
that document exists to resist.

Round 7's recorded result is untouched: its script passes no sequence, so the
model falls back to the uniform scale, and there is a test asserting the two
agree to 1e-12.

Suite 445 → 462 passing; 18 documented numbers reproduce.

## Round 27 — more variants, and what kind of evidence they are (2026-08-06)

Round 20 measured the constraint: 42 variants for a large effect, 98 for a
medium one, against 25 usable. Round 22 then ran into it from the other side,
with six loss-of-function missense variants. This round went looking for more.

### ClinVar does not answer the question directly

It reports **pathogenicity**, not **direction**. A gain-of-function and a
loss-of-function variant are both "Pathogenic".

What makes direction recoverable for PIEZO1 specifically is that its two
diseases have opposite, well-established mechanisms — dehydrated hereditary
stomatocytosis is dominant and acts by slowed inactivation; generalised
lymphatic dysplasia is recessive and acts by loss of function. So the condition
implies the direction.

That is **weaker evidence than measuring the current**, because it assumes the
variant acts by the usual mechanism for its disease. So it is recorded as a
separate evidence level rather than pooled, and `build_analysis_set` defaults to
the measured set: a caller who does not think about evidence strength gets the
smaller, better-supported answer rather than the larger one.

### The ambiguity is bigger than expected

**Eleven of 63 directed records are reported under both diseases at once.**
ClinVar submitters routinely attach the whole gene's disease list to a variant,
so a single entry carries both "Dehydrated hereditary stomatocytosis" and
"Lymphatic malformation 6". Those get no direction. Resolving them by preferring
one disease would have manufactured 11 labels out of nothing.

The wild-type gate also did its job: three records disagree with Q92508 at their
stated position and were dropped, and 117 had a protein change that would not
parse unambiguously.

### An independent check that came free

Nine of the ClinVar variants are already in the curated set from
electrophysiology. The condition-based inference agrees with the measured
direction **8 times out of 9** — which is the only evidence available that the
inference is worth anything at all.

The one disagreement is V598M: we have it as gain-of-function, the inference
says loss. Our own record reads *"increased opening (one report); no change in
another"*, so the literature is mixed and neither side is simply wrong. It is
reported rather than resolved, because picking a side by fiat would hide that
from any test built on the set.

### What it bought

The directional missense set goes from **26 (20 GoF, 6 LoF) to 46 (27 GoF,
19 LoF)**. The loss-of-function class more than triples, which is precisely what
Round 22 was starved of, and the design becomes far better balanced.

Minimum detectable effect: **0.61 → 0.41**. Power at a large effect rises from
0.50 to 0.83.

### What it did not buy, stated because the roadmap asked

**A medium effect is still out of reach.** Power at δ = 0.28 is **0.49**, not
0.80, and reaching it at this ratio would need **104** variants against the 46
now available. The constraint has loosened by a useful amount; it has not
lifted.

The original `variants.json` is untouched. Round 7 and Round 22 reference it,
and growing it underneath a recorded result would invalidate that result while
nothing appeared to change.

Suite 462 → 477 passing.

## Round 28 — the footprint in the area change (2026-08-06)

Round 18 built the nonlinear footprint solver and showed the linearised version
is 3.5× too large at PIEZO1's 63° contact slope. Only `DomeModel` consumed it;
the gating energetics still ran on linear numbers.

### The quantity was wrong before the model was

ΔA is a **change between states**, not an absolute area. Round 3's framing —
"the footprint stores 622 nm² against the dome's 256" — invites treating a
stored area as the gating area, and they are different things. What tension does
work on is the *increase* in projected in-plane area on opening.

So both endpoints had to be measured: 7WLT closed at R_c 9.72 nm and contact
slope 1.992 (63.3°), 7WLU flattened at R_c 18.38 nm and slope 0.839 (40.0°).

### The correction is bigger on the difference than on either endpoint

This is the part I did not anticipate. The closed state sits at 63°, where the
small-slope expansion fails badly; the open state at 40°, where it is much less
wrong. Taking a difference between one badly overestimated number and one mildly
overestimated one amplifies the error rather than cancelling it:

- footprint stored closed: 622 → 179 nm² (3.5×, as Round 18 found)
- footprint stored open: 159 → 108 nm² (1.5×)
- **footprint released on opening: 463 → 71 nm² (6.5×)**

Total ΔA falls from 664 to 272 nm².

### The roadmap's question, answered

It asked to report the change *including if the linear version happened to agree
better*, since a wrong model can fit a right number. It did not: T₅₀ moves from
0.060 to 0.147 mN/m, toward the measured 2.7 ± 0.1 rather than away. There is a
test asserting that direction, because the opposite was a live possibility and
would have been the more interesting result.

### What the round actually established

**The correction does not close the gap, and that is worth more than the
correction.** Improving the membrane physics by a factor of six moved T₅₀ by a
factor of 2.4 and left it about **eighteen times below measurement**. The
structural ΔA is still 34× the functional 8 nm².

So the structural-versus-functional discrepancy that `docs/SCIENCE.md` has
carried since Round 3 is **not a membrane-modelling error**. No further
refinement of the footprint will fix it. The two numbers measure different
things: the functional ΔA is the area change along the gating reaction
coordinate, the structural one is the whole protein-plus-footprint deformation,
and only part of that is coupled to the gate.

Knowing that a candidate explanation has been ruled out is a real result, and it
took building the better model to rule it out.

Suite 477 → 489 passing; 19 documented numbers reproduce.
