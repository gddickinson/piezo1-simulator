# INTERFACE.md — navigation map

The top-level map of this project. **Read this before opening any source file.**
Status markers: ✅ implemented · 🚧 in progress · 📋 planned.

---

## Dependency direction

```
io ──▶ core ──▶ structure ──▶ physics ──▶ analysis
                                  │
                   render ◀───────┴───────▶ ui
```

`physics` and below never import from `render` or `ui`. That keeps the science
testable headlessly and lets the whole engine be driven from a notebook.

---

## Repository layout

| Path | Contents |
|---|---|
| `piezo1/` | The application package |
| `scripts/` | One-shot maintenance scripts (env setup, resource building) |
| `docs/` | Science notes, architecture notes, figures |
| `tests/` | pytest suite |
| `ref/` | Downloaded reference material — **git-ignored, regenerable** |
| `data/` | Fetched and derived artefacts — **git-ignored, regenerable** |

---

## `piezo1/` — module by module

### Top level

| File | Purpose | Key names | Status |
|---|---|---|---|
| `config.py` | All filesystem paths, physical constants, runtime settings. Every module imports paths from here. | `PROJECT_ROOT`, `REF_DIR`, `STRUCTURE_DIR`, `RESOURCE_DIR`, `HUMAN_ACC`, `MOUSE_ACC`, `KT_ROOM`, `RenderSettings`, `AppSettings`, `SETTINGS`, `ensure_dirs()` | ✅ |

### `piezo1/io/` — data acquisition and parsing

| File | Purpose | Key names | Status |
|---|---|---|---|
| `cif_reader.py` | Fast mmCIF/PDB coordinate readers producing numpy arrays. Handles quoting, multi-line text fields and model selection. ~0.6 s for a 34k-atom trimer. | `read_cif_atoms()`, `read_pdb_atoms()`, `read_structure_file()`, `parse_cif_categories()` | ✅ |
| `fetch.py` | Cached downloaders for RCSB mmCIF, AlphaFold DB, UniProt, PubChem. Idempotent. AlphaFold versions discovered from the API, never guessed. | `fetch_pdb()`, `fetch_alphafold()`, `fetch_uniprot()`, `fetch_ligand()`, `fetch_all()`, `DEFAULT_PDB_IDS` | ✅ |
| `registry.py` | Curated catalogue of 21 PIEZO structures: state, gating, resolved range, ligands, citation, and what each is recommended for. | `StructureRecord`, `Registry`, `load_registry()` | ✅ |

### `piezo1/core/` — data model

| File | Purpose | Key names | Status |
|---|---|---|---|
| `structure.py` | The central structure-of-arrays container. Boolean-mask selections, residue-level index, vdW radii, element colours, PDB output. | `Structure`, `AA3TO1`, `ELEMENT_RADII`, `ELEMENT_COLORS` | ✅ |
| `sequence.py` | Sequence alignment and the **human↔mouse residue numbering map**, built from a real global alignment. The only sanctioned way to convert between numbering systems. | `align_global()`, `NumberingMap`, `load_numbering_map()`, `human_to_mouse()`, `mouse_to_human()` | ✅ |
| `annotations.py` | Loads `resources/*.json` into typed records carrying provenance and confidence. Represents "we do not know" explicitly. | `Domain`, `ResidueGroup`, `Variant`, `Annotations`, `load_annotations()` | ✅ |

### `piezo1/structure/` — structural operations

| File | Purpose | Key names | Status |
|---|---|---|---|
| `superpose.py` | Kabsch superposition, RMSD, C3 axis recovery, and **protomer correspondence matching** — deposited chain labels are not a reliable guide to rotational order. | `kabsch()`, `superpose()`, `rmsd()`, `SymmetryAxis`, `detect_c3_axis()`, `match_protomers()`, `ProtomerMatch`, `align_axis_to_z()`, `rotation_matrix()` | ✅ |
| `geometry.py` | **Membrane-dome measurement.** Sphere fitting, radial height profile, dome depth / area / excess area. Reproduces published dome curvature. | `fit_sphere()`, `SphereFit`, `radial_profile()`, `RadialProfile`, `DomeGeometry`, `measure_dome()` | ✅ |
| `hybrid.py` | Assembles the full-length model: experimental core + AlphaFold distal blade, with the seam recorded and renderable. | `build_hybrid_model()`, `HybridModel` | 📋 |
| `morph.py` | Conformational interpolation between endpoints: linear, distance-restrained, and elastic-network-subspace methods, each reporting its own bond-geometry error. | `morph()`, `MorphTrajectory`, `prepare_endpoints()`, `restrained_morph()`, `modal_morph()` | ✅ |
| `pore.py` | Pore-radius profile along the conduction axis: Apollonius clearance maximisation per slice with a leash constraint, bottleneck and constriction detection, per-slice lining residues. | `pore_profile()`, `PoreProfile`, `PoreSlice` | ✅ |

### `piezo1/physics/` — the simulation engine

| File | Purpose | Key names | Status |
|---|---|---|---|
| `anm.py` | Anisotropic network model: sparse Hessian, shift-invert Lanczos modes, C3 irreducible-representation labelling, disconnected-network detection. | `ANM`, `ModeSet`, `build_hessian()`, `SPRING_MODELS` | ✅ |
| `modes.py` | Further mode analysis beyond `ModeSet.overlap` / `.cumulative_overlap` / `.msf` / `.collectivity`, which already live on the mode set. | `hinge_sites()`, `project()` | 📋 |
| `membrane.py` | Linearised Helfrich footprint: exact K₀ solution, second-order-convergent finite-difference solver (coupled second-order form), energy, excess area, decay-length recovery, and an explicit small-slope validity check. | `MembraneParameters`, `FootprintSolution`, `solve_footprint()`, `analytic_profile()`, `analytic_energy()`, `decay_length()` | ✅ |
| `dome.py` | Two-state dome energetics ΔG = ΔG₀ − T·ΔA, open probability, T₅₀, footprint coupling, and a side-by-side comparison of functional vs structural ΔA estimates. | `DomeModel`, `DomeGeometrySummary`, `open_probability()`, `half_activation_tension()`, `PUBLISHED_AREA_ESTIMATES` | ✅ |
| `kinetics.py` | Four-state tension-dependent Markov gating (Young et al. 2023): rate matrix with enforced microscopic reversibility, steady state, step/ramp protocols, Gillespie single-channel simulation, and fold-change-calibrated mutant presets. | `GatingModel`, `GatingResult`, `MUTANT_PRESETS`, `STATE_NAMES` | ✅ |

### `piezo1/analysis/` — interpretation

| File | Purpose | Key names | Status |
|---|---|---|---|
| `measure.py` | Distance/angle/dihedral, radius of gyration, principal axes, helix axis, tilt and crossing angles, Shrake–Rupley SASA, buried interface area, and the pore hydrophobicity profile. | `distance()`, `angle()`, `dihedral()`, `sasa()`, `buried_area()`, `helix_axis()`, `tilt_angle()`, `hydrophobicity_profile()` | ✅ |
| `interactions.py` | Hydrogen bonds, salt bridges, hydrophobic contacts, π-stacking, cation–π and disulfides, with published geometric criteria; cross-selection mode for interfaces; state-to-state comparison. | `detect_interactions()`, `Interaction`, `InteractionSet`, `compare_interactions()`, `CUTOFFS` | ✅ |
| `variant_impact.py` | Predicts a variant's mechanical effect on gating as ½·dᵀ(H_mut − H_wt)·d — the change in elastic cost of the observed gating motion. Volume-based spring perturbation, all protomers, coverage reported. | `VariantImpactModel`, `VariantPrediction`, `spring_scale_from_volume()`, `RESIDUE_VOLUME` | ✅ |
| `variants.py` | Maps variants onto structure, reports domain context, contacts lost/gained. | `map_variants()`, `VariantImpact` | 📋 |
| `ensemble.py` | Builds a cross-species, coverage-matched, protomer-corrected structure ensemble and runs PCA on it; compares principal components with elastic-network modes by overlap, subspace overlap and RWSIP. | `build_ensemble()`, `StructureEnsemble`, `PCAResult`, `rwsip()`, `subspace_overlap()`, `DEFAULT_EXCLUSIONS` | ✅ |
| `allostery.py` | Perturbation response scanning, dynamic cross-correlation, correlation-weighted contact networks, shortest allosteric paths, via-point detour cost and path betweenness. | `perturbation_response()`, `PRSResult`, `cross_correlation()`, `build_network()`, `allosteric_path()`, `detour_cost()`, `path_betweenness()` | ✅ |
| `contacts.py` | Residue contact maps, interface detection, contact changes between states. | `contact_map()`, `interface_residues()` | 📋 |
| `pockets.py` | Delaunay alpha-sphere pocket detection (fpocket construction, reimplemented in numpy) with a burial filter that stops surface percolation, Monte-Carlo union volumes, ray-cast buriedness, and resolved-ligand contact mapping. | `find_pockets()`, `Pocket`, `alpha_spheres()`, `AlphaSpheres`, `ligand_contact_residues()` | ✅ |
| `validation.py` | Non-parametric statistics for the blind test: permutation test with the (r+1)/(n+1) convention, Cliff's delta with a bootstrap CI, and tie-averaged AUROC. Implemented directly so the conventions are visible and testable. | `permutation_test()`, `cliffs_delta()`, `bootstrap_cliffs_delta()`, `auroc()` | ✅ |
| `docking.py` | Optional AutoDock Vina integration; degrades gracefully when absent. | `dock()`, `available()` | 📋 |

### `piezo1/render/` — OpenGL 4.1 renderer

Ray-cast impostor rendering (the PyMOL/ChimeraX technique): spheres and
cylinders are drawn as screen-space quads whose fragment shader solves the
ray–quadric intersection and writes `gl_FragDepth`. This gives pixel-perfect
geometry at a fraction of the triangle count.

| File | Purpose | Key names | Status |
|---|---|---|---|
| `camera.py` | Quaternion trackball camera. `frame()` solves for the exact tight-fit distance rather than using a bounding sphere. | `Camera`, `perspective()`, `look_at()`, `quat_to_matrix()` | ✅ |
| `scene.py` | Named batches, shared uniforms, opaque-then-transparent render order. | `Scene`, `Light` | ✅ |
| `representations.py` | Turns a `Structure` into GPU batches; owns style and colouring state. | `MolecularView`, `Style`, `ColorBy` | ✅ |
| `geometry_builders.py` | Swept tubes, cartoon ribbons with arrowheads, membrane surface of revolution. | `Mesh`, `build_tube()`, `build_cartoon()`, `build_membrane_mesh()`, `build_disc()` | ✅ |
| `colormaps.py` | Chain, domain, secondary-structure, B-factor, pLDDT, element and scalar-value colouring. | `DomainPalette`, `load_domain_palette()`, `domain_colors()`, `value_colors()`, `plddt_colors()` | ✅ |
| `animation.py` | Offscreen frame capture, ease-in-out and ping-pong timing, burnt-in captions and progress bar, GIF (shared adaptive palette) and MP4 writers. | `Animator`, `AnimationSpec`, `ease_in_out()`, `ping_pong()`, `write_gif()`, `write_mp4()` | ✅ |
| `shaders/` | GLSL 4.1 sources, loaded at runtime. | `sphere.vert/frag`, `cylinder.vert/frag`, `mesh.vert/frag` | ✅ |

### `piezo1/ui/` — PyQt6 application

| File | Purpose | Key names | Status |
|---|---|---|---|
| `main_window.py` | Application shell: docks, menus, structure loading, highlighting, click-to-identify. | `MainWindow`, `main()` | ✅ |
| `physics_controller.py` | Dome measurement, threaded mode calculation, mode animation and displacement colouring. | `PhysicsController`, `ModeWorker` | ✅ |
| `morph_controller.py` | Builds and plays back the curved-to-flat morph. | `MorphController` | ✅ |
| `gl_widget.py` | `QOpenGLWidget` hosting the moderngl context; input, picking, animation ticks. | `ViewportWidget`, `configure_surface_format()` | ✅ |
| `panels/` | `structure_panel` (chooser + appearance), `annotation_panel` (domains, sites, variants), `physics_panel` (dome, modes, animation). | `StructurePanel`, `AnnotationPanel`, `PhysicsPanel` | ✅ |

### `piezo1/resources/` — curated data (committed)

| File | Contents | Status |
|---|---|---|
| `uniprot_human.json` | Distilled UniProt Q92508: sequence, 38 TM segments, topology, PTMs, 26 natural variants, disulfide, coiled coil. Built by `scripts/build_uniprot_annotations.py`. | ✅ |
| `uniprot_mouse.json` | Same for mouse E2JF22 (2547 aa). | ✅ |
| `domains.json` | 17 architectural domains with ranges in both numbering systems, provenance (uniprot / derived-by-rule / literature) and confidence. | ✅ |
| `variants.json` | 68 curated variants, every wild-type residue verified against Q92508, each annotated with which structures resolve it. | ✅ |
| `functional_residues.json` | 37 residues in 11 groups: hydrophobic gate, selectivity glutamates, CTD constrictions, Yoda1 pocket, PIP2 cluster, basic patches. | ✅ |
| `numbering_human_mouse.json` | Cached human↔mouse alignment map. | ✅ |
| `ligands.json` | Yoda1, Yoda2, Jedi1/2, Dooku1, GsMTx4 and lipids with chemistry and binding-site residues. | 📋 |
| `structures.json` | Registry of 21 structures with state, resolution, coverage, ligands, citation. | ✅ |

---

## `scripts/`

| File | Purpose | Status |
|---|---|---|
| `create_env.sh` | Creates the `piezo1` conda environment with the full stack. | ✅ |
| `build_uniprot_annotations.py` | Distils UniProt JSON into committed resource files. | ✅ |
| `build_domains.py` | Authors `domains.json`. | ✅ |
| `build_functional_residues.py` | Authors `functional_residues.json`, verifying each residue against the sequence. | ✅ |
| `build_variants.py` | Promotes researched variants into `variants.json` behind a validation gate. | ✅ |
| `build_structure_registry.py` | Authors `structures.json`. | ✅ |
| `render_offscreen.py` | Headless render to PNG; also a renderer regression check. | ✅ |
| `make_figures.py` | All README/doc figures, on shared scale and orientation. | ✅ |
| `run_validation.py` | Executes the pre-registered blind test and writes the result. | ✅ |
| `make_animations.py` | The seven-animation library: gating morph, normal mode, ligand and lipid sites, variant context. | ✅ |
| `build_references.py` | Resolves the bibliography from Europe PMC behind a title-verification gate; downloads open-access full texts. | ✅ |
| `screenshot_app.py` | Drives the real GUI as a smoke test and captures screenshots. | ✅ |
| Data download | Use `python -m piezo1.io.fetch`. | ✅ |

## `tests/`

| File | Covers | Status |
|---|---|---|
| `conftest.py` | Fixtures; skips rather than fails when data is not downloaded. | ✅ |
| `test_cif_reader.py` | Tokenizer whitespace/quoting, column alignment, residue indexing, selections. | ✅ |
| `test_geometry.py` | Sphere fitting on synthetic spheres and caps; dome curvature regression against the published 10.2 nm; curved vs flat separation. | ✅ |
| `test_superpose.py` | Kabsch round-trip, reflection exclusion, C3 exactness, reversed-handedness detection. | ✅ |
| `test_anm.py` | Hessian symmetry, zero modes, disconnected networks, symmetry characters, and the gating-overlap result. | ✅ |
| `test_sequence_and_resources.py` | Ten cross-species equivalences, non-constant offset, resource integrity. | ✅ |
| `test_morph.py` | Endpoint preparation, the chord artefact and its removal, modal capture fraction, mismatched-mode-set rejection. | ✅ |
| `test_pore.py` | Leash enforcement, closed-vs-open bottleneck, and agreement between detected constrictions and the curated gate/CTD residues. | ✅ |
| `test_kinetics.py` | Published rate values, microscopic reversibility, generator validity, half-activation against measured T50, Gillespie-vs-analytic agreement, and mutant direction. | ✅ |
| `test_measure.py` | Geometry on analytic shapes, SASA of an isolated atom against 4πr², determinism, and the pore-helix tilt result. | ✅ |
| `test_interactions.py` | The annotated disulfide, the R2456–E2117 inter-protomer salt bridge, cutoff enforcement, and the donor–donor exclusion. | ✅ |
| `test_membrane.py` | Unit conversion, the κ/γ/λ triple, exact-vs-numerical profile and energy, second-order convergence, small-slope validity, and Cox's T₅₀ round trip. | ✅ |
| `test_ensemble.py` | Shared-basis construction, paralogue exclusion, reversed-protomer detection, PC1-as-gating-coordinate, and A-mode dominance against a random control. | ✅ |
| `test_pockets.py` | Circumsphere geometry against a known tetrahedron, percolation prevention, union-not-sum volumes, gate/anchor recovery, and the Yoda1 groove-versus-cavity result. | ✅ |
| `test_validation.py` | Statistical instruments against known cases, plus a pin on the published Round 7 null result so a predictor change cannot silently move it. | ✅ |
| `test_variant_impact.py` | The quadratic-form identity against an explicit Hessian, sign conventions, all-protomer mutation, and honest coverage reporting. Deliberately contains no phenotype comparison. | ✅ |
| `test_allostery.py` | Correlation-matrix validity, chunking invariance, the anchor on the optimal route, the beam as a near-degenerate alternative, and the invariant that a constrained path can never beat a free one. | ✅ |

## `docs/`

| File | Purpose | Status |
|---|---|---|
| `SCIENCE.md` | The scientific basis: mechanism, parameters, provenance, open gaps. | ✅ |
| `PREREGISTRATION.md` | The frozen hypothesis, statistic and decision rule for the Round 7 blind test, written before any comparison was run. | ✅ |
| `VALIDATION.md` | The Round 7 result: a null result, reported in the pre-registered order with a post-hoc diagnostic of why. | ✅ |
| `REFERENCES.md` | Generated bibliography, 51 verified references. | ✅ |
| `img/` | Generated figures (`make_figures.py`, `screenshot_app.py`). | ✅ |
| `ARCHITECTURE.md` | Why the code is shaped this way; the rendering approach in detail. | 📋 |
