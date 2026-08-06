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
| `membrane.py` | Monge-gauge Helfrich solver for the membrane footprint around the dome. | `MembraneFootprint`, `solve_shape()`, `footprint_energy()` | 📋 |
| `dome.py` | Dome-model energetics: ΔE = −T·ΔA, tension–area coupling, state free energies. | `DomeModel`, `gating_energy()` | 📋 |
| `kinetics.py` | Tension-dependent Markov gating model; stochastic and deterministic solutions; simulated patch-clamp traces. | `GatingModel`, `simulate_trace()`, `open_probability()` | 📋 |

### `piezo1/analysis/` — interpretation

| File | Purpose | Key names | Status |
|---|---|---|---|
| `variants.py` | Maps variants onto structure, reports domain context, contacts lost/gained, predicted mode perturbation. | `map_variants()`, `VariantImpact` | 📋 |
| `contacts.py` | Residue contact maps, interface detection, contact changes between states. | `contact_map()`, `interface_residues()` | 📋 |
| `pockets.py` | Grid-based pocket detection for ligand sites (zero external dependencies). | `find_pockets()`, `Pocket` | 📋 |
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

## `docs/`

| File | Purpose | Status |
|---|---|---|
| `SCIENCE.md` | The scientific basis: mechanism, parameters, provenance, open gaps. | ✅ |
| `img/` | Generated figures (`make_figures.py`, `screenshot_app.py`). | ✅ |
| `ARCHITECTURE.md` | Why the code is shaped this way; the rendering approach in detail. | 📋 |
