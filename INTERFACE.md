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
| `fetch.py` | Cached downloaders for RCSB mmCIF, EMDB maps, AlphaFold DB, UniProt, PubChem. Idempotent; skips files already present. | `fetch_pdb()`, `fetch_alphafold()`, `fetch_uniprot()`, `fetch_ligand()`, `fetch_all()` | 📋 |
| `registry.py` | Curated catalogue of PIEZO structures: state, species, resolution, modelled range, what each is good for. | `StructureRecord`, `REGISTRY`, `lookup()` | 📋 |

### `piezo1/core/` — data model

| File | Purpose | Key names | Status |
|---|---|---|---|
| `structure.py` | The central structure-of-arrays container. Boolean-mask selections, residue-level index, vdW radii, element colours, PDB output. | `Structure`, `AA3TO1`, `ELEMENT_RADII`, `ELEMENT_COLORS` | ✅ |
| `sequence.py` | Sequence alignment and the **human↔mouse residue numbering map**. The only sanctioned way to convert between numbering systems. | `align()`, `NumberingMap`, `human_to_mouse()`, `mouse_to_human()` | 📋 |
| `annotations.py` | Loads `resources/*.json` into typed records; maps annotations onto a loaded `Structure`. | `Domain`, `Variant`, `LigandSite`, `Annotations`, `load_annotations()` | 📋 |

### `piezo1/structure/` — structural operations

| File | Purpose | Key names | Status |
|---|---|---|---|
| `superpose.py` | Kabsch superposition, RMSD, C3 symmetry-axis recovery, axis-to-z alignment. | `kabsch()`, `superpose()`, `rmsd()`, `SymmetryAxis`, `detect_c3_axis()`, `align_axis_to_z()`, `rotation_matrix()` | ✅ |
| `geometry.py` | **Membrane-dome measurement.** Sphere fitting, radial height profile, dome depth / area / excess area. Reproduces published dome curvature. | `fit_sphere()`, `SphereFit`, `radial_profile()`, `RadialProfile`, `DomeGeometry`, `measure_dome()` | ✅ |
| `hybrid.py` | Assembles the full-length model: experimental core + AlphaFold distal blade, with the seam recorded and renderable. | `build_hybrid_model()`, `HybridModel` | 📋 |
| `morph.py` | Conformational interpolation between curved and flattened states with bond-geometry restraints. | `morph()`, `MorphTrajectory` | 📋 |
| `pore.py` | Pore radius profile along the conduction axis (HOLE-equivalent, self-contained). | `pore_profile()`, `PoreProfile` | 📋 |

### `piezo1/physics/` — the simulation engine

| File | Purpose | Key names | Status |
|---|---|---|---|
| `anm.py` | Anisotropic network model: sparse Hessian, low-frequency modes, C3 symmetry adaptation. | `ANM`, `build_hessian()`, `calc_modes()` | 📋 |
| `modes.py` | Mode analysis: overlap with an observed transition, cumulative overlap, deformation projection. | `overlap()`, `cumulative_overlap()`, `project()` | 📋 |
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
| `camera.py` | Trackball camera, projection, screen↔world unprojection for picking. | `Camera` | 📋 |
| `scene.py` | Draw-list assembly, per-frame uniform management, render passes. | `Scene`, `RenderPass` | 📋 |
| `representations.py` | Spheres, sticks, backbone tube, cartoon, surface, membrane slab. | `Representation` and subclasses | 📋 |
| `geometry_builders.py` | CA spline, ribbon frames, cartoon extrusion, membrane mesh. | `spline_ca()`, `build_cartoon()`, `build_membrane()` | 📋 |
| `colormaps.py` | Colouring by chain, domain, variant, pLDDT, conservation, mode amplitude. | `ColorScheme` and subclasses | 📋 |
| `shaders/` | GLSL 4.1 sources, loaded at runtime. | `sphere_impostor.vert/frag`, `cylinder_impostor.*`, `mesh.*` | 📋 |

### `piezo1/ui/` — PyQt6 application

| File | Purpose | Key names | Status |
|---|---|---|---|
| `main_window.py` | Application shell: menus, docks, status bar, session state. | `MainWindow` | 📋 |
| `gl_widget.py` | `QOpenGLWidget` hosting the moderngl context; input handling. | `ViewportWidget` | 📋 |
| `panels/` | Dockable control panels — structure, representation, variants, ligands, physics, kinetics, sequence. | one class per panel | 📋 |

### `piezo1/resources/` — curated data (committed)

| File | Contents | Status |
|---|---|---|
| `uniprot_human.json` | Distilled UniProt Q92508: sequence, 38 TM segments, topology, PTMs, 26 natural variants, disulfide, coiled coil. Built by `scripts/build_uniprot_annotations.py`. | ✅ |
| `uniprot_mouse.json` | Same for mouse E2JF22 (2547 aa). | ✅ |
| `domains.json` | PIEZO1 architectural domains (blade THU repeats, beam, anchor, OH, IH, cap/CED, CTD) with residue ranges in both numbering systems and per-entry provenance. | 📋 |
| `variants.json` | Curated variant table: classification, phenotype, functional effect, PMID. | 📋 |
| `ligands.json` | Yoda1, Yoda2, Jedi1/2, Dooku1, GsMTx4 and lipids with chemistry and binding-site residues. | 📋 |
| `structures.json` | Structure registry metadata. | 📋 |

---

## `scripts/`

| File | Purpose | Status |
|---|---|---|
| `create_env.sh` | Creates the `piezo1` conda environment with the full stack. | ✅ |
| `build_uniprot_annotations.py` | Distils UniProt JSON into committed resource files. | ✅ |
| `fetch_data.py` | Downloads every structure, sequence and ligand the app needs. | 📋 |

## `tests/`

| File | Covers | Status |
|---|---|---|
| `test_cif_reader.py` | Parser correctness on real files, including quoted tokens and HETATM. | 📋 |
| `test_geometry.py` | Sphere fitting on synthetic caps; dome curvature regression against 7WLT. | 📋 |
| `test_superpose.py` | Kabsch round-trip, C3 axis recovery. | 📋 |

## `docs/`

| File | Purpose | Status |
|---|---|---|
| `SCIENCE.md` | The scientific basis: mechanism, parameters, provenance, known controversies. | 📋 |
| `ARCHITECTURE.md` | Why the code is shaped this way; the rendering approach in detail. | 📋 |
