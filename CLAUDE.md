# CLAUDE.md — PIEZO1 Dynamic Structural Simulator

Guidance for any Claude (or human) working in this repository.
**Read `@INTERFACE.md` before opening any source file.**

---

## 1. What this project is

An interactive desktop application that renders the PIEZO1 mechanosensitive ion
channel in 3D and drives it with real physics, so that a user can *see* how a
membrane-tension signal becomes an ionic current — and can test hypotheses
about mutants, lipids and drugs while doing so.

It is deliberately two things at once:

- **A learning instrument.** A student should be able to open it, rotate the
  trimer, click a blade, watch the dome flatten under tension, and understand
  why PIEZO1's shape *is* its mechanism.
- **A research instrument.** A working scientist should be able to load a
  disease variant, compute how it changes the elastic-network mode that couples
  blade motion to the pore, and get a number they can put in a grant or a
  paper — with the provenance of every parameter visible.

## 2. Project aims

**A1. Structural fidelity.** Build the best available full-length model of human
PIEZO1 by combining experimental cryo-EM (which resolves roughly residues
570–2521) with predicted structure for the unresolved distal blade, and keep
the seams honest and visible rather than hidden.

**A2. Physics, not animation.** Every motion shown must come from a model with a
citable basis — elastic network normal modes, Helfrich membrane mechanics, a
tension-dependent Markov gating scheme — never from a hand-tuned keyframe. If
something is interpolated for visual purposes, the UI must say so.

**A3. Provenance everywhere.** Every number in the app carries its source. A
tooltip on a rate constant gives the PMID. A domain boundary states whether it
came from UniProt, from a paper, or from our own structural analysis. Values we
could not verify are labelled `UNVERIFIED` rather than quietly rounded.

**A4. The variant-to-mechanism path.** The core scientific workflow: choose a
variant → see it on the structure → see which domain and which normal mode it
perturbs → see the predicted change in gating energetics → compare with the
measured phenotype. This is the pipeline that could generate a real hypothesis.

**A5. Reproducibility.** A fresh clone plus `scripts/create_env.sh` plus
`python -m piezo1.io.fetch` must reproduce the entire working state. No data
file is sacred; all of it is regenerable.

**A6. Extensibility over multiple iterations.** Small focused modules, a stable
public interface documented in `INTERFACE.md`, and a session log that explains
*why* past decisions were made, so the next iteration does not re-litigate them.

## 3. Non-goals

- Not a general-purpose molecular viewer. PyMOL and ChimeraX exist and are
  better at that. This app is PIEZO-specific and opinionated.
- Not an all-atom MD engine. OpenMM is available for local refinement, but the
  interactive dynamics are coarse-grained by design — that is what makes a
  2500-residue trimer responsive.
- Not a clinical tool. Nothing here is validated for diagnosis.

## 4. Scientific ground truth

The science this project rests on is documented in `docs/SCIENCE.md`, with the
underlying literature dossiers in `ref/research/` (git-ignored, regenerable).
Key anchors, each traceable to a citation:

| Quantity | Value | Source |
|---|---|---|
| Dome radius of curvature (closed) | 10.2 nm | Haselwandter & MacKinnon 2018, eLife |
| Membrane bending modulus | 20–25 k_BT | Haselwandter 2018; Dixit 2025 |
| Footprint decay length | 14 nm | Haselwandter & MacKinnon 2018 |
| Half-activation tension | 2.7–4.7 mN/m | see `docs/SCIENCE.md` |
| Unitary conductance | ~25–30 pS | Shi 2020; Vaisey & MacKinnon 2026 |
| Yoda1 EC50 (human) | 26.6 µM | Syeda 2015, eLife |
| GsMTx4 K_D | 155 nM | Bae 2011, Biochemistry |

Our own structural analysis reproduces the published dome geometry
(9.7 nm from PDB 7WLT), which is the standing regression test that the geometry
pipeline is correct.

## 5. Working conventions

**Residue numbering.** Human PIEZO1 is 2521 aa (UniProt Q92508); mouse Piezo1 is
2547 aa (E2JF22). *Most functional literature uses mouse numbering; most disease
variants use human.* The offset is **not** a constant. Never hard-code a
conversion — always go through `piezo1.core.sequence`. Any residue number in
code, data or UI must state its numbering system.

**Every number is a registered parameter. This is not optional.**

Any number a calculation depends on lives in the registry, not in a literal.
A constant written into a function default is invisible: it cannot be listed,
shown to a user, or traced to a paper — and several numbers this project has
had to correct were invisible in exactly that way.

- Declare it in `scripts/build_parameters.py` with a **unit**, **bounds**, a
  **kind** (`physical` / `empirical` / `method` / `convention`), a
  **description**, and a **citation**. Rebuild with
  `python scripts/build_parameters.py`.
- The citation must be a key that exists in `references.json`. If there is no
  paper, use one of the sentinels (`method_choice`, `measured_here`, `derived`,
  `convention`, `unverified`) — each of which **obliges you to say why** in
  `source_note`. The build refuses to write the resource otherwise.
- Consume it as `_P.value("key")` — via `field(default_factory=...)` for a
  dataclass field, or a `None` default resolved in the body for a function
  argument. Resolve at **call time**, so an override takes effect on the next
  call rather than at import.
- `python -m piezo1.parameter_audit` must stay clean. It scans `physics/`,
  `structure/` and `analysis/` and fails on any numeric literal that is neither
  registered nor listed in `EXEMPT`/`EXEMPT_NAMES` **with a stated reason**.
  Tolerances, seeds, iteration caps and zero-initialised fields are legitimately
  exempt; a physical quantity never is. `tests/test_parameters.py` enforces this
  on every run.
- Anything that records a result must record the parameter set with it.
  Reports carry a banner when a parameter is non-default, and
  `verify_claims` **refuses to run** against a modified registry — the
  documented numbers were produced at the defaults, and comparing against
  anything else would report drift the user caused.

**A checking instrument is a measuring instrument. Calibrate it first.**

The most expensive errors in this project have all had the same shape: an
*alternative* route, written to check the main one, was itself wrong — and it
returned a plausible number rather than an error, so the disagreement looked
like a finding. A spheroid fitter that would have reported 89% model error. A
document checker that could not read the Unicode minus its own documents use. A
parameter probe whose "no effect" came from coordinates too diffuse to form a
single alpha sphere. A conservation cross-check whose bias I first explained
with the wrong mechanism.

So an uncalibrated checker is worse than none: it manufactures findings.

- Before a cross-check, re-derivation, audit or probe is believed against real
  data, run it on a case whose answer is **known independently** — an analytic
  shape, a planted signal, a true null, an enumerable exact value, or a
  deliberately inert input.
- The calibration must be able to **fail**. A check that would pass on a broken
  instrument asserts nothing; if there is no input that makes it say "no", it is
  not a calibration.
- Register it. Every public callable in `analysis/crosscheck.py`,
  `crosscheck_methods.py`, `model_error.py`, `uncertainty.py`, `validation.py`,
  `design.py`, `provenance_chain.py` and `parameter_effect.py` must appear in
  `CALIBRATED` in `tests/test_calibration.py` with the test that calibrates it.
  `test_every_checking_instrument_has_a_calibration` fails otherwise, and
  `test_named_calibrating_tests_exist` fails if the named test does not exist.
- When a checker disagrees with the pipeline, **suspect the checker first**.
  Historically it has been wrong more often than the thing it was checking.

**Code style.**
- Files stay under 500 lines. If one is heading past that, split it first.
- Structure-of-arrays with numpy, not per-atom Python objects. A PIEZO1 trimer
  is 120k atoms; per-object designs will not survive.
- Selections are boolean masks over atom arrays.
- Physics modules must not import anything from `render/` or `ui/`. The
  dependency arrow points one way: `io → core → structure → physics → analysis`,
  with `render` and `ui` consuming all of them.
- Comments explain constraints and units, not what the next line does.

**Data hygiene.**
- `ref/` and `data/` are git-ignored. Nothing downloaded is ever committed.
- Curated, hand-authored annotation lives in `piezo1/resources/` and *is*
  committed, because it is authored content with provenance, not a download.
- Units: Angstrom for coordinates, nm for reported dome geometry (matching the
  literature), k_BT for energies, mN/m for tension. Convert at the boundary and
  say so in the docstring.

**Environment.** Use the `piezo1` conda environment:
`source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate piezo1`

## 6. Documentation duties

When the structure of the project changes, update in the same commit:
- `INTERFACE.md` — the navigation map. Non-negotiable.
- `SESSION_LOG.md` — what was done and, more importantly, *why*.
- `README.md` — if the user-visible feature set changed.
- `docs/SCIENCE.md` — if a scientific parameter or its source changed.

@INTERFACE.md
