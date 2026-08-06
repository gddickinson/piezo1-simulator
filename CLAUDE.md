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
