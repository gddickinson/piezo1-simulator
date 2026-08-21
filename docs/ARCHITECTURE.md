# Architecture — why the code is shaped this way

`INTERFACE.md` says *what* each module contains. This says **why**, and it is
organised around the constraints that forced each decision rather than the
structure that resulted. Almost every rule below exists because something went
wrong first; where that is so, the incident is named, because a rule whose
reason has been forgotten is the next thing to be undone.

---

## 1. The dependency arrow points one way

```
io ──▶ core ──▶ structure ──▶ physics ──▶ analysis
                                  │
                   render ◀───────┴───────▶ ui
```

**The constraint:** the science has to be testable with no display attached.
A trimer is 120,000 atoms and the suite runs on every commit; if computing a
dome radius required a GL context, it could not run in CI, in a notebook, or on
a machine over SSH.

**What forced it:** `analysis` was importing `protomer_blocks` from
`ui/model_utils.py`. Nothing failed — the GUI was always up when anyone
noticed — but the import meant a headless analysis pulled in Qt. The function
moved to `structure/protomers.py` and the old module became a re-export.

The same shape recurred in Round 78 inside `analysis`: `report_tags` imported a
helper from `report` while `report` imported `report_tags` at the bottom of the
file. `import report_tags` failed in a fresh interpreter and the suite never
noticed, because something always imported `report` first. A cycle only bites
the person who reaches for one module directly — a notebook user.
`tests/test_imports.py` now imports every module alone.

**What it costs:** controllers in `ui/` are thin to the point of being dull.
They marshal arguments and hand off. That is deliberate: logic in a controller
is logic that only the GUI can reach, and Round 58 found a scoring function in
exactly that position.

---

## 2. Structure-of-arrays, not objects

`Structure` holds parallel numpy arrays — `xyz`, `element`, `res_seq`, `chain`,
`b_factor` — and a selection is a boolean mask over them.

**The constraint:** 120,000 atoms. A per-atom Python object costs roughly 50–100
bytes of overhead each before any data, and every operation becomes an
interpreted loop. Framing, superposing and measuring all become single vectorised
expressions in this layout, and the pore search — which evaluates clearance at
thousands of probe positions — is only tractable at all because the distance
computation is one call into a KD-tree.

**What it costs:** anything that wants "the atoms of residue 2456" writes a mask
rather than following a pointer, and a mask is easy to get subtly wrong. The
compensation is that masks compose with `&` and `|`, so a selection reads as the
sentence it came from, and `entities.py` exists to make "which atoms are the
channel" a computed answer rather than a chain-label guess.

---

## 3. The renderer draws quadrics, not triangles

A sphere is a four-vertex screen-facing quad whose fragment shader solves the
ray–sphere intersection and writes `gl_FragDepth`. Cylinders likewise.

**The constraint:** a tessellated sphere at a quality that survives zooming is
several hundred triangles. At 120,000 atoms that is tens of millions of
triangles per frame for a picture that is mostly featureless. The impostor
draws two triangles and produces a mathematically exact silhouette at any zoom.
This is the technique PyMOL, VMD and ChimeraX use, for the same reason.

**What it costs:** the shaders own the depth buffer, so anything that mixes with
them has to agree about depth. That is why text and the HUD are drawn by
**sibling Qt widgets** rather than by `QPainter` inside `paintGL` — the two
fight over GL state and the text simply did not appear.

---

## 4. Every number a calculation uses is registered

Not a constant in a function default: an entry in `resources/parameters.json`
with a unit, bounds, a kind and a citation, read through `_P.value("key")` at
**call time**.

**The constraint:** a constant written into a default is invisible. It cannot be
listed, shown to a user, traced to a paper, or swept. Several numbers this
project has had to correct were invisible in exactly that way.

**What forced it:** Round 49 found 26 registered parameters that no code read —
the dialog offered them, an override was recorded, reports carried the
non-default banner, and the computed number did not move. The reverse of the
same problem. `parameter_audit.py` now fails the build on an unregistered
literal, and `provenance_chain.py` reports parameters nothing reads.

**Resolution at call time, not import time**, so an override takes effect on the
next call rather than requiring a restart. The cost is a dictionary lookup per
call, which is why hot loops read the value once into a local first.

---

## 5. A checking instrument is a measuring instrument

Every cross-check, re-derivation, audit and probe is registered in
`tests/test_calibration.py` against the known-answer case that calibrates it,
and a test fails if one is added without one.

**The constraint — and this is the most expensive lesson in the project.** Six
times, an *alternative* route written to check the main one was itself wrong,
and it returned a plausible number rather than an error, so the disagreement
looked like a finding. A spheroid fitter that would have reported 89% model
error. A document checker that could not read the Unicode minus its own
documents use. A parameter probe whose "no effect" came from coordinates too
diffuse to form a single alpha sphere. A dead-code detector that would have
deleted the CLI.

An uncalibrated checker is worse than no checker, because it manufactures
findings. The rule is therefore: before a checker is believed against real data,
run it on a case whose answer is known independently — and the calibration must
be able to *fail*. If no input makes it say "no", it is not a calibration.

---

## 6. Results carry the state that produced them

A report records the structure, the parameter set and the commit. The result
window records them **at compute time**, not at display time.

**What forced it:** the result window is non-modal. A user can compute a pore
profile, open the parameters dialog, change `pore.step`, and still have the old
window on screen. Reading the registry when the window draws would label the old
numbers with the new settings. `verify_claims` goes further and **refuses to
run** against a modified registry, because the documented numbers were produced
at the defaults and any comparison would report drift the user caused.

The rule travels one window further. The exploration window opened from a
result is **handed** the result and its stamp rather than reading either again,
and its charts are built from that dict rather than from a second run — two
windows side by side disagreeing about which run the numbers came from is worse
than one window with no numbers at all. Its simulations pass slider values per
call and write nothing to the registry.

**What forced that:** an override survives the window. A slider that set one
would leave every later report carrying the amber non-default banner and
`verify_claims` refusing to run — correctly, and for a reason the user could
not possibly connect to a curve they moved ten minutes earlier.

---

## 7. Downloads are never committed; authored content is

`ref/` and `data/` are git-ignored and regenerable by `python -m piezo1.io.fetch`.
`piezo1/resources/*.json` is committed.

**The distinction:** a download is somebody else's copyrighted artefact and can
be re-fetched. A curated annotation — 68 variants with every wild-type residue
verified against the reference sequence, 17 domains with provenance and a
confidence — is authored content with no other source.

**What forced the verification:** the size floor on a download is necessary and
not sufficient. Round 60 found an Ensembl endpoint returning HTML; Round 65
found two 127-byte error pages *stored as structures*, where every later step
treated them as data. Downloads now declare what they should be and are checked
before anything is written — because a file already on disk is served from cache
without being re-checked, so one bad write persists indefinitely.

---

## 8. Documents that state numbers are generated or guarded

The parameter table, the bibliography and the notebooks are **generated** behind
gates: a citation that does not resolve stops the parameter build, a title that
does not match Europe PMC stops the reference build, and a notebook cell that
raises stops the notebook build. `CONCLUSION.md` and the README's closing
summary are **guarded**: every number in them must come from the claims
registry, the validation record or the published-interval table.

**The constraint:** prose drifts from code silently and reads as authoritative
while it does. A number in a document that nothing recomputes is a number nobody
is checking. That is also why the notebooks ship without stored outputs.

---

## 9. Files stay under 500 lines, and split at real seams

**The constraint** is reviewability, but the rule that matters is the second
half: a split is at a *seam*, not at a line count. `report_tags.py` separated
from `report.py` because each of its entries carries a caveat about what is
modelled rather than measured. `fusion_pose.py` separated from `fusion.py`
because `fusion` deliberately produces a region and `fusion_pose` deliberately
produces one arbitrary draw from it — the opposite claim, and worth its own
module docstring saying so.

A split made only to satisfy the counter produces two files nobody can name.

---

## 10. A picture presses the control the user would press

Nothing that draws on the 3-D view has two ways in. An exhibit in the
exploration window that offers to draw the pore surface **presses the View menu
entry**, not the controller behind it; the drawn overlays read the analysis
object the panel already computed rather than re-running it.

**What forced it:** the constraint is that a picture is more persuasive than
the number it came from, so two sources of truth for what is drawn is two
answers to "what am I looking at". A controller called directly would leave the
menu entry unticked while its overlay was on screen — and the cost of the rule
is real: the table of which control each exhibit presses can drift when a menu
entry is renamed, so a test resolves every one of them against a real window.

---

## 11. What this architecture is not for

It is **not a general-purpose viewer**. PyMOL and ChimeraX exist and are better
at that; this is PIEZO-specific and opinionated, and that is what lets it ship
a curated structure registry, a domain palette and a variant set rather than
asking the user to supply them.

It is **not an all-atom MD engine**. The interactive dynamics are coarse-grained
by design, because that is what makes a 2,500-residue trimer responsive. OpenMM
is available for offline refinement; `analysis/external_md.py` measured whether
deposited MD could contribute and found 1 of 21 entries in MemProtMD, resolving
none of the three blade basic clusters.

It is **not a clinical tool**, and the variant-effect prediction it was built
for does not work — five pre-registered tests, five nulls, and a measured
demonstration that the data which would settle it cannot be assembled. See
[`CONCLUSION.md`](CONCLUSION.md).
