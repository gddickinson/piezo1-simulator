# SESSION LOG

Running record of what was done and — more importantly — *why*. Newest first.

---

## Round 88 — The tag travels with the morph, and is re-solved rather than carried

Asked for directly: the morph did not include the HaloTags. It could not — the
fusion is drawn as its own batches from the loaded structure, so flattening the
channel left three tags hanging in space beside a dome that had moved out from
under them. They are attached to a C-alpha of the channel; that is the whole
point of the model.

### Carrying it rigidly would have been the obvious thing, and wrong

The tag's position is not a pose. It is the centroid of the region the tag
centre can occupy without clashing, so it depends on the channel's shape and not
only on where the anchor is. Flattening occludes a different part of the
tether's reach:

    accessible volume   7WLT -> 7WLU   242 -> 177 nm3
                        7WLT -> 11ZC   242 -> 152
                        8IXN -> 11ZC   261 -> 154

— a third of it gone, on every pair offered. The anchor-to-centre offset changes
with it, so a tag translated with its anchor lands about **7 A** out at the far
end. So the model is **re-solved on every frame**, and a test measures that
offset change rather than asserting the choice.

### Which needed the envelope to be affordable

One `query_ball_point` — 87k grid points against 32k atoms — is 90% of the cost
of a fusion model. `workers=-1` takes it from **456 ms to 75 ms**, and a whole
41-frame path from 25 seconds to about six. It is the same query threaded, so
`test_performance` asserts the kept set is **identical** on four entries, with
the guard that the query rejects something on each — "identical" is free if
nothing is being rejected.

Nothing was added to playback: the models are solved once when the path is
built, and the drawing costs 0 ms a frame. Scrubbing already costs 783 ms a
frame, all of it `MolecularView.update_coords` rebuilding the cartoon, which is
untouched and unrelated.

### The number that must not be quoted

At the far end the tag-to-pore-exit distance reads **3.92 nm** where the
deposited 7WLU's own model gives **3.59** — 8% apart, in the quantity the
calcium nanodomain work depends on. The path ends on 7WLU's C-alpha positions,
but the pore exit is *whichever atom reaches furthest down the axis*, and at
that frame the atoms are 7WLT's: its side chains, and its residues 7WLU does not
resolve. So the end entry's own model is measured at build time and both numbers
go on the status line with the instruction to quote the deposited one.

This is the same lesson as Round 87 one level up. A morph endpoint is exact
where the interpolation determines it and inherited everywhere else, and every
derived quantity has to be checked for which of the two it is.

### Interpolating the tag would have been cheaper and could put it inside

`show_frame` takes the **nearest solved frame** rather than blending two. Every
stored centre is a position the envelope admitted; a point between two of them
is not, and on a path where the channel closes around the tether it can lie
inside the channel. At 41 frames the step is 2.5% of the path.

### The speedup exposed a defect that had nothing to do with any of this

`test_ui_controls` started failing: a spliced model's `PART PREDICTED` status
had been replaced by `bottleneck 0.93 Å · non-conductive`. Stashing the change
put it back, and stashing *only* `fusion.py` narrowed it to the threading —
a pure timing shift, no logic involved.

What it exposed is real. The analyses run off the GUI thread and a worker cannot
be interrupted, so a pore profile launched on one entry can finish after the
user has loaded another. `load_structure` calls `analysis.reset()`, which clears
the *stored* result — and the run already in flight then landed on top of the
new entry. Not merely recorded: `_on_finished` calls `pore_surface.refresh()`
and `nanodomain.refresh()`, so one entry's bottleneck and calcium source were
**drawn inside another's lumen**. That is the thing the clearing in
`load_structure` says in its own comment that it is there to prevent.

Each run is now stamped with the structure object it was launched for, and a
result whose stamp no longer matches is discarded — reported on the Analysis
panel, not the status line, which belongs to the load that invalidated it.
Hijacking it there is what would have replaced `PART PREDICTED`, the one place a
user is told the model is partly a prediction.

It is in `ui/hazards.py` as the eleventh entry, with the positive control the
register promises: the situation constructed, the guard watched firing. A hazard
that only shows when the scheduler is slow enough is not guarded.

### One expression for the drawn coordinates

Writing the tests found that the coordinates the tag was solved on and the ones
it was drawn on differed by ~0.005 A — `float32 + float32` against `float64`
then rounded. Both now go through `_coords_at`. Small, but it is the same defect
as Round 87 in miniature: two routes to the same coordinates.

---

## Round 87 — The morph never reached the structure it said it reached

Reported from use, not by a guard: *the start and end of the morph do not look
like the structures they are supposed to be*. They did not, and the reason is
worth keeping, because the interpolation itself was never wrong.

`tests/test_morph.py` asserts the trajectory hits both endpoints to 1e-9, and it
does. The controller then draws it — and a morph is drawn as a **displacement
field**, `frame - frames[0]` added to the coordinates already on screen, which
is what carries the side chains, lipids and everything else the C-alpha path
does not describe. A displacement means nothing outside the frame it was
measured in, and the endpoints were read from disk while the viewport shows the
canonical frame.

    trajectory centroid (7WLT deposited)   197.6, 197.6, 257.6
    displayed centroid (canonical)           0.0,   0.0,  -1.5
    rotation between the two frames         180 degrees

So every atom was pushed the wrong way. Measured against the deposited flat
structure as a *shape*, so no assumption about frame enters the comparison:

    drawn endpoint vs 7WLU, before   35.98 A
    drawn endpoint vs 7WLU, after     0.000 A
    the change being interpolated    19.70 A

**Worse than the change itself** — the far end of the slider was not an
approximation of the flat state, it was the curved state driven backwards. The
camera was framed on the trajectory too, 381 A from the model, which is why the
*start* looked wrong as well while being numerically exact.

### Fixed by removing the possibility, not by rotating the field

The start structure now goes on screen first and the path is built from the
coordinates that are displayed. Multiplying the field by the frame rotation
would have worked and would have left the same trap for the next person: two
sources of coordinates that have to be kept in step. There is now one.

### The second fault, which the first was hiding

Atoms at a residue outside the shared basis were given a site by **residue
number**. 7WLT's lipids are numbered 2601-2609 and the last shared residue is
2546, so all 1,407 of them were tied to the C-terminal C-alpha and travelled
with the CTD tip — one of the largest displacements in the motion.

    lipid to the site it was tied to   median 64.8 A, max 98.9 A
    lipid to the nearest site          median  6.2 A

They now take the nearest site in space, which is the helix they are actually
against. A residue *in* the basis still takes its own site, so whole residues
stay together; a test checks both halves, because "nearest" applied to
everything would have quietly dissolved the residues.

### A morph is a stored array, and nothing was clearing it

`load_structure` clears the dome, the contacts, the pore, the pockets, the route
and the calcium field. The morph was not on that list, and it is the only one
holding a base coordinate array of its own — so after loading another entry the
slider stayed live and would have pushed 7WLT's motion onto a different
structure. Added to the list.

### What the last frame still is not

It ends on the target's C-alpha positions; every other atom is carried with its
residue, so it is not the deposited entry and the status line now says so. The
**modal** method does not even do that: it is confined to the elastic-network
subspace, captures 95% of the change on this pair and deliberately stops 6 A
short. One sentence covering both would have been false for one of them, and the
modal case is the one a user would report as this same bug.

---

## Round 85 — The review, and the hole it found

Block R's standing question was whether the discipline held now that results
could come back *positive*, after seventy rounds of establishing what this
project cannot do.

For positives it held, and visibly. Every new instrument in the block was
calibrated on a known answer before it was believed, and every calibration
caught something: an electrostatics constant 10^10 too large, a helix detector
passing 41% of a random walk's windows, a numbering test that read mouse PIEZO2
as human PIEZO1, a gate-radius claim confounded with resolution.

**For negatives it did not, and that is the review.**

Round 84d reported that the lateral conduction pathway "does not separate open
from closed". It recorded that as an honest negative — the sort of thing this
project is proud of — and **pinned it with a test**. It was wrong. It was an
artefact of evaluating both halves of the wetting verdict on a truncated
profile, which collapses the Rao score to zero on every entry. Round 84f found
it, and only because the question was asked directly.

Round 84c has the same shape a step earlier: "17 of 19 refused" was reported as
an honest outcome while the stated reason for it was wrong.

The pattern is worth stating plainly, because it is not obvious from inside.
Every guard here — calibrate the checker, suspect the checker first, report a
null as a null, a checker that cannot say no is not a calibration — is aimed at
**not over-claiming**. A guard that asks "are you sure?" only ever fires at
confidence. A negative from a broken instrument sails through all of them,
because nothing in this repository interrogates a "no".

    A null needs a positive control — an input on which the instrument MUST
    return "yes" — exactly as a checking instrument needs an input on which it
    must say "no".

The pre-registered statistical nulls already satisfy this without calling it
that: a power curve is a positive control. Round 82's planted fluctuation is
one. The engineering negatives of 84c and 84d had neither, and those are the two
that turned out to be wrong.

Applied rather than noted: the test that pinned 84d's negative now keeps the
broken composition as a *demonstration* of the defect, beside the composed
verdict that refuses the same entries, and the separation is checked against Liu
et al.'s Figure 5D ordering rather than against a stored copy of our own
numbers. The rule is in `docs/METHODS_NOTE.md`, which is the file meant to be
useful to somebody else's project.

### The roadmap is empty, and that had to be made a legitimate state

`test_the_open_items_are_all_in_the_live_roadmap` required at least one open
item, on the assumption there is always work. Block R finished it. An empty list
is meaningful — but indistinguishable from a file that lost its contents, so the
guard now accepts an empty roadmap **only if the file says so**.

### One observation about where the corrections came from

Three of this block's corrections came from being asked a question rather than
from a guard firing: why does 8IXO not conduct, why does the lateral route not
separate, does this need further work. Each time the honest answer required
measuring something nobody had measured. That is not a substitute for a guard.
It is evidence that the most valuable input to this project is still somebody
asking why a number is what it is.

---

## Round 86 — Two entries read in the wrong numbering, and now read right

Round 83 built an instrument that identifies which protein and which numbering a
deposited file is in, and it found two register errors and **reported the shift
that would repair each**. Nothing applied it. So five entries were read at
residue numbers that point at the wrong residue, and everything this project
reads by number — transmembrane helices, domain boundaries, variants — was wrong
inside the affected range.

    6LQI                      0.447 -> 1.000   +24 over 765 residues
    8ZU3, 8YFC, 9VMX          0.932 -> 1.000   +22 over 767-857
    8YFG                      0.931 -> 0.999   the same, and correctly short

8YFG stopping at 0.999 is the detail worth keeping: it carries R2456H, a genuine
substitution. A numbering fix that reached 1.000 there would have absorbed a
real residue change, and the test requires it not to.

### A renumberer is a rewriter, so the null is the load-bearing test

Applied to a file that is already right it would silently corrupt it. **8YEZ
resolves the same 767-857 region and is numbered correctly throughout**, and it
has to come back unchanged — by identity, not by equality, so a future edit
cannot start returning a copy. 7WLT, 3JAC and 6B3R likewise.

### Chance agreement splits a register error

`mismatch_blocks` finds *runs* of disagreement, and a residue that happens to
agree at the wrong numbering ends a run. On 8ZU3 that split one 91-residue error
into three blocks — 772-787, 789-834, 839-857 — and correcting only those gave
0.969 where a uniform read of the whole stretch gives 1.000.

So spans sharing a shift are merged across small gaps and then grown outward one
residue at a time, keeping each step **only if the corrected identity does not
fall**. That makes the extension a measurement rather than a guess. It grows past
the data — 8ZU3's span reaches 713-914 — which looked wrong until I checked:
that entry resolves *nothing* between 713-766 or 858-914, so the wider span
moves exactly the 91 residues that exist. The report says the spans are in file
numbering and may cover unresolved gaps, because the next reader will wonder the
same thing.

### Which published numbers move, which was the round's actual question

**The dome moves and the pore does not.**

    8ZU3, 8YFC, 9VMX   12.50 -> 11.91 nm      TM helices 25 -> 26
    8YFG               11.12 -> 10.77 nm      TM helices 25 -> 26
    6LQI                9.35 -> 11.03 nm      TM helices 23 -> 26
    pore bottleneck    unchanged to the Angstrom on all five

The dome is fitted to the **annotated** transmembrane helices, so a register
error inside the blade puts the wrong atoms in the fit — and the correction
recovers helices that were being missed entirely, which is why the count rises.
The pore profile reads no residue number at all and is untouched. That contrast
is the useful statement: it says exactly which quantities the numbering reaches.

**No frozen claim uses any affected entry**, so nothing is superseded — and
there is a test that fails if one ever starts, rather than a note saying it was
checked once.

---

## Round 84 — Nothing computed could leave the application

The oldest open item in Block R, and the smallest: conservation, mechanical
coupling, perturbation response, mode displacement and the wetting score are all
per-residue scalars, and none of them could be handed to another viewer.
`to_pdb` writes coordinates. The standard route — put the scalar in the
**B-factor column** and open the file in PyMOL or ChimeraX — is about twenty
lines and did not exist, so every number this project computes was trapped in it
or in a JSON blob.

The twenty lines were the easy part. The rest is what stops the file lying
somewhere none of this project's guards can reach.

### Unscored is not zero, and the display array is not the data

A residue the analysis could not score must not arrive in PyMOL looking like one
that scored zero. Unscored atoms go out with **occupancy 0.00** and a B-factor
at the floor of the column, so either column selects them and the count is in
the header.

The GUI export needed a second version of that same thought. `view.values` — the
array actually painted on the model — has unmeasured residues **filled to the
map floor** so the colour ramp behaves. Exporting that would have shipped a file
where an unscored residue is indistinguishable from a low-scoring one, with no
flag at all. The controller now remembers which scalar is on screen and the
export writes the **raw residue map** instead.

### The field is six characters, and my constant was wrong

The first sentinel was `-1.0`, chosen because "no real value of anything this
project computes is negative". The wetting energies are negative; the guard
caught it on the first real analysis, which is the guard working.

The replacement was `-999.99`, the floor of what I believed the column held.
That is **seven characters in a six-character field**. It overflowed, shifted
every column after it, and produced a file that still parses and is wrong —
which is precisely the class of error this module exists to prevent, committed
by the module itself. `F6.2` actually holds `-99.99` to `999.99`.

So the limit is no longer remembered. `fits_column` formats the value and
measures the result, and the test checks the **element symbol** still reads out
of columns 77-78 — the thing a B-factor overflow destroys, and the only
assertion that would have caught it.

### Reachable, which was the complaint

`piezo1.cli export 8YEZ --scalar wetting` and **File → Export coloured
structure**, which writes whatever is currently painted. On 8YEZ the wetting
energies cover 1,125 of 31,599 atoms — a partial scalar, which is exactly the
case the occupancy flag exists for.

---

## Round 84f — Which profile decides which half

Round 84d made the conduction pathway selectable and reported that the lateral
route "does not separate open from closed": 8IXO at 53.8 pS, but several curved
entries at 6-12 pS. That was pinned as an honest negative. Asked why, it turned
out not to be a property of the channel at all.

### The lateral verdict was close to noise

The wetting verdict has two halves, and on the full axis they do very different
jobs. **Sterics refuses everything** — including 8IXO and 11ZC — because the cap
apex and the cytoplasmic neck are shut in every deposited structure; that is a
constant "no" that only looked decisive. **Hydrophobicity actually separates**:
8IXO at 0.31 and 11ZC at 0.00 sit below the 0.55 cutoff while 15 of 16 curved
entries are above it.

Round 84d evaluated both halves on whichever profile the pathway produced.
Truncating fixes the sterics and *guts the score*, because the Rao score is a
sum over lining residues and its contributions come almost entirely from narrow
slices — which is exactly what the truncation removes:

    7WLT   1.35 -> 0.13
    6B3R   2.05 -> 0.22
    entries above the cutoff:  13 of 18  ->  0 of 18

With the chemistry disabled the verdict rested entirely on whether a residual
radius cleared 1.5 A. Measured, that residual was **Y2335, the cap gate, in 14
of 18 entries, and the transmembrane gate in none**. So the Round 84d lateral
conductances were not a statement about gating; they are superseded.

### The rule, and the order it was arrived at

    Hydrophobic gate from the COMPLETE AXIAL PROFILE — Rao et al.'s cutoff is
    a sum over a whole channel's lining and is not comparable between profiles
    of different length, so applying it to a truncated one is invalid on its
    own terms.

    Steric occlusion from the CONDUCTION PATH — whether an ion fits is a
    property of the route it takes.

`analysis/conduction.py`. On the axial pathway the two profiles are the same
object, so the answer is bit-identical; that is checked on every downloaded
entry one by one, because a single entry would not exercise a case where they
differ.

**How it was reached is recorded in the module, because it matters.** Three
compositions were tried and the third separated the states. That is the wrong
order. The rule is adopted on the calibration argument, which would hold
whichever way the numbers came out, and the agreement is reported as a check
rather than as the derivation. I also tested the obvious alternative first —
cutting at the cap gate rather than the apex, since the paper says ions enter
*through* the cap gates — and it made things worse: 17 of 18 pass.

### The check

Under the corrected rule, on the lateral route:

    curved (closed)        refused, 15 of 16     their Fig 5D: ~0 Na+/us
    flattened (7WLU)       9.2 pS                             ~10
    intermediate (8IXO)    40.1 pS                            ~20

Same ordering, including the part that is easy to get wrong: the *flattened*
state conducts, in their data and in ours, because its transmembrane gate is
dilated while its cap gate is shut. The one curved entry that slips through is
3JAC — 4.8 A with 346 unnamed residues, so nothing in it is narrow and its
score is 0.06. A coverage artefact rather than a state, and named in a test so
that fixing it cannot look like an improvement in the model.

### One mistake on the way, pinned

The first version handed `solve_pnp` the raw full-axis wetting object, which
re-imposes the axial steric block on a route chosen to avoid it. The sweep
silently returned **0.0 pS for both 8IXO and 7WLU** while the composed verdict
said they conduct — a disagreement between two parts of the same answer, which
is the kind that survives a green test run. It now gets the composed verdict,
and there is a test that would fail if it stopped.

---

## Round 84e — Showing one part of it

The request: display selected components of the assembly with the important
residues and the relevant backbone, several presets including the pore region
of Liu et al.'s Figure 2E, plus opacity control on the drawn pore and more
colouring options including hydrophobicity.

### Components come from the annotation, not from ranges in the viewer

`structure/components.py` defines ten parts entirely as *annotation ids* —
domains from `domains.json` for the backbone, residue groups from
`functional_residues.json` for what to pick out. Nothing in the viewer says
"the cap is 2214-2457". That matters because those files have had corrections
(the cuff added in Round 84, the cap gate and spring linker added in 84d) and a
component inherits every one of them. It also means a component is exactly the
residues the analyses read, rather than a second opinion about where the cap is.

The one worth having is `pore_module`: outer helix, cap, spring linker, inner
helix and CTD, with all four gates in ball-and-stick. On 8IXO that is 8,343 of
32,112 atoms, and it is the only view in which the paper's argument is legible.

A highlighted residue outside its component's domains is **added** to the
backbone rather than left floating in space — the PIP2 lysines sit just past
THU9's end, and a residue picked out as important and then not drawn is the
worst of both.

### It hides; it does not subset

The structure object is untouched. The pore profile, the dome and the modes all
still run on the whole trimer, the status line says so on every switch, and a
test measures the dome either side of a selection and requires it to be
bit-identical. This project's rule is that what is drawn never decides what is
computed, and a component selector is the most tempting possible place to break
it.

### The defect, found in pixels

Setting the residue filter and calling `rebuild()` filtered every atom
representation and **left the cartoon drawing the whole chain**. `traces` is
built once in `__post_init__`; `rebuild()` re-uploads batches from the traces it
already has. At the default style the ribbon is most of what is on screen, so
hiding 97% of the atoms changed the picture by about a tenth — which looks like
a subtle bug and is a completely wrong picture.

Nothing about the batch contents would have caught it. `set_visible_residues`
now rebuilds the traces too, and the test renders the transmembrane gate and
requires it to draw under 30% of the whole assembly's pixels.

### Opacity needed the shader

The sphere impostor wrote `vec4(color, 1.0)` and had no alpha at all — which is
why `geometry_builders` notes that a translucent iso-surface has to be a mesh.
It now takes a per-batch `u_alpha`, and a batch below full opacity declares
itself `transparent` so `Scene.render` moves it into the blended pass with
depth writes off. Drawn in the opaque pass it would write depth and hide the
lining it was meant to reveal.

The test asserts brightness falls while the **pixel count does not change** —
a guard that only counted lit pixels would pass on a batch that had simply
stopped drawing.

### Hydropathy, on a fixed scale

`ColorBy.HYDROPHOBICITY`, Kyte-Doolittle, blue polar to orange apolar, fixed at
+-4.5 rather than auto-ranged. Same argument as the electrostatic colouring: an
auto-ranged hydropathy map paints a uniformly polar loop in full orange, and two
structures coloured that way cannot be compared with each other or with a
published panel. A test checks alanine does not change colour when nothing else
is on screen. Anything the scale does not name is mid-grey, so "not scored"
cannot be read as "neutral".

---

## Round 84d — The paper that says the current does not go down the axis

Follows directly from 84c. Having established that the ion animation refuses
every structure on constrictions Liu et al. 2025 report as bypassed, the ask was
to reproduce their figures and their Figure 5 permeation simulation, and to make
the scientific choices selectable rather than baked in.

### Why 8IXO does not conduct, completely

Round 84c got half of it. Truncating the path at the cytoplasmic neck does not
rescue 8IXO — it stays blocked at 1.03 A. Looking at every sub-2 A slice between
their own pore endpoints, R2295 and E2537, every remaining constriction is
**R2295 itself and its immediate neighbours**, on both 8IXO and 7WLT. Which is
what the paper says: *"The top portion of the cap remains closed above the
residue R2295 position among all the structures, consistent with mutagenesis
studies suggesting that ions might not permeate through the top of the cap."*

So PIEZO1's axial pore is closed at **both** ends, deliberately. Ions enter the
cap vestibule through three lateral cap gates and leave the inner vestibule
through intracellular lateral portals; the axial segment is only the middle.
A bulk-to-bulk 1-D model must traverse both closed ends, and that — not the
radii — is why it refuses a structure whose gate has demonstrably opened.

### The pathway is now a choice

`physics/conduction_path.py` offers `axial` (default, returns the *same profile
object*, so nothing recorded can drift), `lateral_entry`, `lateral_exit` and
`lateral`. On the lateral route **8IXO conducts at 53.8 pS** where the axial
model refuses it.

And the honest half, pinned as a test: **it does not separate open from closed.**
7WLT, 6B3R, 8IMZ, 6BPZ, 11YE and 9VED also conduct once the closed ends are
excluded, at 6-12 pS against 8IXO's 53.8. The right ordering, roughly fivefold,
and not the open/shut contrast their simulations show. Opening the ends is
necessary and not sufficient.

The option deliberately does **not** model the portal. The truncated end slice
becomes the mouth and the Hall access resistance is computed from the pore's
radius there rather than the portal's, which this project does not measure — so
a lateral current is an upper bound and the caveat says so wherever it appears.

### Seven published distances reproduce

All four of their states are deposited and all four are in the catalogue, which
is what made this possible: 7WLT, 8IXN, 8IXO, 7WLU.

| | published | ours |
|---|---|---|
| R2295-E2537 pore axis | 110 -> 100 A | 109.5 -> 96.2 |
| V2476 side-chain diagonal | 7 -> 14 A | 7.7 -> 14.2 |
| A2328-P2382 cap loops | 4.3 -> 16.2 A | 4.8 -> 16.1 |
| D2326-E2383 cap loops | 4.8 -> 12.8 A | 5.7 -> 11.4 |
| Y2464 spring, compressed | 17 A | 16.6 |

The cavity volumes reproduce the *direction* of every change — CV, EV and MV
grow into the intermediate state, IV does not — and not the values, because
ours is a solid of revolution and theirs is not. Filed as an analogue for that
reason.

### One panel disagrees, and it is recorded rather than adjusted

Their Figure 6 curvature radii are ~10-12, 14, 32 and 117 nm across the four
states. Ours are 9.7, 11.2, 16.5 and 18.4. The fitter agrees where it was
calibrated — 9.7 nm against Guo & MacKinnon's 10.2 on the curved state — and
saturates where it was not. Fitting a sphere to a nearly flat surface is
ill-conditioned and under-estimating a large radius is exactly how that fails.
The test pins it under 40 nm and says that if it stops saturating, the panel can
be promoted from analogue.

### Figure 5, and the line it cannot cross

`analysis/liu2025_permeation.py` sweeps their four voltages. On 8IXO the slope
is **40.1 pS against their 20 pS** — twice, consistent with the 1.4x this model
already overestimates by on 11ZC.

Their **5C is refused, not approximated**. It counts Na+ that accessed each
cavity, including the ones that entered and turned back; a one-dimensional
steady state carries the same flux through every slice by construction and
cannot distinguish the cavities at all. A number there would be an artefact of
the discretisation.

`physics/martini.py` prepares the coarse-grained system and does not run it.
The boundary is enforced rather than asserted: everything `prepare` produces is
an input, `MartiniRun` is constructed in exactly one place (a test counts the
occurrences in the source), and `load_results` raises when there is no
trajectory instead of estimating one. When the run happens, its counts convert
to picoamperes through the same arithmetic the continuum current uses, so the
two land on one axis.

### Four defects found on the way, three of them latent for a long time

**`_accessible_area` returns exactly zero** for Cl- (1.81 A) in any pore between
1.38 and 1.81 A, and a zero-area row makes the PNP system singular — LinAlgError,
not an answer. Physically a species the pore excludes has zero flux, since the
slices are in series; now stated. 11ZC unchanged at 40.12 pS.

**Chord conductance at 0 V was 0/0.** A solver flux of 1e-20 A made
`abs(0/current)` exactly zero, the pore dropped out of the series, and the
reported conductance became 1/R_access — **1,586 pS on 8IXO, forty times its own
value at -0.1 V**. Found only because their protocol starts at 0 V and no
earlier caller used one.

**`pore_charge._curated_map` selected by `category == "pore"`.** Adding the
cap-gate and spring-linker residues as curated annotation — all genuinely pore
elements, all correctly categorised — took the curated charge set from 6 to 12
and **flipped the measured selectivity from cation- to anion-selective**. An
annotation edit must not be able to redefine a recorded measurement, so the
group list is written down explicitly now.

**`flux.timebase_for_structure` passed `default_species()` as `wetting`**, the
second positional argument. Harmless only because the caller had already gated
on the same verdict two lines above, so the solver's own blocking check saw a
list, found no `available` attribute and skipped.

### The citation gate caught a wrong PMID, then died instead of saying so

I cited Liu et al. as PMID **39674176** throughout the new modules. The
bibliography's title-verification gate rejected it: that PMID resolves to a
paper about IgG accumulation in adipose tissue. The right one is **39719701**,
found by DOI. Corrected in nine files.

Which is the gate working. What was not working is what happened next —
`build_references.py` referenced a `MANUAL` table on the rejection path that
has not existed for some time, so instead of reporting the rejection it raised
`NameError`. Latent precisely because the path only runs when a citation is
wrong. Fixed; 80/80 references now resolve.

### Two architecture violations, caught by the guard that exists for them

Both mine, both real. `physics/martini.py` and `analysis/liu2025_permeation.py`
imported `render.flux` for the elementary charge, pointing the dependency arrow
backwards; and `physics/conduction_path.py` and `physics/martini.py` imported
`analysis.pore_regions` for the numbering, which crosses the science layers the
wrong way.

Fixed by moving rather than by relaxing the test. `ELEMENTARY_CHARGE` and
`ion_rate` are now `physics/charge.py`, re-exported from `render.flux` so every
existing import still works. And "which numbering is this PIEZO1 entry in" is
now `core.numbering_check.piezo1_numbering`, which is where it belonged: two
different layers need it, and the arrow between them only points one way.
`menus.py` crossed 500 lines on the way and split at the same seam.

### And two claims of my own that the calibrations caught

`gate_numbering` decided human-versus-mouse from the three gate residues' names
and read **mouse PIEZO2 as human PIEZO1** — Ile, Val and Phe are not three
distinguishing observations. It goes through `identify_numbering` now.

And I wrote that 8IXO has the widest gate in the catalogue. It does not: 7WLU
(4.67 A) and 3JAC (4.34 A) are wider and are the two worst-resolved entries in
the set. Gate radius is confounded with resolution and cannot carry the
dilation claim; the side-chain diagonal can, and does. The corrected fact is a
test that fails if the confound ever disappears.

---

## Round 84c — Why nothing conducted, and which constriction actually refused it

The report was that the ion flux animation shows every structure as
non-conducting, flat or curved. Two separate things were true, and the second is
the more interesting one.

### The animation had never drawn an ion

`timebase_for_structure` is right about 11ZC — 2.4 pA, gated by the wetting
verdict, and `test_ion_flux` has pinned that since Round 33. But the *controller*
died on its first frame. Twelve particles enter per displayed second, so at
60 fps the first four frames hold none; `SphereBatch.upload` passed a zero-length
payload to `ctx.buffer`, which refuses one — "the buffer cannot be empty" — and
`ViewportWidget._on_tick` answers an exception from a frame callback by
**unregistering the animation**. So the stream stopped before the first ion
existed, the status line was overwritten with `animation stopped: the buffer
cannot be empty`, and the two conducting entries in the catalogue were
indistinguishable from the seventeen that are refused.

Every existing test of that controller passed throughout. They checked the time
base, the gating and the colour; none drove `_step` against a real scene. The
new ones count lit pixels, which is the same lesson the invisible cylinders
taught in Round 79 and it did not generalise on its own.

Two smaller things came out with it. The stream was drawn on the straight C3
axis, but the pore is measured as a probe whose centre is *leashed* to within
8 Å of that axis rather than pinned to it — on 11ZC the fitted centre sits a
median 0.56 Å off it and at 11 of 125 heights the axis line falls **outside**
the sphere fitted there, so the ions were crossing the wall of the pore they
were meant to be in. And the direction was whatever sign `detect_c3_axis`
happened to return, which fixes a line and not a sign; it now runs towards the
end `cytosolic_end` measures, the same function `selectivity` uses.

### The refusal was true and said the wrong thing

Seventeen of nineteen deposited PIEZO1 entries are refused — and the two that
are not are the two worst-resolved models in the set, 11ZC at 6.0 A with no side
chains at all and 3JAC at 4.8 A with 346 unnamed residues. The message was
`sterically occluded` with a bottleneck radius. That reads as *the gate is
shut*. It is not what was measured in a single one of them.

Locating the transmembrane gate from the curated `hydrophobic_gate` residues
rather than assuming a range, the narrowest point of the axial profile is **at
the gate in none** of the 18 entries where the gate can be located — below it at
the cytoplasmic constriction in 16, above it in the cap in 2. The gate itself
measures **2.4–4.7 Å** in every entry, at or above the 1.5 Å water radius the
steric test uses.

8IXO is the entry that makes this unambiguous, and it is why the paper the user
supplied mattered. It is Liu et al. 2025's intermediate-**open** S2472E
structure. Our profile gives its gate 3.52 Å; the Rao score clears the cutoff at
0.31; and on the paper's own measure — the V2476 side-chain diagonal — we get
**14.2 Å against 7.7 Å on the curved 7WLT**, reproducing their 7 → 14 Å. It is
refused on a **0.98 Å neck at E2537**, which is precisely the constriction that
paper reports as *remaining closed* in this structure, because "the lateral
portals rather than the constriction neck represent the major ion-permeation
routes".

So the model must pass through a constriction the channel goes around. That is a
stated limit of a one-dimensional axial conduction model, not a property of the
structures, and `analysis/pore_regions.py` now says so on the status line: the
gate's radius first, then the narrow point on each side of it. It reports three
regions rather than the global minimum because the two flanking constrictions are
within 0.02 Å on 8IXO and which of them wins flips with the frame the structure
was loaded in — a number that changes when you rotate the molecule is not one to
put in a sentence.

**No verdict changed.** A test asserts that `describe_bottleneck` leaves
`predict_wetting` bit-identical. Whether the conduction model should be allowed
to leave the pore laterally is a real question and a much larger one; it would
move Round 19's wetting verdict, Round 34's null and the frozen claims, and it
is not something to do as a side effect of fixing an animation.

### Two calibrations, and both caught me

The first version of `gate_numbering` decided human-versus-mouse itself, by
checking the three gate residues' names at each numbering's positions. It read
**mouse PIEZO2 (6KG7) as human PIEZO1** — Ile, Val and Phe are not three
distinguishing observations, and any transmembrane helix is full of them. The
job went to `identify_numbering`, which scores every residue against all six
references and gives 1.000 for the right one and under 0.25 for the rest. It
also correctly refuses 6LQI, which is in the Piezo1.1 isoform's numbering: its
"V2476" side chains sit 31 Å apart against 7.7 Å on 7WLT.

The second was a claim, not code. I wrote that 8IXO has the widest gate in the
catalogue. It does not — 7WLU (4.67 Å) and 3JAC (4.34 Å) are wider, and are the
two worst-resolved entries in the set at 6.81 and 4.8 Å. Gate radius is
confounded with resolution and cannot carry the dilation claim; the side-chain
diagonal can, and does. The corrected fact is now a test that fails if the
confound ever disappears.

### A mislabelling found on the way

`PoreSlice.lining` was `sorted(set(numbers))` and `lining_names` was
`sorted(set(names))` — two independently sorted sets, documented as parallel and
zipped by `predict_wetting`. So the residue the application *names* as the worst
dewetted one was paired with an unrelated name: 8YEZ reported GLU2510 where it is
PRO2510, and 11YE reported LEU2427 where it is ILE2296 — a different residue, not
just a different label, because `zip` truncates to the shorter tuple and dropped
residues out of the point list entirely.

Measured before fixing, across all 19 entries: scores move by at most 4.6%
(6B3R 1.96 → 2.05), **no hydrophobic-gate verdict flips**, and both frozen claims
(8YEZ 0.82, 11ZC 0.00) are unchanged. Eight of nineteen entries had the wrong
residue named. The fix is one `np.unique(..., return_index=True)` in `pore.py`,
and the test only counts slices that touch more than one kind of residue —
on a slice of all-leucines the old code was accidentally right.

---

## Round 84b — Reproducing the paper the whole project rests on

The dome model is Guo & MacKinnon 2017. The 10.2 nm radius in `parameters.json`,
the 120 nm² in `PUBLISHED_AREA_ESTIMATES`, the two-state Boltzmann in
`physics/dome.py` — all of it is that paper's Figure 7, cited and never
recomputed. The ask was to be able to replicate all parts of its figures from
our own structures. Doing that turned out to be less about drawing pictures than
about being honest, panel by panel, about which ones we can actually make.

### The registry is the deliverable, and its refusals are half of it

`analysis/guo2017.py` holds all **31 panels** as data: what each shows, whether
it reproduces, and — for the ones that do not — why not. **16 reproduce, 3 have
an analogue that is a different quantity, 12 need experimental data we do not
hold.**

The 12 could simply have been left out of a module called "replicate the
figures". They are in it because a tool that quietly covers the tractable parts
of a paper leaves a reader assuming the rest. An FSC curve needs two half maps.
A local-resolution map needs the map. Four panels are micrographs of
proteoliposomes. One needs P2X and ASIC coordinates, which are *deliberately*
absent: the structure catalogue, the numbering checks and the entity classifier
all assume a PIEZO, and admitting two unrelated channels to make one panel would
weaken every one of those guards. That is a real gap and it is recorded as one.

The three analogues are the dangerous entries, because they produce a picture
that can sit beside the original. A projection of an atomic model is **not** a
2D class average — no CTF, no solvent, and no detergent micelle, which is most
of what Figure 2b's envelope actually is. A screened-Coulomb surface is **not**
APBS. Both carry the caveat structurally rather than in a docstring, and a test
asserts that no path shows either without it.

### Figure 7 reproduces exactly, and that is a check on arithmetic, not on PIEZO1

Every number in Figure 7 and its supplement falls out of two lengths by
closed-form spherical-cap geometry — an 18.8 nm opening, a 6.2 nm depth,
397 nm² of surface, 277 nm² projected, 121 nm² released, 153 k_BT of bending,
42 k_BT of stabilisation. That is a good calibration precisely *because* it is
not a measurement: the answers are known independently of any code here, so an
implementation that reproduces them is doing cap geometry correctly, and one
that does not is broken in a way a real structure would have hidden.

Our own measurement of 6B3R's dome does **not** match it — 568 nm² of surface
against 397 — and `compare_with_measured` reports both and adjusts neither. The
idealised dome is a shape the authors chose to make the energetics tractable
and they say so; the gap between an idealisation and a measurement is a result.

Making Figure 7c quantitative — flattening at constant membrane area — surfaced
something the paper does not spell out. Complete flattening releases the *whole*
152.8 k_BT of bending energy as well as the whole 120.8 nm² of projected area.
With the paper's own (ΔG_prot + ΔG_bend) of 20–40 k_BT, that puts ΔG_prot near
**+180 k_BT**.

### Three instruments, and what each calibration caught

**Electrostatics.** `e²/4πε₀` was written in Joule-Angstrom and then converted
from metres a second time. The Bjerrum length came out at 7×10¹⁰ Å, the Debye
length underflowed to zero, every potential was exactly zero — and the
truncation-error check reported a flawless 0.000% because it was comparing zero
with zero. Nothing raised, and a "coloured surface" would have been uniformly
one colour and might well have been shipped. The single-point-charge case is
what said no. The cutoff's documented error was wrong too: I wrote "under 0.2%"
from the analytic falloff, which bounds *one* charge and not 804; measured on
6B3R's own surface it is 2.2%.

**Planarity.** My first control was one protomer's points replicated three
times — same point count as the trimer, same shape as a protomer. It agreed
with the protomer RMSD to every digit on every structure, because replication
leaves a least-squares plane exactly unchanged. A tautology dressed as evidence.
Replaced with a decomposition that can fail, and checked for closure.

**The helix detector.** Rise and radius alone passed **41%** of the windows of a
synthetic random walk — a walk with a fixed step length looks locally like a
helix on both. The turn criterion, applied to the worst step rather than the
mean, takes that under 5% while admitting every window of an ideal helix. And
using the window's own principal axis biased the estimator to 1.63 Å rise
against the textbook 1.50; sub-window centroids fixed it, which matters because
otherwise the tolerances would have been absorbing a systematic error as well
as a real spread.

### Coverage, again, in a new place

Measured on whatever each entry resolves, 6B3R's non-planarity is 17.2 Å and
6BPZ's is 4.7 Å. Two structures of the same protein apparently disagreeing about
whether it is curved. They do not: 6BPZ resolves 14 transmembrane helices and
6B3R resolves 26, and coverage-matched both give ~4.6 Å. The non-planarity lives
almost entirely in the distal blade.

This is exactly the trap `analysis/paralogue.py` was written after in Round 83,
in a completely different place, and I walked into it again — which is the
argument for `blade_dependence` measuring the split on **one** entry, where no
second entry's coverage can confound it.

### The 4-TM repeat, measured rather than cited

The paper extends six visible units to nine on the strength of hydropathy. That
inference is load-bearing for *this* project — it is why `domains.json` has nine
THUs and why the full-length model grafts a distal blade — so it is now measured
against a shuffled control. Supported in both mammalian PIEZOs (z = 4.5, 5.0),
**not** in PEZO-1 or dPIEZO.

The control has to be **register-maximised** like the statistic. Taking the best
of four registers and comparing against an unmaximised null moved the control
mean from 38 to 0 and would have manufactured about a standard deviation of
significance from nothing.

Separately: PIEZO1's transmembrane helices average +1.22 on the Kyte-Doolittle
scale, below the conventional +1.6 cut, so a threshold call recovers 6 of 38.
The threshold is left at the published value and the whole recall curve is
reported. Tuning it to this protein would have made the agreement a statement
about the tuning.

### The topology diagram, and the one thing it must never do

Asked for as a GUI feature: a monomer's topology in a membrane, with selectable
groups of TMs as in Figure 3b's boxes. `analysis/topology.py` builds it as data
and `ui/topology_view.py` paints it, so the widget, a script and a test all draw
the same diagram and it cannot disagree with an analysis of the same
architecture.

A helix the entry does not model is drawn **dashed, never dropped**. Dropping
one would put TM13 where TM1 belongs and silently renumber every helix after it,
on a picture that still looked entirely reasonable — the single worst thing a
topology figure can do. 6B3R greys out TM1–12 and 7WLT greys out TM1–16, both
read from the coordinates rather than written down.

The boxes are a *selection*, not an annotation: ticking a unit lights the same
residues on the 3-D model. When the ticked units are not adjacent the model
highlights the whole span between the first and last, and the status line says
so rather than letting the picture imply otherwise.

### The three Figure 4 views, and what each could be mistaken for

Asked for after the topology diagram: the micelle density (4b), electrostatic
colouring (4c) and a monomer in a planar membrane (4a), all in the main window.
Each is a picture that can be read as a measurement it is not, so each is built
around the specific thing it could be mistaken for.

**The micelle is the hardest case**, because Figure 4b's envelope is *evidence*
— it is a measurement of the detergent, and it is the paper's direct
demonstration that PIEZO1 bends its surroundings. We have no map. What is drawn
is the iso-surface of the distance to the hydrophobic belt, which has two
properties that make it worth drawing at all: it is calibratable (one atom
gives a sphere of exactly the offset radius, a line a capsule, both checked),
and it has **no free shape parameters**. That second point is what lets the
status line divide the picture honestly: the shell *thickness* is a registered
parameter and carries no information, because an offset surface around a sphere
is a sphere with the radius increased by exactly the offset — measured, not
argued, by building it at 7 Å and 13 Å and checking the fitted curvature is
identical to floating point. The *curvature* is a fit to the belt atoms
themselves, and comes out at 9.8 nm against the paper's 10.2 nm idealisation
and our own 10.8 nm dome fit.

**The potential colouring needed a new ColorBy.** The obvious route was to
reuse `ColorBy.VALUE`, which already carries a per-atom array. It auto-ranges
to the 2nd and 98th percentiles, and for a potential that is not a cosmetic
difference: 6B3R's surface never leaves ±2 k_BT/e, so an auto-ranged map would
paint an almost-neutral protein in full saturated red and blue, and it could
not be compared with Figure 4c or with the same protein in another state.
`ColorBy.POTENTIAL` exists to hold the scale fixed at the panel's own ±5.

**The planar membrane draws its own control.** Any point set has a best-fit
plane, so drawing one proves nothing; `use_trimer` draws the same construction
on the assembly, where the paper says it fails, and the status line reports the
slab each would need — 42 Å for a protomer, 60 Å for the trimer, against a real
bilayer's 36 Å.

### The caveats were making the application unusable

Reported from use: running an analysis grew the window past the edge of the
monitor. The cause is this project's own discipline. Long status text goes into
a plain `QLabel`; a non-wrapping label reports its **full text width** as a size
hint, `QStatusBar` passes that up as a minimum, and `QMainWindow` honours it. So
the window widened in proportion to how careful the caveat was, and the worst
offenders are the newest — the micelle's "not the observed density", the pore
surface's two sentences about what a probe sphere is not.

Measured rather than described, because the fix needed a reason: a plain label
demands a **12,566-pixel** window for a 1500-character message, and 1,698 px for
200. The replacement demands the same width at 0 characters as at 4,000.

The fix separates the display from the text, and the second half is the part
that matters. `StatusMessage.text()` returns the **whole** message — every
guard in the suite that asserts a caveat cannot be omitted reads exactly that,
and had elision reached it, all of them would have started passing vacuously
while the user saw less. What is *painted* is elided; the rest is in the
tooltip and in a scrollable history behind a **⋯** button.

The history was not in the ask and is the part I would keep if only one
survived: several controllers set a status while they work and another when
they finish, so the first was previously lost with no way back to it.

### Two defects found on the way

**The bibliography had 6B3R and 6BPZ swapped.** `guo2017` was recorded against
PDB 6BPZ and `saotome2018` against 6B3R. Both PMIDs were right, so nothing
resolved wrongly and no test could see it; the deposited entries themselves say
which is which. This matters more than a typo, because 6B3R is the entry every
Figure 7 number is measured against.

**Adding the cuff broke the anchor.** Guo & MacKinnon give ranges for the elbow,
base and hairpin, and all three sit *inside* domains the project already had.
`Annotations.domain_at` returns the **smallest** containing domain, so adding
them moved the anchor from first to thirteenth in the allosteric-betweenness
ranking and broke the conservation ranking too. The fix is a `sub_element` flag:
the chain's architecture is a partition and callers rely on it being one; these
are named features inside two of its parts. Two existing tests caught it, which
is what they are for.

**`build_disc` had never worked and nothing had ever called it.** It builds a
flat disc through `build_membrane_mesh` with `n_radial=2`, and that function
takes the normals from `np.gradient(..., edge_order=2)`, which needs three
samples. Every call raised. The dome controller draws its flat projection
through `build_membrane_mesh` directly, so nothing had exercised it; the Figure
4a planes are its first caller. Found by a test that renders to a framebuffer
and counts lit pixels — the same assertion that caught the invisible cylinders,
for the same reason.

The pore-extension helix is the one range the paper never states. It is derived
— the pore-lining segment between TM38 and the hairpin, mouse 2479–2500 — and
the derivation is checked rather than asserted: over that range every C-alpha
sits within 11 Å of the three-fold axis, and it contains E2487, H2490 and M2493,
the three residues the paper itself labels as constricting there. Marked medium
confidence, and the description says the boundaries are ours.

### Where the numbers disagree, and why that is left alone

Our pore radii are 0.62 Å wider than the published HOLE ones at all three named
residues. Same constrictions, same order, same closed verdict. Two independent
pore algorithms differing systematically is expected; the offset is pinned in
**both** directions, because zero would mean the profiler had been fitted to the
paper.

Figure 4—supplement 1 names two salt bridges. E2257–R1762 is reproduced by both
conventions, in all three protomers, entirely domain-swapped — the paper's
actual claim, and `patch_interaction` puts it at −6.18 k_BT with a same-chain
term of −0.001. D2264–R1761 is **not** found: 6.43 Å charge-centroid separation
against a 5.5 Å cutoff, though its closest atoms are 4.58 Å apart. The two
conventions disagree about that contact and the paper does not say which it
used, so neither is called wrong.

---

## Round 83 — PIEZO2, and the answer that the mechanism is the fold's

6KG7 has been downloaded since the beginning, entity-classified, and then
excluded from every ensemble as a paralogue. Correct for a PIEZO1 ensemble.
Not an answer to the question the project had never asked: **how much of this
mechanism is PIEZO1, and how much is the fold?**

### The catalogue was wrong about the entry

The registry note said 6KG7 "resolves residues 8-823, so it is the best
experimental view of the distal blade". It resolves **8 to 2822, in sixteen
segments, 1,817 C-alphas per protomer** — *more resolved residues than any
PIEZO1 entry in the catalogue*, which run 1,223 to 1,502 — including all 38
transmembrane helices.

Nobody had asked the file. The note reads like it was written from the paper's
emphasis, and its effect was to make the one structure that could answer the
generality question look like a fragment not worth loading. Corrected, and
pinned by a test that reads the coordinates rather than the note.

### Which protein, and which numbering, measured

The roadmap said PIEZO2 is 2,752 aa. That is *human* PIEZO2. 6KG7 is **mouse**
Piezo2, 2,822 aa — a third length beside mouse Piezo1's 2,547, with no constant
offset relating any pair.

Rather than trust a label, every entry is scored residue by residue against all
four committed UniProt sequences: if a file's numbering belongs to a sequence,
its residue *names* agree with it at every position. Each entry matches exactly
one reference at **1.000** with the runner-up below 0.25. That is a
known-answer measurement, and it fails as it should on a PIEZO1 entry
renumbered by a constant — which is exactly the mistake it exists to catch.

Two new committed resources, `uniprot_mouse_piezo2.json` and
`uniprot_human_piezo2.json`, built by the same script from the same source as
PIEZO1's, so the two dome measurements cannot differ by how their membrane
surface was defined.

### The naive comparison was a coverage artefact, and that is the result

Measured directly, PIEZO2's dome looks dramatically different: 8.51 nm deep
against 4.92, with 462 nm² of excess area against 256. I nearly had a finding.

The cause is that 6KG7 resolves 38 transmembrane helices where 7WLT resolves
22, so the two surfaces trace different amounts of blade. Restricted to the
helices resolved in both — paired by index, a pairing the global alignment
confirms for **37 of 38** — PIEZO2 gives 5.64 nm and R_c 10.32 against 9.72,
inside the PIEZO1 range on every quantity.

Both rows are reported, because the gap between them is a caveat on this
project's own numbers: **dome depth and excess area scale with how much blade
an entry resolves.** Only the radius of curvature is robust to it, which is
fortunate, since that is the one this project has been quoting against the
published 10.2 nm all along.

### The gating coordinate is not PIEZO1's

With the sites coverage-matched through the alignment (1,236 per protomer), the
protomer correspondence **searched** rather than read off chain labels, and
PIEZO2's mode vectors rotated into PIEZO1's frame by the same superposition
that aligns the sites:

| | |
|---|---|
| Overlap of PIEZO1's lowest A mode with one PIEZO2 A mode | **0.804** |
| Fraction inside PIEZO2's symmetric subspace | **0.925** |
| Shuffled-correspondence control | 0.190 |
| Superposition RMSD, 3,708 C-alphas at 48% identity | 4.36 Å |

The protomer order came out **(2, 0, 1)**. Chain labels would have been wrong,
for the third time in this project's history — and across a paralogue there was
never any reason for them to agree.

So the motion identified as the candidate gating coordinate is a property of
the PIEZO fold. That is a result about generality rather than a failure, and it
cuts both ways: the mechanism is more general than one protein, and nothing in
it distinguishes two proteins whose inactivation kinetics and tissue roles
differ. With one PIEZO2 structure it says the fold *admits* the mechanism, not
that every PIEZO uses it.

### The instrument found two things nobody was looking for

The identification exists to stop a PIEZO1 entry being read with PIEZO2's
transmembrane annotation. It found two live defects on its first pass over the
catalogue, and both matter because this project applies every annotation —
domains, helices, variants, functional residues — **by residue number**.

**6LQI is deposited in the splice isoform's own numbering.** Piezo1.1 lacks
residues 1382–1405, and the file numbers straight across the deletion: 1.000
agreement with canonical mouse Piezo1 before the splice site, 0.058 after, and
1.000 again shifted by **+24**. That is 764 of its 1,301 resolved residues.

**Four human entries carry a block numbered 22 low.** 8ZU3, 8YFC, 9VMX and
8YFG score 0.932 — comfortably past any floor I would have set — and the 7%
they are missing is not spread out. It is residues 767–857, ninety-one of
them, every one disagreeing and every one agreeing again read +22. 8YEZ
resolves the same region and does not have it.

The second is the more instructive. A whole-file identity of 0.932 looks fine;
runs of disagreement are what makes it visible, and reporting only the total
would have passed four entries with a real fault in them. Both are recorded as
Round 86 rather than fixed here, because fixing them means re-reading five
entries everywhere and recomputing every number this project has published for
them.

**And one thing I got wrong first.** 3JAC scored 0.623 and I had it as a third
case. Every single mismatch was a `UNK` — the depositor declining to name a
residue rather than disagreeing about one. `AA3TO1` maps UNK to X, so
membership in it is not the test, which is what my first filter used. Excluding
unassigned residues, 3JAC matches at 1.000 over the 572 it names. Three
findings became two, and the correction is in the module docstring because the
next person to write such a filter will reach for `in AA3TO1` too.

### A definition that existed twice

The dome surface — one point per transmembrane helix per protomer — was defined
inline in the report and again in the claims registry. Survivable while only
PIEZO1 was measured; not survivable the moment a second protein had to be
compared against the same surface, because two definitions is exactly how a
comparison ends up measuring how each side was defined. Extracted to
`structure.geometry.tm_surface_points`, with both callers collapsed onto it and
the 9.72 nm claim unmoved.

### And the catalogue quietly grew a HaloTag

Rebuilding `structures.json` swept **6U32 — the HaloTag crystal structure** —
into a registry documented as a catalogue of PIEZO structures, as an
unclassified entry of "unknown" species. The build globs every `.cif` in the
structure directory, which was harmless until Round 31 downloaded the tag for
the fusion geometry and nobody rebuilt until now.

It took two unrelated tests failing to notice: the ligand audit found the tag's
tetramethylrhodamine and reported it as a possible undocumented PIEZO
modulator. `fusion.load_halotag` reads 6U32 by path and never through the
registry, so the exclusion costs nothing — and every registry entry must now be
identifiable as a PIEZO, which the round's own instrument checks.

Suite 1320 → 1344 passing (1354 collected, 10 skipped for want of
downloaded data); 116 → 121 registered parameters. Reachable as
`piezo1 paralogue <entry>` and under Analysis → PIEZO2 comparison.

---

## Round 82 — the B-factors, read at last, and what they say about the network

Every structure this project loads carries a `b_factor` for every atom and no
analysis had ever read one. That is the first thing anyone does with an elastic
network — does the predicted mean-square fluctuation track the deposited
B-factor? — and this project's central mechanism claim rests on such a network.

**One correction to the premise.** The roadmap said the ANM has no
scalar-fluctuation counterpart. It has: `ModeSet.msf` sums `|v|²/λ` over modes
and both the feature table and the fluctuation colouring already consume it.
What was missing was the *comparison*. That changes what this round is — a
missing validation rather than missing physics — and it is worth writing down,
because "we never built X" and "we built X and never checked it" call for
different responses.

### The column decides more than the network does

An observed B-factor is a temperature factor only in the sense that refinement
put it there. In a cryo-EM map it absorbs local resolution, sharpening,
per-particle scaling and the restraints. So the comparison is gated on the
column before any Hessian is built, and three kinds are refused:

- **uniform** — carries no per-residue information;
- **grouped** — 3JAC and 6BPZ each carry **212 distinct values over ~2,700
  C-alphas**, one per thirteen residues, and cannot resolve what is being asked
  of them;
- **pLDDT** — an AlphaFold model puts a *confidence* in that field, high where
  the model is certain, which is where a real B-factor would be low.

The last gate is decided by provenance, which is a decision rather than a
measurement, so the measurement is made separately: build the network on the
AlphaFold monomer and correlate its own column, and it comes out at **−0.57**.
The gate points the right way, and that is checked rather than asserted.

### The control is what makes the number mean anything

A residue with more neighbours moves less. That is true of any packed solid
with no normal modes in it at all, so a correlation between an ANM and a
B-factor can be entirely burial. Every result therefore carries the same
correlation for **contact number** — no Hessian, no eigenvalues, no gating
coordinate.

That control is what turned a clean-looking result into an honest one:

| | Network | Contact number |
|---|---|---|
| Median Spearman | **0.74** | 0.32 |
| Median Pearson | **0.48** | 0.39 |
| Wins on Spearman | **13 of 15** | — |
| Wins on Pearson | **9 of 15** | — |

**The elastic network orders PIEZO1's residues by mobility much better than
burial does, and predicts how far they move barely better.** Had I reported
only the rank correlation this would have read as a validation; had I reported
only Pearson it would have read as a failure. Both are in the module, in the
result dialog and in `SCIENCE.md`, and a test fails if the asymmetry ever
disappears without the prose changing.

### Three things the survey turned up

**18 of 21 entries can answer**, and the three that cannot say why rather than
being dropped.

**Three entries have a *negative* control** — 8YEZ, 8ZU8 and 6B3R. Their
B-factor *rises* with burial, which no mobility does. On the first two the
network scores 0.10, and the honest reading is that the column is not a
temperature factor rather than that the network failed. They are excluded from
the counts above; counting them as wins, which the naive comparison does, would
have added two spurious successes.

**The two entries the network loses on are named in the test**, so a change in
either reopens the question. One of them is **6KG7 — PIEZO2**, where burial
predicts the column at 0.55 and the network at 0.07. That is a hint about
Round 83's question rather than an answer to it.

### The calibration, and the thing it had to be able to say

A planted fluctuation generated *by* the network is recovered exactly. The same
comparison against a network built on **shuffled coordinates** — same beads,
same cutoff, same spring model, only the neighbours changed — collapses. A
shuffled observation returns nothing, bounded by the null's own standard
deviation rather than by a round number, because at sixty residues a
correlation of 0.25 is two sigma and a fixed threshold would either be vacuous
or flake.

`fluctuations` is registered in `tests/test_calibration.py` as a checking
instrument in its own right, since that is exactly what it is: an instrument
for checking the network, and therefore subject to the rule.

Suite 1304 → 1320 passing (1330 collected, 10 skipped for want of
downloaded data); 115 → 116 registered parameters. Reachable as
`piezo1 fluctuations <entry>` and under Analysis → Fluctuation vs B-factor.

---

## Round 81 — the pore gets its charge, and a sign error falls out

`solve_pnp` has taken a `fixed_charge` argument since Round 33. Its documented
equation carries the ρ_fixed term. **No caller had ever supplied one.** So every
current this project has published was computed for an electrically neutral
pore — a cation channel with nothing in it that could prefer a cation.

### The argument that was already there

Wiring it in is not just passing an array. Fixed charge only means anything
through the *equations*, and the electroneutral limit the solver already uses is
where it belongs: a charge density sets a local Donnan potential, counterions
are enriched against it and coions excluded, and the concentration gradients
that creates carry a diffusion current. That current is what gives a pore a
reversal potential, which is what selectivity is measured as. The same
double-layer overlap that stopped the Poisson coupling converging in Round 33 —
5.7–8.1 Å screening against a 3.3 Å radius — is what makes the uniform-potential
limit the right one here, so the justification was already measured.

### The calibration failed, and it was supposed to be able to

Before believing anything about PIEZO1 I ran three known answers: a cation must
drift downhill, the GHK inversion must return its Nernst and unity limits, and
an uncharged pore must return its ions' mobility ratio rather than one.

The third came out **inverted**: a pore given more negative charge became more
*anion*-selective, monotonically and smoothly. Suspecting the checker first is
the standing rule here, and it was right to — but this time the checker was
fine. The Scharfetter-Gummel drift term had its two Bernoulli factors attached
to the wrong nodes, so cations drifted *up* the potential gradient. The
decisive test needs no discretisation at all: at a field weak enough that SG
reduces to a centred difference, the flux must equal `-D A z c φ'/φ_T`, and it
came back with the right magnitude and the wrong sign.

Nothing could have caught it. Every current the project had computed was
between identical baths, where the concentration term vanishes and reversing
the field only reverses the sign — and the sign was then discarded by
`pore_ohm = abs(voltage / current)`. Two independent checks agreed with each
other (`series_conductance` to 1.5%) because both are magnitudes. Correcting it
moved the recorded 41 pS by one part in 10¹⁴.

### And then the closure was wrong too

With the sign fixed and real charge supplied, the Gummel loop did not converge
and the electroneutrality residual reached **43 times the ionic content**. The
cause is a discretisation mismatch: I had extended the ohmic closure with a
diffusion-current source term, discretised centrally, while the flux is
discretised Scharfetter-Gummel — and where a carboxylate sits the potential
changes by four thermal voltages across one Ångström, at which point the two
disagree by more than the term itself.

The fix is not a better source term, it is the right closure. A charged pore's
potential is set by **local electroneutrality**, not by Ohm's law; the ohmic
operator contains no term the fixed charge could enter through at all. Imposing
electroneutrality directly converges in tens of iterations, drives the residual
to 1e-10, and reproduces Donnan equilibrium exactly — zero current and Boltzmann
partitioning to machine precision at zero applied voltage. The neutral path
keeps the old closure, untouched, which is how an explicitly zero charge still
returns the old numbers bit for bit.

### What the charge is, and what it is not

The residues came out more interesting than the wiring. Positions have to come
from C-alpha, because **11ZC — the only open structure — is the only entry
deposited without side chains**; measuring the others from their real side
chains would make them incomparable. A residue counts if its charge, on a fully
extended side chain, could reach the lumen at its own height, which is
deliberately permissive.

On that criterion **three of the four curated "selectivity glutamates" are not
within reach of the lumen**. E2117 sits 12.9 Å past the wall. That is not a
contradiction of the curation: Coste et al., who identified E2117 (mouse E2133)
by mutagenesis, concluded from function alone that it *"may not lie in the
selectivity filter but could be located close enough to the pore to
allosterically modulate its properties"*. The geometry and the
electrophysiology agree without either having been fitted to the other, which
is the most satisfying thing in the round.

### The result

| Route | Charges | Net | Conductance | P_Cl/P_Na |
|---|---|---|---|---|
| none | 0 | 0 | 40.1 pS | 0.904 |
| curated | 6 | −6 e | 29.6 pS | 0.021 |
| every group reaching the lumen | 46 | **+8 e** | 4.1 pS | 0.207 |
| measured (Coste 2015) | — | — | 25–30 pS | 0.14 |

**The direction is right and the number is not.** Both routes make the model
cation-selective and they bracket the published ratio tenfold apart. Three
things stop that being read as agreement, and all three are reported beside it:
the uncharged pore is *already* cation-selective at 0.904, because chloride is
nearly twice sodium's radius and loses more cross-section at a 3.3 Å bottleneck
than it gains in mobility; the curated route only reaches 0.021 at an in-pore
concentration of **13.9 M**, which is flagged as past any packing a solution
could reach; and the two routes are net negative and net *positive*
respectively, so they disagree in kind rather than in value.

That the curated route also moves the conductance from 41 pS to 29.6 pS, into
the published 25–30 pS band, is recorded and explicitly not claimed. It is one
of several free choices moving a number the module already documents as
spanning 16–94 pS across two unmeasured parameters.

### A guard that earned its keep

The four side-chain reach parameters were reported as **read by nothing** the
first time the suite ran. They were read — through a dictionary of key strings,
which `provenance_chain.resolved_keys` cannot see, because it scans for the call
rather than running it. That is the failure mode the check exists for: such a
parameter still appears in the dialog, still trips the non-default banner, still
stops `verify_claims`, and does nothing. The four calls are now written out, and
a test holds the two routes together so they cannot drift.

Suite 1275 → 1304 passing (29 new tests, all but four of them known-answer
cases run before any PIEZO1 number); 104 → 115 registered parameters; two new
claims so the selectivity numbers cannot drift out of `SCIENCE.md`. Three files
reached the 500-line limit on the way and were split at real seams:
`_pnp_kernels.py` took the discretisation, `claims_structural.py` took the claim
recomputations that need downloaded coordinates, and `help_topics_tags.py` took
the longest help topic.

---

## Rounds 78–80 — three documentation rounds that found three code defects

Block Q's remaining items were all documentation. Each turned up something in
the code on the way, which is the argument for doing them rather than declaring
them cosmetic.

**Round 78 — retiring `HALOTAG_CALCIUM_PLAN.md`.** Rounds 29–32 implemented all
of it. The validation clause was that nothing may be lost, and **two of its five
recorded risks were carried nowhere**: that PNP is a mean-field theory applied
to a pore two ions wide — good for a conductance, useless for single-ion
energetics — and that using the wetting heuristic as an on/off switch is a
stronger claim than AUROC 0.91 over ~200 channels was validated for. Both are
about over-claiming, which is what this project exists to guard, so both moved
to `SCIENCE.md` §8f before the file went.

Checking one of those claims — *is* the wetting verdict reported beside the
conductance rather than multiplied into it? — I could not import
`analysis.report_tags` at all. It imported a helper from `report`, and `report`
imports `report_tags` at the bottom, so a fresh interpreter got a
half-initialised module. The suite never noticed because something always
imports `report` first; a cycle only bites the person who reaches for one
module directly, which is a notebook user. The helper turned out to be a
**duplicate** of `structure.protomers.protomer_blocks` — identical output on
every trimer, differing only in its non-trimer sentinel and in hardcoding the
300-C-alpha floor instead of using `well_resolved_chains`. Deleting it fixed
both.

**And then my own checker broke seven unrelated tests.** `test_imports.py`
clears `sys.modules` to reproduce a fresh interpreter, which is right, and my
first version restored nothing. `PARAMETERS` is a singleton: other test modules
bind it at import time, so leaving fresh module objects behind made the rest of
the session resolve a *different* registry, and overrides silently stopped
taking effect. The failures looked like parameter-wiring bugs and were nothing
of the kind. The restore now saves and replaces the original objects **by
identity**, and two tests pin that.

**Round 79 — `ARCHITECTURE.md`, written rather than deleted.** The easy way to
close it would have been a summary of the module layout, which is what
`INTERFACE.md` already is; a second copy would go stale differently and settle
nothing. So it is organised around constraints, and each section names the
incident that exposed one. The test caps module-shaped rows at fewer than five,
and — more useful — checks the architecture it *describes* is the one the code
still has: the one-way dependency arrow is verified against real imports rather
than asserted.

**Round 80 — how far is the conclusion from where a reader starts?** Measured:
**four of seven entry points could not reach it in one step** — `SCIENCE.md`,
the guided tour's own closing step, the command line, and the notebooks index.
`ROADMAP.md` was a fifth, found by the test itself. The roadmap's clause is that
the surface is wrong rather than the reader, so all five were fixed.

Two of my checks were wrong before the things they checked, again. Reading the
tour's static body missed the count because the tour **computes** it from the
validation record — which is the better design and the reason it cannot go
stale. And the wrong-count guard was too narrow to catch a surface drifting to
a different number, which I found only by mutating one.

Suite 1116 → 1271.

## Round 77, and a codebase review against comparable software

**Round 77 — a fetch that verifies what arrived.** The size floor has been
necessary and insufficient twice: Round 60 found an Ensembl endpoint returning
HTML, Round 65 found two 127-byte error pages *stored as structures*. Both were
caught only because they happened to be tiny. `_download` now takes a `kind` and
refuses to write anything that is not it.

The load-bearing check is that the guard rejects **what the size check
accepts** — a 554-byte HTML error page, comfortably past the 200-byte floor. A
guard that only fired on things the old one already caught would be worth
nothing. Nothing is written before the check, which matters more than it looks:
`_download` serves an existing file from cache *without re-checking it*, so a
bad file written once is served forever. That is what made Round 65's failure
durable.

Writing the "a real mmCIF still passes" half, I served the first 200 kB of 8YEZ
and it was refused. The guard was right and my test was wrong — mmCIF metadata
runs well past 200 kB before the first coordinate — and the behaviour is worth
keeping, so it is now its own test: a connection that drops mid-transfer
otherwise leaves a file that opens, parses, and contains a fraction of the
molecule.

**The review.** Compared against HOLE, CHAP, MOLEonline, ProDy and the
APBS/PDB2PQR route, three gaps are statable as measurements rather than
opinions, and they became Block R:

1. **`solve_pnp` takes a `fixed_charge` argument that no caller has ever
   supplied.** The documented equation carries a ρ_fixed term; every permeation
   number this project has produced treats the pore as electrically neutral.
   Meanwhile `functional_residues.json` has curated four sequence-verified
   acidic residues as "acidic residues setting ion selectivity" since early on —
   twelve charges across the trimer — and `default_species()` offers a generic
   "cation" and "anion" not even named Na⁺, K⁺ and Ca²⁺. **A cation channel
   modelled without charge cannot be selective, and this one is not.** This is
   the strongest item in the block: the machinery, the annotation and the
   equation are all present, and only the wire between them is missing.
2. **`b_factor` is parsed for every atom and no analysis reads one.** Predicted
   fluctuation against observed B-factor is the first validation a ProDy user
   runs, and this project has never run it on the network its central mechanism
   claim rests on. It needs calibrating before it is believed — sharpened
   cryo-EM values, backbone-only entries, and AlphaFold files carrying pLDDT in
   that column are three ways to get a confident wrong correlation.
3. **PIEZO2 (6KG7) is downloaded, classified, and only ever excluded.** Correct
   for a PIEZO1 ensemble, wrong as a final answer: it is the obvious control for
   the question never asked — how much of this mechanism is PIEZO1 and how much
   is the fold?

Plus a fourth from the interop comparison: **nothing computed can leave the
application.** Conservation, coupling, PRS response and mode displacement are
all per-residue scalars, and `to_pdb` writes coordinates only, so the standard
route into PyMOL or ChimeraX — the scalar in the B-factor column — does not
exist.

Block R is deliberately the first block in a long while whose items could each
return a *positive* result, so Round 84c's review carries a standing question:
check the same discipline held.

Suite 1101 → 1116.

---

## The empty `notebooks/` folder — what it was, and filling it

George asked what the empty `notebooks/` directory was for. The answer is
unflattering: it was created on **the first day of the project** and never
filled. Git does not track empty directories, so it was never in a clone at
all — it existed only on this machine, which is why nothing ever noticed it,
and why no audit could have. `docs/NOTEBOOK.md` took its place, documenting the
headless API in prose.

Prose is not a substitute, because prose cannot be executed. Four notebooks
now: what is in a structure and how to frame it; the elastic network model and
the symmetry rule; the pore through wetting to a current; and the variant
workflow, walked to the point where it stops working.

**They are generated, not hand-written.** The cell content lives in
`scripts/notebook_content*.py` as ordinary Python so it can be reviewed and
diffed, and `build_notebooks.py` **executes every code cell, in order, in one
namespace, before writing anything.** That gate earned its place immediately —
it caught six wrong API calls I had written from memory, including
`ANM.from_structure` (does not exist), `measure_dome` missing its surface
argument, and a dimension mismatch from building the network on one residue set
and the displacement on another.

Two of those were more than typos. The gating notebook needed the project's
real route — resample both structures onto shared residues, `match_protomers`
because chain labels do not encode rotational order, and Kabsch before taking
the displacement — and every one of those steps produces a plausible wrong
number if skipped rather than an error. They are now what the notebook teaches.
And my frame check took "the last 200 rows" of the C-alpha array, which is the
exact defect `test_frame` documents: the slice straddles chains, and it read
+7.7 Å where the cytosolic end must be negative.

**No stored outputs.** A committed output is a number nobody recomputes; it
goes stale silently while reading as authoritative, which is the failure this
project spends most of its machinery avoiding. The notebooks `assert` the
numbers they quote instead, so running one checks the science and not only the
syntax.

`tests/test_notebooks.py` executes the **committed JSON** rather than the
content module — the builder already runs its own source, so re-running it
would prove nothing about what a reader downloads. It also pins that a
hand-edited `.ipynb` is caught before the next build discards it, and that
notebook 03 puts the parameter registry back: a leaked override would poison
every later test, and `verify_claims` refuses to run against one.

Jupyter is **not** a build dependency — nbformat 4.4 is emitted directly, and
the cells are executed in-process. `pip install -e ".[notebooks]"` adds
JupyterLab for anyone who wants to open them. `make notebooks` rebuilds; that
needed `.PHONY`, since the target name collides with the directory and make
otherwise decided it was already built.

Suite 1076 → 1101.

---

## A documentation pass — and two claims I invented while making it clearer

George asked for the docs to be reviewed and updated, the README rewritten in
clear English with citations and a reference list, and new figures where they
help.

**What was actually wrong with the README.** Not only the prose. Its *Status*
section listed as "Planned" a Helfrich membrane solver, tension-dependent
Markov kinetics, conformational morphing, pore profiling, pocket detection and
docking — every one of which had shipped, some dozens of rounds earlier. A
reader would have concluded the project did a fraction of what it does. It also
predated the HaloTag work, ion permeation, the nanodomain, the full-length
model, the parameter registry and the guided tour, none of which it mentioned.

Rewritten around what a reader wants to know in the order they want it: what
PIEZO1 is and why its shape *is* the mechanism, what the software does grouped
by task, the numbers it reproduces against the published ones, how to install
and drive it, and the closing record — which is guarded by `test_conclusion.py`
and was left intact. Citations are inline and there is a numbered reference
list with DOIs, all 22 of them traceable to `docs/REFERENCES.md`, which is
itself built behind a title-verification gate.

**And I introduced two false claims while doing it**, both in one paragraph,
both by writing from memory instead of checking. I said E756del is carried by
roughly a third of people of African ancestry — the project's own `SCIENCE.md`
records the gnomAD AFR frequency as 0.166–0.173, so a sixth, which is what the
old README correctly said. And I wrote that it protects against severe malaria,
where `SCIENCE.md` records the association as **contested**: a later study found
OR 0.91, p = 0.19, and the original mouse work tested R2482H, a different
allele. Both are now stated as the project's own data states them, and
`tests/test_readme.py` pins each so the same slip fails rather than ships.

That test also caught something real. Checking that every README figure can be
rebuilt, it found `gating_morph_small.gif` committed with **no script that
produces it** — `docs/anim/` is git-ignored, so anything the README shows must
be committed, which makes it exactly the file that quietly stops matching the
code. It now comes from `make_animations.py --only readme`, the same code path
at a smaller frame. My first version of that check was also too naive to be
trusted: it scanned the scripts for each figure's name and flagged `app_modes`
and `app_pore`, which are built from an f-string. Suspect the checker first —
it was two-thirds wrong and one-third right.

**Two new figures**, both for features that are hard to describe in words:
`hybrid_model.png` (the graft in AlphaFold's confidence colours against the grey
experimental core, seam marked) and `halotag_fold.png` (the real tag fold at its
modelled position). `scripts/make_model_figures.py` generates both. The tag
figure draws the channel in a uniform colour rather than by chain, because the
tag's orange sits 0.10 from the chain palette's orange and a per-protomer
picture makes the modelled tag look like part of the experimental trimer.

**`docs/NOTEBOOK.md`** covered none of ion permeation, the full-length model,
the fusion pose or the nanodomain. All four added, and every snippet was run
against real data before being written down — which caught two wrong signatures
I would otherwise have documented.

Suite 1065 → 1076. All 21 documented numbers still reproduce.

---

## Round 76 — the full-length model reaches the GUI, and the colouring it needs was broken

`structure/hybrid.py` had existed since Round 65 and nothing outside a notebook
could build one — the same exposure gap Round 58 found for the coupling score.
It is now reachable three ways, all through the one shared `ANALYSES` entry so
they cannot diverge: **View → Full-length model** draws it, **Analysis →
Full-length model numbers…** tabulates it, `python -m piezo1.cli hybrid <PDB>`
prints it. The registry's own guard caught me when only the drawing existed —
`test_every_shared_analysis_is_reachable_from_the_gui` failed on the entry I
had just added, which is the mechanism working.

**The validation clause was the interesting part.** Not "it draws" but *the
seam must be visibly rendered*, because a complete-looking PIEZO1 trimer whose
569 distal residues are AlphaFold is exactly the confident-wrong-picture failure
the Round 50 hazard audit exists for. So the experimental core is flat grey —
dull on purpose, and not any of the real colourings, so it reads as background —
the graft is coloured by pLDDT in AlphaFold's own bands, the seam carries a
marker, and the status line states the range, the 48% clearing pLDDT 70, and the
**75 Å** the two models differ by away from a seam that fits to 2.4 Å. That last
number is the one a good local fit hides.

**And the signal the whole design leans on did not work.** `plddt_colors`
applied the bands in the order `PLDDT_BANDS` declares them — highest first — so
the final pass at `>= 0.0` overwrote every atom. **Colour by → AlphaFold pLDDT
painted the entire model one flat orange**, and had since the feature was
written; on PIEZO1's prediction, which spans 22–95 pLDDT with 2,626 atoms above
90, every one came out "very low". A confidence colouring that shows no
variation is worse than none, because it reads as uniformly bad.

Fixed once, in a shared `plddt_band_colors` the controller reuses rather than
copying — the copy is how I found it, since I wrote the same loop into the
controller and the graft came out one colour. Four bands, four colours, each
value taking the highest threshold it clears.

Suite 1053 → 1065.

---

## Round 75 — retiring the roadmap, and what its length was hiding

`ROADMAP.md` had reached 2,702 lines, of which **96% was a record of work
already finished**, duplicating `SESSION_LOG.md`. The question the file exists
to answer — *what is left?* — took scrolling past three hundred ticked boxes to
find.

Split into a 108-line `ROADMAP.md` holding only open work, and
`docs/ROADMAP_COMPLETED.md` holding the finished record verbatim. The rule was
mechanical — a section goes wholesale to one file or the other — and it turned
out every section was entirely open or entirely done, so nothing had to be cut
in half and no completed item could be stranded beside an unfinished one.

**The validation clause asked for a count before and after; it got a
comparison.** A number typed in from memory proves nothing, so
`tests/test_roadmap.py` reads the pre-split file **out of git** at commit
`4c1c61c` and compares item by item: every one of the 358 completed items is
present afterwards, none appears in both files, every open item is in
`ROADMAP.md` and none in the archive. Frozen counts back it up where history is
unavailable, and the completed count is ratcheted so it can only grow. Deleting
three items, filing one in both files, and tidying away the awkward Round 68
entry each fail it.

**What the length was hiding.** Two adjacent headings, both numbered *Round 68*
— one ticked as superseded by Round 63, one still open asking the identical
question. Round 67 had recorded the supersession by *adding* a heading rather
than ticking the original, and nobody saw the duplicate because nobody reads
three hundred ticked boxes. Merged, with the duplication itself written into
the entry rather than tidied away.

My first description of this was wrong and is worth recording as such: I read
the open copy alone and wrote that "Round 68 asked a question Round 63 had
already answered, and nobody noticed". The supersession *was* noticed and
recorded — in a second heading, six lines above. The defect was duplication,
not oversight. Corrected before it reached the archive.

This is the second stale roadmap entry found (Round 53 was the first), and both
survived for the same reason. That is the argument for the split, stated as a
measurement rather than a preference.

---

## Auditing every control, and the one that had been wrong all along

Asked to confirm that the mouse, keyboard and menu controls all work, the only
honest way to answer was to fire them rather than read them — and reading would
have missed the bug, because the handler looked correct.

**Every rotation ended in a pick.** `mouseReleaseEvent` told a click from a drag
by the distance from `_last_pos` — which `mouseMoveEvent` overwrites on every
step to compute the drag delta. By release time it *was* the release position,
so the distance was always zero and a 390-pixel rotation registered as a click.
Harmless while a pick only rewrote the status bar, which is why it survived;
the moment a pick highlights, every turn of the structure repainted the
selection. Fixed by keeping the press position separately (`CLICK_SLOP = 3`,
enough for trackpad tremor).

Everything else checked out: 75 menu actions, all with their own tooltips, no
shortcut collisions, none disabled, and all 55 non-dialog ones fire without
raising against a loaded structure. Left-drag rotates, shift- and middle-drag
pan, right-drag and the wheel zoom, `R`/`O`/`Space`/`+`/`-` all act, and the
viewport takes focus so the keys are reachable at all.

**One of my own checks was wrong first**, which is the standing lesson again in
miniature: I tested panning by watching `camera.pivot`, which `translate` does
not touch — it moves `camera.pan`. Reported as working once measured properly,
rather than filed as a bug.

Four bindings were undocumented — middle-drag, right-drag, `O` and `+`/`-`. The
shortcut list is now complete, and `test_ui_controls` reads the viewport's own
`keyPressEvent` for `Qt.Key.*` names and fails if one is missing from the help,
so a new key cannot be added unannounced.

## A right-click menu, built so it cannot drift from the panels

George asked for right-click to open a context menu of useful actions. Right-
*drag* keeps the zoom and the menu is on the *click*, which is the same
distinction the left button already makes between picking and rotating — a user
dragging to zoom never gets a menu they did not ask for.

Two design rules, both tested rather than asserted:

*Nothing is implemented twice.* The Representation and Colour-by entries set
the **Model panel's combo boxes**, not `view.style` directly. Calling the setter
would change the model while the panel went on displaying the old value, and
the two would disagree the first time anyone used the menu. The test watches
the combo move, not the code.

*Opening the menu is not a selection.* It identifies the residue under the
cursor so the entries can name it, and dismissing it must leave the model as it
was. That needed splitting `_pick_at` into `atom_at`, which answers the
question, and `_pick_at`, which announces it.

`Add to measurement` **arms picking** rather than refusing when it is off. This
menu is exactly where a user who never found the Measure button ends up, and
telling them to go and press it first would repeat the problem it solves.

Shown with `popup()` rather than `exec()`: `exec()` spins its own event loop and
does not return until the user chooses, so the whole path would have been
untestable — the first probe hung for ten minutes and that is how I found out.

Suite 1002 → 1040. `ui/context_menu.py` 240 lines.

---

## Selecting atoms — a feature that existed and could not be found

George reported three things while testing: clicking the structure put a
selection in the status bar but not in the Measure panel, nothing was ever
highlighted on the model, and it was not clear whether selecting required some
other click.

All three were real, and they were separate.

**Picking is armed, and nothing said so.** Clicks reach the Measure panel only
when its button is toggled on. That gate is correct and worth keeping: a click
already means "tell me about this residue", and a measurement tool that quietly
consumed those clicks would break inspection. But the button said `Measure`,
which reads as a verb you press *after* selecting, and the hint under it said
"Pick 2 atoms for a distance" — the goal, not the step, phrased as though
clicking would already work. The button now says **Start picking** and then
**Picking — click atoms**, and the hint states which of the two states you are
in and what to do in it.

**A pick appeared only when the measurement completed.** `_refill` drew
`self.set.measurements`, so selecting one atom of two left the Selection table
empty — exactly what was reported. Pending picks now get their own row,
dimmed in the pick marker's blue, with the number still to go, so the table
shows the selection rather than only the result. Deleting that row abandons the
half-made selection; it sits past the end of `measurements`, which a naive row
index would have got wrong.

**A click marked nothing.** `_on_pick`'s unarmed branch built a status string
and stopped. The machinery to highlight was already there — `_highlight` draws
a gold sphere batch and is used by the annotation and analysis panels — it was
simply never called from a click. It is now, with one change: `_highlight`
takes an optional `chains`, because a residue *number* means all three
protomers and that is what an annotation wants, while a click means the one
copy under the cursor.

Discoverability is also paid for in the status bar, but only three times
(`PICK_HINTS`). A permanent hint would crowd out the residue identification
that is the status bar's actual job.

**Why this survived.** `MeasurementSet` had fourteen tests and the panel that
drives it had none, so the whole wiring layer was uncovered — including the arm
gate, the table refill and the click routing. `tests/test_ui_measure.py` covers
it with the real widgets offscreen; reverting each of the three fixes fails 1,
5 and 1 tests respectively.

`main_window.py` reached 497 lines, so selection, camera focus and the pick
path moved to `ui/selection.py` as `SelectionMixin` — the seam this work had
just drawn, and consistent with the five mixins already there. 373 + 149 lines.
Suite 987 → 1002.

---

## Drawing the real HaloTag fold — a user request, and the finding it produced

Not a roadmap round. George asked, while testing the GUI, whether the three
orange spheres under **View → HaloTag fusion → Show modelled tags** could show
the tag's actual structure. They can: the tag's fold is experimental (6U32) and
already downloaded — what is modelled is where it sits, not what it looks like.

**Why the sphere was there in the first place, and why replacing it needed
care.** `fusion.py` deliberately produces a region rather than a pose, and the
sphere was the visual form of that refusal. A drawn fold is more informative and
more dangerous: it looks like a determined structure. So the fold is placed with
the freedoms named — position and seam direction from the model, the spin about
the linker left free — and `View → HaloTag fusion → Turn tag orientation` turns
it, because a user watching the fold spin while nothing else moves understands
"undetermined" in a way no caption achieves. The guard is a test that no path
can put a fold on screen without the status line saying so; it caught a real
case, where switching dyes on replaced the caveat with the dye count.

**The measurement, which is the reason this was worth doing properly.**
`accessible_volume` treats the tag as a sphere of its radius of gyration and
says in its own docstring that the real fold, reaching 30.0 Å, "clashes where
this says it does not". With real coordinates that becomes a number. Over 36
spins: 7WLT 27 clear, 8YFG 7, 8YEZ 1, 11ZC 0. The two models agree on the
question that matters — 11ZC is the one entry whose sphere clearance (15.7 Å)
falls below the radius of gyration (17.6 Å), and the one where no orientation of
the fold clears — while the sphere is generous about how much room there is.
That is a genuine cross-check of an approximation that had only ever been
asserted.

**And the eleventh instance of the standing lesson.** The first version counted
every touching atom and reported **0 of 36 orientations clear on all four
structures** — the fold apparently contradicting the sphere. It did not. The
persistent contact was the tag's own N-terminal residue against PIEZO1's
C-terminus: the two ends of the linker, which the placement rule deliberately
points at each other. The instrument was reporting its own construction as a
finding, and it returned a plausible number rather than an error, exactly as
CLAUDE.md warns. Excluding the anchor residue is what turned a manufactured
disagreement into the agreement above. `test_the_attachment_residue_is_excluded_and_that_is_what_mattered`
pins it on synthetic coordinates so it holds even if the deposited structures
change, and the three mutations tried against the suite — dropping the
exclusion, placing the fold at the anchor, making the spin a no-op — each fail
four to six tests.

**A claim corrected rather than defended.** The controller said its colours were
"deliberately unlike any colouring the channel uses". Measured, the tag's orange
is 0.10 from the chain palette's orange and the dye's red 0.10 from its red, and
there is nowhere to move: every colour far from those eight hues is too dark to
see against the background. The comment now says what is true — colour is not
the guard, the status line is — and a test pins both halves together so that
finding a free colour cannot become a reason to drop the caveat.

Registered `fusion.pose_contact_distance` (3.4 Å, twice Bondi's carbon radius);
the finding is insensitive to it from 2.0 to 4.0 Å. Suite 955 → 987 passing
(20 new in `test_fusion_pose.py`, 12 in `test_ui_fusion.py`);
`screenshot_app.py --structure 8YEZ` clean; the fold verified through the real
GL path, not only the fake scene.

---

## Session handoff — paused after Round 74

**State: clean, after correcting a regression Round 74 shipped.** Head is
`74680ac`. Full suite 955 passed / 10 skipped in 5:05; `make coldclone` clean in
61 s (628 passed, 0 failed, 337 skipped); no file over 500 lines.

**The regression, because it is worth remembering.** Round 74's
`test_the_script_runs_and_reports` ran `cold_clone_check.py`, which clones the
repository and runs the suite *inside the clone* — and the clone contains that
same test, which runs the script again. Unbounded recursion, timing out at 1800
seconds and making both `pytest` and `make coldclone` unusable.

It passed when written **only because the test file was not yet committed**, so
the clone did not contain it. The bug activated on commit: green locally, broken
for everyone afterwards. I then reported the suite as clean on the strength of
that run, which was wrong and is corrected here.

The repair is not the obvious one. An environment guard breaks the recursion,
but `git clone` copies *committed* state — so an end-to-end test that clones
HEAD exercises the last commit rather than the working tree. It cannot see the
change being made. A test that cannot see the code under test does not belong in
a suite, so the end-to-end run is now `make coldclone` only, and the eleven
remaining tests cover the logic against planted inputs.

**Where the roadmap stands.** Rounds 1–74 are complete or explicitly superseded.
The remaining unchecked items are:

- **Round 75 — retire the roadmap itself.** *Started but not begun in earnest:*
  the only work done was measuring the file, and no edits were made. The task is
  to split `ROADMAP.md` into "what is not done" and a completed-rounds record,
  since the changelog half now duplicates this log. Its validation condition is
  the hard part: **no completed item may be lost**, and a test should count the
  measured results before and after the split.
- **Block Q (Rounds 76–80)**, appended after the Rounds 61–65 review: make the
  hybrid model reachable from the GUI with the seam visibly rendered; a fetch
  that verifies *what* it downloaded rather than only the size; retire
  `HALOTAG_CALCIUM_PLAN.md`; write or delete `ARCHITECTURE.md`; and measure how
  many steps a reader needs to reach the conclusion from each entry point.

**Where the science stands, and why it is finished.** Five pre-registered tests,
five predictor families, five nulls. Round 47 closed the across-position route
(134 variants needed, 59 reachable) and Rounds 54 and 61–62 closed the
within-position one (8 shared positions needed at an implausibly good predictor,
one available, 3–4 reachable). Round 64 recorded in writing that no further test
will be attempted on this variant set, with a ratcheting guard so the question
reopens by itself if the data ever changes. `docs/CONCLUSION.md` states all of
it on one page with every number traceable to the code.

**What remains is engineering and communication, not science.**

**One untidy thing worth knowing.** Commit `835b079` is a scratch `wip: cold
clone fixes` commit, made so that `cold_clone_check.py` — which clones the
repository — could see fixes that were not yet committed. It is real work and it
is pushed; it simply has a worse message than it deserves.

**If the GUI test finds something**, the most likely recent culprits are the two
things Round 74 touched in the interface's vicinity: `View → Ion flux animation`
(Round 33, added late) and the result-window provenance stamp (Round 50). Both
have tests, but `scripts/screenshot_app.py` has caught what tests did not twice
before.


## Round 74 — automating the cold-clone run, which immediately caught me

**Five rounds were already done and I verified each rather than ticking them.**
69 and 71 by Round 65 (hybrid implemented, four rows deleted, `dead_code.audit()`
at zero), 70 by Round 60, 72 by Round 67, 73 by Round 63. Round 74 was genuinely
open.

**The argument for automating it turned out to be stronger than I expected.**
Round 60 found eight tests that failed instead of skipping on an empty clone.
Fourteen rounds later, `scripts/cold_clone_check.py` found **two more on its
first run** — both written by me, in Rounds 65 and 67. The same class of bug
came back within five rounds of being fixed, which is exactly the case for a
command rather than a habit.

**One of the two was better fixed in the code than the test.**
`build_hybrid_model` loaded a 2.4 MB AlphaFold model *before* checking the
structure it had been handed, so a caller's own error surfaced as a missing
file. Validating the argument first is better regardless of clones: a large read
should not stand between a caller and the mistake in their input.

**The checker needed calibrating twice, which is becoming the pattern.** It
first reported a clean run as broken: with `-q` and no failures pytest's final
line is a progress bar, so reading the tail found zero of everything. Fixing
that, the counts stayed at zero — because `pytest.ini` already sets `-q`, and
passing another made it `-qq`, which suppresses the summary line however it is
parsed. The exit code is now the signal and the counts are decoration.

That is the tenth time an instrument built to check something in this project
was itself the thing at fault. The tests for it plant a failing run, a clean
run, and an all-skipped run, because "all skipped" must count as *success* on an
empty clone and that is the assertion most likely to be written backwards.

**What it does not do.** The `--fetch` path exists and is not run by default:
it needs the network, and a check that fails when a server is busy would train
people to ignore it. Round 60's finding was that the *fetch* was broken, so that
half still needs running deliberately.


## Round 67 — the methods note, organised around what each safeguard caught

**Round 66 was already done, and I checked rather than assumed.** Block O
proposed the conclusion document before Round 56 reached it. Verified: the page
exists, states the unprovability result rather than listing nulls, carries both
figures, gives 134 against 59, and its traceability guard passes. Ticked with a
note. Round 68 is likewise superseded by Round 63.

**The note's organising principle was the decision worth making.** The obvious
structure is a list of mechanisms with what each is *for*. That reads as a
manifesto and is unfalsifiable. So every section says what the mechanism
**caught** instead — the three-clause decision rule because Round 41's p was
0.0477 with an interval spanning zero; negative controls because in two rounds a
deliberately meaningless predictor matched or beat every real endpoint;
feasibility because "we need more data" was measurably the wrong conclusion.

A safeguard nobody has watched fire is a claim, not a method.

**Leading with the failure is the point, not a disclaimer**, and the note argues
that explicitly: a pipeline that only ever confirmed things would tell you
nothing about whether its safeguards work. This one produced five nulls and then
showed the claim cannot be settled, which is exactly the condition under which
the machinery is worth copying.

**The section I found hardest to write honestly** is the calibration one,
because the list of incidents is a list of my own errors: a spheroid fitter that
would have reported 89% model error, a document checker that could not read its
own documents' minus sign, two dead-code detectors that would have deleted the
CLI, an audit that missed calibrations named in test names, a graft anchored on
a whole flexible arm. Nine of them. Writing that as the strongest section rather
than the most embarrassing one is the correct framing, and it took a moment to
see it.

**Guarded like the conclusion.** A methods note that drifted from the code would
mislead precisely the reader most likely to act on it, so a test resolves every
module it cites and imports it, checks the quoted figures still match `CLAIMS`
and `paired_positions_required`, and asserts the calibration incidents are
recorded in the session log rather than remembered.

**Two of my own test's assumptions were wrong**, both about notation rather than
facts. It searched for `analysis/x.py` at the repository root when that is the
project's shorthand for `piezo1/analysis/x.py`, and it looked for "unicode
minus" where the log says "U+2212". Both failed on wording rather than on
substance — the same failure mode as Round 55's planted probe finding itself.


## Round 65 — one module implemented, four rows deleted, and a real modelling choice

**All five had been promises, not stubs** — no file existed for any of them. The
decision was one-per-module rather than a policy.

`physics/modes.py` was deleted on the evidence of its own row, which said the
useful parts already live on `ModeSet`. `analysis/contacts.py` and
`analysis/variants.py` went because `interactions`, `allostery`, `annotations`,
`variant_sets` and `variant_impact` already cover what they promised.
`analysis/docking.py` went for a better reason: it contradicts the project's own
stance. `build_ligands.py` *refuses* a `bound_structure` claim and verifies
against the downloads that no bound modulator exists — shipping a docking module
would generate poses the ligand file declines to accept as evidence.

**`structure/hybrid.py` was implemented because nothing else serves aim A1.**
The full-length model — cryo-EM core plus predicted distal blade, seam visible
— is a stated project aim that had gone fifty rounds without an implementation.

**What it reports about itself is the point.** 8YEZ resolves 570–2521, so the
graft is 569 residues of distal blade. Mean pLDDT there is 64.5 with only 48%
above the conventional 70: the prediction is least confident exactly where it is
being relied on. Every atom carries `source` and `plddt` so no analysis can
average across the join silently.

**The modelling decision was measured, not assumed.** My first version anchored
the graft on the whole 1279-residue overlap and reported 19.0 Å RMSD. That is a
plausible number for a flexible multi-domain protein, and I nearly recorded it
as the finding. Anchoring near the seam instead fits to **2.4 Å** — because what
determines where a grafted piece attaches is the local geometry, not the average
over an arm whose conformation differs. Then the far end of the blade sits 75 Å
from the experiment, and that is reported too, because a good local fit hides
exactly that.

**A bug I found by not trusting the first number.** The 19 Å looked large enough
to check, and the check found that `build_hybrid_model` fell back to the whole
trimer when `Structure.select` did not exist — building the residue map over
three chains and keeping whichever came last. That is the `dict(zip(...))`
pattern from an earlier round. Fixing it did **not** change the 19 Å, which I
should say plainly: the bug was real, the number was not caused by it, and I had
assumed otherwise for a few minutes.

**Two error pages were living in `ref/structures/`** — 127-byte XML `NoSuchKey`
responses named as AlphaFold v4 models, from before `_download` gained its size
guard. Nothing read them, but a glob over the structure directory would have.

**And the provenance count is now explained rather than displayed.** 7 of 21
chains complete reads as failure and is not: every break is a claim computed
from a frozen record consuming no parameter, or an analytic claim reading no
structure. Document breaks — the kind that would mean something — are zero, and
a test asserts every break stays benign.


## Round 64 — declining to pre-register, and recording that as a result

**The item was written conditionally and the condition failed.** "Only if
Rounds 61–63 leave a design with adequate power." They leave 8 positions
required at an implausibly good predictor against 1 available, 3–4 reachable,
and 0 added by the engineered variants. So nothing was pre-registered.

**Writing that down was the round.** The easy version of this is to tick the box
with "not applicable" and move on. But a decision not to test is a decision, and
this one rests on numbers from three rounds that a later reader would have to
reassemble. `NOT_PREREGISTERED_ROUND64.md` states them once, in the same shape
as a null result — what was proposed, why it will not be run, and what would
change it.

**The harder question was whether to run it exploratorily.** §2 of the protocol
explicitly permits exploratory work if it is labelled. So the refusal needed a
reason beyond "underpowered", and there are three.

A sign test on one pair has a minimum one-sided p of **0.5** — I checked rather
than asserting, and even four perfect pairs give only 0.0625. There is no
outcome. Second, the one available position is R2456, which every round since 7
has cited as the example that *breaks* the predictor: the informal answer is
already known, so a pre-registration written afterwards would be a formality
dressed as a safeguard. Third, a δ from one pair would be quoted somewhere
without its caveat, which this project has watched happen enough times to treat
as a prediction rather than a worry.

**A refusal is worth nothing unless something makes it stick.** The test does
three things: re-checks the arithmetic, asserts no within-position comparison
exists in the codebase, and — the part that matters — **ratchets the count of
discriminating positions**. If a new one ever appears the suite fails and points
at the document. The question reopens by itself rather than depending on a
future round finding a file.

**And my own guard had a false positive, caught by calibrating it.** The
codebase check flagged `feasibility.py`, which simulates a sign test to compute
the *required* sample size and imports the discriminating positions to compare
against. That is the design analysis, not the comparison — the two look alike to
a keyword scan and are opposites. Sharpened to require variant scores as well.
That is the seventh time a checking instrument in this project has caught itself
first.


## Round 63 — a scientific decision, settled by the project's own annotations

**This was the one remaining question that was not about counting.** Fifteen
engineered variants carry measured functional effects and no analysis set uses
them. Round 54 marked them `blocked` on a genuine scientific question — may a
conductance or selectivity change stand for gain or loss of mechanosensitive
function? — and deliberately declined to answer it while counting routes,
because a costing exercise that answers a scientific question by arithmetic is
how bad decisions get made quietly.

**The answer is no, and the argument did not have to be general.** I expected to
reason from the disease phenotypes: hereditary xerocytosis variants show slowed
inactivation and larger mechanically-evoked currents, lymphatic dysplasia
variants show reduced ones, and neither is a statement about unitary
conductance. That argument is sound but it is the kind a reader has to take on
trust.

Two of the fifteen make it concrete instead. **A2078W**: "Yoda1 sensitivity
severely reduced *while mechanosensitivity to stretch is retained*".
**KKKK2166-**: "selectively removes inactivation *without changing mechanical
sensitivity*". Both dissociate the axes at a single residue, in the project's
own curated annotations. So the refusal rests on measurements already in the
file rather than on a principle I am asserting.

**Refusing is not refusing all fifteen, and that distinction mattered.** Five
are on the right axis — S1335A and S1335V raise the force threshold, A1718W
loses stretch-activated currents, P2113A desensitises mechanically, S2446E
stabilises an open intermediate. A blanket exclusion would have been easier to
write and would have thrown away real evidence.

**Then the measurement settled it.** Admitting all five adds **zero**
discriminating positions. None sits at a position carrying a directional curated
variant. Position 1335 does hold the only engineered pair — and both variants
*raise* the threshold, so it is same-direction, and a within-position test needs
the two to disagree. That near miss is pinned, because it is exactly the kind of
thing a later round could misread as a usable pair.

**A test I nearly wrote badly.** The first version checked each verdict's axis
against my own classification, which is circular — I wrote both. It now checks
that each verdict's stated basis is a real phrase from the curated
`functional_effect`, so a loose paraphrase cannot put a variant on whichever
axis the argument needs.


## Round 62 — the count was already known; the evidence level was not

**Most of this item was answered two rounds early.** It asks how many of the
"40 positions" carry directions — and Round 54 had already established the 40
was wrong and the answer is one. Repeating that would have been busywork.

**What was genuinely unmeasured is the evidence-level split**, and it is not a
detail. `variant_sets` refuses to pool `measured` electrophysiology with
`disease_mechanism` inference, deliberately, so a within-position design has a
*different ceiling at each level* — and only the pooled number had ever been
counted.

**One finding is reassuring.** R2456's four variants are all at `measured`
evidence. The single discriminating position the project has is not a weak one
that a stricter design would throw away; it is its best-evidenced position. I
had half expected the opposite — that the one usable site would turn out to
rest on inference — and it was worth checking rather than assuming.

**The other is not.** A pair is no stronger than its weaker half. M870I is
`disease_mechanism`, so curating M870V — one of Round 54's three named targets —
can never produce a `measured`-level pair. At the level a confirmatory test
would actually need, the ceiling is **3**, not 4.

**So the answer to "does this collapse the Round 61 design" is no: it confirms
it and makes it slightly worse.** Round 61 measured that even a predictor
ordering nine pairs in ten correctly needs 8 shared positions. The best
reachable count is 4 pooled, 3 on measured evidence alone. No level reaches the
requirement.

That closes the last open question about the central claim. Both designs —
across positions and within them — are now costed and closed at every evidence
level, by measurement rather than by discouragement.


## Round 61 — costing the design that was supposed to be the way out

**The premise was already dead, and the question still worth asking.** I wrote
this item after Round 50, when I believed 40 shared positions existed. Round 54
corrected that to one. But "how many would be enough" does not require having
them, and it is the question that decides whether the route Round 47 left open
is a route at all.

**Choosing the statistic was the substantive decision.** A within-position
comparison has no distribution to estimate — there is no sample. So the sign
test is the right instrument: at each position carrying both directions, does
the predictor rank the gain-of-function variant above the loss-of-function one?
Under the null that is a coin flip, and the requirement follows from the
binomial alone. Cliff's delta on a paired ordering is 2p − 1, which keeps the
answer on the same effect scale the other rounds use.

**The answer closes the route.** At the across-position effect (δ = 0.249) the
paired design needs about **102** positions. At δ = 0.5, twenty-six. Even at
δ = 0.8 — the predictor ordering nine pairs in ten correctly, which nothing in
five rounds suggests it can do — it needs **eight**. There is **one**, and Round
54 put the absolute ceiling at four if three named variants could be directed.

**The part that makes this decisive rather than merely discouraging.** Pairing
is *not cheaper at the same effect*: 102 positions against Round 47's 134
variants is the same order. The entire argument for pairing was that it would
**enlarge** the effect by removing the between-position variance that consumed
99.8% of Round 7's predictor. So the route only ever paid off if the paired
effect were much larger — and it needs eight positions even then.

**Calibrated before believed, as the standing rule requires.** The instrument is
a simulation, so it is checked at both ends first: a coin-flip predictor must
return "not detectable at any sample size" rather than the search bound (the
failure mode the kinetics calibration was fixed for), and a perfect predictor
must need a handful. Only then is the middle of the curve trusted.

**And no comparison was run**, as in Round 47. A test asserts the module imports
no statistic that could produce one — because the temptation here is obvious:
there is one discriminating position, R2456, and looking at how the predictor
orders it would be a one-line answer to a question that is not pre-registered.

**One thing the parameter audit improved.** I first wrote the available count as
`SHARED_POSITIONS_AVAILABLE = 1`. The audit flagged the literal, and the right
fix was not to register it but to *derive* it from
`data_routes.discriminating_positions()` — so the requirement and the supply
cannot disagree. A test had already been asserting they matched, which is a sign
the constant should not have existed.


## Round 60 — an empty clone, and the three things a warm cache was hiding

**The point of the exercise is that none of these could be seen here.** Every
defect this round found is invisible on a machine that already has the data,
because the failing code path is never taken. That is the argument for doing it
at all, and it produced three.

**Bug 1: eight tests failed instead of skipping.** `conftest.py` states the
project's rule — skip when data is absent, do not fail — and eight tests did not
follow it. Three of the eight I wrote in the last fifteen rounds. They now skip
*and name the fetch command*, because a skip that does not say what to run
leaves whoever hits it no better off.

**Bug 2: the Ensembl CDS download had been broken.** The URL sent the content
type as a `;content-type=text/x-fasta` query parameter. Ensembl used to honour
that and no longer does — the plain URL answers **415**, and the project's
downloader reported it as a 500. It must be an HTTP header. Nobody with the
file on disk would ever see this, and the file has been on disk here for many
rounds. `_download` now takes headers and both species fetch.

While fixing it I hit a 503 on the second request and briefly took it for the
same fault. It was rate limiting: a retry succeeded immediately. Worth
recording because the two failures look identical in a log and mean opposite
things — one is a broken integration, the other is a busy server.

**Bug 3 is the one that matters, and it is not a crash.** `feasibility.assess()`
calls the literature harvest to compute its ceiling. With no corpus the harvest
returns nothing, so the ceiling fell from **59 to 34** — and the function said
nothing. A reader on a fresh clone would have got a *stronger* conclusion than
the data supports, out of the same code, with no indication. 59 is a documented
number that appears in `CONCLUSION.md`, the README and the tour.

That is precisely "anything that only works because of a stale cache", and it is
worse than the two crashes because it produces an answer. It now reports
`harvest_available = False` and refuses to state a ceiling at all.

**What was not a bug.** The 110 skips on a populated clone are the recorded
validation runs — `data/derived/*.json`, regenerable by their scripts — and the
optional PAE download. `verify_claims` reports those four claims as *skipped,
not drifted*, which is the behaviour `test_reproducibility.py` already pins. The
distinction matters: drift means the code changed, skipped means the input is
not here, and conflating them would make a fresh clone look broken.

**And the test I wrote to pin bug 2 tripped on its own documentation** — the
scan for the old `;content-type=` form found it in the comment explaining the
fix. Comments are stripped now. That is the third time in six rounds that a
checking instrument has caught itself rather than the code, which is at least a
reassuring failure mode.


## Round 59 — the README ended on its licence, not on its result

**The item was stale twice over, and checking which halves were real was most
of the work.** It said the claim "has been tested four times" — it is five —
and it asked for the tour's closing steps to be rewritten, which Round 53
already did. Rather than assume, I checked: no tour step before the ending
presents the claim as open, and the variant step already frames R2456 as *why
this is hard*. That half was genuinely done.

**The README half was not.** Round 56 put a link to the conclusion at the top,
which is where a reader starts — but the document *ended* on data sources and a
licence. Its last scientific statement was "The result that says it works",
which is the structural validation: the mode overlap and the dome curvature.
Both are true and neither is the project's answer. A reader who read to the end
would finish on the part that worked.

**So the README now ends where the science does**: five pre-registered tests
across five predictor families, five nulls, the forest plot — and then the part
that makes it more than a list of failures. 134 variants needed against a
ceiling of 59 across positions; exactly one usable site within them. Therefore a
sixth test should not be run, whatever predictor goes into it. It closes on what
is reusable, which is the apparatus rather than the predictor.

**Guarded, because three surfaces stating one record is exactly how drift
happens.** The conclusion page's traceability guard now covers the README
summary too — every number must come from the claims registry, the validation
record or the published-interval table. And a new cross-surface test asserts the
tour, the README and `CONCLUSION.md` agree on the *count*, so a sixth test
updates all three or fails.

That test exists because these surfaces have already gone out of step. The tour
said "tested twice" after five tests had run, and it survived for several rounds
because each surface was hand-written at a different time and nobody compared
them. Linking them by a test is cheaper than remembering to.

**One thing I deliberately left alone.** "The result that says it works"
overstates nothing — the elastic network really does find the gating transition
through the symmetry channel theory permits, and the dome really does reproduce
the published curvature. Softening a section that is accurate, because a
different claim failed, would be its own kind of dishonesty.


## Round 58 — the illegitimate reading was written into the API

**What I expected to find, and what was actually there.** Round 39 recorded
that the mechanical score has a legitimate use — locating mechanically coupled
positions — and an illegitimate one. I expected to be renaming something
mildly suggestive. The output was a class called `VariantPrediction` carrying a
property called `direction`, whose docstring read *"Predicted direction:
stiffening (LoF-like) or softening (GoF-like)"*. The claim five pre-registered
tests failed to support was not merely implied by the naming; it was the
docstring of a public property.

**So the rename is the deliverable, not cosmetics.** `CouplingScore` with a
`gating_cost_change` field and a `sign` property that says explicitly it is not
a gain/loss mapping, and names R2456 — where H, K and P are gain-of-function, C
is loss, and this model gives all four nearly the same number — as the
demonstration. `ddG` went too: it implies a thermodynamic free energy of
folding, which this quantity is not.

**No alias.** The tempting move is `VariantPrediction = CouplingScore` for
compatibility. An alias is exactly how the old reading comes back, and a test
asserts the old name resolves nowhere.

**A correction to my own Round 50 register.** The hazard
`prediction_read_as_validated` described a user selecting a variant and *seeing
a mechanical ΔΔG*. Checking where the score is actually computed: nowhere in
the GUI, nowhere in the CLI. It is reachable only from a notebook or the
validation scripts. The registered scenario could not happen as written — I had
written down a plausible hazard rather than a measured one, which is the same
failure as the Round 50 review's forty positions. Corrected to the real,
narrower exposure rather than deleted, because a notebook user does still reach
it.

**The frozen record needed the opposite treatment.** `run_validation.py`
regenerates the frozen Round 7 result, whose JSON carries the keys `ddg` and
`ddg_normalised`. The attributes had to be renamed; the stored keys had to not
be, or the script would stop reproducing the record it exists for. Both halves
are pinned: the stored keys stay, and no attribute read uses an old name. It is
a small thing, but renaming through a codebase is exactly where a frozen record
gets quietly rewritten.

**What is kept.** The coupling map — PRS gate response, PRS coupling, dynamic
cross-correlation to the gate, betweenness — is untouched and is the part with
a defensible use. A test asserts the score still computes, because a retirement
that broke the calculation would be a deletion wearing a rename.


## Round 57 — reading all 35 by hand, and what that was worth

**The round asked a precise question and it has a precise answer: five.** Round
45 extracted 35 substitutions the curated set does not have, each with its
sentence. How many carry a direction a human can recover from the sentence
alone? I read all 35. Five do.

**And the five are worth less than five.** Every one is an alanine-scanning
mutant whose sentence says "non-functional" — D1975A, D2034A, L2131A, R2135A,
W2140A in the source numbering. That is loss of *channel* function in a
mutagenesis screen, not the loss-of-function-in-disease the curated set records;
admitting them at the same evidence level would be exactly the pooling
`variant_sets` refuses. More decisively: **none of the five sits at a position
carrying any other variant**, so they create no within-position pair and Round
54's count of one usable position is untouched.

**The other thirty are informative about why the corpus does not help.** Four
give a direction only for chemical agonist response, and two of those are
double mutants whose single-mutant phenotype is never separated. Five are
conductance changes, which Round 54 already marked as a different question.
Seventeen have no phenotype in the sentence at all — construct lists, figure
legends listing plot symbols, or statements that some mutations are "scattered
throughout the channel". Three are a clone's sequencing differences from the
reference rather than mutants anyone tested.

**Two faults in the harvest, found only by reading.** V190P is a **STOML3**
mutation from a paper studying STOML3 and PIEZO1 together. The wild-type gate —
which rejects 23% of raw hits and is the reason to trust the rest — passed it
because position 190 is valine in PIEZO1 as well. No residue-identity check can
catch a substitution that is real but belongs to another protein, and the honest
response is to record the class of error rather than filter this one instance
away.

The second: Round 45 reported two candidates carrying a measurement. Read in
context, both are fragments of a conductance list — `'7 pS, V2132A; 59.'` is the
tail of something like "56.7 pS". So the number of harvested candidates with a
usable measurement is **zero, not two**, and the small encouraging figure in the
Round 45 record was a parsing artefact.

**What this closes.** The harvest was the last route that could add *measured*
directions without new experiments. It is now spent, and the answer is recorded
as data with the phrase each verdict rests on, so nobody has to read 35
sentences again to find out.


## Round 56 — the conclusion, and a guard so it cannot go stale

**The item was written before the result it now states existed.** It asked for a
page collecting "the four nulls, the data limits and the model-error result".
There are five nulls, and more importantly Rounds 47 and 54 changed what the
page can claim: not "the central claim is unproven" but **"the central claim
cannot be settled with data that could exist"**, established by two independent
measured routes — 134 variants needed against a reachable 59 across positions,
and exactly one usable position within them.

That is a result rather than an apology, and the document is organised around
it: what was established, what was not, and then why the second is a finding.

**Why the page is guarded rather than proof-read.** This project has shipped
stale prose twice in the last ten rounds — the tour still said "tested twice"
after five tests had run, and my own Round 50 review counted 40 usable
within-position sites where there is one. Both survived because nobody
re-derived them. A page whose whole purpose is to be trusted without reading
the working is exactly where that failure would matter most.

So `tests/test_conclusion.py` extracts every number from the document and
requires each to be supported by the claims registry, the validation record,
the published-interval table, or an allowlist in which every entry states why
it is exempt. It caught two numbers on its first run — the derived agreement
percentages 0.4% and 0.01% — which are legitimate but are arithmetic on two
sourced values, so they are allowed *with the derivation written down* rather
than silently permitted.

**And the guard is calibrated**, per the standing rule: a test asserts that an
invented number would fail the check. A guard that cannot fail is a rubber
stamp, and this one would have been easy to write that way — the allowlist is
generous enough that a lazy version would pass everything.

**One thing I deliberately did not do.** The obvious way to satisfy "every
number traceable" is to generate the document from the registries. I wrote it
by hand instead, because a generated page would be a table rather than an
argument, and the argument — that two independent routes are closed, so a sixth
test should not be run — is the part a reader needs. The guard gets the safety
without giving up the prose.


## Round 33 (completed) — a rate made visible, and the label that makes it honest

**Why this was left undone and why it was worth finishing.** Round 33 built the
permeation physics and deferred the animation with a one-line note that the
morph clock's discipline applies. It does, and inverted: the morph refuses a
seconds axis because an interpolation between two endpoints is not a
trajectory, whereas a current genuinely *is* a rate — so a time base here is
meaningful and the honest move is to name the factor rather than avoid it.

**The factor is what makes the label non-optional.** A single channel passes
1.5 × 10⁷ ions per second. At 60 fps that is a quarter of a million ions per
frame, so anything a person can follow runs about a millionfold slow. An
unlabelled stream of particles reads as "this is what it looks like", which is
wrong by six orders of magnitude — the same class of error as quoting a
confidence interval where a model spread dominates.

**A shut pore animates nothing, and that is a first-class outcome.** The
permeation result is gated by the wetting verdict and every deposited human
structure is closed. Measured: 8YEZ draws no ions and shows why — score 0.82
against a 0.55 cutoff, bottleneck 0.095 nm, occluded *and* hydrophobically
gated. 11ZC gives 2.44 pA. Drawing a trickle through a closed gate would
contradict the project's own structural result while looking like a
demonstration of it.

**And the rate declares what it inherits.** The solver gives 41 pS against a
published 25–30 — a disagreement this project already records and deliberately
does not tune away. So the animation runs about 1.5× fast, and the HUD says
that as well. A stream calibrated to a number 1.5× too large, shown without
saying so, is precisely the confident-wrong-picture failure Round 50 audited
for; having just built that audit, shipping one would have been poor.

**Two faults in my own work, and the second is the useful one.** The controller
first read `result.current_pA`, a field that does not exist — and
`PermeationResult.current` is in **amperes**, so had the name been right and
the units unconverted, the rate would have been wrong by 10¹² and would have
looked like a plausible animation either way.

The test meant to cover this inspected the controller's *source* for
`predict_wetting`. It passed the entire time the units were wrong, and then
failed the moment I moved the gating into a function where it could actually be
run. A source-inspection test proves reachability, not correctness, and it
fails for the wrong reason exactly when the code improves. The physics now
lives in `render.flux.timebase_for_structure`, Qt-free and exercised on real
8YEZ and 11ZC coordinates.

**Also closed: a stale checkbox from Round 1.** "Expose the pore profile in the
GUI" was deferred, delivered in a later round, and never ticked — the Analysis
dock has had a two-axis pore plot with click-to-locate for some time. Verified
rather than assumed before ticking it.


## Round 55 — retiring what earns nothing, and three attempts at the detector

**The answer to the question asked was "almost nothing", which is worth
recording.** The round assumed the codebase had accumulated scaffolding to
delete. Measured across 109 modules and 568 top-level definitions, every module
is imported, tested or a documented entry point; at definition level exactly
five things earned nothing. So the project is already load-bearing, and the
round's value is the standing guard rather than the deletions.

**The hard part was the detector, and two versions of it would have done real
damage.** The first was a grep over `__all__`: 102 unused public names,
including `format_result` (used inside its own module) and every return-type
dataclass, which is constructed but never named elsewhere. The second used the
AST but collapsed same-file references into a *set*, so a helper called only
from its own module looked unreferenced — 129 dead functions including
`fetch_pdb`, `cmd_list` and `_optimise_slice`. Acting on that list would have
deleted the CLI.

Round 51's rule is the only reason neither was believed. Both outputs were
plausible: a long list of unfamiliar names is exactly what a real finding would
look like, and nothing about it announces itself as wrong. What caught it was
checking the list against names I *knew* were used before reading the rest.

**So the calibration is now part of the instrument.** `dead_code.calibration()`
requires that nine known-used names are not flagged and that a planted
unreferenced name is, and every test runs it before touching the audit's
output. The fix to the counting was to count *occurrences* rather than files,
and to count bare words inside string literals — this project dispatches by
string through `ANALYSES` and the CLI subcommands, so a name reached only that
way is live.

**The calibration probe defeated itself first.** The planted name was written
as a literal, so the string scanner found it in the detector's own source and
reported the instrument as broken. That is the string scanning working exactly
as designed; the probe now builds the token from fragments so it appears
nowhere. A small thing, but it is the third time this round that the checking
apparatus was the thing at fault.

**One deletion needed permission from the record.** `_poisson_newton_step` was
42 lines of dead code whose docstring carried a real numerical finding — the
Gummel loop going −0.37 V, then −171 V, then −2×10¹⁶ V in three iterations,
which is why the solver uses the electroneutral limit. I checked that the
divergence is recorded in the module, in `test_permeation.py` and in
`SCIENCE.md` before removing it. Deleting the last copy of a finding would be
worse than keeping dead code, and a test now pins all three records.


## Round 54 — costing the data limit, and correcting my own review

**The item was stale, like Round 53's.** It named two routes needing no new
experiments: gnomAD constraint and published supplementary tables. Both had
already been run — Round 41 returned null with its negative control
indistinguishable from the predictor, and Round 45 found 35 candidates of which
none carries a direction and only two have any measurement. So "cost them and
do the cheaper one" had no cheaper one left to do, and the useful work was
costing what remains.

**The real result is a correction to something I wrote last round.** The Round
50 review counted 40 positions carrying more than one variant and called that
"a real design, not a curiosity". That was the most encouraging number this
project had produced in a while, and it was wrong for the purpose. A
within-position comparison needs two or more **missense** variants at a
position, each carrying a **direction**, from sources that do not **disagree**.
Applying those filters gives **one** position — R2456 — which is exactly what
Round 48 measured before I talked myself past it.

What the 40 contained: nonsense variants (Q1009\*, a truncation, not a
substitution a structural predictor can score), insertions (E2496ELE), positions
whose second variant has no direction at all, and V598M, which is curated as
gain-of-function and inferred from ClinVar as loss-of-function. That last one is
worth noting for the right reason: `variant_sets.disagreements()` was **already
reporting it**. The machinery was correct and the review was not — I had
counted rows without asking what each row was.

**"Not enough variants" is now a named list.** Exactly three variants would each
unlock one more position: M870V (the position already has M870I, LoF), R1358C
(has R1358P, GoF) and A2020V (has A2020T, GoF). That is the cheapest remaining
route and the only one that is a finite list rather than a search. It is also an
upper bound — two of the three are curated as VUS *precisely because* the
evidence to direct them was not found, so the expected yield is below three and
the ceiling is four positions. Four positions is not a design.

**Why the engineered variants are marked blocked rather than open.** Fifteen
carry a measured functional effect, which is the only untapped measured data in
the project. But the effects are changes in conductance and selectivity, not in
gain or loss of mechanosensitive function. Whether one may stand for the other
is a scientific question and belongs in its own round with its own reasoning,
not in a costing exercise that would quietly answer it by counting.

**What this closes.** Round 47 established the across-position route needs 134
variants against a ceiling of 59. Round 54 establishes the within-position route
needs positions the data does not contain. Both are now recorded as measured
costs rather than impressions, so the next round to reach for either has a
number to argue with.


## Round 53 — ending the tour on the record, and finding the roadmap item stale

**The item was out of date, which is the finding.** It asked for the closing
steps to be rewritten because the tour ended on *two* nulls when there were
*three*. There are now **five**, and the tour's final step still said "it has
been tested twice, both times pre-registered". A student taking the tour was
being told the project had made two attempts on its central claim when it had
made five, all null. That is the exact failure the round exists to fix, and it
had happened to the fix itself.

**So the rewrite is structural, not editorial.** The closing measurements now
read `ALL_PREREGISTERED` and the claims registry rather than restating numbers,
so a sixth test updates the tour automatically. A test asserts that
`_data_limit`'s source contains neither 134 nor 59 — the numbers must arrive
from the registry, because prose is what went stale last time.

**Three closing steps instead of one.** The record (five tests, five
predictors, five nulls), the data limit (134 needed against a ceiling of 59),
and what remains uncertain even where the measurements worked (the dome's
shape spread, ~18× its bootstrap interval). Splitting them matters because they
say different things: the first is "the claim is unsupported", the second is
"and it cannot be supported with reachable data", the third is "and even the
successes are less precise than they look".

**The sharpest line is Round 48's.** A feature computed on the wild-type
structure has exactly zero within-position variance, so R2456H, R2456K, R2456P
and R2456C receive the identical value to every digit. That is worth a student's
attention more than any p-value: it says the *kind* of predictor cannot answer
the question, independent of how well it is fitted.

**Where I kept the scoping rather than simplifying it.** The obvious move was to
append Rounds 41 and 48 to `VALIDATION_RECORD`. I did not: that record is what
the GUI shows beside a variant's ΔΔG score, and pooling in tests of population
constraint and wild-type context would make the caveat describe a number the
user is not looking at. `OTHER_PREREGISTERED` holds those, `ALL_PREREGISTERED`
joins them, and the tour reads the joined view.

**Figures, and why they must be allowed to be absent.** `TourStep` gained an
`image` field, and `body_html()` appends the figure only when the PNG exists.
`docs/img` is regenerable and git-ignored, so a fresh clone would otherwise get
a tour that raised or showed a broken image. Degrading to prose is the correct
behaviour and is pinned. Both figures took two layout passes — the first put the
title through the "no effect" label, the second hid the ceiling annotation
behind a bar label — which is the usual cost of not looking at what you made.


## Round 52 — publishing the interval that answers the question actually asked

**The problem, restated.** Round 38 measured that the dome radius's model
spread is six times its bootstrap interval. Nothing about that was wrong: the
bootstrap correctly reports how well a sphere is determined by 66 surface
points. It just is not the limiting uncertainty, because the open question is
whether a sphere is the right shape. Publishing ±0.9 nm was answering a
question nobody asked.

**The decision, as a rule rather than four judgements.** Publish the *widest*
term and name its kind in the same breath; never call a sensitivity range or a
model spread a confidence interval, because a network cutoff has no sampling
distribution; and where the widest term is a model spread, say it is a **lower
bound**, since two models agreeing does not bound the error from above. Written
as `analysis/published_interval.py` so the rule and the numbers live together
and a fifth quantity cannot be added without stating its terms.

**T₅₀ was the surprise, and it is the most useful result of the round.** It
looked well determined — the matrix exponential and an adaptive ODE integration
agree to 0.6% — so the natural thing to quote is the solver agreement. But
perturbing the published Young 2023 rate constants by ±20% moves it to
[2.584, 2.838], **sixteen times wider**. The number is limited by its inputs,
not its arithmetic, and quoting the solver agreement would have implied
precision it does not have.

This also makes the headline agreement *stronger* rather than weaker. The
measured 2.7 ± 0.1 mN/m lies inside the input-propagated range, so the match
with Lewis & Grandl survives the uncertainty on the rate constants instead of
depending on their exact published values. A test pins that, because it is the
part that would quietly stop being true if a rate were revised.

**A mismatch found on the way, kept rather than repaired silently.** The dome's
model comparison is anchored on the **untrimmed** sphere fit (9.45 nm) while
the published number is **trimmed** (9.72 nm) — the two were never
like-for-like, and nobody had noticed. Before deciding whether it mattered I
measured the trim's own effect: `geometry.sphere_trim` moves the radius by
0.30 nm across 0–0.25, against a 5.54 nm model spread. So the conclusion is
untouched, and the honest thing is to record the mismatch as its own term
rather than delete it by re-anchoring the comparison and saying nothing.

**Why `verify_claims` is untouched, and why that is the right answer.** The
roadmap asked that the documented numbers and their stated uncertainties move
together. They did — by the point estimates deliberately not moving at all.
What changed is only what is claimed *about* them. A claim tolerance exists to
detect code drift; a published interval is a scientific statement about what
the measurement can support. Conflating them would mean widening drift
detection every time an honest uncertainty grew, which would make the drift
guard weaker exactly when the science got less certain.


## Round 51 — calibrating the checkers, including the one doing the audit

**Why this round exists.** Four of the last five rounds found their defect in
the *instrument* rather than the science, and always the same way: an
alternative route built to check the main one was itself wrong, and returned a
plausible number instead of an error. The disagreement then looks like a
finding. This round turns that pattern into a rule and a mechanism.

**The audit's own instrument was uncalibrated, which is the round in
miniature.** My first pass was a keyword scan over the test files, and it
reported twelve instruments as having no calibration. Reading the tests showed
why: it searched the test *bodies*, and several calibrations are named in the
test *name* — `test_auroc_known_cases`, `test_cliffs_delta_extremes`. Corrected,
the true count was four. Had I believed the first run, this round would have
"fixed" eight things that were never broken, and the roadmap would now record
eight false findings. I noticed only because the claim "the statistical
instruments are uncalibrated" contradicted what `INTERFACE.md` said, and
checking which was right is exactly the habit the rule now demands.

**The four genuine gaps, and why each was a gap rather than an oversight.**
`spring_model_error` was checked only by asserting its spread on real
structures was "modest" — a statement about the answer, not the instrument; a
version computing overlaps wrongly would have produced a modest spread too.
`minimum_detectable_effect` was checked only against this project's own
recorded results, which is circular, since the same function supplies the power
statements those rounds are judged by. `sensitivity` was tested on its
*wording* — that it refuses to call itself a confidence interval — and never on
its arithmetic. `permutation_test` had a null case and a real-shift case but no
exactly-known value, though for eight observations all seventy partitions can
simply be enumerated.

**Each calibration has to be able to fail.** A check that would pass on a broken
instrument asserts nothing, so the register carries a test that the new
calibrations reject a degenerate instrument — a `sensitivity` returning the
reference at both ends, or a permutation test calling identical samples
significant.

**A reproducibility finding, and the check that stopped it becoming a scare.**
The spring calibration was flaky, and the cause was not test order: `calc_modes`
returns overlaps between 0.954 and 0.997 on *identical* inputs, because ARPACK
starts from an unseeded random vector and the random test geometry has
near-degenerate low modes. The immediate worry was `anm.gating_overlap`, which
is a documented claim. Measured over four runs it is bit-identical to eight
decimal places — real PIEZO1's low modes are well separated, so the
non-determinism never reaches the science. Reporting "the mode solver is
non-deterministic" without that check would have been a true statement that
implied something false. The test now asserts the separation between spring
models rather than an absolute bound, which is the scientific content anyway.

**Made enforceable rather than aspirational.** `CALIBRATED` maps every public
callable in the eight checking modules to the test that calibrates it;
`test_every_checking_instrument_has_a_calibration` fails when one is added
without a case, and `test_named_calibrating_tests_exist` fails when the named
test does not exist. That second guard immediately caught two test names I had
guessed from memory rather than read — which is the whole argument for having
it.

**The rule is now in `CLAUDE.md`** with the four incidents that motivated it,
including the standing instruction that matters most: when a checker disagrees
with the pipeline, suspect the checker first, because historically it has been
wrong more often than the thing it was checking.


## Round 50 — auditing for wrong numbers rather than missing buttons

**The distinction from Round 33.** That round asked whether every analysis was
*reachable* from the menus. This asks the harder question: given that a user can
reach it, can they be handed a number that is wrong with nothing saying so. A
wrong number that announces itself is a bug; one that looks exactly like a right
one is the failure this project exists to avoid.

**Two of the three named suspects were open.** The cross-species overlay
refusal was real and fires. The other two were not guarded at the point where
the number appears: a result window named **neither** the structure the numbers
came from **nor** the parameter set they were computed under. With companions
displayed there was nothing on the window identifying which structure it was —
precisely the "companion mistaken for the primary" failure — and an overridden
registry produced numbers that looked documented. The status bar did warn about
overrides, but a separate, non-modal window is exactly where that warning fails
to reach.

**The stamp records compute time, not read time.** This mattered enough to test
directly. A stamp recomputed on access would quietly agree with whatever the
registry says later — but the numbers in the window were produced under the
earlier set, so a window that outlives a parameter change must keep saying what
it was computed under and visibly disagree with the status bar.

**A third gap fell out of the audit.** `CAVEATS["interactions"]` was the empty
string — the one tabular analysis shown with no warning at all. Interactions
are measured directly from coordinates, so there is no *modelling* caveat, which
is presumably why it was left blank; but there is a real one, and it is about
state and resolution: the contacts are those of this structure in this
conformation, unresolved side chains cannot contribute a bond, and the criteria
are heavy-atom based because deposited entries carry no hydrogens. A test now
forbids showing any analysis without a caveat.

**Positive controls, not assertions.** Each hazard is exercised by constructing
the dangerous situation and watching the guard fire. That immediately corrected
one of my assumptions: I expected `verify_claims` to return results marked
incomparable, and it *raises* instead — a stronger guard than I had written the
test for, and the test now pins the real behaviour plus the deliberate
`allow_overrides` escape that keeps a refusal from being an obstruction.

**Generated, not duplicated.** The register produces the help topic rather than
being re-written as HTML. A hand-written second copy would have drifted from the
guards the first time one changed — the same reasoning as `prediction_record`
and `claims`.

**The review that came with this round found something the rounds themselves
missed.** Round 48 measured one position carrying more than one variant, which
made the within-position design look impossible. That count was over the 46
directional missense subset. Over the full curated and ClinVar sets it is
**40 positions** (6 curated, 13 ClinVar, 40 combined, 8 with three or more) —
a real design rather than a curiosity, and reachable by curation rather than by
experiments nobody has run. There are also 15 engineered variants, every one
carrying a measured functional effect, excluded from every analysis set because
`engineered` is not `GoF`/`LoF`. Block N is built around both.


## Round 49b — wiring the other twenty-one, and why static proof was not accepted

**The job.** Round 49 found 26 registered parameters that no code read and
wired the five `pore.*` ones. This wired the remaining 21 — 11 modules, 28 call
sites — taking the unwired count to **zero**. Every documented number is
unchanged: 759 tests pass and all 18 claims verify with no drift.

**Why a static check was not enough to close the round.** `provenance_chain`
can show that code *reads* a key. That is not the same as the key mattering: it
can be read into a variable that is shadowed, passed to a function that ignores
it, or resolved on a branch nothing takes. Round 49 had already been misled once
by reading source, so `parameter_effect` settles it the only way that is not an
inference — override, recompute, and see whether the number moved. It checks
both halves, because a parameter that changes the answer but does not restore it
on reset is worse than one that does nothing: it leaves the process unable to
reproduce the documented value.

**The probe immediately caught a defect the wiring introduced.** `value()`
returns a float, so `n_permutations` reached numpy as `10000.0` and was
rejected. Eleven count-valued parameters now cast where they resolve. A static
wiring check would have marked all eleven done — this is the concrete argument
for paying for the empirical version.

**And a probe is worthless until it is calibrated.** Two parameters first read
"no effect" for reasons unrelated to wiring. The pockets probe used coordinates
too diffuse to form a single alpha sphere, so every pockets key looked dead; and
`pockets.r_max` was pushed *upward*, where nothing new can be admitted because
the widest sphere in that geometry is 5.0 Å against a 5.5 Å default. Lowering it
to 4.0 moved 79 → 66. A badly aimed probe produces exactly the reading a broken
wire does, which is the same lesson this project has now learned in four
different forms — an alternative must be checked against a known answer before
its disagreement means anything.

**A real bug fell out of it.** Probing `conservation.min_coverage` showed
`top_conserved` returning 10 residues when only 1 met the coverage requirement.
It sorted the failing ones to the bottom with `−inf` and then sliced `[:n]`,
which returns them anyway — each carrying its true conservation value with
nothing to mark it as unqualified. That is reachable from the CLI as
`conservation --top`, so the project has been able to print a confidently wrong
list. Sorting a failing entry to the bottom is not excluding it. Now fixed to
return at most *n*, all passing, and pinned.

**Nothing was deleted.** The round anticipated that some parameters might be
better removed than wired, on the grounds that a parameter no calculation needs
should not sit in the registry claiming a citation. Every one of the 21 turned
out to have a real call site, so the question did not arise.


## Round 49 — checking the path instead of the number, and what it found

**The distinction the round rests on.** `verify_claims` recomputes a documented
number and compares. That catches drift in the *value* and is silent about
everything else: which structure file was read, which registered parameters were
consumed, which commit produced the document. A number can be perfectly correct
and completely untraceable. `make provenance` walks the five links — document,
code, parameters, data, commit — and reports where each breaks.

**The parameter and data links are measured, not declared.** `record_sources`
wraps the single registry read path and the two file doors for the duration of
one claim, so what comes back is what the computation *actually touched*. That
is the whole reason the round found anything: a static reading of the source
would have agreed with the declarations.

**The checker's first run found a bug in the checker.** It reported Round 22's
δ = −0.2105 as missing from `VALIDATION_ROUND22.md`. It is written there — as
`−0.211` with U+2212 MINUS SIGN, where the pattern matched only the ASCII
hyphen, so the document's negative numbers parsed as positive. I had written in
the docstring that this check "cannot fail when it should pass"; that was wrong
and is now corrected in the docstring rather than quietly patched. The fix
deliberately does *not* treat an en-dash as a minus, because the science
documents use it for ranges like "2.7–4.7 mN/m", and a test pins both.

**Then the real finding: 26 of 101 registered parameters were read by no code.**
Not unused — *inert while advertised*. Such a parameter appears in the dialog
with a unit and a citation, an override on it is recorded, reports carry the
non-default banner, and `verify_claims` refuses to run against it, while the
number it claims to control does not move. That is strictly worse than an
unregistered literal, which is at least honestly invisible.

**Proved rather than inferred.** Setting `pore.step` from 1.0 to 0.25 left the
8YEZ bottleneck at 0.951756 Å — bit-identical. The parameter audit passes this
by construction: `MAPPED` records that a literal *corresponds to* a registered
parameter, and a correspondence is not a wire. The audit was doing exactly what
it was written to do; it was never asked whether the connection existed.

**Fixing it needed both ends.** Wiring `pore_profile` to resolve at call time
was not enough — the number still would not move, because three callers passed
`step=1.0` explicitly and six more had their own literal defaults. A parameter
is only wired when the callee resolves it *and* no caller overrides it with a
constant. That is why the fix touched six modules for five parameters.

**One inconsistency surfaced on the way.** `analysis_pore` sampled at 1.5 Å
while the registry advertised `pore.step = 1.0`. Two different numbers for one
named quantity, in a project whose stated aim is that every number has one
traceable source. It now uses the registered value.

**Scope, and why the rest are deferred honestly.** The other 21 are recorded as
Round 49b with the same proof obligation — override it, show the number moves,
show the default is unchanged — and held by a ratcheting test that fails if the
count grows. Some are probably better *deleted* than wired: a parameter no
calculation needs should not sit in the registry claiming a citation.

**A limitation pinned rather than hidden.** Most loaders here memoise, so only
the first call reads anything. Running the whole registry gave
`hydration.score_11zc` zero parameters; running it alone gave four. An empty
`data_files` therefore means "read nothing during this call", not "depends on
nothing" — documented in the module and pinned by a test, because a silently
order-dependent provenance report is worse than none.


## Round 48 — the LoF gap, and the ceiling measured rather than argued

**Why this was worth running at all after Round 47.** Round 47 had just
concluded that no reachable dataset resolves the effect the mechanical
predictor produces. That conclusion was about a *per-variant* predictor at
δ ≈ 0.25. Round 48 asks a different question — do LoF and GoF variants sit at
structurally different *positions* — where the design has 80% power at
\|δ\| ≥ 0.495. A large positional effect was still detectable, and the feature
table already existed, so the test was cheap. It was pre-registered as
exploratory below that threshold and committed alone (`7ffb008`) before
anything ran.

**The ceiling went in §2, before the hypothesis, on purpose.** Every previous
round put the "this cannot become the predictor we want" caveat at the end. For
this round it is the dominant fact: a feature computed on the wild-type
structure has zero within-position variance, so it cannot distinguish two
variants at one residue. Stating it first meant a positive result could not have
been over-read afterwards — the pre-registration says in advance that it would
license "GoF and LoF variants occur at different positions" and nothing more.

**And then it was measured rather than asserted**, which is what the roadmap
asked for. Between-position share **1.000000**, within-position **0.000000**,
via the same `variance_decomposition` Round 26 was judged by. The demonstration
is R2456: four curated variants, three GoF and one LoF, all valued 0.127326 to
every digit. Set beside 4.9% (Round 7) and 52.5% (Round 26), the progression is
the clearest statement of the confound the project has produced. The test for
this includes a control on the instrument — a synthetic case where the
decomposition *does* return non-zero — because a zero is only meaningful if the
measurement could have said otherwise.

**The result: the fifth null, and the flattest.** Primary δ = +0.036, p = 0.509,
AUROC 0.482, all three decision clauses failing. Nothing in the six-endpoint BH
family survives; the smallest q is 0.930. Distance to the gate — the endpoint
with the clearest mechanical story — separates *exactly* nothing, δ = +0.000.

**The negative control is again what decides the reading.** Distance from the
three-fold axis was pre-registered because no mechanism predicts it, and at
δ = +0.268 it out-performs every mechanistic endpoint. So the spread across
endpoints is noise at this sample size, and any single large result would have
been indistinguishable from it. Round 41 produced the same diagnostic; having it
twice is what makes it a property of the data rather than of one predictor.

**A trap the write-up had to name.** The median LoF position is *more* exposed
than the median GoF position (0.219 vs 0.143), which reads as a reversal of the
hypothesis. Cliff's δ is +0.036 — a rank statistic — so the distributions
overlap almost entirely and the median gap rests on a few positions. Reporting
the medians alone would have described an effect that is not there. The
pre-registered statistic is what prevented that, and a test pins the warning.

**Kept out of the record deliberately.** `prediction_record.VALIDATION_RECORD`
is scoped to the mechanical ΔΔG score — the number the GUI shows next to a
variant — so Rounds 41 and 48 do not belong in it. What *was* stale there is now
fixed: the caveat list said "roughly 130 variants would be needed" (Round 47
made that exactly 134 and added the reachable ceiling of 59), and it did not
mention that two other predictor families had been pre-registered and failed.
A user reading the GUI now sees five tests, five predictors, five nulls.


## Round 47 — a predictor that could survive its own data limit

**The question, and why it is not another test.** Round 26 made the mechanical
predictor substitution-aware (within-position variance 4.9% → 52.5%) and
Round 36 tested it: δ = −0.249, p = 0.405, the fourth null. The roadmap asked
what effect size is now detectable and whether the predictor reaches it — and
said explicitly not to run the comparison without a pre-registration first. So
this round is a **design analysis on already-recorded effects**, and the module
is built so that it cannot quietly become a test.

**The answer is no, and the margin is what matters.** Round 26 genuinely helped:
the effect grew from −0.083 to −0.249 and the required sample fell from over 800
variants to 134. But the reachable maximum is 59 — the 46 directional missense
variants plus Round 45's 35 harvest candidates, times the 74% that survived
Round 36's modelling gate. At 59 the minimum detectable effect is 0.356 against
an observed 0.249, and power is 0.51.

**Why that distinction is worth a module.** "We need more data" invites another
curation round. "The data that could exist is not enough" says a fifth
pre-registered test on this variant set should not be run *whatever predictor
goes into it*, because a null would be uninformative by construction. Only the
second is actionable, and it is the one the numbers support.

**Derived, not restated.** The first draft hard-coded `OBSERVED_EFFECT = -0.249`
and a `split = 0.55`. The parameter audit flagged both, which was the right
call for the wrong reason — the real problem was not that they were
unregistered but that they were *copies*. Everything now reads from
`prediction_record.VALIDATION_RECORD`: the effect, the group sizes, the split
(19/34 = 0.559), and the survival rate (34/46 = 0.739). The module cannot
disagree with the record it argues about. Fixing this moved the ceiling 61 → 59,
because my hand-computed survival rate had used the wrong denominator.

**A guard that was too crude.** The test forbidding a comparison originally
banned the string `cliffs_delta` anywhere in the module — and then failed, because
reading `record.cliffs_delta` off the recorded result is exactly what the module
should do. Rewritten to ban the *imports* that could compute a fresh statistic.
The distinction is the whole point of the round: reading a recorded effect is
allowed, producing a new one is not.

**Also corrected.** A first draft of the document claimed
`analysis.paired_variant` had been built for within-position pairing and found
too few positions. It had not — that module pairs a variant *structure* against
wild-type structures, which is Round 34's n = 1 problem. Checked before
committing rather than after.

**Guarded.** `feasibility.required_n` (134 ± 22) and `feasibility.ceiling`
(59 ± 1) are in the claims registry, so the conclusion fails loudly if a larger
variant source ever appears — which is the one thing that would overturn it.


## Round 46 — the one pair, and the control that interprets it

### The question a single pair can answer
Round 34 left exactly one comparison open: 8YFG (R2456H) is the only deposited
entry that resolves its own mutation and is coordinate-distinct. Two structures
always differ, so "does R2456H differ from wild type?" has no useful answer. The
question that does is **"does it differ by more than wild-type entries differ
among themselves?"**

### Widening the control, not narrowing it
The obvious comparison is 8YFG against 8YEZ — one pair. Instead the wild-type
side is three independent entries (8YEZ, 8ZU3, 8ZU8), which gives a measurable
within-wild-type spread.

The subtlety is which entries count as independent. 8YFC and 9VMX are
byte-identical to 8ZU3 — Round 34's finding — so including them would add two
pairs differing by exactly zero, shrink the wild-type spread, and make the
variant look more distinct. They are excluded **by coordinate fingerprint rather
than by name**, so a future duplicate is caught without anyone maintaining a
list.

### The result
R2456H: bottleneck **0.808 Å**, wetting **0.904**. Wild type: 0.673–0.930 and
0.457–0.986. The variant sits *inside* both ranges, and its largest difference
from any wild-type entry (0.135, 0.446) is smaller than the largest difference
between two wild-type entries (0.257, 0.529).

Tested generously — *any* measure exceeding the wild-type spread would have
counted as distinguishable — and it still is not.

### Why that is not a disappointment
Every deposited human structure is closed, and R2456H's phenotype is *slowed
inactivation*: it changes how long the channel stays open, not how wide the
closed pore is. A closed structure has no obligation to show it.

Saying that is part of the result rather than an excuse for it, and it was worth
measuring precisely because the intuition "a severe gain-of-function variant
should look different" is strong and, here, wrong.

### Testing the control both ways
"Not distinguishable" and "the test cannot distinguish anything" look the same
from the outside. So the control is exercised on synthetic cases: a variant far
outside the wild-type set is detected, one just inside is not, and — the one
that matters — a wild-type set of *identical* structures makes even a 0.001
difference "distinguishable". That last test is the duplicate trap made
concrete, and it shows exactly what including 8YFC and 9VMX would have done.

709 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 45 — the last route to more data, costed

### The premise, and what it was worth
"Harvesting published supplementary tables is the only route to that number that
does not require new experiments." Round 36 needed ~130 directional variants and
had 34, so this was the highest-value round remaining.

It yields **2**, and neither carries a direction.

### The funnel
86 raw substitution matches across 15 open-access papers → **66** pass the
wild-type gate → 66 mappable to human → **35** not already curated → **2** carry
an extractable measurement.

Three things in that chain are worth stating separately.

**The gate removes 23%.** cDNA changes are written in exactly the shape of
protein substitutions — C7366T looks like a cysteine-to-threonine at 7366 — and
without the wild-type check they would enter a curated resource as variants. The
regex is deliberately left loose and the gate does the work, because tightening
the pattern would drop real variants while still admitting these.

**40 of 66 are mouse-numbered**, against 18 human. That is the project's standing
trap made quantitative: a harvest that assumed one numbering would mis-assign the
majority. Conversion goes through the alignment, and a test asserts the offsets
are *not* all equal — otherwise the harvest could pass while accidentally relying
on a constant shift.

**31 of 66 are already curated.** The hand curation was more thorough than I
expected, and that number bounds what any harvest of this corpus could ever have
added. It is the most encouraging thing in the round.

### Why it fails, which is not where I expected
Not the gate. Of the 35 fresh candidates, **33 appear only in prose**. Across all
38 downloaded papers the *tables* contain four substitution strings, two of them
cDNA. The measurements this project needs are in sentences and in
non-open-access supplements — the round's premise about "supplementary tables"
does not describe the corpus that is actually reachable.

### The line I did not cross
`Candidate` has no direction field, and a test enforces its absence. It would
have been easy to regex "slowed inactivation" or "increased current" out of the
context sentence and emit a GoF/LoF label — and that would have put unreviewed
labels into the set four pre-registered blind tests depend on. The harvest
produces candidates *for* curation; it does not curate.

What it does produce is a bounded, reviewable list: 35 substitutions with their
sentence and their source. Round 57 proposes working through them by hand, which
is now a task with a known denominator rather than an open-ended search.

### Block M
Five rounds, five answers about data rather than method. The destination is now
measurably out of reach with the available data, and every route has been tried
and costed — which is itself a result. What the project should not do next is a
fifth predictor on 34 variants.

Two habits from these rounds are worth keeping: every negative result was
measured against a **control**, and each was implemented as a **check that will
change by itself** when the world does.

702 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 44 — reading the confidence we had already downloaded

### What was ignored
`fetch_alphafold(with_pae=True)` has existed for many rounds. Nothing ever read
the matrix. pLDDT — which the project does use, for colouring — says how well a
residue's *local* environment is predicted, and says nothing about whether two
domains are correctly placed relative to each other. That second question is the
only one a hybrid model cares about.

### The trap I nearly fell into
The obvious comparison is mean PAE across the seam versus within each region:
**27.3 Å against 16.1 and 20.7 Å**. That looks like a decisive statement that
AlphaFold has no idea where the distal blade sits.

It conflates two things. Pairs spanning a boundary near one end of a 2521-residue
chain are systematically *further apart in sequence* than pairs within a region,
and PAE grows with separation. Controlled for separation the picture changes:
the penalty peaks at **+4.3 Å on a 31.75 Å scale**, and at 50–150 separation it
**reverses** — cross-seam pairs score 13.25 against 15.82 within.

So **pLDDT agrees with where the seam had to be placed and PAE does not single
it out**, which is a more interesting answer than either alone.

The control is tested both ways, against a planted penalty and against a matrix
with none, because "no penalty found" and "the control does not work" are
otherwise indistinguishable. A third test builds a matrix with *no* seam
penalty and shows the uncontrolled comparison manufacturing one anyway — the
mistake made concrete rather than described.

### The finding that actually matters
PAE is **85% saturated** beyond 800 residues of separation, and — the part worth
stopping on — **80% saturated within the cryo-EM-resolved core alone**, a region
experiment places with confidence.

AlphaFold does not determine PIEZO1's long-range architecture *anywhere*. So for
`hybrid.py` the seam is not the weak point. Wherever the cut is made, the global
arrangement of the two halves is unconstrained by the prediction, which argues
for placing the distal blade by the experimental C3 symmetry and dome geometry
and treating AlphaFold as a source of local fold only.

That is a design constraint on a module that has not been written yet, arrived
at before writing it, which is the cheapest moment to learn it.

### What I did not do
AlphaFold3, Boltz-2 and Chai-1 were named in the round. I did not pursue them:
the PAE result says the limitation is the prediction's global architecture rather
than its vintage, so a newer predictor is not obviously the fix — and any
replacement would need its own confidence readout before it could be trusted more
than this one. Recorded on the roadmap rather than silently dropped.

692 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 43 — the ligands that have no structure

### What the resource is really for
`ligands.json` had been 📋 since the project began, and it would have been easy
to fill it with pocket residues from the literature. The reason not to is that
**no PIEZO structure with a bound small-molecule modulator has been deposited**.
Every "Yoda1 pocket" in the field comes from mutagenesis, docking or geometry —
and once residues are drawn on a structure they look exactly like residues that
were observed there.

So the file's organising idea is a graded `site_evidence` field that travels
with each site, and a build that rejects `bound_structure` outright.

### Verifying rather than asserting
The claim "no bound modulator exists" is the kind that quietly stops being true.
`deposited_modulators()` scans the heteroatoms of all 21 downloaded structures
against the set that is legitimately there — lipid, detergent, glycan, ion — and
reports anything else. Nothing else is present. If a Yoda1-bound entry is ever
deposited and downloaded, the build **fails**, which is the correct outcome
because the resource would then be wrong.

That is the same shape as Round 42's control: an absence is only evidence if the
instrument can detect presence.

### What the six ligands actually support
Only **one of six** carries a residue-level site — Yoda1's 1718/2075/2078, at
`docking_md`. The other five record *why* they have none, which matters because
silence would read as "not looked at" rather than "deliberately absent":

- **Dooku1** competes with Yoda1, which implies a shared site but does not
  locate one. No residues claimed.
- **GsMTx4** partitions into the outer leaflet and acts on the bilayer. Recording
  protein residues would misrepresent the mechanism, not merely overstate it.
- **Yoda2** is assumed to share the Yoda1 pocket by analogy — an assumption, so
  no site.
- **Jedi1/2** act through the blade and beam, mapped by mutation, but no residue
  set is specific enough to call a binding site.

### Gates the build enforces
Chemistry is fetched from PubChem and the returned **InChIKey is compared with
the recorded one**, so a wrong CID cannot pass silently — the same failure mode
that gave this project a bone-marrow-transplantation paper as a PIEZO1 citation
in an earlier round. Every citation must resolve in `references.json`. A site
with residues must carry a citation; a site without must carry an explanation.

### Validated against the project's own anchors
Yoda1 **EC50 26.6 µM** (Syeda 2015) and GsMTx4 **Kd 155 nM** (Bae 2011), both
matching the ground-truth table in `CLAUDE.md`. Those are the two numbers the
project committed to at the start, and they are now in a machine-readable
resource with their provenance rather than only in prose.

683 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 42 — the premise was wrong, and that is the result

### What the round assumed
"Two sources make this cheap." The idea was sound: comparing our geometrically
found lipid sites against occupancies from simulations other people have already
run would be an independent check, and an independent method agreeing is worth
more than a better version of ours.

The data is not there.

### What was measured
**MemProtMD holds 1 of 21** catalogued PIEZO entries — only 3JAC, the 2015
structure. Absent are 7WLT, 7WLU, 6B3R and 8YEZ: every structure this project
actually uses.

That absence is only evidence because the probe has a **control**. 2RH1 and
1M0L return 200 on the identical request, so "PIEZO is not there" is
distinguishable from "the URL is wrong". This project has been caught by that
class of error more than once — a PMID resolving to an unrelated paper, a frame
mismatch reporting a closed channel as conducting — and an absence without a
control is not a finding.

**And the one entry cannot answer the question.** 3JAC resolves 918 of 2,547
residues (36%), and of the 15 curated lipid-associated residues it resolves 4:
the PIP2 cluster in full, and none of the three blade basic clusters. A
simulation of a model that omits the lipid-binding residues cannot report their
contacts.

Zenodo's PIEZO1 records turned out to be microscopy TIFFs and PDFs rather than
trajectories. GPCRmd is GPCR-specific and PIEZO1 is not a GPCR.

### Where I stopped
MemProtMD's site is a single-page app and its analysis is browsable but not
fetchable: ten candidate API paths returned 404 against page URLs that return
200. I stopped guessing endpoints at that point rather than continuing, because
the coverage result already answered the round — the API being inaccessible is a
second obstacle behind a first one that is decisive on its own.

### What was built instead
`analysis/external_md.py` implements the **check**, not the comparison. The
conclusion is therefore reproducible, and it will change by itself: if MemProtMD
ever ingests a modern PIEZO structure the coverage function says so and the test
guarding this null fails.

One test exists specifically for the way this round could have fooled itself: an
offline run must return "not checked", never an empty coverage that reads as a
measured absence. A network failure would otherwise manufacture exactly the
conclusion the round reached.

### The pattern
This is the fourth round in a row where the answer was about data rather than
method — Round 34 (one informative variant structure), Round 36 (34 variants
against the 130 needed), Round 41 (an unconstrained gene), and now this. The
project's models are not the limiting factor and have not been for some time.

674 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 41 — the fourth null, and the clause that caught it

### The gene-level answer arrived before the test
gnomAD's constraint metrics for PIEZO1 came back in the first query, and they
are discouraging: **LOEUF 1.10**, **pLI 1.5e-100**, `oe_mis` **1.45**, `mis_z`
**−11.3**. The gene is not merely unconstrained — missense variation is
*enriched* relative to the mutational model. Loss-of-function tolerance is
biologically unsurprising (E756del is a common African allele), but the missense
enrichment means selection is not removing missense variation from PIEZO1 at
all.

That went into §2 of the pre-registration, before the per-position test, because
it changes what a null means. Looking for regional depletion inside a gene with
no global depletion was always a long shot, and saying so beforehand is the
difference between a predicted null and an excuse.

### The result, and the clause that decided it
Primary Cliff's δ **−0.269**, CI **[−0.595, +0.074]**, **p = 0.0477**, on 18 LoF
versus 24 GoF positions. Medians differ in the hypothesised direction — LoF
2.275 against GoF 2.520 missense per residue, about 10% lower.

**p = 0.0477 is below the conventional threshold.** The pre-registered rule
required three things: p < 0.05, δ negative, *and* the bootstrap interval
excluding zero. The first two pass; the interval contains zero. **Fail to
reject.**

Had the rule been written as "p < 0.05" — which is how most of these are written
— this would have been the project's first positive result, on a p-value 0.0023
under the line, with an effect estimate whose interval comfortably contains no
effect. The conjunction was fixed in advance for precisely this case, and this
is the first time it has been the deciding clause rather than a formality.

### The negative control did more work than the primary
§7.4 pre-registered the raw per-residue count as a control that "should be
dominated by shot noise and show nothing". It gives δ **−0.231, p 0.078** —
statistically indistinguishable from the ±25 smoothed predictor, and from both
the ±10 and ±50 alternatives.

So the windowing that makes this a *regional* constraint estimate contributes
nothing. Whatever weak tendency is there is present in the unsmoothed counts,
and the predictor cannot be told apart from its own negative control. That is a
stronger reason to disbelieve the primary than the interval alone, and it exists
only because the control was pre-registered rather than added afterwards to
explain a result.

The summed allele count ran the *opposite* direction (+0.250): how often a
position is hit says the reverse of how many ways it is hit. Uninformative at
q = 0.716, and reported because the family is reported whole.

### R2456 again
Excluded in writing beforehand as the only position carrying both directions —
R2456H/K/P gain, R2456C loss. A position-level predictor is structurally
incapable of separating them, so the exclusion is both necessary and the
cleanest available statement of this predictor's ceiling. The same residue was
Round 7's diagnostic example for the mechanical predictor. Two different
predictors, the same residue, the same reason.

### Where four nulls leave the project
Four pre-registered tests, four nulls, four different predictors: elastic-network
ΔΔG, FoldX, substitution-aware ΔΔG, and now population constraint. Every
diagnosis has ended in data rather than method.

`docs/SCIENCE.md` §8b now says four rather than three. The honest summary is that
this project has found **no signal that separates gain- from loss-of-function**,
from structure or from population genetics.

I named the biggest methodological weakness for whoever picks this up: observed
counts conflate "few variants seen" with "few expected". gnomAD's own regional
missense constraint model fits expected counts per region and would separate
them — that is a better instrument than the one used here, and it is the fair
next test rather than another window width.

668 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 40 — reproducing a published figure, and getting one of two

### Choosing what to reproduce
Haselwandter & MacKinnon's footprint would largely have tested our own linear
solver against our own parameters. Young et al. 2023 was the better choice
because their **full rate set is published and registered**, so the output can
be checked against papers the model was not built from. That is the difference
between an integration test and a restatement.

### One agreement, one disagreement
**T₅₀ = 2.711 mN/m** against Lewis et al. 2015's measured **2.7 ± 0.1** — a
0.4% difference. Young's rates, this project's solver, a third group's
measurement: three independent things landing on the same number. That is the
strongest cross-source agreement in the project.

**τ_inact = 73.3 ms** against Bae et al. 2013's **8.6 ± 0.4 ms** — **8.5×**
slower. Checked that this is not a fitting artefact: the decay is cleanly
mono-exponential, and a bi-exponential fit sends its second time constant to
2.2e6 ms, i.e. a constant. k₂ carries the timescale — at 8 s⁻¹ it sets ~125 ms
before the rest of the system pulls it to 73 — and reaching 8.6 ms needs a
12.8-fold increase to ~103 s⁻¹.

The agreement is what makes the disagreement informative. If both had failed the
natural conclusion would be that our solver is wrong.

### It justifies a policy that already existed
`kinetics.wt_tau_ms` already carried the note that mutants are calibrated by
**fold change** against the wild-type τ, "never by absolute τ across
preparations". That was written as caution, before anyone measured how far apart
the preparations are. This round supplies the number: 8.5×. The policy was right
and is now quantified rather than merely prudent.

### An API I misread, and why it matters here specifically
`with_modification` takes **fold changes, not absolute rates**. Passing
``k2=8.0`` to a model whose k2 is already 8 gives **64**. I read it as a setter,
measured τ = 13 ms, and briefly had two irreconcilable numbers for the same
model — 13 ms and 73 ms.

What makes this worth recording is that the misreading was wrong by roughly the
same factor the round was measuring. Had I not already had the 73 ms figure from
the default model, 13 ms would have looked like a modest disagreement with Bae
rather than a large one, and the round's conclusion would have flipped. The
reconciliation came from diffing the two model objects rather than from
re-reading the docstring, which is the faster route when two numbers disagree
and both look plausible.

`calibrate_k2_for_tau` came out well: it **refuses** an out-of-reach target and
states the reachable range instead of returning its search bound. Clipping would
have looked like an answer.

### Notes
- The figure is `docs/img/young2023_response.png`: tension steps, the T₅₀ curve
  with Lewis's band, and τ against tension on a log axis with Bae's band, so
  both the agreement and the 8.5× gap are visible at a glance.
- No parameters changed; this round consumed the registry rather than adding to
  it, which is what an integration test should do.

660 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 39 — putting the record where the user meets the claim

### The problem
The mechanical ΔΔG is this project's central claim and it has failed three
pre-registered tests. All of that lived in `docs/` and the CLI. A GUI user could
select a variant, see a classification and a score, and have no way of learning
that the score has never predicted anything.

That is the worst possible arrangement: the machinery is visible and the
evidence about the machinery is not.

### What was built
`analysis/prediction_record.py` holds the record as **data, Qt-free**, so the
GUI, the CLI and the tests read the same numbers and cannot drift. Following
`tour.py`, which solved the same problem for the guided tour.

Three things now reach the user:

**Evidence level, beside every variant.** A classification alone reads as a
fact. For 20 of the 46 directional variants the direction is *inferred* from
which disease the variant causes rather than measured, and the panel now says
which, colour-coded, with the sentence explaining the difference.

**The conflicting variant, flagged where it is met.** V598M is the one variant
where the curated and ClinVar-inferred directions disagree. Round 36 excluded it
in writing beforehand; the panel now shows "sources disagree: curated says GoF,
the disease mechanism implies LoF — this project reports the disagreement rather
than resolving it".

**The record itself**, under *Analysis → Variant prediction record…*: all three
tests, the power statement, and five standing caveats.

### The caveat that took the most thought
Four of the five caveats say what the score cannot do. The fifth says what it
still can: *use it to ask which residues sit in mechanically coupled positions;
do not use it to assign a direction to a variant.*

That distinction is real and worth stating. Round 7's diagnostic was that the
predictor reports **where a residue sits** rather than which substitution
occurred — which is a failure for direction prediction and is exactly what a
coupling map should do. Leaving it out would have made the interface
discouraging rather than accurate, and a caveat nobody believes is a caveat
nobody reads.

### Keeping it honest
`verify_record()` re-reads the stored Round 36 run and fails if the frozen
constants drift; measured agreement is 1.2e-4. Same discipline `analysis.claims`
applies to the documented numbers — a frozen record that drifts from its own run
is worse than no record, because it looks like provenance.

The GUI-reachability guard from Round 34 caught the new analysis the moment it
entered the registry without a menu entry. Second time it has paid for itself,
and both times within a round of being written.

### Block L appended
Five rounds since the last review. The pattern worth naming: **three times in
five rounds an alternative built to check the pipeline was itself defective**,
and each time it produced a plausible number rather than an error — a 15–16 nm
tag distance, an 89% dome model error, 32 pA through a closed channel. What
caught all three was calibrating the alternative against a known answer first.
Round 51 proposes making that a written rule rather than a habit.

651 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 38 — the error the intervals never contained

### What the round is for
Round 29 attached three kinds of spread to the headline numbers and stated on
each that **model error is not captured**. Bootstrapping a sphere fit says how
well a sphere is determined; it cannot say whether a sphere was the right shape.
This estimates that missing term where a second defensible model exists.

`ModelError` is deliberately a separate type from `Sensitivity`. Changing a
spring exponent is a knob within one model; changing a sphere into a spheroid is
a different claim about the object, and letting a "sensitivity range" stand in
for "we do not know the shape" would be the wrong kind of tidy.

### The headline: the dome interval is too narrow by six times
A sphere fitted to the transmembrane surface gives **9.454 nm** with a geometric
rmse of 6.180 Å. An oblate spheroid gives flattening **+0.431** and an apex
curvature of **14.991 nm** with rmse **5.243 Å** — it fits *better*, as it must
with an extra parameter, and the surface is plainly not spherical.

Those two radii differ by **5.54 nm (59%)**. The bootstrap interval on the
sphere is **0.92 nm**. So model error is **6.0×** sampling error, and the
interval this project has been quoting measures how well a sphere is determined
rather than whether a sphere was right.

The published 10.2 nm is itself a sphere-based number, so the sphere stays the
right comparator for the literature. What changes is the honest width.

### The fitter was wrong, and would have been reported as science
The first `fit_spheroid` alternated between the centre and the semi-axes with a
hand-rolled gradient step. On the dome it produced an apex curvature of 18.4 nm
and an "89% model error" — a number I nearly wrote down.

Testing it on a *known* spheroid first is what stopped that. Given a full
surface with true (a, c) = (100, 60) it returned (163, 98) and put the centre
89 Å away. Both axes inflated by the same 1.63×, which is the signature of a
drifting centre rather than a shape error — the ratio was right, so any check on
shape alone would have passed it.

Replaced with the exact linear solution: in the axis frame the implicit
quadric is linear in six coefficients, the null vector of the design matrix
gives them up to scale, and completing the square recovers centre and axes with
no iteration. It now recovers known spheroids to **0.01 Å**.

The lesson is the ordering. An alternative model is an instrument, and an
instrument has to be calibrated on something with a known answer before it is
pointed at the thing you do not know. This is the second time in two rounds that
the alternative route turned out to be the faulty one.

### Two smaller results
**Springs, 5.2%.** Cumulative gating overlap is 0.890 / 0.912 / 0.937 across the
three published spring models. All find the transition; the elastic network is
far less model-sensitive than the dome geometry.

**Pore, exactly zero — with a mechanism.** Apollonius and a uniform probe agree
*exactly* at 1.70 Å, because 7WLT's bottleneck lining is carbon and 1.70 Å is
carbon's radius. Moving the probe off 1.70 shifts the answer by precisely the
offset (0.300 Å at both 1.40 and 2.00), which proves the check is live rather
than silently returning one number. So the per-atom refinement buys nothing at a
carbon-lined constriction — a real null, with a reason.

### What is claimed
Every result says "**lower bound**" on its face. Two models disagreeing bounds
model error from below. Two agreeing does not bound it from above, because both
may be wrong in the same direction — and for a dome fitted only over a cap, that
is a live possibility rather than a formality.

640 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 37 — cross-checking four methods, and being wrong about why one disagreed

### What was built
Round 30 re-derived three headline *physics results* by independent routes. This
does the same for four *algorithms*, each swapping out the part most likely to
be wrong:

* **pore radius** — the pipeline runs a coarse polar grid then a shrinking
  pattern search, a local optimiser on a piecewise-smooth surface whose failure
  mode is a local maximum. The alternative maximises the same clearance by 20k
  uniform random probe centres, which has no basin structure to get stuck in.
* **SASA** — Shrake–Rupley with a fixed golden-spiral point set, against
  Monte-Carlo with independent random directions. A defect in the point set is
  invisible to a test that uses the same point set.
* **conservation** — Needleman–Wunsch anchoring, against exact k-mer seeds with
  no dynamic programming and no gap penalties at all.
* **PCA** — the SVD of the centred coordinate matrix, against power iteration on
  XᵀX, which touches no library eigensolver and never forms the 7389×7389
  covariance.

Each alternative is tested first against a case with an analytic answer — an
isolated sphere, a planted component, a slab with a known gap — because two
routes agreeing means nothing if neither is right.

### Results
PCA agrees **exactly** (0.0%, |cos| = 1.000000). SASA to **0.1%**.

The pore check agrees to 5.2%, and the *sign* is the useful part: brute force
finds **0.9783 Å** where the pattern search finds **0.9300 Å**. A brute force can
only match or beat a local optimiser, so the larger value means the pattern
search stops slightly short. That is under-convergence, not a wrong answer, and
a test now bounds the gap so a real regression would widen it.

### The one I got wrong
Conservation correlates at **0.817**, and I wrote in the module docstring —
before measuring — that the residual would concentrate near indels, since that
is where an alignment does work a seed-based anchor cannot. That was a
plausible story and it was wrong. The eight worst-disagreeing positions all have
coverage 1.00: perfectly aligned.

The real cause is a bias in *my* alternative. Anchoring by **maximum exact
matches** is a selection — given a choice of offsets it prefers the one where
residues agree — so it inflates conservation exactly where a position is
variable. The evidence is unambiguous: the k-mer profile's floor is **0.36
rather than 0**, it agrees at invariant positions (0.993 where the pipeline says
1.00), and it reads **0.653** where the pipeline says below 0.50. A test now
isolates the bias on random unrelated sequences, where it reads them as partly
conserved.

So the k-mer route is the weaker instrument, and the check confirms the pipeline
rather than indicting it — the same shape of conclusion Round 30 reached about
the parabola. The docstring was corrected rather than left with the tidier
explanation.

### Why that matters beyond this round
A cross-check whose disagreement is explained by a guess is not a cross-check.
The value came from measuring *where* the two routes part company, which took
about as long as writing the alternative did, and overturned the explanation I
had already committed to prose.

### Notes
- Three cross-check arguments map onto the pipeline's own registered parameters
  (`pore.leash`, `pore.search`, `sasa.probe_radius`); the two k-mer arguments are
  exempted with the reason that they are properties of the alternative
  instrument, not of PIEZO1.
- Six background research agents were restarted on request; all hit an API
  session limit and none completed.

627 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 36 — the third null, pre-registered first

### Why a third test at all
Two nulls already stand. Re-testing the same hypothesis with the same predictor
on the same data is precisely the drift `NEGATIVE_RESULT_PROTOCOL.md` exists to
prevent. What justified a third was that **both sides of the question had
changed, and neither change was chosen after seeing an outcome**: Round 26 took
the predictor's within-position variance from 4.9% to 52.5% (and deliberately
contained no phenotype comparison), and Round 27 grew the directional set from
26 to 46.

The pre-registration was written and committed **alone**, in `af37a82`, before
anything ran. That separation is the whole mechanism: with the document in
history, the endpoints, the decision rule, the exclusions and the caveats cannot
be adjusted to suit what came back.

### The result
Primary Cliff's δ = **−0.249**, CI **[−0.628, +0.151]**, p = **0.405**,
AUROC 0.625, on 19 GoF vs 15 LoF. The direction is as hypothesised; the
significance is not. **Fail to reject.** Nothing in the six-test secondary family
survives Benjamini–Hochberg (min q = 0.591).

Per the protocol, the predictor is **not adjusted and re-run**, and the record
is written in the prescribed order in `docs/VALIDATION_ROUND36.md`.

### What is genuinely interesting, and what it is not
The substitution-aware predictor gives δ = −0.249 where the volume-only
predictor Round 7 used gives −0.025 on the *same 34 variants*. Tenfold larger,
in the predicted direction. Round 26's improvement is real and shows up here.

And across the three tests the effect has grown monotonically: −0.083, −0.211,
−0.249. It is very tempting to read that as a signal emerging as the method
improves. It is **not evidence**, and the write-up says so in those words. Three
point estimates drifting the right way across designs with 13%, low and 50%
power is exactly what either a real medium effect *or* chance would produce.
Separating them needs about 130 variants; there are 34.

That is the honest reading, and writing it down was the point of having
pre-committed the power section.

### Two defects, both found by doing the round properly
**`design.sample_size_for` returned `max_n` for any positive effect size.**
`power_curve` defaults to `alternative="less"`, so a positive δ injected the
effect *against* the alternative, gave ~zero power at every size, and the
bisection walked to its ceiling — which the docstring described as "not reached
within the search range", i.e. as a statement about the design. Anyone asking
"how many variants would a large effect need?" got a confident wrong answer
dressed as a finding. Fixed to use the magnitude and take the sign from
`alternative`; it now reproduces the protocol's own recorded table.

That surfaced while *writing* a power section rather than reading one, which is
the argument for the protocol requiring the section at all.

**The AUROC came back `nan`.** `auroc()` takes scores and a boolean mask; I
passed it two groups, so the second was cast to all-True, there were no
negatives, and it returned nan. My misuse, not a library defect — caught only
because a nan in the output is loud. Had the arrays been ordinary integers it
would have returned a plausible number instead.

### One endpoint could not be run
FoldX ΔΔG is absent from the offline ProtVar cache for **0 of 34** variants, so
a pre-registered secondary endpoint was untestable. It is recorded as such
rather than dropped: §7 forbids removing a test from the family after the fact,
and "could not be run" is a different statement from "was not significant".

### Where this leaves the destination
Three pre-registered tests, three nulls, three different predictors. Round 34
established the structural side cannot supply more data — one informative
variant structure, all gain-of-function. The remaining route that needs no new
experiments is Block K §45: harvesting T50 and inactivation constants from
published supplementary tables, behind the same wild-type gate the curated set
uses.

Until that exists the project should state the position plainly rather than keep
testing, and `docs/SCIENCE.md` §8b now does.

### Notes
- Also fixed a numbering collision I introduced in Round 34: Block J claimed
  Rounds 36–40, which Block I already had. Block J renumbered to 46–50.
- The Round 36 record is pinned by `tests/test_validation_round36.py`, including
  a guard that the write-up still states what the null does and does not exclude.

618 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 35 — the calcium nanodomain, and a prediction that held

### The result
Unusually for this project, the round's prediction survived. The screened
Green's function, fed the two numbers Rounds 31 and 33 already produce, gives
**113.8 µM at the tag** against a **0.2 µM** sensor Kd — **99.82% occupancy**.
A BAPTA-based sensor on a C-terminal HaloTag is therefore saturated whenever its
own channel opens, and reports opening as a binary event.

The roadmap expected ~200 µM at 4–6 nm. The measured value is half that, because
Round 31 moved the tag to 3.95 nm rather than the 4–6 nm assumed when the
prediction was written. Same order, and far above the Kd either way, so the
conclusion the number existed to support is unchanged — which is worth saying
explicitly rather than quietly reporting a different figure.

λ comes out at **148 nm**, far larger than the tag distance, so the exponential
is essentially 1 where it matters: this nanodomain is set by geometry, not by
buffering. That is worth knowing because it means the two buffer parameters —
which are not separately identifiable anyway, only their product enters — barely
influence the answer.

### Making the claim falsifiable rather than merely robust
A prediction that survives every parameter in its own model is only interesting
if the range swept could have broken it, so the sweep spans 80 combinations of
tag distance, calcium fraction and buffering, and a test asserts that **at least
one combination desaturates the sensor**. Two do — and each needs 20 nm *and*
0.5% calcium *and* ≥1 mM buffer simultaneously.

The falsifiers are stated as numbers someone could go and measure: the tag would
have to sit at 373 nm, or calcium carry 4.4e-5 of the current, or free buffer
reach 0.14 M. All are two to three orders of magnitude from reality.

### Two things found on the way
**The sensor has a floor.** Resting calcium at 100 nM against a 0.2 µM Kd
already holds it 33% occupied, so its dynamic range is 33–100%, not 0–100%.
Asking for an occupancy below that floor has no answer at any distance, and the
solver now returns infinity rather than its search bound. This surfaced as two
test failures that looked like numerical noise and were not.

**A silent frame bug.** The report entry detected the C3 axis on the *unframed*
structure and applied it to the *framed* one. The axis and the coordinates were
then in different frames, so the pore was measured along a line that misses the
pore — and the **closed** 8YEZ came back carrying **32 pA** and making a 1.5 mM
nanodomain. Every number downstream stayed finite and plausible; the only tell
was that a closed structure was conducting. Caught by comparing the CLI output
against the standalone calculation, fixed, and pinned by a test asserting a
closed structure never reports its own current.

That is the second time in three rounds that mixing framed and unframed
coordinates has produced a confident wrong answer. Both times the symptom was a
plausible number rather than an error.

### What it joins up with
Round 32 found that a saturating labelling protocol puts a dye on all three
tags, so a 1:2:3 brightness mixture cannot come from a short incubation. Round 35
finds the sensor saturated whenever the channel opens, so brightness cannot
track calcium amplitude either. Together they say published puncta-brightness
heterogeneity points at **unreactive tags and open probability** — two things
that are measurable — rather than at dye kinetics or graded calcium.

### Notes
- Four references added through the title gate: `stern1992`, `naraghi1997`,
  `allbritton1992`, `tsien1980bapta`. 72 resolve.
- 6 parameters registered; the calcium share of the current is `unverified` and
  is the one the answer scales linearly with, so it is swept rather than trusted.
- The GUI-reachability test written last round earned its keep immediately: it
  failed the moment `nanodomain` entered the registry without a menu entry.

607 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 34 — the comparison that was not available, and a GUI audit

### What the round wanted, and why it could not happen
Compare permeation across the four deposited variant structures and read a
direction of change against the measured phenotype. Three things stopped it, and
each is measured rather than asserted:

**Every deposited human PIEZO1 structure is closed.** Bottleneck radii 0.67–0.93
Å against a 1.38 Å cation, so every conductance is exactly zero. A *difference*
in conductance between them does not exist to be measured.

**Three of the four variant entries do not contain their variant.** A1988 is
unmodelled in both entries named for A1988V; E756 is unmodelled in the E756del
entry. Only 8YFG (R2456H) shows its mutation — histidine there, arginine in
every other entry, which is the control that makes it meaningful.

**Three of the entries share one model.** 8ZU3, 8YFC and 9VMX have byte-identical
protein coordinates: 31,839 atoms, the same hash, 0.000 Å RMSD.

That last one needed ruling out before it could be reported. Far the likeliest
explanation was our own fetch cache writing one file under three names, which
would have made it a fact about this project rather than about the PDB. The
files have different sizes, different md5s, and each identifies itself with its
own `data_` block and its own title, so the identity is in the depositions. A
test pins that check alongside the finding, because the finding is worthless
without it.

### What follows
Four deposited variant entries → one resolves its own mutation → one informative,
against 68 curated variants. And all four are gain-of-function: there is no
deposited loss-of-function structure at all, so this route cannot discriminate
direction even in principle. Round 22 found too few phenotyped variants; this
finds too few structures. Both ends are data-limited.

The null is pinned by tests written to **fail if the situation improves** — if a
human structure ever conducts, or a second variant becomes informative, the
suite says so rather than letting a stale limitation stand.

### Out of band: is everything reachable from the GUI?
Asked to check, and the answer was no. `permeation` and `interactions` were in
the shared `ANALYSES` registry and wired into the CLI but absent from every
menu — computable, and invisible to anyone not using a terminal. The pattern is
easy to fall into: an analysis with no picture to draw has no obvious home in a
3-D viewer.

`ui/result_dialog.py` renders any registry result as formatted text, and
`ui/tabular_analyses.py` adds the menu entries. Deliberately generic over the
result dict: a bespoke panel per analysis would look better and would be one
more thing to drift from the function it displays. Each entry carries the caveat
its numbers need *above* them rather than in a docstring.

The guide had stopped at Round 30 and therefore misdescribed the application —
nothing on the tags, the labelling, the ion current, the canonical framing or
multiple structures. Two topics added, including what those models cannot do.

Three tests keep it from decaying: every registry analysis must have a GUI entry
point, every menu action must carry a tooltip, and the guide must mention both
what was added and its limits. The last one caught the hand-wrapped-HTML trap
this project has hit before — a phrase split across a line break fails a naive
substring test for reasons unrelated to content, so the assertion normalises
whitespace.

### Notes
- The audit exemption list grew by two: a zero-initialised counter and a pore
  step already mapped elsewhere. Both mechanical.
- No new parameters: this round measured structures rather than introducing
  physics.

592 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 33 — ion permeation, and where continuum theory runs out

### What was built
`physics/permeation.py` turns the measured pore radius profile into a current:
steady-state drift-diffusion per species, Scharfetter–Gummel discretised so it
stays stable where drift dominates, with Hall access resistance in series at
both mouths — for a pore this short the mouths are a real share of the total
resistance, not a correction.

### The Poisson half does not converge, and that is the finding
Feeding the solved potential back through Poisson diverged immediately: −0.37 V,
then −171 V, then −2×10¹⁶ V. The obvious diagnosis is a missing restoring term,
so I added the one a proper Newton step needs, ∂ρ/∂φ = −F²/(RT)·Σz²c, which is
negative and makes the operator negative-definite. It still did not converge —
the update plateaued and the potential swung ±1.5 V at every damping I tried.

The reason turned out to be physics rather than arithmetic. In 150 mM the Debye
length is **5.7–8.1 Å**; PIEZO1's open bottleneck radius is **3.3 Å**. The
double layers from opposite walls overlap completely, so the pore has no
electroneutral core for a Gummel map to relax onto. That is not a solver defect
to be worked around; it is the statement that a continuum treatment of this pore
is at the edge of its validity.

So the potential is solved in the **electroneutral limit** — current continuity
with the local conductivity — which converges, and agrees with the independent
closed-form series-resistance formula to **1.5%** (41.0 vs 40.4 pS). The Poisson
machinery is kept and `debye_length()` is reported on every result, because the
honest record is that it was attempted and why it was set aside.

### Validation
**41.0 pS for the open 11ZC against a published 25–30 pS** — high by about half.
Closed structures give exactly zero.

The roadmap asked specifically *which* mechanism blocks which, and getting that
right needed a change: my first `_blocked_reason` returned on the first match,
which made 8YEZ and 7WLU look identical when the whole point is that they are
not. `blocking_mechanisms()` now returns every reason. 8YEZ is shut by **two**
(0.95 Å bottleneck *and* a hydrophobic gate at wetting 0.82); 7WLU by **one**
(0.98 Å, wetting 0.11 — no gate).

### Why the agreement is not a prediction
Sweeping the two unmeasured confinement parameters over plausible ranges —
in-pore diffusivity 0.25–1.0 of bulk, ion radius 1.0–2.0 Å — moves the answer
across **16–94 pS**, a 5.8× span straddling the measurement. The model can be
made to agree by choosing values nobody has measured. Both are registered
`unverified`, and a test pins the span so the point cannot be quietly forgotten.

### Two numerical faults, both of which produced a number before an error
The interior rows of the discretised system are built from `D·A/h` with `A` a
few square angstroms, so they sit around 1e-18 while a Dirichlet row is 1.
LAPACK called that singular. It is not singular, it is badly scaled, and since
the interior rows are homogeneous, dividing each by its own largest coefficient
is free. The Bernoulli function also overflowed at large argument, propagating
inf into the matrix as nan; its asymptotes are now explicit. Both live in
`_pnp_kernels.py` — pure linear algebra with no ions in it, testable alone.

### Out of band: seeing the HaloTag fusion
Asked mid-round how to view the tags, the honest answer was that you could not:
Rounds 31 and 32 computed the fusion and the labelling but drew nothing, and the
numbers were reachable only through the CLI. `ui/fusion_controller.py` now draws
them under **View → HaloTag fusion** — tag bodies, linker seams in a colour the
channel never uses, the accessible-volume cloud, and the bound dyes. Everything
is styled so it reads as a model: sphere-of-gyration bodies rather than the
fold, straight seams rather than a conformation, and the envelope shown
precisely so a single sphere is not mistaken for a determined position.

While doing it I found I had **misreported** where the Round 31/32 options live.
I described them as *View*; they had in fact been added to *Options*. Moved to
*View*, which is both where they belong and what the documentation says.

### Notes
- Three references added through the title gate: `coste2010piezo`,
  `gnanasambandam2015`, `hall1975access`. The gate earned its keep again —
  the Coste entry was first rejected because my `expect` word did not appear in
  the resolved title, and a PMID lookup I tried resolved to an unrelated 2025
  paper.
- 12 parameters registered. Several tabulated diffusivities and Shannon radii
  had to move from `physical` to `convention`: the audit is right that a
  `physical` kind promises a citation, and a community-standard tabulated value
  does not have one in this bibliography.
- `report.py` and `permeation.py` both passed 500 lines and were split at real
  seams — the newer analyses into `report_tags.py`, the linear algebra into
  `_pnp_kernels.py`.
- **Not done:** the particle animation driven by the computed current. The
  physics is in place and the current is available to drive it.

580 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 32 — HaloTag labelling: an import, and what it then says

### The criterion was "exactly", so that is what was tested
The kinetics were not derived here. `analysis/labelling.py` brings across the
three equations from `halotag_binding_sim` — exposure, per-site `p(t)`, and
Binomial(3, p) over the trimer — and `compare_with_source()` re-runs the
original functions to check. Over 241 time points the maximum absolute
difference is **0.0** for `p_site`, `p³`, the occupancy distribution and the
sampled histogram, and the Monte-Carlo populations are identical channel for
channel.

That last part drove a design decision. The sampler reproduces the source's two
`rng.random((n_channels, 3))` draws **in the same order**, reactivity first.
A different order gives a statistically identical population and a numerically
different one — which would have satisfied a "close enough" test while hiding a
real divergence behind sampling noise. Reusing a single uniform draw across
snapshots is also what makes an individual channel's dye count monotonic; the
covalent bond does not reverse, and resampling per snapshot would let a
three-dye channel drop back to two while still giving the right marginals.

The equations are vendored rather than imported at runtime, because the sibling
project is not on a fresh clone's path and A5 says a fresh clone must reproduce
the working state. `compare_with_source()` returns `{"available": False}` there
and the test skips; where the source *is* present it must agree exactly.

### What the model then says, which was not what the round expected
At the standard protocol — 200 nM JF646, 30 min, live cell — labelling is
complete in **54 s** to 99%. Per-site p = 1.0000 and the population is **100%
three-dye**. There is no kinetic dye mixture at any realistic concentration;
producing one needs sub-nanomolar ligand or an incubation under a minute.

This bears on why puncta brightness is heterogeneous, where sub-saturation
labelling is a named candidate. The model separates two things that share that
name. The kinetic route is closed at a saturating concentration. The other route
is not: a population of chemically unreactive tags leaves a mixture at **every**
time, because the ceiling is `active_fraction³`. At 90% reactive that is 72.9%
three-dye and 24.3% two-dye, and no incubation removes it. So an observed 1:2:3
mixture under a saturating protocol argues for unreactive tags, not for a short
incubation.

That conclusion rests on two registered-`unverified` numbers — `k_perm_live`
and `active_fraction` — and both say so in their `source_note`. The first is a
transport estimate, the second an assumption; neither is measured, in this
project or the source.

### What the structure adds
`label_sites()` puts the statistics on the three tag centres from Round 31's
fusion model, so occupancy is a set of places rather than a count. The sites are
treated as **equivalent**, and that is a claim worth stating: the trimer is
C3-symmetric, the three C-termini sit at the same height and radius, and the
ligand reaches them from one cytosolic pool. Geometry therefore decides where a
dye is drawn, not whether it binds.

### Notes
- Three references added through `build_references.py`'s title-verification
  gate rather than by remembered PMID: `los2008halotag`, `grimm2015jf`,
  `bertaccini2025piezo1`. 65 references resolve.
- 8 parameters registered. Two of them are `unverified`, and two more had to be
  reclassified from `empirical` to `method` — the audit correctly refused to let
  an unmeasured model estimate sit in a category that promises a citation.
- `parameter_table.py` passed 500 lines; the tag parameters moved to
  `parameter_table_tags.py`, a real seam rather than an arbitrary cut.
- **Not done:** the brightness *animation*. The histogram, the per-site
  occupancy and the figure are there; nothing is rendered on the trimer over
  time. Recorded on the roadmap rather than quietly dropped.

566 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 31 — HaloTag fusion geometry, and a ghost structure

### What was built
`structure/fusion.py` places a HaloTag at each of PIEZO1's three C-termini.
There is no structure of the fusion, so the module produces an **accessible
volume** — the region the tag centre can occupy without clashing — rather than a
pose, in the manner used for FRET position modelling. Measured inputs from 6U32
(1.8 Å, TMR ligand bound): Rg **17.6 Å**, N-terminus **19.9 Å** from the centre,
ligand 21.8 Å from that N-terminus. A C-terminal fusion attaches to the tag's
*N*-terminus, so it is that offset, not the radius of gyration, that decides
where the body sits.

### The validation: two of three
**C3 symmetry — pass.** 0.0000 Å on all 20 entries. Exact by construction, since
one envelope is solved and rotated, but measured anyway: a placement that
quietly broke the symmetry would look right in a picture and be wrong in every
calculation after it.

**No steric clash — pass on 18 of 20.** 21.5 Å clearance against the tag's
17.6 Å radius on 8YEZ. 3JAC and 11ZC are marginal at 17.6 and 15.7 Å.

**Tag centre 4–6 nm from the pore exit — misses.** 3.95 nm on 8YEZ, and
3.27–4.21 nm across all 20 structures. The interesting part is that this is not
the unverified linker's fault: sweeping it from 1 to 30 residues — a 30× change
in accessible volume — moves the answer only between 3.0 and 4.0 nm, and moves
it *down*, because a longer tether wraps further around the channel and pulls
the centroid back. So the miss is structural, and it is explained: the 4–6 nm
estimate added the tag's ~2 nm anchor-to-centre offset to the anchor's 2.6 nm
from the pore exit, which assumes the tag points straight away from the channel.
Averaged over accessible directions, many of which run sideways along the
membrane, the mean is pulled in. The band is not unreachable — the envelope
spans 1.7–7.9 nm and 51% of it lies inside 4–6 nm — so the honest statement is
that **the window describes an achievable position, not the ensemble mean**.
Round 35 should take its nanodomain distance from the envelope, not the centroid.

### Two sign faults, both of which gave confident wrong answers
The sweep across all 20 structures is what exposed them; a single structure
would have passed silently.

**The pore exit is not the lowest atom.** In a curved trimer that is a distal
blade tip. Restricting to atoms within a registered `fusion.pore_mouth_radius`
of the axis finds the CTD bundle, which is what an ion actually leaves through.

**`SymmetryAxis.direction` has no fixed sign.** It comes from the rotation
operator relating two protomers, and re-detecting the axis on an already-framed
structure returns −z as readily as +z — it does for 7WLT and 8YFG but not for
8YEZ. Trusting it put those two structures' pore exit at the *extracellular*
end. Together the two faults reported the tag 15–16 nm from the pore exit for
7WLT and 8YFG against 3.9 nm for the same construct on 8YEZ. The cytosolic
direction is now taken from the C-terminal anchors, which are intracellular by
topology.

### The ghost structure
Reported mid-round: loading a new structure left the previous one drawn.
`OverlayController.clear()` finishes by rebuilding the primary view to undo
deviation colouring, and `load_structure` called it *after* clearing the old
view — so the rebuild put the old batches straight back. The bug predates this
round; it was survivable only while deposited frames sat 100 Å apart, and became
obvious as soon as canonical framing made structures superimpose. Fixed by
clearing the overlay first and nulling the view before anything can rebuild it.

**And the feature it suggested.** `ui/companions.py` adds deliberate
multi-structure display: off by default, opt-in under *View → Show multiple
structures at once*, each extra structure drawn in its own colour in the shared
frame, with the Structure panel naming what is on screen. Loading with the
option on demotes the previous primary to a companion, so the feature needs no
second way of opening a file. Turning the option off drops the extras, because a
setting reading "one structure" while three are drawn is worse than either
state. It is deliberately **not** `overlay_controller`, which superposes one
nominated structure and reports an RMSD — that measures; this displays. Every
analysis still runs on the primary, whatever else is drawn.

### Notes
- `fusion.linker_residues` is registered with the `unverified` sentinel. No
  source for the construct states a linker, so it is the one assumed input, and
  the module says to vary it rather than trust it. The sweep above is why the
  result survives that.
- `build_parameters.py` passed 500 lines, so the table moved to
  `scripts/parameter_table.py` and the gate stayed behind. 75 parameters now.
- The GL-dependent companion tests skip on the offscreen Qt platform the rest of
  the UI suite forces, so the ordering fault is *also* pinned by a test that
  needs no GL at all.

550 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Structure framing — two faults behind one complaint

### What was reported
Reloading 8YEZ after visiting another entry "looked different", and different
structures loaded in orientations that did not overlap.

### Fault 1 — the camera drifted on every load
`_reset_camera` called `Camera.orbit(0.0, -0.42)`. `orbit` **composes** with the
current rotation, which is right for a mouse drag and wrong for restoring a
standard view. Every structure load therefore added another 0.42 rad of pitch,
so the same entry came back at a different angle on each visit. Added
`Camera.set_orientation`, which sets the rotation absolutely, and used it for the
reset. `test_orbit_accumulates_but_set_orientation_does_not` pins both halves —
including a guard that fails if `orbit` ever stops accumulating, so the test
cannot silently stop testing anything.

### Fault 2 — deposited frames are unrelated
Nothing about a PDB frame is canonical. Measured across the 20 downloaded
entries, pairs sit **29–147 Å** apart before any alignment. New
`piezo1/structure/frame.py` offers three modes, exposed under *View → Structure
alignment*: `deposited` (unchanged), `canonical` (the structure's own C3 axis on
+z, cytosolic side at −z, centred) and `reference` (least-squares onto the first
structure loaded).

**Result.** Canonical framing takes those pairs to 0.9–25 Å, in every case
within about an angstrom of what an unconstrained least-squares superposition
achieves — 19.77 vs 19.70 Å for 7WLU on 7WLT, 0.91 vs 0.91 Å for 8ZU8 on 8YEZ.
The residual is real conformational difference, not framing. The one large
value, 6KG7 at 57.8 Å, is PIEZO2 rather than PIEZO1 and its free superposition
is no better at 55.4 Å.

### Three things that had to be got right, each found by measuring

**The z sign straddled the cap.** The rule "the top 10% of residues by number is
the cytosolic end" is wrong: PIEZO's extracellular cap runs to ~2457 of 2547, so
the top tenth mixes cap and CTD and its mean z depends on which is better
resolved. 7WLU and 11ZC loaded **upside down** while reporting a perfect C3 fit
— the failure was silent, which is the dangerous kind. Checked 0.02 / 0.05 /
0.10 / 0.20 against the last 15 residues over all 20 entries: the first two
agree everywhere, 0.10 fails on two, 0.20 on nineteen. `CTERM_FRACTION = 0.02`,
and the rule now *measures where the C-terminus lands* rather than predicting it.

**Chain labels can run either way round the ring.** 8YFG and 8ZU3 both present
chains A, B, D, but numbered in opposite rotational senses. Taking the labels at
face value costs 60 Å of apparent RMSD against 8YEZ — a bookkeeping error that
reads as a conformational change. All six protomer permutations collapse into
two classes of three (cyclic 72.88 Å, reversed 12.45 Å, identical within each),
so `PERMUTATIONS` holds one representative of each and the search is complete.

**Residue number alone is not a correspondence.** A residue number occurs once
*per protomer*, so `dict(zip(res_seq, xyz))` silently keeps whichever chain was
read last and discards two thirds of the atoms — and if the two structures order
their chains differently, the survivors are not equivalent. Correspondence is now
built on protomer blocks over the shared residue basis.

### Why the roll is solved for, not chosen from three
C3 symmetry leaves the roll about z free, so it can only be pinned against a
reference. The starting rule (protomer 0's centroid on +x) is coverage
dependent: 7WLU resolves from residue 576 and 7WLT from 784, and the extra blade
density swings the azimuth by 47.7° — not a multiple of 120°, so no symmetry-
equivalent choice undoes it. Given a reference, the roll is now the closed-form
one-parameter Procrustes optimum. This is a constrained fit and is documented as
one: the axis and the z sign still come from the molecule, only the genuinely
free degree of freedom is set by the reference.

### Consequences elsewhere
- The load path now hands **transformed** coordinates to the physics, so
  `test_measured_geometry_is_unchanged_by_reframing` pins that curvature, depth
  and areas are invariant (they agree to 2e-8 relative). `pore_profile` shifts by
  ≤0.22 Å because it samples on a fixed step whose origin moves with the frame —
  discretisation, not error.
- `ui/model_utils.py` moved to `structure/protomers.py`. `structure/frame.py`
  needed it, and importing from `ui` would have pointed the dependency arrow
  backwards; `analysis/claims.py` was **already** reaching into `ui` for it,
  which quietly made a headless analysis depend on PyQt. `ui/model_utils.py`
  remains as a re-export shim.
- `_open_file` referenced an undefined `rec` on four lines — a guaranteed
  `NameError` for anyone opening a file outside the catalogue. It also skipped
  the per-structure resets (modes, measurements, overlay) that `load_structure`
  does, so stale state carried over. Fixed both.
- `main_window.py` reached 574 lines; the framing, file-open and camera-reset
  code moved to `ui/alignment.py` (`AlignmentMixin`).
- `CTERM_FRACTION` and `MIN_CA_PER_PROTOMER` added to the audit's `EXEMPT` with
  reasons: both are criteria for reading a deposited file, not measured
  quantities with a literature source.

534 tests pass; `parameter_audit` clean; `screenshot_app.py --structure 8YEZ`
completes and now reports the frame in the status line.

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

## Round 29 — intervals, and three different things people call one (2026-08-06)

Rounds 18 to 28 kept finding the same failure: a number recorded with more
confidence than it had earned. The footprint area wrong by 3.5×. A null result
that could only ever have excluded a large effect. A T₅₀ eighteen times off. A
point estimate invites that, so this round attaches a spread.

### The distinction that shapes the module

Three things get called "the uncertainty" and they mean different things:

- **Bootstrap** — resample the data. This is a genuine confidence interval for
  sampling variability and the only one that deserves the name.
- **Sensitivity** — vary a method choice, such as the elastic-network cutoff.
  There is no sampling distribution here; a cutoff is not a random variable.
  Quoting its spread as a confidence interval would be a second kind of
  overconfidence dressed as a cure for the first.
- **Parameter range** — vary a registered input over its published values.
  Propagated uncertainty from someone else's measurement, not a statement about
  this dataset.

They are separate classes with different names, and the printed summary of a
sensitivity range says *"sensitivity, not a confidence interval"* in so many
words.

### Two headline numbers change meaning

**The dome radius.** 9.73 nm with a 95% CI of **[8.83, 10.34]**, bootstrapped
over 66 transmembrane surface points. The published 10.2 nm is *inside* that
interval.

The project has been reporting this as agreement-with-a-small-gap — "our 9.7
against their 10.2". The honest statement is that the two are **statistically
indistinguishable**: a stronger claim of consistency, and a weaker claim of
precision than the two-decimal figure implied.

**The gating overlap.** 0.705 at the default cutoff, but **0.554–0.723** across
cutoffs from 10 to 20 Å, non-monotonically. The qualitative conclusion — a
substantial overlap carried entirely by A-symmetric modes — holds at every
cutoff, and the E-mode contribution stays negligible throughout. The third digit
was never meaningful, and it has been quoted as though it were since Round 4.

Ensemble PC1 is 0.900 [0.796, 0.972] over ten structures; dropping 11ZC alone
moves it to 0.832, which is the largest single influence but not a dominating
one.

### Two details worth recording

The confidence level is **derived from `stats.alpha`** rather than written as
0.95. Two copies of a significance level drift, and this project already
registers one — setting α to 0.01 now widens every interval automatically.

The parameter sweep restores the registry in a `finally`, including when the
statistic raises. Leaving it modified would make every later number in the
session incomparable with the documentation, which is exactly what the override
tracking from the parameter round exists to prevent.

### What none of it does

**It does not capture model error.** Bootstrapping a sphere fit tells you how
well the sphere is determined; it says nothing about whether a sphere was the
right shape for the dome, or whether springs are the right physics for a
2500-residue trimer. That is stated on every result rather than left for the
reader to remember, because it is precisely the assumption an interval invites
you to stop questioning.

The parameter audit caught the new module on its first run — seven unregistered
literals — which is the rule from the parameter round working as intended.

Suite 489 → 507 passing.

## Round 30 — checking without reusing the derivation (2026-08-06)

Round 18 is the reason this round exists. The elastica solver integrated
Euler–Lagrange equations derived by hand, and a boundary-value solver converges
happily onto wrong equations. What caught the error was evaluating the exact
functional in a different gauge — a check that reused none of the derivation
being checked. A test written from the same understanding as the code shares its
blind spots; if the author misread the physics, the test encodes the misreading.

Three headline results, each re-derived by a route sharing no machinery with the
original. Two of the three taught me something.

### The dome: a disagreement, and the checking route was the wrong one

A parabola through the radial height profile gave **8.12 nm against the sphere
fit's 9.72** — 16.6% apart. That looks like a finding.

It was, but not the one it appeared to be. On synthetic caps of *known* radius,
the parabola is accurate to 0.6% at an 8.6° contact angle and **25.8% low at
63.4°** — and PIEZO1's dome sits at 63.3°. The parabola `h = h₀ − r²/2R` is a
shallow-cap approximation, which is the *same small-slope failure* Round 18
found in the membrane theory, now turning up in a geometry method. The sphere
fit is exact at every slope on the same synthetics.

So I built a route that is actually valid: the exact cap relation
`R = −(h² + r²)/2h`, inverted per point, no fitting and no expansion. It gives
**10.17 nm against 9.72** — 4.5% apart, both inside Round 29's bootstrap
interval and both consistent with the published 10.2.

Both alternative routes are kept. The invalid one is worth keeping precisely
because knowing *why* a check disagrees is the useful part, and a future reader
who reaches for a parabolic fit should find the answer already written down.

### The overlap: superposition-free, and it holds

Distances between sites are invariant to rotation and translation. Comparing the
observed transition and each normal mode in *pairwise distance changes* uses no
Kabsch fit and no protomer matching — so neither can manufacture the result.
**0.641 against 0.705, 9.0% apart.** Given this project has been bitten twice by
protomer correspondence, that is a check worth having.

### T₅₀: a disagreement that was mine

The first alternative used the analytic steady state and returned zero. I nearly
wrote that up as a discrepancy.

The diagnosis: at equilibrium this channel sits **~96% inactivated at every
tension**, so steady-state open occupancy runs only 0.030 to 0.036 and never has
a half-maximum. The route was computing a *different quantity*. T₅₀ is
necessarily a property of the peak transient — which is also what a patch-clamp
measures, so the pipeline's definition is the experimentally right one.

Replaced with adaptive Runge–Kutta integration of the same master equation: same
quantity, different numerics, no matrix exponential. **2.727 against 2.711 —
0.6%.** The steady-state function is kept, with a test asserting it must *not*
reproduce the peak-based number: if it ever does, one of the two has stopped
computing what it claims.

### A small thing the audit prompted

The parameter audit flagged the new module twice. The second flag was
`OPEN_STATE = 1`, an array index. Rather than exempt it, it now derives from
`STATE_NAMES.index("O")` — a hardcoded index would silently read the wrong
occupancy if the state order ever changed, and that is exactly the class of
error this round is about.

Suite 507 → 518 passing. Block I appended after the thirty-round review,
numbered from 36 since Block H already claims 31–35.

---

## Contacts on the model, and the bonds that had never been drawn

Two requests in one: draw the interaction inventory, and fix ball-and-stick,
which was showing balls only.

### The second one was not a missing feature

`Style.BALL_AND_STICK` built its cylinder batch correctly, uploaded it, and drew
nothing. Two independent faults, either of which alone was fatal:

**A sign error in `cylinder.frag`.** The ray–cylinder quadratic used `-oc` where
it needed `oc` when forming the perpendicular component. That negates B, which
negates both roots, so the near hit came out *behind the eye* and the `t < 0`
guard discarded every fragment. **No cylinder has ever been drawn by this
renderer.** Ball-and-stick has drawn balls only since the day it was written,
and the HaloTag linker seam — a feature with its own controller, its own tests
and its own paragraph in `INTERFACE.md` — has drawn nothing at all.

**Face culling.** `Scene.render` enabled `CULL_FACE` once at setup. That is
right for meshes and wrong for impostors: an impostor quad is oriented by the
geometry it stands for, not by the camera, so roughly half of them are
back-facing and were discarded. `render()` now reads a per-batch `cull` flag.

### Why the test suite did not catch it

Because nothing rendered. Every existing renderer test checked that a batch was
built, that it had the right vertex count, that the style switch changed state.
All of that passed on a shader that emitted no fragments. `test_render_impostors`
now renders to a framebuffer and **counts lit pixels**, which is the only check
that could have failed. Both faults are pinned separately, because fixing one
without the other still gives an empty picture and would have looked like the
fix not working.

Three new styles came out of it — Balls, Sticks and the corrected Ball and
stick. The help page had listed sticks as available throughout.

### The contacts: getting the default wrong first

`interaction_controller` draws each contact as a cylinder between the two atoms
`detect_interactions` found it between. It calls the analysis rather than
reimplementing it, so if the picture and the table ever disagree the picture is
wrong.

Two mistakes, both caught by looking at the screenshot rather than the tests:

**Keyed on the wrong strings.** `KIND_COLORS` used `hbond`; the analysis emits
`hydrogen_bond`. So 7,984 of 9,863 contacts — five sixths of everything found —
silently did not draw, while the status line cheerfully reported the rest. The
test now takes the kinds from `detect_interactions`'s own output rather than a
list copied into the test, which is the only version that catches a rename.

**Then it drew them, and that was worse.** With the key fixed the picture was
eight thousand green lines and unreadable. The instinct is to call that "too
many", but the real reason is better than that: most of those hydrogen bonds are
backbone *i*→*i*+4 — they **are** the secondary structure, and the cartoon was
already drawing it. So the default is now the specific contacts only:
disulfides, salt bridges, cation–π, π-stacks. 270 cylinders on 8YEZ, each saying
something the ribbon does not — and the blue salt bridges land exactly where
they should, at the pore and along the beam. What is hidden is counted on the
status line rather than dropped.

The status line carries two caveats that cannot be omitted, the same two the
table carries: no deposited PIEZO entry has hydrogens, so every criterion is
heavy-atom geometry and a drawn hydrogen bond is an inference; and a contact
belongs to *this* structure in *this* state.

---

## Five more overlays, and a checker that could not say "no"

The list this closes is the visualisation shortlist from the end of the last
session: interactions (done then), pore surface, pockets, allosteric path,
calcium nanodomain, predicted fluctuation. All five remaining are in, each as a
View-menu toggle or a Physics-panel button.

They share a rule, which is why they were worth doing together: **each reads the
result the corresponding analysis already produced rather than computing its
own.** The pore surface draws `analysis.pore`; the pockets draw
`analysis.pockets`; the fluctuation colouring reads `ModeSet.msf`. Two rankings
of the same pockets, or two pore profiles from different parameter values, would
both look entirely reasonable on screen and there would be nothing to say which
was right. So the picture and the panel are the same object, asserted by
identity rather than by closeness.

### The calibration found a real defect, in the checker rather than the thing

The allosteric route is the most persuasive picture this application can draw
and the easiest to over-read: a single line from blade to gate reads as *the*
pathway. Round 25 already knew it is not — the beam is a near-degenerate
alternative — but that lives in a test, and nobody looking at the picture reads
tests. So the drawn route carries a measurement: re-run the same search with
this route's own steps removed, and report what the best remaining route costs.

The rule says calibrate a checker before believing it, on a case whose answer is
known. Two graphs: one where a single bridge is the only way across, one lattice
where the routes are interchangeable. The lattice passed. **The bridge did not.**

The first version suppressed the route's correlations to the registry's
`min_correlation` floor. That does not remove an edge — it leaves one costing
−log(0.001) = 6.9, which is finite, so the search walked straight back over the
route it had been told to avoid. On the real trimer it returned 1.0013× and
looked perfect, because on a real network an alternative always exists. No input
could ever have made it answer "unique". A checker that cannot say no is not a
check, and this one would have passed forever on the only data it was ever
pointed at.

Fixed at the source: `allosteric_path` gained an `exclude` set of site pairs
that are deleted from the graph before dijkstra runs. The real answer did not
move — 1.0013× either way — which is the point. The correction was invisible in
the measurement and total in what the measurement was entitled to claim.

### The endpoints were trivially wrong, and the picture is what showed it

The first version handed every blade residue to the search. The result was a
five-step hop from THU9 — the blade unit that happens to sit nearest the pore —
straight to the gate, never going near the beam. Perfectly correct as "the
cheapest blade-to-gate route", and useless as a picture of a lever.

The source is now the most distal blade unit the entry actually *resolves*,
taken from the domain records' own `thu_index`: THU4 on a deposited structure,
THU1 on a full-length model, and reported rather than assumed since it changes
what the route means. On 8YEZ that gives 32 steps through
THU4 → THU5 → THU6 → THU7 → THU8 → THU9 → Anchor → CTD.

### Two surfaces that are correct and unwatchable

The calcium nanodomain's iso-surfaces are the ones carrying Round 32's
conclusion: at 11ZC's 2.43 pA the sensor is still 90% occupied at **119 nm** and
half-occupied at its Kd at **372 nm**, against a channel reaching 14 nm. Drawn,
that is a viewport of shell with a speck inside, and the speck is the protein.

This is the far-field footprint problem again and it takes the same answer: the
numbers go on the status line and the surfaces are not drawn. What *is* drawn is
the near field at the protein's own scale — shells at decade concentrations, out
to a stated multiple of the structure's own extent. The budget is the part the
tests check, in both directions: the occupancy surfaces must measurably exceed it
(if they ever fit, draw them) and something must still fall inside it, because a
filter that admits nothing looks exactly like a bug.

And a shut structure draws **nothing**. `report_tags` borrows 11ZC's current when
the loaded entry is closed, labelled — but a caption cannot carry that weight
against a picture, because a cloud drawn around 8YEZ reads as 8YEZ's. The empty
screen and its sentence are Round 34's result: no deposited human PIEZO1 entry
conducts.

### The pixel count, again

Every controller test here checks that a batch was built with the right
contents. Every one of those passed for the entire life of the renderer while
cylinders drew nothing. So `test_ui_overlays_render` renders each overlay to a
framebuffer and counts lit pixels.

It earned its place immediately. The route's tube is thinner than its node
markers, so a screenshot showing beads is no evidence the tube drew at all — and
the tube is what carries the correlation colouring. Hidden the nodes; the tube
still draws. That could not have been established from the picture, because the
route on 11ZC runs mostly in the membrane plane and an edge-on camera
foreshortens 122 Å into a blob.

### Smaller things

- `build_sphere` is a new mesh primitive. The impostor `SphereBatch` draws
  spheres far better and has no alpha channel, which is fine for an atom and
  fatal for an iso-surface that encloses the protein.
- The pore bands are applied **ascending**, so each slice takes the highest band
  it clears, with the Round 76 pLDDT defect pinned: applied in declaration order
  a wide-open pore comes out entirely red, which is a perfectly plausible
  picture of a shut channel.
- The bottleneck is marked by drawing its lining residues, not by drawing its
  own probe sphere differently — that would have been a lie about a radius.
- The two value colourings share one slot, so each unticks the other with
  signals blocked. Without the block, unchecking fires the other handler and
  resets the colouring that was just set: the button works and the model stays
  grey.
- Colour-by-fluctuation is measured *not* to be a second view of colour-by-
  displacement: r < 0.95 between them, or one of the two buttons should go.
- The overlay help moved into a topic of its own, `help_topics_views.py`, which
  is a real seam — every other topic explains part of the model, this one
  explains the seven things drawn on top of it, and they share a failure mode.

Suite 1399 → 1478 passing, 10 skipped.

### Follow-up: trimming the pore's escaped ends

The drawn pore included the slices past each mouth, where the probe balloons
into bulk solvent — on 11ZC, five spheres up to 12.2 Å across hanging under the
channel like a bulb. Faithful to the profile and misleading as a picture,
because it reads as the pore being longer and wider than it is.

**The first instrument I reached for was the wrong one.** Radial enclosure —
cast rays perpendicular to the axis, keep slices where most of them hit protein
— gives 8YEZ a median of 1.00 and 11ZC 0.75. That difference is almost entirely
that **11ZC is backbone-only**: with no side chains the rays slip between atoms.
It would have trimmed the one open structure hardest, for a reason that has
nothing to do with its pore. Not used.

The criterion that works needs no new measurement at all: it is the method's own
parameter. `pore_profile` tethers the probe centre within `pore.leash` of the
axis, and that tether is what makes the number mean "radius of *the pore*" —
unconstrained the probe escapes to R ≈ 6188 Å. Once the probe's *radius* exceeds
the leash, the tether has stopped localising anything, and that is exactly where
the profile stops describing a lumen.

**Only the ends, and that distinction is the whole design.** My first version
kept the contiguous in-leash run around the bottleneck. On 11ZC that removed the
5 bad slices at the bottom and **71 good ones off the top** — the upper vestibule
through the CED is genuinely over-leash in places and genuinely surrounded by
protein. A wide slice with protein on both sides is a vestibule, not an escape;
the profile alone cannot tell them apart, but their position can. So the leading
and trailing runs go and nothing else does. Measured: 5 slices on 11ZC, 4 on
7WLU, **0 on every closed entry** — 8YEZ, 7WLT and 8ZU3 never escape at all,
which is the check that the rule is not just shortening every pore it meets.

Display only. A trimmed slice is wider than the leash, so it can never have been
the minimum — driven over 200 random profiles rather than argued, because that
is the kind of argument that stays right until the rule changes. The profile
object, the plot and the bottleneck are untouched, and the status line says how
many slices went and that nothing moved.

---

## 2026-08-12 — Round 89: the family, and what a sequence can say about it

**Asked.** Find the other homologues, orthologs, isoforms and related PIEZOs and
add them so they can be compared. Decide whether the application needs BLAST
searches built in, and improve the alignment tool. Bring across what is useful
from the sister project at `../piezo_genes`.

### The family is nine, which is the answer to the BLAST question

One UniProt query — `reviewed:true AND protein_name:piezo` — returns exactly
nine entries and that is the whole reviewed family. This project held six. The
three missing are **rat Piezo1** (Q0KL00, most of the electrophysiology
literature is rat), ***Arabidopsis* PIEZO** (F4IN58) and ***Dictyostelium*
pzoA** (Q54S52). The last two matter most: PIEZO is not a metazoan invention,
and they take the generality question from one vertebrate duplication out to the
root of the eukaryotes.

A search tool is for when the answer set is unknown. It is not, and swapping a
pinned nine-row table for a live query would trade something reproducible for
something that is not. **No BLAST client**, recorded as a decision in
`docs/HOMOLOGY_SEARCH.md` in the form `NOT_PREREGISTERED_ROUND64.md` uses.

### What was missing was BLAST's *statistic*, and the measurement is stark

Aligning each pair, then re-aligning against **composition-matched shuffles of
the same partner**:

- **15 of the 36 pairs fall below Rost's 30% twilight line.**
- PEZO-1 against Arabidopsis PIEZO: **23.8% identical, and a scrambled
  Arabidopsis sequence of the same composition gives 22.5%** — z = 1.5. The
  percentage carries essentially nothing.
- The **local alignment score on the same pair is z = 64.** The homology is not
  marginal at all.

So the conclusion is not that these proteins are distantly related. It is that
**percent identity is the wrong statistic below that line**, which is what the
line means. `analysis.homology` refuses to report either statistic alone, and
`Relationship.verdict` says in words which one a pair is entitled to.

### Two instruments got it wrong before they were calibrated

`homology_sites` asks whether the curated gate, glutamates and PIP2 lysines
survive across the family. Its reliability gate was wrong twice:

1. **It used window identity** — self-contradicting, since the module next to it
   argues identity is the weak statistic. Now the BLOSUM62 window score.
2. **The width was chosen by taste and had no power.** At 31 residues the gate's
   z runs 1.3–2.6 in *every* non-mammalian member — refusing a mapping that is
   visibly right (human I2447/V2450/F2454 land on V2363/V2366/L2371 in the worm
   at a constant −84 register, and on F/L/F at offset 0 in the fly). A power
   scan across widths 21–201 put it at **101**, where the worm, fly and plant
   clear 3σ and Dictyostelium — correctly — still does not. The cost is stated:
   at 101 columns it answers *is this region in register*, not *is this residue*.

Split into `analysis/alignment_windows.py` at the 500-line limit, along a real
seam: that module is about a pair of sequences and has nothing to do with PIEZO.

### What the sites measure

- **The cap does not travel.** All three curated cap groups become unreadable
  outside the vertebrates. A limit on this project's own annotation, since every
  cap-gate distance in `liu2025_panels` is quoted from mouse residue numbers.
- **The gate erodes rather than switching off:** 3/3 identical in the mammals,
  2/3 in PIEZO2, 1/3 in worm/fly/plant. A gate made of a *property* rather than
  a contact would look exactly like that.
- **One group is identical everywhere it can be read, Dictyostelium included,
  and it is not the pore** — the anchor brake, human P2113/F2114.

A group with nothing readable returns `None`, not 0%. Perfect non-conservation
is the most confident possible way of saying nothing.

### The structural comparison, generalised

`paralogue.compare` refused anything that was not one PIEZO1 and one PIEZO2 —
right when PIEZO2 was the only other structure, wrong now four PIEZOs are
deposited. `analysis/homology_structure.py` removes that restriction and keeps
every guard, adding two `paralogue` never needed: it **renumbers first** (9W7X
is in a dPIEZO isoform's own numbering, +3 after residue 1570) and it **refuses
to pair helices by index** where the counts differ — 36, 38, 40, 35, 35, so all
four non-vertebrates. Reporting "2 of 38 agree" there would be a fact about
counting dressed as one about structure.

**And the first version of this reported a cherry-pick.** Comparing 7WLT with
9W7X gives a gating-mode overlap of **0.980** with dPIEZO — a striking sentence
at 30% sequence identity, and it was already written into the decision
document. Comparing **8YEZ** with the same 9W7X gives **0.189**. Same two
proteins, same method. Running the full 3x5 grid before believing it:

| partner | overlap over 3 PIEZO1 entries | beats control | |
|---|---|---|---|
| PIEZO2 (6KG7, 9VEE) | **0.80 – 0.98** | 6 / 6 | stable |
| PEZO-1 (9UOY, 9ZIS) | **0.18 – 0.98** | 5 / 6 | not stable |
| dPIEZO (9W7X) | **0.19 – 0.98** | 2 / 3 | not stable |

So `OverlapSpread` reports a range and a stability verdict, and the report entry
and the CLI show that rather than any single pair. PIEZO2 is the **positive
control** that makes this an instability rather than a broken instrument — the
method can say *stable*, and does, for the paralogue. This is the Round 85 rule
working in the direction it was added for.

The coverage-matched dome radius is 9.25 nm for mouse PIEZO1 against **9.24**
for PEZO-1, which does hold.

**The plant cannot join.** Its only structural representation is an AlphaFold
*monomer*, and Dictyostelium has not even that. So whether the dome is a property
of the fold rather than of animals **cannot be asked from structure at all** with
what exists. Reported as a gap in the world, not worked around.

### Defects this found

- **`fetch_alphafold` took `entries[0]`.** The endpoint returns one entry *per
  isoform*. For human PIEZO2 there is **no canonical model at all** — only
  isoform 3 (**709 aa**) and isoform 2 (2,689) — and `entries[0]` is the
  709-residue one, arriving as a well-formed mmCIF of the right protein that
  nothing downstream could question. For PEZO-1 the endpoint returns **twelve**
  entries and the canonical happened to be first, so this was right by luck.
  Now selected by exact accession and **refused** when there is no canonical
  model, naming the isoforms it was offered.

  This is precisely the lesson `piezo_genes` was built around — its own run
  flagged AlphaFold reporting 709 aa for PIEZO2 where every sequence database
  said ~2,800 — arriving in this codebase on the same protein.

- **`prediction_confidence` picked its model with `sorted(glob)[-1]`.** Right
  only because human PIEZO1 happened to sort last among two files. `AF-Q9H5I5`
  sorts *after* `AF-Q92508`, so adding the family would have silently changed
  the default to a different protein, read in human numbering. Named now.

- **`ModeComparison.summary` had "PIEZO2" as a literal** in both slots. The
  moment it could be handed a dPIEZO entry it printed an overlap with
  *Drosophila* PIEZO as an overlap with PIEZO2. Numbers right, sentence wrong,
  nothing raises.

- **4PKE/4PKX were excluded for the wrong stated reason.** The note said "a
  Piezo domain from a distant organism ... cataloguing it would need a seventh
  reference". Both halves wrong: the RCSB cross-references them to A0A061ACU2,
  so they are C. elegans PEZO-1, the same protein as 9UOY — and adding the
  reference does not help, because they are numbered in the expression
  construct's own coordinates (14-278, scoring 0.081) that no shift repairs.
  `EXCLUDED` now carries a reason per file.

- **`compare_sequences(method="positional")` had no numbering guard.** Harmless
  while the viewer offered only human and mouse; a live hazard the moment it
  offered nine, since pairing human 2447 with plant 2447 would report about two
  thousand confident substitutions. It raises now, and a structure chain carries
  the numbering `identify_numbering` measures rather than defaulting to human.

### The one that looked certain to break and did not

Rat Piezo1 is **94.2% identical to mouse** as a sequence, so adding it to
`REFERENCES` looked certain to collapse the margin `NumberingIdentity.confident`
requires and make eleven catalogue entries unreadable. It scores **0.066**
against a mouse entry. The identification is not a similarity test — it reads
each residue's name at its own *number*, and the twelve-residue length
difference puts everything past the first indel out of register. Measured on
7WLT, 3JAC and 6B3R, whose margins are unchanged to three decimals, and pinned.

### Also

Catalogue 28 → 34 entries. Added the two missing PEZO-1 isoform structures
(9ZIS isoform g, 9UOX isoform k — giving the catalogue its **only replicate
pair**, which is what lets an isoform difference be told from a dataset
difference), the two human PIEZO2 entries, dPIEZO, and four AlphaFold models.
Five new verified references (BLOSUM62, Rost's twilight zone, Smith–Waterman,
neighbour-joining, Felsenstein's bootstrap); 85 total. Eight new registered
parameters. The sequence viewer now offers all nine PIEZOs where it offered two.


---

## 2026-08-12 — Round 89b: the annotation gap the family made visible

**Reported.** Selecting *plant* in the structure filter crashed the application
with `KeyError: 'plant_piezo'` in `colormaps.lookup_table`.

**And it was not mine, though my change is how it was found.** `domains.json`,
`variants.json` and `functional_residues.json` are curated in human and mouse
PIEZO1 numbering and in nothing else. The registry already carried
`worm_piezo`, `fly_piezo`, `human_piezo2` and `mouse_piezo2` **before** Round
89 — since Round 83 — so 6KG7, 9VEE, 9VEF, 9UOY, 9ZIT and 9W7X had all been
uncloseable in the GUI the whole time. Round 89 added the plant entry, and the
plant filter put a new numbering one click away.

**Two faces, and the second is much worse than the crash.**

`DomainPalette` indexed `d[self.species]` straight into the domain records, so
colouring by domain — the default colouring — raised on any non-PIEZO1 entry.
Loud, and at least unmistakable.

`Annotations` built each functional-residue group with `r["human"] if species
== "human" else r["mouse"]`. **Every non-human numbering therefore got mouse
PIEZO1 residue numbers.** Selecting the hydrophobic gate on a *C. elegans*
structure highlighted whatever sits at mouse 2473, 2476 and 2480. Nothing
raised, the picture looked right, and the residues were a different protein's.

**The fix is one line stated in one place.** `ANNOTATED_NUMBERINGS = ("human",
"mouse")` in `core/annotations.py`, with `is_annotated()` and
`annotation_gap()` beside it. An unannotated numbering yields **nothing** —
checked before any file is read, because there is no partial answer to give —
and the reason travels on the object to the status line, which now leads with
**NO ANNOTATION FOR THIS PROTEIN**. That last part is not decoration: an empty
domain list drawn as uniform grey reads as *this protein has no domains*, which
is the opposite of what is true.

**And chasing it found one that was mine.** `piezo1_numbering()` gated on
`PIEZO1_REFERENCES`, and Round 89 added rat to that tuple — correctly, since a
rat entry should be identified as rat rather than mis-read as mouse. But every
caller takes its return value straight to `load_annotations`, so the rat
AlphaFold model would have reached the component selector, the conduction path
and the pore charge map as a numbering they find nothing in, and reported an
empty result rather than a refusal. It gates on `ANNOTATED_NUMBERINGS` now.
*Which protein is this* has nine answers; *can I read annotation into it* has
two, and conflating them is how a refusal turns into a null result.

Verified by loading all six previously-crashing entries through a real
`MainWindow` with a GL context, and pinned in `tests/test_annotation_coverage.py`
by driving `domain_colors` — the exact call `rebuild` makes — over one entry per
catalogued numbering.


---

## 2026-08-12 — Round 89c: a trimer from one protomer, and what it is worth

**Asked.** When a structure has only one protein, offer to display a trimer
assembled from it based on similar structures.

**Why it matters here.** Half this project needs three protomers — the dome
fit, the pore profile, the elastic network, the paralogue comparison all take
three blocks and refuse anything else. Round 89 recorded the cost as a gap in
the world: the only structural representation of a non-animal PIEZO is an
AlphaFold **monomer**, so *is the dome a property of the fold rather than of
animals?* could not be asked at all.

**Built** (`structure/assembly.py`): place the monomer onto each protomer of a
deposited trimer. Correspondence through a real alignment, because the point is
a plant protein on a mouse template and no two PIEZOs share a length. Each
protomer placed independently rather than by rotating one about the template's
C3 axis — that would give a trimer whose C3 deviation is a constructed zero,
indistinguishable on screen from a measured one.

**Two wrong versions before a right one, both caught by the numbers.**

A global superposition placed human PIEZO1's *own* AlphaFold model onto 6B3R at
**19 Å**, and the plant at 25. The fit is dominated by the distal blade — the
lesson `hybrid` already learned. So: fit the rigid core, found by outlier
rejection rather than by naming residues, which matters because the protein this
exists for has no curated helix ranges to name one with.

Then rejecting a fixed *fraction* per cycle drove the core to its 200-residue
floor on every entry — 200 of 2,500 fitted to 1.2 Å, which is not a core but
the 200 that agree best and always exists. A *distance* criterion alone never
starts, because from a 19 Å fit nothing is within 3 Å. Both: descend by
fraction until the distance criterion can bite, then converge on it. Cores are
now 200–948 residues at 1.3–2.9 Å, and `n_core` is a measurement — rat 948 on
6B3R, the plant only 200, which is flagged **AT THE FLOOR**.

**Calibrated on an exact known answer.** Pull chain A out of 6B3R, hand it back
as a monomer, rebuild against 6B3R: **0.00 Å over all 1,502 residues and 8
clashes — 6B3R's own count**. Rebuilt against 7WLT instead: 1.34 Å and 36
clashes. Both halves needed, or 0.00 Å would be equally consistent with a
function that returns its input. The clash counter is calibrated the same way,
and had to be: assemblies score thousands, and that means nothing until a real
trimer's score is known — 6B3R 8, 7WLT 3, 9ZIS 6.

**And measuring what it is worth is the result, which is less than hoped.**
`borrowed_fraction` splits the assembly's departure from planarity using
`structure.planarity`'s existing decomposition: **79% worm, 83% plant, 96%
rat** of it is the template's arrangement. So a dome radius measured on an
assembly is mostly a measurement of the template, and this does **not** answer
the question it was reached for. It narrows "cannot be asked from structure at
all" to "can be asked, and would be 83% about 9ZIS" — a sharper statement of
the same gap rather than a way round it. Said in the caveat, on the status
line, and pinned in a test.

Reachable as **Completeness → Assembled trimer (MODELLED)**, which is the right
home because it is the existing mechanism for "whatever this produces is what
every analysis runs on, without any of them knowing" — and therefore the one
that already understands that provenance cannot be optional. The amber HUD
banner reads *MODELLED TRIMER — arrangement taken from 9ZIS, not measured*.

---

## 2026-08-12 — Round 90: rendering-style controls for the features that had none

**Asked.** More controls for choosing the type of rendering of structural
features — including the ability to change the HaloTag to ribbon etc., and the
other structural elements that were inaccessible.

**What was actually inaccessible.** The primary structure has had a style
selector since the beginning (Model panel) and the superposition overlay since
the Overlay panel existed. Everything else drawn as *structure* was hard-coded
to one representation: the HaloTag fold (half-vdW sphere cloud), the
full-length graft (1.6 Å sphere cloud), companions (backbone), the component
highlight (gold ball-and-stick), and the resolved lipids/ligands (vdW spheres,
with only a visibility toggle).

**The rule that shaped the design: restyling moves no caveat.** The features
being opened up are exactly the ones this project draws most carefully — a
fold whose orientation is undetermined, a graft that is 48% low-confidence
prediction. A cartoon looks *more* like a determined pose than a sphere cloud
does, so the guards matter more after restyling, not less. Concretely:

- `MolecularView` gained `color_override`, a per-atom RGB that wins over every
  palette. The fold's contact-red atoms and the graft's pLDDT bands are the
  visible half of reported numbers; a style change must carry them, not hand
  the colouring to a palette that knows nothing of them.
- The **fold styles** (`FOLD_STYLES`) build a real `Structure` of the three
  placed tags — one chain per copy, so bonds and cartoon traces stay within a
  tag — and drive the same `MolecularView` machinery that styles the channel.
  The default is the pinned sphere cloud, bit-for-bit: `test_ui_fusion`'s
  batch layout is untouched, and the UNDETERMINED status line is asserted
  across **every** style. The radius-of-gyration sphere is not restylable,
  because it is a statement about what the model determined, not a preference.
- The **full-length model** (`HYBRID_STYLES`) can be a backbone tube or trace;
  `HybridModel` now records its C-alphas (`ca`) so the controller does not
  re-derive atom names from the two source files. The ribbon reads its colours
  from the same per-atom array the sphere cloud uses, so grey stays grey and
  the bands stay bands.
- **Companions** share one selectable style (persisted; backbone default) —
  one for all rather than one each, because companions are told apart by
  colour, and a mixture of styles would hand that job to shape as well.
- The **component highlight** offers ball-and-stick, sticks and vdW spheres,
  all gold, and none changes which residues are highlighted.
- **Ligands** get a style combo in the Model panel (spheres / balls /
  ball-and-stick), with their bonds in a separate batch so hiding or clearing
  cannot strand them.

**Where it lives.** `menus.py` was at exactly 500 lines, so the new submenus —
and the HaloTag submenu itself, which now carries its style group — moved to
`ui/menus_styles.py`, the same split-at-the-seam `menus_flux.py` used: every
entry there chooses what is computed, every entry here chooses how a feature
is drawn.

**Tests** (`test_ui_feature_styles.py`, 21): the default fold style leaves the
pinned batches alone; every fold style keeps the caveat and the seam; the
contact-red count survives restyling exactly (`body_contacts × n_tags`); the
hybrid ribbon must show the flat grey *and* more than one pLDDT band; the
highlight styles change radii and nothing about which atoms; unknown keys are
refused everywhere; `color_override` is checked at the layer that implements
it. The existing pinned suites — `test_ui_fusion`, `test_ui_hybrid`,
`test_ui_companions`, `test_ui_components`, `test_ui_controls` (which fires
the new menu actions against a real window), and the pixel-counting render
suites — all pass unchanged, which is the point: nothing that was pinned
moved.

---

## 2026-08-12 — Round 90b: tiered test runs, with the full suite as the occasional check

**Asked.** The full suite takes a long time now — shorter versions for
different situations, keeping the full suite as an occasional check.

**The trap in a shorter suite.** A short run is a *selection*, and a selection
can silently drop things: a test that no situation runs would decay without
anyone deciding that, and the only symptom of a selector that quietly runs
everything is the time the user was trying to save. Both failure modes are
invisible in a green run, which is what shaped the design:

- **The tiers are a partition** (`tests/tiers.py`): five situations — `quick`
  (sanity on every edit: imports, registries, resources, pure logic; no Qt, no
  GL, no heavy computation), `science` (physics/structure/analysis on real
  coordinates), `ui` (offscreen Qt suites), `render` (real OpenGL, judged in
  pixels), `records` (documentation, frozen records, reproducibility) — and
  every test file belongs to exactly one. `test_tiers.py` fails on a file in
  no tier or two, so assigning a new test file is enforced, not remembered.
- **Selection is by deselection** (`--suite` in `conftest.py`), so a short
  run's *real* skips — the ones that mean data is missing — stay visible
  instead of drowning in a forest of tier skips. An unknown tier name is a
  usage error, because a typo that selected nothing would report a green run.
- **The selector is calibrated two-sided**: a science file must vanish from a
  quick collection *and* a quick file from a science collection, or a hook
  that selects everything would pass half the check.

**Measured.** `make test-quick` runs 517 tests in ~50 s (measured while the
full suite was still occupying the machine, so an upper bound). The Makefile
targets: `test-quick`, `test-science`, `test-ui` (both include quick),
`test-render`, `test-records` — and `make test` unchanged as the full,
occasional check. The situational tiers together *are* the full suite, by the
partition, so nothing can fall between them.

---

## 2026-08-12 — Round 90c: a click answers for everything drawn

**Asked.** Select atoms in any structure — as for the main structure: the
HaloTag, ligands, and the rest.

**What a click could and could not do.** The pick source was the primary
structure's atom array and nothing else. The HaloTag, the extra structures
and the full-length graft were mute — worse, a click aimed at them
identified whatever primary atom lay *behind* them, which is a wrong answer
delivered confidently. And a click on a lipid worked but then looked the
lipid's author-assigned residue number up in the **curated protein
annotation**: mouse numbering runs to 2547, the lipids' numbers land inside
it, and the status line named a domain the lipid is not part of.

**Two rules carry the feature.**

- **Nearest wins, whatever drew it.** `nearest_hit` in `gl_widget.py` takes
  named coordinate sets and answers "what did I click" by the one honest
  rule: the thing in front. Pure geometry, Qt-free, calibrated on rays whose
  answer is known by construction. Controllers register their drawn atoms
  (`register_pick_feature`) and unregister when they clear, so a click can
  never identify something not on screen — which is also why `load_structure`
  now clears the full-length overlay: it survived a structure change, drawn
  over the new entry, and was one registration away from *answering clicks*
  for it.
- **A feature identifies as what it is.** The describe text is part of the
  registration, and it must say what the thing *is*, not only which atom: a
  tag atom answers MODELLED with the spin UNDETERMINED, a graft atom answers
  PREDICTED with its pLDDT (its experimental half answers experimental), a
  companion names itself and adds that the analyses run on the primary. A
  tag atom identified like a deposited one would be the confident wrong
  answer the rest of this project spends its guards on.

An armed measure click on a feature **refuses out loud** — swallowing the
click silently would break inspection, and a distance to a modelled position
would be a measurement of a guess. A feature pick marks the atom exactly as
a primary pick marks its residue, one current selection at a time in either
direction.

**Tests** (`test_ui_picking.py`, 15): the geometry both ways round (a nearer
feature beats the primary *and* a nearer primary beats the feature), the
identify-as-what-it-is texts for all three features, the loud refusal, the
marker dropping with its source, and the HETATM fix pinned on a real 7WLT
lipid. The `ui` and `render` tiers pass, `test_ui_controls`' pinned
click-versus-drag behaviour among them — one situational run each, which is
what Round 90b was for.

---

## 2026-08-12 — Round 91: what survives a load, and what a click may claim

**Asked.** Research further improvements and implement them. The roadmap has
no open items, so the research was an audit of the seams the last three
rounds exposed — and it found one bug *class* with five live instances, plus
two picking gaps of the class Round 90c fixed.

**The stale-overlay class.** `load_structure`'s clear list has a history of
omissions: the morph was missing until Round 87, the full-length overlay
until Round 90c. This round found four more — the **micelle**, the **planar
membrane**, the **potential colouring**, and worst, the **ion stream**: left
running across a load it keeps animating the old entry's ions along the old
entry's path over the new structure, with the old current on the HUD — over
8YEZ, the very entry the wetting model refuses to animate. And the
open-a-file path (`_open_file`) had a two-entry copy of the list, so *every*
overlay it lacked survived opening a file; it also carried the documented
overlay-resurrection ordering bug live in a second place, and left
spliced-model state that could put the amber PART PREDICTED banner over a
deposited file.

**The fix is one list** — `_clear_structure_overlays`, shared by both entry
points — **and the guard is not a list**, because lists are what decayed.
`test_ui_load_hygiene` holds both paths to an *equivalence*: loading B after
using A must leave the same scene as loading B fresh, whatever was drawn in
between. Overlays are switched on by discovery (has `show`, and
`clear`/`reset`), with a ratchet on the discovery count so the sweep cannot
quietly find nothing. Calibrated the way this project calibrates checkers:
with `micelle.clear()` deliberately removed, the guard fails naming
`micelle_envelope` as the survivor.

**Picking follows visibility.** The pick source was the full atom array, so
a hidden entity category — or everything a component hides — kept answering
clicks: click where an invisible lipid is, get told about the invisible
lipid instead of the visible helix behind it. `MolecularView.pickable_mask`
now states which atoms are on screen (entity filter, component residues, the
ligand toggle — which silences ligands in ribbon styles and not in atom
styles, because that is what is drawn), every visibility handler refreshes
the viewport's pick mask, and `nearest_hit` takes the mask with the rule
sharpened to **nearest visible wins**. The unmasked control is asserted
beside the masked case, so the test shows the mask doing the work.

**The right-click joins the rule.** It resolved against the primary only, so
right-clicking through a drawn tag named the occluded residue behind it —
the menu's residue entries are annotation read by primary atom index. It now
routes through `hit_at`: a feature in front opens the generic menu rather
than describing an atom the user did not aim at.

**Hygiene.** `fusion_controller.py` (489 lines) split along the
lifecycle-versus-geometry seam: `fold_view.py` owns drawing and labelling
the placed fold, as pure functions of the pose. `_open_file` gained
`show_opened_structure`, split from the dialog so the whole path can be
driven in a test — which is how its gaps were found.

**Not done, recorded rather than dropped:** session persistence for the
Round 90 style controls (fold, companion, hybrid, highlight, ligand styles
are not saved in sessions). A candidate for a future round; it needs
`session.py`'s format guards extended, not just keys added.

---

## 2026-08-13 — Round 92: viewer appearance options, and one rule for the Options menu

**Asked.** Options for the viewer — background colour, GUI style — and all
options consolidated into the Options menu, each tested to make an actual
change.

**The consolidation needed a rule, not a tidy-up.** Four persisted
preferences had accreted under View — the structure alignment mode, the
multi-structure toggle, the companion style and the display-options dialog —
so where to change a remembered behaviour depended on which menu a feature
had happened to be added to. The rule, stated in `menus_options.py` where
the menu is built: **Options holds what is remembered across sessions; View
holds what is shown right now.** Two things that look like options
deliberately stay put, and the test pins them staying: the ion-flux pathway
and voltage change what is *computed* (the reasoning recorded in
`menus_flux.py` when they were placed), and the per-feature style submenus
are choices about a feature, made while looking at it.

**The two new options.**

- **Viewport background** — five steps from midnight to white, no colour
  picker: the useful question is dark room versus manuscript figure. The
  background is also the depth-cue fog colour, and the two share one
  settings object precisely so they cannot disagree — a background change
  that left the fog behind would haze every atom toward the old colour, and
  a test reads the fog uniform beside the clear colour. The "default" entry
  is asserted at import to *be* the `RenderSettings` default, because two
  copies of one colour drift. The scale bar and readouts already carry dark
  halos, which is what keeps them legible on white.
- **Interface theme** — dark (default), light, system. The existing dark
  stylesheet became one template over a per-theme token table, so the
  themes cannot drift apart in structure and a missing token raises at
  import instead of borrowing the other theme's colour. Applied at startup
  before the window builds (no dark flash for a light-theme user), on
  change, and on Restore-defaults — a reset that left a white viewport or
  light chrome behind would not be a reset.

**Tested to the pixel** (`test_ui_options.py`, 9, render tier): the
background choice is verified in the rendered framebuffer's corner pixel
(white > 240, midnight < 40) and in the fog uniform; the theme in the
application palette's lightness and the stylesheet's presence or absence;
every moved option through the state it stores; every simple option through
its observable effect (hint visibility, spin speed, focus mode, layout
memory); and the consolidation itself as a menu-placement guard. The
fixture snapshots and restores every setting and the application style, so
running the suite cannot restyle the user's next session. The scripted GUI
smoke test passes over the refactored menus — the check that exists because
mechanical Qt refactors have broken the app silently twice.

---

## Round 93 — the PIEZO family census, imported and answered back

**What this round is.** A sibling project, `piezo_genes`, has spent twenty-four
sessions on a 194-genome, eukaryote-wide census of the PIEZO family: what its
true range is, that vertebrates have a **third** PIEZO gene the databases
largely missed, and which parts of the protein half a billion years of evolution
has refused to change. It is a sequence project — no coordinates, no physics,
and no way to ask *why* a residue is conserved. This project is the mirror
image. The round brings its results across and asks what this project's
coordinates say back.

**Why an import needed a gate rather than a copy.** The failure mode for
somebody else's result is specific and silent: they keep working, correct a
number, and the copy here becomes a confident quotation of a superseded value.
So `scripts/build_family_findings.py` makes every imported statement declare how
to re-read the number it rests on, re-reads all thirty-two from the census's own
result files, requires claims never reduced to a table to appear literally in
its `FINDINGS.md`, records the source commit — and **refuses to write at all**
when the census is not on disk, because a rebuild with no source would otherwise
re-stamp the resource with a fresh date and nothing behind it.

The gate that earns its place is the fifth one. The census works in three
numbering systems and this project works in two others, and a per-residue track
joined to the wrong sequence produces a colouring that looks entirely plausible
and is off by an indel. So the imported PIEZO1 track's own amino acids are
checked residue by residue against `uniprot_human.json` before anything is
written.

**Four replications, on boundaries the census did not choose.** The point of
re-running rather than quoting is that it can fail. The two projects' domain
partitions put the **anchor 141 residues apart** and the outer helix 120, while
agreeing on the cap to within four — the signature of two ranges taken from
different papers. On ours: anchor 0.832, CTD 0.810, inner helix 0.789 against
THU1's 0.630, so the core-is-pore ordering is a property of the protein. The
pairwise identities recomputed from this project's own alignment put the cap —
the census's one exception, the piece of pore machinery *below* the
whole-protein figure — at **0.404 against their 0.402**. The fourteen pathogenic
pore positions are confirmed on a *different UniProt record* of the same
zebrafish gene from the one the census scored. And two evolutionary measurements
sharing no data and no statistic agree at ρ = 0.88.

**Four places the answer came out differently.** None overturns a headline; two
sharpen one and two say what a finding is really about.

- **The distal-versus-proximal blade gradient is composition, not biology.** The
  census's bands reproduce here exactly, so the import is sound — and the bands
  are a single chain cut containing 29% and **77%** inter-unit linker
  respectively, with linker scoring the same either side (0.517 against 0.515).
  Restricted to the four-TM units the ordering **reverses**. What survives the
  composition is the opposite gradient.
- **The disease enrichment is real and it turns on 120 residues.** Re-tested on
  PIEZO1 alone against gnomAD *population* missense rather than ClinVar benign
  labels — a better control for the ascertainment problem the census names in
  its own caveat — it replicates strongly on their boundaries (OR 3.63,
  P = 0.0033 against their 3.9 and 0.0014) and does not reach significance on
  ours (OR 1.60, P = 0.25). The disputed band is 2057–2176, and it carries six
  pathogenic positions including E2117 and T2127. Both rows are printed, because
  picking the one that agrees would be reporting a boundary choice as a finding.
- **The blades "splaying" is the prediction, not the paralogue.** The census's
  structural result superposes a *predicted* piezo3 monomer on experimental
  mouse Piezo1 and reports the cores agreeing while the blades splay. Generalised
  here to every pair, with the control that decides it: an AlphaFold monomer
  splays **7–9× from an experimental structure of the protein it is a model
  of**, while two experimental structures of *different paralogues* splay
  **0.8–1.2×**. The core agreement survives and is stronger than reported — two
  experimental paralogues superpose at 3.6 Å against the census's 3.9 — and the
  splay is an artefact of comparing a model with an experiment. The 19× at the
  top of that table is PIEZO1's own gating motion: core-conserved and
  periphery-free is what flattening looks like inside one protein.
- **`best_template` picks the wrong template for a new paralogue.** Its rule —
  same protein first, then most residues resolved — falls through for a
  paralogue nobody has a structure of, to the **worm** PEZO-1 entry at 28%
  identity and 13,839 inter-protomer clashes against a PIEZO1 trimer's 44% and
  2,922. Right for the case it was written for. `analysis.piezo3` chooses
  explicitly and prints the comparison rather than overriding it silently.

**The joint question neither project could ask.** Is a residue's evolutionary
constraint predicted by how mechanically coupled it is, or only by how buried it
is? Buried residues are conserved in every protein ever studied for reasons
having nothing to do with mechanotransduction, and burial correlates with almost
every mechanical quantity — so the second half of that sentence is the whole
difficulty. Three controls: the null is a **circular shift, not a permutation**
(both series are strongly autocorrelated along the chain, and a permutation null
is measured to be three times narrower — which is how a comparison invents
significance from structure that was in both series to begin with); burial is
partialled out, rank-first because the confounder is monotone and not linear;
and eight features are corrected together. Measured on 7WLT: PRS response at the
gate ρ = 0.373, holding **0.287** with burial fixed at q = 0.007, against
burial's own 0.369. Five of eight survive all three controls, and the signs are
the census's picture in mechanics — coupled to the gate ⟹ constrained, mobile ⟹
free.

**A defect this found, and one it did not fix.** That result only appeared after
`analysis/features.py` was corrected. It defaulted to *human* annotation
whatever entry it was handed, so on a **mouse** entry — most of the catalogue —
the hydrophobic-gate group, the blade range and the human-anchored conservation
profile were all looked up at human residue numbers against mouse coordinates,
a non-constant offset reaching 26 residues. Measured before and after on 7WLT:
the conservation column scored **ρ = 0.29** against the same profile read at the
correct residue and now scores 1.00; `distance_to_gate` had no residue at zero
and now finds mouse 2473/2476/2480, which is what `functional_residues.json`
says the gate is; `prs_gate_response` against constraint went from **−0.02 to
+0.37**. It survived because `tests/test_features.py` uses a *human* fixture,
where the default is right, and the actual use is mouse.

`docs/PREREGISTRATION_ROUND48.md` records that Round 48's endpoints — burial,
conservation and gate coupling — were built by `build_feature_table` on 7WLT.
Round 48 is a frozen null and is **not revised here**: the standing policy in
`NEGATIVE_RESULT_PROTOCOL.md` is that a recorded result is superseded rather
than edited, and superseding it is a decision about the project's central claim
rather than a side effect of an import. It is recorded in `docs/FAMILY.md` §6
rather than left to be discovered.

**piezo3 as a structure.** The only coordinates the third vertebrate PIEZO has
anywhere are one AlphaFold model of the zebrafish protein — human piezo3 has
been the pseudogene `PIEZO1P2` since before the primate radiation. It is now a
tenth family reference, a catalogue entry, and identifies at 1.000 against a
runner-up of 0.068 **with no other entry's identification moving**, which is
pinned the way rat's addition was. Run through the pipeline it gives dome
R = 10.8 nm and a 0.37 Å bottleneck against 7WLT's 9.7 nm and 0.73 Å by the
identical route — and neither number is evidence about piezo3, because 96% of
the assembled trimer's departure from planarity is the template's arrangement
and the protomer is a prediction whose blades this round has just measured to
sit 33–45 Å from where cryo-EM puts them. What the numbers *can* do is fail, and
they did not: the protomer arranges into a closed trimer with an axis and a
continuous lumen. A negative that survived, not a positive demonstrated.

**Calibration.** Eight new modules join `CHECKING_MODULES`, because none of them
measures a property of PIEZO1 and all of them exist to decide whether something
may be believed. Every public callable has a named calibrating test — and the
register earned its keep immediately, failing on two tests I had named and not
written.

Two calibrations were wrong first and are recorded as such. A planted
correlation on two random walks reached only z = 1.6, which looked like the
instrument failing and is the instrument being correct: two random walks
correlate at |ρ| near 1 by chance, so a shift null built from one is enormous.
Pinned as a bound rather than deleted. And the rank-partialling test first
claimed ranking removes more of a *skewed* confounder, which is false when the
dependence is linear in the raw values — the case that decides it is a monotone
**non-linear** confounder, which is what burial is, and the test now measures
the ranked control against the raw one rather than asserting the choice.

**Six numbers registered, three exempted with reasons.** The census's own band
boundaries (mouse 1300, its human equivalent, and where its proximal band ends)
are exempt rather than registered: they are quoted to reproduce somebody else's
number, and an override would silently stop the reproduction being one.
