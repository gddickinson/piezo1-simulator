"""Cells for the first two notebooks: what is in a structure, and how it moves.

Kept as Python rather than JSON so the content can be reviewed and diffed like
anything else. `build_notebooks.py` runs every code cell before publishing.

Written to be read top to bottom. Each notebook ends with `assert`s on the
numbers it quotes, so running it checks the science and not just the syntax.
"""

from __future__ import annotations

__all__ = ["NOTEBOOKS"]


_FIRST_LOOK = [
    ("markdown", """
# 1 · A first look at PIEZO1

What is actually in a deposited structure, how to put it in a frame you can
reason about, and how to measure the membrane dome from the coordinates.

**You need the data first:**

```bash
python -m piezo1.io.fetch
```

Every cell below runs in a few seconds.
"""),
    ("code", """
from piezo1.config import STRUCTURE_DIR
from piezo1.core.structure import Structure

st = Structure.from_file(STRUCTURE_DIR / "8YEZ.cif")
print(st.n_atoms, "atoms |", st.n_residues, "residues |", st.chains, "chains")
"""),
    ("markdown", """
## What is in the file, and what is the channel

A deposited entry is not only the protein you asked for. It carries lipids,
detergent, glycans, ions, water — and sometimes a whole other subunit. Six
PIEZO entries include three copies of **MDFIC**, a 21-residue auxiliary
protein whose residue numbers (226–247) sit *inside* PIEZO1's own range, so a
selection by residue number alone would silently mix the two.

`classify` sorts every atom, and the analyses use the channel protomers
whatever else is present.
"""),
    ("code", """
from piezo1.core.entities import classify

entities = classify(st)
for name, count in sorted(entities.counts().items(), key=lambda kv: -kv[1]):
    print(f"  {name:16s} {count:6d} atoms")
print("\\nchannel protomer chains:", entities.protomer_chains)
print("auxiliary chains        :", entities.auxiliary_chains or "none")
"""),
    ("markdown", """
## Put it in a frame you can reason about

Deposited entries sit wherever the depositor left them — up to 147 Å apart, and
some of them upside down. Every geometric number in this project (dome
curvature, pore profile, tag distance) is quoted relative to the three-fold
axis, so the structure has to be framed before any of it means anything.

`canonical_transform` finds the channel's own C3 axis, puts it on **+z**, and
puts the cytosolic side at **−z**.
"""),
    ("code", """
from piezo1.structure.frame import apply_frame, canonical_transform

transform = canonical_transform(st)
framed = apply_frame(st, transform)

print("mode              :", transform.mode)
print("C3 axis fit RMSD  :", round(transform.axis_rmsd, 3), "A")
print("note              :", transform.note)

# The cytosolic end must now be at negative z. Every distance quoted against
# the conduction axis depends on that sign being right.
#
# Select by residue NUMBER, not by position in the array. Taking "the last N
# rows" straddles chains and gave the wrong answer on 7WLU and 11ZC while
# still reporting a perfect C3 fit — a structure upside down in the viewport
# is obvious, one upside down inside a calculation is not.
mask = framed.mask_ca() & ~framed.hetero
seq, xyz = framed.res_seq[mask], framed.xyz[mask]
cterm = xyz[seq >= seq.max() - 40][:, 2].mean()
nterm = xyz[seq <= seq.min() + 40][:, 2].mean()
print(f"C-terminal (cytosolic) end: z = {cterm:+6.1f} A")
print(f"N-terminal end            : z = {nterm:+6.1f} A")
assert cterm < 0.0, "the cytosolic end must be at negative z"
"""),
    ("markdown", """
## Measure the membrane dome

PIEZO1's blades bend the membrane into a dome, and the curvature of that dome
is the mechanism (Guo & MacKinnon 2017). The measurement fits a sphere to the
mid-membrane surface and reports the radius of curvature.

The regression case is curved mouse Piezo1 in a bilayer, **7WLT**: the code
should return about 9.7 nm against a published 10.2 nm (Haselwandter &
MacKinnon 2018).
"""),
    ("code", """
# ANALYSES is the shared registry the GUI and the command line both dispatch
# through, so a notebook using it cannot drift from what the application shows.
from piezo1.analysis.report import ANALYSES

curved = Structure.from_file(STRUCTURE_DIR / "7WLT.cif")
dome = ANALYSES["dome"](curved, "mouse")

for key, value in dome.items():
    if isinstance(value, float):
        print(f"  {key:24s} {value:8.2f}")
print(f"\\n  reference: {dome['reference']}")

assert 9.0 < dome["radius_of_curvature_nm"] < 10.5, dome["radius_of_curvature_nm"]
"""),
    ("markdown", """
## Curved against flat

The same measurement on the flattened state separates the two clearly. This is
the transition the whole project is about: tension flattens the dome, and the
blades lever the pore open.
"""),
    ("code", """
for pdb, species, state in (("7WLT", "mouse", "curved, bilayer"),
                            ("7WLU", "mouse", "flattened"),
                            ("11ZC", "mouse", "flat, native vesicle")):
    entry = Structure.from_file(STRUCTURE_DIR / f"{pdb}.cif")
    d = ANALYSES["dome"](entry, species)
    print(f"  {pdb}  {state:22s} R_c {d['radius_of_curvature_nm']:5.1f} nm   "
          f"depth {d['dome_depth_nm']:4.1f} nm")
"""),
    ("markdown", """
## The trap that will bite you: residue numbering

Most PIEZO1 **mechanism** papers number by mouse Piezo1 (2,547 residues). Most
**disease** papers number by human PIEZO1 (2,521 residues). The offset between
them is **not constant** — it runs from 0 to +26 across twelve blocks and
passes through zero twice.

Never subtract a constant. Always go through the alignment.
"""),
    ("code", """
from piezo1.core.sequence import human_to_mouse

for human in (1718, 2456, 2496, 756):
    mouse = human_to_mouse(human)
    print(f"  human {human} -> mouse {mouse}   (offset {mouse - human:+d})")

offsets = {human_to_mouse(h) - h for h in range(600, 2500, 25)
           if human_to_mouse(h) is not None}
print("\\ndistinct offsets across the chain:", sorted(offsets))
assert len(offsets) > 1, "a constant offset would mean the map is not being used"
"""),
    ("markdown", """
## Where to go next

* `02_gating_motion` — the elastic network model, and the symmetry rule that
  says which motions can couple to membrane tension.
* `03_pore_to_current` — is the pore open, would water stay in it, and what
  current would flow.
* `04_variants_and_the_null` — the variant workflow, and the result that did
  not work.

`docs/NOTEBOOK.md` is the full API reference, with a table of the things that
will bite you.
"""),
]


_GATING = [
    ("markdown", """
# 2 · The gating motion, and why symmetry matters

PIEZO1 opens when membrane tension flattens its dome. This notebook builds an
**elastic network model** from a single closed structure and asks whether the
motion the protein makes most cheaply is the motion it actually makes.

It is a real test, because the answer is checkable: two deposited structures
give the observed transition, and the model never sees the second one.
"""),
    ("code", """
import numpy as np

from piezo1.config import STRUCTURE_DIR
from piezo1.core.structure import Structure
from piezo1.physics.anm import ANM
from piezo1.structure.protomers import protomer_blocks

curved = Structure.from_file(STRUCTURE_DIR / "7WLT.cif")
flat = Structure.from_file(STRUCTURE_DIR / "7WLU.cif")

_, curved_res = protomer_blocks(curved)
_, flat_res = protomer_blocks(flat)
common = np.array(sorted(set(curved_res.tolist()) & set(flat_res.tolist())))
print(f"{len(common)} residues resolved in all six protomers")
"""),
    ("markdown", """
## Three things to get right before any of this means anything

Comparing two structures is where most of the mistakes live, and all three of
these produce a plausible wrong number rather than an error.

1. **Compare like with like.** The two entries do not resolve the same
   residues, so both are resampled onto the set they share. Skip this and you
   are comparing different molecules.
2. **Do not trust chain labels.** They do not encode rotational order.
   Overlaying 7WLU on 7WLT by label gives 90.7 Å; searching the protomer
   correspondence gives 12.3 Å. `match_protomers` searches.
3. **Remove rigid-body motion.** A structure that has merely been *moved* would
   otherwise look like a huge conformational change. Kabsch superposition
   takes it out, leaving only the shape change.
"""),
    ("code", """
from piezo1.structure.superpose import kabsch, match_protomers


def resample(structure, residues):
    \"\"\"The same residues, in the same order, from each protomer.\"\"\"
    out = []
    for chain in structure.chains:
        mask = structure.mask_ca() & (structure.chain == chain)
        if mask.sum() < 300:
            continue
        index = {int(r): i for i, r in enumerate(structure.res_seq[mask])}
        xyz = structure.xyz[mask]
        out.append(np.array([xyz[index[r]] for r in residues], dtype=float))
    return out[:3]


curved_blocks = resample(curved, common)
flat_blocks = resample(flat, common)

match = match_protomers(curved_blocks, flat_blocks)
flat_blocks = [flat_blocks[i] for i in match.order]
print("protomer order found:", match.order, " RMSD", round(match.rmsd, 1), "A")

rotation, translation, centroid = kabsch(np.vstack(flat_blocks),
                                         np.vstack(curved_blocks))
fitted = (np.vstack(flat_blocks) - centroid) @ rotation.T + translation
displacement = (fitted - np.vstack(curved_blocks)).ravel()

rmsd = float(np.sqrt((displacement ** 2).reshape(-1, 3).sum(1).mean()))
print(f"shape change after removing rigid motion: {rmsd:.1f} A RMSD")
"""),
    ("markdown", """
## Build the elastic network, and label every mode by symmetry

The channel is a trimer, so every normal mode transforms as one of the
irreducible representations of C3: **A** (symmetric under 120° rotation) or
**E** (a degenerate pair that is not).

This is not bookkeeping. Membrane tension is isotropic, so it is itself
three-fold symmetric — and a symmetric perturbation cannot drive an
antisymmetric motion at first order. **Only A modes can be gating
coordinates.** E modes are forbidden by symmetry, whatever their frequency.

Note the model is built on the **curved** structure alone. It never sees the
flattened one.
"""),
    ("code", """
anm = ANM.from_trimer(curved_blocks, cutoff=15.0,
                      spring="inverse_square").build()
modes = anm.calc_modes(n_modes=30)
anm.label_symmetry(modes)

print(f"{modes.n_modes} modes over {anm.n_sites * 3} degrees of freedom\\n")
for i in range(8):
    print(f"  mode {i:2d}   symmetry {modes.symmetry[i]:2s}   "
          f"frequency {modes.frequencies[i]:.5f}")
"""),
    ("markdown", """
## Compare with the observed transition

The overlap between a mode and the observed displacement is the cosine between
them: 1.0 is a perfect match, 0.0 is unrelated.
"""),
    ("code", """
overlaps = np.abs(np.asarray(modes.overlap(displacement), dtype=float))
order = np.argsort(-overlaps)

print("best modes by overlap with the observed change:")
for i in order[:5]:
    print(f"  mode {i:2d}  symmetry {modes.symmetry[i]:2s}  "
          f"overlap {overlaps[i]:.4f}")
"""),
    ("markdown", """
## The result

A single symmetric mode captures most of the transition, and every
symmetry-forbidden mode scores essentially zero. The model does not merely fit
the change; it finds it through the channel the physics permits.
"""),
    ("code", """
a_modes = [i for i in range(modes.n_modes) if modes.symmetry[i] == "A"]
e_modes = [i for i in range(modes.n_modes) if modes.symmetry[i] == "E"]

best_a = max(overlaps[i] for i in a_modes)
best_e = max(overlaps[i] for i in e_modes)
share = sum(overlaps[i] ** 2 for i in a_modes) / float((overlaps ** 2).sum())

print(f"best A-mode overlap    : {best_a:.4f}")
print(f"best E-mode overlap    : {best_e:.4f}")
print(f"share of overlap in A  : {share:.2%}")
print(f"cumulative over {modes.n_modes} modes: "
      f"{modes.cumulative_overlap(displacement)[-1]:.4f}")

assert best_a > 0.6, "the symmetric mode should capture most of the transition"
assert best_e < 0.05, "a forbidden mode should score essentially zero"
assert share > 0.99, "the overlap should sit almost entirely in A"
"""),
    ("markdown", """
## What this does *not* say

The overlap depends on the elastic-network cutoff: over 10–20 Å it ranges from
0.554 to 0.723. The qualitative result — one symmetric mode, forbidden modes at
zero — survives every cutoff. The third decimal place does not, and
`docs/SCIENCE.md` publishes the range rather than the point estimate.

An elastic network also says nothing about **energetics**. It tells you the
motion is cheap, not that tension is enough to drive it. That calculation is in
`piezo1.physics.dome` and `piezo1.physics.elastica`, and it is where the linear
Helfrich theory usually applied to PIEZO1 turns out to overestimate the
footprint energy by 3.65×.
"""),
    ("code", """
from piezo1.analysis.uncertainty import sensitivity

def best_overlap(cutoff):
    trial = ANM.from_trimer(curved_blocks, cutoff=cutoff,
                            spring="inverse_square").build()
    return float(abs(trial.calc_modes(n_modes=20).overlap(displacement)).max())

spread = sensitivity(best_overlap, [10.0, 13.0, 15.0, 18.0, 20.0],
                     knob="anm.cutoff", what="best mode overlap")
print(spread.summary())
print("\\nThis is a SENSITIVITY range over a method choice.")
print("It is explicitly not a confidence interval, and the class refuses to")
print("call it one.")
"""),
]


NOTEBOOKS = {
    "01_first_look": {
        "title": "A first look at PIEZO1",
        "cells": _FIRST_LOOK,
    },
    "02_gating_motion": {
        "title": "The gating motion, and why symmetry matters",
        "cells": _GATING,
    },
}
