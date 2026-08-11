# Using the engine from a notebook

Nothing below `render` imports from the GUI, so the whole scientific engine
runs headlessly. This is the documented path for doing something the
application does not have a button for.

```python
from piezo1.core import Structure
st = Structure.from_file("ref/structures/8YEZ.cif")
st          # <Structure 8YEZ: 31599 atoms, 3843 residues, chains=ABC>
```

---

## The one thing to get right first

Almost every analysis on a C3 trimer needs the same preparation: equal-length
C-alpha blocks, one per protomer, in the same residue order. Get this wrong and
everything downstream returns a plausible wrong number rather than an error.

```python
import numpy as np

def protomer_blocks(st, min_ca=300):
    """Equal-length, identically ordered CA blocks, one per protomer."""
    chains = []
    for c in st.chains:
        m = st.mask_ca() & (st.chain == c)
        if m.sum() > min_ca:
            chains.append((st.xyz[m], st.res_seq[m]))
    common = set(chains[0][1].tolist())
    for _, seq in chains[1:3]:
        common &= set(seq.tolist())
    residues = np.array(sorted(common))
    blocks = [x[np.searchsorted(s, residues)].astype(float)
              for x, s in chains[:3]]
    return blocks, residues
```

**Never assume chain labels give rotational order.** Deposited entries label
protomers in either sense; 7WLT and 7WLU disagree, and superposing them by
label gives 71 Å RMSD instead of 19.7. Use
`piezo1.structure.superpose.match_protomers`.

---

## Residue numbering

Human PIEZO1 is 2521 aa, mouse 2547, and **the offset between them is not
constant** — it varies from 0 to +26 in twelve blocks and passes through zero
twice. Never add a constant.

```python
from piezo1.core.sequence import load_numbering_map
nm = load_numbering_map()
nm.to_b(2456)      # human -> mouse: 2482
nm.to_a(2496)      # mouse -> human: 2470   (NOT 2496)
```

---

## Membrane dome

```python
from piezo1.structure.geometry import measure_dome
from piezo1.structure.superpose import detect_c3_axis

blocks, residues = protomer_blocks(st)
dome = measure_dome(blocks, tm_surface_points)   # see cli.py for the TM helper
dome.radius_of_curvature / 10   # nm; ~9.7 for curved 7WLT
```

## Elastic network modes

```python
from piezo1.physics.anm import ANM

anm = ANM.from_trimer(blocks, cutoff=15.0, spring="inverse_square").build()
modes = anm.calc_modes(n_modes=30)
anm.label_symmetry(modes)        # 'A' or 'E' per mode

modes.symmetry                   # only A modes couple to isotropic tension
modes.overlap(displacement)      # cosine overlap with an observed change
modes.collectivity(0)
```

`from_trimer` is required rather than optional: it keeps protomer blocks
contiguous and identically ordered, which is what makes the symmetry labelling
valid.

## Does the network describe the molecule?

The standard check, and the one to run before quoting any mode: correlate the
predicted fluctuation against the B-factor the entry was deposited with.

```python
from piezo1.analysis.fluctuations import assess_b_factors, compare_fluctuations

assess_b_factors(st).summary()   # check the COLUMN before the network
res = compare_fluctuations(st)   # 7WLT: r = +0.434, Spearman +0.726
res.spearman_r, res.control_spearman      # 0.726 against a burial control's 0.578
res.beats_control, res.control_inverted   # True, False
res.by_mode_count                # how much the mode truncation matters
```

**Three traps, all of them live.** `assess_b_factors` refuses a uniform column,
a grouped one (3JAC has 212 distinct values over 2,754 residues) and an
AlphaFold model, whose B column holds **pLDDT** — a confidence that
anti-correlates with fluctuation at −0.57 and would otherwise look like a
strong negative result. Read `control_spearman`: contact number uses no network
at all, and if it wins, the agreement is burial. A control that comes out
**negative** means the entry's B-factor rises with burial, so the column is not
a mobility and neither number means anything.

## Pore profile

```python
from piezo1.structure.pore import pore_profile
prof = pore_profile(st, detect_c3_axis(blocks), step=1.0)
prof.bottleneck_radius, prof.is_conductive()
prof.constrictions(threshold=3.0)
```

The probe is tethered near the axis. Without that leash the clearance function
has no interior maximum and the answer runs away to ~6000 Å — a true maximum,
and useless.

## Gating kinetics

```python
from piezo1.physics.kinetics import GatingModel
model = GatingModel()                       # Young et al. 2023, four states
model.half_activation()                     # 2.71 mN/m
trace = model.step(tension=5.0, duration=0.5, n_channels=100)
trace.inactivation_tau()
model.mutant("R2456H")                      # calibrated by fold change
```

## Membrane footprint

```python
from piezo1.physics.membrane import MembraneParameters, solve_footprint
p = MembraneParameters(kappa=20.0, tension=20.0 / 14.0**2)   # lambda = 14 nm
sol = solve_footprint(r0=10.0, slope=0.35, params=p)
sol.energy, sol.excess_area(), sol.validity_note()
```

Check `validity_note()`. PIEZO1's real contact slope is ~2.0, far outside the
small-slope regime the linear theory assumes — at which point you want the
nonlinear solver instead:

```python
from piezo1.physics.elastica import solve_elastica, compare_with_linear
nl = solve_elastica(r0=8.69, slope=1.99, params=p)
nl.energy, nl.excess_area, nl.force_residual   # 25.6 kT, 177 nm^2, ~7e-11

c = compare_with_linear(r0=8.69, slope=1.99, params=p)
print(c.summary())          # linear is 3.64x too large at this slope
```

`force_residual` is the drift of the conserved axial force along the solution.
It is an error estimate that costs nothing — if it is not ~0, do not trust the
energy.

## Ensemble PCA

```python
from piezo1.analysis.ensemble import build_ensemble
ens = build_ensemble(species="mouse", min_common=900)
pca = ens.pca()
pca.variance_explained[0]        # 0.90 — PC1 is the gating coordinate
pca.best_mode_for(modes, pc=0)   # (mode index, overlap)
```

## Allostery

```python
from piezo1.analysis.allostery import (perturbation_response, cross_correlation,
                                       build_network, detour_cost, path_betweenness)
prs = perturbation_response(modes, np.tile(residues, 3))
dcc = cross_correlation(modes)
graph = build_network(np.vstack(blocks), dcc, contact_cutoff=10.0)
detour_cost(graph, blade_sites, gate_sites, anchor_sites)   # penalty ~0 for the anchor
```

Use `detour_cost`, not two separate shortest paths — independent legs can pick
endpoints in different protomers and never join, which produced a "detour"
cheaper than the direct route.

## Conservation, pockets, interactions

```python
from piezo1.analysis.conservation import load_orthologs, conservation_profile
from piezo1.analysis.pockets import find_pockets, ligand_contact_residues
from piezo1.analysis.interactions import detect_interactions

profile = conservation_profile(load_orthologs())
profile.domain_means(annotations)      # anchor is the most conserved, 0.987

find_pockets(st)                       # alpha spheres; ligands excluded by default
detect_interactions(st, min_sequence_separation=3).counts()
```

## Annotations

```python
from piezo1.core.annotations import load_annotations
ann = load_annotations("human")
ann.domain_at(2456)                    # Inner helix
ann.group("hydrophobic_gate").residues
ann.annotate_residue(2456)             # domain, sites and variants at once
```

## Reports

```python
from piezo1.analysis.report import build_report
report = build_report(st, analyses=["dome", "pore", "modes"])
report.to_markdown("report.md")
report.to_json("report.json")
```

Every report carries software version, input file, parameters and library
versions. A number without provenance is not a result.

## Hydrophobic gating

```python
from piezo1.analysis.hydration import load_grid, predict_wetting
from piezo1.structure.pore import pore_profile
from piezo1.structure.superpose import detect_c3_axis

prof = pore_profile(st, detect_c3_axis(blocks), step=1.0)
pred = predict_wetting(st, prof, load_grid())
pred.score            # 0.82 for closed 8YEZ; closed above 0.55
pred.verdict          # 'non-conductive (sterically occluded + hydrophobic gate)'
pred.dewetted[:3]     # the residues carrying the score, worst first
```

Needs `python -m piezo1.io.fetch` for the CHAP grid; without it `pred.available`
is False and the score is NaN rather than an exception.

`hydrophobic_gate` and `sterically_occluded` are **separate**. The heuristic
asks whether water would dewet, not whether water fits — 7WLU has no hydrophobic
gate but a 0.098 nm bottleneck, so it is shut for the other reason.

---

## Ion permeation

Drift-diffusion for each species along the measured pore. Gated by the wetting
verdict, because a current through a pore the model calls shut is meaningless.

```python
from piezo1.physics.permeation import (solve_pnp, default_species,
                                       blocking_mechanisms)

res = solve_pnp(prof, default_species())
res.current * 1e12       # 2.44 pA for open-like 11ZC at the default voltage
res.conductance * 1e12   # 40.7 pS
res.blocked_by           # None here; a string naming the mechanism if shut
res.converged            # the Gummel loop diverges if you push it; check this

blocking_mechanisms(pred, prof.radius, default_species())   # ALL of them, not the first
```

**The answer is dominated by two unmeasured inputs.** The in-pore diffusivity
and the effective ion radius have never been measured for PIEZO1; across their
plausible ranges the conductance spans 16–94 pS against a published 25–30 pS.
The 40.7 pS above is an overestimate and is reported as one. `series_conductance`
is the independent closed-form check on the solver.

---

## The pore's fixed charge, and selectivity

Everything above is for an electrically neutral pore, which is what
`solve_pnp` computes unless you give it a charge. Building one takes the
structure, the profile and the axis, because a residue only counts if it can
reach the lumen at its own height.

```python
from piezo1.physics.pore_charge import pore_charge, cytosolic_end
from piezo1.physics.selectivity import measure_selectivity

charge = pore_charge(st, prof, axis, mode="curated", species="mouse")
charge.summary()             # '6 charged groups (curated), net -6 e, ...'
charge.residue_summary()     # which residue numbers, how many copies, how far

res = solve_pnp(prof, fixed_charge=charge.density)
res.meta["peak_in_pore_M"]           # 13.9 — check this before believing it
res.meta["exceeds_packing_limit"]    # True: past what a solution can reach
res.meta["electroneutrality_residual"]   # ~1e-10 if the closure converged

sel = measure_selectivity(prof, fixed_charge=charge.density,
                          cytosolic_index=cytosolic_end(st, axis))
sel.summary()                # 'P_Cl-/P_Na+ = 0.021 against a published 0.14; ...'
```

**Three things will bite you here.** `cytosolic_index` is a sign, not a
formality — pass the wrong end and a cation-selective pore is reported as
anion-selective with a perfectly plausible number, which is why
`cytosolic_end` measures it. `mode="lining"` gives a *different kind* of
answer from `mode="curated"`, not a refinement of it: on 11ZC the geometric
route is net **positive** where the curated one is net negative. And an
uncharged pore is not an unselective one — it returns 0.9035 here, from size
exclusion alone, so that is the baseline any charged result has to beat.

---

## The full-length model

Experimental core plus the AlphaFold distal blade, with the join kept visible.

```python
from piezo1.structure.hybrid import build_hybrid_model

model = build_hybrid_model(st)          # st = a human PIEZO1 entry
model.seam_residue                      # 570
model.predicted.sum()                   # 4437 atoms over 569 residues
model.confident_prediction.sum() / model.predicted.sum()    # 0.48
model.overlap_rmsd, model.global_rmsd   # 2.39 A at the seam, 75.2 A overall
model.warnings()                        # what a caller must be told
```

**Use `model.experimental_only` unless you mean to include the prediction.**
Every atom carries its `source` and `plddt` precisely so that no analysis can
average across the join. The 2.4 Å seam fit says nothing about the rest of the
graft — the two models differ by 75 Å over the region they share.

---

## HaloTag fusion

There is no structure of the fusion, so this produces a **region**, not a pose.

```python
from piezo1.structure.fusion import build_fusion, load_halotag
from piezo1.structure.fusion_pose import pose_for_display, spin_scan

tag = load_halotag()                    # PDB 6U32, with its TMR ligand
model = build_fusion(st, tag)
model.pore_exit_distances()[0]          # 3.95 nm on 8YEZ
model.volume.volume                     # nm^3 the tag centre can occupy
model.meta["clashes"]                   # sphere-model verdict

pose = pose_for_display(st, model, tag) # the real fold, one orientation
pose.body_contacts                      # excludes the anchor residue - see below
pose.meta["clear_spins"]                # 1 of 36 on 8YEZ, 27 on 7WLT, 0 on 11ZC
```

Two traps. **The spin about the linker is undetermined** — `pose_for_display`
picks the least-contacting draw as a display convention, not a determination.
And when counting contacts, **exclude the anchor residue**: the placement rule
aims the tag's N-terminus at PIEZO1's C-terminus, so counting that contact
reports the rule's own construction and made all four structures look as though
no orientation fits.

Every distance here depends on `fusion.linker_residues`, which no source for
the construct states. Vary it; do not trust it.

---

## Calcium at the tag

```python
from piezo1.physics.nanodomain import Nanodomain

nd = Nanodomain(current_A=res.current, calcium_fraction=0.1,
                distance_m=3.95e-9)
nd.concentration_M * 1e6     # 225.6 uM
nd.occupancy                 # 0.999 against a 0.2 uM sensor Kd
nd.saturated                 # True whenever the channel is open
nd.falsifiers                # what would have to be observed to refute this
```

The tag distance is **modelled**, and the calcium share of the current is
unverified. `nd.sweep(...)` runs the 80-combination parameter sweep that shows
the saturation conclusion survives both.

---

## Things that will bite you

| | |
|---|---|
| Chain labels | Do not encode rotational order. Use `match_protomers`. |
| Numbering | Human↔mouse offset is not constant. Use `load_numbering_map`. |
| Coverage | Most variants of interest are not resolved. Check `modelled_in`. |
| 6KG7 | Is PIEZO2, a paralogue. Excluded from ensembles by default. |
| 3JAC | Has poly-UNK regions with arbitrary numbering; dome fits are meaningless. |
| 6LQI | Splice isoform missing 1382–1411 — a *sequence* outlier that dominates a PC. |
| Units | Coordinates Å, dome geometry nm, energies k_BT, tension mN/m. |
| Footprint | The linear solver is not quantitative at PIEZO1's 63° contact slope — 3.6× too large. Use `physics.elastica`. |
| Wetting | `predict_wetting` finds *hydrophobic* gates, not steric occlusion. Check both flags. |
| Permeation | Two inputs are unmeasured; the answer spans 16-94 pS across their ranges. Never quote it as a prediction. |
| Hybrid model | Use `experimental_only`. A 2.4 A seam fit hides a 75 A global disagreement. |
| Tag orientation | Undetermined. The drawn spin is a display convention, and contact counts must exclude the anchor residue. |
| pLDDT colours | Were applied in the wrong order until Round 76, painting every atom "very low". Use `plddt_band_colors`. |
