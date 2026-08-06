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
