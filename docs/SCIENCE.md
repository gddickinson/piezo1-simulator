# The science behind the simulator

What this application models, what the numbers mean, and where each came from.
Full literature dossiers live in `ref/research/` (git-ignored; regenerable).

---

## 1. The molecule

| | Human | Mouse |
|---|---|---|
| UniProt | Q92508 | E2JF22 |
| Length | 2521 aa | 2547 aa |
| Mass | 286.8 kDa | — |
| Transmembrane helices | 38 per protomer | 38 |
| Assembly | homotrimer, C3 | homotrimer, C3 |
| Identity (global alignment) | 82.5% to mouse | — |

PIEZO1 is one of the largest polytopic membrane proteins known. Each protomer
contributes 38 transmembrane helices: TM1–36 form **nine four-helix
transmembrane helical units (THUs)** that make up a curved blade, TM37 is the
**outer helix**, and TM38 is the **inner helix** that lines the pore.

Note that NCBI's RefSeq summary still describes PIEZO1 as having 36
transmembrane segments and forming a homotetramer. Both are wrong; the cryo-EM
structures settled it.

### Domain architecture (human Q92508 numbering)

| Domain | Human | Mouse | Provenance |
|---|---|---|---|
| THU1 (TM1–4) | 1–138 | 1–138 | UniProt TM features, grouped in fours |
| THU2 | 195–330 | 202–337 | " |
| THU3 | 418–530 | 425–536 | " |
| THU4 | 575–696 | 581–702 | " |
| THU5 | 817–954 | 812–949 | " |
| THU6 | 995–1118 | 990–1113 | " |
| THU7 | 1154–1302 | 1149–1297 | " |
| Beam | 1305–1370 | 1300–1365 | literature (approximate) |
| Beam coiled coil | 1339–1370 | 1334–1365 | UniProt coiled-coil feature |
| Piezo1.1 spliced segment | 1388–1411 | 1382–1405 | PDB 6LQI |
| THU8 | 1657–1801 | 1657–1801 | UniProt TM features |
| THU9 | 1961–2076 | 1977–2092 | " |
| Anchor | 2077–2176 | 2093–2192 | derived: cytoplasmic domain between TM36 and TM37 |
| Outer helix (OH) | 2177–2197 | 2193–2213 | derived: TM37 |
| Cap / CED | 2198–2431 | 2214–2457 | derived: extracellular domain between TM37 and TM38 |
| Inner helix (IH) | 2432–2452 | 2458–2478 | derived: TM38 |
| CTD | 2453–2521 | 2479–2547 | derived: cytoplasmic domain after TM38 |

THU1 is the **pore-distal** end of the blade. The derivation was checked four
independent ways: mouse UniProt's separately curated TM coordinates land where
the alignment predicts for all 38 helices; PDB 4RAX's cap construct (mouse
2214–2457) reproduces the derived cap exactly; the 2025 human structural paper
places E756del at the THU4/THU5 linker and A1988V at the THU9 linker, both
matching; and InterPro PF24874 ("THU9 and anchor domain", 1958–2195) brackets
the derived THU9 + anchor.

**Deliberately omitted:** "clasp" and "latch". Neither survived verification —
"clasp" is not an established PIEZO1 term and appears to conflate the latch with
Guo & MacKinnon's cross-helices, and "latch" is used inconsistently with no
agreed residue range.

---

## 2. Cross-species numbering

This is the single largest source of silent error in the PIEZO literature, and
the project treats it as a first-class problem.

**The human↔mouse offset is not constant.** A global alignment (Biopython
`PairwiseAligner`, BLOSUM62, gap −11/−1; 82.5% identity, 2566 alignment columns)
gives twelve distinct offset blocks:

| Human range | Offset | Contains |
|---|---|---|
| 1–155 | 0 | THU1 |
| 156–488 | +7 | THU2, start of THU3 |
| 490–737 | +6 | end of THU3, THU4 |
| 753–777 | −6 | acidic low-complexity (unreliable) |
| 778–1375 | −5 | THU5–7, beam, coiled coil |
| 1377–1462 | −6 | Piezo1.1 spliced segment |
| 1466–1638 | −9 | disordered (unreliable) |
| 1639–1817 | 0 | THU8, Yoda1 residue A1718 |
| 1818–1910 | +6…+15 | disordered, dense indels (unreliable) |
| 1911–2393 | +16 | THU9, anchor, OH, most of cap, PIP2 cluster |
| 2394–2399 | +17 | |
| 2400–2521 | +26 | end of cap, IH, CTD, R2456, E2470 |

Consequences worth memorising:

- Mouse **A1718 = human A1718** — the same number, by coincidence. Do not
  generalise from it.
- Mouse **E2496 = human E2470**. Human E2496 is a *different residue*.
- Mouse **E2133 = human E2117**; mouse **S2472 = human S2446**; mouse
  **S1330 = human S1335**; mouse **R2482 = human R2456**.
- Mouse **K2182–K2185 = human K2166–K2169** — the same PIP2-binding lysine
  cluster. This equivalence was flagged as unverified in the literature review
  and is settled here by direct alignment, with byte-identical 29-residue
  sequence context.

Everything in the codebase converts through `piezo1.core.sequence`, never by
adding a constant.

---

## 3. The dome model of mechanotransduction

PIEZO1's blades are intrinsically curved, and they bend the surrounding bilayer
into a dome that bulges towards the cytoplasm. Because the dome's *projected*
in-plane area is smaller than its *surface* area, flattening it under tension
increases the projected area of the protein-plus-membrane unit. Membrane tension
T acting through an area change ΔA supplies the gating energy:

    ΔE = −T · ΔA

This is the essential physics: PIEZO1 is a channel gated by an area change, and
the area change comes largely from the *membrane* it deforms, not from the
protein alone.

### Measured and published geometry

| Quantity | Value | Source |
|---|---|---|
| Radius of curvature, closed | 10.2 nm | Haselwandter & MacKinnon 2018, *eLife* |
| Radius of curvature, outside-in | ~11.8 nm | Vaisey & MacKinnon 2026, *Sci Adv* |
| Mid-bilayer dome area | ≈390 nm² | Haselwandter & MacKinnon 2018 |
| Dome is protein/lipid | ~20% / ~80% | " |
| Footprint decay length λ | 14 nm | " |
| Membrane bending modulus κ | 20–25 k_BT | Haselwandter 2018; Dixit 2025 |
| Intrinsic radius R₀ | 42 ± 12 nm | Haselwandter 2022, *PNAS* |
| Protein stiffness K_P | 18 ± 2.1 k_BT | " |
| Flattening force constant | 60 ± 20 pN/nm | Dixit, Noé & Weikl 2025, *eLife* |
| Dome depth | 6–7 nm | Chong et al. 2021 |
| CED displacement on flattening | ~12 Å | Vaisey & MacKinnon 2026 |

**The footprint matters more than the dome.** Haselwandter & MacKinnon's key
result is that the membrane deformation extending *beyond* the protein, decaying
over λ ≈ 14 nm, dominates the tension sensitivity — not the dome itself.

### The footprint, solved

The linearised Helfrich energy in the Monge gauge,

    E = ∫ [ (κ/2)(∇²h)² + (γ/2)|∇h|² ] dA

minimises to κ∇⁴h − γ∇²h = 0, whose axisymmetric decaying solution is
`h(r) = A·K₀(r/λ)` with **λ = √(κ/γ)**. For a fixed contact slope `s` at the
inclusion radius `r₀` the energy has the closed form

    E = π κ s² (r₀/λ) K₀(r₀/λ) / K₁(r₀/λ)

Note the direction of that Bessel ratio — inverting it gives an answer 2.5×
too large at the r₀/λ where PIEZO1 sits.

Our solver reproduces this to second order and recovers λ = 13.998 nm from its
own output against an input of 14.0 nm. Applied to the measured 7WLT dome it
stores **622 nm² of excess area against the dome's own 256 nm²**, i.e. the
surrounding membrane holds about 2.4× as much deformable area as the dome —
the quantitative form of Haselwandter & MacKinnon's claim that the footprint,
not the dome, dominates tension sensitivity.

**Validity caveat, stated because it matters.** PIEZO1's dome meets the bilayer
at a contact slope near 2.0, about 63°. The small-slope expansion behind the
Monge gauge drops terms of order |∇h|², so at that slope the neglected terms
are larger than the ones kept. The solver returns numbers and flags them as
indicative of scale and trend only; quantitative work needs a full nonlinear
Helfrich or Euler–elastica treatment.

### An unresolved discrepancy, stated plainly

Reported area changes for gating differ by more than an order of magnitude
depending on how they are obtained:

- **Functional** ΔA from Boltzmann fits to tension–response curves: **6–20 nm²**
  (e.g. 8 ± 1 nm², Cox et al. 2016).
- **Structural/simulated** nanodome excess area: **~40 nm²** at zero tension
  (Dixit 2025), ~120 nm² projected (Guo & MacKinnon 2017), ~300 nm² in-plane
  expansion (Yang et al. 2022), up to 500 nm² by some measures.

These are measuring *different quantities* and should not be quoted
interchangeably. The dome model makes the consequence explicit by computing the
half-activation tension each would imply, holding ΔG₀ = 9.7 k_BT fixed:

| ΔA source | ΔA (nm²) | Kind | Implied T₅₀ (mN/m) |
|---|---|---|---|
| Cox 2016 | 8 | functional | **4.99** |
| Dixit 2025 | 40 | structural | 1.00 |
| Guo & MacKinnon 2017 | 120 | structural | 0.33 |
| Yang 2022 | 300 | structural | 0.13 |

Measured T₅₀ is 2.7–5.1 mN/m. The functional area reproduces it — which is a
consistency check, since all three of Cox's numbers come from the same fit —
while the structural areas under-predict the threshold by 5–40×. They are not
the gating area.

### What this code measures

Taking the mid-point of each of the 38 transmembrane helices in each protomer as
a sample of the mid-membrane surface, recovering the three-fold axis, and
fitting a sphere:

| Structure | State | R_c | Dome depth | Excess area |
|---|---|---|---|---|
| 7WLT | curved, bilayer | **9.7 nm** | 4.9 nm | 256 nm² |
| 11YE | curved, native vesicle | **10.4 nm** | 4.6 nm | 293 nm² |
| 8YEZ | human apo | 12.0 nm | 5.8 nm | 279 nm² |
| 8ZU3 | human + MDFIC | 12.5 nm | 5.4 nm | 270 nm² |
| 7WLU | flattened | 18.4 nm | 2.5 nm | 379 nm² |
| 11ZC | flat, native vesicle | 21.6 nm | 3.5 nm | 335 nm² |

The curved values sit squarely on the published 10.2 nm. This is the standing
regression test for the geometry pipeline.

*Caveat:* curved and flat entries resolve different residue ranges, so their
footprint radii are not directly comparable without restricting to commonly
resolved residues.

---

## 4. Elastic network models and the symmetry argument

An anisotropic network model treats each C-alpha as a bead joined to its
neighbours by harmonic springs. Its low-frequency normal modes describe the
large collective motions that all-atom simulation cannot reach for a
2500-residue trimer.

### Why symmetry labelling is not decoration

PIEZO1 is a C3 homotrimer, so every normal mode carries an irreducible
representation label: **A** (totally symmetric) or **E** (a degenerate pair).
Isotropic membrane tension is itself C3-symmetric. By the standard
group-theoretic selection rule, **only A modes can couple to it at first
order.** Identifying the lowest A mode therefore identifies the candidate gating
coordinate, and E modes can be excluded on principle rather than by inspection.

The application computes the character ⟨u, Su⟩ where S rotates by 120° and
permutes protomers; it comes out at exactly +1.000 for A modes and −0.500 for E
modes on real structures.

### The validation

Curved 7WLT versus flattened 7WLU, 1274 residues common to all six protomers,
19.7 Å trimer RMSD after superposition:

| | |
|---|---|
| Best single mode overlap | **0.705** (mode 3, symmetry A, collectivity 0.610) |
| Cumulative overlap over 40 modes | **0.964** |
| Best E-mode overlap | **0.0011** |
| Fraction of overlap² in A modes | **100.00%** |

An overlap of 0.7 for a single mode against a 19.7 Å conformational change is
high by the standards of the ENM literature. The complete absence of E-mode
overlap is a strong internal consistency check: the observed transition is
C3-symmetric, and the analysis recovers that without being told.

### The experimental ensemble agrees

The single-transition test above uses one pair of structures. A stronger test
asks whether the elastic network predicts the direction the *whole deposited
record* varies in. Principal component analysis over 10 mouse PIEZO1
structures, placed on a shared 1091-residue basis in human numbering:

| | |
|---|---|
| Variance in PC1 | **90.0%** |
| PC1 overlap with ANM mode 6 (symmetry A) | **0.804** |
| Cumulative overlap over 30 modes | **0.960** |
| RWSIP, 8 components vs 8 modes | **0.555** |
| Random-vector control | 0.001 |

PC1 orders the structures by gating state without ever being told what those
states are — seven curved entries negative, the intermediate at +334,
flattened at +678, flat at +1045. And the top three principal components all
match **A**-symmetric modes, even though E modes outnumber them two to one.

Two exclusions were necessary and both are instructive. **6KG7 is PIEZO2**, a
40%-identity paralogue; it is a different protein and does not belong in this
ensemble. **6LQI is the Piezo1.1 splice isoform** missing residues 1382–1405,
so what distinguishes it is sequence, not conformation — included, it dominates
a component on its own and splits the gating coordinate across PC1 (58%) and
PC2 (36%).

### A trap this exposed

Deposited chain labels do **not** reliably indicate rotational order around the
symmetry axis. 7WLT and 7WLU label their protomers in opposite senses.
Superposing by chain label gave 71.2 Å RMSD instead of 19.7 Å, and any
difference vector built that way is meaningless. Protomer correspondence is now
always determined by superposition
(`piezo1.structure.superpose.match_protomers`).

---

## 4b. The pore, measured from coordinates

The pore-radius profile is the largest sphere that fits at each height along
the three-fold axis without overlapping any van der Waals surface, with the
probe tethered within 8 Å of the axis.

**The leash is not a convenience, it is a correctness requirement.** The
clearance function has no interior maximum — a free probe leaves the pore
sideways and finds bulk solvent, growing without bound. Unconstrained on real
PIEZO1 coordinates it reaches R ≈ 6188 Å, which is a true maximum and a
completely useless answer.

Measured on the closed human structure 8YEZ:

| Feature | z (Å) | Radius (Å) | Lining |
|---|---|---|---|
| Transmembrane hydrophobic gate | −17.7 | **3.01** | 2449, 2450, 2451 |
| CTD constriction 1 | −46.2 | **1.24** | M2467 |
| CTD constriction 2 | −55.2 | **0.76** | 2509, 2514 |

The global bottleneck is 0.76 Å, so the structure is not conductive — correct
for a closed channel. The flat, open-like 11ZC gives 3.25 Å and is conductive.

These constrictions were found from coordinates alone. That they coincide with
the residues independently curated from the literature as the hydrophobic gate
(I2447/V2450/F2454) and the CTD constrictions (M2467/F2468, P2510/E2511) is a
mutual validation of the profiler and the annotation.

Note that HOLE cannot be used here: it has no Apple-Silicon build, and
MDAnalysis's `hole2` wrapper is an empty stub as of 2.10.

## 5. Electrophysiology

| Quantity | Value | Notes |
|---|---|---|
| Unitary conductance, with divalents | 29.1 ± 0.4 pS | Coste et al. 2015 |
| Unitary conductance, divalent-free | 58.6 ± 1.2 pS | same channel |
| Native endothelial | ~25 pS | Shi et al. 2020 |
| In GPMVs | ~27.5 pS | Vaisey & MacKinnon 2026 |
| Permeability P_K:P_Cs:P_Na:P_Li | 1 : 0.88 : 0.82 : 0.71 | |
| P_Ca/P_Cs | 1.21 | |
| P_Cl/P_Na | 0.14 | |
| Inactivation τ at −80 mV | ~34 ms | Romero et al. 2019 |
| Half-activation tension T50 | 2.7 ± 0.1 mN/m (cell-attached) | Lewis & Grandl 2015 |
| | 4.7 ± 0.3 mN/m (inside-out) | " |
| | 5.1 ± 0.2 mN/m | Cox et al. 2016 |
| Gating free energy ΔG₀ | 9.7 ± 1.5 k_BT | Cox et al. 2016 |
| P50, cell-attached | 36.4 ± 3 mmHg | Ridone et al. 2020 |

The 23–130 pS spread reported across the literature resolves once ion identity,
divalent content and current direction are accounted for — it is not a
disagreement. **P_Ca/P_Na has never been directly measured.**

Rectification is a *gating* effect, not a pore effect: the instantaneous
rectification index is 1.13 ± 0.06 versus 5.3 ± 0.6 at peak.

**Clustering does not alter gating.** P50 and open probability are invariant
from 1 to 100 channels/µm² (Lewis & Grandl 2021), so cooperativity should not be
modelled.

---

## 6. Lipids are not optional

The most consequential recent result: **mechanical force is necessary but not
sufficient**. PIEZO1 reconstituted into synthetic lipids (soy PC, or
POPC/DOPS/cholesterol 8:1:1) fails to activate mechanically, while
HEK293-derived plasma-membrane vesicles support activation (Vaisey & MacKinnon,
*Sci Adv* 2026). Cryo-EM shows an unidentified lipid-like density with a
headgroup and two branched acyl chains entering laterally from the inner
leaflet, engaging the conserved lysine cluster — a phosphoinositide or ceramide.

Other quantified lipid effects:

- **PIP2/PI(4)P**: both must be depleted to inhibit; depleting PI(4,5)P₂ alone
  is not significant (Borbiro 2015). Binding site = human K2166–K2169. But
  exogenous PIP2 *suppresses* open probability in brain capillary endothelium
  (Hashad 2025, *PNAS*) — **an unresolved directional conflict**.
- **Cholesterol**: depletion right-shifts P50 from 36.4 ± 3 to 54.2 ± 1.1 mmHg
  and doubles the diffusion coefficient (Ridone 2020).
- **Fatty acids**: margaric acid inhibits, IC50 28.3 ± 3.4 µM. EPA speeds
  inactivation (34 → ~20 ms), DHA slows it (34 → ~53 ms) — opposite directions.
  EPA normalises the R2456H gain-of-function phenotype to wild type.
- **Sphingomyelinase**: ceramide production makes native endothelial PIEZO1
  non-inactivating (Shi 2020).

---

## 7. Pharmacology

| Compound | Type | Potency | Site |
|---|---|---|---|
| Yoda1 | agonist | EC50 17.1 µM (mouse), 26.6 µM (human) | pocket at human A1718/A2075/A2078 |
| Yoda2 | agonist | EC50 0.15 µM | — |
| Jedi1 / Jedi2 | agonist | EC50 ~200 / 158 µM | extracellular distal blade |
| Dooku1 | antagonist | IC50 1.3–1.5 µM | competitive Yoda1 congener, zero efficacy |
| GsMTx4 | inhibitor | K_D 155 nM | amphipathic, acts through the bilayer |

**No agonist co-structure exists.** Every PIEZO PDB entry contains only lipids.
The Yoda1 site is mapped by mutagenesis and simulation, and the application
labels it as *predicted* rather than experimental. Notably, a lipid (PLX)
occupies that pocket in 7WLT — so Yoda1 may act by displacing a lipid.

That GsMTx4's **D-enantiomer is equally effective** is the key evidence that it
acts through the membrane rather than at a stereospecific protein site.

---

## 8. Variants and disease

68 curated variants ship with the application: 22 gain-of-function, 17
loss-of-function, 8 VUS, 6 blood-group, 15 engineered. Every wild-type residue
was verified against Q92508.

- **Hereditary xerocytosis (DHS1)** — gain of function, slowed inactivation.
  Archetype **R2456H**: τ 22.2 ± 2.1 ms versus 8.6 ± 0.4 ms wild type (2.6×).
- **Lymphatic dysplasia (LMPHM6)** — recessive loss of function. Constraint is
  consistent: pLI ≈ 0, LOEUF 1.097.
- **E756del** — a polyglutamate microsatellite contraction, not a point
  deletion, so the "756" assignment is arbitrary within the tract. gnomAD AFR
  allele frequency 0.166–0.173. **Malaria protection is contested**: Thye 2022
  found heterozygote OR 0.91, p = 0.19, and the original mouse work tested
  R2482H, not E756del.
- **Er blood group** — PIEZO1 carries the Er antigens; G2394 is Erᵃ/Erᵇ. All Er
  variants give wild-type Yoda1 currents: antigenic, not channelopathic.

**Coverage caveat, built into the UI.** All six human PIEZO1 structures model
from residue 570 only. E756 is not modelled in 9VMX; A1988 is not modelled in
8ZU8 or 8YFC. Of the 68 variants, **14 are resolved in no human structure at
all**, and only R2456 appears in its own structure (8YFG).

---

## 9. Known gaps

Stated so nobody has to rediscover them:

- No measured P_Ca/P_Na.
- No measured PIP2–PIEZO1 binding affinity.
- No experimental test of PIEZO1 in bilayers of systematically varied thickness.
- No agonist-bound structure.
- Functional and structural ΔA differ by 10–50×, unexplained.
- The PIP2 direction-of-effect conflict between Borbiro 2015 and Hashad 2025.
- Residues 1–570 (human) are unresolved in every experimental structure.
  AlphaFold covers them but its PAE between that block and the rest is 25–29 Å
  against a 31.75 Å maximum — i.e. **AlphaFold cannot place the distal blade
  relative to the core**. PIEZO2 (6KG7, resolving residues 8–823) is the best
  experimental guide.

---

## Key references

Full bibliographies with PMIDs are in `ref/research/`.

- Haselwandter & MacKinnon, *eLife* 2018 — membrane footprint. PMID 30480546
- Guo & MacKinnon, *eLife* 2017 — dome model. PMID 29231809
- Zhao et al., *Nature* 2018 — lever-like transduction. PMID 30089899
- Yang et al., *Nature* 2022 — curved and flattened states in bilayer
- Vaisey & MacKinnon, *Sci Adv* 2026 — lipid cofactor requirement. PMID 42234740
- Dixit, Noé & Weikl, *eLife* 2025 — nanodome elasticity
- Lewis & Grandl, *eLife* 2015 — tension sensitivity
- Coste et al., *Nature* 2015 — pore properties. PMID 26649819
- Buyan et al., *Biophys J* 2020 — lipid interaction sites. PMID 32949489
- Botello-Smith et al., *Nat Commun* 2019 — Yoda1 mechanism. PMID 31582801
- Syeda et al., *eLife* 2015 — Yoda1. PMID 26001275
- Bae, Sachs & Gottlieb, *Biochemistry* 2011 — GsMTx4. PMID 21696149
