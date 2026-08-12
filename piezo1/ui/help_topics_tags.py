"""The help topic for the tag, the labelling and the ion current.

Split out of ``help_content.py`` at the project's length limit, along the same
seam ``report_tags.py`` and ``parameter_table_tags.py`` already use: this
describes the tag and what flows through the pore, not the channel's mechanics.
It is the longest topic because it is the one carrying the most caveats — a
modelled tag position, two unmeasured confinement parameters, and a selectivity
result whose direction is right and whose value is not.
"""

from __future__ import annotations

__all__ = ["TAGS_AND_CURRENT"]

TAGS_AND_CURRENT = """
<h2>HaloTag fusion, labelling and ion current</h2>

<h3>Seeing the tags — View &rarr; HaloTag fusion</h3>
<p>PIEZO1 imaging constructs fuse <b>HaloTag</b> to the cytosolic C-terminus,
one per protomer. <b>There is no structure of the fusion</b>, so everything
drawn here is a model: the linker is a straight seam rather than a
conformation, and the <b>accessible volume</b> is shown as a point cloud
precisely so a single position is not mistaken for a determined one.</p>
<p>The tag body can be drawn two ways. <b>Show modelled tags</b> gives a sphere
of the tag's radius of gyration — the shape that claims exactly what the model
determined, a position and nothing about orientation. <b>Show tag structure</b>
gives the real 6U32 fold, rigidly placed at the same centre and turned so its
N-terminus faces the channel's C-terminus. That is more informative and more
dangerous, because a drawn fold reads as a determined pose. It is not: the
<b>spin about the linker is undetermined</b>, and <b>Turn tag orientation</b>
rotates it 10° at a time so you can see that for yourself — the fold turns and
nothing else moves. Atoms inside the channel are drawn red, and the status line
says how many of the 36 sampled orientations clear it.</p>
<p>Colour is <i>not</i> what keeps a modelled tag from reading as experimental
structure: the tag's orange sits 0.10 from the chain palette's orange, and every
colour genuinely distant from the eight chain hues is too dark to see. The
status line is the guard, which is why the fold cannot be drawn without it.</p>
<p>Drawing the fold measures something the sphere could only assert. Over 36
spins the fold clears the channel in <b>27 of 36</b> orientations on 7WLT, 7 on
8YFG, <b>1</b> on 8YEZ and <b>none</b> on 11ZC — and 11ZC is exactly the entry
whose sphere clearance (15.7 Å) falls below the radius of gyration (17.6 Å). The
two models agree on which structures admit a tag, while the sphere is generous
about how much room there is.</p>
<p>Measured inputs, from PDB <b>6U32</b> (1.8 Å, TMR ligand bound): radius of
gyration <b>17.6 Å</b>, N-terminus <b>19.9 Å</b> from the centre, furthest atom
<b>30.0 Å</b>. A C-terminal fusion attaches to the tag's <i>N</i>-terminus, so
that offset — not the radius of gyration — sets where the body sits.</p>
<p>The tag centre lands <b>3.3–4.2 nm</b> from the pore exit across all twenty
downloaded entries. That is <i>below</i> the 4–6 nm this project first
estimated, and sweeping the unverified linker length from 1 to 30 residues moves
it by under a nanometre, so the miss is structural rather than an artefact of the
assumption. About half the accessible volume does lie in the 4–6 nm band: the
window describes a reachable position, not the ensemble mean.</p>

<h3>Labelling — Analysis &rarr; HaloTag labelling</h3>
<p>The kinetics are <b>imported</b> from the companion
<tt>halotag_binding_sim</tt> project and reproduced to machine precision, so any
divergence would mean the import is wrong rather than that anything was
discovered.</p>
<p>Because all three sites must bind for a channel to be fully labelled, a
per-site shortfall is <b>cubed</b>: p = 0.9 leaves only 0.73 of channels fully
labelled. At the standard protocol (200 nM, 30 min) labelling is complete in
<b>54 s</b>, giving a 100% three-dye population — so at any realistic
concentration there is <i>no</i> kinetic dye mixture. A mixture at a saturating
protocol argues instead for <b>chemically unreactive tags</b>, whose ceiling is
the reactive fraction cubed and which no incubation removes.</p>

<h3>Ion current — Analysis &rarr; Ion permeation</h3>
<p>Drift-diffusion for each ion species along the measured pore, with the
spreading resistance of the pore mouths in series, and <b>gated by the wetting
verdict</b>: a pore wide enough for an ion still carries no current if it has
dewetted, because the hydration shell the ion must shed into is not there.</p>
<p>When the pore is shut the report lists <b>every</b> mechanism, not the first.
That distinction is real: <b>8YEZ</b> is shut both sterically <i>and</i> by a
hydrophobic gate, while <b>7WLU</b> is shut only sterically.</p>
<p>For the open <b>11ZC</b> the model gives <b>41 pS</b> against a published
<b>25–30 pS</b>. Two things temper that. The in-pore diffusivity and the ion
radius are <b>unmeasured</b>, and across their plausible ranges the answer spans
<b>16–94 pS</b> — so the model can be made to agree, but that would be tuning.
And the Debye length in 150 mM (5.7–8.1 Å) <i>exceeds</i> the open pore radius
(3.3 Å), so the double layers overlap and a continuum treatment of this pore is
at the edge of its validity. The potential is therefore solved in the
electroneutral limit, which agrees with an independent closed-form check to
1.5%.</p>
<p><b>Selectivity, added in Round 81.</b> Every current above was computed for
an electrically <i>neutral</i> pore: the solver had taken a fixed-charge
argument since it was written and nothing had ever supplied one. It now does,
from the ionisable residues that can reach the lumen, and the report runs the
published protocol — 150/30 mM NaCl, reversal potential inverted through GHK —
against Coste et al.'s measured <b>P_Cl/P_Na = 0.14</b>.</p>
<p>The charge makes the model <b>cation-selective</b>, which is the direction
the measurement has. It does not reproduce the value: the curated pore residues
give <b>0.021</b> and every group that reaches the lumen gives <b>0.207</b>,
bracketing 0.14 rather than landing on it. Three things are reported beside
that and none should be skipped. The <i>uncharged</i> pore already sits at
0.9035, because chloride is nearly twice sodium's radius and loses more
cross-section at a 3.3 Å bottleneck than it gains in mobility — so part of the
preference for cations is size, not charge. The curated route only reaches
0.021 at an in-pore concentration of <b>13.9 M</b>, past anything a solution
could hold; the result carries that flag rather than being clipped. And the two
routes are net negative and net <i>positive</i> respectively, so they disagree
about what kind of thing the pore is, not merely by how much.</p>
<p>One measurement fell out on the way: <b>three of the four glutamates the
annotation calls selectivity determinants are not within side-chain reach of
the lumen</b> on the open structure — E2117 sits 12.9 Å past the wall. That
agrees with the paper that identified it, which concluded from function alone
that the residue probably modulates the pore rather than lining it.</p>

<h3>Why almost nothing animates — View &rarr; Ion flux animation</h3>
<p>The stream is gated by the same wetting verdict, so <b>17 of the 19
deposited PIEZO1 entries show no ions at all</b>. That is the honest outcome
and not a fault — but the reason matters, and until Round 84c the status line
did not give it.</p>
<p><b>The refusal is almost never about the gate.</b> Locating the
transmembrane gate from the curated <tt>hydrophobic_gate</tt> residues and
asking where the profile is actually pinched: in <b>none</b> of the 18 entries
whose gate can be located is the narrowest point <i>at</i> the gate. It is
below it, at the cytoplasmic constriction, in 16, and above it in the cap in
2. The gate itself measures <b>2.4–4.7 Å</b> everywhere — at or above the
1.5 Å water radius the steric test uses. The status line now names the
constriction and quotes the gate's radius beside it.</p>
<p><b>8IXO is the case that makes the point.</b> It is Liu et al. 2025's
intermediate-<i>open</i> S2472E structure. Its gate is 3.52 Å; its lining
clears the hydrophobicity cutoff at 0.31; and the V2476 side-chain diagonal
measures <b>14.2 Å</b> against <b>7.7 Å</b> on the curved 7WLT, reproducing the
7 &rarr; 14 Å dilation that paper reports. It is still refused — on a
<b>0.98 Å</b> neck at E2537, the vertical constriction the same paper says
<i>remains closed</i>, because <b>the lateral portals carry the current</b>.
Our conduction model is one-dimensional and axial: it has no lateral portals in
it, so it must pass through a constriction the real channel goes around. That
is a stated limit of the model, not a property of the structures.</p>
<p><b>And the two entries that do animate are the two weakest models.</b>
<b>11ZC</b> is deposited at 6.0 Å with backbone atoms only — no side chains
anywhere — which is why nothing in it is narrow, why its gate measures 7.04 Å
and why it runs at 2.4 pA. <b>3JAC</b> is 4.8 Å with 346 unnamed residues and
runs at 0.18 pA, a tenth of the measured conductance rather than the 1.5×
overestimate 11ZC gives. Both numbers carry that comparison on the HUD. The
animation is a picture of a rate, and on these two entries the rate is a
picture of the resolution.</p>

<h3>Choosing the route — View &rarr; Ion flux pathway</h3>
<p>Round 84d made the conduction pathway a <b>choice</b> instead of an
assumption, because the assumption was wrong for this channel.</p>
<p><b>Axial</b> is the default and is what every number this project has
recorded was computed on: bulk solvent, down the three-fold axis, bulk solvent.
It refuses all but two entries, and the reason is that PIEZO1's axis is closed
at <i>both</i> ends on purpose — pinched to about 1 Å at <b>R2295</b> on top
and at the cytoplasmic constriction beneath. Liu et al. 2025 report both:
Na⁺ enters the cap vestibule through <b>three lateral cap gates</b> and leaves
the inner vestibule through <b>intracellular lateral portals</b>, and neither
closed end is ever traversed.</p>
<p>The <b>lateral</b> options exclude those ends, which is the smallest change
that lets the model represent that route. On it, <b>8IXO conducts at 53.8 pS</b>
where the axial model refuses it. Two things must be read with it. The portal
is <b>not modelled</b> — the truncated end slice becomes the mouth and its
radius is the pore's, not the portal's — so the current is an <b>upper
bound</b>. And it does <b>not</b> separate open from closed: 7WLT, 6B3R, 8IMZ
and others also conduct at 6–12 pS once the ends are open, against 8IXO's 53.8.
The right ordering, roughly fivefold, and not the clean contrast their
simulations show. Opening the ends is necessary and not sufficient.</p>
<p><b>View &rarr; Ion flux voltage</b> sweeps the four transmembrane potentials
of their Figure 5A. At 0 V the current is zero and nothing animates, which is
the correct picture rather than a failure. Across the four, the slope on 8IXO
is <b>40.1 pS against their 20 pS</b> — twice, which is the same overestimate
this solver already carries against the measured unitary conductance.</p>
<p>What the animation still is not: their Figure 5 has explicit ions with
trajectories, and this has a concentration field. The particles are a rate made
visible. Their Figure 5C — how many Na⁺ reached each cavity — has <b>no
analogue here at all</b> and the panel registry says so: a one-dimensional
steady state carries one flux through every slice by construction, and their
count includes ions that entered and turned back.</p>

<h3>What the variant structures can support — Analysis &rarr; Variant structures</h3>
<p>A null result, shipped rather than worked around. Of the four deposited
variant entries, <b>one</b> (8YFG, R2456H) actually resolves its own mutation;
A1988 and E756 are unmodelled in the entries named for them. <b>8ZU3, 8YFC and
9VMX share one set of coordinates</b> — separate depositions, different files,
identical atoms — so they cannot distinguish anything from each other. And every
deposited human structure is <b>closed</b>, so no difference in conductance can
be measured between them.</p>
<p>All four are <b>gain-of-function</b>. There is no deposited loss-of-function
structure, so this route cannot discriminate direction even in principle. That
is the same data limit the blind test met from the other side.</p>
"""
