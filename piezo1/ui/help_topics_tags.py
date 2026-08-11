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
