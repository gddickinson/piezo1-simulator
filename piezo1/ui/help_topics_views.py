"""The help topic for everything drawn *on top of* the structure.

Split out of `help_topics_physics.py` at the length limit and along a real
seam. Every other topic explains part of the model; this one explains the
overlays — seven things that put a measurement on screen next to the protein it
was measured from.

They share a failure mode, which is why they share a topic: a drawn surface is
more persuasive than the number it came from, and each one of these can be read
as claiming something it does not. A sphere fitted to the wrong atoms still
looks like a dome. A probe sphere looks like a wall. A pile of alpha spheres
looks like a volume you could add up. A single line looks like *the* pathway.
So each section ends on what the picture must not be read as, and every one of
those sentences is also on the status line, where it cannot be skipped.
"""

from __future__ import annotations

__all__ = ["VIEWS"]

VIEWS = """
<h2>Drawing what was measured</h2>
<p>Everything under the <b>View</b> menu below the display options puts a
computed result on screen beside the structure it came from. None of them
computes anything of its own: each calls the same function the corresponding
panel or table calls, so if a picture and a panel ever disagree, the picture is
wrong and that is a defect rather than a finding.</p>

<h3>The dome — View &rarr; Dome surface</h3>
<p>The dome is the project's central geometric claim and it was four numbers in
a status bar. This draws it. The <b>blue cap</b> is the sphere fitted to the
transmembrane helices on screen, out to the footprint radius — and a picture is
the only check that catches a sphere fitted to the <i>wrong</i> atoms, which
returns a perfectly reasonable radius. The <b>grey disc</b> is that cap's own
flat projection, so the gap between the two surfaces <i>is</i> the excess area
the gating model is built on.</p>
<p><b>What is not drawn, and why.</b> The obvious third surface is the bilayer
relaxing back to flat outside the rim. It is not here. PIEZO1's cap meets the
membrane at a slope near <b>1.9</b> — a 63° contact angle — and the linearised
Helfrich theory is a small-slope expansion: continued from that rim it plunges
158 Å over a 526 Å skirt and overestimates the real footprint <b>3.65×</b>. The
footprint radius is given as a number instead.</p>

<h3>The contacts — View &rarr; Contacts</h3>
<p>The interaction inventory has existed since Round 21 and only ever as a
table. A table can tell you that the R2456–E2117 salt bridge exists; it cannot
tell you that it joins <i>two different protomers</i> at the pore, which is the
part that matters. This draws each contact as a cylinder between the two atoms
the analysis actually found it between, coloured by kind — gold disulfides,
blue salt bridges, green hydrogen bonds, purple π-stacks, pink cation–π, grey
hydrophobic.</p>
<p><b>Two kinds are off until you ask for them, and not because there are a lot
of them.</b> 8YEZ yields 9,863 contacts, of which 7,984 are hydrogen bonds —
and most of those are backbone <i>i</i> to <i>i</i>+4, which <i>are</i> the
secondary structure the cartoon is already drawing. Turning them on repeats the
ribbon in eight thousand green lines. Hydrophobic contacts are off because two
carbons near each other is the weakest claim in the set. What is left is 270
cylinders, each saying something the ribbon does not. Whatever is hidden is
counted on the status line rather than dropped.</p>
<p><b>Two things the picture must not be read as.</b> No deposited PIEZO entry
has hydrogens, so every criterion here is heavy-atom geometry — a drawn hydrogen
bond is an inference, not an observed proton. And a contact belongs to
<i>this</i> structure in <i>this</i> state: a closed entry does not show the
open state's salt bridges, and a residue whose side chain was never resolved
cannot contribute one at all.</p>

<h3>The pore — View &rarr; Pore surface</h3>
<p>The pore radius has been a two-axis plot since Round 12, and a plot answers
<i>how narrow</i> without ever answering <i>where</i>. This draws the profile as
the probe spheres it was measured with: at each height, the largest sphere that
fits without overlapping an atom. Same profile the Analysis panel plots — it is
read, not recomputed, so the picture and the plot cannot be of different runs.
</p>
<p>Three colours, from two registered numbers. <b>Red</b> is narrower than
<tt>pore.ion_radius</tt>, so a bare ion of that size does not fit.
<b>Amber</b> clears that but not <tt>pore.constriction_threshold</tt>, the
conventional hydrated-ion cut. <b>Blue</b> clears both. The pale spheres are the
C-alphas of the residues lining the narrowest slice, in all three protomers —
the bottleneck marked without misstating any radius.</p>
<p><b>These are the spheres that fit, not the pore wall.</b> They are the space
left over; the protein surface is somewhere outside them at an unstated
distance. And a radius does not settle whether the channel conducts: a lumen
wide enough for an ion still dewets if it is hydrophobic, which is the whole
content of Round 19. The wetting verdict is on the status line, and it — not the
colour — is this project's answer to "is it open?".</p>

<h3>The pockets — View &rarr; Pockets</h3>
<p>A pocket was a list of residue numbers, which says what lines a cavity and
nothing about its shape — and shape is why anyone looks. The top five ranked
pockets are drawn as the alpha spheres the detector found them with, one colour
each. Lower-ranked ones are counted on the status line rather than dropped.</p>
<p><b>The spheres overlap, heavily.</b> Counting them is counting the sampling,
not the cavity: summing 4/3πr³ overcounts several fold, which is why
<tt>Pocket.volume</tt> integrates a Monte-Carlo union instead. <b>And a cavity
is not a binding site.</b> No deposited PIEZO entry contains a bound modulator,
and the one residue-level site this project holds is inferred from docking
rather than observed. Detection also excludes ligands first — a resolved lipid
fills the very pocket being looked for — so a drawn pocket may sit exactly on
top of one that is also on screen.</p>

<h3>The blade-to-gate route — View &rarr; Allosteric path</h3>
<p>This is the picture of the project's central mechanical claim: force on the
blade reaching the pore. The tube is the cheapest route through a contact graph
whose edges cost <tt>−log|DCC|</tt>, so a pair that moves together is cheap to
cross. It needs normal modes first, from the Physics panel. The tube is coloured
by each step's own correlation, so the weakest link in the chain is visible
rather than averaged away.</p>
<p>The search starts at <b>the most distal blade unit the entry actually
resolves</b> — THU4 on a deposited structure, THU1 on a full-length model.
Handing it the whole blade makes the answer trivial: the cheapest route then
starts at whichever blade residue happens to sit nearest the pore and never goes
near the beam.</p>
<p><b>A drawn line reads as unique, and it is not.</b> Rather than say so and
leave it, the same search is re-run with this route's own steps deleted from the
graph, and the status line reports what the best remaining route costs. On 8YEZ
it costs <b>1.001×</b> — the line on screen is a representative of a
near-continuum of equivalent routes, not the pathway. The correlations are the
elastic network's own, so none of this is a measured signal; it is a statement
about the structure.</p>

<h3>The calcium — View &rarr; Calcium nanodomain</h3>
<p>Round 32 concluded the calcium sensor on a HaloTag is saturated whenever its
own channel opens, so puncta brightness reports labelling stoichiometry and open
probability rather than calcium amplitude. That is a claim about a distance, and
this draws it: translucent shells at decade concentrations around the measured
cytosolic mouth of the pore. The Green's function is spherically symmetric, so
these are exactly spheres — nothing is idealised for drawing.</p>
<p><b>The two surfaces that carry the conclusion are not drawn.</b> At 11ZC's
2.43 pA the sensor is still 90% occupied out to <b>119 nm</b> and half-occupied
at its Kd out to <b>372 nm</b>, against a channel reaching about 15 nm. Drawn,
they are two enormous shells with a speck inside, and the speck is the protein —
so they are reported as numbers, the same rule the far-field membrane footprint
follows. The whole channel sits inside both, which <i>is</i> the result.</p>
<p><b>A shut structure draws nothing.</b> The analysis report borrows 11ZC's
current when the loaded entry is closed and labels the substitution; a picture
cannot carry that label convincingly, because a cloud drawn around 8YEZ reads as
8YEZ's. So the answer is an empty screen and the reason — which is Round 34's
result, that no deposited human PIEZO1 entry conducts.</p>
<p>And the model has never heard of the protein: a point source in free
solution, screened by a uniform buffer. The shells pass straight through the
channel and the membrane because the equation they come from does too.</p>

<h3>Predicted fluctuation — Physics panel &rarr; Colour by fluctuation</h3>
<p>Beside <b>Colour by displacement</b>, and deliberately one click from it. That
button shows how far <i>one</i> mode moves each residue; this one shows the
whole mode set's mean-square fluctuation, Σ|v|²/λ — the quantity Round 82
validated against the deposited B-factors. The two are routinely confused and
seeing them differ is the cheapest cure.</p>
<p><b>The scale is arbitrary.</b> This network has no fitted spring constant, so
only the <i>ordering</i> means anything. Whether the ordering is right is
measured under <b>Analysis &rarr; Fluctuation vs B-factor</b>, and the answer is
split: across the catalogue the network's median Spearman is <b>0.74</b> against
a burial-only control's <b>0.32</b> — but on Pearson it is <b>0.48</b> against
<b>0.39</b>. It orders residues by mobility much better than burial does and
predicts how far they move barely better. Colour by B-factor from the Model
panel to see the observation next to the prediction.</p>
"""
