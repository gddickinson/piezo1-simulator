"""The physics topic: what the model does, and what is model rather than data.

Split from ``help_content.py`` at the length limit, and along the seam that had
made it the longest topic by far. Everything here is about the *model* — the
dome it measures, the modes it moves along, whether the elastic network
describes the molecule at all, whether any of it is PIEZO1 rather than the
fold, and the full-length models that put prediction into the coordinates
themselves.

That last one is why this kept growing: the Completeness selector changes what
every other number in the application is measured on, so the place it is
explained has to explain what that costs.
"""

from __future__ import annotations

__all__ = ["PHYSICS"]

PHYSICS = """
<h2>Physics panel</h2>

<h3>Dome geometry</h3>
<p>Fits a sphere to the transmembrane surface and reports the radius of
curvature, depth, area and excess area. On the curved structure 7WLT this gives
<b>9.7 nm</b> against the published <b>10.2 nm</b> (Haselwandter &amp;
MacKinnon 2018) — the standing regression test that the geometry pipeline is
correct.</p>

<h3>Normal modes</h3>
<p>An anisotropic elastic network model. Springs connect C-alpha pairs within a
cutoff; the Hessian's low-frequency eigenvectors are the motions the fold makes
most cheaply.</p>
<p>Each mode is labelled <b>A</b> or <b>E</b> by its behaviour under the
three-fold rotation. This is not decoration. Membrane tension is isotropic and
therefore C3-symmetric, so <b>only A modes can couple to it at first order</b>.
The lowest A mode overlaps the observed curved-to-flat transition at
<b>0.705</b>, while E modes contribute <b>0.001</b> — the selection rule falls
out of the calculation rather than being imposed.</p>
<p><b>Colour by displacement</b> shows which parts of the protein a mode moves.
<b>Animate</b> drives the structure along it. The amplitude is a visualisation
choice and is not a physical prediction of how far the protein travels.</p>

<h3>Full-length models — Model panel &rarr; Completeness</h3>
<p><b>This is the one that changes what everything else measures.</b> The
Completeness selector sits with the structure chooser because it decides what
is <i>loaded</i>, not how it is drawn — so the dome, the pore, the normal
modes, the animations and every measurement then run on the spliced model with
no extra step. It works on <b>any</b> entry, so a full-length model is
available in every gating state the catalogue has.</p>
<ul>
<li><b>Deposited only</b> — exactly what the experiment resolved. The default.</li>
<li><b>+ AlphaFold gaps</b> — fills the unresolved stretches <i>inside</i> the
deposited range. Each is anchored on resolved residues at <b>both</b> ends, so
it is <i>interpolated</i>; gaps longer than the registered ceiling are left
empty and counted rather than invented.</li>
<li><b>+ AlphaFold blade</b> — adds the distal blade that no PIEZO1 structure
resolves. Anchored at one end only, so it is <i>extrapolated</i> and the far
end is the least constrained part of the model.</li>
<li><b>+ AlphaFold (full length)</b> — both. Roughly half the residues are
then prediction.</li>
</ul>
<p><b>You can always tell which you are looking at.</b> The structure's name
gains a suffix, so every result window and every exported report is stamped
with it. The status line begins with <tt>DEPOSITED</tt> or <tt>PART
PREDICTED</tt> and the count. An amber banner sits in the viewport whenever
prediction is on screen, and unlike the other overlay elements it cannot be
switched off. Any result computed on a spliced model opens with a caveat saying
so.</p>
<p>The three blades are grafted <b>independently</b>, one per protomer, so how
three-fold related they come out is a <i>measurement</i> — 0.0 Å on 7WLT, and
3.6 Å on the low-resolution 11ZC, which is the model telling you it is less
reliable there.</p>

<h3>The full-length model — View &rarr; Full-length model</h3>
<p>Cryo-EM resolves roughly residues <b>570–2521</b>; the distal blade has
never been resolved in any PIEZO1 structure. This grafts AlphaFold's prediction
of those <b>569 residues</b> onto the experimental core and keeps the join
visible: the measured part is flat grey, the predicted part is coloured by
AlphaFold's own <b>pLDDT</b> confidence bands, and the seam is marked.</p>
<p>It is a toggle, so it draws on top of whatever else is on screen. The
numbers behind it — the fraction of the graft clearing pLDDT 70, and how far
the two models disagree away from the seam — are under <b>Analysis &rarr;
Full-length model numbers</b>, and the status line carries them whenever the
graft is drawn.</p>
<p><b>Two populations, not one structure.</b> Only about <b>48%</b> of the
grafted region clears pLDDT 70, and the two models differ by <b>75 Å</b> over
the region they share even though the seam itself fits to <b>2.4 Å</b> — a good
local fit says nothing about the rest of a long flexible arm. Every atom
carries its source so no analysis can average across the join.</p>
<p>The prediction used matches the entry: the human model for a human
structure, the mouse model for a mouse one. There is no PIEZO2 prediction
downloaded, so a PIEZO2 entry is refused rather than given a PIEZO1 blade —
a wrong graft is not visible on screen, it is a blade in roughly the right
place in the wrong sequence.</p>

<h3>Seeing all of this</h3>
<p>The dome, the contacts, the pore, the pockets, the blade-to-gate route and
the calcium field can each be drawn on top of the structure from the
<b>View</b> menu. They have a topic of their own — <b>Drawing what was
measured</b> — because they share a failure mode: a drawn surface is more
persuasive than the number it came from.</p>

<h3>Is any of this PIEZO1? — Analysis &rarr; PIEZO2 comparison</h3>
<p><b>The catalogue now holds four PIEZOs.</b> PIEZO1 (human and mouse),
PIEZO2 (<b>human</b> and mouse), <i>C. elegans</i> PEZO-1 and <i>Drosophila</i>
PIEZO. Filter the entry list by <b>Protein</b> in the Model panel to reach
them. The invertebrates are neither PIEZO1 nor PIEZO2 — that duplication is
vertebrate — and they do not even share the 38-helix architecture, with 36 and
40 respectively.</p>
<p>Having <b>human</b> PIEZO2 matters: the comparison no longer has to cross a
species boundary as well as the paralogue. Mouse-to-mouse gives an overlap of
0.804 with PIEZO2's seventh symmetric mode; human-to-human gives <b>0.962</b>
with its <b>lowest</b>. The comparison picks the same-species partner
automatically.</p>
<p><b>PIEZO2 is the only control available</b> for the question this project
otherwise never asks: how much of the mechanism is PIEZO1, and how much is the
fold? 6KG7 was downloaded, classified and then excluded from every ensemble as
a paralogue — correct for a PIEZO1 ensemble, and not an answer.</p>
<p><b>Read the two dome blocks together.</b> Measured naively the two proteins
look very different — PIEZO2's dome 8.5 nm deep against 4.9. That is because
6KG7 resolves <b>all 38</b> transmembrane helices where a PIEZO1 entry resolves
22-26, so the two measurements trace different amounts of blade.
Coverage-matched, PIEZO2 gives 5.6 nm and falls inside the PIEZO1 range on
every quantity.</p>
<p><b>The gating mode is the fold's.</b> With the sites matched through a real
alignment and the protomer order searched rather than trusted — it is (2, 0, 1)
here, so chain labels would have been wrong — PIEZO1's lowest symmetric mode
overlaps a single PIEZO2 symmetric mode at <b>0.804</b>, with <b>0.925</b> of
it inside PIEZO2's symmetric subspace, against a shuffled-correspondence
control of 0.190. The candidate gating coordinate is not specific to
PIEZO1.</p>
<p>That is a result about generality, not a failure — and it cuts both ways.
Nothing in this mechanism distinguishes two proteins whose inactivation
kinetics and tissue roles differ. With one PIEZO2 structure it says the fold
<i>admits</i> the mechanism, not that every PIEZO uses it.</p>

<h3>Does the network describe this molecule? — Analysis &rarr; Fluctuation vs B-factor</h3>
<p>The standard check on any elastic network, and one this project had never
run until Round 82: does the predicted mean-square fluctuation track the
<b>B-factor the entry was deposited with</b>?</p>
<p><b>The column is checked before the network is.</b> A cryo-EM B-factor
absorbs local resolution and sharpening as much as motion, so three kinds of
column are refused rather than correlated — a uniform one, a <b>grouped</b> one
(3JAC carries 212 distinct values over 2,754 residues), and an AlphaFold model,
whose B column holds <b>pLDDT</b>, a confidence that runs the other way.</p>
<p><b>Read the control.</b> A buried residue moves less in any packed solid, so
every correlation is shown beside the same correlation for contact number,
which uses no network at all. Across the catalogue the network's median rank
correlation is <b>0.74</b> against the control's <b>0.32</b> and it wins on 13
of 15 entries — but on Pearson it is 0.48 against 0.39 and wins only 9 of 15.
The network orders residues by mobility much better than burial does, and
predicts how much they move barely better.</p>
<p>A <b>negative</b> control is a verdict on the entry: its B-factor rises with
burial, which no mobility does. Three entries do that, and on two of them the
network gets 0.10 — there the column is the problem, not the model.</p>

<h3>Morph</h3>
<p>Interpolates between two conformational endpoints. Three methods, each
reporting its own bond-geometry error: linear (fast, distorts bonds),
distance-restrained, and elastic-network-subspace. A morph is an
<i>interpolation</i>, not a simulated trajectory, and the panel says so.</p>
<p><b>What the two ends of the slider are.</b> Only C-alphas are interpolated.
At the far end the restrained and linear methods put every C-alpha exactly on
the target entry's, but each remaining atom is carried with its own residue and
keeps the starting structure's side-chain geometry — so the last frame is
<i>not</i> the deposited target, and loading that entry instead will not look
identical. The elastic-network method does not even reach the target: it is
confined to the subspace the network supports, which on 7WLT&nbsp;&rarr;&nbsp;7WLU
captures 95% of the change and stops about 6&nbsp;&Aring; short. The status line
states which of the two you are looking at.</p>
<p>The two endpoints are also reduced to the residues both entries resolve and
rigid-body superposed onto the displayed one, so the morph shows the target's
shape where the deposited target sits, not where it was deposited.</p>
"""
