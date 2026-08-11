"""The topic for reproducing a published paper — Guo & MacKinnon 2017.

Split off at the length limit and along a real seam: every other topic explains
part of *this* application, while this one explains what the application can
and cannot say about somebody else's figures. It is also the only topic whose
main content is a list of refusals, and those need room to state their reasons
rather than being compressed into a sentence each.
"""

from __future__ import annotations

__all__ = ["PAPER"]

PAPER = """
<h2>Reproducing Guo &amp; MacKinnon 2017</h2>

<p>The membrane dome mechanism this whole application is built around comes
from one paper — Guo &amp; MacKinnon, <i>Structure-based membrane dome mechanism
for Piezo mechanosensitivity</i>, eLife 2017;6:e33660, PDB <b>6B3R</b>. The
10.2&nbsp;nm radius, the 120&nbsp;nm² of released area, the two-state
energetics: all of it is that paper's Figure 7.</p>

<p><b>Analysis → Guo &amp; MacKinnon 2017 figures…</b> replicates it panel by
panel and reports what it cannot. Of its 31 panels, 16 reproduce from deposited
coordinates, 3 have an <i>analogue</i> that is a different quantity, and 12 need
experimental data this project does not hold.</p>

<h3>What reproduces exactly</h3>
<p><b>Figure 7 and its supplement.</b> Every number follows from two lengths —
a sphere of radius 10.2&nbsp;nm centred 4.0&nbsp;nm above the plane the membrane
returns to — by closed-form spherical-cap geometry. An 18.8&nbsp;nm opening
against the paper's "about 18", a 6.2&nbsp;nm depth against "about 6",
397&nbsp;nm² of surface against 400, 277&nbsp;nm² projected against 280,
121&nbsp;nm² released against 120, and 153&nbsp;k<sub>B</sub>T of bending energy
against "approximately 150". Read that agreement for what it is: a check that
the arithmetic is right, <i>not</i> a measurement of PIEZO1. The idealised dome
is a shape the authors chose to make the energetics tractable, and they say so.
Our own measurement of 6B3R's dome gives different areas, and both are
reported.</p>

<p><b>Figure 4a</b> — the claim that one subunit fits a plane and the trimer
does not. Measured as residuals and split into two terms: how far one protomer
departs from its own best plane, and what is left when each protomer is
<i>made</i> planar and the three are left where the symmetry puts them. On a
curved entry the second term dominates; on the flattened 7WLU it does not,
which is what makes the measurement about curvature rather than about trimers.
The beam comes out at 56° against the paper's "about 60", and opens towards 90°
when the channel flattens.</p>

<p><b>Figure 4—supplement 1</b> — the cap's acidic patch against the loops'
basic patch. Both the hydrogen bonds and the salt bridges are found, in all
three protomers, and every one of them is domain-swapped, which is the paper's
specific claim. E2257–R1762 is reproduced; D2264–R1761 is not, at 6.4&nbsp;Å
charge-centroid separation against a 5.5&nbsp;Å cutoff though its closest atoms
are 4.6&nbsp;Å apart. The two conventions disagree about that contact and the
paper does not say which it used.</p>

<p><b>Figure 3 and its supplements</b> — the topology and the hydropathy behind
it. The 4-TM repeat the authors infer is supported against a shuffled control
in both mammalian PIEZOs and is <i>not</i> in <i>C.&nbsp;elegans</i> PEZO-1 or
<i>Drosophila</i> PIEZO. See <b>Analysis → Topology diagram</b> for the figure
itself.</p>

<h3>Where the numbers differ, and why</h3>
<p><b>Figure 6b, the pore profile.</b> The same three residues constrict, in the
same order, and the closed verdict agrees — but our radii are about 0.6&nbsp;Å
wider than the published ones throughout. Guo &amp; MacKinnon used HOLE; this
project's profiler is an independent implementation. A systematic offset between
two pore algorithms is expected, and it is reported rather than absorbed. If it
ever became zero, the profiler would have been fitted to the paper.</p>

<h3>The three analogues — do not put these beside the originals unread</h3>
<p><b>Figure 2a,b</b> are 2D class averages: thousands of real particle images,
aligned and averaged. What we compute is the projection of the atomic model,
which is the quantity a class average <i>estimates</i>. No contrast transfer
function, no defocus, no solvent, and — decisively for the side view — no
detergent micelle. Figure 2b's dome-shaped envelope is substantially
micelle.</p>

<p><b>Figure 4b</b> shows that micelle directly, from the unsharpened map. What
this application can draw is the fitted dome surface: a model of the same object
arrived at from the protein rather than from the density.</p>

<p><b>Figure 4c</b> was computed with APBS, which solves the Poisson–Boltzmann
equation on a grid with a dielectric boundary between a low-dielectric protein
interior and high-dielectric solvent. Ours is screened Coulomb through a
uniform solvent dielectric — the same approximation PyMOL and Chimera offer as
"coulombic colouring". It gets the sign and the pattern right and
systematically <b>under-estimates the magnitude</b>: on 6B3R nothing reaches
the panel's ±5&nbsp;k<sub>B</sub>T/e saturation, where the published surface
visibly saturates. That is the expected direction of the approximation, not a
disagreement about the protein.</p>

<h3>What cannot be reproduced at all</h3>
<p>Six panels need the cryo-EM map or the half maps (the FSC curve, the
local-resolution colouring, the three local-density panels, the stereo view of
the constrictions). Four need micrographs of proteoliposomes. One needs P2X and
ASIC coordinates, which are deliberately absent: this project's structure
catalogue, its numbering checks and its entity classifier all assume a PIEZO,
and admitting two unrelated channels to make one panel would weaken every one
of those guards. Figure 1 is a drawing.</p>

<p>Each is listed with its reason. A tool that quietly covered the tractable
parts of a paper would leave a reader assuming the rest.</p>

<h3>The three Figure 4 views — View menu</h3>

<p><b>Micelle density (modelled).</b> Figure 4b is the unsharpened cryo-EM map
at 6 sigma, and what it shows around the protein is <i>detergent</i>. This
project holds no map, so what is drawn is the surface a fixed distance outside
the hydrophobic transmembrane belt — apolar side chains of the annotated TM
helices, which is the band the paper itself describes. Two things follow, and
the status line says both. The <b>shell thickness</b> is a registered parameter
and carries no information: an offset surface around a sphere is a sphere with
the radius increased by exactly the offset. The <b>curvature</b> is a sphere
fitted to the belt atoms themselves and <i>is</i> a measurement of the protein
— 9.8&nbsp;nm on 6B3R, against the paper's 10.2&nbsp;nm idealisation and this
project's own 10.8&nbsp;nm dome fit.</p>

<p><b>Planar membrane (one protomer).</b> Figure 4a's two grey rules, fitted to
one protomer's transmembrane band. The paper's point is the <i>contrast</i>, so
the status line reports both residuals and the slab each would need against a
real 36&nbsp;Å bilayer: about 42&nbsp;Å for a protomer and 60&nbsp;Å for the
trimer on 6B3R. Every point set has a best-fit plane, so the number to read is
the residual, not the lines.</p>

<p><b>Colour by electrostatics.</b> Figure 4c's surface, on the same fixed
scale — red at −5&nbsp;k<sub>B</sub>T/e, white at zero, blue at +5. Fixed, not
auto-ranged, and that is deliberate: the application's other value colourings
stretch themselves to the data, which for a potential would paint an almost
neutral protein in full red and blue and make it incomparable with a published
surface. Atoms with no accessible surface have no potential and are painted the
same neutral white as a measured zero; the status line says how many.</p>

<h3>The topology diagram — Analysis &rarr; Topology diagram (Ctrl+Shift+T)</h3>
<p>Figure 3a as a live diagram of whichever entry is loaded: 38 transmembrane
helices in a membrane band, grouped into the nine 4-TM units, with the cap
above and the beam and cuff below. Two things it does that the published figure
cannot.</p>
<p><b>It follows the coordinates.</b> A helix the entry does not model is drawn
<i>dashed</i>, not dropped. That is not a stylistic choice: dropping one would
put TM13 where TM1 belongs and silently renumber every helix after it, on a
picture that still looked entirely reasonable. 6B3R greys out TM1–12 and 7WLT
greys out TM1–16, both read from the model rather than written down.</p>
<p><b>The boxes are a selection.</b> Ticking a unit boxes it exactly as Figure
3b does <i>and</i> selects its residues on the 3-D model, so "which part of the
blade is this" is answered by looking. Shift-click a helix to do the same in one
step. Boxing units that are not adjacent highlights the whole span between the
first and the last, and the status line says so — the 3-D view highlights a
range, and a gapped selection drawn as one range would misrepresent what is
lit.</p>
<p>Every element carries its residue range in the loaded entry's own numbering,
named on the diagram. Export writes a PNG at twice the on-screen size.</p>

<h3>Numbering</h3>
<p>The paper is in <b>mouse</b> numbering throughout (UniProt E2JF22, 2547 aa)
and 6B3R is deposited in it. Most of this project's human entries are not. Every
residue number quoted above is mouse; the human equivalent is not a fixed
offset away, and the conversion goes through <code>piezo1.core.sequence</code>
rather than arithmetic.</p>
"""
