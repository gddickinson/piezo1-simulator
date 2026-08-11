"""The in-application help text.

Kept as data in its own module so the dialog stays a widget and the text stays
reviewable. Every scientific claim here carries the same provenance the rest of
the project demands: where the number came from, and what it does not mean.
"""

from __future__ import annotations

from .help_topics_tags import TAGS_AND_CURRENT

__all__ = ["TOPICS", "DOC_LINKS", "SHORTCUTS", "topic_html"]

#: Documents shipped with the project, opened in the system viewer.
DOC_LINKS = [
    ("What was established", "docs/CONCLUSION.md",
     "The one page: what this project showed, and what it could not"),
    ("README.md", "README.md", "Overview, install, and the headline results"),
    ("Scientific basis", "docs/SCIENCE.md",
     "Every parameter, its source, and the open gaps"),
    ("Navigation map", "INTERFACE.md",
     "What each module contains and how they connect"),
    ("Headless API", "docs/NOTEBOOK.md",
     "Driving the engine from a notebook, with the traps listed"),
    ("Pre-registration", "docs/PREREGISTRATION.md",
     "The frozen hypothesis for the blind variant test"),
    ("Validation result", "docs/VALIDATION.md",
     "The Round 7 null result, reported in the pre-registered order"),
    ("Negative-result protocol", "docs/NEGATIVE_RESULT_PROTOCOL.md",
     "Standing policy on power, multiplicity and what a null may claim"),
    ("Bibliography", "docs/REFERENCES.md", "All 60 verified references"),
    ("Session log", "SESSION_LOG.md", "What was done each round, and why"),
    ("Roadmap", "ROADMAP.md", "What is still open — 12 items across six rounds"),
    ("Completed rounds", "docs/ROADMAP_COMPLETED.md",
     "The result of every finished round, 75 of them"),
]

#: Every binding the application has. The list is complete rather than
#: representative, and a test enforces that: `test_ui_controls` fails if a menu
#: action carries a shortcut this does not mention, or if a key the viewport
#: handles is missing. Four bindings were undocumented until that check existed
#: — middle-drag, right-drag, O and the +/- pair.
SHORTCUTS = [
    # Mouse — the viewport
    ("Drag", "Rotate the model"),
    ("Shift + drag", "Pan"),
    ("Middle-drag", "Pan"),
    ("Right-click", "Open the context menu: select this residue, this chain "
                    "or all three copies, add it to a measurement, centre on "
                    "it, copy its label — plus representation, colouring and "
                    "the view toggles"),
    ("Right-drag", "Zoom — dragging never opens the menu"),
    ("Wheel", "Zoom"),
    ("Click an atom", "Identify it and mark it in gold — chain, residue, "
                      "domain, any annotation. Dragging to rotate never picks, "
                      "however the drag ends. Press Start picking in the "
                      "Measure panel to measure between clicks instead"),
    # Keys — the viewport, which must have focus (click it once)
    ("R", "Reset the camera to frame the model"),
    ("O", "Switch between perspective and orthographic projection"),
    ("Space", "Toggle spin"),
    ("+ / -", "Grow or shrink the drawn atoms"),
    # Menu shortcuts
    ("Ctrl+O", "Open a structure file"),
    ("Ctrl+S / Ctrl+L", "Save / load a session"),
    ("Ctrl+E", "Export an analysis report"),
    ("Ctrl+Shift+S", "Sequence window"),
    ("Ctrl+D", "Display options — what the overlay shows"),
    ("Ctrl+P", "Parameters — every registered number, with its source"),
    ("Ctrl+R", "Reset the panel layout"),
    ("F11", "Presentation mode (full screen)"),
    ("F1", "This help"),
    ("F2", "Guided tour"),
    ("Ctrl+Q", "Quit"),
]

_MODEL = """
<h2>Model panel</h2>
<p>Chooses which structure is loaded and how it is drawn.</p>

<h3>Structure</h3>
<p><b>Four filters narrow the catalogue</b> — protein, species, state and
gating — and each box lists the values the records actually take, so a new
entry appears in it without anyone remembering to add it. The count underneath
says how much is hidden.</p>
<p><b>Protein</b> is the one worth knowing about: the catalogue contains one
PIEZO2 entry, <b>6KG7</b>, and it is filed as <i>mouse</i> like fifteen other
entries, so before this it was reachable only by knowing which one it was. It
is identified by measurement — each file's own residue names scored against
every reference sequence — not by a curated label. <b>Gating</b> is the other
useful one: every deposited human entry is closed, and 11ZC is the only
open-like structure in the set.</p>
<p>Filters never hide a structure something else asked for by name: opening one
with <tt>--structure</tt>, restoring a session or starting a morph clears
whatever filter was in the way rather than silently leaving the previous model
on screen.</p>
<p>Twenty-one curated PIEZO structures. Each entry states its <b>gating
state</b>, resolution, the residue range actually resolved, its <b>numbering
species</b>, bound ligands and the paper it came from. That matters more than it
looks: most functional literature uses <i>mouse</i> numbering and most disease
variants use <i>human</i>, and the offset between them is <b>not constant</b> —
it varies from 0 to +26 across twelve blocks. The application converts through a
real alignment, never by adding a constant.</p>

<h3>Representation</h3>
<ul>
<li><b>Cartoon</b> — secondary structure ribbons. Helices and strands are
assigned from C-alpha geometry (P-SEA), because most of these entries are
C-alpha traces without the backbone atoms DSSP needs.</li>
<li><b>Spheres / sticks / ball-and-stick</b> — atomic detail. Drawn as ray-cast
impostors, so the spheres are mathematically exact rather than tessellated.</li>
<li><b>Backbone</b> — a fast trace for orientation on a 120,000-atom trimer.</li>
</ul>

<h3>Colour by</h3>
<ul>
<li><b>Domain</b> — the 17 architectural domains, each with its provenance
(UniProt, derived by rule, or literature).</li>
<li><b>Chain</b> — the three protomers.</li>
<li><b>Secondary structure</b>, <b>B-factor</b>, <b>pLDDT</b> (AlphaFold
confidence), <b>element</b>.</li>
<li><b>Computed value</b> — whatever the Analysis panel last produced.</li>
</ul>
"""

_ANNOTATION = """
<h2>Annotation panel</h2>
<p>Curated knowledge mapped onto the structure. Selecting anything highlights it
in the viewport and explains it, with the PMID.</p>

<h3>Domains</h3>
<p>17 domains with ranges in both numbering systems. Each carries a
<b>confidence</b> and a source; boundaries we could not verify are labelled
rather than quietly rounded.</p>

<h3>Sites</h3>
<p>37 residues in 11 functional groups — the hydrophobic gate, the selectivity
glutamates, the CTD constrictions, the Yoda1 pocket, the PIP2 cluster, the basic
patches. Every one was verified against the UniProt sequence when the resource
was built.</p>

<h3>Variants</h3>
<p>68 curated variants, each with its wild-type residue checked against Q92508,
its phenotype classification, and which deposited structures resolve it. The
<b>classification</b> field distinguishes gain-of-function (22),
loss-of-function (17), engineered mutations (15), variants of uncertain
significance (8) and blood-group antigens (6). Many variants are
<i>not resolved</i> in any structure — the panel says so rather than showing
nothing.</p>
"""

_PHYSICS = """
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

<h3>Is any of this PIEZO1? — Analysis &rarr; PIEZO2 comparison</h3>
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
"""

_ANALYSIS = """
<h2>Analysis panel</h2>

<h3>Pore</h3>
<p>The radius of the largest sphere that fits at each height along the
conduction axis, with the probe tethered near the axis. Without that leash the
clearance function has no interior maximum and the answer escapes to ~6000 Å —
a true maximum, and useless.</p>
<p>The <b>hydrophobicity trace</b> drawn against it is the point. Radius alone
predicts the conductive state at <b>AUROC 0.59</b>; radius combined with
hydrophobicity reaches <b>0.91</b> (Rao et al. 2019). A pore can be wide enough
for a hydrated ion and still block, because a hydrophobic neck expels liquid
water.</p>
<p>The verdict separates two independent ways to be shut:</p>
<ul>
<li><b>Sterically occluded</b> — narrower than a water molecule (0.15 nm).</li>
<li><b>Hydrophobic gate</b> — wide enough, but the lining would dewet.</li>
</ul>
<p>They are reported separately because PIEZO1 has structures that are one and
not the other: 7WLU and 8IXO are sterically shut with hydrophilic linings.</p>
<p>Click anywhere on the plot to select the residues lining the pore there.</p>

<h3>Pockets</h3>
<p>Delaunay alpha-sphere detection — the fpocket construction, reimplemented in
numpy — with a burial filter. Without that filter the largest "pocket"
percolates over the entire surface and comes out at 408,000 Å&sup3; across 601
residues. Volumes are Monte-Carlo <i>unions</i>, not sums of overlapping
spheres.</p>

<h3>Conservation</h3>
<p>Per-residue Shannon entropy over vertebrate PIEZO1 orthologs, one sequence
per species, anchored on the human reference. Positions where fewer than 70% of
orthologs align are dropped: their value measures the alignment rather than
selection pressure on the residue.</p>

<h3>Coupling to the gate (PRS)</h3>
<p>Perturbation response scanning — push each residue in turn, measure how much
the gate moves. Requires normal modes first. The anchor domain comes out as the
force-transmission hub, which agrees with the independent conservation result
and with the lever model from the literature.</p>
<p>Values are averaged over the three protomers, since the trimer is symmetric
and a per-protomer difference is numerical noise.</p>
"""

_MEASURE = """
<h2>Selecting and measuring</h2>

<h3>What a click does</h3>
<p>Clicking an atom in the 3-D view always does two things: it names the
residue in the status bar — with its domain, any annotated site, and any
variant reported there — and it <b>marks that residue in gold on the model</b>.
The mark is the copy you clicked, not all three protomers: a click means one
specific atom, where selecting a residue from the Annotations panel means the
residue number and therefore all three.</p>
<p>Rotating does not select. A click only registers if the mouse has barely
moved between press and release, so dragging to turn the structure will not
pick anything by accident, however the drag ends.</p>

<h3>Right-click — the context menu</h3>
<p>Right-clicking opens a menu for whatever is under the cursor. On a residue
it offers to select that copy, the same residue number in all three protomers,
or the whole chain; to add the atom to a measurement (<i>arming picking for
you</i>, so this is the short way in); to centre the view on it; and to copy its
label. Any variant reported at that residue is named. On empty background the
residue entries are simply absent.</p>
<p>Either way the menu carries <b>Representation</b>, <b>Colour by</b>, and the
view toggles — reset, spin, orthographic projection. Those drive the Model
panel's own controls rather than a second copy of them, so the panel and the
menu can never disagree about what is on screen.</p>
<p>Opening the menu selects nothing. It identifies the residue so the entries
can name it, and dismissing it leaves the model exactly as it was. Right-<i>drag</i>
still zooms and never opens the menu.</p>

<h3>Measuring — Measure panel</h3>
<p>Measuring needs picking to be <b>armed</b>, because clicks already mean
"tell me about this residue" and a measurement tool that silently consumed
them would break inspection. So:</p>
<ol>
<li>Choose <b>distance</b> (two atoms), <b>angle</b> (three) or
<b>dihedral</b> (four).</li>
<li>Press <b>Start picking</b>. The button changes to
<i>Picking — click atoms</i>.</li>
<li>Click the atoms in the 3-D view. Each one appears in the
<b>Selection</b> table straight away, in blue, with the number still to go —
and is marked in blue on the model.</li>
<li>On the last atom the row resolves to a value and units.</li>
</ol>
<p>Pressing the button again abandons any half-made selection. Selecting the
blue pending row and pressing <b>Delete</b> does the same without disarming;
selecting a completed row and pressing Delete removes that measurement.</p>
<p>Measurements accumulate in the table and export to CSV. The measurement
logic is deliberately free of Qt so it can be tested without a display; the
regression case is the C2411–C2415 disulfide, which must come out at
<b>2.04 Å</b>.</p>
"""

_HONESTY = """
<h2>What this application will not do</h2>
<p>The project's central claim — predicting gain- versus loss-of-function from
structure — <b>has been tested and did not work</b>, and that result is shipped
rather than hidden.</p>

<h3>The blind test</h3>
<p>A pre-registered comparison of an elastic-network ΔΔG against 25 phenotyped
variants returned <b>p = 0.234, Cliff's delta −0.083, AUROC 0.542</b>. The null
hypothesis was not rejected. The diagnostic is the useful part: <b>99.8% of the
predictor's variance is between-position</b>, meaning it reports <i>where a
residue sits</i> rather than <i>which substitution occurred</i>. All four R2456
substitutions score as "softening" although three are gain-of-function and one
is loss-of-function.</p>

<h3>What that null is entitled to claim</h3>
<p>With 16 gain- against 9 loss-of-function variants, the design reached 80%
power only at |Cliff's delta| ≥ 0.55 — beyond "large". Power at the observed
effect was <b>0.13</b>. So the result excludes a <i>large</i> mechanical effect
and says little about a small or medium one. Detecting a medium effect would
need about <b>98</b> phenotyped variants against the 25 available; the binding
constraint is data, not method.</p>

<h3>Numbers that have been corrected</h3>
<p>The linearised membrane footprint was reported at 622 nm² of excess area.
Solved without the small-slope approximation — which PIEZO1's 63° contact slope
badly violates — it is <b>179 nm²</b>. The earlier figure was 3.5× too large and
the conclusion drawn from it reversed. Both are left visible in the session log,
because the correction is the informative part.</p>

<h3>Provenance</h3>
<p>Every parameter in this application carries its source. Values that could not
be verified are labelled <tt>UNVERIFIED</tt> rather than rounded into
plausibility. Anything interpolated for visual purposes says so.</p>
"""

_GETTING_STARTED = """
<h2>Getting started</h2>
<ol>
<li>Pick a structure in the <b>Model</b> panel. <b>8YEZ</b> is the closed human
channel; <b>11ZC</b> is flattened and open-like.</li>
<li>Drag to rotate. Press <b>R</b> to reframe, <b>space</b> to spin.</li>
<li>Click any atom to identify it.</li>
<li>In <b>Physics</b>, press <i>Measure dome</i> — you should get about 9.7 nm
for a curved structure.</li>
<li>In <b>Analysis → Pore</b>, press <i>Compute pore profile</i>. On 8YEZ the
bottleneck is 0.95 Å and the verdict is non-conductive.</li>
<li>Compute normal modes in <b>Physics</b>, then return to <b>Analysis</b> to
scan coupling to the gate.</li>
</ol>

<h3>Panels</h3>
<p>Every panel is a dock. Drag its title bar to move it to any edge, drag it out
of the window to float it as an independent window, or close it and bring it
back from <b>View → Panels</b>. <b>View → Reset layout</b> (Ctrl+R) restores the
arrangement the application ships with.</p>

<h3>If something is missing</h3>
<p>Structures, sequences and the hydrophobic-gating grid are downloaded, not
bundled. Run <tt>python -m piezo1.io.fetch</tt> once. Analyses degrade to
"unavailable" rather than failing when their data is absent.</p>
"""

#: Ordered so the list reads as a tour rather than an index.

_FRAMING = """
<h2>Where structures sit, and showing more than one</h2>

<h3>Structure alignment — View &rarr; Structure alignment</h3>
<p>Nothing about a PDB coordinate frame is canonical: across the twenty
downloaded entries, pairs sit <b>29–147 Å</b> apart before any alignment, which
is why loading one after another used to look like the molecule jumping.</p>
<ul>
<li><b>As deposited</b> — the file's own frame.</li>
<li><b>Canonical</b> — a frame defined by the structure's <i>own</i> three-fold
symmetry: axis vertical, cytosolic side down, centred. Needs no other structure,
so it works across species and for PIEZO2. This brings those same pairs to
<b>0.9–25 Å</b>, within about an angstrom of what a least-squares superposition
achieves — what is left is real conformational difference, not framing.</li>
<li><b>Superpose on the loaded structure</b> — least-squares onto the first
structure loaded. Needs a shared residue numbering, so it falls back to
canonical across species rather than aligning non-equivalent residues.</li>
</ul>

<h3>Showing several structures — View &rarr; Show multiple structures</h3>
<p><b>Off by default</b>, because two entries in the same frame sit on top of
each other and one left behind reads as extra density. With it on, loading keeps
the previous structure as a companion in its own colour, and the Model panel
lists what is drawn. Turning it off drops the extras.</p>
<p>This is <i>display</i>, not measurement. The overlay feature superposes one
nominated structure and reports an RMSD; this simply draws several. Every
analysis runs on the primary structure, whatever else is on screen.</p>
<p><b>This is how to see PIEZO1 and PIEZO2 together.</b> The shared frame is
each structure's own three-fold symmetry, so it needs no correspondence between
them and works across the paralogue. The measuring <i>overlay</i> does not and
refuses the pair: it joins on residue number, and mouse Piezo1 residue 1500 is
not mouse Piezo2 residue 1500. For a measured comparison of the two, use
Analysis &rarr; PIEZO2 comparison, which matches them through a real
alignment.</p>
"""


def _hazards_topic() -> str:
    """Built from the hazard register, so the guide cannot overstate the guards.

    Hand-written HTML here would be a second copy of `ui/hazards.py` and would
    drift from it the first time a guard changed.
    """
    from .hazards import HAZARDS

    rows = []
    for hazard in HAZARDS:
        rows.append(
            f"<li><b>{hazard.scenario}</b><br>"
            f"<i>What would be wrong:</i> {hazard.wrong}<br>"
            f"<i>What stops it ({hazard.status}):</i> {hazard.guard}</li>")
    return ("\n<h2>Ways to get a wrong number, and what stops them</h2>\n"
            "<p>Every analysis here can be run in a situation where its answer "
            "would not mean what it appears to. These are the ones known "
            "about, each with the mechanism that prevents or marks it. "
            "Anything not on this list has not been audited.</p>\n<ul>\n"
            + "\n".join(rows) + "\n</ul>\n")


TOPICS: list[tuple[str, str]] = [
    ("Getting started", _GETTING_STARTED),
    ("Model panel", _MODEL),
    ("Annotation panel", _ANNOTATION),
    ("Physics panel", _PHYSICS),
    ("Analysis panel", _ANALYSIS),
    ("Measure panel — selecting atoms", _MEASURE),
    ("HaloTag and ion current", TAGS_AND_CURRENT),
    ("Framing and multiple structures", _FRAMING),
    ("Limits and honesty", _HONESTY),
    ("Wrong numbers, and what stops them", _hazards_topic()),
]


def topic_html(title: str, body: str) -> str:
    """Wrap a topic in styling that matches the dark theme."""
    return f"""<html><head><style>
      body {{ color:#c8ccd4; font-size:13px; line-height:1.5; }}
      h2 {{ color:#7aa7ff; font-size:16px; margin-bottom:2px; }}
      h3 {{ color:#e8ecf3; font-size:13px; margin-top:14px; margin-bottom:2px; }}
      b {{ color:#f0f3f8; }}
      tt {{ color:#f2a65a; }}
      li {{ margin-bottom:4px; }}
    </style></head><body>{body}</body></html>"""
