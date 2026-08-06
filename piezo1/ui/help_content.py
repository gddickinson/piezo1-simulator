"""The in-application help text.

Kept as data in its own module so the dialog stays a widget and the text stays
reviewable. Every scientific claim here carries the same provenance the rest of
the project demands: where the number came from, and what it does not mean.
"""

from __future__ import annotations

__all__ = ["TOPICS", "DOC_LINKS", "SHORTCUTS", "topic_html"]

#: Documents shipped with the project, opened in the system viewer.
DOC_LINKS = [
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
    ("Roadmap", "ROADMAP.md", "Completed and planned work"),
]

SHORTCUTS = [
    ("Drag", "Rotate the model"),
    ("Shift + drag", "Pan"),
    ("Wheel", "Zoom"),
    ("R", "Reset the camera to frame the model"),
    ("Space", "Toggle spin"),
    ("Click an atom", "Identify it — chain, residue, domain, any annotation"),
    ("Ctrl+O", "Open a structure file"),
    ("Ctrl+S / Ctrl+L", "Save / load a session"),
    ("Ctrl+E", "Export an analysis report"),
    ("Ctrl+R", "Reset the panel layout"),
    ("F1", "This help"),
]

_MODEL = """
<h2>Model panel</h2>
<p>Chooses which structure is loaded and how it is drawn.</p>

<h3>Structure</h3>
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
<h2>Measure panel</h2>
<p>Click atoms in the viewport to measure between them.</p>
<ul>
<li><b>Distance</b> — two atoms.</li>
<li><b>Angle</b> — three.</li>
<li><b>Dihedral</b> — four.</li>
</ul>
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
TOPICS: list[tuple[str, str]] = [
    ("Getting started", _GETTING_STARTED),
    ("Model panel", _MODEL),
    ("Annotation panel", _ANNOTATION),
    ("Physics panel", _PHYSICS),
    ("Analysis panel", _ANALYSIS),
    ("Measure panel", _MEASURE),
    ("Limits and honesty", _HONESTY),
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
