"""The in-application help text.

Kept as data in its own module so the dialog stays a widget and the text stays
reviewable. Every scientific claim here carries the same provenance the rest of
the project demands: where the number came from, and what it does not mean.
"""

from __future__ import annotations

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

_TAGS_AND_CURRENT = """
<h2>HaloTag fusion, labelling and ion current</h2>

<h3>Seeing the tags — View &rarr; HaloTag fusion</h3>
<p>PIEZO1 imaging constructs fuse <b>HaloTag</b> to the cytosolic C-terminus,
one per protomer. <b>There is no structure of the fusion</b>, so everything
drawn here is a model and is drawn to look like one: the tag body is a sphere of
its radius of gyration rather than its fold, the linker is a straight seam in a
colour the channel never uses, and the <b>accessible volume</b> is shown as a
point cloud precisely so a single sphere is not mistaken for a determined
position.</p>
<p>Measured inputs, from PDB <b>6U32</b> (1.8 Å, TMR ligand bound): radius of
gyration <b>17.6 Å</b>, N-terminus <b>19.9 Å</b> from the centre. A C-terminal
fusion attaches to the tag's <i>N</i>-terminus, so that offset — not the radius
of gyration — sets where the body sits.</p>
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
    ("Measure panel", _MEASURE),
    ("HaloTag and ion current", _TAGS_AND_CURRENT),
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
