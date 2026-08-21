"""The guide's topic for the Explore button.

Its own module for the reason the other topics have theirs: ``help_content``
is at the length limit, and this one is about a surface every analysis window
shares rather than about one panel. What it has to get across is not how to
use it — a list and a slider need no explaining — but **what the four kinds of
exhibit are evidence of**, because on screen they look identical.
"""

from __future__ import annotations

__all__ = ["EXPLORE"]

EXPLORE = """
<h2>Explore these findings</h2>

<p>Every analysis under the <b>Analysis</b> menu opens a window of numbers, and
every one of those windows has an <b>Explore these findings…</b> button. A
number is where a piece of reasoning ends; the button is the rest of it — the
figure it is drawn in, the model it came out of, and the same result on the
structure.</p>

<h3>Four kinds, and they are not the same evidence</h3>
<p>Each exhibit says which kind it is, what it rests on, and — in amber at the
bottom — what it must <b>not</b> be read as. That last line is the one to read
first: a picture is more persuasive than the number behind it, and a curve you
have just produced by moving a slider is more persuasive still.</p>
<ul>
<li><b>Figure.</b> A generated picture from <tt>docs/img</tt>. Regenerable, so
a clone that has not built one shows the command that builds it rather than a
broken image.</li>
<li><b>Chart.</b> Drawn from the result <i>already in the window</i> — never
recomputed. The picture and the table are therefore always of one run.</li>
<li><b>Simulation.</b> A model you drive. Where a control is a registered
parameter it says so and starts at the registry's own default, marked on the
plot. Moving it is a <b>sensitivity</b>, not a measurement.</li>
<li><b>On the model.</b> A button that turns on the matching overlay in the 3-D
view. It presses the same menu entry or panel button you would press yourself,
so the View menu can never disagree with what is drawn.</li>
</ul>

<h3>What a slider cannot do</h3>
<p>Nothing in this window writes to the parameter registry. Explore the ion
current across every value of the in-pore diffusivity you like: the reports,
the exported numbers and the claims verifier are untouched, and the application
goes on quoting the documented defaults. To <i>change</i> a parameter, use
<b>Options → Parameters</b>, which marks the override in amber and makes every
report say so.</p>

<h3>Showing it on the structure</h3>
<p>This application exists because PIEZO1's shape is its mechanism, so most
results can put something in the 3-D view: the entry the result is about, a
second structure beside it or superposed on it, one part of the assembly on its
own, the residues a finding is made of, a recolouring, or the gating motion.</p>
<p>Three of those are worth knowing about directly, because they are useful
outside the result windows too:</p>
<ul>
<li><b>Fit on the pore module only</b> (Overlay panel). Puts two channels'
pore modules on top of each other and then <i>measures</i> where the blades
land. It is the only superposition mode that works across paralogues — by
residue number PIEZO1 and PIEZO2 do not correspond at all, and that fit is
refused — and the status line reports the core, the blades and the ratio
between them together. PIEZO1's own curved-to-flat pair gives 19×.</li>
<li><b>Show component</b> (View menu). One named part of the assembly, up
close. It <i>hides</i> rather than subsets: every analysis still runs on the
whole trimer, which is what the status line says on each switch.</li>
<li><b>Mark census positions</b> (View menu). The pathogenic pore positions and
the two positions equivalent to a PIEZO2 disease residue, converted into the
loaded entry's own numbering. These are the piezo_genes census's positions, not
this project's measurements, and an entry whose numbering they cannot be
carried into is refused rather than marked.</li>
</ul>
<p>Anything marked on the structure is marked in <b>this entry's numbering</b>.
Most of the catalogue is mouse-numbered and the human-to-mouse offset is not
constant — it reaches 26 residues — so a human number on a mouse entry lands on
a real, wrong residue rather than on nothing. A variant keeps its published
name and moves its position: R2456H is R2456H everywhere, and on 7WLT it sits
at 2482.</p>

<h3>Some worth starting with</h3>
<ul>
<li><b>Ion permeation.</b> How far the constriction would have to open before
the model passes the published 25–30 pS — with this entry's own bottleneck
marked, and the two unmeasured inputs on sliders that move the answer by a
factor of six.</li>
<li><b>Calcium nanodomain.</b> Free calcium against distance, with the sensor's
K<sub>D</sub> and the 90%-occupancy distance drawn on it. The whole channel
sits inside both, which is the result.</li>
<li><b>Guo &amp; MacKinnon 2017.</b> Figure 7c as arithmetic — flatten the
idealised dome and watch the released area and the bending energy — and Figure
7d's tension response with the area change on a slider.</li>
<li><b>R2456H against wild type.</b> One bar inside a band. The band is three
independent wild-type entries, and it is wide enough to hide a real
difference.</li>
</ul>

<p>An analysis whose exhibits are all charts of a result it could not produce
shows the reason instead of an empty plot: a pore the wetting model calls shut
has no conductance to draw, and saying so is the honest picture.</p>
"""
