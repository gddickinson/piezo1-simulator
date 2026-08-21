"""The topic for results this application did **not** produce.

Split off at ``help_content.py``'s length limit and along a seam none of the
other topics share: every one of them explains something this project measured,
and this one explains something it imported. That difference is the whole
content of the topic, because a number on screen looks the same either way.
"""

from __future__ import annotations

__all__ = ["FAMILY"]

FAMILY = """
<h2>The PIEZO family census — somebody else's results, answered back</h2>

<p>A sibling project, <b>piezo_genes</b>, ran a 194-genome, eukaryote-wide
census of the PIEZO family. It established the family's true range, that
vertebrates have a <b>third</b> PIEZO gene the databases largely missed, and
which parts of the protein half a billion years of evolution has refused to
change. It is a <i>sequence</i> project: no coordinates, no physics, and no way
to ask <i>why</i> a residue is conserved.</p>

<p>This application is the mirror image. <b>Analysis &rarr; The PIEZO family
census</b> shows the thirteen statements that were imported; the four entries
below it are the ones that measure something here. The full account is in
<code>docs/FAMILY.md</code>.</p>

<h3>Nothing in the census window was measured here</h3>

<p>Every statement in it came from the other project and was re-verified against
that project's own result files when this application's resource was built —
thirty-two numbers, re-read on every build. If the census is not on disk the
importer <b>refuses to write at all</b> rather than re-stamping the resource
with a fresh date and nothing behind it. An imported finding whose source has
since been corrected is a confident quotation of a superseded value, and nothing
else in the application would notice.</p>

<h3>Evolutionary constraint</h3>

<p><b>Analysis &rarr; Evolutionary constraint</b>, and <b>View &rarr; Colour by
evolutionary constraint</b> for the picture. The per-residue values are the
census's — Jensen-Shannon divergence over 174 genome-backed PIEZO1 orthologues —
and the partition is ours. On our own domain boundaries the conserved core is
still the pore machinery: anchor 0.83, CTD 0.81, inner helix 0.79, against
THU1's 0.63. The two projects' partitions put the anchor <b>141 residues
apart</b>, so that agreement is a property of the protein rather than of a
boundary choice.</p>

<p>Two things about the colouring. The scale is <b>fixed</b> at 0&ndash;1, so
the same protein is not repainted depending on how much blade the entry
resolved and two entries stay comparable. And an <b>unscored residue is grey,
not dark</b>: the blade tips are where alignment coverage is worst and where low
constraint is exactly the claim being made, so a coverage hole must not be
readable as a finding. A PIEZO2 or invertebrate entry is <b>refused</b> rather
than coloured by whatever sits at those residue numbers in PIEZO1.</p>

<p>One census finding does <b>not</b> survive here, and the window says so. Its
result that the distal blade is more conserved than the proximal one holds on
its own chain-cut bands and <b>reverses</b> on the transmembrane units: the
bands are 29% and 77% inter-unit linker respectively, and linker scores the same
either side (0.517 against 0.515). What is left after the composition is removed
is the opposite gradient.</p>

<h3>Where human disease sits</h3>

<p><b>Analysis &rarr; Where disease sits</b>. The census found pathogenic
missense concentrating in the pore module (odds ratio 3.9, P = 0.0014). This
re-tests it on PIEZO1 alone against <b>gnomAD population missense</b> rather
than ClinVar benign labels — variation that exists in people rather than
variation somebody classified as harmless, which is a better control for the
ascertainment problem the census names in its own caveat.</p>

<p>The result is reported under <b>both</b> domain partitions because the answer
follows: odds ratio 3.63 (P = 0.0033) on the census's boundaries and 1.60
(P = 0.25) on ours. The two disagree about <b>120 residues, 2057&ndash;2176</b>,
and those 120 carry six pathogenic positions including E2117 and T2127. Reading
one row and not the other would be reporting a boundary choice as a finding.</p>

<h3>Core and periphery</h3>

<p><b>Analysis &rarr; Core and periphery</b> superposes a partner on the loaded
entry <b>by the pore module alone</b> and then measures where the blades land.
A fit on everything shared spreads the error and makes a pair with a rigid
common core look the same as a pair with no common core at all; a core-only fit
asks a directional question and can <b>fail</b>, in which case no splay ratio is
reported.</p>

<p>The control is what makes the answer readable: an AlphaFold monomer splays
<b>7.2&ndash;9.1&times;</b> from an experimental structure <i>of the protein it
is a model of</i>, while three experimental cross-paralogue pairs splay
<b>0.8&ndash;2.5&times;</b>. So a large splay measured against a predicted model
is a statement about the model. Within one protein, 7WLT against the flattened
7WLU splays <b>19&times;</b> — core-conserved and periphery-free is what
PIEZO1's own gating motion looks like.</p>

<p>For a PIEZO1-versus-PIEZO2 pair the window also locates the census's most
striking clinical result: PIEZO1 <b>R2456</b> (hereditary xerocytosis) and
PIEZO2 <b>R2686</b> (Gordon syndrome) are the same residue of the same machine,
as are R2488 and R2718. Three things are checked in order — that this
application's own alignment pairs the same residues, the C&nbsp;alpha separation
after the fit, and the part that is actually the evidence: whether the claimed
partner is the <b>nearest residue of the other paralogue</b>. The whole pore
module superposes, so proximity is what every pair gives; only the register
distinguishes a correspondence from an error one residue along.</p>

<h3>piezo3 — the third vertebrate PIEZO</h3>

<p><b>Analysis &rarr; piezo3</b>. Human piezo3 has been the pseudogene
<code>PIEZO1P2</code> since before the primate radiation, so the only
coordinates the paralogue has anywhere are one AlphaFold model of the zebrafish
protein. It carries the identical human residue at all <b>fourteen</b>
pathogenic pore-module positions, checked here against a different UniProt
record of the same gene from the one the census scored.</p>

<p>Run through the pipeline it gives a dome radius of 10.8&nbsp;nm against
7WLT's 9.7&nbsp;nm by the identical route — and <b>neither number is evidence
about piezo3</b>. It is a monomer, so a trimer has to be built from a deposited
template, and 96% of the resulting departure from planarity is that template's
arrangement. What the numbers <i>can</i> do is fail, and they did not: the
protomer arranges into a closed trimer with an axis and a continuous lumen,
which a protein not built like a channel need not have done. That is a negative
that survived, not a positive demonstrated. No current has ever been recorded
from any piezo3.</p>

<h3>The one question neither project could ask alone</h3>

<p>Is a residue's evolutionary constraint predicted by how mechanically coupled
it is, or only by how <b>buried</b> it is? Buried residues are conserved in
every protein ever studied for reasons that have nothing to do with
mechanotransduction, and burial correlates with almost every mechanical quantity
the elastic network produces — so a raw correlation is worth very little.</p>

<p>Burial alone reaches &rho; = 0.37 on 7WLT. With burial held fixed, the
perturbation response <i>at the gate</i> keeps <b>&rho; = 0.29</b>, against a
null built by circularly shifting the constraint track rather than shuffling it
(both series are strongly autocorrelated along the chain, and a permutation null
is measured to be three times too narrow). Five of eight mechanical features
survive that null, the multiple-comparison correction and the burial control.
The signs are the census's picture in mechanics: coupled to the gate means
constrained, mobile means free. The effect is modest, and that is its honest
size.</p>
"""
