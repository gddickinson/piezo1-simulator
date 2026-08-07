"""The curated ligand table — authored content, validated on the way out.

Split from ``build_ligands.py`` so the data can be read and diffed without the
validation machinery, exactly as ``parameter_table.py`` is split from
``build_parameters.py``.

**Every binding site here is inferred.** No PIEZO structure with a bound
small-molecule modulator has been deposited, so ``site_evidence`` can never be
``bound_structure``; the build refuses to write if it is. The levels, weakest
last:

* ``mutagenesis`` — residues whose substitution measurably changes the response
* ``docking_md``  — a computational pose, no experimental contact
* ``geometric``   — a cavity this project found, with no ligand evidence at all
* ``none``        — no site proposed
"""

from __future__ import annotations

__all__ = ["LIGANDS", "SITE_EVIDENCE"]

SITE_EVIDENCE = {
    "mutagenesis": "residues whose substitution measurably changes the response",
    "docking_md": "a computational pose; no experimental contact",
    "geometric": "a cavity found from geometry alone, with no ligand evidence",
    "none": "no binding site proposed",
}

LIGANDS = [
    dict(
        key="yoda1", name="Yoda1", role="activator", kind="small_molecule",
        pubchem_cid=2746822, inchikey="BQNXBSYSQXSXPT-UHFFFAOYSA-N",
        chembl_id="CHEMBL4303374",
        description="The first PIEZO1 chemical activator; a gating modifier "
                    "that shifts the tension-response curve rather than opening "
                    "the channel by itself.",
        potency=dict(value=26.6, unit="uM", measure="EC50", species="human",
                     assay="Ca2+ influx, HEK293T",
                     citation="syeda2015"),
        site_residues=(1718, 2075, 2078), site_numbering="human",
        site_evidence="docking_md",
        site_note="A pocket between the beam and the repeat-A region, proposed "
                  "from MD simulation and supported by the fenestration it "
                  "occupies. No bound structure exists.",
        site_citation="botellosmith2019",
    ),
    dict(
        key="yoda2", name="Yoda2", role="activator", kind="small_molecule",
        pubchem_cid=170908031, inchikey="OQEIWUHZUPFUQU-UHFFFAOYSA-M",
        chembl_id=None,
        description="A more potent and more soluble Yoda1 analogue.",
        potency=None,
        site_residues=(), site_numbering="human", site_evidence="none",
        site_note="Assumed to share the Yoda1 pocket by analogy, which is an "
                  "assumption rather than a measurement; no site is recorded.",
        site_citation=None,
    ),
    dict(
        key="dooku1", name="Dooku1", role="antagonist", kind="small_molecule",
        pubchem_cid=137321150, inchikey="MNPOBXLPCWFONX-UHFFFAOYSA-N",
        chembl_id=None,
        description="A Yoda1 analogue that antagonises Yoda1 without blocking "
                    "mechanical activation — evidence that the Yoda1 site is "
                    "distinct from the pore.",
        potency=dict(value=1.3, unit="uM", measure="IC50_vs_Yoda1",
                     species="human", assay="Ca2+ influx, HEK293T",
                     citation="evans2018dooku"),
        site_residues=(), site_numbering="human", site_evidence="none",
        site_note="Competes with Yoda1, which implies a shared site but does "
                  "not locate it. No residues are claimed.",
        site_citation="evans2018dooku",
    ),
    dict(
        key="jedi1", name="Jedi1", role="activator", kind="small_molecule",
        pubchem_cid=736516, inchikey="VLMNACSEESRUAK-UHFFFAOYSA-N",
        chembl_id=None,
        description="A hydrophilic activator acting on the distal blade rather "
                    "than the Yoda1 pocket.",
        potency=None,
        site_residues=(), site_numbering="human", site_evidence="mutagenesis",
        site_note="Acts through the peripheral blade and the beam: the lever "
                  "mechanism was mapped by mutation, but no residue set is "
                  "specific enough to record as a binding site.",
        site_citation="wang2018jedi",
    ),
    dict(
        key="jedi2", name="Jedi2", role="activator", kind="small_molecule",
        pubchem_cid=2796026, inchikey="YXDCRSXNEOKXDE-UHFFFAOYSA-N",
        chembl_id=None,
        description="Jedi1 analogue, same proposed mechanism.",
        potency=None,
        site_residues=(), site_numbering="human", site_evidence="mutagenesis",
        site_note="As Jedi1.", site_citation="wang2018jedi",
    ),
    dict(
        key="gsmtx4", name="GsMTx4", role="inhibitor", kind="peptide",
        pubchem_cid=None, inchikey=None, chembl_id=None,
        uniprot="Q7YT39",
        description="A 34-residue tarantula peptide that inhibits cationic "
                    "mechanosensitive channels. It partitions into the outer "
                    "leaflet and acts on the bilayer rather than by plugging "
                    "the pore, so it is not selective for PIEZO1.",
        potency=dict(value=155.0, unit="nM", measure="Kd", species="n/a",
                     assay="bilayer partitioning", citation="bae2011"),
        site_residues=(), site_numbering="human", site_evidence="none",
        site_note="Acts on the lipid bilayer, not on a protein site. Recording "
                  "residues would misrepresent the mechanism.",
        site_citation="bae2011",
    ),
]
