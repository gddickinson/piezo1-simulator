"""The bibliography seed — every citation the code and docs rely on.

Split out of ``build_references.py`` to keep both files under the project's
500-line limit, and along the same seam ``parameter_table.py`` uses: this
module is **data**, so the whole bibliography can be read and diffed without
the resolution machinery in the way.

Each entry is ``(key, query, topic, expect)``. ``expect`` is a word that must
appear in the resolved title or abstract, and it exists because several PMIDs
entered from memory resolved cleanly to entirely unrelated papers. It has since
caught two more things: a PMID for the Kyte-Doolittle paper that returned a 2019
cholesterol-trafficking study, and two of this project's own ``expect`` words
that were subtly wrong ("hydropathy" for a title that says "hydropathic",
"helix" for one that says "helical").
"""

from __future__ import annotations

__all__ = ["SEED"]

SEED = [
    # --- structure -------------------------------------------------------
    ("ge2015", "PMID:26390154", "First Piezo1 cryo-EM structure (3JAC)", "piezo1"),
    ("kamajaya2014", "PMID:25242456", "First Piezo structures: C. elegans CED (4PKE/4PKX)", "piezo"),
    # The two PDB IDs below were swapped until Round 84. Both PMIDs were
    # right, so nothing resolved wrongly and no test could see it; the entries
    # themselves say which is which — 6B3R cites Elife 2017 'Structure-based
    # membrane dome mechanism' (PMID 29231809) and 6BPZ cites Nature 2018
    # 'Structure of the mechanically activated ion channel' (PMID 29261642),
    # matching the Data availability section of each paper.
    ("guo2017", "PMID:29231809", "The dome model; PDB 6B3R", "piezo"),
    ("saotome2018", "PMID:29261642", "Piezo1 structure 6BPZ", "piezo1"),
    ("zhao2018", "PMID:29469092", "Lever-like transduction; PDB 5Z10", "piezo1"),
    ("wang2019piezo2", "PMID:31435011", "Piezo2 structure 6KG7", "piezo2"),
    ("geng2020", "PMID:32142647", "Plug-and-latch gating; Piezo1.1 isoform 6LQI", "piezo"),
    ("yang2022", "PMID:35388220", "Curved and flattened mPIEZO1 in bilayer; 7WLT/7WLU", "piezo1"),
    ("zhou2023mdfic", "PMID:37590348", "MDFIC is a PIEZO auxiliary subunit", "piezo"),
    ("liu2025", "PMID:39719701", "Intermediate-open S2472E; PDB 8IXN/8IXO; lateral portals carry the current", "piezo1"),
    ("vaisey2026", "PMID:42234740", "Lipid cofactor required; force alone insufficient", "piezo1"),
    # --- membrane mechanics ----------------------------------------------
    ("haselwandter2018", "PMID:30480546", "Membrane footprint; R_c 10.2 nm, lambda 14 nm", "piezo"),
    ("haselwandter2022", "PMID:36166476", "Elastic properties and shape of the Piezo dome", "piezo"),
    ("haselwandter2022b", "PMID:36166475", "Quantitative prediction and measurement of the membrane footprint", "piezo"),
    ("dixit2025", "DOI:10.7554/eLife.105138.3", "Nanodome excess area and elasticity", "piezo"),
    ("chong2021", "PMID:33582137", "Full-length model; dome depth 6-7 nm", "piezo1"),
    ("devecchis2021", "PMID:33582135", "MD of opening by membrane tension", "piezo1"),
    # --- electrophysiology and kinetics ----------------------------------
    ("coste2015", "PMID:26008989", "Pore properties dictated by the C-terminal region", "piezo1"),
    ("bae2013", "PMID:23487776", "Xerocytosis mutations slow inactivation", "piezo1"),
    ("lewis2015", "PMID:26646186", "Tension sensitivity; T50 2.7 and 4.7 mN/m", "piezo1"),
    ("cox2016", "PMID:26785635", "Bilayer tension gating; dG0 9.7 kT, dA 8 nm2", "piezo1"),
    ("young2023", "PMID:36795747", "Four-state TENSION model - the one implemented", "mechanotransduction"),
    ("lewis2021", "PMID:34711306", "Clustering does not alter gating", "piezo1"),
    # --- lipids and pharmacology -----------------------------------------
    ("borbiro2015", "PMID:25670203", "Phosphoinositide dependence", "piezo"),
    ("ridone2020", "PMID:32582958", "Cholesterol; P50 shift on depletion", "piezo1"),
    ("shi2020", "PMID:33027663", "Sphingomyelinase disables inactivation", "piezo1"),
    ("romero2019", "PMID:30867417", "Dietary fatty acids tune the mechanical response", "piezo1"),
    ("romero2020", "PMID:32561714", "Margaric acid and PIEZO2", "mechanical"),
    ("buyan2020", "PMID:32949489", "PIP2 and cholesterol sites; K2166-K2169", "piezo1"),
    ("buyan2023", "PMID:35927961", "Lipid redistribution in the curved footprint", "piezo1"),
    ("hashad2025", "PMID:41433068", "PIP2 corrects an endothelial channelopathy", "piezo1"),
    ("syeda2015", "PMID:26001275", "Yoda1", "mechanotransduction"),
    ("botellosmith2019", "PMID:31582801", "Yoda1 mechanism and binding pocket", "piezo1"),
    ("wang2018jedi", "PMID:29610524", "Jedi1/2 and the lever transduction pathway", "piezo1"),
    ("evans2018dooku", "PMID:29498036", "Dooku1 antagonises Yoda1", "yoda1"),
    ("bae2011", "PMID:21696149", "GsMTx4 inhibits Piezo1", "piezo1"),
    ("poole2014", "PMID:24662763", "STOML3 tunes the displacement threshold", "piezo"),
    ("qi2015", "PMID:26443885", "STOML3 membrane stiffening", "stoml3"),
    ("wetzel2017", "PMID:27941788", "STOML3 inhibitors reverse mechanical hypersensitivity", "stoml3"),
    # --- genetics and disease --------------------------------------------
    ("zarychanski2012", "PMID:22529292", "PIEZO1 mutations cause hereditary xerocytosis", "piezo1"),
    ("albuisson2013", "PMID:23695678", "DHS1 mutations; M2225R, R2456H", "piezo1"),
    ("andolfo2013", "PMID:23479567", "DHS1; multiple clinical forms", "piezo1"),
    ("fotiou2015", "PMID:26333996", "PIEZO1 loss of function causes lymphatic dysplasia", "piezo1"),
    ("ma2018malaria", "PMID:29576450", "E756del and malaria resistance", "piezo1"),
    ("karamaticcrew2023", "PMID:36723926", "PIEZO1 carries the Er blood group antigens", "er"),
    # --- methods ----------------------------------------------------------
    ("atilgan2001", "PMID:11159421", "Anisotropic network model", "elastic network"),
    ("bahar2010", "PMID:19785456", "Normal mode analysis of membrane proteins - review", "normal mode"),
    ("smart1996hole", "PMID:9195488", "HOLE pore-radius algorithm", "pore"),
    ("rao2019heuristic", "PMID:31235590",
     "Hydrophobic-gating heuristic: the (hydrophobicity, radius) landscape",
     "hydrophobic"),
    ("klesse2019chap", "PMID:31220459",
     "CHAP - pore annotation; source of the MIT-licensed grid", "chap"),
    ("beckstein2003", "PMID:12740433",
     "Liquid-vapour oscillations of water in hydrophobic nanopores", "water"),
    ("aryal2015", "PMID:25106689", "Hydrophobic gating in ion channels - review",
     "hydrophobic"),
    ("labesse1997", "PMID:9183534", "P-SEA: secondary structure from CA geometry", "secondary structure"),
    ("leguilloux2009", "PMID:19486540",
     "fpocket: alpha-sphere pocket detection - source of the 3.0-5.5 A radii",
     "pocket"),
    ("shrake1973", "PMID:4760134",
     "Shrake-Rupley numerical SASA; source of the 1.4 A water probe", "solvent"),
    ("kabsch1976", "DOI:10.1107/S0567739476001873", "Optimal rotation superposition", "rotation"),
    ("jumper2021", "PMID:34265844", "AlphaFold", "protein structure prediction"),
    ("varadi2024", "PMID:37933859", "AlphaFold DB", "alphafold"),
    # --- sequence homology, Round 89 ---------------------------------------
    # The family comparison needs the conventions it is scored against to be
    # cited rather than remembered, because the distant PIEZOs land exactly on
    # the boundary these papers define.
    ("henikoff1992", "PMID:1438297",
     "BLOSUM62 - the substitution matrix every alignment here scores with",
     "amino acid substitution"),
    ("rost1999", "PMID:10195279",
     "The twilight zone: below ~30% identity, a pairwise alignment stops "
     "distinguishing homology from chance", "twilight"),
    ("smith1981", "DOI:10.1016/0022-2836(81)90087-5",
     "Smith-Waterman local alignment - the statistic that survives the "
     "twilight zone where percent identity does not", "molecular subsequences"),
    # --- external variant predictors, reached through ProtVar (CC BY 4.0) ---
    # Free-text queries rather than remembered PMIDs: six citations in Round 8
    # were resolved from memory to entirely unrelated papers, so every entry
    # here is looked up by title and gated on `expect`.
    ("stephenson2024protvar",
     'TITLE:"ProtVar: mapping and contextualizing human missense variation"',
     "ProtVar API - serves the predictors below under CC BY 4.0", "protvar"),
    ("cheng2023alphamissense",
     'TITLE:"Accurate proteome-wide missense variant effect prediction with AlphaMissense"',
     "AlphaMissense pathogenicity", "alphamissense"),
    ("frazer2021eve",
     'TITLE:"Disease variant prediction with deep generative models of evolutionary data"',
     "EVE - unsupervised variant effect from evolutionary data", "variant"),
    ("brandes2023esm1b",
     'TITLE:"Genome-wide prediction of disease variant effects with a deep protein language model"',
     "ESM-1b variant effects", "language model"),
    ("schymkowitz2005foldx",
     'TITLE:"The FoldX web server: an online force field"',
     "FoldX force field - source of the precomputed ddG", "foldx"),
    # --- HaloTag labelling, imported with the kinetics in Round 32 ----------
    # Looked up by title rather than by remembered PMID, for the reason stated
    # above the ProtVar block.
    ("los2008halotag",
     'TITLE:"HaloTag: a novel protein labeling technology for cell imaging and protein analysis"',
     "HaloTag chemistry; source of the covalent on-rate and its irreversibility",
     "halotag"),
    ("grimm2015jf",
     'TITLE:"A general method to improve fluorophores for live-cell and single-molecule microscopy"',
     "Janelia Fluor dyes including JF646 - cell-permeable, so partition ~ 1",
     "dye"),
    # --- ion permeation, Round 33 ------------------------------------------
    ("coste2010piezo",
     'TITLE:"Piezo1 and Piezo2 are essential components of distinct mechanically activated cation channels"',
     "The original PIEZO1 characterisation; single-channel conductance and "
     "cation non-selectivity", "piezo"),
    ("gnanasambandam2015",
     'TITLE:"Ionic Selectivity and Permeation Properties of Human PIEZO1 Channels"',
     "PIEZO1 selectivity and permeation - the direct target of the PNP model",
     "selectivity"),
    ("hall1975access",
     'TITLE:"Access resistance of a small circular pore"',
     "Access resistance of a circular pore mouth - the term that limits a short "
     "wide pore", "access resistance"),
    # --- population constraint, Round 41 ------------------------------------
    ("chen2024gnomad",
     'TITLE:"A genomic mutational constraint map using variation in 76,156 human genomes"',
     "gnomAD constraint: LOEUF, pLI and missense z-scores", "constraint"),
    # --- calcium nanodomain, Round 35 ---------------------------------------
    ("stern1992",
     'TITLE:"Buffering of calcium in the vicinity of a channel pore"',
     "The steady-state buffered-diffusion Green's function this model uses",
     "calcium"),
    ("naraghi1997",
     'TITLE:"Linearized buffered Ca2+ diffusion in microdomains and its implications for calculation of [Ca2+] at the mouth of a calcium channel"',
     "Linearised buffered diffusion; the screening length and its validity",
     "calcium"),
    ("allbritton1992",
     'TITLE:"Range of messenger action of calcium ion and inositol 1,4,5-trisphosphate"',
     "Cytosolic calcium diffusion coefficient and buffering range", "calcium"),
    ("tsien1980bapta",
     'TITLE:"New calcium indicators and buffers with high selectivity against magnesium and protons: design, synthesis, and properties of prototype structures"',
     "BAPTA - the chelator scaffold of the JF646-BAPTA sensor, and its Kd",
     "calcium"),
    ("bertaccini2025piezo1",
     'TITLE:"Visualizing PIEZO1 Localization and Activity in hiPSC-Derived Single Cells and Organoids with HaloTag Technology"',
     "The tagged-PIEZO1 experiment this labelling model describes; three tags "
     "per channel and the multi-level brightness histogram",
     "halotag"),
    # --- replicating Guo & MacKinnon 2017's figures, Round 84 ---------------
    # Every one of these is a method or a number the paper's own figures rest
    # on and this project had no citable source for. Resolved by title: the
    # PMID form of kyte1982 came back as a 2019 cholesterol-trafficking paper,
    # which is exactly the failure the ``expect`` gate exists to catch.
    ("kyte1982",
     'TITLE:"A simple method for displaying the hydropathic character of a protein"',
     "The Kyte-Doolittle hydropathy scale and sliding window, the method behind "
     "Figure 3-figure supplements 1-3", "hydropathic"),
    ("vonheijne1992",
     'TITLE:"Membrane protein structure prediction. Hydrophobicity analysis and the positive-inside rule"',
     "The positive-inside rule the paper checks its charge distribution "
     "against in Figure 4c", "topology"),
    ("rawicz2000",
     'TITLE:"Effect of chain length and unsaturation on elasticity of lipid bilayers"',
     "Lytic tension of a bilayer, ~3.5 k_BT/nm^2 - the scale the paper's 42 "
     "k_BT stabilisation is quoted against", "membrane"),
    ("dolinsky2004",
     'TITLE:"PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics calculations"',
     "APBS/PDB2PQR, the Poisson-Boltzmann route Figure 4c was computed with "
     "and the one our screened-Coulomb surface deliberately is not",
     "electrostatics"),
    ("pauling1951",
     # PubMed renders this title with a semicolon, not the colon the paper
     # itself uses, so an exact-phrase match on the full title fails. The
     # distinctive tail matches either way.
     'TITLE:"two hydrogen-bonded helical configurations of the polypeptide chain"',
     "The alpha helix: 1.5 A rise and 100 degrees of turn per residue, the "
     "geometry the cross-helix detector is calibrated against", "helical"),
    ("perozo2002",
     'TITLE:"Open channel structure of MscL and the gating mechanism of mechanosensitive channels"',
     "MscL's ~20 nm^2 in-plane expansion, the comparison the dome mechanism is "
     "argued against", "mscl"),
]


#: Entries Europe PMC cannot resolve, typed by hand. Kept deliberately small —
#: anything here is a citation whose metadata is not machine-verified.
MANUAL = {
    "kabsch1976": {
        "title": "A solution for the best rotation to relate two sets of vectors",
        "authors": "Kabsch W.",
        "journal": "Acta Crystallographica Section A",
        "year": "1976", "volume": "32", "pages": "922-923",
        "pmid": None, "pmcid": None,
        "doi": "10.1107/S0567739476001873",
        "is_open_access": False, "has_pdf": False,
        "url": "https://doi.org/10.1107/S0567739476001873",
        "note": "Predates PubMed indexing; metadata entered by hand.",
    },
}
