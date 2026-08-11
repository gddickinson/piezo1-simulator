"""Parameters for replicating Guo & MacKinnon 2017's figures.

Split out of ``parameter_table.py`` to keep both under the project's 500-line
limit, and along a real seam: every number here exists because a *published
figure* depends on it. Some are the paper's own idealisations (the 10.2 nm
sphere's centre height, the 3.6 nm bilayer); some are the methods its panels
were computed with (the Kyte-Doolittle window, the 150 mM NaCl of the APBS
surface); some are choices this project had to make to reproduce a panel the
paper describes only in words (which helices count as "cross" helices).

The third group is the one to be careful with. A number invented to reproduce a
picture is a number tuned to an answer, so each of those says in its
``source_note`` what the picture would look like if it were wrong, and each is
probed by :mod:`piezo1.analysis.parameter_effect` like any other.
"""

from __future__ import annotations

__all__ = ["FIGURE_PARAMETERS"]

FIGURE_PARAMETERS = [
    # ------------------------------------- the idealised dome (Figure 7) ----
    dict(key="membrane.thickness", name="Bilayer thickness", value=3.6,
         unit="nm", minimum=2.0, maximum=6.0, kind="physical",
         category="Membrane mechanics", citation="guo2017",
         source_note="Figure 7a: 'a semi-sphere-shaped membrane 3.6 nm thick', "
                     "idealised from the curved micelle density. Sets only "
                     "where the two leaflet surfaces are drawn relative to the "
                     "mid-plane; no area or energy in the dome model depends "
                     "on it, because Helfrich theory works on the mid-plane.",
         description="Thickness of the bilayer the dome carries."),
    dict(key="dome.idealised_center_height",
         name="Idealised dome centre height", value=4.0, unit="nm",
         minimum=0.0, maximum=15.0, kind="empirical",
         category="Dome energetics", citation="guo2017",
         source_note="Figure 7a: the mid-plane semi-sphere 'is centered 4.0 nm "
                     "above the projection plane'. With the 10.2 nm radius it "
                     "fixes every other number in Figure 7-figure supplement "
                     "1 — an 18.8 nm opening, a 6.2 nm depth, 397 nm^2 of "
                     "surface and 277 nm^2 projected. Raising it flattens the "
                     "cap and shrinks the released area towards zero.",
         description="Height of the idealised sphere's centre above the plane."),
    dict(key="dome.lytic_tension", name="Bilayer lytic tension", value=3.5,
         unit="k_BT/nm^2", minimum=0.5, maximum=10.0, kind="physical",
         category="Dome energetics", citation="rawicz2000",
         source_note="The scale Guo & MacKinnon quote their 42 k_BT "
                     "stabilisation against: one tenth of lytic tension acting "
                     "through 120 nm^2. Rawicz et al. 2000 measure rupture "
                     "tensions of 5-10 mN/m across chain lengths; 3.5 "
                     "k_BT/nm^2 is 14.4 mN/m, at the top of that range, so the "
                     "42 k_BT is an optimistic figure and is reported as one.",
         description="Tension at which a bilayer ruptures."),

    # --------------------------------- hydropathy (Figure 3-S1 to S3) ------
    dict(key="hydropathy.window", name="Hydropathy window", value=19,
         unit="residues", minimum=5, maximum=31, kind="method",
         category="Hydropathy", citation="kyte1982",
         source_note="Kyte & Doolittle recommend 19-21 for identifying "
                     "membrane-spanning segments and 7-11 for surface "
                     "exposure. 19 is the shorter of the two membrane values, "
                     "so it resolves the 4-TM repeat rather than smoothing it "
                     "into one broad hydrophobic block.",
         description="Sliding-window width for the hydropathy average."),
    dict(key="hydropathy.tm_threshold", name="Transmembrane call threshold",
         value=1.6, unit="hydropathy", minimum=0.0, maximum=4.5, kind="method",
         category="Hydropathy", citation="kyte1982",
         source_note="Kyte & Doolittle's own figures place the membrane-"
                     "spanning cut at a window average near +1.6 on the same "
                     "scale. It is a threshold on a smoothed curve, not a "
                     "prediction method: raising it loses genuine TM segments "
                     "before it removes false ones, which is measured against "
                     "UniProt's 38 rather than asserted.",
         description="Window hydropathy above which a segment is called TM."),
    dict(key="hydropathy.min_tm_length", name="Shortest TM call", value=15,
         unit="residues", minimum=5, maximum=40, kind="method",
         category="Hydropathy", citation="method_choice",
         source_note="a run shorter than this cannot cross a 3.6 nm bilayer as "
                     "an alpha helix (1.5 A rise per residue gives 22.5 A at "
                     "15 residues, against a 36 A bilayer, so this is already "
                     "generous and admits tilted and interfacial helices)",
         description="Minimum run length above threshold to call a segment."),

    # --------------------------------- surface electrostatics (Figure 4c) --
    dict(key="electrostatics.ionic_strength", name="Ionic strength", value=0.15,
         unit="M", minimum=0.0, maximum=1.0, kind="convention",
         category="Electrostatics", citation="guo2017",
         source_note="Figure 4c states 'aqueous solution containing 150 mM "
                     "NaCl'. It enters only through the Debye length (7.86 A "
                     "at this value), so it sets how quickly the surface "
                     "potential decays, not its sign anywhere.",
         description="Bulk 1:1 salt concentration screening the surface."),
    dict(key="electrostatics.dielectric_solvent", name="Solvent permittivity",
         value=78.5, unit="", minimum=1.0, maximum=100.0, kind="physical",
         category="Electrostatics", citation="dolinsky2004",
         source_note="water at 298 K, the value APBS uses by default and the "
                     "one Figure 4c will have been computed with",
         description="Relative permittivity of the aqueous phase."),
    dict(key="electrostatics.temperature", name="Electrostatics temperature",
         value=298.15, unit="K", minimum=270.0, maximum=320.0,
         kind="convention", category="Electrostatics", citation="convention",
         source_note="25 C, the default APBS calculations are reported at. The "
                     "paper does not state a temperature; the Debye length "
                     "moves by 2% between 25 and 37 C, which is invisible on a "
                     "coloured surface.",
         description="Temperature entering the Debye length and the k_BT unit."),
    dict(key="electrostatics.colour_scale", name="Surface potential scale",
         value=5.0, unit="k_BT/e", minimum=0.5, maximum=50.0,
         kind="convention", category="Electrostatics", citation="guo2017",
         source_note="Figure 4c's colour bar runs -5 to +5 k_BT/e. Reproducing "
                     "the panel means reproducing its saturation points too: a "
                     "wider scale makes any surface look neutral and a "
                     "narrower one makes any surface look charged.",
         description="Potential at which the red/blue surface colouring saturates."),
    dict(key="electrostatics.max_distance", name="Charge interaction cutoff",
         value=30.0, unit="A", minimum=5.0, maximum=200.0, kind="method",
         category="Electrostatics", citation="method_choice",
         source_note="screened Coulomb falls as exp(-r/lambda_D)/r, so a "
                     "cutoff is what makes the sum over 6B3R's 804 formal "
                     "charges and 135k surface points tractable. The "
                     "truncation error is measured rather than argued: "
                     "against the exact all-pairs sum on 6B3R's own surface, "
                     "30 A costs at most 0.089 k_BT/e (2.2% of the peak) and "
                     "no point moves by more than 0.1, while 20 A costs 0.30 "
                     "(7.5%) and moves 19% of points that far. A first draft "
                     "of this note claimed 'under 0.2%' from the analytic "
                     "falloff alone and was an order of magnitude out, "
                     "because the falloff bounds one charge and not 804.",
         description="Range beyond which a point charge is ignored."),

    # --------------------------- simulated projections (Figure 2a and 2b) --
    dict(key="projection.pixel_size", name="Projection pixel size", value=1.3,
         unit="A", minimum=0.2, maximum=10.0, kind="convention",
         category="Projections", citation="guo2017",
         source_note="the calibrated physical pixel size of the K2 data in the "
                     "Methods, so a simulated projection is sampled the way "
                     "the real 2D class averages were",
         description="Sampling of the simulated projection image."),
    dict(key="projection.resolution", name="Projection resolution", value=8.0,
         unit="A", minimum=1.0, maximum=40.0, kind="method",
         category="Projections", citation="method_choice",
         source_note="a 2D class average is not at the map's global "
                     "resolution; 8 A reproduces the blurred, secondary-"
                     "structure-free appearance of Figure 2a,b. It is a "
                     "display choice and the panel is labelled a simulation "
                     "rather than a class average because of it.",
         description="Gaussian blur applied to the simulated projection."),

    # ------------------------------- cross-helices (Figure 7b, in yellow) --
    dict(key="cross_helix.min_tilt_deg", name="Cross-helix minimum tilt",
         value=55.0, unit="degrees", minimum=20.0, maximum=90.0, kind="method",
         category="Architecture", citation="method_choice",
         source_note="Guo & MacKinnon describe the cross-helices as running "
                     "'perpendicular to the TM helices' but give no residue "
                     "ranges anywhere in the paper, so this reproduces their "
                     "yellow helices by the property they are named for. 55 "
                     "degrees from the local membrane normal separates them "
                     "from transmembrane helices, which in this structure tilt "
                     "up to about 40; the count is reported against the "
                     "threshold so a reader can see how sharp the separation "
                     "is rather than trusting the cut.",
         description="Tilt from the membrane normal above which a linker helix "
                     "is called a cross-helix."),
    dict(key="micelle.offset", name="Micelle shell thickness", value=9.0,
         unit="A", minimum=2.0, maximum=25.0, kind="empirical",
         category="Micelle", citation="guo2017",
         source_note="How far outside the hydrophobic belt the modelled "
                     "detergent envelope is drawn. Chosen so the envelope's "
                     "overall width matches the micelle in Figure 4b rather "
                     "than derived from digitonin's dimensions, which the "
                     "paper does not give — so it sets the *thickness* of the "
                     "drawn shell and nothing else. The curvature reported "
                     "with it is fitted to the belt atoms themselves and does "
                     "not depend on this value at all, which is the point: an "
                     "offset surface has the same centre and a radius larger "
                     "by exactly the offset.",
         description="Offset of the modelled micelle surface from the belt."),
    dict(key="micelle.grid_spacing", name="Micelle grid spacing", value=1.6,
         unit="A", minimum=0.5, maximum=6.0, kind="method",
         category="Micelle", citation="method_choice",
         source_note="marching-cubes resolution. Costs cubically and changes "
                     "how faceted the surface looks, not where it is; the "
                     "enclosed volume is stable to under a percent between "
                     "1.2 and 2.4 A, which is measured rather than assumed.",
         description="Grid step for contouring the micelle envelope."),
    dict(key="helix.rise", name="Alpha-helix rise per residue", value=1.5,
         unit="A", minimum=0.5, maximum=4.0, kind="physical",
         category="Architecture", citation="pauling1951",
         source_note="C-alpha advance along the helix axis per residue. With "
                     "the 100-degree turn this gives 3.6 residues and 5.4 A "
                     "per turn, the standard alpha helix. The detector is "
                     "calibrated against an analytically generated helix at "
                     "exactly these values before it is used on coordinates.",
         description="Axial rise per residue of an alpha helix."),
    dict(key="helix.radius", name="Alpha-helix C-alpha radius", value=2.3,
         unit="A", minimum=0.5, maximum=5.0, kind="physical",
         category="Architecture", citation="pauling1951",
         source_note="distance of a C-alpha from the helix axis; with the "
                     "rise it separates an alpha helix from a 3-10 (1.9 A) "
                     "and a pi helix (2.8 A)",
         description="Radial offset of C-alpha from the alpha-helix axis."),
    dict(key="helix.turn", name="Alpha-helix turn per residue", value=100.0,
         unit="degrees", minimum=45.0, maximum=180.0, kind="physical",
         category="Architecture", citation="pauling1951",
         source_note="rotation about the axis between consecutive residues, "
                     "right-handed. This is the criterion that excludes a "
                     "random coil: rise and radius alone passed 41% of the "
                     "windows of a synthetic random walk, because a walk with "
                     "a fixed step length looks locally like a helix on both.",
         description="Rotation per residue about the alpha-helix axis."),
    dict(key="helix.rise_tolerance", name="Helix rise tolerance", value=0.45,
         unit="A", minimum=0.05, maximum=1.5, kind="method",
         category="Architecture", citation="method_choice",
         source_note="wide enough for a real, bent, 3.7-A-resolution helix and "
                     "narrow enough that a 3-10 helix (2.0 A rise) falls "
                     "outside it; the measured margin on 3-10 is 0.05 A, so "
                     "this is the tightest of the three criteria",
         description="Permitted departure from the ideal helix rise."),
    dict(key="helix.radius_tolerance", name="Helix radius tolerance", value=1.1,
         unit="A", minimum=0.1, maximum=3.0, kind="method",
         category="Architecture", citation="method_choice",
         source_note="deliberately loose: at 3.7 A resolution a C-alpha "
                     "position carries most of an Angstrom of uncertainty, and "
                     "the turn criterion is what does the discriminating",
         description="Permitted departure from the ideal C-alpha radius."),
    dict(key="helix.turn_tolerance", name="Helix turn tolerance", value=25.0,
         unit="degrees", minimum=5.0, maximum=90.0, kind="method",
         category="Architecture", citation="method_choice",
         source_note="applied to the *worst* step in a window rather than the "
                     "mean, so one reversal cannot hide inside an average. At "
                     "25 degrees a 3-10 helix (120) and a pi helix (87) are "
                     "both excluded and every window of an ideal alpha helix "
                     "is admitted.",
         description="Permitted departure from the ideal turn, worst step."),
    dict(key="helix.window", name="Helix detection window", value=7,
         unit="residues", minimum=5, maximum=21, kind="method",
         category="Architecture", citation="method_choice",
         source_note="just under two turns, the shortest window in which the "
                     "turn criterion has enough steps to mean anything. It "
                     "slides one residue at a time, so a bent helix is still "
                     "recognised along its length rather than failing as a "
                     "whole.",
         description="Sliding-window width for the helix criteria."),
    dict(key="electrostatics.histidine_charge", name="Histidine formal charge",
         value=0.0, unit="e", minimum=0.0, maximum=1.0, kind="convention",
         category="Electrostatics", citation="convention",
         source_note="At pH 7.4 a typical histidine (pKa ~6.0) is about 10% "
                     "protonated. Zero is used rather than 0.1 because a "
                     "fractional formal charge is not a state any single "
                     "residue is in, and assigning a full +1 would put three "
                     "spurious positive charges per protomer on the surface "
                     "Figure 4c is read from. `compare_conventions` reports "
                     "the effect of both alternatives rather than hiding the "
                     "choice.",
         description="Charge placed on histidine when building formal charges."),
    dict(key="architecture.proximal_first_helix",
         name="First pore-proximal helix", value=25, unit="helix index",
         minimum=1, maximum=38, kind="method", category="Architecture",
         citation="measured_here",
         source_note="TM25 opens THU7, the most distal 4-TM unit every "
                     "deposited PIEZO1 entry in this project's catalogue "
                     "resolves in all three protomers. Measured from the "
                     "downloaded entries rather than chosen: it is the split "
                     "that makes two entries comparable, and 6B3R resolves "
                     "TM13-38 where 6BPZ resolves only TM25-38.",
         description="Helix index dividing the distal blade from the pore "
                     "module for coverage-matched comparisons."),
    dict(key="cross_helix.min_length", name="Cross-helix minimum length",
         value=7, unit="residues", minimum=4, maximum=30, kind="method",
         category="Architecture", citation="method_choice",
         source_note="two turns of alpha helix; shorter runs are single turns "
                     "and loop fragments whose axis direction is not "
                     "meaningful, and admitting them adds orientation noise "
                     "rather than helices",
         description="Shortest helical run admitted as a cross-helix."),
]
