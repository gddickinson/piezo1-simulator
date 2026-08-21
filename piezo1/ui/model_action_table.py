"""The structural displays themselves — one row per thing a result can draw.

Split from :mod:`piezo1.ui.model_actions` at the length limit and along the
seam ``parameter_audit_exemptions`` uses: that module is the *mechanism* (find
the control, press it, say what happened) and this is the *judgement* about
which displays are worth offering and which control each one belongs to —
the half a reader checks by eye.

Every ``target`` here is a real thing elsewhere: a catalogued PDB id, a
component key from :mod:`piezo1.structure.components`, a colouring the Model
panel offers, or a curated residue group. A test resolves all of them against
a real window, because a display naming a structure nobody has downloaded or a
component that has been renamed would be a button that appears to work.
"""

from __future__ import annotations

from .model_actions import ModelActionSpec

__all__ = ["MODEL_ACTIONS"]


MODEL_ACTIONS: dict[str, ModelActionSpec] = {
    # --------------------------------------------- overlays already in View
    "pore_surface": ModelActionSpec("Draw the pore surface", menu="&Pore surface"),
    "ion_flux": ModelActionSpec("Run the ion animation",
                                menu="&Ion flux animation"),
    "contacts": ModelActionSpec("Draw the contacts", menu="&Contacts"),
    "pockets": ModelActionSpec("Draw the pockets", menu="P&ockets"),
    "nanodomain": ModelActionSpec("Draw the calcium field",
                                  menu="&Calcium nanodomain"),
    "halotag": ModelActionSpec("Draw the modelled tags",
                               menu="Show modelled &tags"),
    "full_length": ModelActionSpec("Draw the full-length model",
                                   menu="&Full-length model"),
    "dome_surface": ModelActionSpec("Draw the dome surface", menu="&Dome surface"),
    "micelle": ModelActionSpec("Draw the modelled micelle",
                               menu="&Micelle density (modelled)"),
    "planar_membrane": ModelActionSpec("Draw the planar membrane",
                                       menu="Planar &membrane (one protomer)"),
    "allosteric_path": ModelActionSpec("Draw the blade-to-gate route",
                                       menu="&Allosteric path"),
    "colour_constraint": ModelActionSpec(
        "Colour by evolutionary constraint",
        menu="Colour by evolutionary c&onstraint"),
    "colour_electrostatics": ModelActionSpec(
        "Colour by electrostatics", menu="Colour by &electrostatics"),
    "colour_fluctuation": ModelActionSpec(
        "Colour by predicted fluctuation",
        button="physics_panel.fluctuation_button"),

    # ------------------------------------------------- one part, up close
    "component_pore": ModelActionSpec(
        "Show the pore module", kind="component", target="pore_module"),
    "component_gate": ModelActionSpec(
        "Show the transmembrane gate", kind="component", target="tm_gate"),
    "component_cap": ModelActionSpec(
        "Show the cap and its gates", kind="component", target="cap_and_gates"),
    "component_vestibule": ModelActionSpec(
        "Show the cytoplasmic vestibule", kind="component",
        target="ctd_vestibule"),
    "component_anchor": ModelActionSpec(
        "Show the anchor domain", kind="component", target="anchor"),
    "component_blade": ModelActionSpec(
        "Show the blade", kind="component", target="blade"),
    "component_beam": ModelActionSpec(
        "Show the beam and lateral plug", kind="component",
        target="beam_and_latch"),
    "component_md": ModelActionSpec(
        "Show the simulated construct", kind="component",
        target="md_construct"),
    "component_whole": ModelActionSpec(
        "Show the whole assembly again", kind="component", target="whole"),

    # ------------------------------------------------- another structure
    "load_curved": ModelActionSpec("Load the curved state (7WLT)",
                                   kind="load", target="7WLT"),
    "load_intermediate": ModelActionSpec("Load the intermediate state (8IXN)",
                                         kind="load", target="8IXN"),
    "load_open": ModelActionSpec("Load the intermediate-open state (8IXO)",
                                 kind="load", target="8IXO"),
    "load_flat": ModelActionSpec("Load the flattened state (7WLU)",
                                 kind="load", target="7WLU"),
    "load_guo_entry": ModelActionSpec("Load the paper's own entry (6B3R)",
                                      kind="load", target="6B3R"),
    "load_variant": ModelActionSpec("Load the variant entry (8YFG)",
                                    kind="load", target="8YFG"),
    "load_open_like": ModelActionSpec("Load the open-like entry (11ZC)",
                                      kind="load", target="11ZC"),
    "load_piezo2": ModelActionSpec("Load PIEZO2 (6KG7)",
                                   kind="load", target="6KG7"),
    "load_piezo3_model": ModelActionSpec(
        "Load the piezo3 model (AlphaFold)", kind="load",
        target="AF-A0AB32U1Q1-F1-model_v6"),

    # ------------------------------------------------- side by side
    "companion_flat": ModelActionSpec("Draw the flattened state beside it",
                                      kind="companion", target="7WLU"),
    "companion_open": ModelActionSpec("Draw the intermediate-open state beside it",
                                      kind="companion", target="8IXO"),
    "companion_piezo2": ModelActionSpec("Draw PIEZO2 beside it",
                                        kind="companion", target="6KG7"),

    # ------------------------------------------------- superposition
    "overlay_flat": ModelActionSpec("Superpose the flattened state",
                                    kind="overlay", target="7WLU",
                                    mode="protomer"),
    "overlay_core_piezo2": ModelActionSpec(
        "Superpose PIEZO2 on the pore module", kind="overlay", target="6KG7",
        mode="core"),
    "overlay_core_flat": ModelActionSpec(
        "Superpose the flattened state on the pore module", kind="overlay",
        target="7WLU", mode="core"),
    "overlay_core_worm": ModelActionSpec(
        "Superpose PEZO-1 on the pore module", kind="overlay", target="9ZIS",
        mode="core"),
    "overlay_variant": ModelActionSpec(
        "Superpose the variant entry on the pore module", kind="overlay",
        target="8YFG", mode="core"),

    # ------------------------------------------------- residues
    "highlight_gate": ModelActionSpec("Mark the hydrophobic gate",
                                      kind="highlight", target="site:hydrophobic_gate"),
    "highlight_selectivity": ModelActionSpec(
        "Mark the selectivity glutamates", kind="highlight",
        target="site:selectivity_acidic"),
    "highlight_yoda1": ModelActionSpec("Mark the Yoda1 pocket",
                                       kind="highlight", target="site:yoda1_pocket"),
    "highlight_pip2": ModelActionSpec("Mark the PIP2 cluster",
                                      kind="highlight", target="site:pip2_cluster"),
    "highlight_ctd": ModelActionSpec("Mark the cytoplasmic constrictions",
                                     kind="highlight", target="site:ctd_constriction"),
    "highlight_pathogenic": ModelActionSpec(
        "Mark the pathogenic pore positions", kind="highlight",
        target="family:pathogenic_pore"),
    "highlight_equivalent": ModelActionSpec(
        "Mark the two equivalent positions", kind="highlight",
        target="family:equivalent"),
    "highlight_r2456": ModelActionSpec("Mark R2456", kind="highlight",
                                       target="variant:R2456H"),

    # ------------------------------------------------- colour and motion
    "colour_plddt": ModelActionSpec("Colour by AlphaFold pLDDT", kind="colour",
                                    target="AlphaFold pLDDT"),
    "colour_bfactor": ModelActionSpec("Colour by the deposited B-factor",
                                      kind="colour", target="B-factor"),
    "colour_domain": ModelActionSpec("Colour by domain", kind="colour",
                                     target="Domain"),
    "colour_hydrophobicity": ModelActionSpec(
        "Colour by hydrophobicity", kind="colour",
        target="Hydrophobicity (Kyte-Doolittle)"),
    "morph": ModelActionSpec("Build the curved-to-flat morph", kind="morph"),
}
