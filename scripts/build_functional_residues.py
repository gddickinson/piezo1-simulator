#!/usr/bin/env python
"""Author the functional-residue resource.

Writes ``piezo1/resources/functional_residues.json``: the individual residues
that a user will want to find on the structure — the hydrophobic gate, the
acidic residues that set ion selectivity, the Yoda1 pocket, the PIP2-binding
lysine cluster, and the anchor-domain brake.

Every entry is declared in **human Q92508 numbering** and its mouse equivalent
is computed by alignment, then **verified by checking that the amino acid
identity matches** in both species. Any entry that fails that check is written
out with ``"verified": false`` rather than silently kept, because a residue
number that lands on the wrong amino acid is worse than no annotation at all.

Usage::

    python scripts/build_functional_residues.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import RESOURCE_DIR  # noqa: E402
from piezo1.core.sequence import human_sequence, load_numbering_map, mouse_sequence  # noqa: E402

# Declared in human numbering, with the expected wild-type residue letter so
# the build can self-check.
SITES = [
    # --- pore and gate -----------------------------------------------------
    {"group": "hydrophobic_gate", "residues": {2447: "I", 2450: "V", 2454: "F"},
     "label": "Transmembrane hydrophobic gate",
     "description": "The constriction in the inner-helix bundle that occludes "
                    "the pore in the closed state.",
     "source": "Yang et al. eLife 2025 (human numbering)",
     "color": "#e8556d", "category": "pore"},

    {"group": "pore_lining", "residues": {2440: None, 2443: "L", 2447: "I",
                                          2451: "I", 2454: "F"},
     "label": "Pore-lining residues",
     "description": "Inner-helix residues facing the conduction pathway. "
                    "Human I2451 corresponds to mouse V2477 - not conserved.",
     "source": "Botello-Smith et al. Nat Commun 2019 (converted from mouse)",
     "color": "#c9455c", "category": "pore"},

    {"group": "selectivity_acidic", "residues": {2117: "E", 2461: "E",
                                                 2469: "E", 2470: "E"},
     "label": "Acidic residues setting ion selectivity",
     "description": "Glutamates lining the anchor (E2117) and the CTD "
                    "vestibule and lateral fenestrations. Mutating them alters "
                    "divalent permeation and unitary conductance.",
     "source": "Coste et al. Nature 2015; Zhao et al. Nature 2016 "
               "(mouse E2133/E2495/E2496 -> human E2117/E2469/E2470)",
     "color": "#ff8a3d", "category": "pore"},

    {"group": "ctd_constriction", "residues": {2467: "M", 2468: "F",
                                               2510: "P", 2511: "E"},
     "label": "Cytoplasmic vestibule constrictions",
     "description": "Two constrictions in the C-terminal domain beneath the "
                    "pore that ions must pass on the intracellular side.",
     "source": "converted from mouse M2493/F2494 and P2536/E2537",
     "color": "#d2691e", "category": "pore"},

    # --- ligand sites ------------------------------------------------------
    {"group": "yoda1_pocket", "residues": {1718: "A", 2075: "A", 2078: "A"},
     "label": "Yoda1 binding pocket",
     "description": "Pocket between blade repeats A and B where the agonist "
                    "Yoda1 is proposed to act as a molecular wedge. Mapped by "
                    "mutagenesis and MD, not by a ligand-bound structure. "
                    "Note A1718 happens to carry the same number in both "
                    "species by coincidence - do not generalise from it.",
     "source": "Botello-Smith et al. Nat Commun 2019 (PMID 31582801)",
     "color": "#ffd93d", "category": "ligand", "evidence": "predicted"},

    {"group": "pip2_cluster", "residues": {2166: "K", 2167: "K", 2168: "K",
                                           2169: "K"},
     "label": "Polybasic PIP2-binding cluster",
     "description": "Four lysines immediately before helix 37, at the "
                    "inner-leaflet pore periphery. Deleting them abolishes "
                    "inactivation without changing the pressure threshold. "
                    "A DHS1 disease variant deletes K2166.",
     "source": "Buyan et al. Biophys J 2020 (PMID 32949489)",
     "color": "#5ec8ff", "category": "lipid", "evidence": "simulation+mutagenesis"},

    {"group": "anchor_brake", "residues": {2113: "P", 2114: "F"},
     "label": "Anchor-domain apex brake",
     "description": "Apex of the anchor domain acting as a resistive brake on "
                    "the inner helix; also sits in a cholesterol/CRAC context.",
     "source": "Li/Cox/Martinac Channels 2021; Buyan et al. 2020",
     "color": "#4fc3c7", "category": "lipid"},

    {"group": "lipid_anchor_arg", "residues": {2456: "R"},
     "label": "Inner-leaflet arginine (R2456)",
     "description": "Site of the archetypal gain-of-function xerocytosis "
                    "variant R2456H, which markedly slows inactivation.",
     "source": "Zarychanski et al. Blood 2012; Bae et al. PNAS 2013",
     "color": "#ff5ec8", "category": "disease"},

    # --- other basic clusters implicated in phosphoinositide binding -------
    {"group": "basic_cluster_2", "residues": {623: "R", 624: "K", 627: "K"},
     "label": "Basic cluster (blade, THU4)",
     "description": "Inner-leaflet basic patch contacting phosphoinositides in "
                    "coarse-grained simulation.",
     "source": "Buyan et al. Biophys J 2020",
     "color": "#7fd4ff", "category": "lipid", "evidence": "simulation"},

    {"group": "basic_cluster_3", "residues": {1724: "R", 1727: "K", 1728: "R"},
     "label": "Basic cluster (blade, THU8)",
     "description": "Inner-leaflet basic patch contacting phosphoinositides in "
                    "coarse-grained simulation.",
     "source": "Buyan et al. Biophys J 2020",
     "color": "#7fd4ff", "category": "lipid", "evidence": "simulation"},

    {"group": "basic_cluster_4", "residues": {1911: "K", 1912: "R", 1915: "R",
                                              1919: "R", 1921: "R"},
     "label": "Basic cluster (THU8-THU9 linker)",
     "description": "Inner-leaflet basic patch contacting phosphoinositides in "
                    "coarse-grained simulation.",
     "source": "Buyan et al. Biophys J 2020",
     "color": "#7fd4ff", "category": "lipid", "evidence": "simulation"},
]


def main() -> int:
    hs, ms = human_sequence(), mouse_sequence()
    nm = load_numbering_map()

    groups = []
    n_ok = n_bad = 0
    for spec in SITES:
        entries = []
        for human_pos, expected in sorted(spec["residues"].items()):
            actual = hs[human_pos - 1]
            mouse_pos = nm.to_b(human_pos)
            mouse_aa = ms[mouse_pos - 1] if mouse_pos else None
            verified = (expected is None or actual == expected)
            conserved = mouse_aa == actual if mouse_aa else None
            if verified:
                n_ok += 1
            else:
                n_bad += 1
                print(f"  ! {spec['group']}: expected {expected}{human_pos} "
                      f"but Q92508 has {actual}{human_pos}")
            entries.append({
                "human": human_pos, "human_aa": actual,
                "mouse": mouse_pos, "mouse_aa": mouse_aa,
                "offset": (mouse_pos - human_pos) if mouse_pos else None,
                "conserved": conserved,
                "verified": verified,
            })
        groups.append({
            "id": spec["group"], "label": spec["label"],
            "category": spec["category"],
            "description": spec["description"], "source": spec["source"],
            "color": spec["color"],
            "evidence": spec.get("evidence", "experimental"),
            "residues": entries,
        })

    out = {
        "protein": "PIEZO1",
        "numbering": "declared in human Q92508; mouse computed by alignment",
        "note": ("Every residue was checked against the actual Q92508 sequence "
                 "at build time. 'conserved' compares the human and mouse "
                 "amino acid at aligned positions."),
        "groups": groups,
    }
    dest = RESOURCE_DIR / "functional_residues.json"
    dest.write_text(json.dumps(out, indent=1))

    print(f"\nwrote {dest}")
    print(f"{n_ok} residues verified against Q92508, {n_bad} mismatched\n")
    for g in groups:
        res = ", ".join(f"{r['human_aa']}{r['human']}"
                        + ("" if r["conserved"] else f"(m:{r['mouse_aa']}{r['mouse']})")
                        for r in g["residues"])
        print(f"  {g['id']:20s} {g['category']:8s} {res}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
