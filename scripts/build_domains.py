#!/usr/bin/env python
"""Author the PIEZO1 domain-architecture resource.

Writes ``piezo1/resources/domains.json``: every architectural element of PIEZO1
with residue ranges in **both** human (Q92508) and mouse (E2JF22) numbering, a
provenance record, and a confidence label.

Design note on provenance. Domain boundaries come from three different kinds of
source and the file records which:

``uniprot``     Taken directly from the UniProt feature table — transmembrane
                segments and topological domains. High confidence.
``derived``     Computed here from UniProt features by an explicit rule, e.g.
                "the cap is the extracellular topological domain lying between
                TM37 and TM38". High confidence, and the rule is stated.
``literature``  Read from the structural-biology literature, where boundaries
                are often approximate and stated in mouse numbering. Medium
                confidence; the mouse range is authoritative and the human
                range is obtained by alignment.

Mouse ranges are produced with :mod:`piezo1.core.sequence`, never by adding a
constant — the human/mouse offset varies from 0 to +26 along the chain.

Usage::

    python scripts/build_domains.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import RESOURCE_DIR  # noqa: E402
from piezo1.core.sequence import load_numbering_map  # noqa: E402

# --------------------------------------------------------------------------
# Domains defined from the literature, in MOUSE numbering where the source
# used mouse. Human ranges are derived by alignment.
# --------------------------------------------------------------------------

LITERATURE_DOMAINS = [
    {
        "id": "beam",
        "name": "Beam",
        "category": "lever",
        "numbering": "mouse",
        "start": 1300, "end": 1365,
        "description": (
            "Long intracellular helix running from beneath the distal blade to "
            "the pore module. The lever arm of the lever-like transduction "
            "mechanism: blade motion is transmitted along the beam to the "
            "anchor and thence to the inner helix."
        ),
        "source": "Zhao et al. Nature 2018 (PMID 30089899); Wang et al. Nat Commun 2018 (PMID 29610524)",
        "confidence": "medium",
        "color": "#f2a541",
    },
    {
        "id": "coiled_coil",
        "name": "Beam coiled coil",
        "category": "lever",
        "numbering": "mouse",
        "start": 1334, "end": 1365,
        "description": (
            "Coiled-coil segment at the distal end of the beam, annotated "
            "independently by UniProt in both species."
        ),
        "source": "UniProt Q92508/E2JF22 coiled-coil feature",
        "confidence": "high",
        "color": "#d98032",
    },
    {
        "id": "splice_1_1",
        "name": "Piezo1.1 spliced segment",
        "category": "lever",
        "numbering": "mouse",
        "start": 1382, "end": 1405,
        "description": (
            "The 24-residue segment removed in the Piezo1.1 splice isoform "
            "(PDB 6LQI). Identical in sequence between human and mouse, so the "
            "human deletion is 1388-1411."
        ),
        "source": "Geng et al. Neuron 2020, PDB 6LQI",
        "confidence": "high",
        "color": "#9d7bd8",
    },
]

# Terms deliberately NOT included, and why. Both appear in secondary sources
# but neither survived verification against the primary literature:
#   "clasp" — not an established PIEZO1 term; appears to be a conflation of the
#             latch with Guo & MacKinnon's "cross-helices".
#   "latch" — used inconsistently across papers with no agreed residue range.
# Adding a domain with invented boundaries would be worse than omitting it,
# because the viewer would colour residues with false confidence.

# Elements derived from UniProt features by explicit rules.
DERIVED_SPEC = {
    "anchor": {
        "name": "Anchor domain",
        "category": "lever",
        "rule": "cytoplasmic topological domain between TM36 and TM37",
        "description": (
            "Amphipathic intracellular domain wedged between the blade and the "
            "pore module. Transmits beam motion to the outer helix. Its apex "
            "(human P2113/F2114) acts as a resistive brake on the inner helix; "
            "the same region also sits in a cholesterol/CRAC context."
        ),
        "source": ("UniProt Q92508 topology; Li/Cox/Martinac Channels 2021; "
                   "Buyan et al. Biophys J 2020 (PMID 32949489)"),
        "confidence": "high",
        "color": "#4fc3c7",
    },
    "outer_helix": {
        "name": "Outer helix (OH)",
        "category": "pore",
        "rule": "transmembrane segment 37",
        "description": "Penultimate transmembrane helix; forms the outer wall of the pore module.",
        "source": "UniProt Q92508 transmembrane features",
        "confidence": "high",
        "color": "#5b8def",
    },
    "cap": {
        "name": "Cap / C-terminal extracellular domain (CED)",
        "category": "cap",
        "rule": "extracellular topological domain between TM37 and TM38",
        "description": (
            "Trimeric extracellular cap sitting above the pore, 234 residues per "
            "protomer. Contains the C2411-C2415 disulfide and carries the Er "
            "blood-group antigens. Rotates during gating and controls ion "
            "permeation through lateral portals."
        ),
        "source": "UniProt Q92508 topology; Karamatic Crew et al. Blood 2023",
        "confidence": "high",
        "color": "#8ad35e",
    },
    "inner_helix": {
        "name": "Inner helix (IH, pore-lining)",
        "category": "pore",
        "rule": "transmembrane segment 38",
        "description": (
            "Final transmembrane helix. Three copies line the ion-conduction "
            "pathway and form the transmembrane gate."
        ),
        "source": "UniProt Q92508 transmembrane features",
        "confidence": "high",
        "color": "#e8556d",
    },
    "ctd": {
        "name": "C-terminal domain (CTD)",
        "category": "gate",
        "rule": "cytoplasmic topological domain after TM38",
        "description": (
            "Intracellular domain beneath the pore forming the cytoplasmic "
            "vestibule and the constriction that sets ion selectivity. Hosts the "
            "acidic residues controlling divalent permeation and most "
            "gain-of-function xerocytosis variants."
        ),
        "source": "UniProt Q92508 topology; Coste et al. Nature 2015 (PMID 26649819)",
        "confidence": "high",
        "color": "#e8556d",
    },
}

#: The nine four-helix transmembrane helical units of the blade, numbered from
#: the N-terminus (distal, membrane-remote) to the pore-proximal end.
THU_COUNT = 9
THU_COLORS = ["#1f3f7a", "#26508f", "#2d61a4", "#3472b9", "#3b83ce",
              "#5695d8", "#71a7e1", "#8cb9ea", "#a7cbf3"]


def _tm_by_index(uni: dict) -> list[dict]:
    return sorted(uni["transmembrane"], key=lambda d: d["start"])


def _topology_between(uni: dict, lo: int, hi: int, kind: str) -> dict | None:
    for t in sorted(uni["topology"], key=lambda d: d["start"]):
        if t["start"] > lo and t["end"] < hi and t["description"].lower().startswith(kind):
            return t
    return None


def build() -> list[dict]:
    uni = json.loads((RESOURCE_DIR / "uniprot_human.json").read_text())
    nm = load_numbering_map()
    tms = _tm_by_index(uni)
    if len(tms) != 38:
        raise RuntimeError(f"expected 38 TM segments, found {len(tms)}")

    domains: list[dict] = []

    def add(did, name, category, h_start, h_end, description, source,
            confidence, color, rule=None, extra=None):
        m_start, m_end = nm.convert_range(h_start, h_end, "human")
        rec = {
            "id": did, "name": name, "category": category,
            "human": {"start": h_start, "end": h_end},
            "mouse": {"start": m_start, "end": m_end},
            "description": description, "source": source,
            "confidence": confidence, "color": color,
        }
        if rule:
            rec["rule"] = rule
        if extra:
            rec.update(extra)
        domains.append(rec)

    # --- the blade: nine THUs of four transmembrane helices each -----------
    for i in range(THU_COUNT):
        group = tms[i * 4:(i + 1) * 4]
        start = 1 if i == 0 else group[0]["start"]
        end = group[-1]["end"]
        distal = i < 3
        add(
            f"thu{i + 1}",
            f"THU{i + 1} ({group[0]['name']}-{group[-1]['name']})",
            "blade", start, end,
            (f"Transmembrane helical unit {i + 1} of the mechanosensory blade, "
             f"comprising {group[0]['name']}-{group[-1]['name']}. "
             + ("Distal blade — not resolved in any experimental structure; "
                "modelled from AlphaFold." if distal else
                "Proximal blade — resolved in cryo-EM.")),
            "UniProt Q92508 transmembrane features, grouped in fours from the N-terminus",
            "high", THU_COLORS[i],
            rule=f"transmembrane segments {i * 4 + 1}-{i * 4 + 4}",
            extra={"thu_index": i + 1, "distal": distal,
                   "transmembrane": [t["name"] for t in group]},
        )

    # --- derived elements ---------------------------------------------------
    anchor_topo = _topology_between(uni, tms[35]["end"], tms[36]["start"], "cytoplasmic")
    if anchor_topo:
        s = DERIVED_SPEC["anchor"]
        add("anchor", s["name"], s["category"], anchor_topo["start"], anchor_topo["end"],
            s["description"], s["source"], s["confidence"], s["color"], s["rule"])

    for key, tm_index in (("outer_helix", 36), ("inner_helix", 37)):
        s = DERIVED_SPEC[key]
        tm = tms[tm_index]
        add(key, s["name"], s["category"], tm["start"], tm["end"],
            s["description"], s["source"], s["confidence"], s["color"], s["rule"],
            extra={"transmembrane": [tm["name"]]})

    cap_topo = _topology_between(uni, tms[36]["end"], tms[37]["start"], "extracellular")
    if cap_topo:
        s = DERIVED_SPEC["cap"]
        add("cap", s["name"], s["category"], cap_topo["start"], cap_topo["end"],
            s["description"], s["source"], s["confidence"], s["color"], s["rule"])

    ctd_topo = max((t for t in uni["topology"]
                    if t["start"] > tms[37]["end"]
                    and t["description"].lower().startswith("cytoplasmic")),
                   key=lambda t: t["end"] - t["start"], default=None)
    if ctd_topo:
        s = DERIVED_SPEC["ctd"]
        add("ctd", s["name"], s["category"], ctd_topo["start"], ctd_topo["end"],
            s["description"], s["source"], s["confidence"], s["color"], s["rule"])

    # --- literature elements (declared in mouse numbering) -----------------
    for spec in LITERATURE_DOMAINS:
        if spec["numbering"] != "mouse":
            raise ValueError("literature domains are expected in mouse numbering")
        h_start, h_end = nm.convert_range(spec["start"], spec["end"], "mouse")
        domains.append({
            "id": spec["id"], "name": spec["name"], "category": spec["category"],
            "human": {"start": h_start, "end": h_end},
            "mouse": {"start": spec["start"], "end": spec["end"]},
            "description": spec["description"], "source": spec["source"],
            "confidence": spec["confidence"], "color": spec["color"],
            "declared_numbering": "mouse",
        })

    domains.sort(key=lambda d: (d["human"]["start"] or 0))
    return domains


def main() -> int:
    domains = build()
    nm = load_numbering_map()
    out = {
        "protein": "PIEZO1",
        "numbering": {
            "human": {"accession": "Q92508", "length": 2521},
            "mouse": {"accession": "E2JF22", "length": 2547},
            "identity": round(nm.identity, 4),
            "note": ("The human/mouse offset is NOT constant; it varies from 0 "
                     "to +26 across the chain. Always convert via "
                     "piezo1.core.sequence.load_numbering_map()."),
        },
        "domains": domains,
    }
    dest = RESOURCE_DIR / "domains.json"
    dest.write_text(json.dumps(out, indent=1))

    print(f"wrote {dest} with {len(domains)} domains\n")
    print(f"{'id':14s} {'category':8s} {'human':>12s} {'mouse':>12s}  {'conf':6s} name")
    for d in domains:
        h, m = d["human"], d["mouse"]
        print(f"{d['id']:14s} {d['category']:8s} "
              f"{str(h['start']) + '-' + str(h['end']):>12s} "
              f"{str(m['start']) + '-' + str(m['end']):>12s}  "
              f"{d['confidence']:6s} {d['name']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
