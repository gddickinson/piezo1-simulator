"""A per-residue feature table assembled from everything the project computes.

One row per residue, one column per quantity: mechanical coupling to the gate,
evolutionary constraint, burial, position, domain. It is the common
denominator for three separate uses — colouring a structure by any quantity,
nominating residues for testing, and supplying the structure-derived half of a
variant predictor.

**What this round deliberately does not do.** It computes no phenotype
comparison. The blind test in Round 7 returned a null result and that result
stands as recorded; any use of these features against the variant labels
requires a *new* pre-registration written first (Round 22). Assembling
predictors and evaluating them in the same breath is how a blind test stops
being blind.

A note on what the features mean. Several are strongly position-dependent by
construction — PRS gate response falls off with distance to the gate, and
betweenness is concentrated on whatever lies between the source and target.
That is not a defect, it is what they measure. The Round 7 diagnostic showed
that a *purely* positional score cannot distinguish substitutions at one site,
so these are the ingredients for a predictor, not a predictor.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.annotations import Annotations, load_annotations
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["ResidueFeatures", "build_feature_table", "MAX_ASA", "FEATURE_NOTES"]

#: Theoretical maximum solvent-accessible area per residue, Å² (Tien et al.
#: 2013, empirical Gly-X-Gly values). Used to turn absolute SASA into a
#: relative burial that is comparable between residue types.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

#: What each column means and how to read it.
FEATURE_NOTES = {
    "prs_gate_response": "mean response at the gate to a unit force here; "
                         "falls off with distance to the gate by construction",
    "prs_coupling": "absolute mechanical coupling to the rest of the protein "
                    "(the raw response matrix is symmetric, so effector and "
                    "sensor are the same quantity)",
    "betweenness": "fraction of blade-to-gate shortest paths passing through",
    "dcc_to_gate": "mean |dynamic cross-correlation| with the gate residues",
    "gating_amplitude": "displacement along the lowest symmetric (A) mode",
    "msf": "mean-square fluctuation summed over the low-frequency modes",
    "relative_sasa": "solvent accessibility as a fraction of the Gly-X-Gly "
                     "maximum; low means buried",
    "conservation": "1 - normalised Shannon entropy across vertebrate orthologs",
    "distance_to_gate": "closest heavy-atom approach to a gate residue, A",
    "distance_to_axis": "perpendicular distance from the three-fold axis, A",
    "n_contacts": "C-alpha neighbours within the elastic-network cutoff",
}


def _to_human(residue: int, species: str) -> int | None:
    """A residue number in this entry's numbering, as a human PIEZO1 number."""
    if species == "human":
        return residue
    from ..core.sequence import mouse_to_human
    return mouse_to_human(residue)


def _blade_range(species: str) -> tuple[int, int]:
    """The blade span, in the entry's own numbering.

    Human 570-1302 is where the deposited structures start resolving and where
    THU7 ends. Carried across by the alignment map rather than by an offset,
    for the same reason everything else in this project is.
    """
    human = (570, 1302)
    if species == "human":
        return human
    from ..core.sequence import human_to_mouse
    lo, hi = (human_to_mouse(human[0]), human_to_mouse(human[1]))
    return (lo or human[0], hi or human[1])


@dataclass
class ResidueFeatures:
    """A residue-by-feature table with provenance for each column."""

    residues: np.ndarray
    columns: dict[str, np.ndarray] = field(default_factory=dict)
    domains: dict[int, str] = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.residues)

    @property
    def names(self) -> list[str]:
        return list(self.columns)

    def get(self, residue: int) -> dict | None:
        idx = np.flatnonzero(self.residues == residue)
        if len(idx) == 0:
            return None
        i = int(idx[0])
        out = {"residue": int(residue), "domain": self.domains.get(int(residue))}
        out.update({k: float(v[i]) for k, v in self.columns.items()})
        return out

    def column(self, name: str) -> np.ndarray:
        return self.columns[name]

    def as_dict(self, name: str) -> dict[int, float]:
        """One column as ``{residue: value}``, for the ranking helpers."""
        return {int(r): float(v) for r, v in zip(self.residues, self.columns[name])}

    def rows(self) -> list[dict]:
        return [self.get(int(r)) for r in self.residues]

    def to_csv(self) -> str:
        import csv
        import io
        buf = io.StringIO()
        fields = ["residue", "domain"] + self.names
        writer = csv.DictWriter(buf, fieldnames=fields)
        writer.writeheader()
        for row in self.rows():
            writer.writerow({k: (round(v, 5) if isinstance(v, float) else v)
                             for k, v in row.items()})
        return buf.getvalue()

    def correlations(self) -> tuple[np.ndarray, list[str]]:
        """Pearson correlation between columns, ignoring non-finite entries."""
        names = [n for n in self.names
                 if np.isfinite(self.columns[n]).sum() > 10]
        matrix = np.full((len(names), len(names)), np.nan)
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                x, y = self.columns[a], self.columns[b]
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() > 10 and x[ok].std() > 0 and y[ok].std() > 0:
                    matrix[i, j] = float(np.corrcoef(x[ok], y[ok])[0, 1])
        return matrix, names

    def percentile(self, name: str) -> np.ndarray:
        """A column converted to percentile rank, for combining unlike units."""
        values = self.columns[name]
        ok = np.isfinite(values)
        out = np.full(len(values), np.nan)
        if ok.sum() > 1:
            order = np.argsort(np.argsort(values[ok]))
            out[ok] = order / (ok.sum() - 1)
        return out


def _protomer_blocks(st: Structure, min_ca: int = 300):
    chains = []
    for c in st.chains:
        m = st.mask_ca() & (st.chain == c)
        if m.sum() > min_ca:
            chains.append((st.xyz[m], st.res_seq[m]))
    if len(chains) < 3:
        raise ValueError("need three well-resolved protomers")
    common = set(chains[0][1].tolist())
    for _, seq in chains[1:3]:
        common &= set(seq.tolist())
    arr = np.array(sorted(common))
    return ([x[np.searchsorted(s, arr)].astype(float) for x, s in chains[:3]],
            arr)


def build_feature_table(structure: Structure,
                        annotations: Annotations | None = None,
                        n_modes: int = 30, cutoff: float = 15.0,
                        gate_group: str = "hydrophobic_gate",
                        blade_range: tuple[int, int] | None = None,
                        include_conservation: bool = True,
                        include_sasa: bool = True,
                        sasa_points: float | None = None) -> ResidueFeatures:
    """Compute every per-residue feature for one structure.

    Everything is averaged over the three protomers, because a homotrimer has
    three chemically identical copies of each residue and reporting one would
    make the answer depend on which chain the file happened to list first.

    **Every residue number here is in the entry's own numbering**, and that has
    to be established rather than assumed. Until Round 93 this defaulted to
    human annotation whatever it was handed, so on a *mouse* entry — which is
    most of the catalogue and is the entry the Round 48 validation used — the
    hydrophobic-gate group, the blade range and the conservation profile were
    all looked up at human residue numbers against mouse coordinates. The
    human/mouse offset is not constant and reaches 26 residues, so the gate
    columns pointed at residues a helix turn away and the conservation column
    dropped from rho = 0.81 to 0.29 against an independent measure of the same
    quantity.

    It survived because ``tests/test_features.py`` uses a **human** fixture,
    where the default is right. The numbering now comes from
    :func:`piezo1.core.numbering_check.piezo1_numbering`, and the
    human-anchored conservation profile is converted through
    :mod:`piezo1.core.sequence` rather than by an offset.
    """
    from ..core.numbering_check import piezo1_numbering

    species = piezo1_numbering(structure) or "human"
    if blade_range is None:
        blade_range = _blade_range(species)
    sasa_points = int(_P.value("sasa.n_points_fast")) if sasa_points is None else sasa_points
    from .allostery import (build_network, cross_correlation,
                            path_betweenness, perturbation_response)
    from ..physics.anm import ANM

    ann = annotations or load_annotations(species)
    blocks, residues = _protomer_blocks(structure)
    per = len(residues)
    coords = np.vstack(blocks)
    site_residues = np.tile(residues, 3)

    anm = ANM.from_trimer(blocks, cutoff=cutoff).build()
    modes = anm.calc_modes(n_modes=n_modes)
    anm.label_symmetry(modes)

    gate_residues = set(ann.group(gate_group).residues) if ann.group(gate_group) else set()
    gate_sites = [i for i, r in enumerate(site_residues) if r in gate_residues]
    blade_sites = [i for i, r in enumerate(site_residues)
                   if blade_range[0] <= r <= blade_range[1]]

    # The response matrix is used two ways but yields only ONE per-residue
    # summary. Both matrices are built because the gate-response column needs
    # the normalised one (relative pattern) while the coupling column needs the
    # raw one (absolute magnitude, and the only version where the row mean
    # means anything).
    prs = perturbation_response(modes, site_residues, normalise=True)
    prs_raw = perturbation_response(modes, site_residues, normalise=False)
    dcc = cross_correlation(modes)
    graph = build_network(coords, dcc, contact_cutoff=10.0)

    columns: dict[str, np.ndarray] = {}

    def collapse(per_site: np.ndarray) -> np.ndarray:
        """Average a per-site quantity over the three protomer copies."""
        return per_site.reshape(3, per).mean(axis=0)

    columns["prs_gate_response"] = collapse(prs.response_at(gate_sites)) \
        if gate_sites else np.full(per, np.nan)
    # ONE column, not two. The raw matrix is symmetric, so its row and column
    # means are literally the same numbers - effector and sensor coincide. Row
    # normalisation appears to break that, but the normalised column mean still
    # correlates with the raw row mean at 0.998, so it is the same quantity
    # wearing a different scale. Shipping both would look like two independent
    # lines of evidence and be one.
    columns["prs_coupling"] = collapse(prs_raw.effectiveness)

    bet = path_betweenness(graph, blade_sites, gate_sites, site_residues,
                           max_pairs=200) if gate_sites and blade_sites else {}
    columns["betweenness"] = np.array([bet.get(int(r), 0.0) for r in residues])

    columns["dcc_to_gate"] = (
        collapse(np.abs(dcc[:, gate_sites]).mean(axis=1))
        if gate_sites else np.full(per, np.nan))

    # The lowest three-fold-symmetric mode is the candidate gating coordinate:
    # isotropic tension is itself C3-symmetric, so only A modes can couple.
    symmetry = modes.symmetry if modes.symmetry is not None else []
    a_indices = [i for i, s in enumerate(symmetry) if s == "A"]
    gating_mode = a_indices[0] if a_indices else 0
    columns["gating_amplitude"] = collapse(
        np.linalg.norm(modes.vectors[gating_mode], axis=1))
    columns["msf"] = collapse(modes.msf())

    axis_direction = None
    try:
        from ..structure.superpose import detect_c3_axis
        axis = detect_c3_axis(blocks)
        axis_direction = axis
        columns["distance_to_axis"] = collapse(axis.radial(coords))
    except Exception:
        columns["distance_to_axis"] = np.full(per, np.nan)

    if gate_sites:
        gate_xyz = coords[gate_sites]
        d = np.linalg.norm(coords[:, None, :] - gate_xyz[None, :, :], axis=2).min(axis=1)
        columns["distance_to_gate"] = collapse(d)
    else:
        columns["distance_to_gate"] = np.full(per, np.nan)

    from scipy.spatial import cKDTree
    tree = cKDTree(coords)
    counts = np.array([len(tree.query_ball_point(p, cutoff)) - 1 for p in coords])
    columns["n_contacts"] = collapse(counts.astype(float))

    if include_sasa:
        from .measure import sasa
        mask = structure.mask_protein() & ~structure.hetero
        result = sasa(structure, n_points=sasa_points, mask=mask)
        totals: dict[int, list[float]] = {}
        for res, name, area in zip(result.residue_seq,
                                   structure.subset(mask).residue_name,
                                   result.residue):
            maximum = MAX_ASA.get(str(name))
            if maximum:
                totals.setdefault(int(res), []).append(float(area) / maximum)
        columns["relative_sasa"] = np.array(
            [float(np.mean(totals[int(r)])) if int(r) in totals else np.nan
             for r in residues])

    if include_conservation:
        try:
            from .conservation import conservation_profile, load_orthologs
            profile = conservation_profile(load_orthologs())
            # The profile is anchored on human Q92508 whatever entry this is,
            # so its keys are human residue numbers and the structure's are
            # not. Converted through the alignment map, never by subtraction.
            lookup = {int(r): float(c) for r, c, cov
                      in zip(profile.residues, profile.conservation, profile.coverage)
                      if cov >= 0.5}
            columns["conservation"] = np.array(
                [lookup.get(_to_human(int(r), species), np.nan) for r in residues])
        except FileNotFoundError:
            columns["conservation"] = np.full(per, np.nan)

    domains = {}
    for r in residues:
        d = ann.domain_at(int(r))
        domains[int(r)] = d.id if d else None

    return ResidueFeatures(
        residues=residues, columns=columns, domains=domains,
        meta={"structure": structure.name, "n_residues": per,
              "n_modes": n_modes, "cutoff": cutoff,
              "gating_mode_index": int(gating_mode),
              "gating_mode_symmetry": (modes.symmetry[gating_mode]
                                       if modes.symmetry is not None else "?"),
              "gate_residues": sorted(gate_residues),
              "notes": dict(FEATURE_NOTES),
              "caveat": ("Several columns are positional by construction. "
                         "These are ingredients for a predictor, not a "
                         "predictor; see docs/VALIDATION.md.")})
