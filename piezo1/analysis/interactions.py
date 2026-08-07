"""Detection of non-covalent interactions from coordinates.

Every criterion here is a published geometric definition with its cutoffs
stated, not a single distance threshold applied to everything. The defaults
follow PLIP (Salentin et al. 2015, *Nucleic Acids Res* 43:W443), which is the
most widely used reference implementation, so results are comparable with what
a reviewer would get from that server.

Cutoffs used:

===================  ========================================================
Hydrogen bond        donor–acceptor ≤ 3.5 Å. PLIP uses 4.1 Å, but only in
                     combination with an explicit D–H···A angle ≥ 100°, which
                     needs hydrogens. Without them, 4.1 Å alone admits far too
                     much — it produced 8005 "bonds" in one PIEZO1 trimer,
                     including arginine-pairs 2.3 Å apart. 3.5 Å is the
                     conventional heavy-atom donor–acceptor distance and is
                     what is used here, with the result flagged
                     ``heavy_atom_only``.
Salt bridge          ≤ 5.5 Å between charged group centroids
Hydrophobic contact  ≤ 4.0 Å between apolar carbons
π-stacking           ring centroids ≤ 5.5 Å; offset ≤ 2.0 Å;
                     planes ≤ 30° apart (parallel) or 60–90° (T-shaped)
Cation–π             cation to ring centroid ≤ 6.0 Å, offset ≤ 2.0 Å
Disulfide            S–S ≤ 2.5 Å
===================  ========================================================

The important caveat, stated once here and carried on every result: cryo-EM
structures of PIEZO1 have **no hydrogens**, and many side chains are modelled
at 3–4 Å resolution where rotamers are uncertain. Hydrogen bonds inferred from
heavy atoms alone are plausible contacts, not established bonds.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from scipy.spatial import cKDTree

from ..core.structure import Structure

__all__ = ["Interaction", "InteractionSet", "detect_interactions",
           "compare_interactions", "CUTOFFS"]

CUTOFFS = {
    "hbond_distance": 3.5,
    "hbond_angle": 100.0,
    "salt_bridge": 5.5,
    "hydrophobic": 4.0,
    "pi_stack_centroid": 5.5,
    "pi_stack_offset": 2.0,
    "pi_stack_parallel_angle": 30.0,
    "pi_stack_tshape_angle": 60.0,
    "cation_pi": 6.0,
    "disulfide": 2.5,
}

_DONORS = {"N", "O"}
_ACCEPTORS = {"O", "N"}
#: Charged groups, by residue and the atoms whose centroid carries the charge.
_POSITIVE = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",), "HIS": ("ND1", "NE2")}
_NEGATIVE = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
#: Aromatic rings, by residue and ring atoms.
_RINGS = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
}
_APOLAR_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS",
               "TYR"}


@dataclass(frozen=True)
class Interaction:
    kind: str
    distance: float
    atom_i: int
    atom_j: int
    res_i: int
    res_j: int
    chain_i: str
    chain_j: str
    name_i: str
    name_j: str
    atom_name_i: str = ""
    atom_name_j: str = ""
    angle: float | None = None
    note: str = ""

    def key(self) -> tuple:
        """Chain/residue-level identity, for comparing between structures."""
        a = (self.chain_i, self.res_i)
        b = (self.chain_j, self.res_j)
        return (self.kind,) + tuple(sorted([a, b]))

    def __str__(self) -> str:
        # Atom names are shown, not just residue names. Without them a
        # perfectly ordinary backbone carbonyl-to-amide bond between two
        # arginines reads as "ARG - ARG", which looks like a detection error.
        a = f"{self.name_i}{self.res_i}{self.chain_i}"
        b = f"{self.name_j}{self.res_j}{self.chain_j}"
        if self.atom_name_i:
            a += f".{self.atom_name_i}"
        if self.atom_name_j:
            b += f".{self.atom_name_j}"
        extra = f"  [{self.note}]" if self.note and "heavy_atom" not in self.note else ""
        return f"{self.kind}: {a} – {b}  {self.distance:.2f} Å{extra}"


@dataclass
class InteractionSet:
    interactions: list[Interaction] = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.interactions)

    def __iter__(self):
        return iter(self.interactions)

    def of_kind(self, kind: str) -> list[Interaction]:
        return [i for i in self.interactions if i.kind == kind]

    def counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for i in self.interactions:
            out[i.kind] = out.get(i.kind, 0) + 1
        return out

    def involving(self, residue: int, chain: str | None = None) -> list[Interaction]:
        return [i for i in self.interactions
                if (i.res_i == residue and (chain is None or i.chain_i == chain))
                or (i.res_j == residue and (chain is None or i.chain_j == chain))]

    def keys(self) -> set[tuple]:
        return {i.key() for i in self.interactions}


def _centroids(st: Structure, table: dict[str, tuple[str, ...]],
               mask: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    """Centroids of named atom groups, one per matching residue."""
    pts, info = [], []
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        return np.zeros((0, 3)), []
    keys = np.stack([st.chain[idx], st.res_seq[idx].astype("U8")], axis=1)
    order = np.lexsort((keys[:, 1], keys[:, 0]))
    idx = idx[order]
    start = 0
    for i in range(1, len(idx) + 1):
        same = (i < len(idx) and st.chain[idx[i]] == st.chain[idx[start]]
                and st.res_seq[idx[i]] == st.res_seq[idx[start]])
        if same:
            continue
        block = idx[start:i]
        rname = str(st.res_name[block[0]])
        wanted = table.get(rname)
        if wanted:
            sel = block[np.isin(st.atom_name[block], wanted)]
            if len(sel) >= (3 if table is _RINGS else 1):
                pts.append(st.xyz[sel].mean(axis=0).astype(float))
                info.append({"res": int(st.res_seq[block[0]]),
                             "chain": str(st.chain[block[0]]),
                             "name": rname, "atoms": sel})
        start = i
    return (np.array(pts) if pts else np.zeros((0, 3))), info


def _ring_normal(st: Structure, atoms: np.ndarray) -> np.ndarray:
    pts = st.xyz[atoms].astype(float)
    centred = pts - pts.mean(axis=0)
    _, _, vt = np.linalg.svd(centred)
    return vt[2] / np.linalg.norm(vt[2])


def detect_interactions(structure: Structure,
                        mask_a: np.ndarray | None = None,
                        mask_b: np.ndarray | None = None,
                        kinds: tuple[str, ...] | None = None,
                        min_sequence_separation: int | None = None,
                        ) -> InteractionSet:
    """Detect non-covalent interactions, optionally between two selections.

    With no masks, everything against everything. With both, only interactions
    that cross between them — which is what you want for a protein–ligand or
    protein–protein interface.

    ``min_sequence_separation`` suppresses trivial i,i+1 backbone neighbours
    within a chain; set it to 0 to keep them.
    """
    if min_sequence_separation is None:
        min_sequence_separation = int(_P.value("interactions.min_sequence_separation"))
    st = structure
    all_mask = np.ones(st.n_atoms, bool)
    a = all_mask if mask_a is None else np.asarray(mask_a)
    b = all_mask if mask_b is None else np.asarray(mask_b)
    cross_only = mask_a is not None and mask_b is not None
    want = set(kinds) if kinds else {"hydrogen_bond", "salt_bridge",
                                     "hydrophobic", "pi_stacking",
                                     "cation_pi", "disulfide"}
    found: list[Interaction] = []

    def allowed(i: int, j: int) -> bool:
        if cross_only and not ((a[i] and b[j]) or (b[i] and a[j])):
            return False
        if not (a[i] or b[i]) or not (a[j] or b[j]):
            return False
        if (st.chain[i] == st.chain[j]
                and abs(int(st.res_seq[i]) - int(st.res_seq[j]))
                < min_sequence_separation):
            return False
        return True

    def record(kind: str, i: int, j: int, d: float,
               ang: float | None = None, note: str = "") -> None:
        found.append(Interaction(
            kind=kind, distance=float(d), atom_i=int(i), atom_j=int(j),
            res_i=int(st.res_seq[i]), res_j=int(st.res_seq[j]),
            chain_i=str(st.chain[i]), chain_j=str(st.chain[j]),
            name_i=str(st.res_name[i]), name_j=str(st.res_name[j]),
            atom_name_i=str(st.atom_name[i]), atom_name_j=str(st.atom_name[j]),
            angle=ang, note=note))

    pool = np.flatnonzero(a | b)
    tree = cKDTree(st.xyz[pool].astype(float))

    # --- hydrogen bonds (heavy-atom geometry only) -----------------------
    if "hydrogen_bond" in want:
        polar = pool[np.isin(st.element[pool], list(_DONORS | _ACCEPTORS))]
        ptree = cKDTree(st.xyz[polar].astype(float))
        for p, q in ptree.query_pairs(CUTOFFS["hbond_distance"]):
            i, j = int(polar[p]), int(polar[q])
            if not allowed(i, j):
                continue
            # Without hydrogens, donor/acceptor character has to be inferred
            # from atom identity. Nitrogen in protein is almost always a donor,
            # so an N...N pair is donor-donor and cannot be a hydrogen bond.
            # The exception is histidine, whose ring nitrogens may be
            # unprotonated and therefore acceptors. Admitting all N...N pairs
            # produced arginine-pairs at 2.3 A being reported as bonds.
            if (st.element[i] == "N" and st.element[j] == "N"
                    and "HIS" not in (str(st.res_name[i]), str(st.res_name[j]))):
                continue
            d = float(np.linalg.norm(st.xyz[i] - st.xyz[j]))
            if d < 2.2:                   # too close to be a hydrogen bond
                continue
            record("hydrogen_bond", i, j, d,
                   note="heavy_atom_only: no hydrogens in this model")

    # --- salt bridges ----------------------------------------------------
    if "salt_bridge" in want:
        pos_xyz, pos_info = _centroids(st, _POSITIVE, a | b)
        neg_xyz, neg_info = _centroids(st, _NEGATIVE, a | b)
        if len(pos_xyz) and len(neg_xyz):
            d = np.linalg.norm(pos_xyz[:, None, :] - neg_xyz[None, :, :], axis=2)
            for p, q in zip(*np.where(d <= CUTOFFS["salt_bridge"])):
                i = int(pos_info[p]["atoms"][0])
                j = int(neg_info[q]["atoms"][0])
                if allowed(i, j):
                    record("salt_bridge", i, j, d[p, q])

    # --- hydrophobic contacts -------------------------------------------
    if "hydrophobic" in want:
        apolar = pool[(st.element[pool] == "C")
                      & np.isin(st.res_name[pool], list(_APOLAR_RES))]
        atree = cKDTree(st.xyz[apolar].astype(float))
        seen_pairs: set[tuple] = set()
        for p, q in atree.query_pairs(CUTOFFS["hydrophobic"]):
            i, j = int(apolar[p]), int(apolar[q])
            if not allowed(i, j):
                continue
            # One contact per residue pair, at the closest approach.
            k = tuple(sorted([(str(st.chain[i]), int(st.res_seq[i])),
                              (str(st.chain[j]), int(st.res_seq[j]))]))
            if k in seen_pairs:
                continue
            seen_pairs.add(k)
            record("hydrophobic", i, j,
                   float(np.linalg.norm(st.xyz[i] - st.xyz[j])))

    # --- aromatic stacking and cation-pi ---------------------------------
    ring_xyz, ring_info = _centroids(st, _RINGS, a | b)
    if len(ring_xyz) and {"pi_stacking", "cation_pi"} & want:
        normals = np.array([_ring_normal(st, r["atoms"]) for r in ring_info])

    if "pi_stacking" in want and len(ring_xyz) > 1:
        d = np.linalg.norm(ring_xyz[:, None, :] - ring_xyz[None, :, :], axis=2)
        for p, q in zip(*np.where(np.triu(d <= CUTOFFS["pi_stack_centroid"], 1))):
            i = int(ring_info[p]["atoms"][0])
            j = int(ring_info[q]["atoms"][0])
            if not allowed(i, j):
                continue
            inter = np.degrees(np.arccos(
                abs(np.clip(np.dot(normals[p], normals[q]), -1, 1))))
            vec = ring_xyz[q] - ring_xyz[p]
            offset = np.linalg.norm(vec - normals[p] * np.dot(vec, normals[p]))
            if inter <= CUTOFFS["pi_stack_parallel_angle"] \
                    and offset <= CUTOFFS["pi_stack_offset"]:
                record("pi_stacking", i, j, d[p, q], inter, note="parallel")
            elif inter >= CUTOFFS["pi_stack_tshape_angle"]:
                record("pi_stacking", i, j, d[p, q], inter, note="T-shaped")

    if "cation_pi" in want and len(ring_xyz):
        cat_xyz, cat_info = _centroids(st, _POSITIVE, a | b)
        if len(cat_xyz):
            d = np.linalg.norm(cat_xyz[:, None, :] - ring_xyz[None, :, :], axis=2)
            for p, q in zip(*np.where(d <= CUTOFFS["cation_pi"])):
                i = int(cat_info[p]["atoms"][0])
                j = int(ring_info[q]["atoms"][0])
                if not allowed(i, j):
                    continue
                vec = cat_xyz[p] - ring_xyz[q]
                offset = np.linalg.norm(vec - normals[q] * np.dot(vec, normals[q]))
                if offset <= CUTOFFS["pi_stack_offset"]:
                    record("cation_pi", i, j, d[p, q])

    # --- disulfides -------------------------------------------------------
    if "disulfide" in want:
        sulphur = pool[(st.element[pool] == "S") & (st.res_name[pool] == "CYS")]
        if len(sulphur) > 1:
            stree = cKDTree(st.xyz[sulphur].astype(float))
            for p, q in stree.query_pairs(CUTOFFS["disulfide"]):
                i, j = int(sulphur[p]), int(sulphur[q])
                if allowed(i, j):
                    record("disulfide", i, j,
                           float(np.linalg.norm(st.xyz[i] - st.xyz[j])))

    return InteractionSet(
        interactions=found,
        meta={"cutoffs": dict(CUTOFFS), "structure": st.name,
              "cross_selection_only": cross_only,
              "caveat": ("Hydrogens are absent from cryo-EM models, so "
                         "hydrogen bonds are inferred from heavy-atom geometry "
                         "alone and are plausible contacts, not established "
                         "bonds.")})


def compare_interactions(before: InteractionSet, after: InteractionSet
                         ) -> dict[str, list[Interaction]]:
    """Which interactions break, form, or persist between two states.

    For a mechanosensor this is the central question: the contacts that break
    when the dome flattens are the ones transmitting the conformational change.
    """
    a, b = before.keys(), after.keys()
    lost = [i for i in before if i.key() in (a - b)]
    gained = [i for i in after if i.key() in (b - a)]
    kept = [i for i in after if i.key() in (a & b)]
    return {"lost": lost, "gained": gained, "retained": kept}
