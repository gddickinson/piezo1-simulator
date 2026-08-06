"""Structure-of-arrays container for molecular models.

Everything downstream — renderer, elastic network model, morphing, analysis —
operates on :class:`Structure`.  Atom data live in parallel numpy arrays so
that selections are boolean masks and coordinate maths is vectorised.

A second, residue-level view is built once at load time: ``res_first`` indexes
the first atom of each residue, so per-residue quantities can be gathered with
``np.add.reduceat`` instead of Python loops.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from ..io.cif_reader import read_structure_file

__all__ = ["Structure", "AA3TO1", "ELEMENT_RADII", "ELEMENT_COLORS"]

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M", "UNK": "X",
}

#: van der Waals radii in Angstrom (Bondi 1964 / Rowland & Taylor 1996).
ELEMENT_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90,
    "NA": 2.27, "K": 2.75, "MG": 1.73, "CA": 2.31, "ZN": 1.39, "FE": 2.00,
    "GD": 2.10, "MN": 2.00, "CU": 1.40,
}
_DEFAULT_RADIUS = 1.70

#: CPK-ish palette tuned for a dark background.
ELEMENT_COLORS = {
    "C": (0.60, 0.63, 0.70), "N": (0.30, 0.47, 0.95), "O": (0.94, 0.28, 0.28),
    "S": (0.98, 0.82, 0.24), "P": (0.98, 0.58, 0.20), "H": (0.92, 0.92, 0.95),
    "SE": (1.00, 0.55, 0.35), "CL": (0.35, 0.85, 0.40), "F": (0.55, 0.92, 0.55),
    "NA": (0.55, 0.35, 0.90), "K": (0.45, 0.25, 0.85), "MG": (0.20, 0.80, 0.55),
    "CA": (0.25, 0.85, 0.45), "ZN": (0.55, 0.55, 0.70), "FE": (0.88, 0.45, 0.20),
    "GD": (0.85, 0.40, 0.85),
}
_DEFAULT_COLOR = (0.75, 0.55, 0.95)

_WATER = {"HOH", "DOD", "WAT"}


@dataclass
class Structure:
    """An immutable-by-convention molecular model.

    Attributes
    ----------
    xyz:
        ``(n_atoms, 3)`` float32 coordinates in Angstrom.
    element, atom_name, res_name, chain, alt_loc:
        ``(n_atoms,)`` string arrays.
    res_seq:
        ``(n_atoms,)`` int32 author residue numbers.
    hetero:
        ``(n_atoms,)`` bool, True for HETATM records.
    b_factor, occupancy:
        ``(n_atoms,)`` float32.  For predicted models ``b_factor`` holds pLDDT.
    """

    xyz: np.ndarray
    element: np.ndarray
    atom_name: np.ndarray
    res_name: np.ndarray
    res_seq: np.ndarray
    chain: np.ndarray
    hetero: np.ndarray
    b_factor: np.ndarray
    occupancy: np.ndarray
    alt_loc: np.ndarray
    entity: np.ndarray
    name: str = "structure"
    source: str = ""
    meta: dict = field(default_factory=dict)

    # Residue-level view, built by ``_build_residue_index``.
    res_first: np.ndarray = field(default=None, repr=False)   # type: ignore[assignment]
    res_atom_index: np.ndarray = field(default=None, repr=False)  # type: ignore[assignment]

    # ---------------------------------------------------------------- loading

    @classmethod
    def from_file(
        cls,
        path: str | Path,
        name: str | None = None,
        model: int | None = 1,
        keep_water: bool = False,
        keep_altloc: str = "A",
    ) -> "Structure":
        """Load a ``.cif``/``.pdb`` file (optionally gzipped)."""
        path = Path(path)
        raw = read_structure_file(path, model=model)

        seq = np.array([_int_or(v, -9999) for v in raw["auth_seq_id"]], dtype=np.int32)
        st = cls(
            xyz=np.ascontiguousarray(raw["xyz"], dtype=np.float32),
            element=np.char.upper(raw["type_symbol"].astype("U2")),
            atom_name=raw["label_atom_id"].astype("U6"),
            res_name=raw["auth_comp_id"].astype("U5"),
            res_seq=seq,
            chain=raw["auth_asym_id"].astype("U6"),
            hetero=(raw["group_PDB"] == "HETATM"),
            b_factor=raw["B_iso_or_equiv"].astype(np.float32),
            occupancy=raw["occupancy"].astype(np.float32),
            alt_loc=raw["label_alt_id"].astype("U2"),
            entity=raw["label_entity_id"].astype("U6"),
            name=name or path.stem.split(".")[0].upper(),
            source=str(path),
        )

        keep = np.ones(st.n_atoms, dtype=bool)
        if not keep_water:
            keep &= ~np.isin(st.res_name, list(_WATER))
        if keep_altloc:
            keep &= np.isin(st.alt_loc, (".", "?", "", keep_altloc))
        if not keep.all():
            st = st.subset(keep)
        st._build_residue_index()
        return st

    # -------------------------------------------------------------- basic API

    @property
    def n_atoms(self) -> int:
        return len(self.xyz)

    @property
    def n_residues(self) -> int:
        return 0 if self.res_first is None else len(self.res_first)

    @property
    def chains(self) -> list[str]:
        seen: list[str] = []
        for c in self.chain:
            if c not in seen:
                seen.append(str(c))
        return seen

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (f"<Structure {self.name}: {self.n_atoms} atoms, "
                f"{self.n_residues} residues, chains={''.join(self.chains)}>")

    _ARRAY_FIELDS = ("xyz", "element", "atom_name", "res_name", "res_seq",
                     "chain", "hetero", "b_factor", "occupancy", "alt_loc", "entity")

    def subset(self, mask: np.ndarray, name: str | None = None) -> "Structure":
        """Return a new :class:`Structure` containing only the masked atoms."""
        mask = np.asarray(mask)
        if mask.dtype != bool:
            idx = mask
        else:
            idx = np.flatnonzero(mask)
        kw = {f: getattr(self, f)[idx] for f in self._ARRAY_FIELDS}
        st = Structure(name=name or self.name, source=self.source,
                       meta=dict(self.meta), **kw)
        st._build_residue_index()
        return st

    def copy_with_coords(self, xyz: np.ndarray, name: str | None = None) -> "Structure":
        """Clone this structure, replacing only the coordinates."""
        kw = {f: getattr(self, f) for f in self._ARRAY_FIELDS if f != "xyz"}
        st = Structure(xyz=np.ascontiguousarray(xyz, dtype=np.float32),
                       name=name or self.name, source=self.source,
                       meta=dict(self.meta), **kw)
        st.res_first = self.res_first
        st.res_atom_index = self.res_atom_index
        return st

    # ------------------------------------------------------------- residue map

    def _build_residue_index(self) -> None:
        """Identify residue boundaries; assumes atoms of a residue are adjacent."""
        if self.n_atoms == 0:
            self.res_first = np.zeros(0, dtype=np.int64)
            self.res_atom_index = np.zeros(0, dtype=np.int64)
            return
        changed = np.ones(self.n_atoms, dtype=bool)
        changed[1:] = (
            (self.res_seq[1:] != self.res_seq[:-1])
            | (self.chain[1:] != self.chain[:-1])
            | (self.res_name[1:] != self.res_name[:-1])
        )
        self.res_first = np.flatnonzero(changed)
        self.res_atom_index = np.cumsum(changed) - 1

    # ------------------------------------------------------ residue-level views

    @property
    def residue_seq(self) -> np.ndarray:
        return self.res_seq[self.res_first]

    @property
    def residue_chain(self) -> np.ndarray:
        return self.chain[self.res_first]

    @property
    def residue_name(self) -> np.ndarray:
        return self.res_name[self.res_first]

    @property
    def residue_hetero(self) -> np.ndarray:
        return self.hetero[self.res_first]

    def one_letter_sequence(self, chain: str | None = None) -> str:
        """Single-letter sequence of the polymer residues (``X`` for unknown)."""
        names = self.residue_name
        keep = ~self.residue_hetero
        if chain is not None:
            keep &= self.residue_chain == chain
        return "".join(AA3TO1.get(str(n), "X") for n in names[keep])

    # ------------------------------------------------------------- selections

    def mask_ca(self) -> np.ndarray:
        return (self.atom_name == "CA") & (self.element == "C") & ~self.hetero

    def mask_backbone(self) -> np.ndarray:
        return np.isin(self.atom_name, ("N", "CA", "C", "O")) & ~self.hetero

    def mask_protein(self) -> np.ndarray:
        return np.isin(self.res_name, list(AA3TO1))

    def mask_chain(self, chains: str | Sequence[str]) -> np.ndarray:
        if isinstance(chains, str):
            chains = [chains]
        return np.isin(self.chain, list(chains))

    def mask_resrange(self, lo: int, hi: int, chains: Sequence[str] | None = None) -> np.ndarray:
        m = (self.res_seq >= lo) & (self.res_seq <= hi)
        if chains is not None:
            m &= self.mask_chain(chains)
        return m

    def mask_residues(self, residues: Iterable[int], chains: Sequence[str] | None = None) -> np.ndarray:
        m = np.isin(self.res_seq, np.asarray(list(residues), dtype=np.int32))
        if chains is not None:
            m &= self.mask_chain(chains)
        return m

    def mask_ligands(self, exclude: Sequence[str] = ()) -> np.ndarray:
        m = self.hetero & ~np.isin(self.res_name, list(_WATER) + list(exclude))
        return m

    # ------------------------------------------------------------- coordinates

    def ca_coords(self, chain: str | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(xyz, res_seq)`` for C-alpha atoms, optionally one chain."""
        m = self.mask_ca()
        if chain is not None:
            m &= self.chain == chain
        return self.xyz[m], self.res_seq[m]

    @property
    def center(self) -> np.ndarray:
        return self.xyz.mean(axis=0) if self.n_atoms else np.zeros(3, np.float32)

    def radius_of_gyration(self) -> float:
        d = self.xyz - self.center
        return float(np.sqrt((d * d).sum() / max(self.n_atoms, 1)))

    def bounding_box(self) -> tuple[np.ndarray, np.ndarray]:
        return self.xyz.min(axis=0), self.xyz.max(axis=0)

    def vdw_radii(self) -> np.ndarray:
        out = np.full(self.n_atoms, _DEFAULT_RADIUS, dtype=np.float32)
        for el, r in ELEMENT_RADII.items():
            out[self.element == el] = r
        return out

    def element_colors(self) -> np.ndarray:
        out = np.tile(np.array(_DEFAULT_COLOR, np.float32), (self.n_atoms, 1))
        for el, c in ELEMENT_COLORS.items():
            out[self.element == el] = c
        return out

    # ------------------------------------------------------------------ output

    def to_pdb(self, path: str | Path) -> None:
        """Write a PDB file (chain IDs truncated to one character)."""
        lines = []
        for i in range(self.n_atoms):
            nm = str(self.atom_name[i])
            name = f" {nm:<3}" if len(nm) < 4 else nm
            lines.append(
                f"{'HETATM' if self.hetero[i] else 'ATOM  '}"
                f"{(i + 1) % 100000:5d} {name:<4}{'':1}{str(self.res_name[i]):>3} "
                f"{str(self.chain[i])[:1]:>1}{int(self.res_seq[i]) % 10000:4d}{'':4}"
                f"{self.xyz[i, 0]:8.3f}{self.xyz[i, 1]:8.3f}{self.xyz[i, 2]:8.3f}"
                f"{self.occupancy[i]:6.2f}{self.b_factor[i]:6.2f}"
                f"{'':10}{str(self.element[i]):>2}"
            )
        lines.append("END")
        Path(path).write_text("\n".join(lines) + "\n")


def _int_or(value, default: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default
