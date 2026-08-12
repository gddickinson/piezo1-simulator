"""Getting a computed number out of this application and into another viewer.

Conservation, mechanical coupling, perturbation response, mode displacement and
the wetting score are all per-residue scalars. Every one of them was trapped
here or inside a JSON blob: :meth:`Structure.to_pdb` writes coordinates, and the
standard interoperability route — put the scalar in the **B-factor column** and
open the file in PyMOL or ChimeraX — did not exist.

It is about twenty lines. The other hundred are the three things that make the
file honest:

**Unmeasured is not zero.** A residue the analysis could not score must not
arrive in another viewer looking like a residue that scored zero — that is a
confident wrong number in somebody else's session, where none of this project's
guards can reach it. Unscored atoms get **occupancy 0.00** and a B-factor at the
floor of the column. Select on the occupancy: ``q < 0.5`` in PyMOL,
``@@occupancy<0.5`` in ChimeraX. That is the selector to use because it works
whatever the data are — the first version of this used a *negative* sentinel
and told readers to select on the sign, which the wetting energies broke
immediately by being negative themselves. The count is reported and written
into the header.

**The column quantises.** ``%6.2f`` carries two decimals, so a 0–1 score
arrives with about a hundred distinguishable levels. That is a real loss and
the export states it rather than letting a reader assume the file carries what
the array did. ``scale`` multiplies before writing — a 0–1 score at
``scale=100`` uses the column properly — and the factor goes in the header so
the number can be undone.

**The column has a range, and it is narrower than it looks.** The B-factor
field is six characters in ``F6.2``, so the widest values it can hold are
``999.99`` and ``-99.99`` — seven-character numbers like ``-999.99`` overflow
the field and shift every column after it, producing a file that still parses
and is wrong. This module got that constant wrong on its first attempt and
wrote exactly such a line; the limit is now **derived by formatting the value
and measuring it** rather than asserted, so it cannot be wrong again.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["ScalarExport", "SENTINEL", "COLUMN_STEP", "COLUMN_RANGE",
           "COLUMN_WIDTH", "fits_column", "write_scalar_pdb",
           "read_scalar_pdb"]

#: Written into the B-factor of any atom the analysis did not score: the floor
#: of the column, so it cannot collide with a real value of anything this
#: project computes. Deliberately *not* a small negative number — the wetting
#: energies are negative, and a sentinel inside the data's own range is not a
#: sentinel. The occupancy flag is the selector; this is the fallback for a
#: reader who only has the one column.
SENTINEL = -99.99

#: The B-factor column is ``F6.2`` — six characters including the sign and the
#: point. Properties of the PDB format, not choices made here.
COLUMN_WIDTH = 6
COLUMN_STEP = 0.01
#: Derived, not asserted. `-999.99` is seven characters and overflows the
#: field; writing it shifts every column after it into a file that still parses
#: and is wrong.
COLUMN_RANGE = (-99.99, 999.99)


@dataclass(frozen=True)
class ScalarExport:
    """What was written, and what the format could not carry."""

    path: Path
    n_atoms: int
    n_measured: int
    n_unmeasured: int
    value_range: tuple = (float("nan"), float("nan"))
    scale: float = 1.0
    levels: int = 0                     # distinguishable steps after scaling

    @property
    def fraction_measured(self) -> float:
        return self.n_measured / self.n_atoms if self.n_atoms else 0.0

    def summary(self) -> str:
        low, high = self.value_range
        text = (f"{self.path.name}: {self.n_measured:,} of {self.n_atoms:,} "
                f"atoms carry a value, range {low:.3g} to {high:.3g}")
        if self.scale != 1.0:
            text += f", written x{self.scale:g}"
        text += (f" · {self.levels} distinguishable levels in the B-factor "
                 f"column · {self.n_unmeasured:,} unscored atoms written as "
                 f"b={SENTINEL:g} with occupancy 0.00, so they are selectable "
                 f"and cannot be read as a score of zero")
        return text


def write_scalar_pdb(structure, values, path, scale: float = 1.0,
                     name: str = "scalar") -> ScalarExport:
    """Write ``structure`` with ``values`` in the B-factor column.

    ``values`` is either a mapping of residue number to value, or a per-atom
    array. A mapping is the usual case: every per-residue analysis in this
    project keys on residue number, and expanding it here keeps that expansion
    in one place rather than in each caller.

    Raises rather than truncating when a value will not fit the column. A
    number written as ``******`` or silently clipped is exactly the failure this
    module exists to prevent, one file further downstream.
    """
    path = Path(path)
    per_atom, measured = _expand(structure, values)

    finite = per_atom[measured]
    if not len(finite):
        raise ValueError("nothing to export — no atom carries a value")

    scaled = finite * float(scale)
    low, high = float(scaled.min()), float(scaled.max())
    for edge in (low, high):
        if not fits_column(edge):
            raise ValueError(
                f"values span {low:.4g} to {high:.4g} after scaling, which the "
                f"PDB B-factor column ({COLUMN_RANGE[0]} to {COLUMN_RANGE[1]}, "
                f"{COLUMN_WIDTH} characters) cannot hold — pass a smaller "
                f"`scale`")
    if low <= SENTINEL:
        raise ValueError(
            f"a scaled value reaches {low:.4g}, at or below the {SENTINEL} "
            f"sentinel that marks unscored atoms — they would be "
            f"indistinguishable. Offset or rescale the values first.")

    written = np.full(structure.n_atoms, SENTINEL, dtype=float)
    written[measured] = per_atom[measured] * float(scale)
    occupancy = measured.astype(float)

    lines = [
        f"REMARK   1 B-FACTOR COLUMN CARRIES {name.upper()}, NOT A B-FACTOR",
        f"REMARK   1 WRITTEN BY piezo1; scale x{scale:g}; "
        f"unscored atoms b={SENTINEL:g} occupancy 0.00",
        f"REMARK   1 {int(measured.sum())} of {structure.n_atoms} atoms scored",
    ]
    for i in range(structure.n_atoms):
        atom = str(structure.atom_name[i])
        label = f" {atom:<3}" if len(atom) < 4 else atom
        lines.append(
            f"{'HETATM' if structure.hetero[i] else 'ATOM  '}"
            f"{(i + 1) % 100000:5d} {label:<4}{'':1}"
            f"{str(structure.res_name[i]):>3} "
            f"{str(structure.chain[i])[:1]:>1}"
            f"{int(structure.res_seq[i]) % 10000:4d}{'':4}"
            f"{structure.xyz[i, 0]:8.3f}{structure.xyz[i, 1]:8.3f}"
            f"{structure.xyz[i, 2]:8.3f}"
            f"{occupancy[i]:6.2f}{written[i]:6.2f}"
            f"{'':10}{str(structure.element[i]):>2}")
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")

    span = high - low
    return ScalarExport(
        path=path, n_atoms=int(structure.n_atoms),
        n_measured=int(measured.sum()),
        n_unmeasured=int(structure.n_atoms - measured.sum()),
        value_range=(float(finite.min()), float(finite.max())),
        scale=float(scale),
        levels=int(round(span / COLUMN_STEP)) + 1)


def fits_column(value: float) -> bool:
    """Whether ``value`` fits the six-character B-factor field.

    Measured by formatting it, not by comparing against a remembered limit —
    the remembered limit was wrong.
    """
    return len(f"{float(value):{COLUMN_WIDTH}.2f}") == COLUMN_WIDTH


def read_scalar_pdb(path):
    """Read back the two columns, for checking the round trip numerically.

    Deliberately parses the fixed columns here rather than going through
    :mod:`piezo1.io.cif_reader`: the point of the check is that the *file* says
    what it should, and reading it with the same code that wrote it would only
    prove the two agree.
    """
    occupancies, b_factors = [], []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        occupancies.append(float(line[54:60]))
        b_factors.append(float(line[60:66]))
    return np.array(occupancies), np.array(b_factors)


def _expand(structure, values):
    """``values`` -> (per-atom array, per-atom measured mask)."""
    n = int(structure.n_atoms)
    if isinstance(values, dict):
        per_atom = np.full(n, np.nan, dtype=float)
        residues = structure.res_seq
        for residue, value in values.items():
            if value is None or not np.isfinite(float(value)):
                continue
            per_atom[residues == int(residue)] = float(value)
    else:
        per_atom = np.asarray(values, dtype=float)
        if per_atom.shape != (n,):
            raise ValueError(
                f"expected one value per atom ({n}) or a residue mapping, "
                f"got shape {per_atom.shape}")
    return per_atom, np.isfinite(per_atom)
