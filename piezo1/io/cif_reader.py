"""Minimal, fast mmCIF/PDB reader specialised for large membrane-protein models.

BioPython's parser is correct but allocates a Python object per atom, which is
painfully slow for a 120 000-atom PIEZO1 trimer.  This module walks the file
once and produces contiguous numpy arrays, which is what the renderer and the
physics engine both want anyway.

Only the parts of the mmCIF grammar that appear in wwPDB coordinate files are
supported: ``data_`` blocks, key/value pairs, ``loop_`` tables, single/double
quoted strings and semicolon-delimited multi-line text fields.
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import Iterator

import numpy as np

__all__ = ["read_cif_atoms", "read_pdb_atoms", "read_structure_file", "parse_cif_categories"]

# Columns we lift out of _atom_site, mapped to (dtype, default).
_ATOM_FIELDS = {
    "group_PDB": ("U6", "ATOM"),
    "id": ("i4", 0),
    "type_symbol": ("U2", "C"),
    "label_atom_id": ("U6", ""),
    "label_alt_id": ("U2", "."),
    "label_comp_id": ("U5", "UNK"),
    "label_asym_id": ("U6", "A"),
    "label_entity_id": ("U6", "1"),
    "label_seq_id": ("U8", "."),
    "auth_seq_id": ("U8", "0"),
    "auth_comp_id": ("U5", ""),
    "auth_asym_id": ("U6", ""),
    "pdbx_PDB_ins_code": ("U2", "?"),
    "Cartn_x": ("f4", 0.0),
    "Cartn_y": ("f4", 0.0),
    "Cartn_z": ("f4", 0.0),
    "occupancy": ("f4", 1.0),
    "B_iso_or_equiv": ("f4", 0.0),
    "pdbx_PDB_model_num": ("U4", "1"),
}


def _open_text(path: Path) -> io.TextIOBase:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def _tokenize(line: str) -> list[str]:
    """Split one mmCIF data line honouring single and double quotes."""
    out: list[str] = []
    ws = " \t\r\n\x0b\x0c"
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in ws:
            i += 1
            continue
        if c == "#":
            break
        if c in "'\"":
            quote = c
            i += 1
            start = i
            # A quote only terminates when followed by whitespace or EOL.
            while i < n:
                if line[i] == quote and (i + 1 >= n or line[i + 1] in ws):
                    break
                i += 1
            out.append(line[start:i])
            i += 1
        else:
            start = i
            while i < n and line[i] not in ws:
                i += 1
            out.append(line[start:i])
    return out


def _iter_values(handle: Iterator[str], first: list[str]) -> Iterator[str]:
    """Yield tokens from a loop body, gluing multi-line ``;`` blocks together."""
    pending = list(first)
    while True:
        while pending:
            yield pending.pop(0)
        try:
            line = next(handle)
        except StopIteration:
            return
        if line.startswith(";"):
            chunk = [line[1:].rstrip("\n")]
            for cont in handle:
                if cont.startswith(";"):
                    break
                chunk.append(cont.rstrip("\n"))
            yield "\n".join(chunk)
            continue
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("_") or stripped.startswith("loop_") or stripped.startswith("data_"):
            # Caller detects the end of the loop by inspecting this sentinel.
            yield "\x00" + line
            return
        pending = _tokenize(line)


def read_cif_atoms(path: str | Path, model: int | None = 1) -> dict[str, np.ndarray]:
    """Read ``_atom_site`` from an mmCIF file into a dict of numpy arrays.

    Parameters
    ----------
    path:
        ``.cif`` or ``.cif.gz`` file.
    model:
        Keep only this ``pdbx_PDB_model_num``.  ``None`` keeps every model.

    Returns
    -------
    dict
        Keys are the ``_ATOM_FIELDS`` names plus ``"xyz"`` (an ``(N, 3)``
        float32 array).  Absent columns are filled with their defaults.
    """
    path = Path(path)
    headers: list[str] = []
    rows: list[list[str]] = []

    with _open_text(path) as fh:
        lines = iter(fh)
        in_loop = False
        for line in lines:
            s = line.strip()
            if s == "loop_":
                in_loop = True
                headers = []
                continue
            if in_loop and s.startswith("_atom_site."):
                headers.append(s.split(".", 1)[1].split()[0])
                continue
            if headers and in_loop and s and not s.startswith("_"):
                ncol = len(headers)
                buf: list[str] = []
                for tok in _iter_values(lines, _tokenize(line)):
                    if tok.startswith("\x00"):
                        break
                    buf.append(tok)
                    if len(buf) == ncol:
                        rows.append(buf)
                        buf = []
                break
            if s.startswith("_") or s.startswith("data_"):
                in_loop = False
                headers = []

    if not rows:
        raise ValueError(f"no _atom_site records found in {path}")

    idx = {name: i for i, name in enumerate(headers)}
    table = np.array(rows, dtype=object)
    out: dict[str, np.ndarray] = {}
    for name, (dtype, default) in _ATOM_FIELDS.items():
        if name in idx:
            col = table[:, idx[name]].astype("U16")
            if dtype.startswith("f") or dtype.startswith("i"):
                col = np.where(np.isin(col, (".", "?", "")), str(default), col)
                out[name] = col.astype(dtype)
            else:
                out[name] = col.astype(dtype)
        else:
            out[name] = np.full(len(rows), default, dtype=dtype)

    if model is not None and "pdbx_PDB_model_num" in idx:
        keep = out["pdbx_PDB_model_num"] == str(model)
        if keep.any():
            out = {k: v[keep] for k, v in out.items()}

    # auth_* fields are optional; fall back to label_* so downstream code can
    # rely on them always being populated.
    for auth, label in (("auth_asym_id", "label_asym_id"), ("auth_comp_id", "label_comp_id")):
        blank = (out[auth] == "") | (out[auth] == ".")
        out[auth] = np.where(blank, out[label], out[auth])

    out["xyz"] = np.stack([out["Cartn_x"], out["Cartn_y"], out["Cartn_z"]], axis=1)
    return out


def read_pdb_atoms(path: str | Path, model: int | None = 1) -> dict[str, np.ndarray]:
    """Read ATOM/HETATM records from a legacy PDB file."""
    path = Path(path)
    recs: list[tuple] = []
    current_model = 1
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                try:
                    current_model = int(line[10:14])
                except ValueError:
                    current_model += 1
                continue
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if model is not None and current_model != model:
                continue
            try:
                recs.append((
                    line[:6].strip(), int(line[6:11]), line[12:16].strip(),
                    line[16:17].strip() or ".", line[17:20].strip(),
                    line[21:22].strip() or "A", line[22:26].strip(),
                    line[26:27].strip() or "?",
                    float(line[30:38]), float(line[38:46]), float(line[46:54]),
                    float(line[54:60] or 1.0), float(line[60:66] or 0.0),
                    (line[76:78].strip() or line[12:16].strip()[:1]),
                    current_model,
                ))
            except ValueError:
                continue
    if not recs:
        raise ValueError(f"no ATOM records found in {path}")

    cols = list(zip(*recs))
    out = {
        "group_PDB": np.array(cols[0], dtype="U6"),
        "id": np.array(cols[1], dtype="i4"),
        "label_atom_id": np.array(cols[2], dtype="U6"),
        "label_alt_id": np.array(cols[3], dtype="U2"),
        "label_comp_id": np.array(cols[4], dtype="U5"),
        "label_asym_id": np.array(cols[5], dtype="U6"),
        "auth_asym_id": np.array(cols[5], dtype="U6"),
        "auth_seq_id": np.array(cols[6], dtype="U8"),
        "label_seq_id": np.array(cols[6], dtype="U8"),
        "pdbx_PDB_ins_code": np.array(cols[7], dtype="U2"),
        "Cartn_x": np.array(cols[8], dtype="f4"),
        "Cartn_y": np.array(cols[9], dtype="f4"),
        "Cartn_z": np.array(cols[10], dtype="f4"),
        "occupancy": np.array(cols[11], dtype="f4"),
        "B_iso_or_equiv": np.array(cols[12], dtype="f4"),
        "type_symbol": np.array([c.capitalize() for c in cols[13]], dtype="U2"),
        "pdbx_PDB_model_num": np.array(cols[14], dtype="U4"),
    }
    out["auth_comp_id"] = out["label_comp_id"]
    out["label_entity_id"] = np.full(len(recs), "1", dtype="U6")
    out["xyz"] = np.stack([out["Cartn_x"], out["Cartn_y"], out["Cartn_z"]], axis=1)
    return out


def read_structure_file(path: str | Path, model: int | None = 1) -> dict[str, np.ndarray]:
    """Dispatch to the mmCIF or PDB reader based on the file extension."""
    name = str(path).lower()
    if name.endswith((".cif", ".cif.gz", ".mmcif", ".mmcif.gz")):
        return read_cif_atoms(path, model=model)
    return read_pdb_atoms(path, model=model)


def parse_cif_categories(path: str | Path, categories: set[str]) -> dict[str, dict]:
    """Extract simple key/value items for the requested mmCIF categories.

    Loop tables are returned as dicts of lists; single records as dicts of
    strings.  ``_atom_site`` is deliberately ignored — use
    :func:`read_cif_atoms` for that.
    """
    wanted = {c.rstrip(".") for c in categories}
    result: dict[str, dict] = {}
    with _open_text(Path(path)) as fh:
        lines = iter(fh)
        for line in lines:
            s = line.strip()
            if s == "loop_":
                headers = []
                for nxt in lines:
                    t = nxt.strip()
                    if t.startswith("_"):
                        headers.append(t.split()[0])
                        continue
                    if not headers:
                        break
                    cat = headers[0].split(".")[0].lstrip("_")
                    if cat not in wanted:
                        break
                    keys = [h.split(".", 1)[1] for h in headers]
                    store = result.setdefault(cat, {k: [] for k in keys})
                    buf: list[str] = []
                    for tok in _iter_values(lines, _tokenize(nxt)):
                        if tok.startswith("\x00"):
                            break
                        buf.append(tok)
                        if len(buf) == len(keys):
                            for k, v in zip(keys, buf):
                                store[k].append(v)
                            buf = []
                    break
                continue
            if s.startswith("_") and "." in s:
                cat = s.split(".")[0].lstrip("_")
                if cat not in wanted:
                    continue
                parts = _tokenize(s)
                key = parts[0].split(".", 1)[1]
                if len(parts) > 1:
                    result.setdefault(cat, {})[key] = parts[1]
                else:
                    val_lines: list[str] = []
                    for nxt in lines:
                        if nxt.startswith(";"):
                            val_lines.append(nxt[1:].rstrip("\n"))
                            for cont in lines:
                                if cont.startswith(";"):
                                    break
                                val_lines.append(cont.rstrip("\n"))
                            break
                        if nxt.strip():
                            val_lines.append(nxt.strip())
                            break
                    result.setdefault(cat, {})[key] = "\n".join(val_lines).strip()
    return result
