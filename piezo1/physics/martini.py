"""Preparing a coarse-grained permeation simulation — and not running one.

Liu et al. 2025 answered the question this project cannot: where the current
actually goes. They ran coarse-grained molecular dynamics on a **truncated**
PIEZO1 construct at four transmembrane potentials, watched Na+ cross the
transmembrane gate into the inner vestibule, and then watched it leave
*sideways* through the intracellular lateral portals. Nothing in a
one-dimensional continuum model can produce that, and
:mod:`piezo1.analysis.liu2025_permeation` says so rather than approximating it.

This module prepares the system that would answer it here. It **does not run
anything**, and the separation is deliberate and enforced:

- :func:`prepare` builds the truncated construct and writes it with a manifest
  describing the box, the membrane, the ion concentration and the voltages.
  Everything it produces is an *input*.
- :class:`MartiniRun` is the only type that carries results, and the only way
  to obtain one is :func:`load_results`, which reads a file some real run
  wrote. There is no code path that manufactures one.
- :func:`results_available` is what every consumer must ask before quoting an
  MD number. Until a trajectory exists it returns ``False`` and the figures
  keep showing the continuum analogue, labelled as such.

**Why prepare at all, if it cannot run?** Because the preparation is the part
that encodes somebody else's methods, and it is the part that goes stale
silently. Their construct is not the deposited trimer: it is the central pore
module (residues 1,956-2,547), the beam (1,315-1,365) and the intracellular
lateral plug gate (1,401-1,421), in mouse numbering. Writing that down against
our own structures, now, with the numbering checked, is what makes the later
run a matter of compute rather than of re-reading the paper.

**What this is not.** It is not a force field, a topology or a parameter set.
Martinising, the membrane build, equilibration and the production run all live
in tools this project does not vendor, and the manifest names them rather than
pretending to be them.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["SEGMENTS", "SYSTEM", "MartiniSystem", "MartiniRun", "prepare",
           "write_inputs", "load_results", "results_available",
           "TOOLCHAIN"]

#: Their construct, in **mouse** numbering (all four of their entries are
#: mouse). A transcription of the STAR methods, not a choice made here.
SEGMENTS = {
    "pore_module": (1956, 2547),
    "beam": (1315, 1365),
    "lateral_plug_gate": (1401, 1421),
}

#: The conditions the manifest declares. Each is somebody else's protocol or a
#: stated convention; none is fitted here, and none is read back as a result.
SYSTEM = {
    "voltages_V": (0.0, -0.1, -0.25, -0.5),
    "production_ns": 1000,
    "extended_production_ns": 10000,
    "salt": "150 mM NaCl",
    "membrane": "POPC bilayer",
    "model": "coarse-grained (Martini-class)",
    "padding_A": 30.0,
    "source": ("Liu et al. Neuron 2025;113:590-604 (PMID 39719701), "
               "STAR Methods: coarse-grained system preparation, "
               "equilibration and production simulation"),
}

#: The tools a run needs, named rather than vendored. Nothing here invokes any
#: of them; this is what a reader has to install before the manifest is useful.
TOOLCHAIN = (
    "martinize2 (vermouth) — atomistic model to coarse-grained topology",
    "insane or TS2CG — bilayer and solvent box",
    "GROMACS or OpenMM with a Martini force field — equilibration, production",
)


@dataclass
class MartiniSystem:
    """A prepared input. Carries no results and cannot be made to."""

    structure: object                        # the truncated Structure
    kept: dict = field(default_factory=dict)  # segment -> residues retained
    missing: dict = field(default_factory=dict)  # segment -> residues absent
    numbering: str = ""
    source_entry: str = ""
    box_A: tuple = ()

    @property
    def n_residues(self) -> int:
        return int(len(np.unique(self.structure.res_seq)))

    @property
    def n_atoms(self) -> int:
        return int(self.structure.n_atoms)

    def manifest(self) -> dict:
        """Everything a run needs, and an explicit statement that it is input."""
        return {
            "is_input_only": True,
            "results": None,
            "note": ("A prepared system. No simulation has been run and no "
                     "number here is a result. See piezo1.physics.martini."),
            "source_entry": self.source_entry,
            "numbering": self.numbering,
            "segments": {k: list(v) for k, v in SEGMENTS.items()},
            "residues_kept": {k: len(v) for k, v in self.kept.items()},
            "residues_missing": {k: len(v) for k, v in self.missing.items()},
            "n_residues": self.n_residues,
            "n_atoms": self.n_atoms,
            "box_A": list(self.box_A),
            "conditions": dict(SYSTEM),
            "toolchain": list(TOOLCHAIN),
        }

    def summary(self) -> str:
        missing = sum(len(v) for v in self.missing.values())
        return (f"{self.source_entry}: {self.n_residues} residues in "
                f"{len(self.kept)} segments ({self.numbering} numbering), "
                f"{missing} requested residues not resolved — INPUT ONLY, "
                f"nothing has been simulated")


@dataclass(frozen=True)
class MartiniRun:
    """Results from a real trajectory. Only :func:`load_results` makes one."""

    voltages_V: tuple
    permeated: tuple                  # ions crossing, per voltage
    duration_ns: float
    source: str                       # where the trajectory came from
    structure: str = ""
    note: str = ""

    def currents_pA(self, valence: int = 1) -> np.ndarray:
        """Current implied by the counts, picoamperes.

        The same arithmetic :mod:`piezo1.render.flux` uses in reverse, so an
        MD count and a continuum current are compared in one unit rather than
        two.
        """
        from .charge import ELEMENTARY_CHARGE

        counts = np.asarray(self.permeated, dtype=float)
        seconds = self.duration_ns * 1e-9
        return counts * valence * ELEMENTARY_CHARGE / seconds * 1e12


def prepare(structure, entry: str = "") -> MartiniSystem:
    """Cut ``structure`` down to the construct they simulated.

    Refuses an entry whose numbering cannot be read from its own coordinates,
    because every range here is a residue number and applying mouse ranges to a
    human entry would silently keep the wrong 600 residues.
    """
    from ..core.numbering_check import piezo1_numbering
    from ..core.sequence import mouse_to_human

    numbering = piezo1_numbering(structure)
    if numbering is None:
        raise ValueError(
            "cannot read the numbering from these coordinates — the construct "
            "is defined by residue number and a wrong reading would keep the "
            "wrong segments")

    kept, missing, mask = {}, {}, np.zeros(structure.n_atoms, dtype=bool)
    present = set(int(v) for v in structure.res_seq.tolist())
    for name, (first, last) in SEGMENTS.items():
        wanted = []
        for mouse in range(first, last + 1):
            number = mouse if numbering == "mouse" else mouse_to_human(mouse)
            if number is not None:
                wanted.append(int(number))
        here = [n for n in wanted if n in present]
        kept[name] = here
        missing[name] = [n for n in wanted if n not in present]
        mask |= np.isin(structure.res_seq, here)

    truncated = structure.subset(mask)
    xyz = truncated.xyz.astype(float)
    pad = SYSTEM["padding_A"]
    box = tuple(float(v) for v in (xyz.max(axis=0) - xyz.min(axis=0) + 2 * pad))
    return MartiniSystem(structure=truncated, kept=kept, missing=missing,
                         numbering=numbering,
                         source_entry=entry or getattr(structure, "name", ""),
                         box_A=box)


def write_inputs(system: MartiniSystem, directory) -> Path:
    """Write the coordinates and the manifest. Returns the directory.

    The manifest is written with ``is_input_only`` and a null ``results`` field
    so that a later reader — or :func:`load_results` — cannot mistake a
    prepared directory for a finished one.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    pdb = directory / f"{system.source_entry or 'construct'}_md.pdb"
    system.structure.to_pdb(pdb)
    manifest = directory / "manifest.json"
    payload = system.manifest()
    payload["coordinates"] = pdb.name
    manifest.write_text(json.dumps(payload, indent=2))
    return directory


def results_available(directory) -> bool:
    """Whether a real trajectory's results sit in ``directory``.

    What every consumer asks before quoting an MD number. A prepared directory
    returns ``False``: the manifest it contains says ``results: null``, and a
    prepared system is not a simulated one.
    """
    path = Path(directory) / "results.json"
    if not path.exists():
        return False
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError):
        return False
    return bool(payload.get("permeated")) and bool(payload.get("duration_ns"))


def load_results(directory) -> MartiniRun:
    """Read a finished run's results. Raises when there are none.

    Deliberately has no fallback. A function that returned an estimate when the
    trajectory was missing would put a continuum number where a reader expects
    a simulated one, which is the single thing this module exists to prevent.
    """
    path = Path(directory) / "results.json"
    if not path.exists():
        raise FileNotFoundError(
            f"no results in {directory} — the system has been prepared, not "
            f"run. Nothing here estimates what the simulation would say.")
    payload = json.loads(path.read_text())
    for required in ("voltages_V", "permeated", "duration_ns", "source"):
        if required not in payload:
            raise ValueError(f"results.json is missing {required!r}; a results "
                             f"file must say what was run and for how long")
    if len(payload["voltages_V"]) != len(payload["permeated"]):
        raise ValueError("one permeation count is needed per voltage")
    return MartiniRun(
        voltages_V=tuple(float(v) for v in payload["voltages_V"]),
        permeated=tuple(int(v) for v in payload["permeated"]),
        duration_ns=float(payload["duration_ns"]),
        source=str(payload["source"]),
        structure=str(payload.get("structure", "")),
        note=str(payload.get("note", "")))
