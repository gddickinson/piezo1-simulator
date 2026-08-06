"""Saving and restoring a working session.

A session records *what you were looking at*, not the data itself: which
structure, how it was drawn, what was selected, which analyses had been run and
with what parameters. Reopening it re-derives everything from the same inputs,
so a session file stays small and never goes stale against a re-downloaded
structure.

It deliberately does **not** store coordinates or results. A file that carried
its own copy of the numbers would let a session drift silently out of step with
the code that produced them, which is the opposite of reproducibility.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from .. import __version__

__all__ = ["Session", "save_session", "load_session", "SESSION_FORMAT"]

#: Bumped whenever the stored fields change incompatibly.
SESSION_FORMAT = 1


@dataclass
class Session:
    """The state needed to put the application back where it was."""

    structure: str = ""
    species: str = "human"
    style: str = "cartoon"
    color_by: str = "domain"
    show_ligands: bool = True
    radius_scale: float = 1.0

    #: Camera: quaternion, pivot, distance.
    camera_rotation: list[float] = field(default_factory=lambda: [1.0, 0.0, 0.0, 0.0])
    camera_pivot: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    camera_distance: float = 300.0

    selected_residues: list[int] = field(default_factory=list)
    selection_label: str = ""

    #: Analyses that had been run, with the parameters they used.
    analyses: dict = field(default_factory=dict)
    morph: dict = field(default_factory=dict)
    notes: str = ""

    format_version: int = SESSION_FORMAT
    software_version: str = __version__
    saved_at: str = ""

    def as_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "Session":
        version = data.get("format_version", 0)
        if version > SESSION_FORMAT:
            raise ValueError(
                f"session format {version} is newer than this build "
                f"understands ({SESSION_FORMAT}); upgrade piezo1 to open it")
        known = {f for f in cls.__dataclass_fields__}
        # Unknown keys are dropped rather than raising, so a session written by
        # a newer minor version still opens with what it can.
        return cls(**{k: v for k, v in data.items() if k in known})

    def describe(self) -> str:
        bits = [f"{self.structure or 'no structure'} ({self.species})",
                f"{self.style}/{self.color_by}"]
        if self.selected_residues:
            bits.append(f"{len(self.selected_residues)} residues selected"
                        + (f" — {self.selection_label}" if self.selection_label else ""))
        if self.analyses:
            bits.append("analyses: " + ", ".join(sorted(self.analyses)))
        return " · ".join(bits)


def save_session(session: Session, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    session.saved_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    session.software_version = __version__
    session.format_version = SESSION_FORMAT
    path.write_text(json.dumps(session.as_dict(), indent=1))
    return path


def load_session(path: str | Path) -> Session:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"no session at {path}")
    return Session.from_dict(json.loads(path.read_text()))
