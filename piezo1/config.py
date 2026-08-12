"""Project-wide paths, constants and runtime settings.

Every other module imports paths from here rather than computing its own, so
that relocating the data directories is a one-line change.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# Filesystem layout
# --------------------------------------------------------------------------

PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_DIR.parent

#: Curated, hand-authored annotation files that ship with the package.
RESOURCE_DIR = PACKAGE_DIR / "resources"

#: External material downloaded from public databases. Git-ignored; fully
#: regenerable with ``python -m piezo1.io.fetch``.
REF_DIR = Path(os.environ.get("PIEZO1_REF_DIR", PROJECT_ROOT / "ref"))
STRUCTURE_DIR = REF_DIR / "structures"
SEQUENCE_DIR = REF_DIR / "sequences"
LIGAND_DIR = REF_DIR / "ligands"
PAPER_DIR = REF_DIR / "papers"
RESEARCH_DIR = REF_DIR / "research"

#: Derived artefacts computed by the application (normal modes, morphs, ...).
DATA_DIR = Path(os.environ.get("PIEZO1_DATA_DIR", PROJECT_ROOT / "data"))
CACHE_DIR = DATA_DIR / "cache"
DERIVED_DIR = DATA_DIR / "derived"

_ALL_DIRS = (
    STRUCTURE_DIR, SEQUENCE_DIR, LIGAND_DIR, PAPER_DIR, RESEARCH_DIR,
    CACHE_DIR, DERIVED_DIR,
)


def ensure_dirs() -> None:
    """Create every directory the application writes to."""
    for d in _ALL_DIRS:
        d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Reference molecules
# --------------------------------------------------------------------------

#: UniProt accession for canonical human PIEZO1 (2521 aa).
HUMAN_ACC = "Q92508"
#: UniProt accession for canonical mouse Piezo1 (2547 aa).
MOUSE_ACC = "E2JF22"
#: UniProt accession for canonical human PIEZO2 (2752 aa).
HUMAN_PIEZO2_ACC = "Q9H5I5"
#: The two invertebrate PIEZOs with deposited structures. They are the most
#: distant homologues available and therefore the strongest generality control:
#: PIEZO2 answers "is this PIEZO1 or the fold?" within mammals, and these ask
#: it across half a billion years.
WORM_PIEZO_ACC = "A0A061ACU2"        # C. elegans PEZO-1, 2442 aa
FLY_PIEZO_ACC = "M9MSG8"             # Drosophila Piezo, 2551 aa

#: UniProt accession for canonical mouse Piezo2 (2822 aa). 6KG7 is deposited in
#: this numbering, so the paralogue comparison needs it rather than the human
#: entry — the two differ by 70 residues and no constant offset relates them.
MOUSE_PIEZO2_ACC = "Q8CD54"

#: The three reviewed PIEZOs with no deposited structure. They exist here for
#: the *sequence* comparison, which is a different question from which
#: numbering a coordinate file is in — see ``core.numbering_check``.
#:
#: Rat is the third mammalian PIEZO1 and much of the electrophysiology
#: literature is rat. The other two are why the family is worth enumerating at
#: all: PIEZO is not a metazoan invention. A plant has one and an amoeba has
#: one, which pushes the generality question this project asks with PIEZO2
#: (mammals) and PEZO-1/dPIEZO (~800 Myr) out to the root of the eukaryotes.
RAT_ACC = "Q0KL00"                   # Rattus norvegicus Piezo1, 2535 aa
PLANT_PIEZO_ACC = "F4IN58"           # Arabidopsis thaliana PIEZO, 2462 aa
DICTY_PIEZO_ACC = "Q54S52"           # Dictyostelium discoideum pzoA, 3080 aa

HUMAN_LENGTH = 2521
MOUSE_LENGTH = 2547

#: PIEZO channels are homotrimers with C3 symmetry about the pore axis.
N_PROTOMERS = 3


# --------------------------------------------------------------------------
# Physical constants (SI unless noted)
# --------------------------------------------------------------------------

KB = 1.380649e-23           # Boltzmann constant, J/K
T_ROOM = 298.15             # K
KT_ROOM = KB * T_ROOM       # J
KT_IN_KCAL_MOL = 0.5921     # kcal/mol at 298.15 K
ANGSTROM = 1e-10            # m
NM = 1e-9                   # m
#: 1 mN/m of membrane tension expressed in kT per nm^2 at 298 K.
MNM_PER_KT_NM2 = 1e-3 * NM * NM / KT_ROOM


@dataclass
class RenderSettings:
    """Defaults for the OpenGL viewport."""

    gl_major: int = 4
    gl_minor: int = 1
    samples: int = 4
    background: tuple[float, float, float, float] = (0.055, 0.063, 0.086, 1.0)
    ambient_occlusion: bool = True
    depth_cue: bool = True
    fov_degrees: float = 35.0
    target_fps: int = 60


@dataclass
class AppSettings:
    """Top-level runtime settings, mutable from the GUI."""

    render: RenderSettings = field(default_factory=RenderSettings)
    default_structure: str = "8YEZ"
    #: Membrane tension shown on first load, in mN/m.
    default_tension: float = 0.0
    temperature: float = T_ROOM
    max_anm_modes: int = 40
    anm_cutoff: float = 15.0        # Angstrom, Ca-Ca contact cutoff
    anm_gamma: float = 1.0          # spring constant, arbitrary units


SETTINGS = AppSettings()
