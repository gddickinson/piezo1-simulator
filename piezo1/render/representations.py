"""Turning a :class:`~piezo1.core.structure.Structure` into GPU batches.

A :class:`MolecularView` owns every batch belonging to one loaded model and
knows how to rebuild them when the style or colouring changes. Keeping that
knowledge here means the Qt layer only has to say "cartoon, coloured by
domain" and never touches an OpenGL buffer.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from ..core.structure import Structure
from . import colormaps
from .geometry_builders import Mesh, build_cartoon, build_tube
from .scene import Scene
from .spline import assign_secondary_structure

__all__ = ["Style", "ColorBy", "MolecularView"]


#: Ball radius in ball-and-stick, Angstrom. Roughly a quarter of a carbon's
#: van der Waals radius, the usual convention.
BALL_RADIUS = 0.42
#: Bond cylinder radius in ball-and-stick, Angstrom.
BOND_RADIUS = 0.20
#: Radius for the sticks-only style, used for both the bonds and the joint
#: spheres at each atom so the two are indistinguishable.
STICK_RADIUS = 0.28


class Style(str, Enum):
    CARTOON = "cartoon"
    TUBE = "tube"
    BACKBONE = "backbone"
    SPHERES = "spheres"
    BALLS = "balls"
    STICKS = "sticks"
    BALL_AND_STICK = "ball_and_stick"
    SURFACE_POINTS = "surface_points"


class ColorBy(str, Enum):
    DOMAIN = "domain"
    CHAIN = "chain"
    SECONDARY = "secondary"
    BFACTOR = "bfactor"
    PLDDT = "plddt"
    ELEMENT = "element"
    UNIFORM = "uniform"
    VALUE = "value"
    #: Kyte-Doolittle hydropathy on a *fixed* -4.5..+4.5 scale, for the same
    #: reason POTENTIAL is fixed: an auto-ranged hydropathy map cannot be
    #: compared with another structure or with a published panel.
    HYDROPHOBICITY = "hydrophobicity"
    #: Electrostatic potential, on a *fixed* diverging scale. Separate from
    #: VALUE because VALUE auto-ranges, and an auto-ranged potential map paints
    #: an almost-neutral surface as violently charged.
    POTENTIAL = "potential"


SS_COLORS = np.array([
    [0.55, 0.58, 0.66],   # coil
    [0.94, 0.42, 0.42],   # helix
    [0.98, 0.82, 0.35],   # strand
], dtype=np.float32)


@dataclass
class ChainTrace:
    """Per-chain C-alpha data cached so restyling does not re-scan the model."""

    chain: str
    indices: np.ndarray        # atom indices of the CA atoms
    xyz: np.ndarray
    res_seq: np.ndarray
    ss: np.ndarray


@dataclass
class MolecularView:
    """All batches belonging to one structure, plus the styling state."""

    scene: Scene
    structure: Structure
    name: str = "model"
    style: Style = Style.CARTOON
    color_by: ColorBy = ColorBy.DOMAIN
    ligands_as_spheres: bool = True
    #: Entity categories to draw. Empty means "everything", so existing callers
    #: are unaffected. See :mod:`piezo1.core.entities` — a deposited PIEZO1 file
    #: also contains lipids, detergent, glycan and, in six entries, the MDFIC
    #: auxiliary subunit, and which of those to show is a choice.
    visible_entities: frozenset = frozenset()
    #: Residue numbers to draw, or None for all. Set by the component
    #: selector. Applied in `_entity_filter`, which every representation goes
    #: through — cartoon, spheres, bonds and the traces — so a component
    #: cannot be half-applied.
    visible_residues: frozenset | None = None
    values: np.ndarray | None = None          # per-atom scalar for ColorBy.VALUE
    highlight: np.ndarray | None = None       # per-atom bool
    #: Colour used by :attr:`ColorBy.UNIFORM`. Per-view rather than global, so
    #: several structures shown together can be told apart.
    uniform_color: tuple = (0.55, 0.62, 0.75)
    traces: list[ChainTrace] = field(default_factory=list)
    _species: str = "human"

    # ------------------------------------------------------------------ init

    def __post_init__(self) -> None:
        self._entities = None
        self._palette = colormaps.load_domain_palette(self._species)
        self._build_traces()

    def _build_traces(self) -> None:
        st = self.structure
        ca = self._entity_filter(st.mask_ca())
        self.traces = []
        for chain in st.chains:
            m = ca & (st.chain == chain)
            n = int(m.sum())
            if n < 4:
                continue
            idx = np.flatnonzero(m)
            xyz = st.xyz[idx].astype(np.float64)
            seq = st.res_seq[idx]
            # Break the trace where the chain is discontinuous, otherwise the
            # spline draws a straight bar across an unmodelled loop. PIEZO1
            # cryo-EM models are full of these gaps.
            self.traces.append(ChainTrace(chain=chain, indices=idx, xyz=xyz,
                                          res_seq=seq,
                                          ss=assign_secondary_structure(xyz)))

    # ---------------------------------------------------------------- colours

    def atom_colors(self) -> np.ndarray:
        st = self.structure
        if self.color_by is ColorBy.DOMAIN:
            return colormaps.domain_colors(st, self._palette)
        if self.color_by is ColorBy.CHAIN:
            return colormaps.chain_colors(st)
        if self.color_by is ColorBy.BFACTOR:
            return colormaps.bfactor_colors(st)
        if self.color_by is ColorBy.PLDDT:
            return colormaps.plddt_colors(st)
        if self.color_by is ColorBy.HYDROPHOBICITY:
            return colormaps.hydrophobicity_colors(st)
        if self.color_by is ColorBy.ELEMENT:
            return st.element_colors()
        if self.color_by is ColorBy.POTENTIAL and self.values is not None:
            return colormaps.potential_colors(self.values)
        if self.color_by is ColorBy.VALUE and self.values is not None:
            return colormaps.value_colors(self.values)
        if self.color_by is ColorBy.SECONDARY:
            out = np.tile(SS_COLORS[0], (st.n_atoms, 1))
            for tr in self.traces:
                out[tr.indices] = SS_COLORS[tr.ss]
            return out
        return colormaps.uniform_color(st, self.uniform_color)

    def set_species(self, species: str) -> None:
        self._species = species
        self._palette = colormaps.load_domain_palette(species)

    # ------------------------------------------------------------- rebuilding

    def set_visible_residues(self, residues) -> None:
        """Restrict what is drawn to these residue numbers, or None for all.

        Rebuilds the **traces** as well as the batches. `traces` is built once
        in `__post_init__`, so setting the field and calling `rebuild()` filtered
        the atom representations and left the cartoon drawing the whole chain —
        which is most of what is on screen at the default style, so hiding 97%
        of the atoms changed the picture by about a tenth.
        """
        self.visible_residues = None if residues is None else frozenset(residues)
        self._build_traces()
        self.rebuild()

    def rebuild(self) -> None:
        """Drop and recreate every batch for the current style."""
        self.clear()
        colors = self.atom_colors()
        if self.style in (Style.CARTOON, Style.TUBE, Style.BACKBONE):
            self._build_ribbons(colors)
        if self.style in (Style.SPHERES, Style.BALLS, Style.STICKS,
                          Style.BALL_AND_STICK, Style.SURFACE_POINTS):
            self._build_atoms(colors)
        if self.ligands_as_spheres:
            self._build_ligands()

    def _segments(self, tr: ChainTrace, max_gap: int = 1):
        """Split a chain trace at unmodelled gaps in residue numbering."""
        breaks = np.flatnonzero(np.diff(tr.res_seq) > max_gap) + 1
        return np.split(np.arange(len(tr.xyz)), breaks)

    def _build_ribbons(self, colors: np.ndarray) -> None:
        mesh = Mesh.empty()
        for tr in self.traces:
            for seg in self._segments(tr):
                if len(seg) < 4:
                    continue
                xyz = tr.xyz[seg]
                col = colors[tr.indices[seg]].astype(np.float64)
                if self.style is Style.CARTOON:
                    part = build_cartoon(xyz, tr.ss[seg], col)
                elif self.style is Style.TUBE:
                    part = build_tube(xyz, col, radius=0.85, sides=10)
                else:
                    part = build_tube(xyz, col, radius=0.35, sides=6,
                                      subdivisions=3)
                mesh = mesh.concat(part)
        if mesh.n_vertices:
            batch = self.scene.mesh(f"{self.name}:ribbon", two_sided=True)
            batch.upload(mesh.positions, mesh.normals, mesh.colors, mesh.indices)

    def _build_atoms(self, colors: np.ndarray) -> None:
        st = self.structure
        if self.style is Style.SURFACE_POINTS:
            mask = ~st.hetero
            radii = st.vdw_radii()
        elif self.style is Style.SPHERES:
            mask = np.ones(st.n_atoms, dtype=bool)
            radii = st.vdw_radii()
        else:
            mask = np.ones(st.n_atoms, dtype=bool)
            radii = np.full(st.n_atoms, BALL_RADIUS, dtype=np.float32)
        mask = self._entity_filter(mask)
        flags = np.zeros(st.n_atoms, dtype=np.float32)
        if self.highlight is not None:
            flags[self.highlight] = 1.0

        if self.style is Style.STICKS:
            # Sticks still need a sphere at every atom, or each bond ends in a
            # flat disc and a branch point shows the seam between two of them.
            # They are drawn at the stick's own radius so they read as joints
            # rather than as balls, which is the whole difference between this
            # style and ball-and-stick.
            radii = np.full(st.n_atoms, STICK_RADIUS, dtype=np.float32)

        batch = self.scene.spheres(f"{self.name}:atoms")
        batch.upload(st.xyz[mask], radii[mask], colors[mask], flags[mask])

        if self.style in (Style.STICKS, Style.BALL_AND_STICK):
            self._build_bonds(colors, mask=mask)

    def _build_bonds(self, colors: np.ndarray, cutoff: float = 1.95,
                     mask: np.ndarray | None = None) -> None:
        """Covalent bonds as cylinders, coloured from each end's atom.

        ``mask`` is the atom selection actually drawn. Bonds used to be built
        from every atom regardless, so switching a category off in the Model
        panel hid its spheres and left its bonds hanging in space.
        """
        from scipy.spatial import cKDTree
        st = self.structure
        tree = cKDTree(st.xyz)
        pairs = np.asarray(list(tree.query_pairs(cutoff)), dtype=np.int64)
        if len(pairs) == 0:
            return
        # Skip anything that would bond across chains at this distance.
        keep = st.chain[pairs[:, 0]] == st.chain[pairs[:, 1]]
        if mask is not None:
            keep &= mask[pairs[:, 0]] & mask[pairs[:, 1]]
        pairs = pairs[keep]
        if len(pairs) == 0:
            return
        radius = STICK_RADIUS if self.style is Style.STICKS else BOND_RADIUS
        batch = self.scene.cylinders(f"{self.name}:bonds")
        batch.upload(st.xyz[pairs[:, 0]], st.xyz[pairs[:, 1]],
                     np.full(len(pairs), radius, np.float32),
                     colors[pairs[:, 0]], colors[pairs[:, 1]])

    def entity_map(self):
        """Cached per-atom entity classification for this structure."""
        if getattr(self, "_entities", None) is None:
            from ..core.entities import classify
            self._entities = classify(self.structure)
        return self._entities

    def _entity_filter(self, mask: np.ndarray) -> np.ndarray:
        """Drop atoms the user has switched off, by category or by component."""
        if self.visible_residues is not None:
            mask = mask & np.isin(self.structure.res_seq,
                                  np.fromiter(self.visible_residues, dtype=int))
        if not self.visible_entities:
            return mask
        allowed = self.entity_map().mask(*self.visible_entities)
        return mask & allowed

    def _build_ligands(self) -> None:
        st = self.structure
        mask = self._entity_filter(st.mask_ligands())
        if not mask.any():
            return
        batch = self.scene.spheres(f"{self.name}:ligands")
        batch.upload(st.xyz[mask], st.vdw_radii()[mask] * 0.85,
                     st.element_colors()[mask],
                     np.zeros(int(mask.sum()), np.float32))

    # ------------------------------------------------------------- animation

    def update_coords(self, xyz: np.ndarray) -> None:
        """Push new coordinates. Rebuilds meshes; spheres take the fast path."""
        # Coordinates change but the composition does not, so the cached
        # classification stays valid through a morph or a mode animation.
        self.structure = self.structure.copy_with_coords(xyz)
        for tr in self.traces:
            tr.xyz = self.structure.xyz[tr.indices].astype(np.float64)
        if self.style in (Style.CARTOON, Style.TUBE, Style.BACKBONE):
            self.rebuild()
        else:
            batch = self.scene.get(f"{self.name}:atoms")
            if batch is not None and batch.count:
                mask = (np.ones(self.structure.n_atoms, dtype=bool)
                        if self.style is not Style.SURFACE_POINTS
                        else ~self.structure.hetero)
                batch.update_centers(self.structure.xyz[mask])
            lig = self.scene.get(f"{self.name}:ligands")
            if lig is not None and lig.count:
                lig.update_centers(self.structure.xyz[self.structure.mask_ligands()])

    # ------------------------------------------------------------------ misc

    def clear(self) -> None:
        for suffix in ("ribbon", "atoms", "bonds", "ligands"):
            self.scene.remove(f"{self.name}:{suffix}")

    def set_visible(self, visible: bool) -> None:
        for suffix in ("ribbon", "atoms", "bonds", "ligands"):
            self.scene.set_visible(f"{self.name}:{suffix}", visible)
