"""A modelled detergent micelle — the envelope Figure 4b shows as density.

Figure 4b of Guo & MacKinnon 2017 is the unsharpened cryo-EM map contoured at
6 sigma, and what it shows around the protein is **detergent**: a digitonin
micelle wrapping the transmembrane belt and curved into the dome the whole
paper is about. That map is not in this project and never will be.

What this module builds is a *model* of where that micelle would sit, and the
distinction is not a formality — the published envelope is an observation and
this is a geometric construction, so :class:`MicelleEnvelope` carries
``is_observed = False`` and every consumer states it.

**The construction, and why it is this one.** The paper describes the belt it
is looking at: "the hydrophobic residues on the TM helices (flanked by charged
amino acids) form a clearly curved band on the trimer surface, matching the
micelle density". So the belt is defined by that description — apolar side
chains of transmembrane residues — and the envelope is the surface at a fixed
offset outside it, i.e. an iso-surface of the distance-to-belt field. That is a
rolled-ball surface with a large ball, which is exactly what a detergent shell
is: acyl chains packing against the hydrophobic surface out to a roughly
constant thickness.

It has two properties that matter. It is **calibratable** — around one atom it
must be a sphere of exactly the offset radius, around a line a capsule — and it
has **no free shape parameters**: one offset, and the shape is then the
protein's. A blob fitted to look like the published density would have neither.

What it is not: a simulation. Real detergent has a headgroup layer, a
disordered core and a curvature preference of its own, none of which is here.
The number worth taking from it is the **curvature of the belt** — which is a
measurement of the protein — rather than the thickness of the shell, which is
a parameter.

Angstrom throughout.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from ..config import RESOURCE_DIR
from ..parameters import PARAMETERS as _P

from .geometry import SphereFit, fit_sphere

__all__ = ["MicelleEnvelope", "belt_atoms", "distance_field",
           "build_micelle", "APOLAR_RESIDUES"]


#: Residues counted as hydrophobic for the belt. The standard apolar set plus
#: the two aromatics that sit at the interface; Tyr and Trp are included
#: because they are the classic interfacial anchors and excluding them leaves
#: gaps in the band at exactly the leaflet boundaries the micelle follows.
APOLAR_RESIDUES = frozenset({"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO",
                             "CYS", "TRP", "TYR", "GLY"})


@dataclass
class MicelleEnvelope:
    """The modelled micelle surface, as a mesh plus what it measures."""

    vertices: np.ndarray            # (n, 3) Angstrom
    faces: np.ndarray               # (m, 3) triangle indices
    normals: np.ndarray             # (n, 3) outward unit normals
    #: Sphere fitted to the envelope's own mid-surface — the curvature the
    #: micelle density is read for, and the thing to compare with the dome.
    sphere: SphereFit | None
    n_belt_atoms: int
    offset: float                   # Angstrom, shell thickness used
    grid_spacing: float
    #: Always False. A construction from coordinates, not a density map.
    is_observed: bool = False
    meta: dict = field(default_factory=dict)

    @property
    def n_vertices(self) -> int:
        return len(self.vertices)

    def area(self) -> float:
        """Surface area of the envelope, Angstrom^2, by triangle sum."""
        if not len(self.faces):
            return 0.0
        a, b, c = (self.vertices[self.faces[:, i]] for i in range(3))
        return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())

    def enclosed_volume(self) -> float:
        """Volume enclosed, Angstrom^3, by the divergence theorem.

        Signed, then absolute: marching cubes gives a consistently oriented
        closed surface, but which way is outward depends on the sign
        convention of the field, and a negative volume would be a units bug
        rather than a shape.
        """
        if not len(self.faces):
            return 0.0
        a, b, c = (self.vertices[self.faces[:, i]] for i in range(3))
        return float(abs(np.einsum("ij,ij->i", a, np.cross(b, c)).sum()) / 6.0)

    @property
    def caveat(self) -> str:
        return ("A MODELLED micelle, not the observed density. Figure 4b is "
                "the unsharpened cryo-EM map at 6 sigma; this is the surface "
                "at a fixed offset outside the hydrophobic belt, so its "
                "thickness is a parameter and only its curvature is a "
                "measurement of the protein.")

    def summary(self) -> str:
        curvature = ("" if self.sphere is None else
                     f"envelope curvature R = {self.sphere.radius / 10:.1f} nm "
                     f"(RMSE {self.sphere.rmse:.1f} A) · ")
        return (f"modelled micelle: {self.n_belt_atoms} belt atoms, "
                f"{self.offset:.0f} A shell · {curvature}"
                f"area {self.area() / 100:.0f} nm^2, "
                f"volume {self.enclosed_volume() / 1000:.0f} nm^3")


def belt_atoms(structure, reference: str,
               mask: np.ndarray | None = None) -> np.ndarray:
    """Boolean mask of the hydrophobic transmembrane belt.

    Apolar residues inside the annotated transmembrane segments, side-chain
    atoms only — the backbone runs through the helix core and contributes
    nothing to the surface the detergent packs against.

    ``reference`` names the committed UniProt resource, so a human entry is not
    read with mouse transmembrane ranges. Getting that wrong shifts the band by
    up to 26 residues and produces a micelle that looks entirely reasonable.
    """
    helices = json.loads(
        (RESOURCE_DIR / f"uniprot_{reference}.json").read_text())["transmembrane"]
    in_helix = np.zeros(structure.n_atoms, dtype=bool)
    for helix in helices:
        in_helix |= ((structure.res_seq >= helix["start"])
                     & (structure.res_seq <= helix["end"]))
    backbone = np.isin(structure.atom_name, ("N", "C", "O"))
    keep = (in_helix
            & np.isin(structure.res_name, tuple(APOLAR_RESIDUES))
            & structure.mask_protein() & (~structure.hetero)
            & (~backbone))
    if mask is not None:
        keep &= np.asarray(mask)
    return keep


def distance_field(points: np.ndarray, spacing: float, pad: float
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Distance from each grid node to the nearest point, and the grid origin.

    The field whose level sets are the envelope. Returned rather than
    thresholded here so a caller can contour it at more than one offset
    without paying for the tree twice.
    """
    points = np.asarray(points, dtype=np.float64)
    if len(points) == 0:
        raise ValueError("no points to build a field around")
    low = points.min(axis=0) - pad
    high = points.max(axis=0) + pad
    shape = np.maximum(np.ceil((high - low) / spacing).astype(int) + 1, 2)
    axes = [low[i] + spacing * np.arange(shape[i]) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
    tree = cKDTree(points)
    distance, _ = tree.query(grid.reshape(-1, 3), k=1)
    return distance.reshape(shape), low


def build_micelle(structure, reference: str, offset: float | None = None,
                  spacing: float | None = None,
                  mask: np.ndarray | None = None) -> MicelleEnvelope:
    """Build the modelled micelle envelope around a structure's belt.

    Raises when the belt is too sparse to enclose anything — a structure whose
    transmembrane helices are not resolved has no band for detergent to follow,
    and a surface built from a handful of atoms would be a set of disconnected
    blobs that still renders.
    """
    from skimage import measure

    if offset is None:
        offset = _P.value("micelle.offset")
    if spacing is None:
        spacing = _P.value("micelle.grid_spacing")

    belt = belt_atoms(structure, reference, mask)
    n_belt = int(belt.sum())
    if n_belt < 100:
        raise ValueError(
            f"only {n_belt} hydrophobic transmembrane atoms found — this entry "
            f"does not resolve enough of the belt for a micelle to follow. "
            f"Is {reference!r} the right numbering for it?")

    points = structure.xyz[belt].astype(np.float64)
    field, low = distance_field(points, spacing, pad=offset + 3.0 * spacing)
    vertices, faces, normals, _ = measure.marching_cubes(
        field, level=float(offset), spacing=(spacing, spacing, spacing))
    vertices = vertices + low
    # marching_cubes returns gradient normals of the *distance* field, which
    # point towards increasing distance, i.e. outward. Normalised here so the
    # renderer does not have to.
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = normals / np.maximum(lengths, 1e-9)

    # Curvature of the belt itself rather than of the offset surface: the
    # offset adds a constant to the radius of any sphere, so fitting the
    # envelope and subtracting is the same number with the parameter removed.
    try:
        sphere = fit_sphere(points, iterations=4,
                            trim=_P.value("geometry.sphere_trim"))
    except ValueError:
        sphere = None

    return MicelleEnvelope(
        vertices=vertices.astype(np.float32), faces=faces.astype(np.int32),
        normals=normals.astype(np.float32), sphere=sphere,
        n_belt_atoms=n_belt, offset=float(offset), grid_spacing=float(spacing),
        meta={"reference": reference,
              "grid_shape": list(field.shape),
              "apolar_residues": sorted(APOLAR_RESIDUES),
              "definition": ("iso-surface of the distance to the apolar "
                             "side-chain atoms of the annotated transmembrane "
                             "helices"),
              "citation": "guo2017"})
