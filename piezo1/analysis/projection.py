"""Simulated 2-D projections — the counterpart of Guo & MacKinnon's Figure 2a,b.

Figure 2a and 2b are **2D class averages**: thousands of real particle images,
aligned and averaged. Nothing computed from an atomic model is a class average,
and this module does not pretend otherwise. What it computes is the thing a
class average is an estimate *of* — the projection of the molecule's density
along a viewing direction — so the two can be put side by side and the question
"does the model account for what the micrograph shows" can be asked.

The distinction matters enough to be structural: the result carries
``is_experimental = False`` and a ``caveat`` that every consumer prints, and the
panel registry files 2a and 2b as *analogues* rather than replications.

What is modelled and what is not:

* **modelled** — atomic scattering as a Gaussian per atom, weighted by atomic
  number, summed along the view direction and sampled at the paper's own
  1.3 A pixel;
* **not modelled** — the contrast transfer function, defocus, solvent
  subtraction, the detergent micelle (which is *visible* in Figure 2b and is
  most of why the side view looks like a curved wedge), radiation damage, and
  alignment error.

The micelle omission is the important one. Figure 2b's dome-shaped envelope is
substantially micelle density, and no projection of a protein model will
reproduce it. :func:`project` therefore has nothing to say about the micelle
and the registry records that panel as only partly replicable.

Angstrom throughout. Returns images as plain arrays so nothing here needs a
plotting library.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["Projection", "project", "project_views", "ATOMIC_NUMBER",
           "STANDARD_VIEWS"]


#: Atomic numbers, as scattering weights. Hydrogen is absent from cryo-EM
#: models and from this table; carbon is the default because a protein is
#: mostly carbon and an unrecognised element is far more likely to be an
#: unusual carbon than an unusual metal.
ATOMIC_NUMBER = {"H": 1, "C": 6, "N": 7, "O": 8, "P": 15, "S": 16,
                 "NA": 11, "MG": 12, "K": 19, "CA": 20, "MN": 25, "FE": 26,
                 "CU": 29, "ZN": 30, "SE": 34, "CL": 17, "F": 9, "BR": 35,
                 "I": 53, "GD": 64}
_DEFAULT_Z = 6

#: The two views Figure 2 shows: down the three-fold axis, and perpendicular
#: to it. Named by what they show rather than by an axis letter, because which
#: Cartesian axis is the three-fold depends on how the structure was framed.
STANDARD_VIEWS = ("top", "side")


@dataclass
class Projection:
    """A simulated projection image."""

    image: np.ndarray                # (n, n), arbitrary density units
    pixel_size: float                # A per pixel
    resolution: float                # A, the Gaussian applied
    view: str
    direction: np.ndarray            # unit vector projected along
    n_atoms: int
    #: Always False. A projection of a model is not an experimental average,
    #: and the flag exists so that nothing downstream can lose track of that.
    is_experimental: bool = False
    meta: dict = field(default_factory=dict)

    @property
    def extent_A(self) -> float:
        """Field of view, Angstrom."""
        return float(self.image.shape[0] * self.pixel_size)

    def scale_bar_pixels(self, nanometres: float = 10.0) -> float:
        """Length of a scale bar in pixels. Figure 2's bar is 10 nm."""
        return float(nanometres * 10.0 / self.pixel_size)

    @property
    def caveat(self) -> str:
        return ("A projection of an atomic model, not a 2D class average. No "
                "CTF, no defocus, no solvent, and no detergent micelle — and "
                "the micelle is much of what Figure 2b's envelope shows.")

    def summary(self) -> str:
        return (f"{self.view} view, {self.image.shape[0]}x{self.image.shape[1]} "
                f"px at {self.pixel_size:g} A ({self.extent_A / 10:.0f} nm "
                f"across), {self.resolution:g} A blur, {self.n_atoms} atoms")


def _view_direction(structure, view: str, axis=None) -> np.ndarray:
    """Unit vector to project along, for a named view.

    The three-fold axis is *measured* rather than assumed to be z. A structure
    straight out of a deposited file is in whatever frame the depositors used,
    and projecting a trimer down the wrong axis produces a perfectly plausible
    picture of nothing in particular.
    """
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis

    if axis is None:
        blocks, _ = protomer_blocks(structure)
        if len(blocks) < 3:
            raise ValueError(
                "need three protomers to find the three-fold axis; pass "
                "`axis` explicitly for a structure that is not a trimer")
        axis = detect_c3_axis(blocks)
    direction = np.asarray(axis.direction, dtype=np.float64)
    if view == "top":
        return direction / np.linalg.norm(direction)
    if view == "side":
        # Any direction perpendicular to the axis. Take the one closest to a
        # global axis so repeated runs of the same structure agree.
        seed = np.eye(3)[int(np.argmin(np.abs(direction)))]
        perpendicular = seed - direction * np.dot(seed, direction)
        return perpendicular / np.linalg.norm(perpendicular)
    raise ValueError(f"unknown view {view!r}; expected one of {STANDARD_VIEWS}")


def project(structure, view: str = "top", pixel_size: float | None = None,
            resolution: float | None = None, size: int | None = None,
            axis=None, mask: np.ndarray | None = None) -> Projection:
    """Project a structure's density along a viewing direction.

    Each atom contributes a Gaussian of integrated weight equal to its atomic
    number. The Gaussian's standard deviation is ``resolution / 2.355`` — the
    full width at half maximum is the stated resolution, which is the
    convention a map's resolution is quoted in.

    The projection is done by binning atoms into pixels and convolving once,
    rather than by evaluating a Gaussian per atom: for 35,718 atoms the two
    agree to rounding and the second is thousands of times slower.
    """
    if pixel_size is None:
        pixel_size = _P.value("projection.pixel_size")
    if resolution is None:
        resolution = _P.value("projection.resolution")

    st = structure
    sel = np.ones(st.n_atoms, bool) if mask is None else np.asarray(mask)
    xyz = st.xyz[sel].astype(np.float64)
    if len(xyz) == 0:
        raise ValueError("no atoms selected")
    weight = np.array([ATOMIC_NUMBER.get(str(e).upper(), _DEFAULT_Z)
                       for e in st.element[sel]], dtype=np.float64)

    direction = _view_direction(st, view, axis)
    # Two axes spanning the image plane.
    seed = np.eye(3)[int(np.argmin(np.abs(direction)))]
    e1 = seed - direction * np.dot(seed, direction)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(direction, e1)

    centred = xyz - xyz.mean(axis=0)
    u, v = centred @ e1, centred @ e2

    if size is None:
        span = 2.0 * max(np.abs(u).max(), np.abs(v).max()) + 6.0 * resolution
        size = int(np.ceil(span / pixel_size))
        size += size % 2                       # keep the centre on a pixel edge
    half = size * pixel_size / 2.0

    bins = np.linspace(-half, half, size + 1)
    image, _, _ = np.histogram2d(u, v, bins=[bins, bins], weights=weight)

    sigma_px = (resolution / 2.3548200450309493) / pixel_size
    if sigma_px > 0:
        from scipy.ndimage import gaussian_filter
        image = gaussian_filter(image, sigma_px, mode="constant")

    return Projection(
        image=image, pixel_size=float(pixel_size), resolution=float(resolution),
        view=view, direction=direction, n_atoms=int(len(xyz)),
        meta={"sigma_px": float(sigma_px),
              "total_weight": float(weight.sum()),
              "structure": st.name,
              "scale_bar_10nm_px": float(10.0 * 10.0 / pixel_size),
              "citation": "guo2017",
              "not_a_class_average": True})


def project_views(structure, views: tuple[str, ...] = STANDARD_VIEWS,
                  **kw) -> dict[str, Projection]:
    """Both of Figure 2's views, sharing one measured three-fold axis.

    Sharing the axis matters: recovering it twice can pick opposite senses on a
    near-symmetric structure, and the two panels would then not be of the same
    orientation.
    """
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis

    axis = kw.pop("axis", None)
    if axis is None:
        blocks, _ = protomer_blocks(structure)
        if len(blocks) >= 3:
            axis = detect_c3_axis(blocks)
    return {view: project(structure, view=view, axis=axis, **kw)
            for view in views}
