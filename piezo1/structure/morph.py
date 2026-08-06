"""Conformational morphing between two experimental states.

Watching PIEZO1 flatten is the single most instructive thing this application
can show, but it has to be done honestly. A morph is an *interpolation*, not a
simulated trajectory: it shows a geometrically plausible path between two
observed endpoints, and says nothing about the barrier between them or the
order in which parts move. Every trajectory produced here records that.

Three methods are provided:

``linear``
    Straight-line interpolation of each atom. Fast, and correct at both
    endpoints, but the middle frames are physically wrong — atoms take chords
    through space, so C-alpha–C-alpha distances contract wherever the local
    motion is rotational. For PIEZO1's blades, which swing through large arcs,
    that contraction is severe.

``restrained`` (default)
    Linear interpolation followed by iterative restoration of C-alpha–C-alpha
    distances toward their interpolated targets. Cheap, robust, and it removes
    the chord artefact, which is what makes the middle frames look and measure
    like a real protein.

``modal``
    Project the endpoint difference onto the low-frequency elastic-network
    modes of the start structure and follow that subspace. The most physically
    motivated option: the path is constrained to motions the elastic network
    actually supports.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .superpose import match_protomers, superpose

__all__ = ["MorphTrajectory", "morph", "interpolate_linear",
           "restrained_morph", "modal_morph", "prepare_endpoints"]


@dataclass
class MorphTrajectory:
    """A sequence of coordinate frames between two endpoints."""

    frames: np.ndarray                    # (n_frames, n_sites, 3)
    method: str
    endpoint_rmsd: float
    #: Per-frame maximum deviation of a C-alpha–C-alpha distance from its
    #: interpolated target, in Angstrom. This is the honesty metric: a large
    #: value means the frame is geometrically strained.
    bond_error: np.ndarray = field(default_factory=lambda: np.zeros(0))
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.frames)

    @property
    def n_sites(self) -> int:
        return self.frames.shape[1]

    def at(self, fraction: float) -> np.ndarray:
        """Frame at a fractional position along the path, with interpolation."""
        f = float(np.clip(fraction, 0.0, 1.0)) * (len(self.frames) - 1)
        lo = int(np.floor(f))
        hi = min(lo + 1, len(self.frames) - 1)
        w = f - lo
        return (1.0 - w) * self.frames[lo] + w * self.frames[hi]

    def summary(self) -> str:
        return (f"{len(self)} frames, {self.method} method, "
                f"endpoint RMSD {self.endpoint_rmsd:.1f} A, "
                f"worst bond error {self.bond_error.max():.2f} A"
                if len(self.bond_error) else
                f"{len(self)} frames, {self.method} method")


# --------------------------------------------------------------------------
# Endpoint preparation
# --------------------------------------------------------------------------

def prepare_endpoints(start_blocks: list[np.ndarray], start_res: np.ndarray,
                      end_blocks: list[np.ndarray], end_res: np.ndarray
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Reduce two trimers to a common, correctly matched, superposed basis.

    Returns ``(start, end, common_residues, info)``.

    Two things must happen before any interpolation, and getting either wrong
    silently invalidates everything downstream: the two models must be reduced
    to residues resolved in *all six* protomers, and the protomers must be put
    into the right rotational correspondence — deposited chain labels are not a
    reliable guide to that.
    """
    common = np.array(sorted(set(start_res.tolist()) & set(end_res.tolist())))
    if len(common) < 50:
        raise ValueError(f"only {len(common)} residues common to both endpoints")

    a = [blk[np.searchsorted(start_res, common)] for blk in start_blocks]
    b = [blk[np.searchsorted(end_res, common)] for blk in end_blocks]

    match = match_protomers(b, a)
    start = np.vstack(a)
    end_unfitted = np.vstack([b[i] for i in match.order])
    end, err = superpose(end_unfitted, start)

    info = {
        "n_common_residues": len(common),
        "protomer_order": match.order,
        "handedness_flipped": match.handedness_flipped,
        "protomer_rmsd_by_order": match.all_rmsd,
        "endpoint_rmsd": err,
    }
    return start, end, common, info


# --------------------------------------------------------------------------
# Interpolation methods
# --------------------------------------------------------------------------

def interpolate_linear(start: np.ndarray, end: np.ndarray,
                       n_frames: int = 30) -> np.ndarray:
    t = np.linspace(0.0, 1.0, n_frames)[:, None, None]
    return (1.0 - t) * start[None] + t * end[None]


def _chain_bonds(n_sites: int, n_protomers: int = 3) -> np.ndarray:
    """Consecutive C-alpha pairs within each protomer block."""
    per = n_sites // n_protomers
    pairs = []
    for p in range(n_protomers):
        base = p * per
        idx = np.arange(base, base + per - 1)
        pairs.append(np.stack([idx, idx + 1], axis=1))
    return np.vstack(pairs)


def _bond_error(coords: np.ndarray, pairs: np.ndarray,
                target: np.ndarray) -> np.ndarray:
    d = np.linalg.norm(coords[pairs[:, 1]] - coords[pairs[:, 0]], axis=1)
    return np.abs(d - target)


def restrained_morph(start: np.ndarray, end: np.ndarray, n_frames: int = 30,
                     iterations: int = 60, n_protomers: int = 3,
                     max_bond_jump: float = 6.0) -> tuple[np.ndarray, np.ndarray]:
    """Linear interpolation with C-alpha–C-alpha distances restored.

    Each frame's target bond lengths are themselves interpolated between the
    two endpoints, then a symmetric shake-style projection pulls the frame onto
    them. Long "bonds" spanning unmodelled gaps are excluded, because forcing a
    chain break to a fixed length would drag two genuinely distant parts of the
    model together.
    """
    pairs = _chain_bonds(len(start), n_protomers)
    d0 = np.linalg.norm(start[pairs[:, 1]] - start[pairs[:, 0]], axis=1)
    d1 = np.linalg.norm(end[pairs[:, 1]] - end[pairs[:, 0]], axis=1)
    # Keep only genuine peptide bonds: ~3.8 A in both endpoints.
    real = (d0 < max_bond_jump) & (d1 < max_bond_jump)
    pairs, d0, d1 = pairs[real], d0[real], d1[real]

    frames = interpolate_linear(start, end, n_frames)
    errors = np.zeros(n_frames)

    for f in range(n_frames):
        t = f / max(n_frames - 1, 1)
        target = (1.0 - t) * d0 + t * d1
        x = frames[f].copy()
        if 0 < f < n_frames - 1:      # endpoints are experimental; leave them
            for _ in range(iterations):
                delta = x[pairs[:, 1]] - x[pairs[:, 0]]
                d = np.linalg.norm(delta, axis=1)
                bad = d > 1e-6
                corr = np.zeros_like(delta)
                corr[bad] = (delta[bad] / d[bad, None]) * \
                            (0.5 * (d[bad] - target[bad]))[:, None]
                # Symmetric: each partner moves half way.
                np.add.at(x, pairs[:, 0], corr)
                np.subtract.at(x, pairs[:, 1], corr)
            frames[f] = x
        errors[f] = _bond_error(frames[f], pairs, target).max()
    return frames, errors


def modal_morph(start: np.ndarray, end: np.ndarray, modes, n_frames: int = 30,
                n_use: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Follow the elastic-network subspace between the two endpoints.

    The endpoint difference is projected onto the low-frequency modes of the
    start structure; the path then moves only along directions the network
    supports. What is *not* captured by the retained modes is reported, so the
    user can see how much of the real change the subspace explains.
    """
    if modes.n_sites != len(start):
        raise ValueError(
            f"mode set covers {modes.n_sites} sites but the morph endpoints "
            f"have {len(start)}. The elastic network must be built on the same "
            f"reduced residue set as the endpoints — use the `start` array "
            f"returned by prepare_endpoints(), not the full structure.")
    disp = (end - start).reshape(-1)
    k = modes.n_modes if n_use is None else min(n_use, modes.n_modes)
    basis = modes.vectors[:k].reshape(k, -1)
    basis = basis / np.linalg.norm(basis, axis=1, keepdims=True)

    coeff = basis @ disp
    projected = (coeff @ basis).reshape(start.shape)
    captured = np.linalg.norm(projected) / max(np.linalg.norm(disp), 1e-9)

    t = np.linspace(0.0, 1.0, n_frames)[:, None, None]
    frames = start[None] + t * projected[None]
    pairs = _chain_bonds(len(start))
    d0 = np.linalg.norm(start[pairs[:, 1]] - start[pairs[:, 0]], axis=1)
    real = d0 < 6.0
    errors = np.array([_bond_error(f, pairs[real], d0[real]).max()
                       for f in frames])
    return frames, errors, captured


# --------------------------------------------------------------------------
# Front door
# --------------------------------------------------------------------------

def morph(start: np.ndarray, end: np.ndarray, n_frames: int = 30,
          method: str = "restrained", modes=None,
          n_protomers: int = 3) -> MorphTrajectory:
    """Build a :class:`MorphTrajectory` between two superposed coordinate sets."""
    start = np.ascontiguousarray(start, dtype=np.float64)
    end = np.ascontiguousarray(end, dtype=np.float64)
    if start.shape != end.shape:
        raise ValueError(f"endpoint shapes differ: {start.shape} vs {end.shape}")

    from .superpose import rmsd
    endpoint_rmsd = rmsd(start, end)
    meta: dict = {"n_frames": n_frames}

    if method == "linear":
        frames = interpolate_linear(start, end, n_frames)
        pairs = _chain_bonds(len(start), n_protomers)
        d0 = np.linalg.norm(start[pairs[:, 1]] - start[pairs[:, 0]], axis=1)
        real = d0 < 6.0
        errors = np.array([_bond_error(f, pairs[real], d0[real]).max()
                           for f in frames])
    elif method == "restrained":
        frames, errors = restrained_morph(start, end, n_frames,
                                          n_protomers=n_protomers)
    elif method == "modal":
        if modes is None:
            raise ValueError("the modal method needs a ModeSet")
        frames, errors, captured = modal_morph(start, end, modes, n_frames)
        meta["fraction_captured_by_modes"] = float(captured)
    else:
        raise ValueError(f"unknown morph method {method!r}")

    meta["note"] = ("This is an interpolation between two observed states, not "
                    "a simulated trajectory. It shows a plausible path; it says "
                    "nothing about the energy barrier or the order of events.")
    return MorphTrajectory(frames=frames, method=method,
                           endpoint_rmsd=endpoint_rmsd,
                           bond_error=errors, meta=meta)
