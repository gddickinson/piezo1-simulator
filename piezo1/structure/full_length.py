"""A full-length trimer: the experimental core with a predicted blade on each protomer.

:mod:`piezo1.structure.hybrid` grafts the AlphaFold distal blade onto **one**
protomer and returns coordinates. That is the right unit for asking how well the
graft fits, and it is the wrong unit for everything else: the dome needs a
three-fold axis, the elastic network needs three identically ordered blocks, the
pore needs a lumen, and none of those exist in a monomer. So nothing in this
project could run physics on a full-length model, and the distal blade — which
is *half the protein by residue count* — was absent from every measurement.

This builds a real :class:`~piezo1.core.Structure` instead, which is the only
thing the rest of the codebase consumes. Once it exists, every analysis, every
animation and every measurement works on it with no further plumbing, because
none of them can tell the difference. That is the point, and it is also the
danger, which is why the provenance below is not optional.

**Each protomer is grafted independently, and that removes an assumption.** The
obvious construction is to graft once and replicate by the measured C3. This
does not: it runs the same seam-local fit against each protomer's own
coordinates, so no symmetry is imposed on the prediction at all. What the C3
then *becomes* is a measurement — :attr:`FullLengthModel.blade_c3_deviation`
reports how far the three independently placed blades are from being
three-fold related, and a large value means the graft is unreliable rather than
that the protein is asymmetric.

**Provenance is derived, not carried.** A boolean mask stored beside the
structure would go stale the first time anything took a subset of it. The
grafted residues are exactly those below each chain's seam — the residues the
experiment does not resolve — so :func:`predicted_mask` recomputes the answer
from the coordinates and the recorded seams every time it is asked.

**What this model is not.** It is two populations. Roughly half the residues
have never been observed by any experiment; on human PIEZO1 only about 48% of
them clear pLDDT 70, and the predicted and experimental blades differ by 75 Å
over the region where both exist. A dome measured on it covers all 38
transmembrane helices instead of 22 — which is more of the molecule and less of
the evidence, and Round 83 measured what that difference alone does to the
numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .hybrid import build_hybrid_model, predicted_model_for

__all__ = ["FullLengthModel", "build_full_length", "predicted_mask",
           "SEAM_KEY", "GAP_KEY", "is_full_length", "FILL_MODES",
           "resolved_gaps"]

#: What to splice in, as (key, label, what it does). The default is ``none``:
#: a model that silently contained prediction would be the worst outcome here,
#: so filling is always something the user chose.
#:
#: The two halves are separate because they are different claims. The blade is
#: a region no PIEZO1 structure has ever resolved, placed by a fit at one end
#: only — it is cantilevered, and the further out it goes the less it is
#: constrained. An internal gap is bracketed by resolved residues on *both*
#: sides, so its fill is interpolated rather than extrapolated and is a much
#: safer thing to look at.
FILL_MODES = (
    ("none", "Deposited only",
     "Exactly what the experiment resolved, and nothing else."),
    ("gaps", "+ AlphaFold gaps",
     "Fill the unresolved stretches *inside* the deposited range. Each is "
     "anchored on resolved residues at both ends, so it is interpolated."),
    ("blade", "+ AlphaFold blade",
     "Add the distal N-terminal blade, which no PIEZO1 structure resolves. "
     "Anchored at one end only, so it is extrapolated."),
    ("full", "+ AlphaFold (full length)",
     "Both: the complete chain, roughly half of it predicted."),
)

#: Where the filled internal gaps live on the built structure, as
#: ``chain -> [[first, last], ...]`` in residue numbers.
GAP_KEY = "full_length_gaps"

#: Where the per-chain seam residues live on the built structure. Anything that
#: needs to know which atoms are predicted derives it from this rather than
#: from a mask that a subset would silently invalidate.
SEAM_KEY = "full_length_seams"


@dataclass
class FullLengthModel:
    """A trimer the analyses can run on, and the reasons to distrust half of it."""

    structure: Structure
    seams: dict                      # chain -> first experimentally resolved residue
    blade_c3_deviation: float        # A, between independently placed blades
    seam_rmsd: dict = field(default_factory=dict)      # chain -> A
    global_rmsd: dict = field(default_factory=dict)    # chain -> A
    meta: dict = field(default_factory=dict)

    @property
    def predicted(self) -> np.ndarray:
        return predicted_mask(self.structure)

    @property
    def n_predicted_atoms(self) -> int:
        return int(self.predicted.sum())

    @property
    def n_predicted_residues(self) -> int:
        mask = self.predicted & (self.structure.atom_name == "CA")
        return int(mask.sum())

    @property
    def confident_fraction(self) -> float:
        """Share of grafted atoms clearing AlphaFold's own pLDDT threshold."""
        from .hybrid import _plddt_confident

        predicted = self.predicted
        if not predicted.any():
            return float("nan")
        return float(
            (self.structure.b_factor[predicted] >= _plddt_confident()).mean())

    def summary(self) -> str:
        from .hybrid import _plddt_confident

        return (f"{self.structure.name}: {self.n_predicted_residues} predicted "
                f"residues grafted onto {len(self.seams)} protomers "
                f"({self.confident_fraction:.0%} above pLDDT "
                f"{_plddt_confident():g}); blades are three-fold related to "
                f"{self.blade_c3_deviation:.1f} A")

    def warnings(self) -> list:
        """What must reach a user before a number from this model does."""
        out = [f"HALF THIS MODEL IS A PREDICTION: {self.n_predicted_residues} "
               f"of its residues are AlphaFold, not experiment, and every "
               f"measurement below includes them."]
        confident = self.confident_fraction
        if np.isfinite(confident) and confident < 0.7:
            from .hybrid import _plddt_confident
            out.append(f"only {confident:.0%} of the grafted atoms clear pLDDT "
                       f"{_plddt_confident():g}; the prediction is least "
                       f"confident exactly where it is being relied on")
        worst = max(self.global_rmsd.values(), default=float("nan"))
        if np.isfinite(worst) and worst > 5.0:
            out.append(f"where both models exist they differ by up to "
                       f"{worst:.0f} A, so the seam fitting well says nothing "
                       f"about the rest of the blade")
        if np.isfinite(self.blade_c3_deviation) and self.blade_c3_deviation > 5.0:
            out.append(f"the three blades were placed independently and came "
                       f"out {self.blade_c3_deviation:.1f} A from three-fold "
                       f"related; treat the symmetry of this model as an "
                       f"output, not an input")
        return out



def is_full_length(structure: Structure) -> bool:
    """Whether this structure carries a grafted prediction."""
    return bool(getattr(structure, "meta", {}).get(SEAM_KEY))


def predicted_mask(structure: Structure) -> np.ndarray:
    """Which atoms are prediction rather than experiment.

    Derived from the recorded per-chain seams and the residue numbers, so it
    survives every operation that produces a new structure from this one —
    a subset, a reframe, a morph frame. A stored mask would not.
    """
    meta = getattr(structure, "meta", {})
    out = np.zeros(structure.n_atoms, dtype=bool)
    for chain, seam in (meta.get(SEAM_KEY) or {}).items():
        out |= (structure.chain == chain) & (structure.res_seq < int(seam))
    for chain, gaps in (meta.get(GAP_KEY) or {}).items():
        for first, last in gaps:
            out |= ((structure.chain == chain)
                    & (structure.res_seq >= int(first))
                    & (structure.res_seq <= int(last)))
    return out


def resolved_gaps(structure: Structure, chain: str) -> list:
    """Unresolved stretches inside one chain's resolved range, as (first, last).

    Every deposited PIEZO structure has them — 7WLT resolves 1,353 residues
    across a 784-2547 span in twelve pieces — and nothing could see them until
    the gap filling needed them. They are the safer half of the model to fill, because unlike the
    blade they have resolved residues on both sides to be anchored against.
    """
    mask = structure.mask_ca() & (structure.chain == chain)
    numbers = np.unique(structure.res_seq[mask])
    if len(numbers) < 2:
        return []
    breaks = np.flatnonzero(np.diff(numbers) > 1)
    return [(int(numbers[i]) + 1, int(numbers[i + 1]) - 1) for i in breaks]


def _concat(parts: list[Structure], name: str, meta: dict) -> Structure:
    """One structure from several, ordered as a deposited file would be.

    The sort is not tidiness. Atoms arrive experiment-first and then blade and
    gap fills, so residue numbers within a chain come out unsorted — and
    ``protomer_blocks`` locates residues with ``searchsorted``, which on
    unsorted input does not fail, it returns the wrong atoms. Here it happened
    to raise; a slightly different residue set would have produced a silently
    scrambled protomer and an elastic network built from it.

    So the built model is ordered by chain, then residue, then arrival — which
    makes it indistinguishable from a deposited file to everything downstream,
    which is the whole requirement.
    """
    arrays = {f: np.concatenate([getattr(p, f) for p in parts])
              for f in Structure._ARRAY_FIELDS}
    order = np.lexsort((np.arange(len(arrays["res_seq"])), arrays["res_seq"],
                        arrays["chain"]))
    arrays = {f: v[order] for f, v in arrays.items()}
    built = Structure(name=name, source="experimental+predicted", meta=meta,
                      **arrays)
    built._build_residue_index()
    return built


def _fill_one_gap(experimental: Structure, predicted: Structure, chain: str,
                  first: int, last: int, window: int):
    """Place the predicted residues of one gap, anchored on both flanks.

    Superposed on resolved residues *either side* of the gap, which is what
    makes an internal fill different in kind from the blade: the fit is
    constrained at both ends, so an error in it shows up as a bad flank RMSD
    rather than as a lever arm. Returns ``(atoms, rmsd)`` or ``None`` when
    there is not enough on both sides to anchor against.
    """
    from .superpose import kabsch

    exp_ca = experimental.mask_ca() & (experimental.chain == chain)
    exp_res = experimental.res_seq[exp_ca]
    pred_ca = predicted.mask_ca()
    pred_res = predicted.res_seq[pred_ca]

    flank = ([r for r in range(first - window, first) if r in set(exp_res.tolist())]
             + [r for r in range(last + 1, last + 1 + window)
                if r in set(exp_res.tolist())])
    shared = [r for r in flank if r in set(pred_res.tolist())]
    if len(shared) < 6 or not any(r < first for r in shared) \
            or not any(r > last for r in shared):
        return None

    exp_index = {int(r): i for i, r in enumerate(exp_res)}
    pred_index = {int(r): i for i, r in enumerate(pred_res)}
    target = experimental.xyz[exp_ca][[exp_index[r] for r in shared]].astype(np.float64)
    moving = predicted.xyz[pred_ca][[pred_index[r] for r in shared]].astype(np.float64)

    rotation, translation, centroid = kabsch(moving, target)
    fitted = (moving - centroid) @ rotation.T + translation
    rmsd = float(np.sqrt(((fitted - target) ** 2).sum(axis=1).mean()))

    inside = (predicted.res_seq >= first) & (predicted.res_seq <= last)
    if not inside.any():
        return None
    atoms = predicted.subset(inside)
    moved = (atoms.xyz.astype(np.float64) - centroid) @ rotation.T + translation
    atoms = atoms.copy_with_coords(moved)
    atoms.chain = np.full(atoms.n_atoms, chain, dtype=experimental.chain.dtype)
    atoms._build_residue_index()
    return atoms, rmsd


def build_full_length(experimental: Structure, mode: str = "full",
                      predicted: Structure | None = None) -> FullLengthModel:
    """Splice prediction into every protomer and return one structure.

    ``mode`` is one of :data:`FILL_MODES`. Raises for anything the graft cannot
    honestly be done on — a PIEZO2 entry, or a fragment with no protomer —
    rather than returning a partial model, because a full-length model missing
    a blade looks exactly like one that has them all.
    """
    from ..config import STRUCTURE_DIR
    from .protomers import well_resolved_chains

    if mode not in dict((k, v) for k, v, _ in
                        ((a, b, c) for a, b, c in FILL_MODES)):
        raise ValueError(f"mode must be one of {[k for k, _, _ in FILL_MODES]}")
    if predicted is None:
        predicted = Structure.from_file(
            STRUCTURE_DIR / predicted_model_for(experimental))

    chains = list(well_resolved_chains(experimental))
    if not chains:
        raise ValueError("no protomer to graft onto")

    window = int(_P.value("full_length.gap_anchor_window"))
    longest = int(_P.value("full_length.max_gap"))

    parts = [experimental]
    seams: dict = {}
    seam_rmsd: dict = {}
    global_rmsd: dict = {}
    gaps: dict = {}
    gap_rmsd: list = []
    skipped: list = []
    blades: list = []

    for chain in chains:
        if mode in ("blade", "full"):
            model = build_hybrid_model(experimental, predicted, chain=chain)
            rotation, translation, centroid = model.meta["transform"]
            seam = int(model.seam_residue)
            blade = predicted.subset(predicted.res_seq < seam)
            moved = ((blade.xyz.astype(np.float64) - centroid)
                     @ rotation.T + translation)
            blade = blade.copy_with_coords(moved)
            blade.chain = np.full(blade.n_atoms, chain,
                                  dtype=experimental.chain.dtype)
            blade._build_residue_index()
            parts.append(blade)
            blades.append(blade)
            seams[str(chain)] = seam
            seam_rmsd[str(chain)] = float(model.overlap_rmsd)
            global_rmsd[str(chain)] = float(model.global_rmsd)

        if mode in ("gaps", "full"):
            filled = []
            for first, last in resolved_gaps(experimental, chain):
                if last - first + 1 > longest:
                    skipped.append((str(chain), first, last))
                    continue
                got = _fill_one_gap(experimental, predicted, chain, first,
                                    last, window)
                if got is None:
                    skipped.append((str(chain), first, last))
                    continue
                atoms, rmsd = got
                parts.append(atoms)
                filled.append([first, last])
                gap_rmsd.append(rmsd)
            gaps[str(chain)] = filled

    deviation = _blade_c3_deviation(experimental, blades)
    suffix = {"none": "", "gaps": "+gaps", "blade": "+AF", "full": "+AF-full"}
    name = f"{experimental.name}{suffix[mode]}"
    meta = dict(experimental.meta)
    meta.update({SEAM_KEY: seams, GAP_KEY: gaps, "full_length": mode != "none",
                 "fill_mode": mode,
                 "experimental_entry": experimental.name,
                 "predicted_model": predicted.name,
                 "blade_c3_deviation_A": deviation})
    built = _concat(parts, name, meta)

    return FullLengthModel(
        structure=built, seams=seams, blade_c3_deviation=deviation,
        seam_rmsd=seam_rmsd, global_rmsd=global_rmsd,
        meta={"n_chains": len(chains), "predicted_model": predicted.name,
              "mode": mode,
              "gaps_filled": sum(len(v) for v in gaps.values()),
              "gaps_skipped": len(skipped), "skipped": skipped,
              "gap_flank_rmsd_max": (max(gap_rmsd) if gap_rmsd else float("nan")),
              "note": "two populations in one file; predicted atoms carry "
                      "pLDDT in the B-factor column"})


def _blade_c3_deviation(experimental: Structure, blades: list) -> float:
    """How far the independently placed blades are from three-fold related.

    Each blade was fitted to its own protomer with no symmetry imposed, so this
    is a *result*: rotate the first blade onto the second by the experimental
    core's own C3 and measure what is left. A small number says the prediction
    and the three protomers agree; a large one says the graft is unreliable and
    the model should not be used for anything symmetry-dependent.
    """
    from .protomers import protomer_blocks
    from .superpose import detect_c3_axis, rotation_matrix

    if len(blades) < 2:
        return float("nan")
    blocks, _ = protomer_blocks(experimental)
    if not blocks:
        return float("nan")
    axis = detect_c3_axis(blocks)
    direction = axis.direction / np.linalg.norm(axis.direction)

    first = blades[0]
    common = min(b.n_atoms for b in blades)
    worst = 0.0
    for turn, other in enumerate(blades[1:], start=1):
        # Both senses, because which way round the ring the chains run is not
        # something to assume — this project has been wrong about it three
        # times, and the blade is the part furthest from the axis.
        candidates = []
        for angle in (2.0 * np.pi * turn / len(blades),
                      -2.0 * np.pi * turn / len(blades)):
            rotated = ((first.xyz[:common].astype(np.float64) - axis.point)
                       @ rotation_matrix(direction, angle).T + axis.point)
            candidates.append(float(np.sqrt(
                ((rotated - other.xyz[:common]) ** 2).sum(axis=1).mean())))
        worst = max(worst, min(candidates))
    return worst
