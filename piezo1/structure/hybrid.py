"""The full-length model: experimental core plus predicted distal blade.

Aim A1 asks for the best available full-length PIEZO1 by combining cryo-EM with
prediction for the unresolved blade, and — the part that matters — **keeping the
seam honest and visible rather than hidden**. A spliced model that does not say
where the splice is invites every downstream number to be read as experimental.

**What is being joined.** 8YEZ resolves residues 570–2521; the AlphaFold model
covers 1–2521. So the graft is the N-terminal 569 residues of the distal blade,
which no experiment has resolved.

**And it is a weak graft, which is the honest headline.** Mean pLDDT over the
grafted region is **64.5**, with only 48% of residues above the conventional 70
confidence threshold — against 74.2 over the region the experiment already
covers. The prediction is *least* confident exactly where it is being relied on.
:class:`HybridModel` reports that rather than leaving a caller to discover it.

**Every atom carries its provenance.** ``source`` is ``"experimental"`` or
``"predicted"`` per atom, so a renderer can colour the seam, an analysis can
exclude the predicted part, and no measurement can silently average across it.
:meth:`HybridModel.experimental_only` is the selection an analysis should use
unless it has a reason not to.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["HybridModel", "build_hybrid_model", "DEFAULT_PREDICTED_MODEL"]

def _plddt_confident() -> float:
    """AlphaFold's own confident/low boundary, resolved at call time."""
    return float(_P.value("hybrid.plddt_confident"))

#: Human PIEZO1, the model the fetcher stores. Version discovered from the API,
#: never guessed — see `io.fetch.fetch_alphafold`.
DEFAULT_PREDICTED_MODEL = "AF-Q92508-F1-model_v6.cif"


@dataclass
class HybridModel:
    """An experimental structure extended by prediction, with the seam kept."""

    xyz: np.ndarray
    res_seq: np.ndarray
    source: np.ndarray                 # "experimental" | "predicted", per atom
    plddt: np.ndarray                  # NaN where experimental
    seam_residue: int | None = None
    overlap_rmsd: float = float("nan")   # A, over the seam-local anchor
    global_rmsd: float = float("nan")    # A, over the whole shared region
    meta: dict = field(default_factory=dict)

    @property
    def predicted(self) -> np.ndarray:
        return self.source == "predicted"

    @property
    def experimental_only(self) -> np.ndarray:
        """The selection an analysis should use unless it means otherwise."""
        return ~self.predicted

    @property
    def confident_prediction(self) -> np.ndarray:
        """Predicted atoms whose pLDDT clears the conventional threshold."""
        return self.predicted & (self.plddt >= _plddt_confident())

    def summary(self) -> str:
        n_pred = int(self.predicted.sum())
        fraction = float(self.confident_prediction.sum()) / max(n_pred, 1)
        return (f"seam at residue {self.seam_residue}: "
                f"{int(self.experimental_only.sum())} experimental atoms, "
                f"{n_pred} predicted ({fraction:.0%} above pLDDT "
                f"{_plddt_confident():g}); the seam fits to "
                f"{self.overlap_rmsd:.2f} A, but the two models differ by "
                f"{self.global_rmsd:.1f} A overall")

    def warnings(self) -> list:
        """What a caller must be told before using this model.

        Returned rather than logged, so a report can carry them and a test can
        assert they are present.
        """
        out = []
        n_pred = int(self.predicted.sum())
        if n_pred:
            fraction = float(self.confident_prediction.sum()) / n_pred
            if fraction < 0.7:
                out.append(
                    f"only {fraction:.0%} of the grafted region clears pLDDT "
                    f"{_plddt_confident():g}; the prediction is least confident "
                    f"where it is being relied on")
        if np.isfinite(self.global_rmsd) and self.global_rmsd > 5.0:
            out.append(
                f"the two models are not the same conformation: they differ by "
                f"{self.global_rmsd:.1f} A over the whole shared region, even "
                f"though the seam itself fits to {self.overlap_rmsd:.1f} A. The "
                f"graft is placed, not validated")
        out.append("the predicted region is a MODEL: no experiment resolves "
                   "these residues, and analyses should use "
                   "`experimental_only` unless they mean to include it")
        return out


def _ca_by_residue(structure) -> dict:
    mask = structure.mask_ca()
    return {int(r): xyz for r, xyz in
            zip(structure.res_seq[mask], structure.xyz[mask])}


def build_hybrid_model(experimental, predicted=None, chain: str | None = None,
                       anchor_window: int | None = None) -> HybridModel:
    """Graft the predicted distal blade onto one experimental protomer.

    One protomer, not the trimer: the AlphaFold model is a monomer, so grafting
    onto all three would require assuming the blade's placement relative to the
    C3 axis, which the prediction does not determine. A caller wanting a trimer
    should apply the experimental structure's own symmetry afterwards, and know
    that it is doing so.
    """
    from ..config import STRUCTURE_DIR
    from ..core.structure import Structure
    from .superpose import kabsch

    if predicted is None:
        predicted = Structure.from_file(STRUCTURE_DIR / DEFAULT_PREDICTED_MODEL)

    if chain is None:
        chains = [c for c in experimental.chains
                  if (experimental.mask_ca()
                      & (experimental.chain == c)).sum() > 300]
        if not chains:
            raise ValueError("no chain with enough C-alphas to graft onto")
        chain = chains[0]

    keep = experimental.chain == chain
    exp_res = experimental.res_seq[keep & experimental.mask_ca()]
    if len(exp_res) == 0:
        raise ValueError(f"chain {chain} has no C-alphas")
    seam = int(exp_res.min())

    # Anchor on residues NEAR the seam rather than on the whole overlap.
    # Fitting all 1279 shared residues gives 19.0 A RMSD, because the AlphaFold
    # and cryo-EM blades are different conformations of a long flexible arm —
    # and spreading that error into the join misplaces the graft. A 200-residue
    # window uses 110 residues and fits to 2.4 A. Aligning on the graft itself
    # would be circular, so only shared residues are used either way.
    # `subset`, and no fallback. The first version guarded with `hasattr` and
    # fell back to the whole structure — which for a trimer builds the residue
    # map over all three chains, keeping whichever came last, and measured the
    # overlap against a mixture. It reported a plausible 19 A.
    exp_ca = _ca_by_residue(experimental.subset(keep))
    pred_ca = _ca_by_residue(predicted)
    shared = sorted(set(exp_ca) & set(pred_ca))
    if len(shared) < 50:
        raise ValueError(f"only {len(shared)} shared residues; cannot superpose")

    if anchor_window is None:
        anchor_window = int(_P.value("hybrid.anchor_window"))
    anchor = [r for r in shared if r < seam + anchor_window] or shared

    moving = np.array([pred_ca[r] for r in anchor], dtype=float)
    target = np.array([exp_ca[r] for r in anchor], dtype=float)
    rotation, translation, centroid = kabsch(moving, target)
    fitted = (moving - centroid) @ rotation.T + translation
    rmsd = float(np.sqrt(((fitted - target) ** 2).sum(axis=1).mean()))

    # The global disagreement is reported too: it is the honest caveat that the
    # two models are not the same conformation, which a good local fit hides.
    all_moving = np.array([pred_ca[r] for r in shared], dtype=float)
    all_target = np.array([exp_ca[r] for r in shared], dtype=float)
    all_fitted = (all_moving - centroid) @ rotation.T + translation
    global_rmsd = float(
        np.sqrt(((all_fitted - all_target) ** 2).sum(axis=1).mean()))

    # Graft only what the experiment does not resolve.
    graft = predicted.res_seq < seam
    graft_xyz = ((predicted.xyz[graft] - centroid) @ rotation.T + translation)

    exp_xyz = experimental.xyz[keep]
    xyz = np.vstack([exp_xyz, graft_xyz])
    res_seq = np.concatenate([experimental.res_seq[keep],
                              predicted.res_seq[graft]])
    source = np.array(["experimental"] * int(keep.sum())
                      + ["predicted"] * int(graft.sum()))
    plddt = np.concatenate([np.full(int(keep.sum()), np.nan),
                            predicted.b_factor[graft]])

    return HybridModel(
        xyz=xyz, res_seq=res_seq, source=source, plddt=plddt,
        seam_residue=seam, overlap_rmsd=rmsd, global_rmsd=global_rmsd,
        meta={"experimental": getattr(experimental, "name", ""),
              "chain": chain, "predicted_model": DEFAULT_PREDICTED_MODEL,
              "anchor_window": anchor_window, "anchor_residues": len(anchor),
              "shared_residues": len(shared),
              "grafted_residues": int(len(set(predicted.res_seq[graft])))})
