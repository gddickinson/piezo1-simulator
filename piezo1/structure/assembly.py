"""A trimer built from one protomer, and everything that makes it a model.

PIEZO is a homotrimer, and half this project's machinery says so: the dome fit,
the pore profile, the elastic network and the paralogue comparison all take
three protomer blocks and refuse anything else. That refusal is right — a dome
radius measured on one blade is not a dome radius — and it has a cost that
Round 89 made concrete. **The only structural representation of a non-animal
PIEZO is an AlphaFold monomer**, so the question the family was added to ask —
is the dome a property of the fold rather than of animals? — could not be asked
at all.

This builds the missing trimer, by superposing the monomer onto each protomer
of a deposited one. Which raises the obvious objection, and the objection is
correct: **the arrangement is the template's, not a measurement of this
protein.** Everything here is shaped around not letting that be forgotten.

**And measuring it settles what the feature is worth, which is less than it
first looked.** :func:`borrowed_fraction` splits the assembled trimer's
departure from planarity into the monomer's own shape and the template's
arrangement, and the arrangement is **79% for the worm, 83% for the plant and
96% for the rat**. So a dome radius measured on an assembly is overwhelmingly a
measurement of the template, and this does **not** answer the question it was
reached for — whether the dome is a property of the fold rather than of
animals. That still needs a non-animal PIEZO trimer, and none exists.

What it does do is make the monomers *drawable and usable*: three protomers on
screen instead of a single arm, and every analysis that refuses a monomer able
to run — with each result readable against the fraction above. Round 89
recorded "cannot be asked from structure at all" as a gap in the world; this
narrows it to "can be asked, and the answer would be 83% about 9ZIS", which is
a more precise statement of the same gap rather than a way round it.

**What is measured and what is borrowed, separately.** ``structure.planarity``
already decomposes a trimer's departure from a plane into a *within-protomer*
term and an *arrangement* term. On an assembled trimer the first is the
monomer's own shape and the second is inherited wholesale, so the split is not
a caveat in prose but a number this project can compute —
:func:`borrowed_fraction` runs it and returns what share of the departure from
planarity came from the template. Anything measured on an assembly that depends
on the arrangement — the dome radius above all — should be read against it.

**Correspondence goes through an alignment, never a residue number.** The whole
point is assembling a *plant* protein onto a *mouse* template, and no two
PIEZOs in the family share a length. Residue 2447 of one is not residue 2447 of
the other, and pairing them by number would superpose the wrong helices onto
each other and still produce a confident, symmetric, entirely wrong trimer.

**Each protomer is placed independently.** The alternative — place one, then
rotate it twice about the template's C3 axis — gives a trimer whose C3 symmetry
is exactly perfect by construction, and this project measures C3 deviation and
reports it. A constructed zero would be indistinguishable on screen from a
measured one. Placing three times against the template's three real chains
leaves the template's own small asymmetry in, so ``detect_c3_axis`` keeps
returning something meaningful, and it hands over a free diagnostic: three
superposition RMSDs that must agree, because a monomer that fits one template
protomer well and another badly has not been placed, it has been forced.

**The interface is invented, so the clashes are the honest symptom.** Nothing
here models an inter-protomer contact. If the monomer's blade is shaped
differently from the template's, the assembly will bury atoms in each other,
and that count is reported rather than relaxed away — a clashing assembly is
telling you the template is wrong for this protein.

**The fit is on the rigid core, found by outlier rejection rather than named.**
The first version superposed every corresponding residue and reported a
placement RMSD of **19 Å for human PIEZO1's own AlphaFold model onto 6B3R, and
25 Å for the plant** — numbers that would have produced a visibly absurd trimer.
The cause is the one ``hybrid`` already ran into: a global fit is dominated by
the distal blade, which is exactly the part that differs most between a
prediction and a deposited entry, and the two models are known to differ by 75 Å
overall while their cores agree to 2.4 Å. So the fit iterates — superpose, drop
the worst-deviating residues, refit — which finds the rigid core without needing
annotation for it. That matters here more than anywhere: the protein this exists
to serve is a plant PIEZO for which no curated helix ranges exist, so any
core defined by naming residues could not be applied to it at all.

Both numbers are reported. The core RMSD says how well the placement is
determined; the **full** RMSD over every corresponding residue says how much of
the molecule is not following it, which is the blade, and is the reason an
assembled trimer's outer rim should not be measured.

Calibrated on the case with a known answer: a protomer *taken out of* a real
trimer and reassembled against that same trimer must reproduce it to ~0 Å, and
against a different template must not.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .clashes import count_clashes

__all__ = ["TrimerAssembly", "assemble_trimer", "best_template",
           "is_monomer", "borrowed_fraction", "ASSEMBLY_CHAINS"]

#: Chain labels the assembled protomers are given. Deliberately not A/B/C —
#: an assembled file that looks like a deposited one is the thing to avoid.
ASSEMBLY_CHAINS = ("X", "Y", "Z")


def is_monomer(structure: Structure) -> bool:
    """One well-resolved chain, so everything needing a C3 axis will refuse it."""
    from .protomers import well_resolved_chains

    return len(well_resolved_chains(structure)) < 3


@dataclass
class TrimerAssembly:
    """A modelled trimer, with what it is a model *of* attached to it."""

    structure: Structure | None
    template: str = ""
    n_corresponding: int = 0
    #: One per placed protomer, over the rigid core the fit converged on.
    #: They must agree — see the module docstring.
    placement_rmsd: tuple = ()
    #: The same fits scored over **every** corresponding residue. Much larger,
    #: and the gap is the distal blade declining to follow the core.
    full_rmsd: tuple = ()
    #: Residues the core fit kept, of those available.
    n_core: int = 0
    #: Heavy-atom pairs closer than ``assembly.clash_distance`` between
    #: different assembled protomers.
    clashes: int = 0
    alignment_identity: float = float("nan")
    refusal: str = ""

    meta: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.structure is not None and not self.refusal

    @property
    def worst_placement(self) -> float:
        return max(self.placement_rmsd) if self.placement_rmsd else float("nan")

    @property
    def placement_spread(self) -> float:
        """How differently the three template protomers fit the same monomer.

        Near zero for a template whose protomers are equivalent, which they are
        in a real trimer. A large spread means the placement is being driven by
        something other than the fold.
        """
        if len(self.placement_rmsd) < 2:
            return float("nan")
        return max(self.placement_rmsd) - min(self.placement_rmsd)

    @property
    def caveat(self) -> str:
        """The sentence that may not be omitted wherever this is shown."""
        return (
            f"MODELLED TRIMER, not a structure. One protomer placed three "
            f"times onto {self.template}'s protomers, so the ARRANGEMENT is "
            f"{self.template}'s and only the within-protomer shape is this "
            f"protein's. No inter-protomer contact is modelled; "
            f"{self.clashes:,} heavy-atom clashes remain against a "
            f"deposited trimer's 3-8. Anything measured "
            f"on this that depends on how the three protomers sit — the dome "
            f"radius, the pore, the C3 axis — is MOSTLY a measurement of "
            f"{self.template}: across the catalogue 79-96% of an assembly's "
            f"departure from planarity is the template's arrangement. Call "
            f"borrowed_fraction() for this one.")

    @property
    def worst_full(self) -> float:
        return max(self.full_rmsd) if self.full_rmsd else float("nan")

    @property
    def at_floor(self) -> bool:
        """The core stopped at its minimum, so the template barely fits.

        Reported because it is the difference between "this protomer follows
        the template over 948 residues" (rat on 6B3R) and "fewer than 200 of
        it does, and the search stopped because it was not allowed to go
        lower" (the plant on 9ZIS). The second is still drawable and should
        not be read as a model of anything.
        """
        return self.ok and self.n_core <= int(
            _P.value("assembly.min_corresponding"))

    @property
    def clash_reference(self) -> str:
        """What a real trimer scores on the same counter, for scale.

        Measured rather than asserted: 6B3R gives 8, 7WLT 3 and 9ZIS 6, so a
        four-figure count here is the interface being wrong rather than the
        counter being noisy.
        """
        return "a deposited trimer scores 3-8 on this counter"

    def summary(self) -> str:
        if not self.ok:
            return f"not assembled: {self.refusal}"
        return (f"trimer assembled on {self.template} (sequence identity "
                f"{self.alignment_identity:.3f}): {self.n_core} of "
                f"{self.n_corresponding} corresponding residues form the "
                f"rigid core, fitted to {self.worst_placement:.2f} A"
                + (f" (AT THE FLOOR — the template accounts for almost none "
                   f"of this protomer)" if self.at_floor else "")
                + (f"; over all {self.n_corresponding} it is "
                   f"{self.worst_full:.1f} A, and that gap is the distal blade"
                   if self.worst_full > self.worst_placement + 0.5 else "")
                + f". {self.clashes:,} inter-protomer clashes "
                f"({self.clash_reference})")


def borrowed_fraction(assembly: "TrimerAssembly", reference: str = "") -> dict:
    """How much of this assembly's non-planarity is the template's arrangement.

    ``structure.planarity`` splits a trimer's departure from a plane into a
    within-protomer term, which on an assembly is the monomer's own shape, and
    an arrangement term, which is inherited entirely. The ratio is the honest
    scale factor to read a dome radius against: at 0.9 the shape on screen is
    almost all template.

    Returns the reason instead of a number when it cannot be computed, which
    is most often that the protein has no annotated transmembrane band —
    :mod:`piezo1.core.annotations` holds those for human and mouse PIEZO1 only.
    """
    from ..core.numbering_check import identify_numbering
    from .planarity import planarity

    if not assembly.ok:
        return {"error": assembly.refusal}
    reference = reference or identify_numbering(assembly.structure).reference
    try:
        split = planarity(assembly.structure, reference)
    except Exception as exc:                              # noqa: BLE001
        return {"error": f"cannot decompose: {type(exc).__name__}: {exc}"}
    within, arrangement = split.protomer_rmsd, split.arrangement_rmsd
    total = (within ** 2 + arrangement ** 2) ** 0.5
    return {
        "within_protomer_A": within,
        "arrangement_A": arrangement,
        "borrowed_fraction": (arrangement / total) if total else float("nan"),
        "template": assembly.template,
        "note": ("the within-protomer term is this protein's own shape; the "
                 "arrangement term is the template's, in full"),
    }


def best_template(structure: Structure) -> str:
    """Which deposited trimer to build against.

    Same protein first, then the most residues resolved — the template supplies
    the arrangement, so the one that resolves most of the molecule constrains
    most of it. Falls back to the best-resolved PIEZO trimer of any protein,
    because a template of the wrong protein and a stated identity is more
    useful than a refusal, and the identity is what tells a reader how much to
    believe.
    """
    from ..core.numbering_check import identify_numbering
    from ..io.registry import load_registry

    try:
        protein = identify_numbering(structure).protein
    except Exception:                                     # noqa: BLE001
        protein = ""

    candidates = []
    for record in load_registry().available():
        if record.state == "predicted" or record.n_protomers < 3:
            continue
        n_ca = max((c["n_ca"] for c in record.protomer_chains), default=0)
        candidates.append((record.protein == protein, n_ca, record.pdb))
    if not candidates:
        return ""
    candidates.sort(reverse=True)
    return candidates[0][2]


def _correspondence(monomer: Structure, template: Structure):
    """``(monomer CA indices, template residue numbers, identity)``.

    Through a real global alignment of the two reference sequences, so a plant
    monomer can be placed on a mouse template. Same-numbering pairs go through
    the same path rather than a shortcut, because a shortcut is a second route
    that has to be kept in step.
    """
    # Built here from ``core`` rather than through ``analysis.paralogue``,
    # which offers the same two lines: ``structure`` importing ``analysis``
    # points the dependency arrow backwards, and `test_architecture` says so.
    from ..core.numbering_check import identify_numbering, reference_entry
    from ..core.sequence import NumberingMap

    first = identify_numbering(monomer).reference
    second = identify_numbering(template).reference
    if not first or not second:
        return None, None, float("nan")
    if first == second:
        mapping, identity = None, 1.0
    else:
        numbering = NumberingMap.from_sequences(
            reference_entry(first)["sequence"],
            reference_entry(second)["sequence"], first, second)
        mapping, identity = numbering.a_to_b, numbering.identity

    mask = monomer.mask_ca()
    if monomer.chains:
        mask = mask & (monomer.chain == monomer.chains[0])
    residues = monomer.res_seq[mask]
    paired = [(int(r), int(r) if mapping is None else mapping.get(int(r)))
              for r in residues]
    paired = [(a, b) for a, b in paired if b is not None]
    return paired, identity, mask


def assemble_trimer(monomer: Structure, template: str | None = None,
                    ) -> TrimerAssembly:
    """Place ``monomer`` onto each protomer of a deposited trimer.

    Returns a :class:`TrimerAssembly` whose ``structure`` carries all three
    copies, or one whose ``refusal`` says why it could not be built. It never
    returns a two-chain or otherwise partial assembly: every consumer of this
    checks for three protomers, and something that passes that check while
    being wrong is the failure mode worth preventing.
    """
    from ..io.registry import load_registry
    from .protomers import well_resolved_chains
    from .superpose import kabsch

    if not is_monomer(monomer):
        return TrimerAssembly(
            structure=None,
            refusal=(f"{monomer.name} already models "
                     f"{len(well_resolved_chains(monomer))} protomers; "
                     f"assembling one would replace measured coordinates with "
                     f"a model of them"))

    template = template or best_template(monomer)
    record = load_registry().get(template) if template else None
    if record is None or not record.available:
        return TrimerAssembly(structure=None,
                              refusal=f"no template trimer available "
                                      f"({template or 'none found'})")
    reference = Structure.from_file(record.path)

    paired, identity, mask = _correspondence(monomer, reference)
    if not paired:
        return TrimerAssembly(
            structure=None, template=template,
            refusal=(f"no residue correspondence between {monomer.name} and "
                     f"{template}; both numberings must be identifiable"))

    minimum = int(_P.value("assembly.min_corresponding"))
    chains = well_resolved_chains(reference)[:3]
    monomer_residues = monomer.res_seq[mask]
    monomer_xyz = monomer.xyz[mask]

    placements, rmsds, full = [], [], []
    n_core = 0
    for chain in chains:
        chain_mask = reference.mask_ca() & (reference.chain == chain)
        chain_residues = reference.res_seq[chain_mask]
        chain_xyz = reference.xyz[chain_mask]
        available = {int(r): i for i, r in enumerate(chain_residues)}

        mobile_index, target_index = [], []
        for source, target in paired:
            if target in available:
                hit = np.where(monomer_residues == source)[0]
                if hit.size:
                    mobile_index.append(int(hit[0]))
                    target_index.append(available[target])
        if len(mobile_index) < minimum:
            return TrimerAssembly(
                structure=None, template=template,
                n_corresponding=len(mobile_index),
                alignment_identity=identity,
                refusal=(f"only {len(mobile_index)} corresponding residues "
                         f"with {template} chain {chain}, below the "
                         f"{minimum} needed to place a protomer"))

        mobile = monomer_xyz[mobile_index]
        target = chain_xyz[target_index]
        rotation, translation, centroid, core, core_rmsd = _core_fit(
            mobile, target)
        fitted = (mobile - centroid) @ rotation.T + translation
        rmsds.append(core_rmsd)
        full.append(float(np.sqrt(((fitted - target) ** 2).sum(1).mean())))
        n_core = int(core.sum())
        placements.append((monomer.xyz - centroid) @ rotation.T + translation)

    assembled = _stack(monomer, placements, template)
    clashes = count_clashes(assembled)
    return TrimerAssembly(
        structure=assembled, template=template,
        n_corresponding=len(paired), placement_rmsd=tuple(rmsds),
        full_rmsd=tuple(full), n_core=n_core,
        clashes=clashes, alignment_identity=identity,
        meta={"template_chains": list(chains),
              "source": monomer.name,
              "note": "arrangement inherited from the template; only the "
                      "within-protomer shape is this protein's"})


def _core_fit(mobile: np.ndarray, target: np.ndarray):
    """Superpose, drop the worst-deviating residues, refit, repeat.

    Returns ``(rotation, translation, centroid, core mask, core RMSD)``. The
    core is found rather than named, because the protein this feature exists
    for has no curated helix ranges to name it with.
    """
    from .superpose import kabsch

    cycles = int(_P.value("assembly.refit_cycles"))
    cutoff = float(_P.value("assembly.core_cutoff"))
    floor = int(_P.value("assembly.min_corresponding"))
    core = np.ones(len(mobile), dtype=bool)
    rotation = translation = centroid = None
    fraction = float(_P.value("assembly.core_fraction"))
    for _ in range(cycles):
        rotation, translation, centroid = kabsch(mobile[core], target[core])
        fitted = (mobile - centroid) @ rotation.T + translation
        deviation = np.linalg.norm(fitted - target, axis=1)

        # Two criteria, in this order, and both are needed.
        #
        # The *distance* one is what the core should be defined by, and what
        # makes `n_core` a measurement rather than a constant: however many
        # residues genuinely follow the template, that is the core. Used alone
        # it fails at the first step — starting from a global fit at 19 A,
        # nothing at all is within 3 A, so the loop breaks immediately and
        # returns the global fit it was supposed to improve on.
        #
        # The *fraction* is only how the search descends to somewhere the
        # distance criterion can be applied. Used alone it drove the core to
        # its floor on every entry in the catalogue — 200 of 2,500 residues
        # fitted to 1.2 A, which is not a core but the 200 that agree best,
        # and always exists.
        within = deviation <= cutoff
        if within.sum() >= floor:
            proposed = within
        else:
            keep = max(int(core.sum() * fraction), floor)
            if keep >= core.sum():
                break
            cut = np.partition(deviation[core], keep - 1)[keep - 1]
            proposed = core & (deviation <= cut)
        if proposed.sum() < floor or np.array_equal(proposed, core):
            break
        core = proposed
    fitted = (mobile[core] - centroid) @ rotation.T + translation
    rmsd = float(np.sqrt(((fitted - target[core]) ** 2).sum(1).mean()))
    return rotation, translation, centroid, core, rmsd


def _stack(monomer: Structure, placements, template: str) -> Structure:
    """Three placed copies as one Structure, labelled as a model."""
    import dataclasses

    n = monomer.n_atoms
    fields = ("element", "atom_name", "res_name", "res_seq", "hetero",
              "b_factor", "occupancy", "alt_loc", "entity")
    tiled = {f: np.concatenate([getattr(monomer, f)] * 3) for f in fields}
    chain = np.concatenate([np.full(n, c) for c in ASSEMBLY_CHAINS])
    return dataclasses.replace(
        monomer, xyz=np.concatenate(placements), chain=chain,
        name=f"{monomer.name}+trimer({template})",
        res_first=None, res_atom_index=None,
        meta={**monomer.meta, "assembled_from": monomer.name,
              "assembly_template": template, "is_observed": False},
        **tiled)
