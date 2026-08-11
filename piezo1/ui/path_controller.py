"""Draw the blade-to-gate route the allostery analysis finds.

This is the picture of the project's central mechanical claim. The lever model
says force on the distal blade reaches the pore through the beam and the
anchor; `allosteric_path` returns that route as a list of residue numbers, and
a list of residue numbers is unreadable as a route. Whether it goes *through*
the beam is the entire question, and only a picture answers it directly.

**The route is the cheapest path in a correlation-weighted contact graph.**
Edges join residues within the registered contact cutoff and cost
``−log|DCC_ij|``, so a pair that moves together is cheap to cross. Nothing here
is a measured signalling pathway: the correlations come from the elastic
network's own modes, which come from the coordinates. It is a statement about
the structure, in the structure's own terms.

**A drawn line reads as unique, and it is not.** That is the one misreading
this feature can produce that a caption cannot fix, so it is measured rather
than disclaimed: the same search is re-run with the drawn route's own edges
made maximally expensive, and the status line reports what the *best remaining*
route costs. If the alternative is barely worse, the drawn route is one of
many and the line on screen is a representative rather than the answer.
Round 25 already measured the beam as a near-degenerate alternative; this puts
the number in front of whoever is looking at the picture.

**PIEZO1 only, and refused rather than approximated otherwise.** The endpoints
are annotation — blade domains and the curated hydrophobic gate — read by
residue number. A PIEZO2 entry read with PIEZO1 annotation would produce a
confident route between two arbitrary places, which is precisely the failure
`numbering_check` exists to prevent.
"""

from __future__ import annotations

import numpy as np

__all__ = ["AllostericPathController", "path_endpoints", "alternative_cost",
           "NAME", "PATH_COLOR", "WEAK_COLOR", "NODE_COLOR", "TUBE_RADIUS",
           "NODE_RADIUS", "MIN_SOURCE_SITES"]

NAME = "allosteric_path"

#: The route where its steps are well correlated.
PATH_COLOR = (0.35, 0.90, 0.75)
#: The route where a step is weakly correlated — the expensive links, which are
#: where the route is least believable.
WEAK_COLOR = (0.95, 0.40, 0.35)
#: The residues the route passes through.
NODE_COLOR = (0.98, 0.98, 0.70)

TUBE_RADIUS = 1.1
NODE_RADIUS = 1.6


#: A blade unit needs at least this many network sites to be usable as the
#: source. Below it the "most distal resolved unit" is a handful of stray
#: residues at the edge of the map rather than a structural unit.
MIN_SOURCE_SITES = 20


def path_endpoints(residues, annotations,
                   min_sites: int = MIN_SOURCE_SITES
                   ) -> tuple[list[int], list[int], str]:
    """Site indices for the source and target sets, and what the source is.

    ``residues`` is the network's residue number per site — three protomers
    tiled, so a residue number appears three times and every copy is a valid
    endpoint. Both sets come from the annotation for the numbering the file is
    actually in, which is why the caller must have established that the file is
    PIEZO1 first.

    **The source is the most distal blade unit this entry actually resolves,
    not the whole blade.** Handing every blade residue to the search makes the
    answer trivial and uninteresting: the cheapest blade-to-gate route then
    starts at whichever blade residue happens to sit nearest the pore, which on
    8YEZ is THU9 and gives a five-step hop that never goes near the beam. The
    lever claim is about force arriving from far out, so the search starts as
    far out as the coordinates allow — THU1 on a full-length model, THU4 on a
    deposited entry that resolves nothing before 570. Which unit was used is
    returned rather than assumed, because it changes what the route means.
    """
    residues = np.asarray(residues)
    blade = [d for d in annotations.domains
             if d.category == "blade" and d.start is not None]
    blade.sort(key=lambda d: d.extra.get("thu_index", 99))

    source: list[int] = []
    label = ""
    for domain in blade:
        mask = (residues >= domain.start) & (residues <= domain.end)
        if int(mask.sum()) >= min_sites:
            source = np.flatnonzero(mask).tolist()
            label = domain.name
            break

    group = annotations.group("hydrophobic_gate")
    gate_residues = sorted(group.residues) if group else []
    gate = (np.isin(residues, np.asarray(gate_residues, dtype=residues.dtype))
            if gate_residues else np.zeros(residues.shape, dtype=bool))
    return source, np.flatnonzero(gate).tolist(), label


def alternative_cost(coords, dcc, source, target, residues, sites,
                     **kwargs) -> float:
    """Cost of the best route that does not use ``sites``' own steps.

    The instrument that stops a single drawn line reading as *the* pathway.
    Returns ``inf`` when deleting those steps disconnects the two sets — the
    answer meaning "genuinely the only way across" — and a number close to the
    original cost when the route is one of many.

    The steps are **deleted from the graph**, not made expensive. The first
    version suppressed their correlations to the registry floor instead, which
    leaves an edge costing ``−log(floor)``: finite, so the search walked
    straight back over the route it was meant to avoid and no input could ever
    return "unique". The calibration in `test_ui_path` is what found it, which
    is the point of calibrating a checker on a case whose answer is known.
    """
    from ..analysis.allostery import allosteric_path

    edges = list(zip(sites, sites[1:]))
    if not edges:
        return float("nan")
    try:
        return float(allosteric_path(coords, dcc, source, target, residues,
                                     exclude=edges, **kwargs).cost)
    except ValueError:                     # nothing left connected
        return float("inf")


class AllostericPathController:
    """Owns the drawn allosteric route under View -> Allosteric path."""

    def __init__(self, window) -> None:
        self.win = window
        self.result: dict | None = None
        self._names: list[str] = []
        self.pending = False

    # ----------------------------------------------------------------- state

    @property
    def visible(self) -> bool:
        return bool(self._names)

    def show(self, on: bool) -> None:
        self.pending = False
        if not on:
            self.clear()
            return
        if self.win.structure is None or self.win.viewport.scene is None:
            self.win._set_status("load a structure first")
            return
        refusal = self.refusal()
        if refusal:
            self.win._set_status(refusal)
            return
        if self.result is None:
            self.pending = True
            self.win._set_status("finding the blade-to-gate route…")
            self.win.analysis.compute_path()
            return
        self._draw()

    def annotations(self):
        """Annotation in the numbering the loaded file is actually in.

        Deliberately not `win.annotations`, which is loaded once as human and
        is what the context menu reads. A mouse entry's residue 1000 is not
        human residue 1000, and the endpoints of this route are annotation.
        """
        from ..core.annotations import load_annotations

        record = self.win.record
        species = record.numbering_species if record else "human"
        return load_annotations(species)

    def refusal(self) -> str:
        """Why the route cannot be drawn on this structure, or an empty string.

        Two refusals, and neither degrades into an approximation. Without
        modes there is no correlation matrix at all. With the wrong protein
        the endpoints are annotation read at residue numbers that mean
        something else, and the route between them would look exactly as
        convincing as a real one.
        """
        record = self.win.record
        if record is not None and getattr(record, "protein", "PIEZO1") not in ("PIEZO1", "unknown"):
            return (f"allosteric path needs PIEZO1 annotation; this entry is "
                    f"{record.protein}, whose residue numbers point at "
                    f"different residues")
        if self.win.modes is None:
            return ("allosteric path needs normal modes — compute them in the "
                    "Physics panel first")
        return ""

    def clear(self) -> None:
        scene = self.win.viewport.scene
        if scene is not None:
            for name in self._names:
                scene.remove(name)
        self._names = []
        self.win.viewport.update()

    def refresh(self, result: dict | None = None) -> None:
        """Draw a freshly computed route. Called when the path run finishes."""
        if result is not None:
            self.result = result
        if not (self.pending or self.visible):
            return
        self.pending = False
        if self.result is None or self.win.viewport.scene is None:
            return
        self._draw()

    def reset(self) -> None:
        self.result = None
        self.clear()

    # -------------------------------------------------------------- building

    def step_colors(self) -> np.ndarray:
        """One colour per step, red where the step's correlation is weakest.

        The cost of an edge is ``−log|DCC|``, so a weak step is an expensive
        one and the route is least believable exactly there. Colouring the
        tube uniformly would hide the one thing worth knowing about a shortest
        path: which link it barely made.
        """
        correlations = np.asarray(self.result["correlations"], dtype=float)
        if correlations.size == 0:
            return np.zeros((0, 3), dtype=np.float32)
        # Interpolate between the weak and strong colours over the observed
        # range, so the scale is the route's own rather than an absolute cut
        # that would paint every real route one colour.
        lo, hi = float(correlations.min()), float(correlations.max())
        t = np.zeros_like(correlations) if hi <= lo else (correlations - lo) / (hi - lo)
        weak = np.asarray(WEAK_COLOR, dtype=float)
        strong = np.asarray(PATH_COLOR, dtype=float)
        return (weak[None, :] + t[:, None] * (strong - weak)[None, :]).astype(np.float32)

    def _draw(self) -> None:
        from ..render.geometry_builders import build_tube

        self.clear()
        scene = self.win.viewport.scene
        coords = np.asarray(self.result["coords"], dtype=float)
        if len(coords) < 2:
            self.win._set_status("route too short to draw")
            return

        # Per-step colours belong to edges; the tube wants one per node, so
        # each node takes the colour of the step leaving it and the last node
        # repeats the one arriving.
        step = self.step_colors()
        node_colors = np.vstack([step, step[-1:]]) if len(step) else \
            np.tile(np.array(PATH_COLOR, np.float32), (len(coords), 1))

        mesh = build_tube(coords, node_colors.astype(np.float64),
                          radius=TUBE_RADIUS, smoothing=0)
        batch = scene.mesh(f"{NAME}:tube", two_sided=True)
        batch.upload(mesh.positions, mesh.normals, mesh.colors, mesh.indices,
                     1.0)
        self._names.append(f"{NAME}:tube")

        batch = scene.spheres(f"{NAME}:nodes")
        batch.upload(coords.astype(np.float32),
                     np.full(len(coords), NODE_RADIUS, np.float32),
                     np.tile(np.array(NODE_COLOR, np.float32), (len(coords), 1)))
        self._names.append(f"{NAME}:nodes")

        self.win.viewport.update()
        self.win._set_status(self.status_line())

    # ------------------------------------------------------------- reporting

    def domains_crossed(self) -> list[str]:
        """Which annotated domains the route passes through, in order."""
        annotations = self.annotations()
        seen: list[str] = []
        for residue in self.result["residues"]:
            domain = annotations.domain_at(int(residue))
            name = domain.name if domain else "unassigned"
            if not seen or seen[-1] != name:
                seen.append(name)
        return seen

    def status_line(self) -> str:
        """The route, and the measurement that stops it reading as unique."""
        if self.result is None:
            return "no allosteric path computed"
        residues = self.result["residues"]
        cost = float(self.result["cost"])
        source = self.result.get("source_name") or "the blade"
        alternative = self.result.get("alternative_cost")
        if alternative is None or not np.isfinite(alternative):
            degeneracy = (" · no route survives when this one's own edges are "
                          "removed, so it is the only way across this network")
        else:
            degeneracy = (f" · the best route AVOIDING these edges costs "
                          f"{alternative / cost:.3f}x this one, so the line "
                          f"drawn is a REPRESENTATIVE of many near-identical "
                          f"routes rather than the pathway")
        return (f"allosteric path: {len(residues)} steps from {source} "
                f"residue {residues[0]} — the most distal blade unit this "
                f"entry resolves — to gate residue {residues[-1]}, cost "
                f"{cost:.2f} (summed -log|DCC|) · through "
                f"{' → '.join(self.domains_crossed())}{degeneracy} · the "
                f"correlations are the elastic network's, not a measured "
                f"signal: this is a statement about the structure")
