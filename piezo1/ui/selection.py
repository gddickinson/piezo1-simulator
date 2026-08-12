"""What is selected, and what a click in the viewport does with it.

Split out of `main_window.py` when it reached the project's 500-line limit, at
the seam a user complaint had just exposed: selection, camera focus and
click-to-measure are one concern, and the window shell around them is another.

The design decision worth keeping in view is that **picking is armed**. A click
in the viewport means "tell me about this residue", and a measurement tool that
silently consumed those clicks would break inspection — so measuring is a mode
you enter. The cost is discoverability, which is paid for here in two ways: a
click always marks what it hit (before this, a click changed only the status
bar, so there was no way to tell it had registered), and the first few unarmed
clicks say where the measuring tool lives.
"""

from __future__ import annotations

import numpy as np

__all__ = ["SelectionMixin", "PICK_HINTS"]

#: How many times an unarmed click says where the measuring tool lives. The
#: status bar's job is to identify the residue clicked, so a permanent nag
#: would crowd out the thing the user asked for; but silence is what left the
#: Measure panel undiscoverable, so it is said a few times and then dropped.
PICK_HINTS = 3


class SelectionMixin:
    """Highlighting, camera focus and the pick-to-measure path."""

    # ------------------------------------------------------------ highlight

    def _highlight(self, residues, label: str, chains=None) -> None:
        """Mark residues on the model. ``chains`` restricts it to some of them.

        Annotations select a residue *number*, which in a homotrimer means all
        three copies — that is what a user asking for "the gate" wants. A click
        means one specific copy, the one under the cursor, so it passes the
        chain it landed on.
        """
        if self.view is None or self.structure is None:
            return
        self.selected_residues = [int(r) for r in (residues or [])]
        self.selection_label = label
        st = self.structure
        if not residues:
            self.view.highlight = None
        else:
            self.view.highlight = np.isin(
                st.res_seq, np.asarray(list(residues), dtype=np.int32))
            if chains is not None:
                self.view.highlight &= np.isin(st.chain, list(chains))
        # Highlight is drawn by the sphere shader, so show the selected atoms
        # as spheres on top of whatever representation is active.
        self.viewport.scene.remove(f"{self.view.name}:selection")
        # One current selection at a time: marking a residue on the primary
        # structure replaces any marked feature atom, and vice versa.
        self.viewport.scene.remove("feature:selection")
        if residues:
            mask = self.view.highlight
            n = int(mask.sum())
            if n:
                batch = self.viewport.scene.spheres(f"{self.view.name}:selection")
                batch.upload(st.xyz[mask], st.vdw_radii()[mask] * 1.05,
                             np.tile(np.array([1.0, 0.83, 0.2], np.float32), (n, 1)),
                             np.ones(n, np.float32))
                self._set_status(f"{label}: {n} atoms highlighted")
            else:
                self._set_status(f"{label}: no atoms — residues not modelled "
                                 f"in {self.record.pdb if self.record else '?'}")
        self.viewport.update()

    def _focus_residues(self, residues) -> None:
        """Move the camera to a selection, if the user has asked for that."""
        mode = self.focus_mode()
        if (mode == "none" or self.structure is None
                or self.viewport.scene is None or not residues):
            return
        mask = np.isin(self.structure.res_seq,
                       np.asarray(list(residues), dtype=np.int32))
        if not mask.any():
            return
        camera = self.viewport.scene.camera
        if mode == "frame":
            # Keep the orientation the user has set; only the pivot and the
            # distance change. Reframing rotation as well would throw away the
            # view they had chosen, which is the complaint this option exists
            # to answer.
            camera.frame(self.structure.xyz[mask])
        else:
            camera.pivot = self.structure.xyz[mask].mean(axis=0)
        self.viewport.update()

    def _set_measure_mode(self, on: bool) -> None:
        self.viewport.measure_mode = on
        self._refresh_measurements()
        self._set_status("measure mode: click atoms to pick" if on
                         else "measure mode off")

    def _refresh_measurements(self) -> None:
        """Push labels and picked-atom markers into the viewport."""
        if self.viewport.scene is None or self.structure is None:
            return
        self.viewport.set_overlay_labels(self.measure_panel.overlay_labels())
        name = "measure:picks"
        self.viewport.scene.remove(name)
        picked = self.measure_panel.highlighted_atoms()
        if picked:
            idx = np.asarray(sorted(set(picked)), dtype=int)
            idx = idx[idx < self.structure.n_atoms]
            if len(idx):
                batch = self.viewport.scene.spheres(name)
                batch.upload(self.structure.xyz[idx],
                             self.structure.vdw_radii()[idx] * 1.15,
                             np.tile(np.array([0.47, 0.78, 1.0], np.float32),
                                     (len(idx), 1)),
                             np.ones(len(idx), np.float32))
        self.viewport.update()

    def _on_pick(self, index: int) -> None:
        if index < 0 or self.structure is None:
            self._set_status("nothing under the cursor")
            return
        st = self.structure
        if self.measure_panel.armed:
            label = (f"{st.res_name[index]}{int(st.res_seq[index])}"
                     f"{st.chain[index]}.{st.atom_name[index]}")
            self.measure_panel.add_pick(int(index), st.xyz[index].astype(float),
                                        label)
            return
        res = int(st.res_seq[index])
        chain = str(st.chain[index])
        bits = [f"{st.res_name[index]}{res} chain {chain} "
                f"atom {st.atom_name[index]}"]
        if st.hetero[index]:
            # The curated domains, sites and variants are protein annotation,
            # looked up by residue number. A lipid's author-assigned number can
            # land inside the protein's range, so looking it up would confidently
            # name a domain the lipid is not part of.
            bits.append("HETATM — a resolved ligand, lipid, glycan or ion; "
                        "modelled density, not a docked pose")
        else:
            info = self.annotations.annotate_residue(res)
            if info["domain"]:
                bits.append(f"domain: {info['domain']}")
            if info["groups"]:
                bits.append("sites: " + ", ".join(info["groups"]))
            for v in info["variants"]:
                bits.append(f"variant {v['label']} ({v['classification']})")

        # Mark it on the model. Without this a click changed only the status
        # bar, so there was no way to tell which atom had been hit — or that
        # the click had registered at all.
        self._highlight([res], f"{st.res_name[index]}{res}{chain}",
                        chains=[chain])
        self._pick_hints += 1
        if not self.measure_panel.armed and self._pick_hints <= PICK_HINTS:
            bits.append("to measure, press Start picking in the Measure panel")
        self._set_status("   ·   ".join(bits))

    # -------------------------------------------------------- feature picks

    #: Marker drawn on a picked feature atom, Angstrom. Fixed rather than van
    #: der Waals because a feature source is bare coordinates — the tag's
    #: sphere centre has no element to look a radius up by.
    FEATURE_MARK_RADIUS = 2.0

    def _feature_registry(self) -> dict:
        # Lazy, so the mixin needs nothing from __init__.
        if not hasattr(self, "_pick_feature_map"):
            self._pick_feature_map: dict = {}
        return self._pick_feature_map

    def register_pick_feature(self, name: str, coords, describe) -> None:
        """Make a drawn feature's atoms answer clicks like the model's do.

        ``describe(index) -> str`` must say **what the thing is**, not only
        which atom — a modelled tag, an extra structure, a predicted residue —
        because the click's whole job is identification, and a feature
        identified as if it were the loaded structure would be the confident
        wrong answer this project spends most of its guards on.
        """
        coords = np.asarray(coords, np.float64).reshape(-1, 3)
        self._feature_registry()[name] = (coords, describe)
        self.viewport.set_feature_pick_source(name, coords)

    def unregister_pick_feature(self, name: str) -> None:
        """Drop a feature's atoms from picking; safe to call when absent."""
        self._feature_registry().pop(name, None)
        self.viewport.set_feature_pick_source(name, None)
        scene = self.viewport.scene
        if scene is not None:
            # Its marker may be the current selection; a marker outliving the
            # thing it marked would float in empty space.
            scene.remove("feature:selection")

    def _on_feature_pick(self, name: str, index: int) -> None:
        entry = self._feature_registry().get(name)
        if entry is None:
            return
        coords, describe = entry
        if not (0 <= index < len(coords)):
            return
        text = str(describe(int(index)))
        if self.measure_panel.armed:
            # Refuse out loud rather than swallow the click: measurements are
            # atom indices into the primary structure, and a distance to a
            # modelled position would be a measurement of a guess.
            self._set_status("measurements are taken on the primary structure"
                             f"   ·   this click hit {text}")
            return
        # Mark it, exactly as a primary click marks its residue — replacing
        # the previous selection, whichever kind it was.
        self._highlight([], "")
        scene = self.viewport.scene
        if scene is not None:
            batch = scene.spheres("feature:selection")
            batch.upload(coords[[index]].astype(np.float32),
                         np.full(1, self.FEATURE_MARK_RADIUS, np.float32),
                         np.array([[1.0, 0.83, 0.2]], np.float32),
                         np.ones(1, np.float32))
        self._set_status(text)
        self.viewport.update()
