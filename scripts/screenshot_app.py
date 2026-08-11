#!/usr/bin/env python
"""Launch the GUI, drive it through a short script, and save screenshots.

Used both as a smoke test (does the whole application actually start and
render?) and to produce the figures in the README.

Usage::

    python scripts/screenshot_app.py
    python scripts/screenshot_app.py --structure 7WLT --modes
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from PyQt6.QtCore import QTimer  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.config import SETTINGS  # noqa: E402
from piezo1.ui.gl_widget import configure_surface_format  # noqa: E402
from piezo1.ui.main_window import MainWindow  # noqa: E402
from piezo1.ui.theme import apply_dark_theme  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structure", default="8YEZ")
    ap.add_argument("--outdir", type=Path, default=Path("docs/img"))
    ap.add_argument("--modes", action="store_true",
                    help="also compute normal modes and screenshot the result")
    ap.add_argument("--analysis", action="store_true",
                    help="also run the pore analysis and screenshot the panel")
    ap.add_argument("--hold", type=int, default=1200,
                    help="milliseconds between scripted steps")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    configure_surface_format(SETTINGS.render)
    app = QApplication(sys.argv)
    apply_dark_theme(app)
    win = MainWindow()
    win.show()

    steps: list = []
    failures: list[str] = []

    def shot(name: str) -> None:
        path = args.outdir / f"app_{name}.png"
        win.grab().save(str(path))
        print(f"  saved {path}")

    def step_select() -> None:
        win.structure_panel.select(args.structure)

    def step_shot_overview() -> None:
        print(f"status: {win.status_label.text()}")
        shot("overview")

    def step_dome() -> None:
        win.physics_panel.measure_dome_requested.emit()

    def step_dome_surface() -> None:
        """Draw the fitted cap and the Helfrich skirt, then photograph them.

        Added when the surface was: a picture is the only check that catches a
        sphere fitted to the wrong atoms, which is the failure the feature
        exists to expose and which no number would show.
        """
        if len(win._mode_blocks) < 3:
            print("  dome surface skipped: needs three protomers")
            return
        win.dome_surface.show(True)
        app.processEvents()
        if not win.dome_surface.visible:
            failures.append("dome surface drew nothing")
            return
        print("dome surface:", win.status_label.text()[:150])
        shot("dome_surface")
        win.dome_surface.show(False)

    def step_shot_dome() -> None:
        # A dome needs three protomers. 4RAX is a single isolated domain and
        # the AlphaFold models are monomers; reporting that as a failure would
        # make correct behaviour look broken.
        if len(win._mode_blocks) < 3:
            print("  skipped: dome needs three protomers, this entry has "
                  f"{len(win._mode_blocks)}")
            shot("dome")
            return
        text = win.physics_panel.dome_label.text()
        print("dome:", text.replace("<br>", " | ").replace("<b>", "")
              .replace("</b>", "")[:200])
        if "radius of curvature" not in text:
            failures.append("dome measurement produced no result")
        shot("dome")

    def step_variant() -> None:
        panel = win.annotation_panel
        panel.tabs.setCurrentIndex(2)
        for row in range(panel.variant_table.rowCount()):
            if panel.variant_table.item(row, 0).text().startswith("R2456"):
                panel.variant_table.selectRow(row)
                break

    def step_shot_variant() -> None:
        print("variant status:", win.status_label.text())
        shot("variant")

    def step_measure() -> None:
        """Measure the curated disulfide — where the structure actually has one.

        Three legitimate reasons a structure cannot answer, all distinguished
        from failure because conflating them would make correct behaviour look
        broken:
          * the model is backbone-only, so there are no sulfur atoms (11ZC);
          * it is a paralogue, where PIEZO1 numbering does not transfer — 6KG7
            is PIEZO2 and has isoleucine at the second position;
          * the residues are cysteine but the bond is not modelled as formed,
            which is the case in both mouse entries at 5.2 A.
        """
        import numpy as np

        from piezo1.core.sequence import human_to_mouse

        st = win.structure
        mp = win.measure_panel
        species = win.record.numbering_species if win.record else "human"
        pair = ((2411, 2415) if species == "human"
                else tuple(human_to_mouse(r) for r in (2411, 2415)))
        chain = st.chains[0]
        print(f"  disulfide residues in {species} numbering: {pair}")

        names, sulfurs = [], []
        for residue in pair:
            here = (st.chain == chain) & (st.res_seq == residue)
            names.append(next(iter(set(st.res_name[here].tolist())), "-"))
            sulfurs.append(np.flatnonzero(here & (st.atom_name == "SG")))

        if any(n != "CYS" for n in names):
            print(f"  skipped: residues are {names}, not both cysteine "
                  f"(paralogue or different numbering)")
            return
        if any(len(s) == 0 for s in sulfurs):
            print("  skipped: no SG atoms — backbone-only model")
            return

        mp.arm_button.setChecked(True)
        for indices in sulfurs:
            win._on_pick(int(indices[0]))
        if not mp.set.measurements:
            failures.append("picking two resolved SG atoms produced no result")
        else:
            value = mp.set.measurements[0].value
            bonded = 1.8 < value < 2.5
            print(f"  S-S distance {value:.2f} A "
                  f"({'disulfide' if bonded else 'not modelled as bonded'})")
            if not (1.0 < value < 12.0):
                failures.append(f"S-S distance {value:.2f} A is implausible")
            elif species == "human" and not bonded:
                failures.append(
                    f"human disulfide measured {value:.2f} A, expected ~2.04")
        mp.arm_button.setChecked(False)

    def step_modes() -> None:
        win.physics_panel.compute_modes_requested.emit(
            {"cutoff": 15.0, "spring": "inverse_square", "n_modes": 20})

    def step_shot_modes() -> None:
        print("modes status:", win.status_label.text())
        if win.modes is None:
            failures.append("normal modes were not computed")
        else:
            win.physics_panel.color_button.setChecked(True)
        shot("modes")

    def step_pore() -> None:
        win.analysis_panel.tabs.setCurrentIndex(0)
        win.analysis_panel.pore_requested.emit()

    def step_shot_pore() -> None:
        text = win.analysis_panel.pore_label.text()
        print("pore:", text.replace("<br>", " | ").replace("<b>", "")
              .replace("</b>", "")[:200])
        if "bottleneck" not in text:
            failures.append("pore profile produced no result")
        if not win.analysis_panel.pore_plot.traces:
            failures.append("pore plot has no traces")
        elif len(win.analysis_panel.pore_plot.traces) < 2:
            failures.append("hydrophobicity trace missing from the pore plot")
        shot("pore")

    def step_focus_mode() -> None:
        """The reported bug: clicking a list entry moved the whole model."""
        import numpy as np
        if win.viewport.scene is None or win.structure is None:
            failures.append("no GL scene for the focus-mode check")
            return
        camera = win.viewport.scene.camera
        residues = [int(win.structure.res_seq[100])]

        win._set_focus_mode("none")
        before = np.array(camera.pivot, dtype=float).copy()
        win._focus_residues(residues)
        still = np.allclose(camera.pivot, before)

        win._set_focus_mode("centre")
        win._focus_residues(residues)
        moved = not np.allclose(camera.pivot, before)

        win._set_focus_mode("none")
        camera.pivot = before
        print(f"focus mode: off -> view still {still}; centre -> view moved {moved}")
        if not still:
            failures.append("focus off still moved the camera")
        if not moved:
            failures.append("focus centre did not move the camera")

    def step_layout() -> None:
        docks = win.docks.docks
        docks["measure"].hide()
        docks["model"].setFloating(True)
        win._reset_layout()
        ok = docks["measure"].isVisible() and not docks["model"].isFloating()
        print(f"layout reset restored hidden and floating panels: {ok}")
        if not ok:
            failures.append("reset layout did not restore the docks")

    def step_tour() -> None:
        """Walk the tour, checking each step actually measures something."""
        win._start_tour()
        empty = []
        for i in range(win.tour_panel.steps.count()):
            win.tour_panel.steps.setCurrentRow(i)
            step = win.tour_panel.current_step()
            text = win.tour_panel.measurement.text()
            if step.measure is not None and (not text or "running" in text):
                empty.append(step.key)
        print(f"tour: {win.tour_panel.steps.count()} steps walked, "
              f"{len(empty)} without a measurement {empty}")
        if len(empty) > 4:
            failures.append(f"tour steps reported nothing: {empty}")

    def step_session() -> None:
        """Round-trip a session through the controller, without a file dialog."""
        import tempfile
        from pathlib import Path as _Path
        from piezo1.io.session import load_session, save_session
        captured = win.session._capture()
        with tempfile.TemporaryDirectory() as tmp:
            path = _Path(tmp) / "smoke_session.json"
            save_session(captured, path)
            restored = load_session(path)
        if restored.structure != captured.structure:
            failures.append("session round-trip lost the structure")
        if restored.selected_residues != captured.selected_residues:
            failures.append("session round-trip lost the selection")
        print(f"session: {restored.structure}, "
              f"{len(restored.selected_residues)} residues, "
              f"analyses {sorted(restored.analyses)}")

    steps += [step_select, step_shot_overview, step_dome, step_shot_dome,
              step_dome_surface, step_variant, step_shot_variant, step_measure]
    if args.analysis:
        steps += [step_pore, None, None, None, None, None, step_shot_pore,
                  step_focus_mode, step_layout, step_tour, step_session]
    if args.modes:
        steps += [step_modes, None, None, None, step_shot_modes]

    def run_next() -> None:
        if not steps:
            if failures:
                print("\nFAILURES:")
                for f in failures:
                    print("  -", f)
            else:
                print("\nall scripted steps completed")
            app.quit()
            return
        fn = steps.pop(0)
        if fn is not None:
            try:
                fn()
            except Exception as exc:
                failures.append(f"{fn.__name__}: {type(exc).__name__}: {exc}")
                print(f"  ! {fn.__name__} raised {type(exc).__name__}: {exc}")
        QTimer.singleShot(args.hold, run_next)

    QTimer.singleShot(1500, run_next)
    app.exec()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
