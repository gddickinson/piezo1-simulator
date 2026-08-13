"""Application launcher: argument parsing and the Qt event loop.

Split out of :mod:`piezo1.ui.main_window` when that file reached the project's
500-line limit. The window is the application's structure; this is how it gets
started, which is a separate concern and the one that changes when someone
needs a different startup geometry.
"""

from __future__ import annotations

from PyQt6.QtWidgets import QApplication

from ..config import SETTINGS
from .main_window import MainWindow
from .menus import make_settings
from .theme import apply_theme

__all__ = ["main"]


def main(argv: list[str] | None = None) -> int:
    import argparse
    import sys

    from .gl_widget import configure_surface_format

    parser = argparse.ArgumentParser(
        prog="python -m piezo1",
        description="PIEZO1 Dynamic Structural Simulator")
    parser.add_argument("--geometry", metavar="WxH",
                        help="window size, e.g. 1280x800. Defaults to the "
                             "smaller of 1680x1000 and your screen.")
    parser.add_argument("--structure", metavar="PDB",
                        help="structure to load at startup, e.g. 8YEZ")
    parser.add_argument("--maximised", "--maximized", action="store_true",
                        dest="maximised", help="start maximised")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    configure_surface_format(SETTINGS.render)
    app = QApplication(sys.argv[:1])
    # The stored choice, not the default: applied before the window exists so
    # a light-theme user does not watch the chrome flash dark while it builds.
    apply_theme(app, make_settings().value("options/ui_theme", "dark",
                                           type=str))
    win = MainWindow()

    if args.geometry:
        try:
            width, height = (int(v) for v in args.geometry.lower().split("x"))
        except ValueError:
            parser.error(f"--geometry expects WxH, got {args.geometry!r}")
        win.resize(width, height)
    if args.maximised:
        win.showMaximized()
    else:
        win.show()
    if args.structure:
        win.structure_panel.select(args.structure.upper())
    return app.exec()
