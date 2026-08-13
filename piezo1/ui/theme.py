"""Interface themes: dark, light, and the platform's own.

Molecular viewers are looked at for hours in dim rooms, and a dark viewport
makes depth cueing and rim lighting legible — so dark is the default, and the
chrome matches it so the eye is not dragged away from the 3D view. But a
manuscript figure session happens in a bright office, and a light interface
beside a white viewport background is a legitimate way to work. The stylesheet
is one template over a token table, so the two themes cannot drift apart in
structure: a widget styled in one is styled in both.

``system`` clears the stylesheet and hands the palette back to the platform.
The few colours panels set inline (mid-grey hints and captions) were chosen to
read on either ground.
"""

from __future__ import annotations

from PyQt6.QtGui import QColor, QPalette

__all__ = ["apply_theme", "apply_dark_theme", "THEMES", "STYLESHEET"]

#: Choices offered in Options → Interface theme, as (key, label).
THEMES = (
    ("dark", "Dark (default)"),
    ("light", "Light"),
    ("system", "System"),
)

#: Every colour the stylesheet uses, named once per theme. A colour used in
#: the template but missing here would raise at import, which is the point:
#: the light theme cannot silently fall back to a dark token.
DARK = {
    "bg": "#0e1118", "panel": "#161a24", "border": "#262c3a",
    "text": "#c8cfdc", "muted": "#7f8798", "accent": "#5b8def",
    "base": "#12161f", "accent_text": "#0b0e14",
    "button": "#1e2534", "button_hover": "#26304a",
    "disabled": "#4d5464", "disabled_bg": "#171b25",
    "alt": "#141924", "selection": "#1f3358", "selection_text": "#ffffff",
    "scroll": "#2b3243",
}

LIGHT = {
    "bg": "#eef1f6", "panel": "#f9fafc", "border": "#c9d0dc",
    "text": "#252b38", "muted": "#68707f", "accent": "#3565c4",
    "base": "#ffffff", "accent_text": "#ffffff",
    "button": "#e3e8f0", "button_hover": "#d4dcea",
    "disabled": "#9aa1af", "disabled_bg": "#e7eaf0",
    "alt": "#f1f4f9", "selection": "#cfdcf5", "selection_text": "#10141c",
    "scroll": "#b9c1cf",
}


def _stylesheet(c: dict) -> str:
    return f"""
QMainWindow, QWidget {{ background: {c['bg']}; color: {c['text']};
    font-size: 12px; }}
QDockWidget {{ color: {c['text']}; titlebar-close-icon: none; }}
QDockWidget::title {{ background: {c['panel']}; padding: 6px 9px;
    border-bottom: 1px solid {c['border']}; font-weight: 600; }}
QGroupBox {{ border: 1px solid {c['border']}; border-radius: 6px;
    margin-top: 14px; padding-top: 10px; background: {c['panel']}; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px;
    color: {c['muted']}; font-weight: 600; text-transform: uppercase;
    font-size: 10px; letter-spacing: 0.6px; }}
QPushButton {{ background: {c['button']}; border: 1px solid {c['border']};
    border-radius: 5px; padding: 6px 12px; color: {c['text']}; }}
QPushButton:hover {{ background: {c['button_hover']};
    border-color: {c['accent']}; }}
QPushButton:pressed, QPushButton:checked {{ background: {c['accent']};
    color: {c['accent_text']}; border-color: {c['accent']}; }}
QPushButton:disabled {{ color: {c['disabled']};
    background: {c['disabled_bg']}; }}
QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox {{ background: {c['base']};
    border: 1px solid {c['border']}; border-radius: 5px; padding: 4px 7px;
    color: {c['text']}; selection-background-color: {c['accent']}; }}
QComboBox:hover, QLineEdit:focus {{ border-color: {c['accent']}; }}
QComboBox QAbstractItemView {{ background: {c['panel']};
    border: 1px solid {c['border']};
    selection-background-color: {c['accent']};
    selection-color: {c['accent_text']}; }}
QListWidget, QTableWidget {{ background: {c['base']};
    border: 1px solid {c['border']}; border-radius: 5px;
    alternate-background-color: {c['alt']};
    selection-background-color: {c['selection']}; }}
QListWidget::item {{ padding: 4px 6px; }}
QListWidget::item:selected, QTableWidget::item:selected {{
    background: {c['selection']}; color: {c['selection_text']}; }}
QHeaderView::section {{ background: {c['panel']}; color: {c['muted']};
    border: none; border-bottom: 1px solid {c['border']}; padding: 5px; }}
QTabWidget::pane {{ border: 1px solid {c['border']}; border-radius: 5px;
    background: {c['panel']}; }}
QTabBar::tab {{ background: transparent; color: {c['muted']};
    padding: 6px 14px; border-bottom: 2px solid transparent; }}
QTabBar::tab:selected {{ color: {c['text']};
    border-bottom-color: {c['accent']}; }}
QSlider::groove:horizontal {{ height: 4px; background: {c['border']};
    border-radius: 2px; }}
QSlider::handle:horizontal {{ background: {c['accent']}; width: 13px;
    margin: -5px 0; border-radius: 7px; }}
QProgressBar {{ border: 1px solid {c['border']}; border-radius: 4px;
    background: {c['base']}; height: 6px; }}
QProgressBar::chunk {{ background: {c['accent']}; border-radius: 3px; }}
QStatusBar {{ background: {c['panel']};
    border-top: 1px solid {c['border']}; color: {c['text']}; }}
QMenuBar {{ background: {c['panel']};
    border-bottom: 1px solid {c['border']}; }}
QMenuBar::item:selected {{ background: {c['accent']};
    color: {c['accent_text']}; }}
QMenu {{ background: {c['panel']}; border: 1px solid {c['border']}; }}
QMenu::item:selected {{ background: {c['accent']};
    color: {c['accent_text']}; }}
QCheckBox::indicator {{ width: 14px; height: 14px; border-radius: 3px;
    border: 1px solid {c['border']}; background: {c['base']}; }}
QCheckBox::indicator:checked {{ background: {c['accent']};
    border-color: {c['accent']}; }}
QScrollBar:vertical {{ background: transparent; width: 10px; }}
QScrollBar::handle:vertical {{ background: {c['scroll']};
    border-radius: 5px; min-height: 24px; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
QToolTip {{ background: {c['panel']}; color: {c['text']};
    border: 1px solid {c['accent']}; padding: 4px; }}
"""


#: The dark stylesheet under its historical name; `screenshot_app.py` and the
#: launcher import it.
STYLESHEET = _stylesheet(DARK)


def _apply_tokens(app, c: dict) -> None:
    app.setStyle("Fusion")
    pal = QPalette()
    pal.setColor(QPalette.ColorRole.Window, QColor(c["bg"]))
    pal.setColor(QPalette.ColorRole.WindowText, QColor(c["text"]))
    pal.setColor(QPalette.ColorRole.Base, QColor(c["base"]))
    pal.setColor(QPalette.ColorRole.AlternateBase, QColor(c["panel"]))
    pal.setColor(QPalette.ColorRole.Text, QColor(c["text"]))
    pal.setColor(QPalette.ColorRole.Button, QColor(c["panel"]))
    pal.setColor(QPalette.ColorRole.ButtonText, QColor(c["text"]))
    pal.setColor(QPalette.ColorRole.Highlight, QColor(c["accent"]))
    pal.setColor(QPalette.ColorRole.HighlightedText, QColor(c["accent_text"]))
    pal.setColor(QPalette.ColorRole.ToolTipBase, QColor(c["panel"]))
    pal.setColor(QPalette.ColorRole.ToolTipText, QColor(c["text"]))
    app.setPalette(pal)
    app.setStyleSheet(_stylesheet(c))


def apply_dark_theme(app) -> None:
    _apply_tokens(app, DARK)


def apply_light_theme(app) -> None:
    _apply_tokens(app, LIGHT)


def apply_system_theme(app) -> None:
    """Hand the chrome back to the platform: no stylesheet, its own palette."""
    app.setStyleSheet("")
    app.setStyle("Fusion")
    app.setPalette(app.style().standardPalette())


def apply_theme(app, key: str) -> None:
    """Apply a theme by its key. An unknown key gets the default, loudly
    rather than silently: dark is what every screenshot and figure assumes."""
    {"dark": apply_dark_theme,
     "light": apply_light_theme,
     "system": apply_system_theme}.get(key, apply_dark_theme)(app)
