"""Dark theme for the application.

Molecular viewers are looked at for hours in dim rooms, and a dark viewport
makes depth cueing and rim lighting legible. The chrome matches it so the eye
is not dragged away from the 3D view.
"""

from __future__ import annotations

from PyQt6.QtGui import QColor, QPalette

__all__ = ["apply_dark_theme", "STYLESHEET"]

BG = "#0e1118"
PANEL = "#161a24"
BORDER = "#262c3a"
TEXT = "#c8cfdc"
MUTED = "#7f8798"
ACCENT = "#5b8def"

STYLESHEET = f"""
QMainWindow, QWidget {{ background: {BG}; color: {TEXT};
    font-size: 12px; }}
QDockWidget {{ color: {TEXT}; titlebar-close-icon: none; }}
QDockWidget::title {{ background: {PANEL}; padding: 6px 9px;
    border-bottom: 1px solid {BORDER}; font-weight: 600; }}
QGroupBox {{ border: 1px solid {BORDER}; border-radius: 6px;
    margin-top: 14px; padding-top: 10px; background: {PANEL}; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px;
    color: {MUTED}; font-weight: 600; text-transform: uppercase;
    font-size: 10px; letter-spacing: 0.6px; }}
QPushButton {{ background: #1e2534; border: 1px solid {BORDER};
    border-radius: 5px; padding: 6px 12px; color: {TEXT}; }}
QPushButton:hover {{ background: #26304a; border-color: {ACCENT}; }}
QPushButton:pressed, QPushButton:checked {{ background: {ACCENT};
    color: #0b0e14; border-color: {ACCENT}; }}
QPushButton:disabled {{ color: #4d5464; background: #171b25; }}
QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox {{ background: #12161f;
    border: 1px solid {BORDER}; border-radius: 5px; padding: 4px 7px;
    color: {TEXT}; selection-background-color: {ACCENT}; }}
QComboBox:hover, QLineEdit:focus {{ border-color: {ACCENT}; }}
QComboBox QAbstractItemView {{ background: {PANEL}; border: 1px solid {BORDER};
    selection-background-color: {ACCENT}; selection-color: #0b0e14; }}
QListWidget, QTableWidget {{ background: #12161f; border: 1px solid {BORDER};
    border-radius: 5px; alternate-background-color: #141924;
    selection-background-color: #1f3358; }}
QListWidget::item {{ padding: 4px 6px; }}
QListWidget::item:selected, QTableWidget::item:selected {{
    background: #1f3358; color: #ffffff; }}
QHeaderView::section {{ background: {PANEL}; color: {MUTED};
    border: none; border-bottom: 1px solid {BORDER}; padding: 5px; }}
QTabWidget::pane {{ border: 1px solid {BORDER}; border-radius: 5px;
    background: {PANEL}; }}
QTabBar::tab {{ background: transparent; color: {MUTED};
    padding: 6px 14px; border-bottom: 2px solid transparent; }}
QTabBar::tab:selected {{ color: {TEXT}; border-bottom-color: {ACCENT}; }}
QSlider::groove:horizontal {{ height: 4px; background: {BORDER};
    border-radius: 2px; }}
QSlider::handle:horizontal {{ background: {ACCENT}; width: 13px;
    margin: -5px 0; border-radius: 7px; }}
QProgressBar {{ border: 1px solid {BORDER}; border-radius: 4px;
    background: #12161f; height: 6px; }}
QProgressBar::chunk {{ background: {ACCENT}; border-radius: 3px; }}
QStatusBar {{ background: {PANEL}; border-top: 1px solid {BORDER};
    color: {TEXT}; }}
QMenuBar {{ background: {PANEL}; border-bottom: 1px solid {BORDER}; }}
QMenuBar::item:selected {{ background: {ACCENT}; color: #0b0e14; }}
QMenu {{ background: {PANEL}; border: 1px solid {BORDER}; }}
QMenu::item:selected {{ background: {ACCENT}; color: #0b0e14; }}
QCheckBox::indicator {{ width: 14px; height: 14px; border-radius: 3px;
    border: 1px solid {BORDER}; background: #12161f; }}
QCheckBox::indicator:checked {{ background: {ACCENT}; border-color: {ACCENT}; }}
QScrollBar:vertical {{ background: transparent; width: 10px; }}
QScrollBar::handle:vertical {{ background: #2b3243; border-radius: 5px;
    min-height: 24px; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
QToolTip {{ background: {PANEL}; color: {TEXT}; border: 1px solid {ACCENT};
    padding: 4px; }}
"""


def apply_dark_theme(app) -> None:
    app.setStyle("Fusion")
    pal = QPalette()
    pal.setColor(QPalette.ColorRole.Window, QColor(BG))
    pal.setColor(QPalette.ColorRole.WindowText, QColor(TEXT))
    pal.setColor(QPalette.ColorRole.Base, QColor("#12161f"))
    pal.setColor(QPalette.ColorRole.AlternateBase, QColor(PANEL))
    pal.setColor(QPalette.ColorRole.Text, QColor(TEXT))
    pal.setColor(QPalette.ColorRole.Button, QColor(PANEL))
    pal.setColor(QPalette.ColorRole.ButtonText, QColor(TEXT))
    pal.setColor(QPalette.ColorRole.Highlight, QColor(ACCENT))
    pal.setColor(QPalette.ColorRole.HighlightedText, QColor("#0b0e14"))
    pal.setColor(QPalette.ColorRole.ToolTipBase, QColor(PANEL))
    pal.setColor(QPalette.ColorRole.ToolTipText, QColor(TEXT))
    app.setPalette(pal)
    app.setStyleSheet(STYLESHEET)
