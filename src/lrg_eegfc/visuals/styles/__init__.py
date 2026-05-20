"""Project-wide matplotlib style.

Single canonical mplstyle (`lrg_eegfc.mplstyle`) tuned for compact
publication figures (font sizes, tick widths, TrueType-embedded PDF
fonts). Activate it once at the top of a script with::

    from lrg_eegfc.visuals.styles import use_lrg_style
    use_lrg_style()
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

__all__ = ["LRG_STYLE_PATH", "use_lrg_style"]

LRG_STYLE_PATH: Path = Path(__file__).with_name("lrg_eegfc.mplstyle")


def use_lrg_style() -> None:
    """Activate the project's canonical mplstyle."""
    plt.style.use(str(LRG_STYLE_PATH))
