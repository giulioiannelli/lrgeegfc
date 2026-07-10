"""Project-wide matplotlib style and the canonical band colour palette.

Activate the style once at the top of a figure script::

    from lrg_eegfc.visuals.styles import use_lrg_style, band_color
    use_lrg_style()

**Band colours are a spectrum.** Slow bands are red, fast bands are blue, with
the rainbow in between (``delta`` -> red, ..., ``high_gamma`` -> blue) — the
visible-light analogy, where low frequency reads as red light. NEVER hardcode a
band colour inside a figure: call :func:`band_color` (or index
:data:`BAND_COLORS`) so the entire codebase restyles from the single constant
:data:`BAND_CMAP_SPAN` / :data:`BAND_CMAP_NAME` below.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, to_hex, to_rgb

from ...config.const import BRAIN_BANDS_NAMES

__all__ = [
    "LRG_STYLE_PATH",
    "use_lrg_style",
    "BAND_CMAP_NAME",
    "BAND_CMAP_SPAN",
    "BAND_COLORS",
    "band_color",
    "band_colors_list",
    "band_cmap",
    "BAND_MARKERS",
    "band_marker",
]

LRG_STYLE_PATH: Path = Path(__file__).with_name("lrg_eegfc.mplstyle")


def use_lrg_style() -> None:
    """Activate the project's canonical mplstyle."""
    plt.style.use(str(LRG_STYLE_PATH))


# --------------------------------------------------------------------------- #
# Canonical band palette — the ONE swap point for every band-coloured figure.
#
# ``turbo`` sampled from its warm end to its cool end is a perceptual rainbow
# that reads as a frequency spectrum: delta(slow) = red -> high_gamma(fast) =
# blue, green/yellow in the middle. To restyle every band figure at once,
# change ``BAND_CMAP_NAME`` or ``BAND_CMAP_SPAN`` here and nothing else.
# --------------------------------------------------------------------------- #
BAND_CMAP_NAME: str = "turbo"
#: (slow-band sample, fast-band sample) in the colormap's [0, 1] domain.
#: turbo(0.92) is red (delta) and turbo(0.08) is blue (high_gamma).
BAND_CMAP_SPAN: Tuple[float, float] = (0.92, 0.08)


def _build_band_colors() -> Dict[str, str]:
    cmap = plt.get_cmap(BAND_CMAP_NAME)
    lo, hi = BAND_CMAP_SPAN
    n = len(BRAIN_BANDS_NAMES)
    return {
        band: to_hex(cmap(lo + (hi - lo) * (i / (n - 1))))
        for i, band in enumerate(BRAIN_BANDS_NAMES)
    }


#: Canonical hex colour per band (``delta`` … ``high_gamma``), slow=red→fast=blue.
BAND_COLORS: Dict[str, str] = _build_band_colors()


def band_color(band: str, shade: float = 1.0) -> str:
    """Canonical hex colour for a frequency band (``delta`` … ``high_gamma``).

    ``shade`` scales luminance for legibility control while keeping the hue:
    ``shade < 1`` darkens (use for text callouts / thin lines on white, since
    the turbo midtones — yellow-green ``alpha``, green ``beta`` — are bright),
    ``shade == 1`` returns the canonical hue (markers/fills with edges).
    """
    if shade == 1.0:
        return BAND_COLORS[band]
    r, g, b = to_rgb(BAND_COLORS[band])
    return to_hex(tuple(min(1.0, max(0.0, c * shade)) for c in (r, g, b)))


def band_colors_list() -> List[str]:
    """Band colours in canonical (slow→fast) order."""
    return [BAND_COLORS[b] for b in BRAIN_BANDS_NAMES]


def band_cmap() -> Colormap:
    """The underlying colormap the band palette samples (for continuous fills)."""
    return plt.get_cmap(BAND_CMAP_NAME)


# --------------------------------------------------------------------------- #
# Canonical per-band marker glyphs — so a band keeps ONE symbol across every
# scatter/forest/legend and the reader never confounds two bands. Filled shapes
# with distinct silhouettes that stay legible at small sizes.
# --------------------------------------------------------------------------- #
_BAND_MARKER_CYCLE: Tuple[str, ...] = ("o", "s", "^", "D", "v", "P")

#: Canonical matplotlib marker per band (``delta`` … ``high_gamma``), in the same
#: slow→fast order as :data:`BAND_COLORS`.
BAND_MARKERS: Dict[str, str] = {
    band: _BAND_MARKER_CYCLE[i % len(_BAND_MARKER_CYCLE)]
    for i, band in enumerate(BRAIN_BANDS_NAMES)
}


def band_marker(band: str) -> str:
    """Canonical marker glyph for a frequency band (``delta`` … ``high_gamma``)."""
    return BAND_MARKERS.get(band, "o")
