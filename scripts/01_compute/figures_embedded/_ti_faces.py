"""Shared helper — load + place the four square guest headshots into the TI cold-open
figures (relation chain + race bracket). Faces live beside the figures in the talk
folder (fb/ag/dg/sm.jpg); missing files degrade gracefully to the initial hint."""
from pathlib import Path

import numpy as np
from PIL import Image

from lrg_eegfc.config.paths import FIGURES_ROOT

FACE_FILE = {"FB": "fb.jpg", "AG": "ag.jpg", "DG": "dg.jpg", "SM": "sm.jpg"}


def load_face(init):
    """Square RGB array for a guest initial, or None if the file is absent."""
    p = Path(FIGURES_ROOT) / "talk" / FACE_FILE[init]
    if not p.exists():
        return None
    im = Image.open(p).convert("RGB")
    w, h = im.size
    s = min(w, h)                        # already square; trim any odd pixel
    l, t = (w - s) // 2, (h - s) // 2
    return np.asarray(im.crop((l, t, l + s, t + s)))


def place_face(ax, img, cx, cy, size, clip_patch, zorder=1.6, margin=0.06):
    """Draw a square face centred at (cx, cy), overfilling `size` by `margin` then
    clipping to `clip_patch` (the rounded box) so it gets rounded corners. Keep the
    face BELOW the border patch's zorder so the coloured border stays on top.

    aspect="equal" is REQUIRED: aspect="auto" makes imshow reset the axes aspect to
    "auto", which stretches the coordinate system and rectangles the (square) boxes
    and faces. "equal" keeps the axes 1:1 so a square extent + square image stays 1:1.
    """
    half = size / 2 + margin
    im = ax.imshow(img, extent=[cx - half, cx + half, cy - half, cy + half],
                   origin="upper", zorder=zorder, aspect="equal",
                   interpolation="antialiased")
    im.set_clip_path(clip_patch)
    return im
