#!/usr/bin/env python3
r"""Talk slide 06 — "from detecting a change to characterizing it".

ONE conceptual figure, the thesis of the slide in a single image: a single-scale
(pairwise) read gives a SCALAR — one number per process — so several genuinely
different processes collapse to the SAME verdict ("something changed"). Put a SCALE
axis under them and the same three processes fan into visibly different PROFILES:
one flat (scale-invariant — no characteristic scale), one peaked (a single
characteristic scale), one rising (emergent only after you coarse-grain). Detection
cannot tell them apart; characterization can.

Deliberately GENERIC / illustrative — NOT our data. No band labels, no p-values, no
result previewed (this is the motivation slide; results start later). The three
archetypes are the *kinds* of thing a scale axis can distinguish, drawn so all three
share the same mean → the same single number on the left.

House rules: use_lrg_style(), transparent PDF (+ .png sibling for Canva), no suptitle,
no rasterisation, plt.close. Light deck → dark ink.

Usage:
    /home/giulio/Documents/miniconda3/envs/lapbrain/bin/python \
        scripts/07_figures/talk_fig_detect_vs_characterize.py
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyArrowPatch

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import FIGURES_ROOT

OUTDIR = FIGURES_ROOT / "talk"

INK = "#1f2a37"        # dark ink for a light deck
MUTE = "#8a929c"       # muted grey (the "one number" reference)
# three distinct, non-band, colour-blind-safe archetype colours (NOT the band palette,
# so nothing here reads as δ/α/β — this is generic).
CFLAT = "#178a80"      # teal   — scale-invariant
CPEAK = "#e08a2e"      # amber  — single characteristic scale
CRISE = "#7a5aa6"      # violet — emergent when coarse-grained
COMMON = 0.5           # the shared mean = the single scalar all three collapse to


def _archetype_curves(s: np.ndarray):
    """Three profiles that share the SAME mean (=COMMON) but different SHAPES."""
    # flat: scale-invariant
    flat = np.full_like(s, COMMON)
    # peak: a single characteristic (meso) scale — gaussian bump, re-centred to mean COMMON
    g = np.exp(-((s - 0.5) / 0.17) ** 2)
    peak = COMMON + 0.50 * (g - g.mean())
    # rise: emergent — invisible when fine, appears only after coarse-graining
    sig = 1.0 / (1.0 + np.exp(-(s - 0.52) / 0.11))
    rise = COMMON + 0.62 * (sig - sig.mean())
    return flat, peak, rise


def main() -> None:
    use_lrg_style()
    s = np.linspace(0.0, 1.0, 400)
    flat, peak, rise = _archetype_curves(s)

    fig = plt.figure(figsize=(8.6, 3.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 3.5], wspace=0.42,
                           left=0.055, right=0.985, top=0.83, bottom=0.16)
    axL = fig.add_subplot(gs[0])
    axR = fig.add_subplot(gs[1])

    # ---- LEFT: DETECT — one number, three processes coincide -------------------
    xs = [0, 1, 2]
    cols = [CFLAT, CPEAK, CRISE]
    axL.bar(xs, [COMMON] * 3, width=0.72, color=cols, edgecolor=INK,
            linewidth=1.0, zorder=3)
    axL.axhline(COMMON, color=MUTE, lw=1.0, ls=(0, (4, 3)), zorder=2)
    axL.text(1.0, COMMON + 0.055, "same score", ha="center", va="bottom",
             fontsize=9.5, style="italic", color=INK)
    axL.set_xlim(-0.7, 2.7)
    axL.set_ylim(0.0, 1.0)
    axL.set_xticks([])
    axL.set_yticks([0.0, COMMON, 1.0])
    axL.set_yticklabels(["0", "", "1"])
    axL.set_ylabel("signal (a.u.)")
    axL.set_title("DETECT\none number", fontsize=11, color=INK, pad=8, linespacing=1.25)
    axL.text(1.0, -0.11, "did it change?", ha="center", va="top", fontsize=9.0,
             color=INK, transform=axL.get_xaxis_transform())
    for sp in ("top", "right"):
        axL.spines[sp].set_visible(False)

    # ---- RIGHT: CHARACTERIZE — a profile across scales -------------------------
    # (no reference line here — the equal bars on the left already carry "one number";
    #  drawing it again would sit on top of the flat scale-invariant curve.)
    axR.plot(s, flat, color=CFLAT, lw=3.0, solid_capstyle="round", zorder=4)
    axR.plot(s, peak, color=CPEAK, lw=3.0, solid_capstyle="round", zorder=4)
    axR.plot(s, rise, color=CRISE, lw=3.0, solid_capstyle="round", zorder=4)

    axR.annotate("scale-invariant", xy=(0.20, flat[80]), xytext=(0.20, 0.72),
                 ha="center", fontsize=9.5, color=CFLAT, fontweight="bold",
                 arrowprops=dict(arrowstyle="-", color=CFLAT, lw=1.0, alpha=0.6))
    axR.annotate("a characteristic scale", xy=(0.44, peak[176]), xytext=(0.46, 0.94),
                 ha="center", fontsize=9.5, color=CPEAK, fontweight="bold",
                 arrowprops=dict(arrowstyle="-", color=CPEAK, lw=1.0, alpha=0.6))
    axR.annotate("emerges only\nwhen you zoom out", xy=(0.70, rise[280]),
                 xytext=(0.83, 0.135), ha="center", fontsize=9.5, color=CRISE,
                 fontweight="bold", linespacing=1.15,
                 arrowprops=dict(arrowstyle="-", color=CRISE, lw=1.0, alpha=0.6))

    axR.set_xlim(0.0, 1.0)
    axR.set_ylim(0.0, 1.0)
    axR.set_xticks([0.0, 1.0])
    axR.set_xticklabels(["fine", "coarse"])
    axR.set_yticks([0.0, COMMON, 1.0])
    axR.set_yticklabels(["0", "", "1"])
    axR.set_xlabel(r"coarse-graining scale  $\longrightarrow$")
    axR.set_title("CHARACTERIZE\na profile across scales", fontsize=11, color=INK,
                  pad=8, linespacing=1.25)
    for sp in ("top", "right"):
        axR.spines[sp].set_visible(False)

    # ---- the arrow between the panels: detection -> resolved across scale ------
    arr = FancyArrowPatch((0.205, 0.46), (0.275, 0.46), transform=fig.transFigure,
                          arrowstyle="-|>", mutation_scale=18, lw=2.2, color=INK)
    fig.add_artist(arr)
    fig.text(0.24, 0.52, "add a\nscale axis", ha="center", va="bottom",
             fontsize=8.5, color=INK, linespacing=1.1)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    pdf = OUTDIR / "slide06_detect_vs_characterize.pdf"
    fig.savefig(pdf, transparent=True, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(OUTDIR / "slide06_detect_vs_characterize.png", transparent=True,
                bbox_inches="tight", pad_inches=0.04, dpi=300)
    plt.close(fig)
    print(f"wrote {pdf}  (+ .png sibling)", flush=True)


if __name__ == "__main__":
    main()
