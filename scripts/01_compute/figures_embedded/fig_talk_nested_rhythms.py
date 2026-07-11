"""Talk slide 6 illustration — cognition across *time*-scales (nested rhythms).

A minimal, synthetic cross-frequency-coupling cartoon: a slow rhythm (theta)
whose phase organises bursts of a fast rhythm (high-gamma) on each crest. A
zoom inset makes the "multiscale" point literal — the fast rhythm is nested
inside the slow one. Paired on the slide with the grid-cell figure (cognition
across *space*-scales), it carries the "why multiscale matters" beat without any
structural anatomy: the signature *is* the cross-scale organisation.

Not a data figure — a schematic. Colours come from the project band palette
(theta = slow/warm, high_gamma = fast/cool) so the cartoon matches every other
band-coded artist in the deck.

Output: data/outputs/figures/talk/slide06_nested_rhythms.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

F_THETA = 5.0      # slow rhythm (Hz, arbitrary time base)
F_GAMMA = 42.0     # fast rhythm
SHARPNESS = 4.0    # how tightly gamma bursts hug the theta crest


def _signals(n: int = 3000):
    t = np.linspace(0.0, 1.0, n)
    phase = 2 * np.pi * F_THETA * t
    theta = np.sin(phase)
    env = np.exp(SHARPNESS * (np.sin(phase) - 1.0))       # peaks at theta crest
    gamma = 0.60 * env * np.sin(2 * np.pi * F_GAMMA * t)
    return t, theta, theta + gamma


def main() -> None:
    use_lrg_style()
    c_slow = band_color("theta")
    c_fast = band_color("high_gamma")
    t, theta, comp = _signals()

    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    ax.plot(t, theta, color=c_slow, lw=5, alpha=0.40, solid_capstyle="round")
    ax.plot(t, comp, color=c_fast, lw=1.15)
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.55, 2.0)
    ax.axis("off")
    ax.text(0.012, -1.36, r"slow rhythm  ($\theta$)", color=c_slow,
            fontsize=12, fontweight="bold")
    ax.text(0.012, 1.68, r"fast bursts  ($\gamma$)", color=c_fast,
            fontsize=12, fontweight="bold")

    # zoom inset -> the fast scale nested inside the slow one
    peaks = t[np.r_[False, (theta[1:-1] > theta[:-2]) & (theta[1:-1] > theta[2:]), False]]
    pc = peaks[np.argmin(np.abs(peaks - 0.5))]
    axi = ax.inset_axes([0.60, 0.58, 0.36, 0.42])
    axi.plot(t, comp, color=c_fast, lw=1.6)
    axi.plot(t, theta, color=c_slow, lw=3, alpha=0.40)
    axi.set_xlim(pc - 0.05, pc + 0.05)
    axi.set_ylim(-0.6, 1.9)
    axi.set_xticks([])
    axi.set_yticks([])
    for s in axi.spines.values():
        s.set_edgecolor("#9aa4b2")
        s.set_linewidth(1.0)
    ax.indicate_inset_zoom(axi, edgecolor="#9aa4b2", alpha=0.9, lw=1.0)

    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / "slide06_nested_rhythms"
    fig.savefig(f"{stem}.pdf", transparent=True, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {stem}.pdf and {stem}.png")


if __name__ == "__main__":
    main()
