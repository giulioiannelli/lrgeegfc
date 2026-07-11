#!/usr/bin/env python3
r"""fig:epi (§3 compound) — one operator, a second read-out: the propagator localizes the SOZ.

Assembles the §3 panels by importing the sibling single-panel scripts and delegating to their
``draw`` helpers (single source of truth -- no forked panel code). The only raster is the
Plotly 3-D seed-spread brain, brought in as its rendered PNG via imshow (exactly as fig_trace1
does for its 3-D panel).

Two page-fitting floats (five full-width analytical strips cannot share one portrait page
legibly -- the same reason §1/§2 each use two floats):

  fig_epi1  — the finding + the detector (R3.1--R3.2)
    (a) seed-spread brain   heat on known seizure seeds diffuses to the unseeded, off-shaft
                            seizure contacts (mechanism; Pat_08 delta).
    (b) relational marker   strength-residual heat-kernel affinity AUC per band.
    (c) calibrated detector 6-band fusion lifts precision@5 34->60%; calibrated P(SOZ).

  fig_epi2  — the two populations + the honest ceiling (R3.3--R3.4)
    (a) two populations     per-patient AUC: community vs hub, and the regime switch.
    (b) scope & ceiling     precision@k vs base rate, gray/white tissue, occult sites.

Run the single-panel scripts first (they write the standalone PDFs and, for the brain, the PNG):
  fig_epi_seed_spread_3d.py  fig_epi_a_relational_marker.py  fig_epi_b_calibrated_detector.py
  fig_epi_c_two_populations.py  fig_epi_d_scope_ceiling.py

Reads : data/preprint/figures/results_section3/fig_epi_seed_spread_3d.png (+ the panel CSVs)
Writes: data/preprint/figures/results_section3/{fig_epi1,fig_epi2}.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

# sibling single-panel scripts (figures_embedded is not a package) importable by name
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig_epi_a_relational_marker as p_a      # noqa: E402
import fig_epi_b_calibrated_detector as p_b     # noqa: E402
import fig_epi_c_two_populations as p_c         # noqa: E402
import fig_epi_d_scope_ceiling as p_d           # noqa: E402

OUTDIR = ROOT / "data/preprint/figures/results_section3"
BRAIN_PNG = OUTDIR / "fig_epi_seed_spread_3d.png"
OUT1 = OUTDIR / "fig_epi1.pdf"
OUT2 = OUTDIR / "fig_epi2.pdf"


def _crop_to_content(img, pad_frac=0.01):
    """Trim the transparent margin of an RGBA array to its content bbox (+ small pad)."""
    if img.shape[-1] < 4:
        return img
    alpha = img[..., 3] > 0.02
    rows = np.where(alpha.any(1))[0]
    cols = np.where(alpha.any(0))[0]
    if not len(rows) or not len(cols):
        return img
    r0, r1, c0, c1 = rows[0], rows[-1], cols[0], cols[-1]
    ph, pw = int((r1 - r0) * pad_frac), int((c1 - c0) * pad_frac)
    r0, r1 = max(0, r0 - ph), min(img.shape[0], r1 + ph + 1)
    c0, c1 = max(0, c0 - pw), min(img.shape[1], c1 + pw + 1)
    return img[r0:r1, c0:c1]


def _letter(subfig, s, x=0.006, y=0.985):
    subfig.text(x, y, s, fontsize=21, fontweight="bold", va="top", ha="left")


def build_fig1():
    """The finding + detector: brain (a), relational band AUC (b), calibrated detector (c)."""
    if not BRAIN_PNG.exists():
        raise SystemExit(f"missing {BRAIN_PNG} — run fig_epi_seed_spread_3d.py first")

    fig = plt.figure(figsize=(12.6, 15.0))
    sfs = fig.subfigures(3, 1, height_ratios=[5.5, 4.9, 4.6], hspace=0.02)

    brain = _crop_to_content(mpimg.imread(str(BRAIN_PNG)))
    axa = sfs[0].add_axes([0.02, 0.0, 0.96, 1.0])
    axa.imshow(brain, interpolation="lanczos")
    axa.axis("off")
    _letter(sfs[0], r"$\mathbf{a}$")

    p_a.draw(sfs[1]);                 _letter(sfs[1], r"$\mathbf{b}$")
    p_b.draw(sfs[2], letters=False);  _letter(sfs[2], r"$\mathbf{c}$")

    OUT1.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT1, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {OUT1}")


def build_fig2():
    """The scope: two populations (a), precision/ceiling/tissue/occult (b)."""
    fig = plt.figure(figsize=(12.6, 10.4))
    sfs = fig.subfigures(2, 1, height_ratios=[5.2, 4.9], hspace=0.03)

    p_c.draw(sfs[0], letters=False);  _letter(sfs[0], r"$\mathbf{a}$")
    p_d.draw(sfs[1], letters=False);  _letter(sfs[1], r"$\mathbf{b}$")

    OUT2.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT2, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {OUT2}")


def main():
    build_fig1()
    build_fig2()


if __name__ == "__main__":
    main()
