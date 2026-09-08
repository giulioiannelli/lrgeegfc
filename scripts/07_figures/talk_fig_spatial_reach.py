#!/usr/bin/env python3
r"""talk_fig_spatial_reach -- slide 14b inset: the physical scale the diffusion resolves,
from the contact pitch to the implant span, as the diffusion scale s = tau*lambda_max grows.

Metric (STARTS AT THE PITCH, unlike the max-diameter reach): at scale s we threshold the
heat kernel K = e^{-tau L} on the mst@0.20 backbone at the uniform-equilibrium share 1/N,
take the connected components (functional groups), and report

    ell(s) = mean over contacts of  { extent of the contact's group } ,
             extent = max pairwise physical distance in the group  (a singleton = the pitch)

So at fine s every contact is its own unit -> ell = pitch (~3.5 mm); as s grows the groups
merge and physically enlarge, up to the whole implant -> ell = span (~9 cm). This mirrors
the micro -> meso -> macro coarse-graining panel and, unlike the max-diameter reach, it does
NOT floor at ~10 mm: the finest resolved scale is the single contact, exactly as it should be.

Cohort median +/- IQR per band, rest_post, mst@0.20. Descriptive geometry -- no null.

Writes: data/outputs/figures/talk/fig_spatial_reach.png   (transparent, Canva-ready)
"""
from __future__ import annotations
import sys

import numpy as np
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import connected_components

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS, load_phase          # noqa: E402
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction        # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig          # noqa: E402
from lrg_eegfc.utils.io.regions import load_channel_regions           # noqa: E402
from lrg_eegfc.config.paths import FIGURES_ROOT                       # noqa: E402
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT                # noqa: E402
from lrg_eegfc.visuals.styles import use_lrg_style, band_color        # noqa: E402

use_lrg_style()

FRAC = 0.20
PHASE = "rest_post"
SGRID = np.logspace(np.log10(0.04), np.log10(180.0), 28)   # s = tau*lambda_max (computed)
XLIM = (0.15, 3.0)          # display window: trim the long pitch/span plateaus, keep the rise
OUT = FIGURES_ROOT / "talk" / "fig_spatial_reach.png"
FIG_W, FIG_H = 7.6, 5.1     # inches (for square-pixel brain insets)
# 3D brain panels hung in the white areas: (tag, centre x, centre y, height) in FIGURE coords.
# micro low on the bottom-left; meso higher and centre-right (before the curve's rise);
# macro in the big right-side white area.
BRAIN_INSETS = [("micro", 0.215, 0.34, 0.32),
                ("meso",  0.40,  0.66, 0.32),
                ("macro", 0.82,  0.52, 0.34)]


def _crop_alpha(img, pad=3):
    a = img[..., 3] > 0.02
    ys, xs = np.where(a)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + 1 + pad, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + 1 + pad, img.shape[1])
    return img[y0:y1, x0:x1]


def add_brain_inset(fig, tag, cx, cy, h):
    """Place cropped brain panel `_scale_brain_{tag}.png` centred at (cx,cy) fig-coords."""
    img = _crop_alpha(plt.imread(FIGURES_ROOT / "talk" / f"_scale_brain_{tag}.png"))
    ih, iw = img.shape[:2]
    w = h * (iw / ih) * (FIG_H / FIG_W)                      # keep pixels square
    axb = fig.add_axes([cx - w / 2, cy - h / 2, w, h], zorder=5)
    axb.imshow(img)
    axb.axis("off")
    axb.patch.set_alpha(0.0)


def phys(pat):
    rdf = load_channel_regions(pat)
    xyz = rdf[["x", "y", "z"]].to_numpy(float) / 1000.0      # um -> mm
    mask = np.isfinite(xyz).all(axis=1)
    return xyz, mask


CACHE = FIGURES_ROOT / "talk" / "_reach_cohort_cache.npz"


def compute_cohort():
    """{band: (n_pat, n_s) ell array}, cohort pitch/span. Cached to CACHE."""
    if CACHE.exists():
        z = np.load(CACHE)
        acc = {b: z[b] for b in BANDS}
        return acc, float(z["pitch_c"]), float(z["span_c"])
    ncols = len(SGRID)
    acc = {b: [np.full(ncols, np.nan) for _ in COHORT] for b in BANDS}
    pitches, spans = [], []
    for ip, pat in enumerate(COHORT):
        xyz, mask = phys(pat)
        p = xyz[mask]
        D = cdist(p, p)
        Dnn = D.copy(); np.fill_diagonal(Dnn, np.inf)
        pitch = float(np.median(Dnn.min(1)))
        span = float(D.max())
        pitches.append(pitch); spans.append(span)
        midx = np.where(mask)[0]
        for band in BANDS:
            W = load_phase(pat, PHASE, band)
            if W.shape[0] != xyz.shape[0]:
                print(f"  [skip] {pat} {band}: N mismatch", flush=True)
                continue
            ev, V = laplacian_eig(mst_union_top_fraction(W, FRAC))
            thr = 1.0 / V.shape[0]
            ell = np.full(ncols, np.nan)
            for k, s in enumerate(SGRID):
                K = (V * np.exp(-(s / ev[-1]) * ev)) @ V.T
                comm = K[np.ix_(midx, midx)] >= thr
                np.fill_diagonal(comm, False)
                ncomp, lab = connected_components(comm, directed=False)
                ext = np.empty(len(midx))
                for c in range(ncomp):
                    g = np.where(lab == c)[0]
                    ext[g] = pitch if len(g) == 1 else D[np.ix_(g, g)].max()
                ell[k] = ext.mean()
            acc[band][ip] = ell
        print(f"[{ip+1}/{len(COHORT)}] {pat} pitch={pitch:.1f} span={span:.0f}", flush=True)
    acc = {b: np.vstack(acc[b]) for b in BANDS}
    pitch_c, span_c = float(np.median(pitches)), float(np.median(spans))
    np.savez(CACHE, pitch_c=pitch_c, span_c=span_c, **acc)
    return acc, pitch_c, span_c


def main():
    acc, pitch_c, span_c = compute_cohort()

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    for band in BANDS:                                       # clean median lines (spread muddy w/ 6 bands)
        med = np.nanmedian(acc[band], axis=0)
        ax.plot(SGRID, med, color=band_color(band), lw=2.8,
                label=BRAIN_BAND_TEX_DICT.get(band, band))

    ax.axhline(pitch_c, ls=":", lw=1.3, color="0.45")
    ax.axhline(span_c, ls="-.", lw=1.3, color="0.45")
    ax.axvline(1.0, ls="--", lw=1.2, color="0.6")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*XLIM)
    ax.set_xlabel(r"diffusion scale   $s=\tau\lambda_{\max}$")
    ax.set_ylabel(r"physical size of the functional group   $\ell(s)$  [mm]")
    ax.set_ylim(pitch_c * 0.8, span_c * 1.25)
    ax.legend(title="band", frameon=False, ncol=2, fontsize=10, loc="lower right")
    fig.tight_layout()
    for tag, cx, cy, h in BRAIN_INSETS:                      # 3D brains in the white corners
        add_brain_inset(fig, tag, cx, cy, h)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, dpi=200)
    plt.close(fig)

    print(f"\ncohort: pitch~{pitch_c:.1f}mm  span~{span_c:.0f}mm")
    for band in BANDS:
        med = np.nanmedian(acc[band], axis=0)
        print(f"  {band:>10s}: ell(finest)~{med[0]:.1f}mm  ell(coarse)~{med[-1]:.0f}mm")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
