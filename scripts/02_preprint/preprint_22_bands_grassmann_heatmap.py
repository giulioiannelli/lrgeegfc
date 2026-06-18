#!/usr/bin/env python3
"""All-bands Grassmann principal-angle residual heatmap — 2×3 compound.

The per-band figures (`preprint_09_beta_grassmann_heatmap.py`) render one
full-page composite per band: principal-angle residual heatmap + per-mode
profiles + T_G(k) line + (now-retracted) glass-brain inset. This script
**strips out only the main heatmap panel** and tiles all six bands in a
2×3 grid (top row δ θ α; bottom row β γ_l γ_h), the global-multiscale
counterpart of the per-pair joint-density compound
(`preprint_18_bands_joint_density_rawfc.py`). One shared colorbar.

The plotted quantity per panel is the matched-strength residual

    H[i, k] = Δθ_i(k)_obs_cohort_median − Δθ_i(k)_surr_cohort_median   (rad)

weighted by (a) the per-k patient-bootstrap evidence
    w(k) = clip(−log10 p_boot(k) / −log10 0.05, 0, 1)
and (b) a band-level verdict scale. **Red on the i ≤ k triangle = trace**
(rsPost subspace closer to taskT than rsPre, in excess of matched-strength
expectation); blue = anti-trace; white = no excess / empty (i > k).

Verdict-scale source — THE LOCKED LEDGER, not the CSV
-----------------------------------------------------
The band-level mute factor encodes the locked Grassmann verdict tier
(`.agents/preprint/locked/VERDICT_LEDGER.md`, Decision 12, 2026-05-28):
strong → 1.00, weak → 0.50, no_trace → 0.15. We read the tier from the
ledger, NOT from `grassmann_cluster_extent/cohort_summary.csv`: that CSV's
`verdict_cluster_extent` column is stale (pre-Decision-12 — it still tags
δ "strong" off the mechanical `cluster_p_mass < 0.01` rule). Decision 12
demoted δ to **weak** because full-data LOO max p_mass = 0.055 (Pat_08)
fails the < 0.05 LOO precondition. The per-band `_visual_bootstrap` PDFs
inherit the CSV staleness and over-weight δ; this compound corrects it.

Locked Grassmann band verdicts (VERDICT_LEDGER.md):
    β    strong trace, both probes (LOO-robust)
    γ_l  strong trace, only Grassmann (LOO-robust)
    δ    weak  trace, only Grassmann (LOO binds at Pat_08; C5 epi-X resolves)
    γ_h  no trace (cluster p_mass = 0.060)
    θ    no trace
    α    no trace on Grassmann (strong only on D_coph)

Inputs (all cached — instant)
-----------------------------
data/preprint/cache/grassmann_principal_angles_obs_rest_pre_A_<band>_full_k.csv
data/preprint/cache/grassmann_principal_angles_surr_cohort_<band>_full_k.csv
data/preprint/cache/grassmann_bootstrap_p_<band>_full_B2000_seed20260528.csv
(loaded via the cached compute helpers in preprint_09_beta_grassmann_heatmap)

Output (PDF only, no PNG sibling)
---------------------------------
data/preprint/figures/all_bands/fig_bands_grassmann_heatmap.pdf
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()

# Reuse the cached compute from the per-band script (no duplication of the
# principal-angle / residual-heatmap / bootstrap machinery). Importing it runs
# that module's own `use_lrg_style()` + LABEL_SCALE rcParams block; we re-assert
# our own compact compound style immediately afterwards.
sys.path.insert(0, str(Path(__file__).resolve().parent))
g09 = importlib.import_module("preprint_09_beta_grassmann_heatmap")

use_lrg_style()

K_MIN, K_MAX, I_MAX = g09.K_MIN, g09.K_MAX, g09.I_MAX
BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

# Locked Grassmann verdict tier per band (VERDICT_LEDGER.md, Decision 12).
GRASS_VERDICT_TIER = {
    "beta": "strong",
    "low_gamma": "strong",
    "delta": "weak",
    "theta": "no_trace",
    "alpha": "no_trace",
    "high_gamma": "no_trace",
}
BAND_SCALE = {"strong": 1.00, "weak": 0.50, "no_trace": 0.15}

# Fixed residual colour scale — same as the per-band figure so panel
# intensities are directly comparable across the two figures.
VMAX_FIXED = 0.07

# Diverging trace colormap (deep blue = anti-trace, white = none,
# deep red = trace) — identical to the per-band heatmap.
TRACE_POS_COLOR = "#b8001f"
TRACE_NEG_COLOR = "#001fb8"
TRACE_CMAP = LinearSegmentedColormap.from_list(
    "trace_continuous",
    [(0.00, TRACE_NEG_COLOR), (0.50, "#ffffff"), (1.00, TRACE_POS_COLOR)],
)
TRACE_CMAP.set_bad("white")  # i > k triangle


def weighted_residual(band: str) -> tuple[np.ndarray, np.ndarray]:
    """H_weighted[i, k] and the k axis for *band* (bootstrap × verdict scale)."""
    pa = g09.compute_principal_angles_full_k(band)
    pa_surr = g09.compute_surrogate_principal_angles_cohort(band)
    H, ks = g09.build_heatmap(pa, pa_surr)

    # Per-k patient-bootstrap weight (self-protects against single-patient
    # leverage; same recipe as the per-band `_visual_bootstrap` figure).
    ks_full, _s_obs, p_full = g09.cohort_strip_bootstrap(band, source="full")
    log_05 = -np.log10(0.05)
    neglog_p = -np.log10(np.maximum(p_full, 1e-12))
    weight_per_k = np.clip(neglog_p / log_05, 0.0, 1.0)
    if len(weight_per_k) != H.shape[1]:
        weight_per_k = np.interp(ks, ks_full, weight_per_k)

    # Band-level verdict scale — from the LOCKED ledger (see module docstring).
    weight_per_k = BAND_SCALE[GRASS_VERDICT_TIER[band]] * weight_per_k
    H_weighted = H * np.tile(weight_per_k, (I_MAX, 1))
    return H_weighted, ks


def plot_band(ax, band: str, H_weighted: np.ndarray, ks: np.ndarray) -> object:
    H_masked = np.ma.masked_invalid(H_weighted)
    im = ax.imshow(
        H_masked, aspect="auto", cmap=TRACE_CMAP,
        vmin=-VMAX_FIXED, vmax=VMAX_FIXED,
        origin="upper", interpolation="nearest",
        extent=(ks.min() - 0.5, ks.max() + 0.5, I_MAX + 0.5, 0.5),
    )
    # i = k diagonal (boundary of the populated upper-right triangle).
    ax.plot([ks.min() - 0.5, ks.max() + 0.5],
            [ks.min() - 0.5, ks.max() + 0.5],
            color="0.40", lw=0.7, ls="--", zorder=4)
    ax.set_xlim(ks.min() - 0.5, ks.max() + 0.5)
    ax.set_ylim(I_MAX + 0.5, 0.5)
    ticks = [2, 20, 40, 60, 80, 100, 112]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_yticklabels([str(t) for t in ticks])
    ax.tick_params(labelsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_box_aspect(1.0)

    # Band symbol in the empty lower-left (large i, small k) triangle.
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    ax.text(0.06, 0.06, band_tex, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=24, fontweight="bold",
            color="0.18", zorder=6)
    return im


def main() -> Path:
    fields = {}
    for band in BAND_ORDER:
        H_w, ks = weighted_residual(band)
        fields[band] = (H_w, ks)
        tier = GRASS_VERDICT_TIER[band]
        finite = H_w[np.isfinite(H_w)]
        pos_frac = float((finite > 0).mean()) if finite.size else float("nan")
        print(f"  {band:11s} tier={tier:8s} scale={BAND_SCALE[tier]:.2f}  "
              f"max|H_w|={np.nanmax(np.abs(H_w)):.3f}  pos_frac={pos_frac:.2f}")

    fig = plt.figure(figsize=(15.0, 10.0))
    outer = fig.add_gridspec(2, 3, wspace=0.22, hspace=0.24,
                             left=0.06, right=0.91, top=0.95, bottom=0.08)
    last_im = None
    for idx, band in enumerate(BAND_ORDER):
        rr, cc = idx // 3, idx % 3
        ax = fig.add_subplot(outer[rr, cc])
        H_w, ks = fields[band]
        last_im = plot_band(ax, band, H_w, ks)
        # Shared axis labels only on the outer edges.
        if cc == 0:
            ax.set_ylabel(r"mode index $i$")
        if rr == 1:
            ax.set_xlabel(r"subspace cutoff $k$")

    cbar_ax = fig.add_axes([0.925, 0.18, 0.014, 0.60])
    cb = fig.colorbar(last_im, cax=cbar_ax, orientation="vertical",
                      extend="both")
    cb.set_label(r"$w(k)\,\Delta\theta_i(k)$  (rad)   "
                 r"(red = trace, blue = anti-trace)",
                 rotation=270, labelpad=24)
    cb.set_ticks([-VMAX_FIXED, -0.5 * VMAX_FIXED, 0.0,
                  0.5 * VMAX_FIXED, VMAX_FIXED])

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_grassmann_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
