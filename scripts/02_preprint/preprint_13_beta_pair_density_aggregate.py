#!/usr/bin/env python3
r"""β cohort-aggregate per-pair density scatter (fig5 aesthetic).

Single-panel figure replicating the
``fig5_distance_class_examples.pdf`` style for the LRG ρ_split^coph
per-pair layer.

Each point is one upper-triangular pair ``(i, j)`` from one of the
ten cohort patients; coordinates are the **within-patient rank** of
the cophenet-distance differences

.. math::

    \Delta_{\mathrm{task}}(i,j) &= D^{\mathrm{coph}}_{\rm test}(i,j)
                                 - D^{\mathrm{coph}}_{\rm preA}(i,j) \\
    \Delta_{\mathrm{rest}}(i,j) &= D^{\mathrm{coph}}_{\rm post}(i,j)
                                 - D^{\mathrm{coph}}_{\rm preB}(i,j)

per the audit_63 split-baseline recipe (§5.3 of the preprint).  The
OLS regression of ``v`` on ``u`` has slope identically the pooled
Spearman ``ρ_split^coph`` correlation; the cohort-median (Wilcoxon)
appears as a stats box matching the fig5 layout.

The visual contract: ~6.4k pairs × 10 patients ≈ 60–70k points fill
a 60×60 [0,1]² bin grid (cividis log(count+1)); the slight diagonal
tilt of the cloud is the cohort ρ ≈ +0.2 that §5.3 reports.  The
matched-strength noise floor ``ρ_surr^{p95}`` is overlaid as a
second, fainter line so the +0.2 reads above the null at the eye.

Output (PDF only)
-----------------
``data/preprint/figures/beta/per_pair_trace/fig_beta_pair_density_aggregate.pdf``
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import (FixedLocator, NullFormatter)
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

LRG_CTM_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
PAIR_SPLIT_DIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"

CMAP_DENSITY = "cividis"
N_BINS = 60
IDENTITY_COLOR = "white"
FIT_COLOR = "#ff5252"
NULL_FIT_COLOR = "#ffd166"   # warm amber — null-floor reference


# ---------------------------------------------------------------------------
def load_pair_data(pat: str, band: str) -> tuple[np.ndarray, np.ndarray]:
    npz_path = PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    return (np.asarray(d["dD_task"], dtype=float),
            np.asarray(d["dD_rest"], dtype=float))


def pool_ranks(band: str) -> tuple[np.ndarray, np.ndarray,
                                    np.ndarray, np.ndarray, list[int]]:
    """Pool within-patient rank pairs across the cohort.

    Returns
    -------
    u, v
        Concatenated within-patient ranks in [0, 1], same length.
    dt_raw, dr_raw
        Concatenated raw signed Δ's (used for the auxiliary stats box).
    sizes
        Per-patient pair count (for caption / annotation).
    """
    us, vs, dts, drs, sizes = [], [], [], [], []
    for pat in COHORT:
        dt, dr = load_pair_data(pat, band)
        n = len(dt)
        if n < 2:
            sizes.append(0)
            continue
        u = (rankdata(dt, method="ordinal") - 0.5) / n
        v = (rankdata(dr, method="ordinal") - 0.5) / n
        us.append(u); vs.append(v)
        dts.append(dt); drs.append(dr)
        sizes.append(n)
    return (np.concatenate(us), np.concatenate(vs),
            np.concatenate(dts), np.concatenate(drs), sizes)


def per_patient_rhos(band: str) -> np.ndarray:
    rhos = []
    for pat in COHORT:
        dt, dr = load_pair_data(pat, band)
        if len(dt) < 2:
            continue
        r, _ = spearmanr(dt, dr)
        rhos.append(float(r))
    return np.asarray(rhos)


# ---------------------------------------------------------------------------
def _apply_unit_ticks(axis):
    """Major ticks at 0 / ½ / 1 (matching the rank-rank scale)."""
    axis.set_major_locator(FixedLocator([0.0, 0.5, 1.0]))
    axis.set_minor_locator(FixedLocator(np.linspace(0.05, 0.95, 19)))
    axis.set_minor_formatter(NullFormatter())


def main(band: str = "beta") -> Path:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")

    cohort = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")
    cohort_row = cohort[cohort.band == band].iloc[0]
    rho_obs = float(cohort_row["obs_median_rho"])
    rho_surr_p95 = float(cohort_row["surr_median_rho_p95"])
    p_wilcoxon = float(cohort_row["paired_wilcoxon_p"])
    n_above_str = str(cohort_row["n_above_surrogate"])

    u, v, dt_raw, dr_raw, sizes = pool_ranks(band)
    n_total = u.size
    n_pat = int(np.sum(np.asarray(sizes) > 0))
    rhos = per_patient_rhos(band)

    # OLS slope of v on u; for centred ranks this equals the cohort
    # Spearman correlation of the pooled (u, v).
    pooled_rho, _ = spearmanr(u, v)

    # ------------------------------------------------------------------
    # Figure layout — single square panel.  Inset legend mirrors fig5.
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(7.4, 7.0))
    gs = fig.add_gridspec(
        1, 1, left=0.11, right=0.985, top=0.93, bottom=0.10,
    )
    ax = fig.add_subplot(gs[0, 0])
    cmap = plt.get_cmap(CMAP_DENSITY)
    ax.set_facecolor(cmap(0.0))

    # 2D histogram on [0,1]², log(count+1) shading — same recipe as
    # fig5 (linear bins here since rank marginals are uniform).
    bins = np.linspace(0.0, 1.0, N_BINS + 1)
    H, xedges, yedges = np.histogram2d(u, v, bins=[bins, bins])
    H_log = np.log10(H + 1.0)
    pcm = ax.pcolormesh(xedges, yedges, H_log.T, cmap=cmap,
                        shading="auto", zorder=1)

    # Identity diagonal (ρ=1) and anti-diagonal (ρ=-1) references.
    ax.plot([0, 1], [0, 1], color=IDENTITY_COLOR, lw=1.4, ls="--",
            alpha=0.85, zorder=3, label="identity ($\\rho=+1$)")
    ax.plot([0, 1], [1, 0], color=IDENTITY_COLOR, lw=0.9, ls=":",
            alpha=0.55, zorder=3, label="anti ($\\rho=-1$)")

    # OLS fit of v on u — slope equals pooled Spearman correlation.
    # We anchor the line at the centre (½, ½) since centred ranks have
    # zero mean by construction.
    slope = float(pooled_rho)
    x_line = np.linspace(0.0, 1.0, 100)
    y_line = 0.5 + slope * (x_line - 0.5)
    ax.plot(x_line, y_line, color=FIT_COLOR, lw=2.0, zorder=5,
            label=rf"cohort fit (slope = $\rho^{{coph}}_{{split}}$"
                  rf" = {pooled_rho:+.3f})")

    # Matched-strength noise-floor reference — what a ρ_surr^p95
    # cohort would look like under this lens.  Drawn fainter so the
    # observed +0.2 reads as visually above noise.
    y_null = 0.5 + rho_surr_p95 * (x_line - 0.5)
    ax.plot(x_line, y_null, color=NULL_FIT_COLOR, lw=1.4, ls="-",
            alpha=0.92, zorder=4,
            label=rf"matched-strength $p_{{95}}$ "
                  rf"(slope = {rho_surr_p95:+.3f})")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    _apply_unit_ticks(ax.xaxis); _apply_unit_ticks(ax.yaxis)
    ax.set_xticklabels(["0", r"$\frac{1}{2}$", "1"])
    ax.set_yticklabels(["0", r"$\frac{1}{2}$", "1"])
    ax.set_xlabel(r"within-patient rank of "
                  r"$\Delta_{\mathrm{task}}(i,j)$  ($u$)",
                  fontsize=10.5, labelpad=4)
    ax.set_ylabel(r"within-patient rank of "
                  r"$\Delta_{\mathrm{rest}}(i,j)$  ($v$)",
                  fontsize=10.5, labelpad=4)
    ax.tick_params(labelsize=9.5, pad=2)

    # Colorbar at the right of the main axes — log(count+1) scale.
    cbar_ax = fig.add_axes([1.0, 0.10, 0.022, 0.83])
    cb = fig.colorbar(pcm, cax=cbar_ax, orientation="vertical")
    cb.set_label(r"$\log_{10}(\mathrm{count} + 1)$  per bin",
                 rotation=270, labelpad=14, fontsize=9.5)
    cb.ax.tick_params(labelsize=8.5)

    # Legend (fig5-style: upper-left, slim).
    leg = ax.legend(loc="upper left", fontsize=8.0, frameon=True,
                    facecolor="white", edgecolor="#bbb",
                    framealpha=0.92, borderpad=0.35,
                    handlelength=1.8, handletextpad=0.45)
    leg.set_zorder(8)

    # Stats annotation box — bottom-right, fig5 style.
    rho_mean = float(np.mean(rhos))
    rho_med = float(np.median(rhos))
    annot = (
        rf"cohort  $n_{{\rm pat}} = {n_pat}$" "\n"
        rf"$n_{{\rm pairs}} = {n_total:,}$" "\n"
        rf"$\rho^{{coph}}_{{split}}$ pooled $= {pooled_rho:+.3f}$" "\n"
        rf"per-pat. median $= {rho_med:+.3f}$" "\n"
        rf"Wilcoxon $p = {p_wilcoxon:.4f}$" "\n"
        rf"$n_{{\rm above\,surr.}} = {n_above_str}$"
    )
    ax.text(0.965, 0.04, annot,
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.6, fontweight="normal", color="0.10",
            bbox=dict(boxstyle="round,pad=0.32",
                      facecolor="white", edgecolor="#999",
                      alpha=0.92), zorder=8)

    fig.text(0.04, 0.965, rf"{band_tex} pair density   "
             rf"$\Delta_{{\rm task}}$ vs $\Delta_{{\rm rest}}$ "
             rf"(within-patient ranks, pooled $n=10$)",
             ha="left", va="top", fontsize=11.5, fontweight="bold")

    out_dir = ROOT / "data" / "preprint" / "figures" / band / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"fig_{band}_pair_density_aggregate.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  pooled n_pairs = {n_total:,}  ({n_pat} patients)")
    print(f"  cohort ρ_obs = {rho_obs:+.4f}  pooled ρ = {pooled_rho:+.4f}")
    print(f"  surr_p95     = {rho_surr_p95:+.4f}    Wilcoxon p = {p_wilcoxon:.4g}")
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", default="beta",
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    args = parser.parse_args()
    main(args.band)
