#!/usr/bin/env python3
"""Substantiation figures for the partition-resolution-vs-fractional-depth claim.

One figure per load-bearing claim from the 2026-04-27 reconciliation report:

  fig1_tree_shape_inventory.pdf   — dmax is constant; dmin / h_rel_min wild
  fig2_gamma_h_rake_dendrograms   — Pat_05/10/14 γ_h rake topologies vs Pat_06
  fig3_hrel_vs_k_curves           — h_rel(k) per-patient overlay, 6 bands
  fig4_three_axis_dvi             — integer-k vs linspace-h_rel vs logspace-h_rel
                                     Δ_VI side-by-side
  fig5_split_baseline_strips      — d_VI_full vs drift_dVI per patient per ridge

All figures rasterised, PDF only, no fig.suptitle. Sidecar reading rules
in `figures/SUBSTANTIATION.md`.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE
from lrg_eegfc.utils.metrics.tree import dmax_from_Z
from lrg_eegfc.workflow.lrg import load_lrg_result


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES = ("rest_pre", "task_test", "rest_post")
RIDGES = {
    "delta":      [(20, 32)],
    "alpha":      [(20, 26)],
    "beta":       [(6, 7)],
    "high_gamma": [(21, 23), (30, 35)],
}
K_RANGE = list(range(2, 50))

OUT_DIR = ROOT / "data" / "audit" / "dvi_split_baseline" / "figures"


def _Z(pat, phase, band):
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return None if r is None else np.asarray(r.linkage_matrix)


# ─────────────────────────── Figure 1 — tree-shape inventory ───────────────


def fig1_tree_shape_inventory(inv: pd.DataFrame):
    fig, axes = plt.subplots(
        1, 3, figsize=(13.5, 4.0), dpi=160,
        gridspec_kw={"wspace": 0.32, "width_ratios": [1.0, 1.0, 1.6]},
    )

    # Panel A — dmax across all trees
    ax = axes[0]
    band_x = {b: i for i, b in enumerate(BRAIN_BANDS_NAMES)}
    rng = np.random.default_rng(0)
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        sub = inv[inv["band"] == band]
        x = ib + rng.uniform(-0.18, 0.18, size=len(sub))
        ax.scatter(x, sub["dmax"].to_numpy(), s=14, alpha=0.65,
                   color="0.2", linewidths=0)
    ax.set_xticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=9)
    ax.set_ylabel("dmax (root cophenetic height)", fontsize=10)
    ax.set_title("A — dmax: identical across every tree (= 1 / 1.01 ≈ 0.9901)",
                 fontsize=10, loc="left")
    ax.axhline(0.9901, color="red", lw=0.7, ls="--")
    ax.set_ylim(0.985, 0.995)
    ax.grid(axis="y", alpha=0.3)

    # Panel B — dmin per band (boxplot, all phases pooled)
    ax = axes[1]
    bands_data = [inv[inv["band"] == b]["dmin"].dropna().to_numpy()
                  for b in BRAIN_BANDS_NAMES]
    bp = ax.boxplot(bands_data, positions=range(len(BRAIN_BANDS_NAMES)),
                    widths=0.55, showfliers=True, patch_artist=True,
                    flierprops={"marker": "o", "markersize": 3,
                                "markerfacecolor": "red",
                                "markeredgecolor": "red", "alpha": 0.6})
    for patch, b in zip(bp["boxes"], BRAIN_BANDS_NAMES):
        patch.set_facecolor("0.85" if b != "high_gamma" else "#fde0a8")
        patch.set_edgecolor("0.2")
    for med in bp["medians"]:
        med.set_color("0.1")
    ax.set_xticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=9)
    ax.set_ylabel("dmin (smallest non-zero merge height)", fontsize=10)
    ax.set_title("B — dmin: wild variance, especially in $\\gamma_{\\mathrm{h}}$ (highlighted)",
                 fontsize=10, loc="left")
    ax.grid(axis="y", alpha=0.3)

    # Panel C — h_rel_min heatmap per (patient, band) pooling phases
    ax = axes[2]
    h_pivot = (inv.groupby(["patient", "band"])["h_rel_min"].max()
                  .unstack("band").reindex(index=COHORT_N10,
                                            columns=BRAIN_BANDS_NAMES))
    cmap = plt.get_cmap("YlOrRd")
    im = ax.imshow(h_pivot.values, aspect="auto", cmap=cmap,
                   vmin=0.0, vmax=0.75, interpolation="nearest",
                   rasterized=True)
    ax.set_xticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=9)
    ax.set_yticks(range(len(COHORT_N10)))
    ax.set_yticklabels(COHORT_N10, fontsize=8)
    for ip in range(len(COHORT_N10)):
        for ib in range(len(BRAIN_BANDS_NAMES)):
            v = h_pivot.values[ip, ib]
            if not np.isfinite(v):
                continue
            ax.text(ib, ip, f"{v:.2f}",
                    ha="center", va="center", fontsize=7,
                    color="white" if v > 0.45 else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("max $h_{rel}^{min}$ across phases", fontsize=9)
    ax.set_title("C — h_rel_min per (patient, band): $\\gamma_{\\mathrm{h}}$ Pat_05/10/14 are the rakes",
                 fontsize=10, loc="left")

    out = OUT_DIR / "fig1_tree_shape_inventory.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# ─────────────────────────── Figure 2 — γ_h rake dendrograms ───────────────


def fig2_gamma_h_rake_dendrograms():
    """4 patients × 1 phase (rest_pre) γ_h dendrograms.

    Pat_06 = balanced control; Pat_05 / Pat_10 / Pat_14 = rakes (sorted by
    h_rel_min severity). Same x-extent (per-tree), same y-range (log).
    """
    pats = ["Pat_06", "Pat_05", "Pat_10", "Pat_14"]
    fig, axes = plt.subplots(1, 4, figsize=(15.0, 4.0), dpi=160,
                              sharey=True,
                              gridspec_kw={"wspace": 0.10})
    for ax, pat in zip(axes, pats):
        Z = _Z(pat, "rest_pre", "high_gamma")
        if Z is None:
            ax.text(0.5, 0.5, "no Z", ha="center", va="center")
            continue
        d_min = float(Z[0, 2])
        d_max = float(Z[-1, 2])
        h_rel_min = d_min / d_max
        ddata = dendrogram(Z, no_plot=True, no_labels=True)
        # plot manually with set rasterization
        for xs, ys in zip(ddata["icoord"], ddata["dcoord"]):
            ax.plot(xs, ys, color="0.2", lw=0.5, rasterized=True)
        ax.set_yscale("log")
        merge_heights = np.sort(Z[:, 2])
        tmin = merge_heights[merge_heights > 0].min() * 0.8
        tmax = merge_heights[-1] * 1.05
        ax.set_ylim(tmin, tmax)
        ax.set_xticks([])
        # annotate dmin / h_rel_min
        ax.axhline(d_min, color="red", lw=0.8, ls="--", alpha=0.8)
        ax.axhline(d_max, color="0.4", lw=0.6, ls=":")
        title = (f"{pat} $\\gamma_{{\\mathrm{{h}}}}$ rest_pre\n"
                 f"$h_{{\\mathrm{{rel}}}}^{{\\min}}$ = {h_rel_min:.3f}, "
                 f"dmin = {d_min:.3f}")
        ax.set_title(title, fontsize=10, loc="left")
    axes[0].set_ylabel("merge height (log)", fontsize=10)
    out = OUT_DIR / "fig2_gamma_h_rake_dendrograms.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# ─────────────────────────── Figure 3 — h_rel(k) curves ───────────────────


def fig3_hrel_vs_k_curves():
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0), dpi=160,
                              sharex=True, sharey=True,
                              gridspec_kw={"hspace": 0.30, "wspace": 0.10})
    cmap = plt.get_cmap("tab10")
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes.flat[ib]
        for ip, pat in enumerate(COHORT_N10):
            Z = _Z(pat, "rest_pre", band)
            if Z is None:
                continue
            n_leaves = Z.shape[0] + 1
            d_max = dmax_from_Z(Z)
            h_rels = []
            for k in K_RANGE:
                if k >= n_leaves:
                    h_rels.append(0.0)
                else:
                    h_rels.append(Z[-k, 2] / d_max)
            ax.plot(K_RANGE, h_rels, color=cmap(ip % 10), lw=1.0,
                    alpha=0.85, label=pat)
        # ridge highlights
        for (k_lo, k_hi) in RIDGES.get(band, []):
            ax.axvspan(k_lo - 0.5, k_hi + 0.5, color="orange", alpha=0.18,
                       zorder=0)
        ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]}", fontsize=11, loc="left")
        ax.set_ylim(0.0, 1.0)
        ax.grid(alpha=0.25)
        if ib == 0:
            ax.legend(fontsize=7, ncol=2, loc="lower left", frameon=False)
    for ax in axes[1]:
        ax.set_xlabel("k (integer-k cut)", fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel("$h_{\\mathrm{rel}}$ at cut", fontsize=10)
    out = OUT_DIR / "fig3_hrel_vs_k_curves.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# ─────────────────────────── Figure 4 — three-axis Δ_VI ─────────────────────


def _cmap_div():
    return LinearSegmentedColormap.from_list(
        "trace_div",
        [(0.05, 0.20, 0.55), (1.0, 1.0, 1.0), (0.70, 0.05, 0.05)],
        N=256,
    )


def _draw_grid_panel(ax, df, contrast, axis_values, axis_label, log_x):
    B = len(BRAIN_BANDS_NAMES)
    n = len(axis_values)
    mean = np.full((B, n), np.nan)
    npos = np.zeros((B, n), dtype=int)
    is_int = isinstance(axis_values[0], (int, np.integer))
    col = "k" if is_int else "h_rel"
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        for jx, val in enumerate(axis_values):
            if is_int:
                sub = df[(df["band"] == band) & (df["k"] == int(val))]
            else:
                sub = df[(df["band"] == band)
                         & np.isclose(df["h_rel"], val, atol=1e-9)]
            x = sub[contrast].to_numpy()
            x = x[np.isfinite(x)]
            if x.size:
                mean[ib, jx] = x.mean()
                npos[ib, jx] = int((x > 0).sum())

    cmap = _cmap_div()
    vmax = np.nanmax(np.abs(mean)) if np.isfinite(mean).any() else 1.0
    norm = plt.Normalize(vmin=-vmax, vmax=+vmax)

    for ib in range(B):
        for jx in range(n):
            v = mean[ib, jx]
            if not np.isfinite(v):
                continue
            if jx == 0:
                lo = axis_values[0] - (axis_values[1] - axis_values[0]) / 2
            else:
                lo = (axis_values[jx - 1] + axis_values[jx]) / 2
            if jx == n - 1:
                hi = axis_values[-1] + (axis_values[-1] - axis_values[-2]) / 2
            else:
                hi = (axis_values[jx] + axis_values[jx + 1]) / 2
            ax.add_patch(Rectangle((lo, ib - 0.5), hi - lo, 1,
                                    facecolor=cmap(norm(v)),
                                    edgecolor="none", zorder=1))
            if npos[ib, jx] >= 8:
                ax.add_patch(Rectangle((lo, ib - 0.5), hi - lo, 1,
                                        facecolor="none", edgecolor="black",
                                        lw=1.4, zorder=3))
            elif npos[ib, jx] >= 6:
                ax.add_patch(Rectangle((lo, ib - 0.5), hi - lo, 1,
                                        facecolor="none", edgecolor="0.4",
                                        lw=0.8, ls=":", zorder=3))
    if log_x:
        ax.set_xscale("log")
    if is_int:
        ax.set_xlim(axis_values[0] - 0.5, axis_values[-1] + 0.5)
    else:
        ax.set_xlim(axis_values[0] * (0.95 if log_x else 1.0),
                    axis_values[-1] + (axis_values[-1] - axis_values[-2]) / 2)
    ax.set_ylim(B - 0.5, -0.5)
    ax.set_yticks(range(B))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=10)
    ax.set_xlabel(axis_label, fontsize=10)
    return cmap, norm, vmax


def fig4_three_axis_dvi():
    df_int = pd.read_csv(ROOT / "data" / "reports" / "imcoh_vi"
                          / "h2_partition_multiscale_raw.csv")
    df_h = pd.read_csv(ROOT / "data" / "audit" / "dvi_split_baseline"
                       / "dvi_hrel_n10_imcoh_abs.csv")
    df_lin = df_h[df_h["grid"] == "linspace"]
    df_log = df_h[df_h["grid"] == "logspace"]
    grid_lin = np.sort(df_lin["h_rel"].unique())
    grid_log = np.sort(df_log["h_rel"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.5), dpi=160,
                              gridspec_kw={"wspace": 0.20})

    cm, nm, vmax_int = _draw_grid_panel(
        axes[0], df_int, "d_VI", K_RANGE, "k (integer cut)", log_x=False,
    )
    axes[0].set_title("integer k — headline (L=13 ridge in $\\delta$ k=20–32, 10/10 at k=28)",
                      fontsize=10, loc="left")

    _draw_grid_panel(
        axes[1], df_lin, "d_VI", grid_lin,
        r"$h_{\mathrm{rel}}$ (linspace)", log_x=False,
    )
    axes[1].set_title("linspace $h_{\\mathrm{rel}}$ — same operator, different axis",
                      fontsize=10, loc="left")

    _draw_grid_panel(
        axes[2], df_log, "d_VI", grid_log,
        r"$h_{\mathrm{rel}}$ (logspace)", log_x=True,
    )
    axes[2].set_title("logspace $h_{\\mathrm{rel}}$",
                      fontsize=10, loc="left")

    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=nm, cmap=cm),
                        ax=axes.tolist(), shrink=0.85, pad=0.02,
                        fraction=0.025)
    cbar.set_label("$\\Delta_{\\mathrm{VI}}$ patient-mean", fontsize=10)
    out = OUT_DIR / "fig4_three_axis_dvi.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# ─────────────────────────── Figure 5 — split-baseline strips ───────────────


def fig5_split_baseline_strips():
    df = pd.read_csv(ROOT / "data" / "audit" / "dvi_split_baseline"
                     / "dvi_split_baseline_n10_imcoh_abs.csv")
    ridges = []
    for band, lst in RIDGES.items():
        for (k_lo, k_hi) in lst:
            ridges.append((band, k_lo, k_hi))
    n_r = len(ridges)
    fig, axes = plt.subplots(1, n_r, figsize=(3.4 * n_r, 4.0), dpi=160,
                              sharey=False,
                              gridspec_kw={"wspace": 0.22})
    if n_r == 1:
        axes = [axes]
    contrast_colors = {
        "d_VI_full": "#7f1d1d",
        "d_VI_A":    "#9c6644",
        "d_VI_B":    "#e09f3e",
        "drift_dVI": "#1f77b4",
    }
    for ax, (band, k_lo, k_hi) in zip(axes, ridges):
        sub = df[(df["band"] == band)
                  & (df["k"] >= k_lo) & (df["k"] <= k_hi)].copy()
        pp = sub.groupby("patient")[
            ["d_VI_full", "d_VI_A", "d_VI_B", "drift_dVI"]].mean()
        pp = pp.reindex(COHORT_N10)
        x = np.arange(len(COHORT_N10))
        w = 0.18
        for offset, col in zip([-1.5 * w, -0.5 * w, 0.5 * w, 1.5 * w],
                               ["d_VI_full", "d_VI_A", "d_VI_B", "drift_dVI"]):
            ax.bar(x + offset, pp[col].values, width=w, color=contrast_colors[col],
                   edgecolor="0.1", linewidth=0.4, label=col, alpha=0.9)
        ax.axhline(0, color="0.1", lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(COHORT_N10, rotation=70, fontsize=7)
        ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} k={k_lo}–{k_hi}",
                     fontsize=10, loc="left")
        ax.grid(axis="y", alpha=0.3)
        ax.set_ylabel(r"$\Delta_{\mathrm{VI}}$" if ax is axes[0] else "",
                      fontsize=9)
    axes[0].legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
    out = OUT_DIR / "fig5_split_baseline_strips.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# ─────────────────────────── main ──────────────────────────────────────────


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    inv = pd.read_csv(ROOT / "data" / "audit" / "dvi_split_baseline"
                      / "dmin_dmax_inventory.csv")

    print("Generating substantiation figures...")
    fig1_tree_shape_inventory(inv)
    fig2_gamma_h_rake_dendrograms()
    fig3_hrel_vs_k_curves()
    fig4_three_axis_dvi()
    fig5_split_baseline_strips()

    md = OUT_DIR / "SUBSTANTIATION.md"
    md.write_text("""---
name: substantiation-figures
type: figure-index
era: COHORT_N10
created: 2026-04-27
---

# Substantiation figures — one figure per claim from the 2026-04-27 reconciliation

| # | Claim being substantiated | Figure |
|---|---|---|
| 1 | dmax is constant by LRG normalization (= 1/1.01 ≈ 0.9901); dmin is wild | `fig1_tree_shape_inventory.pdf` |
| 2 | Pat_05 / Pat_10 / Pat_14 have rake-like $\\gamma_{\\mathrm{h}}$ trees with the first merge near the root | `fig2_gamma_h_rake_dendrograms.pdf` |
| 3 | Per-patient `h_rel(k)` curves diverge cohort-wide; the same integer k corresponds to very different fractional depths | `fig3_hrel_vs_k_curves.pdf` |
| 4 | The integer-k δ ridge does NOT translate to fixed-h_rel space — visual side-by-side | `fig4_three_axis_dvi.pdf` |
| 5 | Δ_VI split-baseline result per ridge (d_VI_full vs d_VI_A vs d_VI_B vs drift_dVI per patient) | `fig5_split_baseline_strips.pdf` |
""", encoding="utf-8")
    print(f"wrote {md}")


if __name__ == "__main__":
    main()
