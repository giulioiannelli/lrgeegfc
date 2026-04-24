#!/usr/bin/env python3
"""Section 2 Figure D: probe enrichment visualization (redesigned).

D1 — Paired dot plot with per-patient MSC→ImCoh transitions + kernel density.
     Each patient is a connected pair showing the bias reduction at a glance.
D2 — Bias-reduction heatmap + summary curves showing the scale-dependent
     improvement of ImCoh over MSC (ratio of enrichment ratios).

Run:
  python scripts/10_notes_imcoh/fig_D_enrichment.py [-v]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde

from lrg_eegfc.utils.metrics.hypothesis import (
    ALL_PATIENTS, BANDS, BAND_TEX, BRAIN_BAND_TEX_DICT,
    PATIENT_COLORS, N_COMMUNITIES,
    CLR_MSC, CLR_IMCOH,
    load_channel_labels, apply_pub_style, save_fig, SECTION2_ROOT,
    compute_probe_weight_ratio, extract_probe_labels,
    compute_enrichment_vs_scale,
    SECTION2_METHODS as METHODS,
    SECTION2_METHOD_LABELS as METHOD_LABELS
)
from lrg_eegfc.utils.scripting import load_all_fc, load_all_lrg



# ---------------------------------------------------------------------------
# Figure D1 — Paired dot plot with distribution
# ---------------------------------------------------------------------------

def fig_d1_paired_dots(
    patients: list[str],
    phase: str = "rest_pre",
    output_dir: Path = SECTION2_ROOT / "fig_D",
    verbose: bool = False,
):
    """Paired scatter with marginal distributions.

    Main square panel
      - x = MSC enrichment ratio (log scale)
      - y = ImCoh enrichment ratio (log scale)
      - diagonal y = x: points below it = ImCoh reduced the bias
      - dashed lines at x = 1 and y = 1 = unbiased reference
      - each dot is one (patient, band) pair, colored by patient

    Top marginal
      - distribution of MSC values (= x projection)
    Right marginal
      - distribution of ImCoh values (= y projection)

    The reader can see at a glance that the cloud is shifted below the
    diagonal and concentrated near y = 1 while x ≫ 1.
    """
    from matplotlib.lines import Line2D
    from matplotlib import gridspec

    rows = []
    for method in METHODS:
        fc_data = load_all_fc(
            patients=patients, bands=BANDS, phases=[phase],
            fc_method=method, verbose=False,
        )
        for (pat, ph, band), W in fc_data.items():
            ch = load_channel_labels(pat)
            ratio = compute_probe_weight_ratio(W, ch)
            rows.append({"patient": pat, "band": band, "method": method, "ratio": ratio})
    df = pd.DataFrame(rows)
    if df.empty:
        print("  no data")
        return df

    msc_map = {(r.patient, r.band): r.ratio for _, r in df[df.method == "msc"].iterrows()}
    imcoh_map = {(r.patient, r.band): r.ratio for _, r in df[df.method == "imcoh_abs"].iterrows()}
    paired = [(pat, band, msc_map[(pat, band)], imcoh_map[(pat, band)])
              for (pat, band) in msc_map if (pat, band) in imcoh_map]
    if not paired:
        return df

    msc_vals = np.array([p[2] for p in paired])
    imcoh_vals = np.array([p[3] for p in paired])

    # Shared log-range for both axes and diagonal
    lo = min(msc_vals.min(), imcoh_vals.min()) * 0.6
    hi = max(msc_vals.max(), imcoh_vals.max()) * 1.4

    # Figure with gridspec: scatter bottom-left, top/right marginals
    fig = plt.figure(figsize=(9, 9))
    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        width_ratios=[4, 1], height_ratios=[1, 4],
        wspace=0.04, hspace=0.04,
    )
    ax_scatter = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_scatter)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_scatter)

    # --- main scatter ---
    # Diagonal y = x
    diag = np.logspace(np.log10(lo), np.log10(hi), 50)
    ax_scatter.plot(diag, diag, color="black", lw=1.0, ls="--", alpha=0.7,
                    zorder=1, label="y = x  (no change)")
    # Unbiased reference lines
    ax_scatter.axhline(1.0, color=CLR_IMCOH, lw=1.0, ls=":", alpha=0.7)
    ax_scatter.axvline(1.0, color=CLR_MSC, lw=1.0, ls=":", alpha=0.7)
    # Unbiased band on ImCoh axis (horizontal)
    ax_scatter.axhspan(0.7, 1.3, color="#4CAF50", alpha=0.08, zorder=0)

    for pat, band, mv, iv in paired:
        ax_scatter.scatter(mv, iv, s=55, color=PATIENT_COLORS[pat],
                           edgecolors="white", linewidths=0.6,
                           alpha=0.9, zorder=5)

    ax_scatter.set_xscale("log")
    ax_scatter.set_yscale("log")
    ax_scatter.set_xlim(lo, hi)
    ax_scatter.set_ylim(lo, hi)
    ax_scatter.set_xlabel("MSC same-probe / cross-probe ratio",
                          fontsize=12, fontweight="bold", color=CLR_MSC)
    ax_scatter.set_ylabel(r"$|\mathrm{ImCoh}|$ same-probe / cross-probe ratio",
                          fontsize=12, fontweight="bold", color=CLR_IMCOH)
    ax_scatter.grid(which="both", alpha=0.15)
    ax_scatter.set_aspect("equal", adjustable="box")

    # --- top marginal: MSC distribution (histogram + KDE) ---
    log_m = np.log10(msc_vals)
    bins = np.logspace(np.log10(lo), np.log10(hi), 22)
    ax_top.hist(msc_vals, bins=bins, color=CLR_MSC, alpha=0.45,
                edgecolor="white", linewidth=0.5)
    if len(log_m) >= 3:
        kde = gaussian_kde(log_m, bw_method=0.3)
        xs = np.linspace(np.log10(lo), np.log10(hi), 200)
        dens = kde(xs)
        # scale to match hist count
        bin_width = (np.log10(hi) - np.log10(lo)) / 22
        ax_top.plot(10**xs, dens * len(log_m) * bin_width,
                    color=CLR_MSC, lw=2)
    med_m = np.median(msc_vals)
    ax_top.axvline(med_m, color=CLR_MSC, lw=2, ls="--")
    ax_top.axvline(1.0, color="gray", lw=1, ls=":")
    ax_top.text(med_m, ax_top.get_ylim()[1] * 0.9,
                f" median MSC\n {med_m:.2f}×",
                color=CLR_MSC, fontsize=10, fontweight="bold",
                ha="left", va="top")
    ax_top.set_xscale("log")
    ax_top.tick_params(labelbottom=False, left=False, labelleft=False)
    for spine in ("top", "right", "left"):
        ax_top.spines[spine].set_visible(False)
    ax_top.set_title(f"Probe-bias shift — {phase}  ({len(paired)} patient×band cells)",
                     fontsize=14, fontweight="bold", loc="left")

    # --- right marginal: ImCoh distribution ---
    log_i = np.log10(imcoh_vals)
    ax_right.hist(imcoh_vals, bins=bins, orientation="horizontal",
                  color=CLR_IMCOH, alpha=0.45,
                  edgecolor="white", linewidth=0.5)
    if len(log_i) >= 3:
        kde = gaussian_kde(log_i, bw_method=0.3)
        ys = np.linspace(np.log10(lo), np.log10(hi), 200)
        dens = kde(ys)
        bin_width = (np.log10(hi) - np.log10(lo)) / 22
        ax_right.plot(dens * len(log_i) * bin_width, 10**ys,
                      color=CLR_IMCOH, lw=2)
    med_i = np.median(imcoh_vals)
    ax_right.axhline(med_i, color=CLR_IMCOH, lw=2, ls="--")
    ax_right.axhline(1.0, color="gray", lw=1, ls=":")
    ax_right.text(ax_right.get_xlim()[1] * 0.9, med_i,
                  f" median ImCoh\n {med_i:.2f}×",
                  color=CLR_IMCOH, fontsize=10, fontweight="bold",
                  ha="right", va="bottom")
    ax_right.set_yscale("log")
    ax_right.tick_params(labelleft=False, bottom=False, labelbottom=False)
    for spine in ("top", "right", "bottom"):
        ax_right.spines[spine].set_visible(False)

    # Patient legend (placed in empty top-right corner)
    pat_handles = [Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=PATIENT_COLORS[p], markersize=9,
                          label=p) for p in patients]
    fig.legend(handles=pat_handles, loc="upper right",
               fontsize=9, title="Patient",
               bbox_to_anchor=(0.99, 0.99), frameon=True)

    # Annotation in the main panel
    ax_scatter.text(0.02, 0.97,
                    "points below y = x: ImCoh reduced the probe bias",
                    transform=ax_scatter.transAxes,
                    fontsize=9, ha="left", va="top",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white",
                              ec="0.6", alpha=0.9))

    save_fig(fig, output_dir / f"fig_D1_paired_dots_{phase}")
    return df


# ---------------------------------------------------------------------------
# Figure D2 — Bias-reduction heatmap at the LRG level
# ---------------------------------------------------------------------------

def fig_d2_bias_reduction(
    patients: list[str],
    band: str = "beta",
    phase: str = "rest_pre",
    output_dir: Path = SECTION2_ROOT / "fig_D",
    verbose: bool = False,
):
    """Main panel: heatmap of log2(MSC_enrichment / ImCoh_enrichment) per
    (patient, scale n). Blue = ImCoh much better, near-zero = tied, red = ImCoh worse.

    Side panel: mean enrichment vs scale per method, with ±1SD shading.
    """
    # Collect enrichment at each scale, for each patient and method
    import math
    records = {m: {pat: {} for pat in patients} for m in METHODS}

    for method in METHODS:
        for pat in patients:
            ch = load_channel_labels(pat)
            probes = extract_probe_labels(ch)
            lrg = load_all_lrg(
                patients=[pat], bands=[band], phases=[phase],
                fc_method=method, verbose=False,
            )
            key = (pat, phase, band)
            if key not in lrg:
                continue
            enrich = compute_enrichment_vs_scale(
                lrg[key].linkage_matrix, probes,
                n_communities=N_COMMUNITIES, n_nodes=lrg[key].n_nodes,
            )
            records[method][pat] = enrich

    # Build MSC matrix (patients x n_scales) and ImCoh matrix
    n_vals = N_COMMUNITIES
    msc_mat = np.full((len(patients), len(n_vals)), np.nan)
    im_mat = np.full((len(patients), len(n_vals)), np.nan)
    for i, pat in enumerate(patients):
        for j, n in enumerate(n_vals):
            if n in records["msc"].get(pat, {}):
                msc_mat[i, j] = records["msc"][pat][n]
            if n in records["imcoh_abs"].get(pat, {}):
                im_mat[i, j] = records["imcoh_abs"][pat][n]

    # Bias-reduction ratio: log2(msc / imcoh) — positive = ImCoh cleaner
    with np.errstate(divide="ignore", invalid="ignore"):
        reduction = np.log2(msc_mat / im_mat)

    fig = plt.figure(figsize=(16, 6.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[2.5, 1.8, 1.8], wspace=0.35)

    # --- Main panel: bias-reduction heatmap ---
    ax = fig.add_subplot(gs[0, 0])
    cmap = plt.get_cmap("RdBu")
    vmax_abs = np.nanmax(np.abs(reduction)) if np.isfinite(reduction).any() else 3.0
    im = ax.imshow(reduction, aspect="auto", cmap=cmap,
                   vmin=-vmax_abs, vmax=vmax_abs, interpolation="nearest")
    for i in range(len(patients)):
        for j in range(len(n_vals)):
            v = reduction[i, j]
            if np.isfinite(v):
                clr = "white" if abs(v) > vmax_abs * 0.6 else "black"
                ax.text(j, i, f"{2**v:.1f}×", ha="center", va="center",
                        fontsize=9, color=clr, fontweight="bold")
    ax.set_xticks(range(len(n_vals)))
    ax.set_xticklabels([str(n) for n in n_vals], fontsize=11)
    ax.set_yticks(range(len(patients)))
    ax.set_yticklabels(patients, fontsize=11)
    ax.set_xlabel("Number of LRG communities (n)", fontsize=12)
    ax.set_title(r"Bias reduction factor: MSC enrichment / $|\mathrm{ImCoh}|$ enrichment"
                 "\n"
                 r"(values > 1 → $|\mathrm{ImCoh}|$ is less probe-biased)",
                 fontsize=12, fontweight="bold")
    cbar = plt.colorbar(im, ax=ax, shrink=0.85,
                        label="log₂ ratio")

    # --- Middle panel: mean ± SD curves ---
    ax2 = fig.add_subplot(gs[0, 1])
    for method, clr, label in [("msc", CLR_MSC, "MSC"),
                                ("imcoh_abs", CLR_IMCOH, r"$|\mathrm{ImCoh}|$")]:
        mat = msc_mat if method == "msc" else im_mat
        mu = np.nanmean(mat, axis=0)
        sd = np.nanstd(mat, axis=0)
        ax2.plot(n_vals, mu, "o-", color=clr, lw=2.5, markersize=7,
                 label=label)
        ax2.fill_between(n_vals, mu - sd, mu + sd, color=clr, alpha=0.18)
    ax2.axhline(1.0, color="black", ls=":", lw=1.2)
    ax2.axhspan(0.8, 1.2, color="#4CAF50", alpha=0.1, zorder=0)
    ax2.set_xscale("log")
    ax2.set_xticks(n_vals)
    ax2.set_xticklabels([str(n) for n in n_vals])
    ax2.set_xlabel("n", fontsize=11)
    ax2.set_ylabel("community-probe enrichment", fontsize=11)
    ax2.set_title("Mean ± SD across patients", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=10, loc="upper left")
    ax2.grid(alpha=0.25, which="both")

    # --- Right panel: distribution of reduction factors per scale ---
    ax3 = fig.add_subplot(gs[0, 2])
    for j, n in enumerate(n_vals):
        col = reduction[:, j]
        col = col[np.isfinite(col)]
        if col.size == 0:
            continue
        # Points
        x_pos = j + np.random.default_rng(j).uniform(-0.12, 0.12, size=col.size)
        ax3.scatter(x_pos, 2**col, s=30, alpha=0.7, color="#555",
                    edgecolors="white", linewidths=0.5)
        # Median bar
        ax3.hlines(2**np.median(col), j - 0.25, j + 0.25,
                   color="black", lw=2.5)
    ax3.axhline(1.0, color="red", ls=":", lw=1)
    ax3.set_xticks(range(len(n_vals)))
    ax3.set_xticklabels([str(n) for n in n_vals])
    ax3.set_yscale("log")
    ax3.set_xlabel("n", fontsize=11)
    ax3.set_ylabel(r"MSC / $|\mathrm{ImCoh}|$ enrichment ratio", fontsize=11)
    ax3.set_title("Per-scale bias-reduction distribution",
                  fontsize=11, fontweight="bold")
    ax3.grid(axis="y", which="both", alpha=0.25)

    band_tex = BRAIN_BAND_TEX_DICT[band]
    save_fig(fig, output_dir / f"fig_D2_bias_reduction_{band}_{phase}")


def main():
    parser = argparse.ArgumentParser(description="Section 2 enrichment figures.")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument("--output-dir", type=Path, default=SECTION2_ROOT)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()
    out = args.output_dir / "fig_D"

    # D1 for representative phases (variety: rest + task)
    for phase in ["rest_pre", "task_learn", "rest_post"]:
        print(f"--- Figure D1: Paired dots ({phase}) ---")
        fig_d1_paired_dots(ALL_PATIENTS, phase, out, args.verbose)

    # D2 for multiple bands and phases
    for band in ["alpha", "beta"]:
        for phase in ["rest_pre", "task_learn"]:
            print(f"\n--- Figure D2: Bias reduction ({band}, {phase}) ---")
            fig_d2_bias_reduction(ALL_PATIENTS, band, phase, out, args.verbose)

    print("\nDone.")


if __name__ == "__main__":
    main()
