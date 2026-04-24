#!/usr/bin/env python3
"""Multiscale H2 task-trace profiles.

For each (patient, band, k):
  H2a(k) = VI(Pre,Post; k) - VI(TT,Post; k)
  H2b(k) = VI(Pre,TT; k)  - VI(TT,Post; k)

  Positive = task persists in rest_post at scale k.
  Negative = brain recovers to pre-task at scale k.

We plot these as continuous curves over k, one per band, with:
  - Individual patient traces (thin lines)
  - Mean across patients (thick line)
  - Shading where ALL patients agree on sign (unanimous region)

This reveals WHERE in the hierarchy each band's task trace lives,
without picking arbitrary scale cutoffs.

Output:
  data/figures/metric_exploration/multiscale_h2_profiles/
    h2a_profiles.pdf    — H2a: Pre↔Post vs TT↔Post
    h2b_profiles.pdf    — H2b: Pre↔TT vs TT↔Post
    h2_combined.pdf     — Both H2a and H2b, band × scale heatmap of mean contrast
    h2_unanimity_map.pdf — Band × k heatmap where color = sign × unanimity
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.paths import FIGURES_ROOT

# ── Load ───────────────────────────────────────────────────────────────
CSV = FIGURES_ROOT / "metric_exploration" / "partition_multiscale" / "results.csv"
df = pd.read_csv(CSV)

OUTDIR = FIGURES_ROOT / "metric_exploration" / "multiscale_h2_profiles"
OUTDIR.mkdir(parents=True, exist_ok=True)

from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS = PATIENTS_4PHASE
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_LABELS = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}

K_VALUES = [2, 3, 4, 5, 6, 8, 10, 15, 20]
K_COLS = [f"VI_k{k}" for k in K_VALUES]


def get_vi_vector(patient: str, band: str, pair: str) -> np.ndarray:
    """Return VI values at all k for a specific (patient, band, pair)."""
    row = df[(df.patient == patient) & (df.band == band) & (df.pair == pair)]
    if len(row) == 0:
        return np.full(len(K_VALUES), np.nan)
    return row[K_COLS].values[0]


def compute_h2_curves():
    """Compute H2a(k) and H2b(k) for all (patient, band)."""
    records = []
    for pat in PATIENTS:
        for band in BANDS:
            vi_pre_post = get_vi_vector(pat, band, "rest_pre-rest_post")
            vi_tt_post = get_vi_vector(pat, band, "task_test-rest_post")
            vi_pre_tt = get_vi_vector(pat, band, "rest_pre-task_test")

            h2a = vi_pre_post - vi_tt_post  # positive = task persists
            h2b = vi_pre_tt - vi_tt_post    # positive = task persists

            for ik, k in enumerate(K_VALUES):
                records.append({
                    "patient": pat, "band": band, "k": k,
                    "h2a": h2a[ik], "h2b": h2b[ik],
                })
    return pd.DataFrame(records)


h2 = compute_h2_curves()


# ── Figure 1 & 2: H2a and H2b profile curves ─────────────────────────
for contrast, label, fname in [
    ("h2a", "H2a: VI(Pre,Post) − VI(TT,Post)", "h2a_profiles.pdf"),
    ("h2b", "H2b: VI(Pre,TT) − VI(TT,Post)", "h2b_profiles.pdf"),
]:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
    axes = axes.ravel()

    for ib, band in enumerate(BANDS):
        ax = axes[ib]
        color = BAND_COLORS[band]

        # Individual patient traces
        for pat in PATIENTS:
            sub = h2[(h2.patient == pat) & (h2.band == band)]
            ax.plot(sub.k, sub[contrast], color=color, alpha=0.25,
                    linewidth=1, marker=".", markersize=4)

        # Mean across patients
        mean_curve = h2[h2.band == band].groupby("k")[contrast].mean()
        ax.plot(mean_curve.index, mean_curve.values, color=color,
                linewidth=2.5, marker="o", markersize=6, zorder=5)

        # Shade unanimous regions
        for k in K_VALUES:
            vals = h2[(h2.band == band) & (h2.k == k)][contrast].values
            if len(vals) == 4:
                if all(v > 0 for v in vals):
                    ax.axvspan(
                        k - 0.4, k + 0.4, color="#4CAF50", alpha=0.15, zorder=0,
                    )
                elif all(v < 0 for v in vals):
                    ax.axvspan(
                        k - 0.4, k + 0.4, color="#F44336", alpha=0.15, zorder=0,
                    )

        ax.axhline(0, color="k", linewidth=0.8, linestyle="--")
        ax.set_title(BAND_LABELS[band], fontsize=13, fontweight="bold")
        ax.set_xlabel("k (number of communities)" if ib >= 3 else "")
        ax.set_ylabel("Contrast (+ = task persists)" if ib % 3 == 0 else "")
        ax.set_xticks(K_VALUES)
        ax.set_xticklabels([str(k) for k in K_VALUES], fontsize=8)
        ax.grid(alpha=0.3)

    # Legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color="gray", linewidth=2.5, marker="o", label="Mean (4 patients)"),
        Line2D([0], [0], color="gray", linewidth=1, alpha=0.4, label="Individual patients"),
        Patch(facecolor="#4CAF50", alpha=0.3, label="Unanimous: task persists"),
        Patch(facecolor="#F44336", alpha=0.3, label="Unanimous: brain recovers"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"{label}\nPositive = rest_post retains task structure | "
                 f"Green shading = all 4 patients agree (+) | Red = all agree (−)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTDIR / fname, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {fname}")


# ── Figure 3: Combined heatmap — mean contrast as function of (band, k) ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ic, (contrast, title) in enumerate([
    ("h2a", "H2a: VI(Pre,Post) − VI(TT,Post)"),
    ("h2b", "H2b: VI(Pre,TT) − VI(TT,Post)"),
]):
    ax = axes[ic]
    mat = np.zeros((len(BANDS), len(K_VALUES)))
    for ib, band in enumerate(BANDS):
        for ik, k in enumerate(K_VALUES):
            vals = h2[(h2.band == band) & (h2.k == k)][contrast].values
            mat[ib, ik] = np.mean(vals)

    vmax = max(0.3, np.abs(mat).max())
    im = ax.imshow(mat, cmap="RdBu", vmin=-vmax, vmax=vmax, aspect="auto")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Annotate with mean values
    for ib in range(len(BANDS)):
        for ik in range(len(K_VALUES)):
            val = mat[ib, ik]
            ax.text(ik, ib, f"{val:+.2f}", ha="center", va="center",
                    fontsize=6, color="white" if abs(val) > vmax * 0.6 else "black")

    ax.set_xticks(range(len(K_VALUES)))
    ax.set_xticklabels([str(k) for k in K_VALUES])
    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BAND_LABELS[b] for b in BANDS])
    ax.set_xlabel("k (communities)")
    ax.set_title(title, fontsize=10)

fig.suptitle("Mean H2 contrast across 4 patients\nBlue = task persists | Red = brain recovers",
             fontsize=12, y=1.05)
fig.tight_layout()
fig.savefig(OUTDIR / "h2_combined_heatmap.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved h2_combined_heatmap.pdf")


# ── Figure 4: Unanimity map — sign × count ───────────────────────────
# Value = number of patients with positive contrast (0-4)
# So: 4 = all persist, 0 = all recover, 2 = split
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ic, (contrast, title) in enumerate([
    ("h2a", "H2a: task persistence in rest_post"),
    ("h2b", "H2b: task persistence in rest_post"),
]):
    ax = axes[ic]
    mat = np.zeros((len(BANDS), len(K_VALUES)))
    for ib, band in enumerate(BANDS):
        for ik, k in enumerate(K_VALUES):
            vals = h2[(h2.band == band) & (h2.k == k)][contrast].values
            n_pos = np.sum(vals > 0)
            mat[ib, ik] = n_pos

    # Custom colormap: 0=red (all recover), 2=white (split), 4=blue (all persist)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("unanimity",
        ["#D32F2F", "#EF9A9A", "white", "#90CAF9", "#1565C0"], N=5)

    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=4, aspect="auto")
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, ticks=[0, 1, 2, 3, 4])
    cbar.set_ticklabels(["0/4\nrecovers", "1/4", "2/4\nsplit", "3/4", "4/4\npersists"])

    # Annotate: bold stars for unanimous
    for ib in range(len(BANDS)):
        for ik in range(len(K_VALUES)):
            val = int(mat[ib, ik])
            sym = "★" if val in (0, 4) else str(val)
            ax.text(ik, ib, sym, ha="center", va="center",
                    fontsize=9 if val in (0, 4) else 7,
                    fontweight="bold" if val in (0, 4) else "normal",
                    color="white" if val in (0, 4) else "black")

    ax.set_xticks(range(len(K_VALUES)))
    ax.set_xticklabels([str(k) for k in K_VALUES])
    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BAND_LABELS[b] for b in BANDS])
    ax.set_xlabel("k (communities)")
    ax.set_title(title, fontsize=10)

fig.suptitle("Unanimity map: how many patients show task persistence at each (band, k)\n"
             "★ = all 4 patients agree",
             fontsize=12, y=1.05)
fig.tight_layout()
fig.savefig(OUTDIR / "h2_unanimity_map.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved h2_unanimity_map.pdf")


# ── Figure 5: The key figure — continuous H2a profiles, all bands on same axes ──
# This is the "money plot": one panel, all 6 bands, showing where each
# band's task trace lives in the hierarchy
fig, ax = plt.subplots(figsize=(10, 6))

for band in BANDS:
    mean_curve = h2[h2.band == band].groupby("k")["h2a"].mean()
    ax.plot(mean_curve.index, mean_curve.values, color=BAND_COLORS[band],
            linewidth=2.5, marker="o", markersize=7, label=BAND_LABELS[band],
            zorder=5)

    # Error band: min-max across patients
    for k in K_VALUES:
        vals = h2[(h2.band == band) & (h2.k == k)]["h2a"].values
        ax.fill_between([k - 0.15, k + 0.15],
                        [np.min(vals)] * 2, [np.max(vals)] * 2,
                        color=BAND_COLORS[band], alpha=0.1)

ax.axhline(0, color="k", linewidth=1, linestyle="--", zorder=0)
ax.set_xlabel("k (number of communities)", fontsize=12)
ax.set_ylabel("H2a contrast: VI(Pre,Post) − VI(TT,Post)", fontsize=12)
ax.set_title("Task trace across hierarchical scales\n"
             "Positive = rest_post retains task structure | Negative = brain recovers",
             fontsize=13)
ax.set_xticks(K_VALUES)
ax.set_xticklabels([str(k) for k in K_VALUES])
ax.legend(fontsize=10, loc="upper left")
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUTDIR / "h2a_all_bands_overlay.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved h2a_all_bands_overlay.pdf")


# ── Figure 6: Same but for H2b ───────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))

for band in BANDS:
    mean_curve = h2[h2.band == band].groupby("k")["h2b"].mean()
    ax.plot(mean_curve.index, mean_curve.values, color=BAND_COLORS[band],
            linewidth=2.5, marker="o", markersize=7, label=BAND_LABELS[band],
            zorder=5)

    for k in K_VALUES:
        vals = h2[(h2.band == band) & (h2.k == k)]["h2b"].values
        ax.fill_between([k - 0.15, k + 0.15],
                        [np.min(vals)] * 2, [np.max(vals)] * 2,
                        color=BAND_COLORS[band], alpha=0.1)

ax.axhline(0, color="k", linewidth=1, linestyle="--", zorder=0)
ax.set_xlabel("k (number of communities)", fontsize=12)
ax.set_ylabel("H2b contrast: VI(Pre,TT) − VI(TT,Post)", fontsize=12)
ax.set_title("Task trace across hierarchical scales (H2b)\n"
             "Positive = rest_post retains task structure | Negative = brain recovers",
             fontsize=13)
ax.set_xticks(K_VALUES)
ax.set_xticklabels([str(k) for k in K_VALUES])
ax.legend(fontsize=10, loc="upper left")
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUTDIR / "h2b_all_bands_overlay.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved h2b_all_bands_overlay.pdf")


# ── Print summary ─────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("MULTISCALE H2 PROFILE SUMMARY")
print("=" * 70)
print("\nUnanimous (★) cells in H2a (all 4 patients agree on sign):")
for band in BANDS:
    stars = []
    for k in K_VALUES:
        vals = h2[(h2.band == band) & (h2.k == k)]["h2a"].values
        if all(v > 0 for v in vals):
            stars.append(f"k={k}(+)")
        elif all(v < 0 for v in vals):
            stars.append(f"k={k}(−)")
    if stars:
        print(f"  {BAND_LABELS[band]:>12s}: {', '.join(stars)}")
    else:
        print(f"  {BAND_LABELS[band]:>12s}: none")

print("\nUnanimous (★) cells in H2b:")
for band in BANDS:
    stars = []
    for k in K_VALUES:
        vals = h2[(h2.band == band) & (h2.k == k)]["h2b"].values
        if all(v > 0 for v in vals):
            stars.append(f"k={k}(+)")
        elif all(v < 0 for v in vals):
            stars.append(f"k={k}(−)")
    if stars:
        print(f"  {BAND_LABELS[band]:>12s}: {', '.join(stars)}")
    else:
        print(f"  {BAND_LABELS[band]:>12s}: none")

print("\nDone!")
