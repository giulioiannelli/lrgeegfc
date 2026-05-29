#!/usr/bin/env python3
"""Step 1: Raw exploration of mean_VI distances.

One figure per band showing the 6 pairwise distances for each of the 4
patients with complete phases.  No aggregation, no z-scoring — just the
raw numbers so we can see where agreement exists and where outliers are.

Produces:
  data/figures/metric_exploration/raw_vi_exploration/
    raw_vi_{band}.pdf           — per-band: 4 patients, 6 pairs
    raw_vi_hypothesis_check.pdf — full summary table
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

# ── Load data ──────────────────────────────────────────────────────────
CSV = FIGURES_ROOT / "metric_exploration" / "partition_multiscale" / "results.csv"
df = pd.read_csv(CSV)

OUTDIR = FIGURES_ROOT / "metric_exploration" / "raw_vi_exploration"
OUTDIR.mkdir(parents=True, exist_ok=True)

PATIENTS = PATIENTS_4PHASE
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_LABELS = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}

# Phase pair ordering: within pairs first, then cross
PAIRS_ORDERED = [
    "task_learn-task_test",   # within (task)
    "rest_pre-rest_post",         # within (rest)
    "rest_pre-task_learn",      # cross
    "rest_pre-task_test",       # cross
    "task_learn-rest_post",     # cross
    "task_test-rest_post",      # cross (the task-trace pair)
]
PAIR_SHORT = {
    "task_learn-task_test": "TL↔TT",
    "rest_pre-rest_post": "Pre↔Post",
    "rest_pre-task_learn": "Pre↔TL",
    "rest_pre-task_test": "Pre↔TT",
    "task_learn-rest_post": "TL↔Post",
    "task_test-rest_post": "TT↔Post",
}
PAIR_COLORS = {
    "task_learn-task_test": "#2196F3",    # blue  — within task
    "rest_pre-rest_post": "#4CAF50",          # green — within rest
    "rest_pre-task_learn": "#9E9E9E",       # grey  — cross
    "rest_pre-task_test": "#BDBDBD",        # light grey
    "task_learn-rest_post": "#FF9800",      # orange — cross (task→rest)
    "task_test-rest_post": "#E91E63",       # pink — the key task-trace pair
}

PAT_MARKERS = {"Pat_02": "o", "Pat_03": "s", "Pat_05": "D", "Pat_08": "^"}

# ── Figure 1: One panel per band ──────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True)
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    sub = df[(df.band == band) & (df.patient.isin(PATIENTS))]

    x_positions = np.arange(len(PAIRS_ORDERED))

    for ip, pair in enumerate(PAIRS_ORDERED):
        for jp, pat in enumerate(PATIENTS):
            row = sub[(sub.patient == pat) & (sub.pair == pair)]
            if len(row) == 0:
                continue
            val = row.mean_VI.values[0]
            offset = (jp - 1.5) * 0.12
            ax.plot(
                ip + offset, val,
                marker=PAT_MARKERS[pat],
                color=PAIR_COLORS[pair],
                markersize=8, markeredgecolor="k", markeredgewidth=0.5,
                zorder=3,
            )

    # Shading for within vs cross
    ax.axvspan(-0.5, 1.5, color="#E3F2FD", alpha=0.3, zorder=0)
    ax.axvspan(1.5, 5.5, color="#FFF3E0", alpha=0.3, zorder=0)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([PAIR_SHORT[p] for p in PAIRS_ORDERED], rotation=45, ha="right", fontsize=8)
    ax.set_title(BAND_LABELS[band], fontsize=14)
    ax.set_ylabel("mean VI (distance)" if idx % 3 == 0 else "")
    ax.set_ylim(0, max(1.2, sub.mean_VI.max() * 1.1))
    ax.grid(axis="y", alpha=0.3)

    # H1 check: is TL-TT the minimum for each patient?
    for pat in PATIENTS:
        pat_sub = sub[sub.patient == pat].set_index("pair")
        if "task_learn-task_test" not in pat_sub.index:
            continue
        tl_tt = pat_sub.loc["task_learn-task_test", "mean_VI"]
        is_min = tl_tt == pat_sub.mean_VI.min()
        if not is_min:
            rank = sorted(pat_sub.mean_VI.values.tolist()).index(tl_tt) + 1
            # Mark the outlier
            offset = (PATIENTS.index(pat) - 1.5) * 0.12
            ax.annotate(
                f"{pat[-2:]}:R{rank}",
                (0 + offset, tl_tt), fontsize=6,
                xytext=(0, 8), textcoords="offset points",
                ha="center", color="red",
            )

# Legend
from matplotlib.lines import Line2D
handles = []
for pat in PATIENTS:
    handles.append(Line2D([0], [0], marker=PAT_MARKERS[pat], color="gray",
                          label=pat, markersize=7, linestyle="none"))
handles.append(Line2D([0], [0], color="none", label=""))
handles.append(Line2D([0], [0], marker="o", color="#2196F3", label="Within (task)", markersize=7, linestyle="none"))
handles.append(Line2D([0], [0], marker="o", color="#4CAF50", label="Within (rest)", markersize=7, linestyle="none"))
handles.append(Line2D([0], [0], marker="o", color="#9E9E9E", label="Cross", markersize=7, linestyle="none"))
handles.append(Line2D([0], [0], marker="o", color="#E91E63", label="TT↔Post (task trace)", markersize=7, linestyle="none"))

fig.legend(handles=handles, loc="lower center", ncol=8, fontsize=8,
           bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Raw mean VI distances — lower = more similar\nBlue shading = within pairs, Orange shading = cross pairs",
             fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(OUTDIR / "raw_vi_all_bands.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved raw_vi_all_bands.pdf")

# ── Figure 2: Hypothesis check — horizontal bar plot per patient ──────
# For each (patient, band): show the 3 key contrasts as signed bars
# H1: mean(cross) - TL↔TT  (positive = H1 holds)
# H3: mean(cross) - mean(within)  (positive = H3 holds)
# H2a: Pre↔Post - TT↔Post  (positive = task persists in rest_post)
# H2b: Pre↔TT - TT↔Post    (positive = task persists)

fig, axes = plt.subplots(4, 1, figsize=(14, 16), sharex=True)

for ip, pat in enumerate(PATIENTS):
    ax = axes[ip]
    sub = df[(df.patient == pat)].set_index(["band", "pair"])

    y_positions = []
    colors_h1, colors_h3, colors_h2a, colors_h2b = [], [], [], []
    vals_h1, vals_h3, vals_h2a, vals_h2b = [], [], [], []

    for ib, band in enumerate(BANDS):
        b = sub.loc[band]
        tl_tt = b.loc["task_learn-task_test", "mean_VI"]
        pre_post = b.loc["rest_pre-rest_post", "mean_VI"]
        tt_post = b.loc["task_test-rest_post", "mean_VI"]
        pre_tt = b.loc["rest_pre-task_test", "mean_VI"]
        pre_tl = b.loc["rest_pre-task_learn", "mean_VI"]
        tl_post = b.loc["task_learn-rest_post", "mean_VI"]

        cross_mean = (pre_tl + pre_tt + tl_post + tt_post) / 4
        within_mean = (tl_tt + pre_post) / 2

        # H1: task stability. cross_mean - tl_tt > 0 means TL-TT is closer than average cross
        h1 = cross_mean - tl_tt
        # H3: within < cross → gap = cross_mean - within_mean > 0
        h3 = cross_mean - within_mean
        # H2a: pre_post - tt_post > 0 → rest_post closer to task than to rest_pre
        h2a = pre_post - tt_post
        # H2b: pre_tt - tt_post > 0 → same idea from different angle
        h2b = pre_tt - tt_post

        vals_h1.append(h1)
        vals_h3.append(h3)
        vals_h2a.append(h2a)
        vals_h2b.append(h2b)

    y = np.arange(len(BANDS))
    bar_h = 0.2
    ax.barh(y + 1.5*bar_h, vals_h3, bar_h, label="H3: within < cross", color="#2196F3", alpha=0.8)
    ax.barh(y + 0.5*bar_h, vals_h1, bar_h, label="H1: TL↔TT most similar", color="#4CAF50", alpha=0.8)
    ax.barh(y - 0.5*bar_h, vals_h2a, bar_h, label="H2a: Pre↔Post − TT↔Post", color="#FF9800", alpha=0.8)
    ax.barh(y - 1.5*bar_h, vals_h2b, bar_h, label="H2b: Pre↔TT − TT↔Post", color="#E91E63", alpha=0.8)

    ax.axvline(0, color="k", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([BAND_LABELS[b] for b in BANDS], fontsize=11)
    ax.set_title(pat, fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)
    if ip == 0:
        ax.legend(loc="upper right", fontsize=8)

axes[-1].set_xlabel("Contrast value (positive = hypothesis supported)", fontsize=10)
fig.suptitle("Hypothesis contrasts per patient — raw mean VI\nPositive bars = hypothesis holds", fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(OUTDIR / "hypothesis_contrasts_per_patient.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved hypothesis_contrasts_per_patient.pdf")


# ── Figure 3: The key 3 pairs across bands — one panel per pair ───────
# Show the actual distances, not contrasts. Connect same patient across bands.
KEY_PAIRS = ["task_learn-task_test", "rest_pre-rest_post", "task_test-rest_post"]
KEY_LABELS = ["TL↔TT (task stability)", "Pre↔Post (rest change)", "TT↔Post (task trace)"]

fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

for ip, (pair, label) in enumerate(zip(KEY_PAIRS, KEY_LABELS)):
    ax = axes[ip]
    for pat in PATIENTS:
        vals = []
        for band in BANDS:
            row = df[(df.patient == pat) & (df.band == band) & (df.pair == pair)]
            vals.append(row.mean_VI.values[0] if len(row) else np.nan)
        ax.plot(range(len(BANDS)), vals, marker=PAT_MARKERS[pat],
                label=pat, linewidth=1.5, markersize=7, alpha=0.8)

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BAND_LABELS[b] for b in BANDS], fontsize=10)
    ax.set_title(label, fontsize=11, fontweight="bold")
    ax.set_ylabel("mean VI distance" if ip == 0 else "")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

fig.suptitle("Key pairs across frequency bands — raw mean VI distances\n(lower = more similar)",
             fontsize=12, y=1.03)
fig.tight_layout()
fig.savefig(OUTDIR / "key_pairs_across_bands.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved key_pairs_across_bands.pdf")


# ── Figure 4: Phase×phase heatmap per patient per band ────────────────
# This is the most complete view: 4 patients × 6 bands = 24 small heatmaps
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]

fig, axes = plt.subplots(4, 6, figsize=(18, 12))

for ip, pat in enumerate(PATIENTS):
    for ib, band in enumerate(BANDS):
        ax = axes[ip, ib]
        mat = np.zeros((4, 4))
        sub = df[(df.patient == pat) & (df.band == band)]
        for _, row in sub.iterrows():
            i = PHASES.index(row.phase1)
            j = PHASES.index(row.phase2)
            mat[i, j] = row.mean_VI
            mat[j, i] = row.mean_VI

        im = ax.imshow(mat, cmap="RdYlGn_r", vmin=0, vmax=1.2, aspect="equal")
        # Annotate
        for i in range(4):
            for j in range(4):
                if i != j:
                    ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", fontsize=6,
                            color="white" if mat[i,j] > 0.7 else "black")

        if ip == 0:
            ax.set_title(BAND_LABELS[band], fontsize=11)
        if ib == 0:
            ax.set_ylabel(pat, fontsize=10, fontweight="bold")

        phase_short = ["Pre", "TL", "TT", "Post"]
        ax.set_xticks(range(4))
        ax.set_xticklabels(phase_short if ip == 3 else [], fontsize=7)
        ax.set_yticks(range(4))
        ax.set_yticklabels(phase_short if ib == 0 else [], fontsize=7)

fig.suptitle("Phase × Phase distance matrices (mean VI) — lower = more similar\nGreen = similar, Red = different",
             fontsize=13, y=1.02)
fig.tight_layout()
fig.savefig(OUTDIR / "phase_heatmaps_all.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved phase_heatmaps_all.pdf")


# ── Print summary table ───────────────────────────────────────────────
print("\n" + "="*90)
print("SUMMARY: Which hypotheses hold per band? (counting patients)")
print("="*90)
print(f"{'Band':>12s}  {'H1: TL-TT min':>16s}  {'H3: within<cross':>16s}  {'H2a: task trace':>16s}  {'H2b: task trace':>16s}")
print("-"*90)

for band in BANDS:
    h1_count = 0
    h3_count = 0
    h2a_pos = 0
    h2a_neg = 0
    h2b_pos = 0
    h2b_neg = 0
    for pat in PATIENTS:
        sub = df[(df.patient == pat) & (df.band == band)].set_index("pair")
        tl_tt = sub.loc["task_learn-task_test", "mean_VI"]
        pre_post = sub.loc["rest_pre-rest_post", "mean_VI"]
        tt_post = sub.loc["task_test-rest_post", "mean_VI"]
        pre_tt = sub.loc["rest_pre-task_test", "mean_VI"]
        pre_tl = sub.loc["rest_pre-task_learn", "mean_VI"]
        tl_post = sub.loc["task_learn-rest_post", "mean_VI"]

        # H1
        if tl_tt == min(tl_tt, pre_post, pre_tl, pre_tt, tl_post, tt_post):
            h1_count += 1
        # H3
        cross_m = (pre_tl + pre_tt + tl_post + tt_post) / 4
        within_m = (tl_tt + pre_post) / 2
        if within_m < cross_m:
            h3_count += 1
        # H2a
        if pre_post > tt_post:
            h2a_pos += 1
        else:
            h2a_neg += 1
        # H2b
        if pre_tt > tt_post:
            h2b_pos += 1
        else:
            h2b_neg += 1

    h1_str = f"{h1_count}/4"
    h3_str = f"{h3_count}/4"
    h2a_str = f"+{h2a_pos}/-{h2a_neg}"
    h2b_str = f"+{h2b_pos}/-{h2b_neg}"
    print(f"{BAND_LABELS[band]:>12s}  {h1_str:>16s}  {h3_str:>16s}  {h2a_str:>16s}  {h2b_str:>16s}")

print()
print("H2a: Pre↔Post − TT↔Post > 0 means task persists in rest_post")
print("H2b: Pre↔TT  − TT↔Post > 0 means task persists in rest_post")
print("+ = persists, - = recovers")
