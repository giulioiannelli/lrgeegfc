#!/usr/bin/env python3
"""Multiscale hypothesis testing: H1–H4 as functions of k.

For each integer k from 2 to 30, cut every dendrogram at k clusters
and compute VI between phase pairs. Then test each hypothesis at each k.

x-axis is always k (number of communities) — directly interpretable and
identical across dendrograms (unlike threshold height h).

Produces one figure per hypothesis:
  H1: Task stability — is TL↔TT the most similar pair?
  H2: Task trace — does rest_post carry the task structure?
  H3: Within vs cross — are same-type phases more similar?
  H4: Frequency gradient — which bands reorganize most?

Output: data/figures/metric_exploration/multiscale_hypotheses/
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.cluster.hierarchy import fcluster

import argparse
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.paths import FIGURES_ROOT, lrg_cache_for

_p = argparse.ArgumentParser(add_help=False)
_p.add_argument("--fc-method", default="imcoh_abs",
                choices=["msc", "corr", "imcoh_abs", "imcoh_sq"])
_FC_ARGS, _ = _p.parse_known_args()
FC_METHOD = _FC_ARGS.fc_method
CACHE_ROOT = lrg_cache_for(FC_METHOD)

OUTDIR = FIGURES_ROOT / "metric_exploration" / "multiscale_hypotheses" / FC_METHOD
OUTDIR.mkdir(parents=True, exist_ok=True)

from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS = PATIENTS_4PHASE
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}
PAT_MARKERS = {
    "Pat_02": "o", "Pat_03": "s", "Pat_05": "D", "Pat_06": "v",
    "Pat_07": "^", "Pat_08": "<", "Pat_10": ">", "Pat_13": "p",
    "Pat_14": "h", "Pat_15": "*",
}

PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
ALL_PAIRS = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("rest_pre", "rest_post"),
    ("task_learn", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]
PAIR_LABELS = {
    ("rest_pre", "task_learn"): "Pre↔TL",
    ("rest_pre", "task_test"): "Pre↔TT",
    ("rest_pre", "rest_post"): "Pre↔Post",
    ("task_learn", "task_test"): "TL↔TT",
    ("task_learn", "rest_post"): "TL↔Post",
    ("task_test", "rest_post"): "TT↔Post",
}
WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [p for p in ALL_PAIRS if p not in WITHIN_PAIRS]

K_RANGE = list(range(2, 31))  # k = 2, 3, ..., 30


from lrg_eegfc.utils.metrics import compute_vi  # canonical implementation


# ── Load all linkage matrices ─────────────────────────────────────────
print("Loading LRG results...")
linkages = {}
for pat in PATIENTS:
    for phase in PHASES:
        for band in BANDS:
            try:
                res = load_lrg_result(pat, phase, band, fc_method=FC_METHOD,
                                      cache_root=CACHE_ROOT)
                linkages[(pat, band, phase)] = res.linkage_matrix
            except Exception:
                pass
print(f"Loaded {len(linkages)} linkages")


# ── Compute VI(k) for all (patient, band, pair, k) ───────────────────
print("Computing VI at each k...")
# vi_data[pat][band][(p1,p2)][k_idx] = VI value
vi_data = {pat: {band: {} for band in BANDS} for pat in PATIENTS}

for pat in PATIENTS:
    for band in BANDS:
        for p1, p2 in ALL_PAIRS:
            k1 = (pat, band, p1)
            k2 = (pat, band, p2)
            if k1 not in linkages or k2 not in linkages:
                continue
            Z1, Z2 = linkages[k1], linkages[k2]
            # Patient may have different channel counts across phases if a
            # phase's recording was partial (e.g. Pat_10 resting=113, task=116).
            # VI across differently-sized label vectors is undefined; skip.
            n1 = Z1.shape[0] + 1
            n2 = Z2.shape[0] + 1
            if n1 != n2:
                continue
            vis = np.zeros(len(K_RANGE))
            for ik, k in enumerate(K_RANGE):
                lab1 = fcluster(Z1, k, criterion="maxclust")
                lab2 = fcluster(Z2, k, criterion="maxclust")
                vis[ik] = compute_vi(lab1, lab2)
            vi_data[pat][band][(p1, p2)] = vis
    print(f"  {pat} done")

K_ARR = np.array(K_RANGE)


# ═══════════════════════════════════════════════════════════════════════
# H1: TASK STABILITY — Is TL↔TT the most similar pair at each k?
# ═══════════════════════════════════════════════════════════════════════
print("\nGenerating H1 figure...")

fig = plt.figure(figsize=(17, 10))
gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.3)

for ib, band in enumerate(BANDS):
    ax = fig.add_subplot(gs[ib // 3, ib % 3])

    # For each patient, plot all 6 pairs with TL↔TT highlighted
    for pat in PATIENTS:
        tl_tt = vi_data[pat][band].get(("task_learn", "task_test"))
        if tl_tt is None:
            continue
        # Plot other pairs in light grey
        for pair in ALL_PAIRS:
            if pair == ("task_learn", "task_test"):
                continue
            v = vi_data[pat][band].get(pair)
            if v is not None:
                ax.plot(K_ARR, v, color="#BDBDBD", linewidth=0.5, alpha=0.4)
        # Plot TL↔TT in blue, bold
        ax.plot(K_ARR, tl_tt, color="#1565C0", linewidth=1.2, alpha=0.5,
                marker=PAT_MARKERS[pat], markersize=3)

    # Mean across patients: TL↔TT vs mean of other 5 pairs.
    # Track the subset of patients that actually contribute (some new patients
    # are missing one phase, so their task_learn/task_test pair may be absent).
    tl_tt_all = []
    others_all = []
    contributing_patients = []
    for pat in PATIENTS:
        if ("task_learn", "task_test") in vi_data[pat][band]:
            tl_tt_all.append(vi_data[pat][band][("task_learn", "task_test")])
            pat_others = []
            for pair in ALL_PAIRS:
                if pair != ("task_learn", "task_test") and pair in vi_data[pat][band]:
                    pat_others.append(vi_data[pat][band][pair])
            if pat_others:
                others_all.append(np.mean(pat_others, axis=0))
            contributing_patients.append(pat)

    if tl_tt_all:
        mean_tl_tt = np.mean(tl_tt_all, axis=0)
        mean_others = np.mean(others_all, axis=0)
        ax.plot(K_ARR, mean_tl_tt, color="#1565C0", linewidth=3, label="TL↔TT (mean)", zorder=5)
        ax.plot(K_ARR, mean_others, color="#757575", linewidth=2.5, linestyle="--",
                label="Other pairs (mean)", zorder=5)

        # Shade where TL↔TT < all others for all contributing patients
        tl_tt_arr = np.array(tl_tt_all)
        unanimous_min = np.ones(len(K_RANGE), dtype=bool)
        for ip, pat in enumerate(contributing_patients):
            all_pair_vals = []
            for pair in ALL_PAIRS:
                v = vi_data[pat][band].get(pair)
                if v is not None:
                    all_pair_vals.append(v)
            if all_pair_vals:
                all_pair_arr = np.array(all_pair_vals)
                pat_min = np.min(all_pair_arr, axis=0)
                unanimous_min &= (tl_tt_all[ip] <= pat_min + 1e-10)

        ax.fill_between(K_ARR, 0, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 2,
                         where=unanimous_min, color="#1565C0", alpha=0.08, zorder=0)

    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_xlabel("k" if ib >= 3 else "")
    ax.set_ylabel("VI (distance)" if ib % 3 == 0 else "")
    ax.set_xlim(2, 30)
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.legend(fontsize=7, loc="upper left")

fig.suptitle(
    "H1: Task stability — Is TL↔TT the most similar pair?\n"
    "Blue = TL↔TT distance, Grey = other pairs | "
    "Blue shading = TL↔TT is minimum for all 4 patients",
    fontsize=12, y=1.01,
)
fig.savefig(OUTDIR / "H1_task_stability.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved H1_task_stability.pdf")


# ═══════════════════════════════════════════════════════════════════════
# H2: TASK TRACE — Does rest_post carry the task structure?
# ═══════════════════════════════════════════════════════════════════════
print("Generating H2 figure...")

# Compute H2a(k) = VI(Pre,Post;k) - VI(TT,Post;k) per (patient, band)
h2a = {pat: {} for pat in PATIENTS}
for pat in PATIENTS:
    for band in BANDS:
        v_pre_post = vi_data[pat][band].get(("rest_pre", "rest_post"))
        v_tt_post = vi_data[pat][band].get(("task_test", "rest_post"))
        if v_pre_post is not None and v_tt_post is not None:
            h2a[pat][band] = v_pre_post - v_tt_post

fig = plt.figure(figsize=(17, 10))
gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.3)

for ib, band in enumerate(BANDS):
    ax = fig.add_subplot(gs[ib // 3, ib % 3])
    color = BAND_COLORS[band]

    # Individual patient traces
    curves = []
    for pat in PATIENTS:
        if band in h2a[pat]:
            c = h2a[pat][band]
            curves.append(c)
            ax.plot(K_ARR, c, color=color, alpha=0.25, linewidth=0.8,
                    marker=PAT_MARKERS[pat], markersize=3)

    if curves:
        curves_arr = np.array(curves)
        mean_c = np.mean(curves_arr, axis=0)
        ax.plot(K_ARR, mean_c, color=color, linewidth=2.5, zorder=5)

        # Shade unanimous regions
        n_pos = np.sum(curves_arr > 0, axis=0)
        n_neg = np.sum(curves_arr < 0, axis=0)
        ylim = max(0.4, np.abs(curves_arr).max() * 1.1)
        # 4/4
        ax.fill_between(K_ARR, -ylim, ylim, where=(n_pos == 4),
                         color="#1565C0", alpha=0.15, zorder=0)
        ax.fill_between(K_ARR, -ylim, ylim, where=(n_neg == 4),
                         color="#B71C1C", alpha=0.15, zorder=0)
        # 3/4 (lighter)
        ax.fill_between(K_ARR, -ylim, ylim, where=(n_pos == 3),
                         color="#1565C0", alpha=0.06, zorder=0)
        ax.fill_between(K_ARR, -ylim, ylim, where=(n_neg == 3),
                         color="#B71C1C", alpha=0.06, zorder=0)
        ax.set_ylim(-ylim, ylim)

    ax.axhline(0, color="k", linewidth=0.8, linestyle="--")
    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_xlabel("k" if ib >= 3 else "")
    ax.set_ylabel("H2a contrast" if ib % 3 == 0 else "")
    ax.set_xlim(2, 30)
    ax.grid(alpha=0.3)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0], [0], color="gray", linewidth=2.5, label="Mean"),
    Patch(facecolor="#1565C0", alpha=0.3, label="4/4 persist"),
    Patch(facecolor="#1565C0", alpha=0.12, label="3/4 persist"),
    Patch(facecolor="#B71C1C", alpha=0.3, label="4/4 recover"),
    Patch(facecolor="#B71C1C", alpha=0.12, label="3/4 recover"),
]
for pat in PATIENTS:
    legend_handles.append(Line2D([0], [0], marker=PAT_MARKERS[pat], color="gray",
                                  linestyle="none", markersize=5, label=pat))
fig.legend(handles=legend_handles, loc="lower center", ncol=9, fontsize=8,
           bbox_to_anchor=(0.5, -0.03))

fig.suptitle(
    "H2: Task trace — VI(Pre,Post) − VI(TT,Post)\n"
    "Positive (blue) = rest_post retains task structure | "
    "Negative (red) = brain recovers to pre-task",
    fontsize=12, y=1.01,
)
fig.savefig(OUTDIR / "H2_task_trace.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved H2_task_trace.pdf")


# ═══════════════════════════════════════════════════════════════════════
# H3: WITHIN < CROSS — Are same-type phases more similar?
# ═══════════════════════════════════════════════════════════════════════
print("Generating H3 figure...")

fig = plt.figure(figsize=(17, 10))
gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.3)

for ib, band in enumerate(BANDS):
    ax = fig.add_subplot(gs[ib // 3, ib % 3])

    # Per patient: mean within vs mean cross
    gap_curves = []
    for pat in PATIENTS:
        within_vals = []
        for p in WITHIN_PAIRS:
            v = vi_data[pat][band].get(p)
            if v is not None:
                within_vals.append(v)
        cross_vals = []
        for p in CROSS_PAIRS:
            v = vi_data[pat][band].get(p)
            if v is not None:
                cross_vals.append(v)
        if within_vals and cross_vals:
            within_mean = np.mean(within_vals, axis=0)
            cross_mean = np.mean(cross_vals, axis=0)
            gap = cross_mean - within_mean  # positive = H3 holds
            gap_curves.append(gap)

            # Thin patient traces
            ax.plot(K_ARR, cross_mean, color="#EF5350", alpha=0.2, linewidth=0.5)
            ax.plot(K_ARR, within_mean, color="#42A5F5", alpha=0.2, linewidth=0.5)

    if gap_curves:
        gap_arr = np.array(gap_curves)

        # Mean within and cross
        all_within = []
        all_cross = []
        for pat in PATIENTS:
            wv = [vi_data[pat][band][p] for p in WITHIN_PAIRS if p in vi_data[pat][band]]
            cv = [vi_data[pat][band][p] for p in CROSS_PAIRS if p in vi_data[pat][band]]
            if wv:
                all_within.append(np.mean(wv, axis=0))
            if cv:
                all_cross.append(np.mean(cv, axis=0))

        mean_within = np.mean(all_within, axis=0)
        mean_cross = np.mean(all_cross, axis=0)

        ax.plot(K_ARR, mean_cross, color="#EF5350", linewidth=2.5, label="Cross (mean)", zorder=5)
        ax.plot(K_ARR, mean_within, color="#42A5F5", linewidth=2.5, label="Within (mean)", zorder=5)

        # Fill between
        ax.fill_between(K_ARR, mean_within, mean_cross,
                         where=mean_cross > mean_within,
                         color="#4CAF50", alpha=0.15, label="Gap (H3 holds)")
        ax.fill_between(K_ARR, mean_within, mean_cross,
                         where=mean_cross <= mean_within,
                         color="#FF9800", alpha=0.15, label="Gap (H3 fails)")

        # Mark where gap is unanimous
        all_pos = np.all(gap_arr > 0, axis=0)
        for ik in range(len(K_RANGE)):
            if all_pos[ik]:
                ax.plot(K_ARR[ik], mean_cross[ik] + 0.02, marker="v", color="#2E7D32",
                        markersize=4, zorder=6)

    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_xlabel("k" if ib >= 3 else "")
    ax.set_ylabel("VI (distance)" if ib % 3 == 0 else "")
    ax.set_xlim(2, 30)
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.legend(fontsize=7, loc="upper left")

fig.suptitle(
    "H3: Within-type vs cross-type similarity\n"
    "Blue = mean within pairs (TL↔TT, Pre↔Post) | "
    "Red = mean cross pairs | Green fill = H3 holds | ▼ = 4/4 unanimous",
    fontsize=12, y=1.01,
)
fig.savefig(OUTDIR / "H3_within_vs_cross.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved H3_within_vs_cross.pdf")


# ═══════════════════════════════════════════════════════════════════════
# H4: FREQUENCY GRADIENT — Which bands reorganize most at each k?
# ═══════════════════════════════════════════════════════════════════════
print("Generating H4 figure...")

fig = plt.figure(figsize=(17, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

# Panel A: Mean cross-pair VI per band (reorganization strength), one line per band
ax = fig.add_subplot(gs[0, 0])
for band in BANDS:
    all_cross_mean = []
    for pat in PATIENTS:
        cv = [vi_data[pat][band][p] for p in CROSS_PAIRS if p in vi_data[pat][band]]
        if cv:
            all_cross_mean.append(np.mean(cv, axis=0))
    if all_cross_mean:
        mean_c = np.mean(all_cross_mean, axis=0)
        ax.plot(K_ARR, mean_c, color=BAND_COLORS[band], linewidth=2.5,
                label=BAND_TEX[band])
        lo = np.min(all_cross_mean, axis=0)
        hi = np.max(all_cross_mean, axis=0)
        ax.fill_between(K_ARR, lo, hi, color=BAND_COLORS[band], alpha=0.08)

ax.set_xlabel("k")
ax.set_ylabel("Mean cross-pair VI")
ax.set_title("A. Reorganization strength per band", fontsize=11, fontweight="bold")
ax.legend(fontsize=8, loc="upper left")
ax.set_xlim(2, 30)
ax.grid(alpha=0.3)

# Panel B: Per-patient rank of bands at each k
ax = fig.add_subplot(gs[0, 1])
# At each k, rank bands by cross-pair mean VI (1=most reorganized)
rank_mat = {band: np.zeros(len(K_RANGE)) for band in BANDS}
for ik, k in enumerate(K_RANGE):
    band_scores = {}
    for band in BANDS:
        scores = []
        for pat in PATIENTS:
            cv = [vi_data[pat][band][p] for p in CROSS_PAIRS if p in vi_data[pat][band]]
            if cv:
                scores.append(np.mean([v[ik] for v in cv]))
        if scores:
            band_scores[band] = np.mean(scores)
    # Rank (higher VI = more reorganized = rank 1)
    sorted_bands = sorted(band_scores.keys(), key=lambda b: band_scores[b], reverse=True)
    for rank, b in enumerate(sorted_bands):
        rank_mat[b][ik] = rank + 1

for band in BANDS:
    ax.plot(K_ARR, rank_mat[band], color=BAND_COLORS[band], linewidth=2,
            label=BAND_TEX[band])
ax.set_xlabel("k")
ax.set_ylabel("Rank (1 = most reorganized)")
ax.set_title("B. Band ranking across scales", fontsize=11, fontweight="bold")
ax.set_ylim(6.5, 0.5)
ax.set_xlim(2, 30)
ax.legend(fontsize=8, loc="center right")
ax.grid(alpha=0.3)

# Panel C: Per-patient profiles at representative scales
for ic, (k_lo, k_hi, title_str) in enumerate([
    (2, 5, "C. Coarse (k=2–5)"),
    (8, 15, "D. Fine (k=8–15)"),
]):
    ax = fig.add_subplot(gs[1, ic])
    k_mask = (K_ARR >= k_lo) & (K_ARR <= k_hi)

    for pat in PATIENTS:
        band_means = []
        for band in BANDS:
            cv = [vi_data[pat][band][p] for p in CROSS_PAIRS if p in vi_data[pat][band]]
            if cv:
                # Average over pairs, then average over k in range
                mean_over_pairs = np.mean(cv, axis=0)
                band_means.append(np.mean(mean_over_pairs[k_mask]))
            else:
                band_means.append(np.nan)
        ax.plot(range(len(BANDS)), band_means, marker=PAT_MARKERS[pat],
                linewidth=1.5, markersize=7, alpha=0.6, label=pat)

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=10)
    ax.set_ylabel("Mean cross VI")
    ax.set_title(title_str, fontsize=11, fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

fig.suptitle(
    "H4: Frequency gradient of reorganization strength\n"
    "Is there a consistent band ordering of reorganization across patients?",
    fontsize=12, y=1.01,
)
fig.savefig(OUTDIR / "H4_frequency_gradient.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved H4_frequency_gradient.pdf")


# ═══════════════════════════════════════════════════════════════════════
# COMBINED SUMMARY: one compact figure showing all 4 hypotheses
# ═══════════════════════════════════════════════════════════════════════
print("Generating combined summary...")

fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(4, 6, hspace=0.45, wspace=0.35,
                        height_ratios=[1, 1, 1, 0.8])

# --- Row 1: H1 per band (TL↔TT rank among 6 pairs) ---
for ib, band in enumerate(BANDS):
    ax = fig.add_subplot(gs[0, ib])
    # For each patient, compute rank of TL↔TT at each k
    rank_curves = []
    for pat in PATIENTS:
        tl_tt = vi_data[pat][band].get(("task_learn", "task_test"))
        if tl_tt is None:
            continue
        ranks = np.zeros(len(K_RANGE))
        for ik in range(len(K_RANGE)):
            all_v = []
            for pair in ALL_PAIRS:
                v = vi_data[pat][band].get(pair)
                if v is not None:
                    all_v.append(v[ik])
            # Rank: 1 = smallest (most similar)
            sorted_v = sorted(all_v)
            ranks[ik] = sorted_v.index(tl_tt[ik]) + 1
        rank_curves.append(ranks)
        ax.plot(K_ARR, ranks, alpha=0.3, marker=PAT_MARKERS[pat], markersize=2,
                color="#1565C0", linewidth=0.7)

    if rank_curves:
        mean_rank = np.mean(rank_curves, axis=0)
        ax.plot(K_ARR, mean_rank, color="#1565C0", linewidth=2.5, zorder=5)
        # Shade where rank=1 for all patients
        all_rank1 = np.all(np.array(rank_curves) == 1, axis=0)
        ax.fill_between(K_ARR, 0.5, 6.5, where=all_rank1,
                         color="#1565C0", alpha=0.1, zorder=0)

    ax.axhline(1, color="#2E7D32", linewidth=1, linestyle=":", alpha=0.5)
    ax.set_ylim(6.5, 0.5)
    ax.set_xlim(2, 30)
    ax.set_title(BAND_TEX[band], fontsize=11, fontweight="bold")
    if ib == 0:
        ax.set_ylabel("TL↔TT rank\n(1=best)")
    ax.set_xticks([5, 10, 15, 20, 25, 30])
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.2)

# --- Row 2: H2 per band (H2a contrast) ---
for ib, band in enumerate(BANDS):
    ax = fig.add_subplot(gs[1, ib])
    color = BAND_COLORS[band]
    curves = []
    for pat in PATIENTS:
        if band in h2a[pat]:
            c = h2a[pat][band]
            curves.append(c)
            ax.plot(K_ARR, c, color=color, alpha=0.25, linewidth=0.7,
                    marker=PAT_MARKERS[pat], markersize=2)
    if curves:
        ca = np.array(curves)
        mean_c = np.mean(ca, axis=0)
        ax.plot(K_ARR, mean_c, color=color, linewidth=2.5, zorder=5)
        n_pos = np.sum(ca > 0, axis=0)
        n_neg = np.sum(ca < 0, axis=0)
        ylim = max(0.3, np.abs(ca).max() * 1.1)
        ax.fill_between(K_ARR, -ylim, ylim, where=(n_pos >= 4),
                         color="#1565C0", alpha=0.15, zorder=0)
        ax.fill_between(K_ARR, -ylim, ylim, where=(n_neg >= 4),
                         color="#B71C1C", alpha=0.15, zorder=0)
        ax.fill_between(K_ARR, -ylim, ylim, where=(n_pos == 3),
                         color="#1565C0", alpha=0.06, zorder=0)
        ax.fill_between(K_ARR, -ylim, ylim, where=(n_neg == 3),
                         color="#B71C1C", alpha=0.06, zorder=0)
        ax.set_ylim(-ylim, ylim)

    ax.axhline(0, color="k", linewidth=0.8, linestyle="--")
    ax.set_xlim(2, 30)
    ax.set_title(BAND_TEX[band], fontsize=11, fontweight="bold")
    if ib == 0:
        ax.set_ylabel("H2a contrast\n(+persist/−recover)")
    ax.set_xticks([5, 10, 15, 20, 25, 30])
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.2)

# --- Row 3: H3 per band (within−cross gap) ---
for ib, band in enumerate(BANDS):
    ax = fig.add_subplot(gs[2, ib])
    gap_curves = []
    for pat in PATIENTS:
        wv = [vi_data[pat][band][p] for p in WITHIN_PAIRS if p in vi_data[pat][band]]
        cv = [vi_data[pat][band][p] for p in CROSS_PAIRS if p in vi_data[pat][band]]
        if wv and cv:
            gap = np.mean(cv, axis=0) - np.mean(wv, axis=0)
            gap_curves.append(gap)
            ax.plot(K_ARR, gap, color="#4CAF50", alpha=0.25, linewidth=0.7,
                    marker=PAT_MARKERS[pat], markersize=2)
    if gap_curves:
        ga = np.array(gap_curves)
        mean_g = np.mean(ga, axis=0)
        ax.plot(K_ARR, mean_g, color="#4CAF50", linewidth=2.5, zorder=5)
        all_pos = np.all(ga > 0, axis=0)
        ylim = max(0.2, np.abs(ga).max() * 1.1)
        ax.fill_between(K_ARR, -ylim, ylim, where=all_pos,
                         color="#4CAF50", alpha=0.1, zorder=0)
        ax.set_ylim(-ylim, ylim)

    ax.axhline(0, color="k", linewidth=0.8, linestyle="--")
    ax.set_xlim(2, 30)
    ax.set_title(BAND_TEX[band], fontsize=11, fontweight="bold")
    if ib == 0:
        ax.set_ylabel("Cross − Within\n(+H3 holds)")
    ax.set_xticks([5, 10, 15, 20, 25, 30])
    ax.set_xlabel("k")
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.2)

# --- Row 4: H4 — band ranking at each k (compact) ---
# Use a single wide panel
ax = fig.add_subplot(gs[3, :])
for band in BANDS:
    all_cross_mean = []
    for pat in PATIENTS:
        cv = [vi_data[pat][band][p] for p in CROSS_PAIRS if p in vi_data[pat][band]]
        if cv:
            all_cross_mean.append(np.mean(cv, axis=0))
    if all_cross_mean:
        mean_c = np.mean(all_cross_mean, axis=0)
        ax.plot(K_ARR, mean_c, color=BAND_COLORS[band], linewidth=2.5,
                label=BAND_TEX[band])
        lo = np.min(all_cross_mean, axis=0)
        hi = np.max(all_cross_mean, axis=0)
        ax.fill_between(K_ARR, lo, hi, color=BAND_COLORS[band], alpha=0.08)

ax.set_xlabel("k (number of communities)", fontsize=10)
ax.set_ylabel("Mean cross VI\n(reorganization strength)")
ax.set_title("H4: Reorganization strength per band", fontsize=11, fontweight="bold")
ax.legend(fontsize=8, ncol=6, loc="upper center")
ax.set_xlim(2, 30)
ax.grid(alpha=0.3)

# Row labels
fig.text(0.01, 0.88, "H1", fontsize=14, fontweight="bold", color="#1565C0",
         va="center", rotation=0)
fig.text(0.01, 0.64, "H2", fontsize=14, fontweight="bold", color="#E91E63",
         va="center", rotation=0)
fig.text(0.01, 0.40, "H3", fontsize=14, fontweight="bold", color="#4CAF50",
         va="center", rotation=0)
fig.text(0.01, 0.12, "H4", fontsize=14, fontweight="bold", color="#FF9800",
         va="center", rotation=0)

fig.savefig(OUTDIR / "all_hypotheses_multiscale.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved all_hypotheses_multiscale.pdf")


# ── Print summary ─────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("MULTISCALE HYPOTHESIS SUMMARY")
print("=" * 80)

print("\nH1 — TL↔TT rank=1 (most similar) for all 4 patients at k=...")
for band in BANDS:
    rank_curves = []
    for pat in PATIENTS:
        tl_tt = vi_data[pat][band].get(("task_learn", "task_test"))
        if tl_tt is None:
            continue
        ranks = []
        for ik in range(len(K_RANGE)):
            all_v = [vi_data[pat][band][p][ik] for p in ALL_PAIRS if p in vi_data[pat][band]]
            sorted_v = sorted(all_v)
            ranks.append(sorted_v.index(tl_tt[ik]) + 1)
        rank_curves.append(ranks)
    if rank_curves:
        all_r1 = np.all(np.array(rank_curves) == 1, axis=0)
        ks = [K_RANGE[i] for i in range(len(K_RANGE)) if all_r1[i]]
        print(f"  {BAND_TEX[band]:>12s}: {ks if ks else 'none'}")

print("\nH2a — unanimous regions:")
for band in BANDS:
    curves = [h2a[pat][band] for pat in PATIENTS if band in h2a[pat]]
    if not curves:
        continue
    ca = np.array(curves)
    pos_ks = [K_RANGE[i] for i in range(len(K_RANGE)) if np.all(ca[:, i] > 0)]
    neg_ks = [K_RANGE[i] for i in range(len(K_RANGE)) if np.all(ca[:, i] < 0)]
    parts = []
    if pos_ks:
        parts.append(f"PERSIST at k={pos_ks}")
    if neg_ks:
        parts.append(f"RECOVER at k={neg_ks}")
    print(f"  {BAND_TEX[band]:>12s}: {' | '.join(parts) if parts else 'none'}")

print("\nH3 — within<cross unanimous at k=...")
for band in BANDS:
    gaps = []
    for pat in PATIENTS:
        wv = [vi_data[pat][band][p] for p in WITHIN_PAIRS if p in vi_data[pat][band]]
        cv = [vi_data[pat][band][p] for p in CROSS_PAIRS if p in vi_data[pat][band]]
        if wv and cv:
            gaps.append(np.mean(cv, axis=0) - np.mean(wv, axis=0))
    if gaps:
        ga = np.array(gaps)
        all_pos = np.all(ga > 0, axis=0)
        ks = [K_RANGE[i] for i in range(len(K_RANGE)) if all_pos[i]]
        print(f"  {BAND_TEX[band]:>12s}: {len(ks)}/{len(K_RANGE)} k values ({ks[:5]}{'...' if len(ks) > 5 else ''})")

print("\nDone!")
