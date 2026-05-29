#!/usr/bin/env python3
"""WP-viz — VI(k) profile visualization candidates.

Options A–G for better visualization of multiscale VI profiles.

Run: python scripts/wp1/wp1_viz.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

# ── Config ────────────────────────────────────────────────────────────
CSV = ROOT / "data" / "wp0_metric_exploration" / "task3_multiscale" / "vi_raw_profiles.csv"
OUT = ROOT / "data" / "wp1_viz_investigation"
OUT.mkdir(parents=True, exist_ok=True)

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {b: BRAIN_BAND_TEX_DICT[b] for b in BANDS}
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
ALL_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"), ("rest_pre", "rest_post"),
    ("task_learn", "task_test"), ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
PAIR_SHORT = {
    ("rest_pre", "rest_post"): "Pre–Post",
    ("task_learn", "task_test"): "TL–TT",
    ("rest_pre", "task_learn"): "Pre–TL",
    ("rest_pre", "task_test"): "Pre–TT",
    ("task_learn", "rest_post"): "TL–Post",
    ("task_test", "rest_post"): "TT–Post",
}
PAIR_COLORS = {
    ("task_learn", "task_test"): "#FF9800",   # orange — H1 focal pair
    ("rest_pre", "rest_post"):      "#2196F3",    # blue — H2 focal pair
    ("rest_pre", "task_learn"):   "#E91E63",    # pink
    ("rest_pre", "task_test"):    "#9C27B0",    # purple
    ("task_learn", "rest_post"):  "#F44336",    # red
    ("task_test", "rest_post"):   "#795548",    # brown
}
PAIR_ZORDER = {("task_learn", "task_test"): 10, ("rest_pre", "rest_post"): 9}

PAT_COLORS = {
    "Pat_02": "#1f77b4", "Pat_03": "#ff7f0e", "Pat_05": "#2ca02c",
    "Pat_07": "#d62728", "Pat_08": "#9467bd",
}


def save_fig(fig, path, *, what, proves, how_to_read, dpi=300):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path}")
    md = path.with_suffix(".md")
    md.write_text(
        f"# {path.stem}\n\n"
        f"## What the figure shows\n{what}\n\n"
        f"## What it proves\n{proves}\n\n"
        f"## How to read it\n{how_to_read}\n"
    )


# ── Load data ─────────────────────────────────────────────────────────
print("Loading VI profiles...")
df = pd.read_csv(CSV)
df["pair_key"] = list(zip(df["phase_a"], df["phase_b"]))
print(f"  {len(df)} rows")

# Common k range across all patients (intersection)
k_common = None
for pat in PATIENTS:
    for band in BANDS:
        ks = set(df[(df["patient"] == pat) & (df["band"] == band)]["k"].unique())
        if k_common is None:
            k_common = ks
        else:
            k_common &= ks
k_range = sorted(k_common)
print(f"  Common k range: {min(k_range)}–{max(k_range)} ({len(k_range)} values)")


# ── Helper: build (band, k, pair) → array of patient VI values ────────
def build_vi_array():
    """Returns dict: (band, pair_key) → DataFrame with columns k + one col per patient."""
    result = {}
    for band in BANDS:
        for pair in ALL_PAIRS:
            sub = df[(df["band"] == band) & (df["pair_key"].apply(lambda x: x == pair))]
            piv = sub.pivot_table(index="k", columns="patient", values="vi")
            piv = piv.reindex(k_range).dropna(how="all")
            result[(band, pair)] = piv
    return result

vi_arr = build_vi_array()


# ======================================================================
# OPTION A: Patient-averaged profiles with ±1 std bands
# ======================================================================
def option_A():
    print("\nOption A: Averaged profiles...")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    axes = axes.flatten()

    for idx, band in enumerate(BANDS):
        ax = axes[idx]
        for pair in ALL_PAIRS:
            piv = vi_arr[(band, pair)]
            if piv.empty:
                continue
            mean = piv.mean(axis=1)
            std = piv.std(axis=1)
            ks = mean.index
            color = PAIR_COLORS[pair]
            zorder = PAIR_ZORDER.get(pair, 2)
            lw = 2.5 if pair in PAIR_ZORDER else 1.0
            alpha = 1.0 if pair in PAIR_ZORDER else 0.5
            ax.plot(ks, mean, color=color, lw=lw, zorder=zorder, alpha=alpha,
                    label=PAIR_SHORT[pair])
            ax.fill_between(ks, mean - std, mean + std, color=color,
                            alpha=0.12, zorder=zorder - 1)

        ax.set_title(BAND_TEX[band], fontsize=14)
        ax.set_ylabel("VI (nats)" if idx % 3 == 0 else "")
        if idx >= 3:
            ax.set_xlabel("k")
        ax.set_ylim(bottom=0)

    axes[0].legend(fontsize=7, loc="upper left", ncol=2)
    fig.suptitle("VI(k) profiles — patient-averaged (N=5) ± 1 std", fontsize=14)
    fig.tight_layout()
    save_fig(fig, OUT / "option_A_averaged_profiles.pdf",
        what="Mean VI(k) ± 1 std across 5 patients, per phase pair. One panel per band. "
             "Orange = TL–TT, Blue = Pre–Post. Other cross-type pairs in muted colors.",
        proves="H1: if orange sits below all others, task phases are most similar. "
               "H2: if blue is elevated in alpha, rest_post moved away from rest_pre.",
        how_to_read="Orange line below everything = H1. Blue overlapping with cross-type = H2. "
                    "Shaded bands show patient variability — narrow bands = consistent.",
    )

    # Normalized version
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    axes = axes.flatten()

    for idx, band in enumerate(BANDS):
        ax = axes[idx]
        for pair in ALL_PAIRS:
            piv = vi_arr[(band, pair)]
            if piv.empty:
                continue
            # Normalize each patient by their max VI across all pairs at each k
            norm_piv = piv.copy()
            for pat in norm_piv.columns:
                # Max VI for this patient/band across all pairs
                all_vi = pd.concat([vi_arr[(band, p)][pat] for p in ALL_PAIRS
                                    if pat in vi_arr[(band, p)].columns], axis=1)
                max_vi = all_vi.max(axis=1)
                norm_piv[pat] = norm_piv[pat] / max_vi.clip(lower=1e-10)

            mean = norm_piv.mean(axis=1)
            std = norm_piv.std(axis=1)
            ks = mean.index
            color = PAIR_COLORS[pair]
            zorder = PAIR_ZORDER.get(pair, 2)
            lw = 2.5 if pair in PAIR_ZORDER else 1.0
            alpha = 1.0 if pair in PAIR_ZORDER else 0.5
            ax.plot(ks, mean, color=color, lw=lw, zorder=zorder, alpha=alpha,
                    label=PAIR_SHORT[pair])
            ax.fill_between(ks, mean - std, mean + std, color=color,
                            alpha=0.12, zorder=zorder - 1)

        ax.set_title(BAND_TEX[band], fontsize=14)
        ax.set_ylabel("Normalized VI" if idx % 3 == 0 else "")
        if idx >= 3:
            ax.set_xlabel("k")
        ax.set_ylim(0, 1.1)

    axes[0].legend(fontsize=7, loc="upper left", ncol=2)
    fig.suptitle("VI(k) profiles — normalized per patient, averaged (N=5)", fontsize=14)
    fig.tight_layout()
    save_fig(fig, OUT / "option_A_normalized_profiles.pdf",
        what="Same as Option A but each patient's VI normalized by their max VI before averaging. "
             "Removes absolute scale differences between patients.",
        proves="Whether the RELATIVE ordering of pairs is consistent across patients, "
               "independent of each patient's overall VI scale.",
        how_to_read="Same as Option A. Values now in [0, 1]. If orange is consistently lowest "
                    "even after normalization, the signal is robust.",
    )


# ======================================================================
# OPTION B: Contrast profiles (H1-focused)
# ======================================================================
def option_B():
    print("\nOption B: Contrast profiles...")

    # B1: plain contrast
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    axes = axes.flatten()

    for idx, band in enumerate(BANDS):
        ax = axes[idx]
        piv_tt = vi_arr[(band, ("task_learn", "task_test"))]
        if piv_tt.empty:
            continue

        # Mean of other 5 pairs
        other_pivs = []
        for pair in ALL_PAIRS:
            if pair == ("task_learn", "task_test"):
                continue
            other_pivs.append(vi_arr[(band, pair)])

        contrasts_per_pat = {}
        for pat in PATIENTS:
            if pat not in piv_tt.columns:
                continue
            tt_vals = piv_tt[pat]
            other_vals = pd.concat([p[pat] for p in other_pivs if pat in p.columns], axis=1)
            other_mean = other_vals.mean(axis=1)
            contrast = tt_vals - other_mean
            contrasts_per_pat[pat] = contrast
            ax.plot(contrast.index, contrast.values, color=PAT_COLORS[pat],
                    lw=1.0, alpha=0.6, label=pat)

        # Patient average
        if contrasts_per_pat:
            all_c = pd.DataFrame(contrasts_per_pat)
            mean_c = all_c.mean(axis=1)
            ax.plot(mean_c.index, mean_c.values, color="black", lw=2.5,
                    zorder=10, label="Mean")

        ax.axhline(0, color="gray", ls=":", lw=0.8)
        ax.set_title(BAND_TEX[band], fontsize=14)
        ax.set_ylabel("Δ VI (TL-TT − others)" if idx % 3 == 0 else "")
        if idx >= 3:
            ax.set_xlabel("k")

    axes[0].legend(fontsize=7, loc="lower left", ncol=2)
    fig.suptitle("H1 contrast: VI(TL,TT) − mean(other pairs)", fontsize=14)
    fig.tight_layout()
    save_fig(fig, OUT / "option_B_contrast_profiles.pdf",
        what="H1 contrast per patient: VI(taskL,taskT) minus mean of other 5 pairs. "
             "One thin line per patient, thick black = mean. Dashed line at 0.",
        proves="H1: negative values mean TL-TT is most similar. All lines below 0 = unanimous.",
        how_to_read="Below zero = H1 supported for that patient at that k. "
                    "If ALL colored lines are below zero, H1 is unanimous. "
                    "Black line well below 0 = strong average effect.",
    )

    # B2: with unanimity bar
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    axes = axes.flatten()

    for idx, band in enumerate(BANDS):
        ax = axes[idx]
        piv_tt = vi_arr[(band, ("task_learn", "task_test"))]
        if piv_tt.empty:
            continue

        other_pivs = [vi_arr[(band, p)] for p in ALL_PAIRS
                      if p != ("task_learn", "task_test")]

        contrasts_per_pat = {}
        for pat in PATIENTS:
            if pat not in piv_tt.columns:
                continue
            tt_vals = piv_tt[pat]
            other_vals = pd.concat([p[pat] for p in other_pivs if pat in p.columns], axis=1)
            other_mean = other_vals.mean(axis=1)
            contrasts_per_pat[pat] = tt_vals - other_mean
            ax.plot((tt_vals - other_mean).index, (tt_vals - other_mean).values,
                    color=PAT_COLORS[pat], lw=1.0, alpha=0.6, label=pat)

        if contrasts_per_pat:
            all_c = pd.DataFrame(contrasts_per_pat)
            mean_c = all_c.mean(axis=1)
            ax.plot(mean_c.index, mean_c.values, color="black", lw=2.5, zorder=10)

            # Unanimity bar at bottom
            ymin = ax.get_ylim()[0]
            bar_h = (ax.get_ylim()[1] - ymin) * 0.04
            for k in all_c.index:
                row = all_c.loc[k].dropna()
                if len(row) == len(PATIENTS) and (row < 0).all():
                    ax.axvspan(k - 0.5, k + 0.5, ymin=0, ymax=0.04,
                               color="#2196F3", alpha=0.5, transform=ax.get_xaxis_transform(),
                               zorder=0)

        ax.axhline(0, color="gray", ls=":", lw=0.8)
        ax.set_title(BAND_TEX[band], fontsize=14)
        ax.set_ylabel("Δ VI" if idx % 3 == 0 else "")
        if idx >= 3:
            ax.set_xlabel("k")

    axes[0].legend(fontsize=7, loc="lower left", ncol=2)
    fig.suptitle("H1 contrast + unanimity bar (blue = 5/5 negative)", fontsize=14)
    fig.tight_layout()
    save_fig(fig, OUT / "option_B_contrast_with_unanimity.pdf",
        what="Same as Option B, with blue bar at bottom marking k values where ALL 5 patients "
             "have negative contrast (5/5 unanimity).",
        proves="Direct visual link between contrast curves and unanimity test.",
        how_to_read="Blue bar = unanimous H1 at that k. Width of blue region = robustness.",
    )


# ======================================================================
# OPTION C: Rank profiles
# ======================================================================
def option_C():
    print("\nOption C: Rank profiles...")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
    axes = axes.flatten()

    for idx, band in enumerate(BANDS):
        ax = axes[idx]

        ranks_per_pat = {}
        for pat in PATIENTS:
            rank_at_k = {}
            for k in k_range:
                vals = {}
                for pair in ALL_PAIRS:
                    piv = vi_arr[(band, pair)]
                    if pat in piv.columns and k in piv.index:
                        vals[pair] = piv.loc[k, pat]
                if len(vals) == 6:
                    sorted_pairs = sorted(vals.items(), key=lambda x: x[1])
                    for rank, (p, v) in enumerate(sorted_pairs, 1):
                        if p == ("task_learn", "task_test"):
                            rank_at_k[k] = rank
                            break
            if rank_at_k:
                ranks_per_pat[pat] = pd.Series(rank_at_k)
                ax.plot(list(rank_at_k.keys()), list(rank_at_k.values()),
                        color=PAT_COLORS[pat], lw=0.8, alpha=0.5, marker=".",
                        markersize=2)

        if ranks_per_pat:
            all_ranks = pd.DataFrame(ranks_per_pat)
            mean_rank = all_ranks.mean(axis=1)
            ax.plot(mean_rank.index, mean_rank.values, color="#FF9800", lw=3,
                    zorder=10, label="Mean rank (TL-TT)")

        ax.axhline(1, color="green", ls="--", lw=0.8, alpha=0.5, label="Rank 1 (best)")
        ax.axhline(3.5, color="red", ls=":", lw=0.8, alpha=0.5, label="Chance (3.5)")
        ax.set_ylim(0.5, 6.5)
        ax.invert_yaxis()
        ax.set_title(BAND_TEX[band], fontsize=14)
        ax.set_ylabel("Rank of TL-TT" if idx % 3 == 0 else "")
        if idx >= 3:
            ax.set_xlabel("k")

    axes[0].legend(fontsize=7, loc="lower left")
    fig.suptitle("Rank of TL-TT among 6 phase pairs (1=most similar)", fontsize=14)
    fig.tight_layout()
    save_fig(fig, OUT / "option_C_rank_profiles.pdf",
        what="Rank of taskL-taskT among 6 pairs at each k. Thin dots = individual patients. "
             "Thick orange = mean rank. Green dashed = rank 1. Red dotted = chance (3.5).",
        proves="H1: if orange stays near rank 1, task phases are consistently most similar.",
        how_to_read="Orange near top (rank 1) = H1 holds. Dipping toward 3.5 = signal lost. "
                    "Patient dots spread = inconsistency.",
    )


# ======================================================================
# OPTION D: Paired panels (profiles + unanimity bar)
# ======================================================================
def option_D():
    print("\nOption D: Paired panels...")
    fig, axes = plt.subplots(6, 2, figsize=(16, 20),
                              gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05})

    for idx, band in enumerate(BANDS):
        ax_prof = axes[idx, 0]
        ax_unan = axes[idx, 1]

        # Left: averaged profiles
        for pair in ALL_PAIRS:
            piv = vi_arr[(band, pair)]
            if piv.empty:
                continue
            mean = piv.mean(axis=1)
            std = piv.std(axis=1)
            color = PAIR_COLORS[pair]
            zorder = PAIR_ZORDER.get(pair, 2)
            lw = 2.5 if pair in PAIR_ZORDER else 1.0
            alpha = 1.0 if pair in PAIR_ZORDER else 0.5
            ax_prof.plot(mean.index, mean, color=color, lw=lw, zorder=zorder,
                         alpha=alpha, label=PAIR_SHORT[pair])
            ax_prof.fill_between(mean.index, mean - std, mean + std,
                                 color=color, alpha=0.1, zorder=zorder - 1)

        ax_prof.set_title(BAND_TEX[band], fontsize=14, loc="left")
        ax_prof.set_ylabel("VI (nats)")
        ax_prof.set_ylim(bottom=0)
        if idx == 0:
            ax_prof.legend(fontsize=6, loc="upper left", ncol=3)
        if idx < 5:
            ax_prof.set_xticklabels([])

        # Right: unanimity strip for H1
        piv_tt = vi_arr[(band, ("task_learn", "task_test"))]
        other_pivs = [vi_arr[(band, p)] for p in ALL_PAIRS
                      if p != ("task_learn", "task_test")]

        unan_map = np.full(len(k_range), np.nan)
        for ki, k in enumerate(k_range):
            signs = []
            for pat in PATIENTS:
                if pat not in piv_tt.columns or k not in piv_tt.index:
                    continue
                tt_v = piv_tt.loc[k, pat] if k in piv_tt.index else np.nan
                others = [p.loc[k, pat] for p in other_pivs
                          if pat in p.columns and k in p.index]
                if np.isfinite(tt_v) and len(others) == 5:
                    contrast = tt_v - np.mean(others)
                    signs.append(np.sign(contrast))
            if len(signs) == len(PATIENTS):
                unan_map[ki] = np.mean(signs)

        ax_unan.imshow(unan_map.reshape(-1, 1), aspect="auto", cmap="RdBu",
                        vmin=-1, vmax=1, interpolation="nearest",
                        extent=[0, 1, k_range[-1], k_range[0]])
        ax_unan.set_xticks([])
        ax_unan.set_ylabel("")
        ax_unan.yaxis.tick_right()
        if idx == 0:
            ax_unan.set_title("H1\nunan.", fontsize=9)
        if idx < 5:
            ax_unan.set_xticklabels([])

    axes[-1, 0].set_xlabel("k (number of clusters)")
    fig.suptitle("VI(k) profiles + H1 unanimity strips", fontsize=14, y=1.01)
    fig.tight_layout()
    save_fig(fig, OUT / "option_D_paired_panels.pdf",
        what="Left column: patient-averaged VI(k) profiles (same as Option A). "
             "Right column: H1 unanimity strip — blue = 5/5 patients support H1, "
             "red = 5/5 disagree, white = split. One row per band.",
        proves="Direct visual link between VI profile shapes and where H1 passes.",
        how_to_read="Read left to right: where orange drops below others (left), "
                    "the strip should be blue (right). Width of blue region = robustness.",
    )


# ======================================================================
# OPTION E: Heatmap of mean contrast
# ======================================================================
def option_E():
    print("\nOption E: Mean contrast heatmap...")
    fig, axes = plt.subplots(1, 4, figsize=(20, 4),
                              gridspec_kw={"width_ratios": [1, 1, 1, 0.04]})

    hyp_defs = [
        ("H1: TL-TT − others", lambda pdf, k: _h1_contrast(pdf, k)),
        ("H2a: Pre-Post − TT-Post", lambda pdf, k: _h2a_contrast(pdf, k)),
        ("H3: cross − within", lambda pdf, k: _h3_contrast(pdf, k)),
    ]

    for hi, (title, cfn) in enumerate(hyp_defs):
        ax = axes[hi]
        mat = np.full((len(BANDS), len(k_range)), np.nan)
        unan_mask = np.zeros_like(mat, dtype=bool)

        for bi, band in enumerate(BANDS):
            for ki, k in enumerate(k_range):
                vals = []
                for pat in PATIENTS:
                    c = cfn(df[(df["patient"] == pat) & (df["band"] == band) & (df["k"] == k)], k)
                    if c is not None and np.isfinite(c):
                        vals.append(c)
                if len(vals) == len(PATIENTS):
                    mat[bi, ki] = np.mean(vals)
                    if all(v < 0 for v in vals) or all(v > 0 for v in vals):
                        unan_mask[bi, ki] = True

        vmax = np.nanmax(np.abs(mat)) if np.any(np.isfinite(mat)) else 1
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                        interpolation="nearest",
                        extent=[k_range[0], k_range[-1], len(BANDS) - 0.5, -0.5])

        # Black borders on unanimous cells
        for bi in range(len(BANDS)):
            for ki in range(len(k_range)):
                if unan_mask[bi, ki]:
                    k_val = k_range[ki]
                    ax.add_patch(plt.Rectangle(
                        (k_val - 0.5, bi - 0.5), 1, 1,
                        fill=False, edgecolor="black", lw=0.5, zorder=5))

        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels([BAND_TEX[b] for b in BANDS], fontsize=10)
        ax.set_xlabel("k")
        ax.set_title(title, fontsize=11)

    cax = axes[3]
    fig.colorbar(im, cax=cax, label="Mean contrast (N=5)")

    fig.suptitle("Mean hypothesis contrasts (black border = 5/5 unanimous)", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT / "option_E_mean_contrast_heatmap.pdf",
        what="Three heatmaps (H1, H2a, H3): band × k, color = patient-averaged contrast. "
             "Blue = hypothesis supported (negative contrast), Red = opposed. "
             "Black borders = 5/5 unanimous.",
        proves="Continuous version of unanimity maps: shows effect STRENGTH not just sign.",
        how_to_read="Blue with black borders = strong unanimous support. "
                    "Compare across panels to see which bands support which hypotheses.",
    )


def _h1_contrast(pdf, k):
    if pdf.empty:
        return None
    pair_vals = {}
    for _, r in pdf.iterrows():
        pair_vals[(r["phase_a"], r["phase_b"])] = r["vi"]
    tt = pair_vals.get(("task_learn", "task_test"))
    others = [v for p, v in pair_vals.items() if p != ("task_learn", "task_test")]
    if tt is None or len(others) < 5:
        return None
    return tt - np.mean(others)


def _h2a_contrast(pdf, k):
    if pdf.empty:
        return None
    pair_vals = {}
    for _, r in pdf.iterrows():
        pair_vals[(r["phase_a"], r["phase_b"])] = r["vi"]
    pre_post = pair_vals.get(("rest_pre", "rest_post"))
    tt_post = pair_vals.get(("task_test", "rest_post"))
    if pre_post is None or tt_post is None:
        return None
    return pre_post - tt_post  # positive = trace


def _h3_contrast(pdf, k):
    if pdf.empty:
        return None
    pair_vals = {}
    for _, r in pdf.iterrows():
        pair_vals[(r["phase_a"], r["phase_b"])] = r["vi"]
    within = [pair_vals.get(p) for p in [("rest_pre", "rest_post"), ("task_learn", "task_test")]]
    cross = [pair_vals.get(p) for p in [("rest_pre", "task_learn"), ("rest_pre", "task_test"),
                                         ("task_learn", "rest_post"), ("task_test", "rest_post")]]
    within = [v for v in within if v is not None]
    cross = [v for v in cross if v is not None]
    if len(within) < 2 or len(cross) < 4:
        return None
    return np.mean(within) - np.mean(cross)  # negative = H3 supported


# ======================================================================
# OPTION F: Three key bands only (larger panels)
# ======================================================================
def option_F():
    print("\nOption F: Three key bands...")
    focus_bands = ["beta", "alpha", "theta"]
    focus_labels = ["$\\beta$ (strongest H1)", "$\\alpha$ (H2 trace)",
                    "$\\theta$ (weakest)"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    for idx, (band, label) in enumerate(zip(focus_bands, focus_labels)):
        ax = axes[idx]
        for pair in ALL_PAIRS:
            piv = vi_arr[(band, pair)]
            if piv.empty:
                continue
            mean = piv.mean(axis=1)
            std = piv.std(axis=1)
            color = PAIR_COLORS[pair]
            zorder = PAIR_ZORDER.get(pair, 2)
            lw = 2.5 if pair in PAIR_ZORDER else 1.0
            alpha = 1.0 if pair in PAIR_ZORDER else 0.5
            ax.plot(mean.index, mean, color=color, lw=lw, zorder=zorder,
                    alpha=alpha, label=PAIR_SHORT[pair])
            ax.fill_between(mean.index, mean - std, mean + std,
                            color=color, alpha=0.1, zorder=zorder - 1)

        ax.set_title(label, fontsize=13)
        ax.set_xlabel("k")
        ax.set_ylim(bottom=0)
        if idx == 0:
            ax.set_ylabel("VI (nats)")

    axes[-1].legend(fontsize=8, loc="upper right")
    fig.suptitle("VI(k) — three contrasting bands (averaged, N=5)", fontsize=14)
    fig.tight_layout()
    save_fig(fig, OUT / "option_F_three_bands.pdf",
        what="Averaged VI(k) profiles for beta (strongest H1), alpha (H2 trace), "
             "theta (weakest). Larger panels for readability. Orange = TL-TT, Blue = Pre-Post.",
        proves="The three scientific points in one figure: (1) beta orange clearly lowest, "
               "(2) alpha blue elevated, (3) theta all curves overlap.",
        how_to_read="Beta: orange clearly separated below → H1. "
                    "Alpha: blue rises above cross-pairs → H2. "
                    "Theta: no clear separation → no signal.",
    )


# ======================================================================
# OPTION G: Contrast + rank combo (custom)
# ======================================================================
def option_G():
    """Combined figure: top row = H1 contrast per band, bottom = H2a contrast for alpha only."""
    print("\nOption G: Custom combo...")

    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 6, hspace=0.35, wspace=0.3,
                          height_ratios=[1, 1])

    # Top row: H1 contrast for all 6 bands
    for idx, band in enumerate(BANDS):
        ax = fig.add_subplot(gs[0, idx])
        piv_tt = vi_arr[(band, ("task_learn", "task_test"))]
        if piv_tt.empty:
            continue
        other_pivs = [vi_arr[(band, p)] for p in ALL_PAIRS
                      if p != ("task_learn", "task_test")]

        for pat in PATIENTS:
            if pat not in piv_tt.columns:
                continue
            tt_vals = piv_tt[pat]
            others = pd.concat([p[pat] for p in other_pivs if pat in p.columns], axis=1)
            contrast = tt_vals - others.mean(axis=1)
            ax.plot(contrast.index, contrast.values, color=PAT_COLORS[pat],
                    lw=0.8, alpha=0.5)

        # Mean + unanimity fill
        all_contrasts = {}
        for pat in PATIENTS:
            if pat not in piv_tt.columns:
                continue
            tt_vals = piv_tt[pat]
            others = pd.concat([vi_arr[(band, p)][pat] for p in ALL_PAIRS
                                if p != ("task_learn", "task_test") and
                                pat in vi_arr[(band, p)].columns], axis=1)
            all_contrasts[pat] = tt_vals - others.mean(axis=1)

        if all_contrasts:
            cdf = pd.DataFrame(all_contrasts)
            mean_c = cdf.mean(axis=1)
            ax.plot(mean_c.index, mean_c.values, "k-", lw=2, zorder=10)

            # Fill where unanimous
            for k in cdf.index:
                row = cdf.loc[k].dropna()
                if len(row) == len(PATIENTS) and (row < 0).all():
                    ax.axvspan(k - 0.5, k + 0.5, color="#2196F3", alpha=0.15, zorder=0)

        ax.axhline(0, color="gray", ls=":", lw=0.8)
        ax.set_title(BAND_TEX[band], fontsize=11)
        if idx == 0:
            ax.set_ylabel("H1 contrast\n(TL-TT − others)")

    # Bottom row: H2a contrast for alpha (big panel) + Pre-Post trace comparison
    ax_h2 = fig.add_subplot(gs[1, :3])
    band = "alpha"
    for pat in PATIENTS:
        contrasts = {}
        for k in k_range:
            sub = df[(df["patient"] == pat) & (df["band"] == band) & (df["k"] == k)]
            c = _h2a_contrast(sub, k)
            if c is not None:
                contrasts[k] = c
        if contrasts:
            ks, cs = zip(*sorted(contrasts.items()))
            ax_h2.plot(ks, cs, color=PAT_COLORS[pat], lw=1.2, alpha=0.7, label=pat)

    # Mean
    mean_h2 = {}
    for k in k_range:
        vals = []
        for pat in PATIENTS:
            sub = df[(df["patient"] == pat) & (df["band"] == band) & (df["k"] == k)]
            c = _h2a_contrast(sub, k)
            if c is not None:
                vals.append(c)
        if len(vals) == len(PATIENTS):
            mean_h2[k] = np.mean(vals)
            if all(v > 0 for v in vals):
                ax_h2.axvspan(k - 0.5, k + 0.5, color="#FF9800", alpha=0.15, zorder=0)

    if mean_h2:
        ks, ms = zip(*sorted(mean_h2.items()))
        ax_h2.plot(ks, ms, "k-", lw=2.5, zorder=10, label="Mean")

    ax_h2.axhline(0, color="gray", ls=":", lw=0.8)
    ax_h2.set_xlabel("k")
    ax_h2.set_ylabel("H2a contrast\n(Pre-Post − TT-Post)")
    ax_h2.set_title(f"H2a task trace — {BAND_TEX[band]} band", fontsize=12)
    ax_h2.legend(fontsize=7, loc="upper right", ncol=3)

    # Bottom right: same for beta (anti-trace)
    ax_h2b = fig.add_subplot(gs[1, 3:])
    band = "beta"
    for pat in PATIENTS:
        contrasts = {}
        for k in k_range:
            sub = df[(df["patient"] == pat) & (df["band"] == band) & (df["k"] == k)]
            c = _h2a_contrast(sub, k)
            if c is not None:
                contrasts[k] = c
        if contrasts:
            ks, cs = zip(*sorted(contrasts.items()))
            ax_h2b.plot(ks, cs, color=PAT_COLORS[pat], lw=1.2, alpha=0.7, label=pat)

    mean_h2b = {}
    for k in k_range:
        vals = []
        for pat in PATIENTS:
            sub = df[(df["patient"] == pat) & (df["band"] == band) & (df["k"] == k)]
            c = _h2a_contrast(sub, k)
            if c is not None:
                vals.append(c)
        if len(vals) == len(PATIENTS):
            mean_h2b[k] = np.mean(vals)

    if mean_h2b:
        ks, ms = zip(*sorted(mean_h2b.items()))
        ax_h2b.plot(ks, ms, "k-", lw=2.5, zorder=10, label="Mean")

    ax_h2b.axhline(0, color="gray", ls=":", lw=0.8)
    ax_h2b.set_xlabel("k")
    ax_h2b.set_ylabel("H2a contrast")
    ax_h2b.set_title(f"H2a task trace — {BAND_TEX['beta']} band (anti-trace)", fontsize=12)
    ax_h2b.legend(fontsize=7, loc="upper right", ncol=3)

    fig.suptitle("H1 (top) + H2a alpha vs beta (bottom)", fontsize=14, y=1.01)
    fig.tight_layout()
    save_fig(fig, OUT / "option_G_custom.pdf",
        what="Top: H1 contrast per band (6 panels), with blue shading = unanimous. "
             "Bottom-left: H2a contrast for alpha (orange shading = unanimous trace). "
             "Bottom-right: H2a contrast for beta (no trace / anti-trace). "
             "Each thin line = one patient, thick black = mean.",
        proves="Combines H1 (top) and H2 (bottom) in a single figure. "
               "The alpha-beta dissociation is directly visible in the bottom row.",
        how_to_read="Top: blue regions = H1 unanimous. Beta has the widest. "
                    "Bottom-left: orange regions = alpha trace unanimous. "
                    "Bottom-right: beta stays near/above zero = no trace.",
    )


# ======================================================================
# Recommendation
# ======================================================================
def write_recommendation():
    print("\nWriting recommendation...")
    md = [
        "# Visualization Recommendation",
        "",
        "## Main figure for paper: Option G (custom combo)",
        "",
        "Option G combines H1 and H2 in a single figure with direct contrast",
        "visualization. The top row shows H1 across all bands (6 small panels",
        "with unanimity shading). The bottom row shows the alpha-beta H2",
        "dissociation as a paired comparison. This is the most information-dense",
        "option that still communicates all three scientific points:",
        "1. H1 holds broadly (top row, blue shading in all bands)",
        "2. H2 is alpha-specific (bottom-left, orange shading)",
        "3. Beta is the anti-trace (bottom-right, no shading)",
        "",
        "## Strong alternatives",
        "",
        "- **Option B (contrast + unanimity)**: best for H1 alone. The unanimity",
        "  bar makes the 5/5 criterion visually explicit. Use as supplementary.",
        "- **Option E (heatmap)**: most compact summary — all hypotheses in one",
        "  figure. Good for supplementary or if space is tight.",
        "- **Option F (three bands)**: clearest for a general audience but loses",
        "  the full cross-band comparison. Good for talks/posters.",
        "",
        "## Supplementary material",
        "",
        "- Option A (averaged profiles): complete reference for all bands/pairs",
        "- Option C (rank profiles): alternative perspective on H1 robustness",
        "- Option D (paired panels): links profiles to unanimity directly",
        "",
        "## What failed",
        "",
        "- **Option A normalized**: the normalization compresses the scale but",
        "  doesn't add information beyond the raw version. The ±std bands are",
        "  actually more informative in absolute units.",
        "- **Option C alone**: the rank view loses magnitude information. Useful",
        "  as a supplement but not as a main figure.",
    ]
    (OUT / "recommendation.md").write_text("\n".join(md))
    print(f"  Saved: {OUT / 'recommendation.md'}")


# ── Main ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    option_A()
    option_B()
    option_C()
    option_D()
    option_E()
    option_F()
    option_G()
    write_recommendation()
    print(f"\nAll done. Outputs in {OUT}")
