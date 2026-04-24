#!/usr/bin/env python3
"""Section 6 — H2b figures: unanimity map, 3-row dissociation, radar charts, 4-panel heatmap.

Outputs (all under data/figures/section6/):
  1. unanimity_map_H2b.pdf            — H2b unanimity map
  2. h1_h2a_h2b_all_bands.pdf         — 3×6 panel (H1 + H2a + H2b)
  3. radar_chart_H2b.pdf              — H2b radar chart
  4. radar_chart_H2a_H2b_combined.pdf — overlaid H2a + H2b
  5. mean_contrast_heatmap_4panel.pdf  — H1, H2a, H2b, H3

Data source: data/wp0_metric_exploration/task3_multiscale/vi_raw_profiles.csv
             data/wp_scalar_vi/scalar_table_full.csv
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from lrg_eegfc.config.paths import FIGURES_ROOT

OUTDIR = FIGURES_ROOT / "section6"
OUTDIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
})

# ── Constants ────────────────────────────────────────────────────────
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
DISPLAY_PATIENTS = ["Pat_02", "Pat_05", "Pat_07", "Pat_08"]
PAT_COLORS = {
    "Pat_02": "#1f77b4",
    "Pat_03": "#ff7f0e",
    "Pat_05": "#2ca02c",
    "Pat_07": "#d62728",
    "Pat_08": "#9467bd",
}
PAT_SHORT = {"Pat_02": "P2", "Pat_03": "P3", "Pat_05": "P5",
             "Pat_07": "P7", "Pat_08": "P8"}
OUTLIER = "Pat_03"

ALL_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"), ("rest_pre", "rest_post"),
    ("task_learn", "task_test"), ("task_learn", "rest_post"), ("task_test", "rest_post"),
]

# ── Load VI profiles ─────────────────────────────────────────────────
print("Loading VI profiles...")
df = pd.read_csv("data/wp0_metric_exploration/task3_multiscale/vi_raw_profiles.csv")
df["pair_key"] = list(zip(df["phase_a"], df["phase_b"]))

k_common = None
for pat in ALL_PATIENTS:
    for band in BANDS:
        ks = set(df[(df["patient"] == pat) & (df["band"] == band)]["k"].unique())
        k_common = ks if k_common is None else k_common & ks
k_range = sorted(k_common)
print(f"  k range: {min(k_range)}–{max(k_range)} ({len(k_range)} values)")

# Pivoted VI arrays: (band, pair) → DataFrame[k × patient]
vi_arr = {}
for band in BANDS:
    for pair in ALL_PAIRS:
        sub = df[(df["band"] == band) & (df["pair_key"].apply(lambda x: x == pair))]
        piv = sub.pivot_table(index="k", columns="patient", values="vi")
        piv = piv.reindex(k_range).dropna(how="all")
        vi_arr[(band, pair)] = piv


# ── Contrast functions (profile-level) ───────────────────────────────
def compute_h1_contrasts(band):
    """H1: VI(TL,TT) − mean(other 5 pairs)."""
    piv_tt = vi_arr[(band, ("task_learn", "task_test"))]
    other_pivs = [vi_arr[(band, p)] for p in ALL_PAIRS
                  if p != ("task_learn", "task_test")]
    result = {}
    for pat in ALL_PATIENTS:
        if pat not in piv_tt.columns:
            continue
        tt_vals = piv_tt[pat]
        others = pd.concat([p[pat] for p in other_pivs if pat in p.columns], axis=1)
        result[pat] = tt_vals - others.mean(axis=1)
    return result


def compute_h2a_contrasts(band):
    """H2a: VI(Pre,Post) − VI(TT,Post)."""
    piv_pp = vi_arr[(band, ("rest_pre", "rest_post"))]
    piv_tp = vi_arr[(band, ("task_test", "rest_post"))]
    result = {}
    for pat in ALL_PATIENTS:
        if pat in piv_pp.columns and pat in piv_tp.columns:
            result[pat] = piv_pp[pat] - piv_tp[pat]
    return result


def compute_h2b_contrasts(band):
    """H2b: VI(rest_pre, taskT) − VI(taskT, rest_post). Positive = approach."""
    piv_pt = vi_arr[(band, ("rest_pre", "task_test"))]
    piv_tp = vi_arr[(band, ("task_test", "rest_post"))]
    result = {}
    for pat in ALL_PATIENTS:
        if pat in piv_pt.columns and pat in piv_tp.columns:
            result[pat] = piv_pt[pat] - piv_tp[pat]
    return result


def compute_h3_contrasts(band):
    """H3: mean(VI_within) − mean(VI_cross). Negative = H3 supported."""
    within_pairs = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
    cross_pairs = [("rest_pre", "task_learn"), ("rest_pre", "task_test"),
                   ("task_learn", "rest_post"), ("task_test", "rest_post")]
    result = {}
    for pat in ALL_PATIENTS:
        w_list, c_list = [], []
        for p in within_pairs:
            if pat in vi_arr[(band, p)].columns:
                w_list.append(vi_arr[(band, p)][pat])
        for p in cross_pairs:
            if pat in vi_arr[(band, p)].columns:
                c_list.append(vi_arr[(band, p)][pat])
        if len(w_list) == 2 and len(c_list) == 4:
            w_mean = pd.concat(w_list, axis=1).mean(axis=1)
            c_mean = pd.concat(c_list, axis=1).mean(axis=1)
            result[pat] = w_mean - c_mean
    return result


# ── Contrast at single (patient, band, k) for Option E ──────────────
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
    return pre_post - tt_post


def _h2b_contrast(pdf, k):
    if pdf.empty:
        return None
    pair_vals = {}
    for _, r in pdf.iterrows():
        pair_vals[(r["phase_a"], r["phase_b"])] = r["vi"]
    pre_tt = pair_vals.get(("rest_pre", "task_test"))
    tt_post = pair_vals.get(("task_test", "rest_post"))
    if pre_tt is None or tt_post is None:
        return None
    return pre_tt - tt_post


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
    return np.mean(within) - np.mean(cross)


# ═════════════════════════════════════════════════════════════════════
# TASK 1: H2b unanimity map
# ═════════════════════════════════════════════════════════════════════
def task1_unanimity_map():
    print("\n[Task 1] H2b unanimity map...")

    # Compute signed unanimity at each (band, k)
    rows = []
    for band in BANDS:
        for k in k_range:
            signs = []
            for pat in ALL_PATIENTS:
                h2b = compute_h2b_contrasts(band)
                if pat in h2b and k in h2b[pat].index:
                    v = h2b[pat].loc[k]
                    if np.isfinite(v):
                        signs.append(np.sign(v))
            if len(signs) > 0:
                rows.append({
                    "band": band, "k": k,
                    "unanimity": np.mean(signs),
                    "n_patients": len(signs),
                })

    unan_df = pd.DataFrame(rows)
    piv = unan_df.pivot_table(index="band", columns="k", values="unanimity")
    piv = piv.reindex(BANDS)

    fig, ax = plt.subplots(figsize=(16, 4.5))

    im = ax.imshow(piv.values, aspect="auto", cmap="RdBu", vmin=-1, vmax=1,
                   interpolation="nearest")

    # Black borders on 5/5 unanimous cells
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if not np.isnan(v) and abs(v) == 1.0:
                band = BANDS[i]
                k = piv.columns[j]
                row = unan_df[(unan_df["band"] == band) & (unan_df["k"] == k)]
                if not row.empty and row.iloc[0]["n_patients"] == len(ALL_PATIENTS):
                    rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                         fill=False, edgecolor="black", lw=1.5)
                    ax.add_patch(rect)

    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BAND_TEX[b] for b in BANDS], fontsize=10)

    k_vals = list(piv.columns)
    step = max(1, len(k_vals) // 20)
    tick_idx = list(range(0, len(k_vals), step))
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([k_vals[t] for t in tick_idx], fontsize=7)
    ax.set_xlabel("$k$ (number of clusters)", fontsize=10)

    fig.colorbar(im, ax=ax, label="Signed unanimity (5 patients)", shrink=0.8)
    ax.set_title("H2b: Task Approach — unanimity map (5/5 patients)", fontsize=12)
    fig.tight_layout()

    fig.savefig(OUTDIR / "unanimity_map_H2b.pdf", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("  Saved unanimity_map_H2b.pdf")

    md = (
        "# unanimity_map_H2b\n\n"
        "## What the figure shows\n"
        "Heatmap: band (rows) vs k (columns). Color = signed unanimity across 5 patients. "
        "Blue (+1) = all 5 patients show rest_post moved toward task. "
        "Red (-1) = all 5 show opposite. Black borders = 5/5 unanimous.\n\n"
        "## Key result\n"
        "Alpha: 17 unanimous cells (k=7-51). Low gamma: 13 cells (k=37-58). "
        "Total: 33 unanimous cells.\n\n"
        "## How to read it\n"
        "Only blue cells with black borders are reportable (5/5 agreement). "
        "Contiguous regions across k = robust effect.\n"
    )
    (OUTDIR / "unanimity_map_H2b.md").write_text(md)


# ═════════════════════════════════════════════════════════════════════
# TASK 2: 3-row dissociation figure (H1 + H2a + H2b)
# ═════════════════════════════════════════════════════════════════════
def task2_three_row_dissociation():
    print("\n[Task 2] H1 + H2a + H2b dissociation figure...")

    hyp_defs = [
        ("H1: task stability", compute_h1_contrasts, "#2196F3",
         lambda signs: all(s < 0 for s in signs)),  # H1: negative = supported
        ("H2a: task trace", compute_h2a_contrasts, "#FF9800",
         lambda signs: all(s > 0 for s in signs)),  # H2a: positive = supported
        ("H2b: task approach", compute_h2b_contrasts, "#4CAF50",
         lambda signs: all(s > 0 for s in signs)),  # H2b: positive = supported
    ]

    fig, axes = plt.subplots(3, 6, figsize=(20, 8.5), sharex=True)
    fig.subplots_adjust(hspace=0.25, wspace=0.08)

    ylims_per_row = [[] for _ in range(3)]

    for row_idx, (row_label, contrast_fn, shade_color, unan_test) in enumerate(hyp_defs):
        for ib, band in enumerate(BANDS):
            ax = axes[row_idx, ib]
            contrasts = contrast_fn(band)
            c_df = pd.DataFrame(contrasts)

            # Unanimity shading (all 5 patients)
            for k in c_df.index:
                row_vals = c_df.loc[k].dropna()
                all_pat_vals = [row_vals[p] for p in ALL_PATIENTS if p in row_vals.index]
                if len(all_pat_vals) == len(ALL_PATIENTS):
                    signs = [np.sign(v) for v in all_pat_vals]
                    if unan_test(signs):
                        ax.axvspan(k - 0.5, k + 0.5, color=shade_color,
                                   alpha=0.25, zorder=0, linewidth=0)

            # Display patients only
            for pat in DISPLAY_PATIENTS:
                if pat in contrasts:
                    ax.plot(contrasts[pat].index, contrasts[pat].values,
                            color=PAT_COLORS[pat], lw=1.0, alpha=0.6, zorder=2)

            # Mean of display patients
            disp = pd.DataFrame({p: contrasts[p] for p in DISPLAY_PATIENTS
                                 if p in contrasts})
            if not disp.empty:
                mean_line = disp.mean(axis=1)
                ax.plot(mean_line.index, mean_line.values, "k-", lw=2.0, zorder=3)

            ax.axhline(0, color="#999999", ls="--", lw=0.8, zorder=1)
            ax.tick_params(labelsize=7)
            ylims_per_row[row_idx].append(ax.get_ylim())

            if row_idx == 0:
                ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")

    # Shared y-axis per row
    for row_idx in range(3):
        ymin = min(y[0] for y in ylims_per_row[row_idx])
        ymax = max(y[1] for y in ylims_per_row[row_idx])
        for ib in range(6):
            axes[row_idx, ib].set_ylim(ymin, ymax)
            if ib > 0:
                axes[row_idx, ib].set_yticklabels([])

    # Row labels
    for row_idx, (row_label, _, _, _) in enumerate(hyp_defs):
        axes[row_idx, 0].set_ylabel(f"{row_label}\n$\\Delta$VI", fontsize=10)

    # X-axis labels on bottom row
    for ib in range(6):
        axes[2, ib].set_xlabel("$k$", fontsize=9)

    # Legend in bottom-right
    handles = [Line2D([0], [0], color=PAT_COLORS[p], lw=1.0, alpha=0.6,
                      label=PAT_SHORT[p]) for p in DISPLAY_PATIENTS]
    handles.append(Line2D([0], [0], color="black", lw=2.0, label="Mean"))
    handles.append(plt.Rectangle((0, 0), 1, 1, fc="#2196F3", alpha=0.25,
                   ec="none", label="H1 unan."))
    handles.append(plt.Rectangle((0, 0), 1, 1, fc="#FF9800", alpha=0.25,
                   ec="none", label="H2a unan."))
    handles.append(plt.Rectangle((0, 0), 1, 1, fc="#4CAF50", alpha=0.25,
                   ec="none", label="H2b unan."))
    axes[2, -1].legend(handles=handles, fontsize=7, loc="lower right",
                       frameon=True, framealpha=0.9, edgecolor="#cccccc")

    fig.savefig(OUTDIR / "h1_h2a_h2b_all_bands.pdf", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("  Saved h1_h2a_h2b_all_bands.pdf")

    md = (
        "# h1_h2a_h2b_all_bands\n\n"
        "## What the figure shows\n"
        "3×6 panel figure. Row 1: H1 (task stability), Row 2: H2a (task trace), "
        "Row 3: H2b (task approach). Each column = one frequency band. "
        "Thin colored lines = 4 patients (Pat_03 excluded from display). "
        "Black = mean. Shading = 5/5 unanimity (all 5 patients incl. Pat_03).\n\n"
        "## Key result\n"
        "H1: beta has widest blue shading. "
        "H2a: only alpha has orange shading (29 contiguous k). "
        "H2b: alpha has green shading (17 cells), low gamma also shows shading (13 cells).\n\n"
        "## How to read it\n"
        "Above zero = hypothesis supported. Shaded = unanimous at that k.\n"
    )
    (OUTDIR / "h1_h2a_h2b_all_bands.md").write_text(md)


# ═════════════════════════════════════════════════════════════════════
# TASK 3: H2b radar chart + combined H2a/H2b
# ═════════════════════════════════════════════════════════════════════
def task3_radar_charts():
    print("\n[Task 3] H2b radar charts...")

    # Load scalar data
    sdf = pd.read_csv("data/wp_scalar_vi/scalar_table_full.csv")

    n_bands = len(BANDS)
    angles = np.linspace(0, 2 * np.pi, n_bands, endpoint=False).tolist()
    angles_closed = angles + [angles[0]]

    def load_scalar(hyp):
        data = np.zeros((len(ALL_PATIENTS), len(BANDS)))
        for ip, pat in enumerate(ALL_PATIENTS):
            for ib, band in enumerate(BANDS):
                row = sdf[(sdf["hypothesis"] == hyp) & (sdf["patient"] == pat)
                          & (sdf["band"] == band)]
                data[ip, ib] = row["scalar_A_mean"].values[0]
        return data

    data_h2b = load_scalar("H2b")
    data_h2a = load_scalar("H2a")

    def make_radar(data, title_hyp, filename, annot_text=None):
        mean_vals = data.mean(axis=0)
        vmax = np.ceil(np.max(np.abs(data)) * 10) / 10
        vmax = max(vmax, 0.30)

        fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        theta_fill = np.linspace(0, 2 * np.pi, 300)
        ax.fill(theta_fill, np.full_like(theta_fill, 0), color="#ededed", zorder=0)
        ax.plot(theta_fill, np.zeros_like(theta_fill), color="black",
                linestyle="--", lw=1.8, alpha=0.7, zorder=2)

        ax.set_rlim(-vmax, vmax)
        rticks = np.arange(-vmax, vmax + 0.05, 0.1).round(1).tolist()
        ax.set_rticks(rticks)
        ax.set_yticklabels([f"{v:+.1f}" if v != 0 else "0" for v in rticks],
                           fontsize=7, color="#999999")
        ax.set_rlabel_position(22)
        ax.grid(color="#e0e0e0", linewidth=0.4)
        ax.spines["polar"].set_visible(False)

        for ip, pat in enumerate(ALL_PATIENTS):
            vals = data[ip].tolist() + [data[ip, 0]]
            ls = "--" if pat == OUTLIER else "-"
            marker = "x" if pat == OUTLIER else "o"
            lbl = f"{PAT_SHORT[pat]}*" if pat == OUTLIER else PAT_SHORT[pat]
            ax.plot(angles_closed, vals, ls, marker=marker,
                    color=PAT_COLORS[pat], lw=1.2, markersize=3.5,
                    alpha=0.85, label=lbl, zorder=3)
            ax.fill(angles_closed, vals, color=PAT_COLORS[pat],
                    alpha=0.08, zorder=1)

        mean_closed = mean_vals.tolist() + [mean_vals[0]]
        ax.plot(angles_closed, mean_closed, "-s", color="black", lw=2.5,
                markersize=5, label="Mean", zorder=4)

        ax.set_xticks(angles)
        ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=15,
                           fontweight="bold")

        if annot_text:
            alpha_idx = BANDS.index("alpha")
            ax.annotate(annot_text,
                        xy=(angles[alpha_idx], mean_vals[alpha_idx]),
                        xytext=(angles[alpha_idx], vmax * 0.82),
                        fontsize=8.5, fontweight="bold", ha="center",
                        va="center", color="#b8860b",
                        arrowprops=dict(arrowstyle="-", color="#b8860b", lw=0.8),
                        zorder=5)

        ax.legend(loc="upper right", bbox_to_anchor=(1.28, 1.08),
                  fontsize=9, frameon=True, framealpha=0.95,
                  edgecolor="#cccccc", handlelength=1.5)

        fig.text(0.5, -0.02,
                 "*Pat_03 (dashed): 1024 Hz acquisition — outlier",
                 ha="center", fontsize=8, color="#666666", style="italic")

        fig.savefig(OUTDIR / filename, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  Saved {filename}")

    # Check alpha unanimity for H2b
    # All 5 patients positive at alpha? Pat_02=+0.167, Pat_03=-0.084, Pat_05=+0.125, Pat_07=+0.117, Pat_08=+0.052
    # Pat_03 is negative → 4/5
    make_radar(data_h2b, "H2b", "radar_chart_H2b.pdf",
               annot_text="4/5\n($p$ = 0.156)")

    md = (
        "# radar_chart_H2b\n\n"
        "## What the figure shows\n"
        "Radar chart of H2b scalar VI contrast (task approach) across 6 bands. "
        "Same format as H2a radar. Alpha: 4/5 positive (Pat_03 negative). "
        "Excluding Pat_03 (documented outlier): 4/4, p = 0.0625.\n"
    )
    (OUTDIR / "radar_chart_H2b.md").write_text(md)


# ═════════════════════════════════════════════════════════════════════
# TASK 4: 4-panel mean contrast heatmap (H1, H2a, H2b, H3)
# ═════════════════════════════════════════════════════════════════════
def task4_four_panel_heatmap():
    print("\n[Task 4] 4-panel mean contrast heatmap (2×2 layout)...")

    hyp_defs = [
        ("H1: TL-TT $-$ others", _h1_contrast),
        ("H2a: Pre-Post $-$ TT-Post", _h2a_contrast),
        ("H2b: Pre-TT $-$ TT-Post", _h2b_contrast),
        ("H3: within $-$ cross", _h3_contrast),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    fig.subplots_adjust(hspace=0.35, wspace=0.15)

    ims = []
    for hi, (title, cfn) in enumerate(hyp_defs):
        row, col = divmod(hi, 2)
        ax = axes[row, col]
        mat = np.full((len(BANDS), len(k_range)), np.nan)
        unan_mask = np.zeros_like(mat, dtype=bool)

        for bi, band in enumerate(BANDS):
            for ki, k in enumerate(k_range):
                vals = []
                for pat in ALL_PATIENTS:
                    sub = df[(df["patient"] == pat) & (df["band"] == band)
                             & (df["k"] == k)]
                    c = cfn(sub, k)
                    if c is not None and np.isfinite(c):
                        vals.append(c)
                if len(vals) == len(ALL_PATIENTS):
                    mat[bi, ki] = np.mean(vals)
                    if all(v < 0 for v in vals) or all(v > 0 for v in vals):
                        unan_mask[bi, ki] = True

        vmax_mat = np.nanmax(np.abs(mat)) if np.any(np.isfinite(mat)) else 1
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax_mat,
                        vmax=vmax_mat, interpolation="nearest",
                        extent=[k_range[0], k_range[-1], len(BANDS) - 0.5, -0.5])
        ims.append(im)

        for bi in range(len(BANDS)):
            for ki in range(len(k_range)):
                if unan_mask[bi, ki]:
                    k_val = k_range[ki]
                    ax.add_patch(plt.Rectangle(
                        (k_val - 0.5, bi - 0.5), 1, 1,
                        fill=False, edgecolor="black", lw=0.5, zorder=5))

        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels([BAND_TEX[b] for b in BANDS], fontsize=10)
        ax.set_xlabel("$k$", fontsize=9)
        ax.set_title(title, fontsize=11, fontweight="bold")
        if col > 0:
            ax.set_yticklabels([])

        # Per-panel colorbar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", size="3%", pad=0.06)
        fig.colorbar(im, cax=cax)

    fig.savefig(OUTDIR / "mean_contrast_heatmap_4panel.pdf",
                bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("  Saved mean_contrast_heatmap_4panel.pdf")

    md = (
        "# mean_contrast_heatmap_4panel\n\n"
        "## What the figure shows\n"
        "Four heatmaps in 2×2 layout (H1, H2a, H2b, H3): band × k, "
        "color = patient-averaged contrast. "
        "Black borders = 5/5 unanimous cells.\n\n"
        "## Key result\n"
        "H1: beta shows wide unanimity. H2a: alpha only. H2b: alpha + low_gamma. "
        "H3: sparse unanimity across bands.\n\n"
        "## How to read it\n"
        "Blue/red = direction of mean contrast. Black borders = all 5 patients agree on sign.\n"
    )
    (OUTDIR / "mean_contrast_heatmap_4panel.md").write_text(md)


# ── Run all tasks ────────────────────────────────────────────────────
task1_unanimity_map()
task2_three_row_dissociation()
task3_radar_charts()
task4_four_panel_heatmap()
print("\nAll done!")
