#!/usr/bin/env python3
"""WP0 Task 3 (rewrite) — Multiscale hypothesis testing with VI and unanimity.

For each (patient, band, phase-pair), compute VI(k) and NVI(k) across the
full range of dendrogram cut levels k = 2, 3, ..., N/2. Then test four
hypotheses AT EACH SCALE using signed unanimity across all 5 patients.

Hypotheses:
  H1 (task stability):  VI(taskL, taskT) < mean(other pairs)
  H2a (task trace):     VI(rest_pre, rest_post) > VI(taskT, rest_post)
  H2b (task approach):  VI(taskT, rest_post) < VI(rest_pre, taskT)
  H3 (within < cross):  mean(VI_within) < mean(VI_cross)

Unanimity: with N=5, we require 5/5 patients to agree on the sign of each
contrast at each (band, k). This is the ONLY valid aggregation for N=5
(Wilcoxon minimum p = 0.031 = P(5/5) under H0).

Outputs (under data/wp0_metric_exploration/):
  task3_multiscale/
    vi_profiles_by_band.pdf + .md
    unanimity_map_H1.pdf + .md
    unanimity_map_H2a.pdf + .md
    unanimity_map_H2b.pdf + .md
    unanimity_map_H3.pdf + .md
    unanimity_summary.md
    vi_raw_profiles.csv
    raw_fc_comparison.md

Run: python scripts/wp0/task3_multiscale_hypotheses.py
"""
from __future__ import annotations

import sys
import time
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "wp0"))
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.cluster.hierarchy import fcluster

from _common import (
    ALL_PAIRS,
    BANDS,
    FC_METHOD,
    LRG_CACHE,
    OUT_ROOT,
    PATIENTS,
    PHASES,
    classify_pair,
    load_all_lrg,
    load_all_fc,
    save_fig,
)
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

OUT = OUT_ROOT / "task3_multiscale"
BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER]

# Phase pair definitions
WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]

PAIR_SHORT = {
    ("rest_pre", "rest_post"): "Pre-Post",
    ("task_learn", "task_test"): "TL-TT",
    ("rest_pre", "task_learn"): "Pre-TL",
    ("rest_pre", "task_test"): "Pre-TT",
    ("task_learn", "rest_post"): "TL-Post",
    ("task_test", "rest_post"): "TT-Post",
}

PAIR_COLORS = {
    ("rest_pre", "rest_post"): "#2196F3",
    ("task_learn", "task_test"): "#FF9800",
    ("rest_pre", "task_learn"): "#E91E63",
    ("rest_pre", "task_test"): "#9C27B0",
    ("task_learn", "rest_post"): "#F44336",
    ("task_test", "rest_post"): "#795548",
}


# ── VI computation ────────────────────────────────────────────────────────
def compute_vi(labels1: np.ndarray, labels2: np.ndarray) -> float:
    """Variation of Information: VI = H(P) + H(Q) - 2*MI(P,Q). Lower = more similar."""
    n = len(labels1)
    if n == 0:
        return 0.0
    classes1, classes2 = np.unique(labels1), np.unique(labels2)

    def _entropy(labels):
        _, counts = np.unique(labels, return_counts=True)
        p = counts / n
        return -np.sum(p * np.log(p + 1e-30))

    h1, h2 = _entropy(labels1), _entropy(labels2)

    mi = 0.0
    for c1 in classes1:
        m1 = labels1 == c1
        n1 = m1.sum()
        for c2 in classes2:
            nij = (m1 & (labels2 == c2)).sum()
            if nij > 0:
                mi += (nij / n) * np.log((nij * n) / (n1 * (labels2 == c2).sum()) + 1e-30)

    return max(h1 + h2 - 2 * mi, 0.0)


# ── Step 1: Compute VI(k) for all (patient, band, pair, k) ───────────────
def compute_all_vi_profiles(lrg_data: dict) -> pd.DataFrame:
    """Compute VI at every k from 2 to N/2 for all combinations."""
    print("Computing VI profiles...")
    rows = []

    for pat in PATIENTS:
        for band in BAND_ORDER:
            # Get all phase results
            results = {}
            for phase in PHASES:
                key = (pat, phase, band)
                if key in lrg_data:
                    results[phase] = lrg_data[key]

            if len(results) < 4:
                print(f"  SKIP {pat} {band}: only {len(results)} phases")
                continue

            # Max k = min(n_nodes)/2 across phases
            n_nodes = min(r.n_nodes for r in results.values())
            max_k = max(n_nodes // 2, 3)
            k_range = list(range(2, max_k + 1))

            for pa, pb in ALL_PAIRS:
                if pa not in results or pb not in results:
                    continue
                ra, rb = results[pa], results[pb]

                # Require same n_nodes for partition comparison
                if ra.n_nodes != rb.n_nodes:
                    continue

                Za, Zb = ra.linkage_matrix, rb.linkage_matrix

                for k in k_range:
                    la = fcluster(Za, k, criterion="maxclust")
                    lb = fcluster(Zb, k, criterion="maxclust")
                    vi = compute_vi(la, lb)

                    # NVI = VI / ln(k) — normalized for cross-scale comparison
                    nvi = vi / np.log(max(k, 2))

                    rows.append({
                        "patient": pat,
                        "band": band,
                        "phase_a": pa,
                        "phase_b": pb,
                        "pair": PAIR_SHORT.get((pa, pb), f"{pa}-{pb}"),
                        "pair_type": classify_pair(pa, pb),
                        "k": k,
                        "vi": vi,
                        "nvi": nvi,
                    })

        print(f"  {pat} done")

    df = pd.DataFrame(rows)
    print(f"  Total: {len(df)} profile points")
    return df


# ── Step 2: Hypothesis testing at each (band, k) ─────────────────────────
def compute_unanimity_maps(df: pd.DataFrame) -> dict:
    """Test H1, H2a, H2b, H3 at every (band, k) with signed unanimity.

    Returns dict of hypothesis name -> DataFrame with columns:
      band, k, unanimity (float in {-1, -0.x, 0, +0.x, +1}), n_agree, n_patients
    """
    results = {}

    for hyp_name, contrast_fn in [
        ("H1_task_stability", _contrast_H1),
        ("H2a_task_trace", _contrast_H2a),
        ("H2b_task_approach", _contrast_H2b),
        ("H3_within_vs_cross", _contrast_H3),
    ]:
        rows = []
        for band in BAND_ORDER:
            bdf = df[df["band"] == band]
            k_values = sorted(bdf["k"].unique())

            for k in k_values:
                kdf = bdf[bdf["k"] == k]
                signs = []

                for pat in PATIENTS:
                    pdf = kdf[kdf["patient"] == pat]
                    if pdf.empty:
                        continue
                    contrast = contrast_fn(pdf)
                    if contrast is not None and not np.isnan(contrast):
                        signs.append(np.sign(contrast))

                n_patients = len(signs)
                if n_patients == 0:
                    continue

                unanimity = np.mean(signs)  # ±1 if all agree
                n_positive = sum(1 for s in signs if s > 0)
                n_negative = sum(1 for s in signs if s < 0)

                rows.append({
                    "band": band,
                    "k": k,
                    "unanimity": unanimity,
                    "n_positive": n_positive,
                    "n_negative": n_negative,
                    "n_patients": n_patients,
                    "is_unanimous": abs(unanimity) == 1.0 and n_patients == len(PATIENTS),
                })

        results[hyp_name] = pd.DataFrame(rows)
        n_unan = sum(1 for _, r in pd.DataFrame(rows).iterrows() if r.get("is_unanimous", False))
        print(f"  {hyp_name}: {n_unan} unanimous (band, k) cells")

    return results


def _get_vi(pdf: pd.DataFrame, pa: str, pb: str) -> float:
    """Get VI value for a specific pair from a patient-k slice."""
    row = pdf[((pdf["phase_a"] == pa) & (pdf["phase_b"] == pb)) |
              ((pdf["phase_a"] == pb) & (pdf["phase_b"] == pa))]
    if row.empty:
        return np.nan
    return row["vi"].iloc[0]


def _contrast_H1(pdf: pd.DataFrame) -> float:
    """H1: task stability. POSITIVE = H1 supported.
    Contrast: mean(other pairs VI) - VI(taskL, taskT)."""
    vi_tt = _get_vi(pdf, "task_learn", "task_test")
    other_vis = []
    for pa, pb in ALL_PAIRS:
        if (pa, pb) == ("task_learn", "task_test"):
            continue
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            other_vis.append(v)
    if np.isnan(vi_tt) or len(other_vis) == 0:
        return np.nan
    return np.mean(other_vis) - vi_tt


def _contrast_H2a(pdf: pd.DataFrame) -> float:
    """H2a: task trace — rest_post closer to task than rest_pre is.
    POSITIVE = trace detected.
    Contrast: VI(rest_pre, rest_post) - VI(taskT, rest_post)."""
    return _get_vi(pdf, "rest_pre", "rest_post") - _get_vi(pdf, "task_test", "rest_post")


def _contrast_H2b(pdf: pd.DataFrame) -> float:
    """H2b: task approach — rest_post moved toward task.
    POSITIVE = approach detected.
    Contrast: VI(rest_pre, taskT) - VI(taskT, rest_post)."""
    return _get_vi(pdf, "rest_pre", "task_test") - _get_vi(pdf, "task_test", "rest_post")


def _contrast_H3(pdf: pd.DataFrame) -> float:
    """H3: within-type < cross-type.
    POSITIVE = H3 supported.
    Contrast: mean(cross VI) - mean(within VI)."""
    within_vis = []
    for pa, pb in WITHIN_PAIRS:
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            within_vis.append(v)
    cross_vis = []
    for pa, pb in CROSS_PAIRS:
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            cross_vis.append(v)
    if not within_vis or not cross_vis:
        return np.nan
    return np.mean(cross_vis) - np.mean(within_vis)


# ── Step 3: Figures ───────────────────────────────────────────────────────
def plot_vi_profiles(df: pd.DataFrame):
    """VI(k) profiles faceted by band, one line per (patient, pair)."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for idx, band in enumerate(BAND_ORDER):
        ax = axes[idx]
        bdf = df[df["band"] == band]

        for (pat, pa, pb), grp in bdf.groupby(["patient", "phase_a", "phase_b"]):
            grp = grp.sort_values("k")
            pair = (pa, pb)
            color = PAIR_COLORS.get(pair, "gray")
            ptype = classify_pair(pa, pb)
            lw = 1.8 if ptype == "within" else 0.7
            alpha = 0.7 if ptype == "within" else 0.3
            ax.plot(grp["k"], grp["vi"], color=color, alpha=alpha, lw=lw)

        ax.set_xlabel("k (number of clusters)", fontsize=9)
        ax.set_ylabel("VI (nats)", fontsize=9)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=12)

    # Legend in first panel
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=PAIR_COLORS[p], lw=2,
               label=PAIR_SHORT[p]) for p in ALL_PAIRS if p in PAIR_COLORS]
    axes[0].legend(handles=handles, fontsize=6, loc="upper left", ncol=2)

    fig.suptitle("VI(k) profiles — all patients, all pairs", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT / "vi_profiles_by_band.pdf",
        what="VI(k) curves for every (patient, phase-pair) combination, faceted by band. "
             "X = number of clusters k. Y = Variation of Information (nats). "
             "Thick colored lines = within-type pairs (Pre-Post blue, TL-TT orange). "
             "Thin lines = cross-type pairs (pink/purple/red/brown). Each patient is a separate line.",
        proves="Shows the FULL multiscale landscape. Within-type pairs (especially TL-TT) should have "
               "lower VI at most scales. If separation only exists at certain k ranges, those are the "
               "informative scales. If lines cross, the hypothesis is scale-dependent.",
        how_to_read="Lower = more similar. TL-TT (orange) should be lowest if H1 holds. "
                    "Pre-Post (blue) height relative to cross-pairs tests H2/H3. "
                    "Look for k ranges where within-type curves are clearly separated from cross-type.",
    )


def plot_unanimity_map(unan_df: pd.DataFrame, hyp_name: str, title: str, description: str):
    """Heatmap: band × k, colored by signed unanimity."""
    fig, ax = plt.subplots(figsize=(16, 4.5))

    # Pivot to matrix: rows = bands, columns = k
    piv = unan_df.pivot_table(index="band", columns="k", values="unanimity")
    piv = piv.reindex(BAND_ORDER)

    # Colormap: red=-1 (all disagree), white=0 (split), blue=+1 (all agree)
    cmap = plt.cm.RdBu
    im = ax.imshow(piv.values, aspect="auto", cmap=cmap, vmin=-1, vmax=1,
                   interpolation="nearest")

    # Mark unanimous cells with black border
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if not np.isnan(v) and abs(v) == 1.0:
                # Check this is truly 5/5
                band = BAND_ORDER[i]
                k = piv.columns[j]
                row = unan_df[(unan_df["band"] == band) & (unan_df["k"] == k)]
                if not row.empty and row.iloc[0]["n_patients"] == len(PATIENTS):
                    rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                         fill=False, edgecolor="black", lw=1.5)
                    ax.add_patch(rect)

    ax.set_yticks(range(len(BAND_ORDER)))
    ax.set_yticklabels(BAND_TEX, fontsize=10)

    # Sparse k labels
    k_vals = list(piv.columns)
    step = max(1, len(k_vals) // 20)
    tick_idx = list(range(0, len(k_vals), step))
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([k_vals[t] for t in tick_idx], fontsize=7)
    ax.set_xlabel("k (number of clusters)", fontsize=10)

    fig.colorbar(im, ax=ax, label="Signed unanimity (5 patients)", shrink=0.8)
    ax.set_title(f"{title} — unanimity map (5/5 patients)", fontsize=12)
    fig.tight_layout()

    save_fig(fig, OUT / f"unanimity_map_{hyp_name}.pdf",
        what=f"Heatmap: band (rows) vs k (columns). Color = signed unanimity across 5 patients. "
             f"Blue (+1) = all 5 patients support the hypothesis at this (band, k). "
             f"Red (-1) = all 5 patients show the opposite. White = split. "
             f"Black borders mark cells with |unanimity| = 1.0 and N = 5.",
        proves=f"{description}",
        how_to_read="ONLY blue cells with black borders are reportable results (5/5 agreement). "
                    "Contiguous blue regions across k = robust effect. "
                    "Band-specific patterns (blue in some bands, white/red in others) = band-dependent effect. "
                    "If no blue regions exist, the hypothesis is not supported.",
    )


def write_unanimity_summary(unan_maps: dict, vi_df: pd.DataFrame):
    """Write detailed summary of where hypotheses hold (strict 5/5)."""
    lines = [
        "# WP0 — Multiscale Hypothesis Testing Summary",
        "",
        f"**N = {len(PATIENTS)} patients**: {', '.join(PATIENTS)}",
        f"**Criterion**: Strict unanimity — 5/5 patients must agree on the sign of the contrast",
        f"**Metric**: Variation of Information (VI), computed at each k = 2, ..., N/2",
        "",
        "With N=5, 5/5 agreement gives Wilcoxon one-sided p = 0.031 (minimum achievable).",
        "",
        "---",
        "",
    ]

    for hyp_name, title, description in [
        ("H1_task_stability", "H1: Task Stability",
         "VI(TL,TT) < mean(other pairs) → task phases are most similar"),
        ("H2a_task_trace", "H2a: Task Trace",
         "VI(rest_pre,rest_post) > VI(taskT,rest_post) → rest_post retains task structure"),
        ("H2b_task_approach", "H2b: Task Approach",
         "VI(rest_pre,taskT) > VI(taskT,rest_post) → rest_post moved toward task"),
        ("H3_within_vs_cross", "H3: Within < Cross",
         "mean(VI_within) < mean(VI_cross) → condition types are discriminable"),
    ]:
        udf = unan_maps[hyp_name]
        lines.append(f"## {title}")
        lines.append(f"*{description}*")
        lines.append("")

        # Find all unanimous (band, k) cells
        unan_cells = udf[udf["is_unanimous"] == True]
        pos_cells = unan_cells[unan_cells["unanimity"] > 0]
        neg_cells = unan_cells[unan_cells["unanimity"] < 0]

        if pos_cells.empty and neg_cells.empty:
            lines.append("**No unanimous (band, k) cells found.** Hypothesis not supported at any scale.")
            lines.append("")
            # Still report best-performing bands
            best = udf.groupby("band")["unanimity"].mean().sort_values(ascending=False)
            lines.append("Best-performing bands (by mean unanimity, not significant):")
            for band, u in best.items():
                lines.append(f"  - {band}: mean unanimity = {u:.3f}")
            lines.append("")
        else:
            if not pos_cells.empty:
                lines.append(f"### Supported (+1 unanimity): {len(pos_cells)} (band, k) cells")
                lines.append("")
                lines.append("| Band | k range | # cells |")
                lines.append("|------|---------|---------|")
                for band in BAND_ORDER:
                    bc = pos_cells[pos_cells["band"] == band]
                    if bc.empty:
                        continue
                    k_min, k_max = bc["k"].min(), bc["k"].max()
                    lines.append(f"| {band} | {k_min}–{k_max} | {len(bc)} |")
                lines.append("")

            if not neg_cells.empty:
                lines.append(f"### Reversed (-1 unanimity): {len(neg_cells)} (band, k) cells")
                lines.append("")
                lines.append("| Band | k range | # cells |")
                lines.append("|------|---------|---------|")
                for band in BAND_ORDER:
                    bc = neg_cells[neg_cells["band"] == band]
                    if bc.empty:
                        continue
                    k_min, k_max = bc["k"].min(), bc["k"].max()
                    lines.append(f"| {band} | {k_min}–{k_max} | {len(bc)} |")
                lines.append("")

        lines.append("---")
        lines.append("")

    # Raw FC comparison section
    lines.append("## Raw FC vs LRG")
    lines.append("")
    lines.append("Raw FC metrics (Frobenius, Spearman, etc.) compare edge-level FC matrices.")
    lines.append("They CANNOT provide scale-resolved information. The question is not")
    lines.append("'does LRG detect effects better?' but 'does LRG reveal WHERE in the")
    lines.append("hierarchy the effects live?'")
    lines.append("")
    lines.append("If unanimous (band, k) cells exist above, they demonstrate that the task")
    lines.append("effect is localized at specific hierarchical scales — information that raw")
    lines.append("FC comparison cannot provide.")
    lines.append("")

    (OUT / "unanimity_summary.md").write_text("\n".join(lines))
    print(f"  Summary written to {OUT / 'unanimity_summary.md'}")


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("Loading LRG data...")
    lrg_data = load_all_lrg()

    vi_df = compute_all_vi_profiles(lrg_data)

    # Save raw profiles
    csv_path = OUT / "vi_raw_profiles.csv"
    vi_df.to_csv(csv_path, index=False)
    print(f"  Saved {csv_path}")

    print("\nTesting hypotheses with unanimity...")
    unan_maps = compute_unanimity_maps(vi_df)

    print("\nGenerating figures...")
    plot_vi_profiles(vi_df)

    for hyp_name, title, desc in [
        ("H1_task_stability", "H1: Task Stability",
         "Tests whether task_learn–task_test VI is lower than all other pairs — "
         "i.e., the two task phases share a common community architecture."),
        ("H2a_task_trace", "H2a: Task Trace (persistence)",
         "Tests whether VI(rest_pre,rest_post) > VI(taskT,rest_post) — "
         "i.e., rest_post is closer to the task state than rest_pre was, meaning "
         "the task left a structural trace in post-rest."),
        ("H2b_task_approach", "H2b: Task Approach",
         "Tests whether VI(rest_pre,taskT) > VI(taskT,rest_post) — "
         "i.e., rest_post moved toward the task-test partition."),
        ("H3_within_vs_cross", "H3: Within < Cross",
         "Tests whether within-type pairs (rest-rest + task-task) have lower VI "
         "than cross-type pairs — i.e., cognitive state shapes community structure."),
    ]:
        plot_unanimity_map(unan_maps[hyp_name], hyp_name, title, desc)

    write_unanimity_summary(unan_maps, vi_df)

    elapsed = time.time() - t0
    print(f"\nDone ({elapsed:.1f}s). Outputs in {OUT}")


if __name__ == "__main__":
    main()
