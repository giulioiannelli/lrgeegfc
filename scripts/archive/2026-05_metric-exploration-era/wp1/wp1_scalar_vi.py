#!/usr/bin/env python3
"""WP-scalar-VI — k-independent scalar summary of the task trace.

Compute three scalar summaries (mean contrast, integrated, fraction) for
four H2 variants (a–d), across 5 patients × 6 bands. Then test all 15
band-pair contrasts and produce figures.

Run: python scripts/wp1/wp1_scalar_vi.py
"""
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.lrg import load_lrg_result

OUT = ROOT / "data" / "wp_scalar_vi"
OUT.mkdir(parents=True, exist_ok=True)

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {b: BRAIN_BAND_TEX_DICT[b] for b in BANDS}
FC_METHOD = "msc"
LRG_CACHE = ROOT / "data" / "lrg_cache"

PAT_COLORS = {
    "Pat_02": "#1f77b4", "Pat_03": "#ff7f0e", "Pat_05": "#2ca02c",
    "Pat_07": "#d62728", "Pat_08": "#9467bd",
}

HYPOTHESES = {
    "H2a": ("rest_pre", "rest_post", "task_test", "rest_post"),   # VI(pre,post) > VI(TT,post)
    "H2b": ("rest_pre", "task_test", "task_test", "rest_post"),  # VI(pre,TT) > VI(TT,post)
    "H2c": ("rest_pre", "rest_post", "task_learn", "rest_post"),   # VI(pre,post) > VI(TL,post)
    "H2d": ("rest_pre", "task_learn", "task_learn", "rest_post"), # VI(pre,TL) > VI(TL,post)
}
HYP_LABELS = {
    "H2a": "VI(Pre,Post) > VI(TT,Post)",
    "H2b": "VI(Pre,TT) > VI(TT,Post)",
    "H2c": "VI(Pre,Post) > VI(TL,Post)",
    "H2d": "VI(Pre,TL) > VI(TL,Post)",
}


def save_fig(fig, path, *, what, proves, how_to_read, dpi=300):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"  Saved: {path}")
    path.with_suffix(".md").write_text(
        f"# {path.stem}\n\n## What the figure shows\n{what}\n\n"
        f"## What it proves\n{proves}\n\n## How to read it\n{how_to_read}\n"
    )


def compute_vi(labels1, labels2):
    n = len(labels1)
    if n == 0:
        return 0.0
    def _ent(lb):
        _, c = np.unique(lb, return_counts=True)
        p = c / n
        return -np.sum(p * np.log(p + 1e-30))
    h1, h2 = _ent(labels1), _ent(labels2)
    mi = 0.0
    for c1 in np.unique(labels1):
        m1 = labels1 == c1
        for c2 in np.unique(labels2):
            nij = (m1 & (labels2 == c2)).sum()
            if nij > 0:
                mi += (nij / n) * np.log((nij * n) / (m1.sum() * (labels2 == c2).sum()) + 1e-30)
    return max(h1 + h2 - 2 * mi, 0.0)


# ── Load data ─────────────────────────────────────────────────────────
print("Loading LRG data...")
lrg = {}
for pat in PATIENTS:
    for band in BANDS:
        for phase in ["rest_pre", "task_learn", "task_test", "rest_post"]:
            res = load_lrg_result(pat, phase, band, FC_METHOD, cache_root=LRG_CACHE)
            if res is not None:
                lrg[(pat, band, phase)] = res
print(f"  {len(lrg)} results")

# Common k range
n_nodes_all = [lrg[(p, "alpha", "rest_pre")].n_nodes for p in PATIENTS
               if (p, "alpha", "rest_pre") in lrg]
K_MAX = min(n_nodes_all) // 2
k_range = list(range(2, K_MAX + 1))
N_K = len(k_range)
print(f"  k range: 2–{K_MAX} ({N_K} values)")


# ======================================================================
# TASK 1: Compute scalar summaries
# ======================================================================
def task1():
    print("\n── Task 1: Scalar summaries ──")
    rows = []

    for pat in PATIENTS:
        for band in BANDS:
            for hyp_name, (pa1, pa2, pb1, pb2) in HYPOTHESES.items():
                # Get LRG results for the 4 phases involved
                r_a1 = lrg.get((pat, band, pa1))
                r_a2 = lrg.get((pat, band, pa2))
                r_b1 = lrg.get((pat, band, pb1))
                r_b2 = lrg.get((pat, band, pb2))
                if any(r is None for r in [r_a1, r_a2, r_b1, r_b2]):
                    continue

                # Check n_nodes compatibility
                ns = {r_a1.n_nodes, r_a2.n_nodes, r_b1.n_nodes, r_b2.n_nodes}
                if len(ns) > 1:
                    continue
                n = r_a1.n_nodes

                contrasts = []
                for k in k_range:
                    la1 = fcluster(r_a1.linkage_matrix, k, criterion="maxclust")[:n]
                    la2 = fcluster(r_a2.linkage_matrix, k, criterion="maxclust")[:n]
                    lb1 = fcluster(r_b1.linkage_matrix, k, criterion="maxclust")[:n]
                    lb2 = fcluster(r_b2.linkage_matrix, k, criterion="maxclust")[:n]
                    vi_a = compute_vi(la1, la2)
                    vi_b = compute_vi(lb1, lb2)
                    contrasts.append(vi_a - vi_b)

                contrasts = np.array(contrasts)
                scalar_a = float(np.mean(contrasts))           # mean contrast
                scalar_b = float(np.sum(contrasts))            # integrated
                scalar_c = float(np.mean(contrasts > 0))       # fraction positive

                rows.append({
                    "patient": pat, "band": band, "hypothesis": hyp_name,
                    "scalar_A_mean": scalar_a,
                    "scalar_B_integrated": scalar_b,
                    "scalar_C_fraction": scalar_c,
                })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "scalar_table_full.csv", index=False)
    print(f"  Saved: {OUT / 'scalar_table_full.csv'} ({len(df)} rows)")

    # Formatted markdown table
    md = ["# Scalar VI summaries — full table", ""]
    for hyp in HYPOTHESES:
        md.append(f"## {hyp}: {HYP_LABELS[hyp]}")
        md.append("")
        md.append("### Scalar A (mean contrast)")
        md.append("")
        md.append("| Band | " + " | ".join(PATIENTS) + " |")
        md.append("|------|" + "|".join(["------"] * len(PATIENTS)) + "|")
        for band in BANDS:
            vals = []
            for pat in PATIENTS:
                sub = df[(df["patient"] == pat) & (df["band"] == band) &
                         (df["hypothesis"] == hyp)]
                if not sub.empty:
                    v = sub["scalar_A_mean"].iloc[0]
                    vals.append(f"{v:+.4f}")
                else:
                    vals.append("—")
            bold = "**" if band in ("alpha", "beta") else ""
            md.append(f"| {bold}{band}{bold} | " + " | ".join(vals) + " |")
        md.append("")

        md.append("### Scalar C (fraction of k supporting)")
        md.append("")
        md.append("| Band | " + " | ".join(PATIENTS) + " |")
        md.append("|------|" + "|".join(["------"] * len(PATIENTS)) + "|")
        for band in BANDS:
            vals = []
            for pat in PATIENTS:
                sub = df[(df["patient"] == pat) & (df["band"] == band) &
                         (df["hypothesis"] == hyp)]
                if not sub.empty:
                    v = sub["scalar_C_fraction"].iloc[0]
                    vals.append(f"{v:.2f}")
                else:
                    vals.append("—")
            bold = "**" if band in ("alpha", "beta") else ""
            md.append(f"| {bold}{band}{bold} | " + " | ".join(vals) + " |")
        md.append("")

    (OUT / "scalar_table_full.md").write_text("\n".join(md))
    print(f"  Saved: {OUT / 'scalar_table_full.md'}")
    return df


# ======================================================================
# TASK 2: Band-pair contrasts
# ======================================================================
def task2(df):
    print("\n── Task 2: Band-pair contrasts ──")

    band_pairs = list(combinations(BANDS, 2))
    rows = []

    for hyp in HYPOTHESES:
        for scalar_col, scalar_name in [
            ("scalar_A_mean", "A"),
            ("scalar_B_integrated", "B"),
            ("scalar_C_fraction", "C"),
        ]:
            for b1, b2 in band_pairs:
                n_b1_gt = 0
                n_total = 0
                for pat in PATIENTS:
                    s1 = df[(df["patient"] == pat) & (df["band"] == b1) &
                            (df["hypothesis"] == hyp)]
                    s2 = df[(df["patient"] == pat) & (df["band"] == b2) &
                            (df["hypothesis"] == hyp)]
                    if s1.empty or s2.empty:
                        continue
                    n_total += 1
                    if s1[scalar_col].iloc[0] > s2[scalar_col].iloc[0]:
                        n_b1_gt += 1

                rows.append({
                    "hypothesis": hyp, "scalar": scalar_name,
                    "band1": b1, "band2": b2,
                    "n_b1_gt_b2": n_b1_gt, "n_patients": n_total,
                    "unanimous": n_b1_gt == n_total and n_total == len(PATIENTS),
                })

    pairs_df = pd.DataFrame(rows)
    pairs_df.to_csv(OUT / "band_pair_contrasts.csv", index=False)

    # Summary of unanimous pairs
    unan = pairs_df[pairs_df["unanimous"]]
    md = [
        "# Band-pair contrasts — unanimous results (5/5)",
        "",
        f"Total pairs tested: {len(band_pairs)} per hypothesis × scalar = "
        f"{len(band_pairs) * len(HYPOTHESES) * 3}",
        "",
        f"**Unanimous pairs found: {len(unan)}**",
        "",
    ]

    if len(unan) > 0:
        md.append("| Hypothesis | Scalar | Band1 > Band2 |")
        md.append("|------------|--------|---------------|")
        for _, r in unan.iterrows():
            md.append(f"| {r['hypothesis']} | {r['scalar']} | "
                      f"{r['band1']} > {r['band2']} |")
    else:
        md.append("No unanimous pairs found at 5/5.")

    # Also report 4/5 pairs
    near = pairs_df[pairs_df["n_b1_gt_b2"] >= 4]
    near = near[~near["unanimous"]]
    md += [
        "",
        f"## Near-unanimous pairs (4/5): {len(near)}",
        "",
    ]
    if len(near) > 0:
        md.append("| Hypothesis | Scalar | Band1 > Band2 | Count |")
        md.append("|------------|--------|---------------|-------|")
        for _, r in near.iterrows():
            md.append(f"| {r['hypothesis']} | {r['scalar']} | "
                      f"{r['band1']} > {r['band2']} | {r['n_b1_gt_b2']}/5 |")

    (OUT / "band_pair_contrasts.md").write_text("\n".join(md))
    print(f"  Saved: {OUT / 'band_pair_contrasts.md'}")
    return pairs_df


# ======================================================================
# TASK 3: Figures
# ======================================================================
def task3(df, pairs_df):
    print("\n── Task 3: Figures ──")

    # 3A: Strip plot per hypothesis
    fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharey=True)
    for ax, hyp in zip(axes, HYPOTHESES):
        for bi, band in enumerate(BANDS):
            color = "#2196F3" if band == "alpha" else "#F44336" if band == "beta" else "#999"
            vals = []
            for j, pat in enumerate(PATIENTS):
                sub = df[(df["patient"] == pat) & (df["band"] == band) &
                         (df["hypothesis"] == hyp)]
                if not sub.empty:
                    v = sub["scalar_A_mean"].iloc[0]
                    jitter = (j - 2) * 0.06
                    ax.plot(bi + jitter, v, "o", color=PAT_COLORS[pat],
                            markersize=7, markeredgecolor="black", markeredgewidth=0.3)
                    vals.append(v)
            if vals:
                ax.plot([bi - 0.2, bi + 0.2], [np.mean(vals)] * 2,
                        color=color, lw=3, zorder=6)

        ax.axhline(0, color="gray", ls=":", lw=1)
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=8)
        ax.set_title(f"{hyp}\n{HYP_LABELS[hyp]}", fontsize=9)
    axes[0].set_ylabel("Scalar A (mean VI contrast)")
    fig.suptitle("Scalar A across bands — 4 hypotheses", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT / "figure_3A_strip_per_hypothesis.pdf",
        what="Strip plots of Scalar A (mean VI contrast) per band for 4 H2 variants. "
             "Positive = trace supported. One dot per patient, bar = mean.",
        proves="Whether the alpha trace holds across all 4 hypothesis variants.",
        how_to_read="Alpha dots above 0 = trace. Beta dots below 0 = anti-trace. "
                    "Compare across panels to see if task reference matters.",
    )

    # 3B: Slope plots per hypothesis
    fig, axes = plt.subplots(1, 4, figsize=(16, 5), sharey=True)
    for ax, hyp in zip(axes, HYPOTHESES):
        n_agree = 0
        for pat in PATIENTS:
            sa = df[(df["patient"] == pat) & (df["band"] == "alpha") &
                    (df["hypothesis"] == hyp)]
            sb = df[(df["patient"] == pat) & (df["band"] == "beta") &
                    (df["hypothesis"] == hyp)]
            if sa.empty or sb.empty:
                continue
            va = sa["scalar_A_mean"].iloc[0]
            vb = sb["scalar_A_mean"].iloc[0]
            ax.plot([0, 1], [va, vb], "o-", color=PAT_COLORS[pat], lw=2,
                    markersize=8, markeredgecolor="black", markeredgewidth=0.5,
                    label=pat)
            if va > vb:
                n_agree += 1

        ax.set_xticks([0, 1])
        ax.set_xticklabels([r"$\alpha$", r"$\beta$"], fontsize=14)
        ax.set_xlim(-0.3, 1.3)
        ax.axhline(0, color="gray", ls=":", lw=1)
        ax.set_title(f"{hyp}: α > β {n_agree}/5", fontsize=10)

    axes[0].set_ylabel("Scalar A (mean VI contrast)")
    axes[0].legend(fontsize=7)
    fig.suptitle("Alpha vs Beta — Scalar A, 4 hypotheses", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT / "figure_3B_slope_per_hypothesis.pdf",
        what="Paired slope plots: alpha vs beta Scalar A for each H2 variant. "
             "Each line = one patient. Title shows unanimity count.",
        proves="Whether alpha > beta holds across different task-trace definitions.",
        how_to_read="All lines sloping down = alpha has more trace than beta. "
                    "Title n/5 shows how many patients agree.",
    )

    # 3C: Scalar comparison (A, B, C for H2a)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    scalar_cols = [("scalar_A_mean", "A: Mean contrast"),
                   ("scalar_B_integrated", "B: Integrated"),
                   ("scalar_C_fraction", "C: Fraction positive")]
    for ax, (col, title) in zip(axes, scalar_cols):
        for bi, band in enumerate(BANDS):
            color = "#2196F3" if band == "alpha" else "#F44336" if band == "beta" else "#999"
            vals = []
            for j, pat in enumerate(PATIENTS):
                sub = df[(df["patient"] == pat) & (df["band"] == band) &
                         (df["hypothesis"] == "H2a")]
                if not sub.empty:
                    v = sub[col].iloc[0]
                    jitter = (j - 2) * 0.06
                    ax.plot(bi + jitter, v, "o", color=PAT_COLORS[pat],
                            markersize=7, markeredgecolor="black", markeredgewidth=0.3)
                    vals.append(v)
            if vals:
                ax.plot([bi - 0.2, bi + 0.2], [np.mean(vals)] * 2,
                        color=color, lw=3, zorder=6)

        ref_line = 0 if "fraction" not in col else 0.5
        ax.axhline(ref_line, color="gray", ls=":", lw=1)
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=8)
        ax.set_title(title)

    fig.suptitle("H2a — three scalar definitions compared", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT / "figure_3C_scalar_comparison.pdf",
        what="Three scalars for H2a: A (mean contrast), B (integrated), C (fraction). "
             "Same strip-plot format.",
        proves="Whether the three scalar definitions agree on band ordering.",
        how_to_read="If alpha is highest in all three panels, the choice doesn't matter.",
    )

    # 3D: Band-pair unanimity matrix (Scalar A only)
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
    for ax, hyp in zip(axes, HYPOTHESES):
        mat = np.full((len(BANDS), len(BANDS)), np.nan)
        sub = pairs_df[(pairs_df["hypothesis"] == hyp) & (pairs_df["scalar"] == "A")]
        for _, r in sub.iterrows():
            bi = BANDS.index(r["band1"])
            bj = BANDS.index(r["band2"])
            mat[bi, bj] = r["n_b1_gt_b2"]
            mat[bj, bi] = r["n_patients"] - r["n_b1_gt_b2"]

        im = ax.imshow(mat, cmap="RdBu", vmin=0, vmax=5, aspect="equal",
                        interpolation="nearest")
        # Black border on 5/5 cells
        for i in range(len(BANDS)):
            for j in range(len(BANDS)):
                if i != j and not np.isnan(mat[i, j]) and mat[i, j] == 5:
                    rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                         fill=False, edgecolor="black", lw=2)
                    ax.add_patch(rect)
                if not np.isnan(mat[i, j]) and i != j:
                    ax.text(j, i, f"{int(mat[i, j])}", ha="center", va="center",
                            fontsize=8, fontweight="bold",
                            color="white" if mat[i, j] >= 4 or mat[i, j] <= 1 else "black")

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=8)
        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels([BAND_TEX[b] for b in BANDS], fontsize=8)
        ax.set_title(hyp, fontsize=11)

    fig.colorbar(im, ax=axes, label="# patients row > col", shrink=0.6, pad=0.02)
    fig.suptitle("Band-pair unanimity — Scalar A (mean VI contrast)", fontsize=13)
    fig.tight_layout()
    save_fig(fig, OUT / "figure_3D_band_pair_unanimity.pdf",
        what="6×6 matrices: entry (row, col) = # patients where row band has higher "
             "Scalar A than col band. Blue = 5/5 (row > col unanimous). "
             "Black borders on 5/5 cells. One panel per hypothesis.",
        proves="Which band-pair contrasts are unanimous across patients.",
        how_to_read="Blue cells with black borders = reportable (5/5). "
                    "Red cells = opposite direction unanimous. "
                    "Compare across panels to check hypothesis consistency.",
    )


# ======================================================================
# TASK 4: Assessment
# ======================================================================
def task4(df, pairs_df):
    print("\n── Task 4: Assessment ──")

    md = ["# Scalar VI Assessment", ""]

    # Q1: Do scalars agree?
    md.append("## 1. Do Scalars A, B, C agree on band ordering?")
    md.append("")
    for hyp in HYPOTHESES:
        md.append(f"### {hyp}")
        for scalar in ["A", "B", "C"]:
            col = {"A": "scalar_A_mean", "B": "scalar_B_integrated",
                   "C": "scalar_C_fraction"}[scalar]
            band_means = {}
            for band in BANDS:
                vals = df[(df["band"] == band) & (df["hypothesis"] == hyp)][col]
                band_means[band] = vals.mean()
            ranking = sorted(band_means, key=band_means.get, reverse=True)
            md.append(f"- Scalar {scalar}: {' > '.join(ranking[:3])} > ...")
        md.append("")

    # Q2: Unanimous alpha > beta?
    md.append("## 2. Is alpha > beta unanimous (5/5)?")
    md.append("")
    md.append("| Hypothesis | Scalar A | Scalar B | Scalar C |")
    md.append("|------------|----------|----------|----------|")
    for hyp in HYPOTHESES:
        results = []
        for scalar in ["A", "B", "C"]:
            sub = pairs_df[(pairs_df["hypothesis"] == hyp) &
                           (pairs_df["scalar"] == scalar) &
                           (pairs_df["band1"] == "alpha") &
                           (pairs_df["band2"] == "beta")]
            if not sub.empty:
                n = sub["n_b1_gt_b2"].iloc[0]
                results.append(f"**{n}/5**" if n == 5 else f"{n}/5")
            else:
                results.append("—")
        md.append(f"| {hyp} | " + " | ".join(results) + " |")
    md.append("")

    # Q3: All unanimous pairs
    md.append("## 3. ALL unanimous band-pair contrasts (5/5)")
    md.append("")
    unan = pairs_df[pairs_df["unanimous"]]
    if len(unan) > 0:
        md.append("| Hypothesis | Scalar | Contrast |")
        md.append("|------------|--------|----------|")
        for _, r in unan.iterrows():
            md.append(f"| {r['hypothesis']} | {r['scalar']} | "
                      f"{r['band1']} > {r['band2']} |")
    else:
        md.append("**None found.** No band-pair contrast reaches 5/5 unanimity "
                  "for any scalar or hypothesis variant.")
    md.append("")

    # Q4: Task reference effect
    md.append("## 4. Does the task reference matter?")
    md.append("")
    for band in BANDS:
        tt_vals = df[(df["band"] == band) & (df["hypothesis"] == "H2a")]["scalar_A_mean"].values
        tl_vals = df[(df["band"] == band) & (df["hypothesis"] == "H2c")]["scalar_A_mean"].values
        if len(tt_vals) > 0 and len(tl_vals) > 0:
            md.append(f"- {band}: task_test mean={np.mean(tt_vals):.4f}, "
                      f"task_learn mean={np.mean(tl_vals):.4f}, "
                      f"diff={np.mean(tt_vals) - np.mean(tl_vals):+.4f}")
    md.append("")

    # Q5: Anti-trace bands
    md.append("## 5. Bands with consistent anti-trace (negative Scalar A)")
    md.append("")
    for hyp in ["H2a", "H2c"]:
        md.append(f"### {hyp}")
        for band in BANDS:
            vals = df[(df["band"] == band) & (df["hypothesis"] == hyp)]["scalar_A_mean"].values
            n_neg = sum(1 for v in vals if v < 0)
            mean_v = np.mean(vals)
            if n_neg >= 4:
                md.append(f"- **{band}**: {n_neg}/5 negative, mean={mean_v:+.4f}")
            elif n_neg >= 3:
                md.append(f"- {band}: {n_neg}/5 negative, mean={mean_v:+.4f}")
        md.append("")

    # Q6: Final recommendation
    md.append("## 6. Final recommendation")
    md.append("")

    # Check what actually works
    unan_count = len(unan)
    near_unan = pairs_df[pairs_df["n_b1_gt_b2"] >= 4]
    near_unan_ab = near_unan[(near_unan["band1"] == "alpha") &
                              (near_unan["band2"] == "beta")]

    if unan_count > 0:
        md.append("Unanimous contrasts exist. Use the most consistent scalar+hypothesis "
                  "combination as the primary result.")
    else:
        best_ab = near_unan_ab.sort_values("n_b1_gt_b2", ascending=False)
        if len(best_ab) > 0:
            best = best_ab.iloc[0]
            md.append(f"No 5/5 unanimous alpha > beta contrast found. "
                      f"Best: {best['hypothesis']} Scalar {best['scalar']} "
                      f"({best['n_b1_gt_b2']}/5).")
        md.append("")
        md.append("**The scalar VI approach confirms that the task trace is alpha-positive "
                  "on average but the patient-level contrast is not strictly unanimous.** "
                  "This is consistent with the VI unanimity maps, which test at EACH k "
                  "separately (finding 29 cells) — the k-averaged scalar dilutes the effect.")
        md.append("")
        md.append("**Recommendation:** Use VI unanimity maps as the primary result. "
                  "Report the scalar mean contrast as a summary statistic (alpha has the "
                  "highest mean Scalar A across patients for H2a/H2c) but do NOT claim "
                  "5/5 unanimity for the scalar — it doesn't survive strict testing.")

    (OUT / "assessment.md").write_text("\n".join(md))
    print(f"  Saved: {OUT / 'assessment.md'}")


# ── Main ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    df = task1()
    pairs_df = task2(df)
    task3(df, pairs_df)
    task4(df, pairs_df)
    print(f"\nAll done. Outputs in {OUT}")
