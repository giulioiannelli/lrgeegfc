#!/usr/bin/env python3
"""Section 3 Figure J — 15-metric concordance MSC vs |ImCoh|.

Produces **two independent figures** (one per fc_method), each replicating
the wp1 ``metric_concordance_15x15_avg.pdf`` recipe (15×15 patient-averaged
Spearman heatmap + dendrogram strip on top, per-figure hierarchical
clustering order — metrics are NOT held in fixed positions across figures).

Outputs:
  data/outputs/figures/section3/fig_J/fig_J1_metric_concordance_MSC.pdf
  data/outputs/figures/section3/fig_J/fig_J2_metric_concordance_ImCoh.pdf
  data/outputs/figures/section3/fig_J/fig_J_concordance_stats.md

Caches:
  data/reports/metric_exploration/task2_full_results.csv       (MSC, existing)
  data/cache/metric_concordance_imcoh/task2_full_results_imcoh_abs.csv (new)
"""
from __future__ import annotations

import argparse
import shutil
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _shared import apply_pub_style, save_fig

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PHASE_LABELS
from lrg_eegfc.config.paths import (
    FIGURES_ROOT,
    CACHE_ROOT,
    REPORTS_ROOT,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.utils.metrics.reorganization import (
    build_metric_specs,
    compute_cluster_labels,
    cluster_swap_coefficient,
)


SECTION3_ROOT = FIGURES_ROOT / "section3"
FIG_J_DIR = SECTION3_ROOT / "fig_J"
WRITING_DIR = SECTION3_ROOT / "for_writing_agent"

MSC_LONG_CSV = REPORTS_ROOT / "metric_exploration" / "task2_full_results.csv"
IMCOH_CONC_CACHE = CACHE_ROOT / "metric_concordance_imcoh"
IMCOH_LONG_CSV = IMCOH_CONC_CACHE / "task2_full_results_imcoh_abs.csv"

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
PHASES = list(PHASE_LABELS)
BANDS = list(BRAIN_BANDS_NAMES)
ALL_PAIRS = list(combinations(PHASES, 2))

FC_METHODS = [
    ("msc",       r"$\mathrm{MSC}$",     "MSC",   MSC_LONG_CSV),
    ("imcoh_abs", r"$|\mathrm{ImCoh}|$", "ImCoh", IMCOH_LONG_CSV),
]

ALL_METRICS = [
    "matrix_distance", "scaled_distance", "rank_distance", "quantile_rmse",
    "permutation_robust",
    "tree_robinson_foulds", "tree_cophenetic_corr", "tree_baker_gamma",
    "tree_fowlkes_mallows",
    "ari", "cluster_swap",
    "fc_frobenius", "fc_scaled_frobenius", "fc_rank_distance", "fc_mean_abs_diff",
]

METRIC_DISPLAY = {
    "matrix_distance":      "Frobenius",
    "scaled_distance":      "Scaled (log)",
    "rank_distance":        "Rank (1−ρ)",
    "quantile_rmse":        "Quantile RMSE",
    "permutation_robust":   "Perm-robust",
    "tree_robinson_foulds": "Robinson–Foulds",
    "tree_cophenetic_corr": "Cophenetic corr.",
    "tree_baker_gamma":     "Baker γ",
    "tree_fowlkes_mallows": "Fowlkes–Mallows",
    "ari":                  "ARI",
    "cluster_swap":         "Cluster swap",
    "fc_frobenius":         "FC Frobenius",
    "fc_scaled_frobenius":  "FC scaled Frob.",
    "fc_rank_distance":     "FC rank (1−ρ)",
    "fc_mean_abs_diff":     "FC mean |Δ|",
}

METRIC_FAMILY = {}
for m in ["matrix_distance", "scaled_distance", "rank_distance",
          "quantile_rmse", "permutation_robust"]:
    METRIC_FAMILY[m] = "Ultrametric"
for m in ["tree_robinson_foulds", "tree_cophenetic_corr",
          "tree_baker_gamma", "tree_fowlkes_mallows"]:
    METRIC_FAMILY[m] = "Tree"
for m in ["ari", "cluster_swap"]:
    METRIC_FAMILY[m] = "Partition"
for m in ["fc_frobenius", "fc_scaled_frobenius",
          "fc_rank_distance", "fc_mean_abs_diff"]:
    METRIC_FAMILY[m] = "Raw FC"

FAMILY_COLORS = {
    "Ultrametric": "#1f77b4",
    "Tree":        "#2ca02c",
    "Partition":   "#ff7f0e",
    "Raw FC":      "#d62728",
}

SIMILARITY_METRICS = {
    "tree_cophenetic_corr", "tree_baker_gamma", "tree_fowlkes_mallows",
    "ari", "cluster_swap",
}


# --- FC-baseline metrics (replicated from scripts/wp0/_common.py) ----------
def _upper(W):
    return W[np.triu_indices_from(W, k=1)]


def fc_frobenius(W1, W2):
    return float(np.linalg.norm(_upper(W1) - _upper(W2)))


def fc_scaled_frobenius(W1, W2):
    v1, v2 = _upper(W1), _upper(W2)
    diff = np.linalg.norm(v1 - v2)
    denom = np.linalg.norm(v1) + np.linalg.norm(v2)
    return float(diff / max(denom, 1e-15))


def fc_rank_distance(W1, W2):
    v1, v2 = _upper(W1), _upper(W2)
    corr, _ = spearmanr(v1, v2)
    return float(1.0 - corr)


def fc_mean_abs_diff(W1, W2):
    return float(np.mean(np.abs(_upper(W1) - _upper(W2))))


FC_BASELINES = {
    "fc_frobenius": fc_frobenius,
    "fc_scaled_frobenius": fc_scaled_frobenius,
    "fc_rank_distance": fc_rank_distance,
    "fc_mean_abs_diff": fc_mean_abs_diff,
}


def compute_long_table(fc_method: str, verbose: bool = False) -> pd.DataFrame:
    """Compute long-form (patient, band, phase_a, phase_b, metric, value) table."""
    if verbose:
        print(f"  Computing 15-metric table for fc_method={fc_method}")
    metric_specs = build_metric_specs()
    rows = []

    for pat in PATIENTS:
        if verbose:
            print(f"    {pat}")
        # Load all (band, phase) data for this patient
        lrg_data, fc_data = {}, {}
        for band in BANDS:
            for phase in PHASES:
                lr = load_lrg_result(pat, phase, band, fc_method)
                if lr is not None:
                    lrg_data[(band, phase)] = lr
                W = load_fc_matrix(pat, phase, band, fc_method)
                if W is not None:
                    fc_data[(band, phase)] = W

        for band in BANDS:
            for p1, p2 in ALL_PAIRS:
                # LRG-derived ultrametric / tree metrics
                lr1 = lrg_data.get((band, p1))
                lr2 = lrg_data.get((band, p2))
                if lr1 is not None and lr2 is not None and lr1.n_nodes == lr2.n_nodes:
                    U1 = squareform(lr1.ultrametric_matrix) \
                        if lr1.ultrametric_matrix.ndim == 1 else lr1.ultrametric_matrix
                    U2 = squareform(lr2.ultrametric_matrix) \
                        if lr2.ultrametric_matrix.ndim == 1 else lr2.ultrametric_matrix
                    Z1, Z2 = lr1.linkage_matrix, lr2.linkage_matrix
                    for mname, mspec in metric_specs.items():
                        try:
                            v = float(mspec["fn"](U1, U2, Z1, Z2))
                        except Exception:
                            v = np.nan
                        rows.append({
                            "patient": pat, "band": band,
                            "phase_a": p1, "phase_b": p2,
                            "metric_name": mname, "value": v,
                        })
                    # ARI + cluster_swap at optimal threshold
                    try:
                        la = compute_cluster_labels(Z1, lr1.optimal_threshold)
                        lb = compute_cluster_labels(Z2, lr2.optimal_threshold)
                        ari = float(adjusted_rand_score(la, lb))
                        swap = float(cluster_swap_coefficient(la, lb))
                    except Exception:
                        ari = swap = np.nan
                    for mname, v in [("ari", ari), ("cluster_swap", swap)]:
                        rows.append({
                            "patient": pat, "band": band,
                            "phase_a": p1, "phase_b": p2,
                            "metric_name": mname, "value": v,
                        })

                # Raw-FC baselines
                W1 = fc_data.get((band, p1))
                W2 = fc_data.get((band, p2))
                if W1 is not None and W2 is not None and W1.shape == W2.shape:
                    for mname, fn in FC_BASELINES.items():
                        try:
                            v = float(fn(W1, W2))
                        except Exception:
                            v = np.nan
                        rows.append({
                            "patient": pat, "band": band,
                            "phase_a": p1, "phase_b": p2,
                            "metric_name": mname, "value": v,
                        })

    return pd.DataFrame(rows)


def load_or_compute_long(fc_method: str, csv_path: Path,
                         verbose: bool = False) -> pd.DataFrame:
    if csv_path.exists():
        if verbose:
            print(f"  [{fc_method}] cached long-form: {csv_path}")
        return pd.read_csv(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df = compute_long_table(fc_method, verbose=verbose)
    df.to_csv(csv_path, index=False)
    if verbose:
        print(f"  [{fc_method}] wrote: {csv_path} ({len(df)} rows)")
    return df


def patient_average_corr(df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Per-patient 15×15 Spearman, then average across patients."""
    # Sign-flip similarity metrics so all are distance-convention
    df = df.copy()
    df["dist_value"] = np.where(
        df["metric_name"].isin(SIMILARITY_METRICS),
        -df["value"], df["value"],
    )

    metrics_present = [m for m in ALL_METRICS
                       if m in df["metric_name"].unique()]
    n = len(metrics_present)

    per_pat = []
    for pat in PATIENTS:
        pdf = df[df["patient"] == pat]
        pivot = pdf.pivot_table(
            index=["band", "phase_a", "phase_b"],
            columns="metric_name", values="dist_value",
        )
        cols = [m for m in metrics_present if m in pivot.columns]
        pivot = pivot[cols]
        corr = np.full((n, n), np.nan)
        for i, mi in enumerate(metrics_present):
            if mi not in pivot.columns:
                continue
            for j, mj in enumerate(metrics_present):
                if mj not in pivot.columns:
                    continue
                vi = pivot[mi].values
                vj = pivot[mj].values
                valid = np.isfinite(vi) & np.isfinite(vj)
                if valid.sum() >= 3:
                    corr[i, j], _ = spearmanr(vi[valid], vj[valid])
        per_pat.append(corr)

    with np.errstate(invalid="ignore"):
        avg = np.nanmean(np.stack(per_pat, axis=0), axis=0)
    return avg, metrics_present


def render_concordance_figure(avg_corr: np.ndarray, metrics: list[str],
                              fc_label_short: str, fc_label_tex: str,
                              out_pdf: Path) -> dict:
    n = len(metrics)

    # Hierarchical clustering on |1 - |ρ||  (per-figure order)
    dist = 1.0 - np.abs(avg_corr)
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)
    dist = np.clip(np.nan_to_num(dist, nan=1.0), 0, None)
    Z = linkage(squareform(dist), method="average")
    leaf_order = leaves_list(Z)

    corr_R = avg_corr[leaf_order][:, leaf_order]
    labels_R = [metrics[i] for i in leaf_order]
    display_R = [METRIC_DISPLAY.get(m, m) for m in labels_R]
    family_R = [METRIC_FAMILY.get(m, "?") for m in labels_R]

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(corr_R, cmap="RdBu_r", vmin=-1, vmax=1,
                   aspect="equal", interpolation="nearest")

    for i in range(n):
        for j in range(n):
            v = corr_R[i, j]
            if np.isfinite(v):
                color = "white" if abs(v) > 0.6 else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=5.5, color=color)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    xls = ax.set_xticklabels(display_R, rotation=45, ha="right", fontsize=8)
    yls = ax.set_yticklabels(display_R, fontsize=8)
    for i, (xl, yl) in enumerate(zip(xls, yls)):
        c = FAMILY_COLORS.get(family_R[i], "black")
        xl.set_color(c)
        yl.set_color(c)

    # Dendrogram strip overhead
    dn = dendrogram(Z, no_plot=True)
    y_top = -0.5
    dend_max_h = max(val for seg in dn["dcoord"] for val in seg)
    dend_scale = 2.5 / max(dend_max_h, 1e-10)
    for xc, yc in zip(dn["icoord"], dn["dcoord"]):
        xs = [(x - 5) / 10 for x in xc]
        ys = [y_top - d * dend_scale for d in yc]
        ax.plot(xs, ys, color="0.3", linewidth=1.0, clip_on=False, zorder=5)
    ax.set_ylim(n - 0.5, y_top - 2.5 - 0.3)

    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02,
                 label=r"Spearman $r_S$")
    handles = [Patch(facecolor=c, label=f) for f, c in FAMILY_COLORS.items()]
    ax.legend(handles=handles, loc="lower left", fontsize=7,
              title="Metric family", title_fontsize=8)

    save_fig(fig, out_pdf)
    return {
        "leaf_order": leaf_order,
        "labels_reordered": labels_R,
        "corr_reordered": corr_R,
    }


def write_stats(per_method_corr: dict, per_method_metrics: dict,
                out_md: Path) -> None:
    lines = []
    lines.append("# Figure J — 15-metric concordance MSC vs |ImCoh|\n")
    lines.append(
        f"Patients: {', '.join(PATIENTS)}. 15 metrics × (6 phase-pairs × "
        "6 bands = 36 obs/metric per patient) → 15×15 Spearman per patient, "
        "averaged across patients. Each panel uses its **own** hierarchical "
        "clustering order (metrics are not pinned).\n"
    )

    # Mean off-diagonal per family pair
    def family_means(corr, metrics):
        out = {}
        for fa in FAMILY_COLORS:
            for fb in FAMILY_COLORS:
                if fa > fb:
                    continue
                vals = []
                for i, mi in enumerate(metrics):
                    for j, mj in enumerate(metrics):
                        if i >= j:
                            continue
                        if METRIC_FAMILY.get(mi) == fa and METRIC_FAMILY.get(mj) == fb:
                            v = corr[i, j]
                            if np.isfinite(v):
                                vals.append(v)
                        elif METRIC_FAMILY.get(mi) == fb and METRIC_FAMILY.get(mj) == fa:
                            v = corr[i, j]
                            if np.isfinite(v):
                                vals.append(v)
                if vals:
                    out[(fa, fb)] = float(np.mean(vals))
        return out

    lines.append("## Mean off-diagonal concordance\n")
    lines.append("| fc_method | mean off-diag $r_S$ |")
    lines.append("|---|---|")
    for fc_method, _, label_short, _ in FC_METHODS:
        C = per_method_corr[fc_method]
        n = C.shape[0]
        iu = np.triu_indices(n, k=1)
        v = C[iu]
        v = v[np.isfinite(v)]
        m = float(np.mean(v)) if v.size else float("nan")
        lines.append(f"| {label_short} | {m:+.3f} |")
    lines.append("")

    for fc_method, _, label_short, _ in FC_METHODS:
        C = per_method_corr[fc_method]
        metrics = per_method_metrics[fc_method]
        fam = family_means(C, metrics)
        lines.append(f"## {label_short} — within / cross-family means\n")
        lines.append("| family A | family B | mean $r_S$ |")
        lines.append("|---|---|---|")
        for (fa, fb), m in sorted(fam.items()):
            lines.append(f"| {fa} | {fb} | {m:+.3f} |")
        lines.append("")

    # Top-10 redundant pairs per fc_method
    for fc_method, _, label_short, _ in FC_METHODS:
        C = per_method_corr[fc_method]
        metrics = per_method_metrics[fc_method]
        n = C.shape[0]
        pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                v = C[i, j]
                if np.isfinite(v):
                    pairs.append((abs(v), v, metrics[i], metrics[j]))
        pairs.sort(reverse=True)
        lines.append(f"## {label_short} — top-10 redundant pairs (by |ρ|)\n")
        lines.append("| metric A | metric B | ρ |")
        lines.append("|---|---|---|")
        for _, v, a, b in pairs[:10]:
            lines.append(
                f"| {METRIC_DISPLAY.get(a,a)} | {METRIC_DISPLAY.get(b,b)} "
                f"| {v:+.2f} |"
            )
        lines.append("")

    lines.append(
        "**Note:** Pat_03 was recorded at 1024 Hz (others 2048 Hz). It is "
        "kept in the average as a documented outlier (per CLAUDE.md "
        "negative-control invariant).\n"
    )

    out_md.write_text("\n".join(lines))
    print(f"  Saved: {out_md}")


def main(verbose: bool = False) -> None:
    apply_pub_style()
    FIG_J_DIR.mkdir(parents=True, exist_ok=True)
    WRITING_DIR.mkdir(parents=True, exist_ok=True)

    per_method_corr, per_method_metrics = {}, {}
    pdf_paths = {}

    for fc_method, label_tex, label_short, csv_path in FC_METHODS:
        df = load_or_compute_long(fc_method, csv_path, verbose=verbose)
        avg, metrics = patient_average_corr(df)
        per_method_corr[fc_method] = avg
        per_method_metrics[fc_method] = metrics

        out_pdf = FIG_J_DIR / (
            "fig_J1_metric_concordance_MSC.pdf" if fc_method == "msc"
            else "fig_J2_metric_concordance_ImCoh.pdf"
        )
        render_concordance_figure(avg, metrics, label_short, label_tex, out_pdf)
        shutil.copy(out_pdf, WRITING_DIR / out_pdf.name)
        pdf_paths[fc_method] = out_pdf

    out_md = FIG_J_DIR / "fig_J_concordance_stats.md"
    write_stats(per_method_corr, per_method_metrics, out_md)
    shutil.copy(out_md, WRITING_DIR / out_md.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(verbose=args.verbose)
