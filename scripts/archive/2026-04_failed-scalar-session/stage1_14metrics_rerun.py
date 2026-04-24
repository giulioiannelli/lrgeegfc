#!/usr/bin/env python3
"""Stage 1 — rerun the 14 tree-comparison metrics under imcoh_abs.

For each (patient, band, metric) under FC=imcoh_abs, compute the H2a
directional contrast:
    similarity metric:  s = sim(task_test, rest_post) - sim(rest_pre, rest_post)
    distance metric:    s = dist(rest_pre, rest_post) - dist(task_test, rest_post)
(Positive s → rest_post is more similar to task_test than to rest_pre.)

Cohort test per band: 1-sided Wilcoxon ``s > 0`` across 9 patients;
FDR-BH within each metric across 6 bands (m=6); rank-biserial effect
size; 10k-bootstrap 95 % CI of the mean.

Pre-registered pass criterion per metric (all four must hold):
    β AND γ_l: q < 0.05, ≥ 7/9 patients positive;
    θ AND α: q > 0.1, null;
    rb(β) - rb(θ) ≥ 0.5.

Outputs under data/reports/imcoh_vi/:
    stage1_14metrics_raw.csv        per (pat, band, phase_pair, metric)
    stage1_14metrics_contrast.csv   per (pat, band, metric)
    stage1_14metrics_stats.csv      per (metric, band)
    stage1_14metrics_verdict.csv    per metric
    stage1_14metrics_rerun.md       narrative + per-metric tables
    figures/stage1_14metrics_heatmap.pdf
    figures/stage1_14metrics_profiles.pdf
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE,
)
from lrg_eegfc.config.paths import REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result


# --- shared stats helpers (Wilcoxon, FDR-BH, rank-biserial, bootstrap) ---
_shared_path = ROOT / "scripts" / "01_compute" / "_shared.py"
_spec = importlib.util.spec_from_file_location("_shared_stage1", _shared_path)
_shared = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_shared)
wilcoxon_z = _shared.wilcoxon_z
rank_biserial = _shared.rank_biserial
boot_ci_mean = _shared.boot_ci_mean
bh_fdr = _shared.bh_fdr

# --- 14-metric registry imported from the MSC-era sweep script ---
_sweep_path = ROOT / "scripts" / "06_metric_sweep" / "sweep_alternative_metrics.py"
_spec = importlib.util.spec_from_file_location("sweep_alt_14", _sweep_path)
_sweep = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_sweep)
METRICS = _sweep.METRICS  # dict: name -> (fn(d1,d2), is_similarity_flag)


FC_METHOD = "imcoh_abs"
PATIENTS = [p for p in PATIENTS_4PHASE if p != "Pat_14"]
PHASE_PAIRS_NEEDED = [
    ("rest_pre", "rest_post"),
    ("task_test", "rest_post"),
    ("rest_pre", "task_test"),
]
OUT_DIR = REPORTS_ROOT / "imcoh_vi"
FIG_DIR = OUT_DIR / "figures"


def load_all_data() -> dict:
    """Return ``{(patient, phase, band): {Z, U, c_cond, n_nodes}}``."""
    data = {}
    needed_phases = sorted({p for pair in PHASE_PAIRS_NEEDED for p in pair})
    for pat in PATIENTS:
        for band in BRAIN_BANDS_NAMES:
            for phase in needed_phases:
                lrg = load_lrg_result(pat, phase, band, FC_METHOD)
                if lrg is None:
                    print(f"  MISSING: {pat} {phase} {band}")
                    continue
                Z = lrg.linkage_matrix
                U = lrg.ultrametric_matrix
                c_cond = cophenet(Z)
                data[(pat, phase, band)] = {
                    "Z": Z, "U": U, "c_cond": c_cond, "n_nodes": lrg.n_nodes,
                }
    return data


def compute_pairs(data: dict) -> pd.DataFrame:
    rows = []
    for pat in PATIENTS:
        for band in BRAIN_BANDS_NAMES:
            for p1, p2 in PHASE_PAIRS_NEEDED:
                k1, k2 = (pat, p1, band), (pat, p2, band)
                if k1 not in data or k2 not in data:
                    continue
                d1, d2 = data[k1], data[k2]
                row = {"patient": pat, "band": band, "phase1": p1, "phase2": p2}
                for name, (fn, _) in METRICS.items():
                    try:
                        v = fn(d1, d2)
                        row[name] = float(v) if v is not None else np.nan
                    except Exception as e:
                        print(f"  WARN: {name} {pat} {band} {p1}-{p2}: {e}")
                        row[name] = np.nan
                rows.append(row)
    return pd.DataFrame(rows)


def build_contrasts(pairs_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pat in PATIENTS:
        for band in BRAIN_BANDS_NAMES:
            pb = pairs_df[(pairs_df["patient"] == pat) & (pairs_df["band"] == band)]
            pre_post = pb[(pb["phase1"] == "rest_pre") & (pb["phase2"] == "rest_post")]
            tt_post  = pb[(pb["phase1"] == "task_test") & (pb["phase2"] == "rest_post")]
            if pre_post.empty or tt_post.empty:
                continue
            for name, (_, is_sim) in METRICS.items():
                v_pp = float(pre_post[name].values[0])
                v_tp = float(tt_post[name].values[0])
                if not (np.isfinite(v_pp) and np.isfinite(v_tp)):
                    continue
                s = (v_tp - v_pp) if is_sim else (v_pp - v_tp)
                rows.append({
                    "metric": name, "patient": pat, "band": band,
                    "v_pre_post": v_pp, "v_tt_post": v_tp, "s": float(s),
                    "is_similarity": bool(is_sim),
                })
    return pd.DataFrame(rows)


def cohort_stats(contrast_df: pd.DataFrame) -> pd.DataFrame:
    """Per (metric, band): Wilcoxon/p, FDR-BH q within metric, rb, n_pos."""
    rows = []
    for metric in METRICS:
        per_band = []
        for band in BRAIN_BANDS_NAMES:
            mb = contrast_df[(contrast_df["metric"] == metric) &
                             (contrast_df["band"] == band)]
            s = mb["s"].values.astype(float)
            s = s[np.isfinite(s)]
            n = int(len(s))
            n_pos = int((s > 0).sum())
            z, p = wilcoxon_z(s)
            rb = rank_biserial(s)
            mean, lo, hi = boot_ci_mean(s)
            per_band.append({
                "metric": metric, "band": band, "n": n, "n_positive": n_pos,
                "median": float(np.median(s)) if n else np.nan,
                "mean": mean, "ci_lo": lo, "ci_hi": hi,
                "z": z, "p": p, "rb": rb,
            })
        raw_p = [r["p"] if np.isfinite(r["p"]) else 1.0 for r in per_band]
        q_vals = bh_fdr(raw_p)
        for r, q in zip(per_band, q_vals):
            r["q_fdr"] = float(q)
        rows.extend(per_band)
    return pd.DataFrame(rows)


def verdict(stats_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in METRICS:
        md = stats_df[stats_df["metric"] == metric].set_index("band")

        def g(band, col):
            return md.at[band, col] if band in md.index else np.nan

        beta_q, gl_q = g("beta", "q_fdr"), g("low_gamma", "q_fdr")
        beta_n, gl_n = g("beta", "n_positive"), g("low_gamma", "n_positive")
        theta_q, alpha_q = g("theta", "q_fdr"), g("alpha", "q_fdr")
        theta_med, alpha_med = g("theta", "median"), g("alpha", "median")
        rb_beta, rb_theta = g("beta", "rb"), g("theta", "rb")

        crit = {
            "c1_bg_q<0.05": bool(np.isfinite(beta_q) and np.isfinite(gl_q)
                                 and (beta_q < 0.05) and (gl_q < 0.05)),
            "c2_bg_n>=7": bool(np.isfinite(beta_n) and np.isfinite(gl_n)
                               and (beta_n >= 7) and (gl_n >= 7)),
            "c3_ta_null": bool(np.isfinite(theta_q) and np.isfinite(alpha_q)
                               and (theta_q > 0.1) and (alpha_q > 0.1)),
            "c4_rb_sep>=0.5": bool(np.isfinite(rb_beta) and np.isfinite(rb_theta)
                                    and (rb_beta - rb_theta >= 0.5)),
        }
        rows.append({
            "metric": metric,
            "beta_q": beta_q, "gl_q": gl_q, "beta_n": beta_n, "gl_n": gl_n,
            "theta_q": theta_q, "alpha_q": alpha_q,
            "rb_beta": rb_beta, "rb_theta": rb_theta,
            **crit,
            "n_criteria_met": sum(v for v in crit.values()),
        })
    return pd.DataFrame(rows).sort_values(
        ["n_criteria_met", "c1_bg_q<0.05", "c2_bg_n>=7"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def plot_heatmap(stats_df: pd.DataFrame, out: Path) -> None:
    metric_order = list(METRICS.keys())
    fig, ax = plt.subplots(figsize=(10, 8))
    rb_mat = np.full((len(metric_order), len(BRAIN_BANDS_NAMES)), np.nan)
    q_mat = np.full_like(rb_mat, np.nan)
    for i, m in enumerate(metric_order):
        for j, b in enumerate(BRAIN_BANDS_NAMES):
            r = stats_df[(stats_df["metric"] == m) & (stats_df["band"] == b)]
            if not r.empty:
                rb_mat[i, j] = float(r["rb"].values[0])
                q_mat[i, j] = float(r["q_fdr"].values[0])

    im = ax.imshow(rb_mat, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES])
    ax.set_yticks(range(len(metric_order)))
    ax.set_yticklabels(metric_order, fontsize=9)
    for i in range(rb_mat.shape[0]):
        for j in range(rb_mat.shape[1]):
            rb = rb_mat[i, j]
            q = q_mat[i, j]
            if not np.isfinite(rb):
                continue
            txt = f"{rb:+.2f}"
            if np.isfinite(q) and q < 0.05:
                txt += "★"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7,
                    color="white" if abs(rb) > 0.5 else "black")
    fig.colorbar(im, ax=ax, shrink=0.7, label="rank-biserial (H2a sign)")
    ax.set_title(f"Stage 1 — 14 metrics × 6 bands, H2a contrast\n"
                 f"|ImCoh|, n={len(PATIENTS)}; ★ = q<0.05 (FDR-BH within metric, m=6)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_profiles(contrast_df: pd.DataFrame, out: Path) -> None:
    metric_order = list(METRICS.keys())
    n_cols = 4
    n_rows = (len(metric_order) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 3.0 * n_rows),
                              squeeze=False)
    for k, m in enumerate(metric_order):
        ax = axes[k // n_cols, k % n_cols]
        for pat in PATIENTS:
            vals = []
            for b in BRAIN_BANDS_NAMES:
                v = contrast_df[(contrast_df["metric"] == m) &
                                (contrast_df["patient"] == pat) &
                                (contrast_df["band"] == b)]["s"].values
                vals.append(float(v[0]) if len(v) else np.nan)
            ax.plot(range(len(BRAIN_BANDS_NAMES)), vals, "o-", markersize=3,
                    linewidth=1, label=pat, alpha=0.7)
        ax.axhline(0, color="k", lw=0.6, ls="--")
        ax.set_xticks(range(len(BRAIN_BANDS_NAMES)))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                            fontsize=8)
        ax.set_title(m, fontsize=10)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7, ncol=2, loc="best")
    for k in range(len(metric_order), n_rows * n_cols):
        axes[k // n_cols, k % n_cols].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def write_report(stats_df: pd.DataFrame, verdict_df: pd.DataFrame, out: Path) -> None:
    lines = [
        "# Stage 1 — 14-metric rerun under |ImCoh|",
        "",
        f"Date: 2026-04-24. FC: `{FC_METHOD}`. Patients: n={len(PATIENTS)} "
        "(Pat_14 excluded). 6 bands × 14 metrics.",
        "",
        "## Pre-registered pass criterion (per metric; all four must hold)",
        "",
        "1. β AND γ_l at FDR-BH q<0.05 (within metric across 6 bands, m=6).",
        "2. β AND γ_l with ≥ 7/9 patients having H2a contrast s > 0.",
        "3. θ AND α null (q > 0.1).",
        "4. rank-biserial(β) − rank-biserial(θ) ≥ 0.5.",
        "",
        "## Verdict (sorted by criteria met)",
        "",
        verdict_df[["metric", "n_criteria_met",
                    "c1_bg_q<0.05", "c2_bg_n>=7", "c3_ta_null", "c4_rb_sep>=0.5",
                    "beta_q", "gl_q", "beta_n", "gl_n",
                    "theta_q", "alpha_q", "rb_beta", "rb_theta"]]
            .round(4).to_markdown(index=False),
        "",
        "## Per-metric statistics",
        "",
    ]
    for metric in METRICS:
        md = stats_df[stats_df["metric"] == metric].set_index("band") \
                    .reindex(BRAIN_BANDS_NAMES)
        lines.append(f"### {metric}")
        lines.append("")
        cols = ["n", "n_positive", "median", "mean", "ci_lo", "ci_hi",
                "z", "p", "q_fdr", "rb"]
        lines.append(md[cols].round(4).to_markdown())
        lines.append("")
    out.write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] Loading {FC_METHOD} LRG caches for {len(PATIENTS)} patients...")
    data = load_all_data()
    print(f"  Loaded {len(data)} entries "
          f"(expected {len(PATIENTS)}×3×{len(BRAIN_BANDS_NAMES)}="
          f"{len(PATIENTS)*3*len(BRAIN_BANDS_NAMES)})")

    print("[2/5] Computing per-pair metrics…")
    pairs_df = compute_pairs(data)
    raw_csv = OUT_DIR / "stage1_14metrics_raw.csv"
    pairs_df.to_csv(raw_csv, index=False)
    print(f"  {len(pairs_df)} rows → {raw_csv}")

    print("[3/5] Building H2a contrasts…")
    contrast_df = build_contrasts(pairs_df)
    con_csv = OUT_DIR / "stage1_14metrics_contrast.csv"
    contrast_df.to_csv(con_csv, index=False)
    print(f"  {len(contrast_df)} rows → {con_csv}")

    print("[4/5] Cohort stats…")
    stats_df = cohort_stats(contrast_df)
    stats_csv = OUT_DIR / "stage1_14metrics_stats.csv"
    stats_df.to_csv(stats_csv, index=False)
    print(f"  {len(stats_df)} rows → {stats_csv}")
    verdict_df = verdict(stats_df)
    verdict_csv = OUT_DIR / "stage1_14metrics_verdict.csv"
    verdict_df.to_csv(verdict_csv, index=False)
    print(f"  verdict → {verdict_csv}")

    print("[5/5] Figures + report…")
    plot_heatmap(stats_df, FIG_DIR / "stage1_14metrics_heatmap.pdf")
    plot_profiles(contrast_df, FIG_DIR / "stage1_14metrics_profiles.pdf")
    write_report(stats_df, verdict_df, OUT_DIR / "stage1_14metrics_rerun.md")

    print("\n" + "=" * 80)
    print("VERDICT (sorted by criteria met)")
    print("=" * 80)
    print(verdict_df.to_string(index=False))
    print("=" * 80)


if __name__ == "__main__":
    main()
