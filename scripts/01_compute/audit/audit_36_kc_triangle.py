#!/usr/bin/env python3
"""Audit 36 — KC (Kendall-Colijn) lambda-blend tree-distance triangle.

For each (patient, band, phase) load Z = LRGResult.linkage_matrix (already
cached). Compute, per phase pair (a, b) and per lambda in {0, 0.25, 0.5,
0.75, 1.0}:

    d_KC(lam; T^a, T^b)  via lrg_eegfc.utils.metrics.tree_distance.kc_distance

Triangle:

    T_KC(p, b, lam) = d_KC(lam; T^TT, T^RPost) - d_KC(lam; T^RPre, T^TT)

Negative = trace direction. lambda = 0 is topology-only; lambda = 1 is
height-only; intermediate values blend.

Outputs
-------
``data/audit/kc_triangle/Td_per_patient_per_band_lambda.csv``
``data/audit/kc_triangle/cohort_summary_lambda.csv``
``data/outputs/figures/section_5_lrg_trace/kc_triangle/boxplots_per_lambda.pdf``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.utils.metrics.tree_distance import kc_distance
from lrg_eegfc.workflow.lrg import load_lrg_result

LAMBDAS = [0.0, 0.25, 0.5, 0.75, 1.0]
OUT_DIR = ROOT / "data" / "audit" / "kc_triangle"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "kc_triangle"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            try:
                Z = {}
                for phi in ("rest_pre", "task_test", "rest_post"):
                    res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                    if res is None:
                        raise FileNotFoundError(f"missing LRG cache: {pat} {phi} {band}")
                    Z[phi] = res.linkage_matrix
                for lam in LAMBDAS:
                    d_pre_tt = kc_distance(Z["rest_pre"], Z["task_test"], lam=lam)
                    d_tt_post = kc_distance(Z["task_test"], Z["rest_post"], lam=lam)
                    rows.append({
                        "patient": pat, "band": band, "lam": lam,
                        "d_pre_tt": d_pre_tt, "d_tt_post": d_tt_post,
                        "T_KC": d_tt_post - d_pre_tt,
                    })
            except Exception as e:
                print(f"[audit_36] WARN {pat} {band}: {e}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "Td_per_patient_per_band_lambda.csv", index=False)

    # Cohort summary
    summary_rows = []
    for b in BRAIN_BANDS_NAMES:
        for lam in LAMBDAS:
            sub = df[(df["band"] == b) & (df["lam"] == lam)]
            v = sub["T_KC"].dropna().values
            n_trace = int((v < 0).sum())
            try:
                _, p = wilcoxon(v, alternative="less")
            except Exception:
                p = float("nan")
            summary_rows.append({
                "band": b, "lam": lam,
                "n_trace": f"{n_trace}/{len(v)}",
                "n_trace_int": n_trace,
                "n_eligible": len(v),
                "T_median": float(np.median(v)),
                "T_iqr": float(np.subtract(*np.percentile(v, [75, 25]))),
                "wilcoxon_one_sided_p": float(p),
            })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "cohort_summary_lambda.csv", index=False)
    print(summary.to_string(index=False))

    # Figure: 5-panel boxplots, one per lambda, bands on x
    fig, axes = plt.subplots(1, len(LAMBDAS), figsize=(15, 3.4), sharey=False)
    for ax, lam in zip(axes, LAMBDAS):
        data = [df[(df["band"] == b) & (df["lam"] == lam)]["T_KC"].dropna().values
                for b in BRAIN_BANDS_NAMES]
        ax.boxplot(data, labels=[BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                   showmeans=True, meanline=True, patch_artist=True,
                   boxprops=dict(facecolor="#e8e8e8"))
        for i, b in enumerate(BRAIN_BANDS_NAMES):
            v = df[(df["band"] == b) & (df["lam"] == lam)]["T_KC"].dropna().values
            jitter = (np.random.RandomState(7).rand(len(v)) - 0.5) * 0.12
            ax.scatter(np.full(len(v), i + 1) + jitter, v, s=10, color="#1f77b4", alpha=0.7, zorder=3)
        ax.axhline(0, color="0.4", lw=0.7, ls="--")
        kind = {0.0: "topology", 0.25: "0.25", 0.5: "balanced", 0.75: "0.75", 1.0: "heights"}[lam]
        ax.set_title(rf"$\lambda = {lam}$ ({kind})")
    axes[0].set_ylabel(r"$T_{KC}(p, b, \lambda)$ -- negative = trace")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "boxplots_per_lambda.pdf")
    plt.close(fig)
    print(f"[audit_36] outputs at {OUT_DIR}; figure at {FIG_DIR}")


if __name__ == "__main__":
    main()
