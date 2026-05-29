#!/usr/bin/env python3
"""Audit 35 — D(tau) phase-distance triangle at the LRG layer.

Direct LRG analogue of substrate audit_25. For each (patient, band, phase),
load the LRG ultrametric distance matrix D = 1/rho via load_lrg_result(...)
.ultrametric_matrix (condensed upper-triangle, length N(N-1)/2). Compute, per
phase pair (a, b):

    d_P_LRG(a, b) = 1 - Pearson(triu D^a, triu D^b)
    d_S_LRG(a, b) = 1 - Spearman(triu D^a, triu D^b)
    d_F_LRG(a, b) = ||D^a - D^b||_F / sqrt(||D^a||_F * ||D^b||_F)

Triangle: T_d_LRG(p, b) = d(rest_pre, task_test) - d(task_test, rest_post).
**Sign convention: T_d > 0 = trace (rsPost closer to task than rsPre),
T_d < 0 = anti-trace, T_d = 0 = no trace.** Mirrors substrate's column
schema in ``data/audit/raw_fc_phase_distance/Td_per_patient_per_band.csv``
(F/P/S).

Outputs
-------
``data/audit/lrg_phase_distance/Td_per_patient_per_band.csv``
``data/audit/lrg_phase_distance/cohort_summary.csv``
``data/outputs/figures/section_5_lrg_trace/d_rank_triangle/`` (boxplots)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.workflow.lrg import load_lrg_result

OUT_DIR = ROOT / "data" / "audit" / "lrg_phase_distance"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "d_rank_triangle"


def triu_vec(D_or_vec):
    """Return condensed upper-triangle vector from a square matrix or condensed input."""
    arr = np.asarray(D_or_vec)
    if arr.ndim == 1:
        return arr
    return squareform(arr, checks=False) if arr.shape[0] != arr.shape[1] else arr[np.triu_indices_from(arr, k=1)]


def distance_pair(va, vb):
    rs, _ = spearmanr(va, vb)
    rp, _ = pearsonr(va, vb)
    fro = np.linalg.norm(va - vb) / np.sqrt(np.linalg.norm(va) * np.linalg.norm(vb))
    return 1.0 - rp, 1.0 - rs, fro


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            try:
                D = {}
                for phi in ("rest_pre", "task_test", "rest_post"):
                    res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                    if res is None:
                        raise FileNotFoundError(f"missing LRG cache: {pat} {phi} {band}")
                    D[phi] = triu_vec(res.ultrametric_matrix)
                # Phase-pair distances
                d_pre_tt = distance_pair(D["rest_pre"], D["task_test"])
                d_tt_post = distance_pair(D["task_test"], D["rest_post"])
                d_pre_post = distance_pair(D["rest_pre"], D["rest_post"])
                # Triangle scalar: d(pre,tt) - d(tt,post)
                # T_d > 0 = trace (rsPost closer to task than rsPre).
                T_P, T_S, T_F = (a - b for a, b in zip(d_pre_tt, d_tt_post))
                rows.append({
                    "patient": pat,
                    "band": band,
                    "d_P_pre_tt": d_pre_tt[0], "d_S_pre_tt": d_pre_tt[1], "d_F_pre_tt": d_pre_tt[2],
                    "d_P_tt_post": d_tt_post[0], "d_S_tt_post": d_tt_post[1], "d_F_tt_post": d_tt_post[2],
                    "d_P_pre_post": d_pre_post[0], "d_S_pre_post": d_pre_post[1], "d_F_pre_post": d_pre_post[2],
                    "T_P": T_P, "T_S": T_S, "T_F": T_F,
                })
            except Exception as e:
                print(f"[audit_35] WARN {pat} {band}: {e}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "Td_per_patient_per_band.csv", index=False)

    # Per-band cohort summary
    bands_ord = BRAIN_BANDS_NAMES
    summary_rows = []
    for b in bands_ord:
        sub = df[df["band"] == b]
        for col, dist_name in (("T_S", "d_S"), ("T_P", "d_P"), ("T_F", "d_F")):
            v = sub[col].dropna().values
            n_trace = int((v > 0).sum())
            from scipy.stats import wilcoxon
            try:
                _, p_one = wilcoxon(v, alternative="greater")
            except Exception:
                p_one = float("nan")
            summary_rows.append({
                "band": b,
                "distance": dist_name,
                "n_trace": f"{n_trace}/{len(v)}",
                "n_trace_int": n_trace,
                "n_eligible": len(v),
                "T_median": float(np.median(v)),
                "T_iqr": float(np.subtract(*np.percentile(v, [75, 25]))),
                "wilcoxon_one_sided_p": float(p_one),
            })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "cohort_summary.csv", index=False)
    print(summary.to_string(index=False))

    # Boxplots: 3 distances x 6 bands, grouped
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4), sharey=False)
    for ax, col, name in zip(axes, ("T_S", "T_P", "T_F"), ("$T_d^{d_S}$", "$T_d^{d_P}$", "$T_d^{d_F}$")):
        data = [df[df["band"] == b][col].dropna().values for b in bands_ord]
        bp = ax.boxplot(data, labels=[BRAIN_BAND_TEX_DICT[b] for b in bands_ord],
                        showmeans=True, meanline=True, patch_artist=True,
                        boxprops=dict(facecolor="#e8e8e8"))
        ax.axhline(0, color="0.4", lw=0.7, ls="--")
        # Patient dots
        for i, b in enumerate(bands_ord):
            v = df[df["band"] == b][col].dropna().values
            jitter = (np.random.RandomState(7).rand(len(v)) - 0.5) * 0.12
            ax.scatter(np.full(len(v), i + 1) + jitter, v, s=10, color="#1f77b4", alpha=0.7, zorder=3)
        ax.set_title(name)
        ax.set_ylabel("triangle scalar (positive = trace)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "boxplots_per_band.pdf")
    plt.close(fig)
    print(f"[audit_35] outputs at {OUT_DIR}; figure at {FIG_DIR}")


if __name__ == "__main__":
    main()
