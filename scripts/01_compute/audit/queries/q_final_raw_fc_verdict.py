"""Build the final raw-FC verdict: one table + one figure.

Inputs (existing artefacts):
- data/audit/fc_phase_geometry/cohort_summary.csv   (audit_26: 4-phase median+IQR)
- data/audit/raw_fc_phase_distance/Td_per_patient_per_band.csv (Q4: T_d per pat/band/dist)
- data/audit/raw_fc_phase_distance/per_patient_with_scale.csv  (Q2: d_obs + null_med + ratio)
- data/audit/raw_fc_phase_distance/Pat_06/drift_triangle.csv   (Q3 sample drift R²)

Outputs:
- data/audit/raw_fc_phase_distance/final_verdict_table.csv
- data/audit/raw_fc_phase_distance/final_verdict.pdf
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import DATA_ROOT

OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"
N_COHORT = len(PATIENTS_4PHASE)  # single source of truth — never hardcode

# ---- Load ------------------------------------------------------------
geom = pd.read_csv(OUT_BASE / "cohort_geometry_4phase_summary.csv")
td = pd.read_csv(OUT_BASE / "Td_per_patient_per_band.csv")
# Pull all patients' 4-phase distances
distances_all = []
for p in sorted(OUT_BASE.glob("Pat_*/distance_4phase.csv")):
    sub = pd.read_csv(p)
    distances_all.append(sub)
distances_all = pd.concat(distances_all, ignore_index=True)

# ---- Build verdict table (per band, distance = d_S) ------------------
def get_geom(band, d_label, A, B):
    row = geom[(geom.band == band) & (geom.distance == d_label) &
               (((geom.phase_A == A) & (geom.phase_B == B)) |
                ((geom.phase_A == B) & (geom.phase_B == A)))]
    if row.empty:
        return float("nan"), float("nan"), float("nan")
    r = row.iloc[0]
    return float(r["median"]), float(r["q1"]), float(r["q3"])


def n_positive(td_long, band, d_label):
    """Count patients with T_d > 0 (trace direction; full cohort; Pat_03 included)."""
    sub = td_long[td_long.band == band]
    if d_label not in sub.columns:
        return 0
    n = int((sub[d_label] > 0).sum())
    return n


rows = []
for band in BRAIN_BANDS_NAMES:
    rpre_tt_med, rpre_tt_q1, rpre_tt_q3 = get_geom(band, "S", "rest_pre", "task_test")
    tt_post_med, tt_post_q1, tt_post_q3 = get_geom(band, "S", "task_test", "rest_post")
    rpre_post_med, _, _ = get_geom(band, "S", "rest_pre", "rest_post")
    tl_tt_med, _, _ = get_geom(band, "S", "task_learn", "task_test")  # within-task floor
    n_pers = n_positive(td, band, "S")
    n_pool = int(td[td.band == band].shape[0])
    rows.append({
        "band": band,
        "d_S_RPre_TT": round(rpre_tt_med, 3),
        "d_S_TT_RPost": round(tt_post_med, 3),
        "d_S_RPre_RPost": round(rpre_post_med, 3),
        "d_S_TL_TT_within_task": round(tl_tt_med, 3),
        # T_d = d(RPre, TT) - d(TT, RPost); positive = trace (RPost closer to TT than RPre).
        "T_d_median": round(rpre_tt_med - tt_post_med, 3),
        "n_persist_d_S": f"{n_pers}/{n_pool}",
        "n_persist_d_F": f"{n_positive(td, band, 'F')}/{n_pool}",
        "verdict": (
            "persistence" if (n_pers >= 6 and tt_post_med < rpre_tt_med)
            else ("drift-only" if n_pers <= 3 else "borderline")
        ),
    })
verdict_df = pd.DataFrame(rows)
verdict_df.to_csv(OUT_BASE / "final_verdict_table.csv", index=False)

print(f"=== Final raw-FC verdict (cohort median, d_S, n={N_COHORT}) ===")
print(verdict_df.to_string(index=False))

# ---- Build verdict PDF -----------------------------------------------
PHASES_4 = ("rest_pre", "task_learn", "task_test", "rest_post")
PHASE_SHORT = {"rest_pre": "RPre", "task_learn": "TL",
               "task_test": "TT", "rest_post": "RPost"}

with PdfPages(OUT_BASE / "final_verdict.pdf") as pdf:
    # ---------- Page 1: per-band persistence scatter (d_S) ----------
    fig, axes = plt.subplots(1, 6, figsize=(15, 2.8), squeeze=False)
    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        # x = d(RPre, TT), y = d(TT, RPost) per patient
        sub = distances_all[distances_all.distance == "S"]
        x_df = sub[(sub.band == band) & (sub.phase_A == "rest_pre")
                   & (sub.phase_B == "task_test")][["patient", "d_obs"]]
        y_df = sub[(sub.band == band) & (sub.phase_A == "task_test")
                   & (sub.phase_B == "rest_post")][["patient", "d_obs"]]
        if x_df.empty or y_df.empty:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10); continue
        merged = x_df.merge(y_df, on="patient", suffixes=("_x", "_y"))
        in_pool = merged  # n=10
        outl = merged[merged.patient == "Pat_03"]
        in_pool_excl = merged[merged.patient != "Pat_03"]
        # Identity line
        lo = float(min(merged.d_obs_x.min(), merged.d_obs_y.min())) * 0.9
        hi = float(max(merged.d_obs_x.max(), merged.d_obs_y.max())) * 1.05
        ax.plot([lo, hi], [lo, hi], color="k", lw=0.7, ls="--", alpha=0.6)
        # 9 in-pool blue circles + Pat_03 orange triangle = n=10
        ax.scatter(in_pool_excl.d_obs_x, in_pool_excl.d_obs_y, c="tab:blue",
                   edgecolor="k", linewidth=0.4, s=36, zorder=3)
        # Pat_03 distinct
        if not outl.empty:
            ax.scatter(outl.d_obs_x, outl.d_obs_y, c="tab:orange",
                       marker="^", edgecolor="k", linewidth=0.4, s=42,
                       zorder=4, label="Pat_03 (1024 Hz)")
        n_below = int((in_pool.d_obs_y < in_pool.d_obs_x).sum())
        n_total = in_pool.shape[0]
        ax.text(0.05, 0.95, f"persistence: {n_below}/{n_total}",
                transform=ax.transAxes, fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="gray", alpha=0.85))
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
        ax.set_xlabel(r"$d_S(\mathrm{RPre, TT})$", fontsize=8)
        if c == 0:
            ax.set_ylabel(r"$d_S(\mathrm{TT, RPost})$", fontsize=9)
        ax.tick_params(labelsize=6)
    fig.text(0.99, 0.02, "raw-FC verdict P1 — per-band persistence scatter",
             ha="right", fontsize=6, color="gray")
    fig.tight_layout()
    pdf.savefig(fig, dpi=200); plt.close(fig)

    # ---------- Page 2: cohort cluster bars (4 reference distances) ----------
    bar_pairs = [
        ("task_learn", "task_test", "TL↔TT\n(within-task)", "tab:purple"),
        ("rest_pre", "rest_post", "RPre↔RPost\n(within-rest)", "tab:gray"),
        ("rest_pre", "task_test", "RPre↔TT\n(pre→task)", "tab:red"),
        ("task_test", "rest_post", "TT↔RPost\n(task→post)", "tab:green"),
    ]
    fig, axes = plt.subplots(1, 6, figsize=(15, 3.0), squeeze=False, sharey=True)
    pooled = distances_all[distances_all.distance == "S"]
    # n=10; Pat_03 marked separately as orange triangle in the figures
    for c, band in enumerate(BRAIN_BANDS_NAMES):
        ax = axes[0, c]
        xs = np.arange(len(bar_pairs))
        meds, q1s, q3s, cols = [], [], [], []
        for A, B, lab, col in bar_pairs:
            v = pooled[(pooled.band == band)
                       & (((pooled.phase_A == A) & (pooled.phase_B == B))
                          | ((pooled.phase_A == B) & (pooled.phase_B == A)))]
            arr = v.d_obs.dropna().values
            if arr.size:
                meds.append(float(np.median(arr)))
                q1s.append(float(np.percentile(arr, 25)))
                q3s.append(float(np.percentile(arr, 75)))
            else:
                meds.append(0); q1s.append(0); q3s.append(0)
            cols.append(col)
        meds = np.asarray(meds)
        ax.bar(xs, meds, yerr=[meds - np.asarray(q1s), np.asarray(q3s) - meds],
               color=cols, edgecolor="k", linewidth=0.4, capsize=2, alpha=0.85)
        ax.set_xticks(xs)
        ax.set_xticklabels([p[2] for p in bar_pairs], fontsize=6)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
        if c == 0:
            ax.set_ylabel(r"$d_S$ median (Q1, Q3)" f"\nn = {N_COHORT} patients",
                          fontsize=9)
    fig.text(0.99, 0.02, "raw-FC verdict P2 — within-task vs within-rest vs cross",
             ha="right", fontsize=6, color="gray")
    fig.tight_layout()
    pdf.savefig(fig, dpi=200); plt.close(fig)

print(f"\nwrote {OUT_BASE / 'final_verdict_table.csv'}")
print(f"wrote {OUT_BASE / 'final_verdict.pdf'}")
