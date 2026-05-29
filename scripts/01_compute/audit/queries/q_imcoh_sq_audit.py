"""Parallel substrate audit on `imcoh_sq` for sensitivity-check.

Mirrors the `imcoh_abs` audit (audit_26 + the per-band T_d / Wilcoxon /
structural-vs-drift queries) but on the squared imaginary coherence
substrate. Produces the CSVs the writing-bundle figure script needs to
generate the parallel set of plots, into

    data/audit/raw_fc_phase_distance_imcoh_sq/

Outputs:
- per-patient `Pat_NN/distance_4phase.csv`
- cohort `cohort_geometry_4phase_summary.csv`
- cohort `Td_per_patient_per_band.csv`
- cohort `Td_dS_vs_dF_band_contrast.csv`
- cohort `Td_significance_tests.csv`
- cohort `final_verdict_table.csv`

NOT computed (would require re-running the within-rest_pre split-half
null from audit_25, ~10-30 min): `cohort_summary.csv`, `null.npz`,
`per_patient_with_scale.csv`. The bundle figure script handles those
gracefully (figS1 + figS3 only generated for the substrate that has
the null artefacts).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import networkx as nx
import pandas as pd
from scipy.stats import pearsonr, spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import DATA_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr

FC_METHOD = "imcoh_sq"
OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance_imcoh_sq"
OUT_BASE.mkdir(parents=True, exist_ok=True)

PHASES_4 = ("rest_pre", "task_learn", "task_test", "rest_post")
DISTANCES = ("P", "S", "F")
PHASE_PAIRS = (
    ("rest_pre", "task_test"),
    ("rest_pre", "rest_post"),
    ("task_test", "rest_post"),
)


def _common_giant_indices(adj_phases):
    sets = []
    for adj in adj_phases.values():
        G = nx.from_numpy_array(np.abs(adj))
        comps = list(nx.connected_components(G))
        sets.append(set(max(comps, key=len)) if comps else set())
    common = set.intersection(*sets) if sets else set()
    return np.array(sorted(common), dtype=int)


def _triu_offdiag(A):
    return A[np.triu_indices_from(A, k=1)]


def _all_distances(A_a, A_b):
    ut_a = _triu_offdiag(A_a)
    ut_b = _triu_offdiag(A_b)
    d_P = float(1.0 - pearsonr(ut_a, ut_b).statistic)
    d_S = float(1.0 - spearmanr(ut_a, ut_b).statistic)
    d_F = float(np.linalg.norm(A_a - A_b, ord="fro") /
                np.sqrt(np.linalg.norm(A_a, ord="fro")
                        * np.linalg.norm(A_b, ord="fro")))
    return {"P": d_P, "S": d_S, "F": d_F}


# ---------------------------------------------------------------------------
# Per-patient: 6 phase pairs × 3 distances × 6 bands
# ---------------------------------------------------------------------------
print(f"=== {FC_METHOD} audit · n = {len(PATIENTS_4PHASE)} patients ===")
all_dist_long = []
for p in PATIENTS_4PHASE:
    out_dir = OUT_BASE / p
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for band in BRAIN_BANDS_NAMES:
        adj = {}
        ok = True
        for ph in PHASES_4:
            A = load_fc_matrix(p, ph, band, FC_METHOD)
            if A is None:
                ok = False
                break
            adj[ph] = np.asarray(A, dtype=np.float64)
        if not ok:
            print(f"[{p}/{band}] FC missing, skipping")
            continue
        common_idx = _common_giant_indices(adj)
        if common_idx.size < 10:
            print(f"[{p}/{band}] |V*|={common_idx.size}, skipping")
            continue
        A_phase = {ph: adj[ph][np.ix_(common_idx, common_idx)]
                   for ph in PHASES_4}
        for i, phA in enumerate(PHASES_4):
            for j, phB in enumerate(PHASES_4):
                if j <= i:
                    continue
                d = _all_distances(A_phase[phA], A_phase[phB])
                for d_lbl in DISTANCES:
                    rows.append({
                        "patient": p, "band": band, "distance": d_lbl,
                        "phase_A": phA, "phase_B": phB,
                        "d_obs": d[d_lbl],
                        "n_common": int(common_idx.size),
                    })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "distance_4phase.csv", index=False)
    all_dist_long.append(df)
    print(f"[{p}] wrote {out_dir / 'distance_4phase.csv'}")

distances_long = pd.concat(all_dist_long, ignore_index=True)


# ---------------------------------------------------------------------------
# Cohort 4-phase geometry summary
# ---------------------------------------------------------------------------
geom_rows = []
all_pairs = [(p[0], p[1]) for p in
             [(PHASES_4[i], PHASES_4[j])
              for i in range(4) for j in range(i + 1, 4)]]
for band in BRAIN_BANDS_NAMES:
    for d_lbl in DISTANCES:
        for phA, phB in all_pairs:
            sub = distances_long[
                (distances_long.band == band)
                & (distances_long.distance == d_lbl)
                & (distances_long.phase_A == phA)
                & (distances_long.phase_B == phB)
            ]
            if sub.empty:
                continue
            arr = sub.d_obs.values
            geom_rows.append({
                "band": band, "distance": d_lbl,
                "phase_A": phA, "phase_B": phB,
                "n": int(arr.size),
                "median": float(np.median(arr)),
                "q1": float(np.percentile(arr, 25)),
                "q3": float(np.percentile(arr, 75)),
                "mean": float(np.mean(arr)),
                "std": float(np.std(arr, ddof=1) if arr.size > 1 else 0.0),
            })
geom_df = pd.DataFrame(geom_rows)
geom_df.to_csv(OUT_BASE / "cohort_geometry_4phase_summary.csv", index=False)
print(f"wrote {OUT_BASE / 'cohort_geometry_4phase_summary.csv'}")


# ---------------------------------------------------------------------------
# T_d per patient per band per distance
# ---------------------------------------------------------------------------
td_rows = []
for p in PATIENTS_4PHASE:
    for band in BRAIN_BANDS_NAMES:
        td_row = {"patient": p, "band": band}
        for d_lbl in DISTANCES:
            d_pre_tt = distances_long[
                (distances_long.patient == p) & (distances_long.band == band)
                & (distances_long.distance == d_lbl)
                & (distances_long.phase_A == "rest_pre")
                & (distances_long.phase_B == "task_test")
            ]
            d_tt_post = distances_long[
                (distances_long.patient == p) & (distances_long.band == band)
                & (distances_long.distance == d_lbl)
                & (distances_long.phase_A == "task_test")
                & (distances_long.phase_B == "rest_post")
            ]
            if d_pre_tt.empty or d_tt_post.empty:
                td_row[d_lbl] = np.nan
            else:
                td_row[d_lbl] = (float(d_tt_post.iloc[0].d_obs)
                                 - float(d_pre_tt.iloc[0].d_obs))
        td_rows.append(td_row)
td_df = pd.DataFrame(td_rows, columns=["patient", "band", "F", "P", "S"])
td_df.to_csv(OUT_BASE / "Td_per_patient_per_band.csv", index=False)
print(f"wrote {OUT_BASE / 'Td_per_patient_per_band.csv'}")


# ---------------------------------------------------------------------------
# d_S vs d_F structural-vs-drift band contrast
# ---------------------------------------------------------------------------
contrast_rows = []
for band in BRAIN_BANDS_NAMES:
    sub = td_df[td_df.band == band].dropna(subset=["S", "F"])
    if sub.shape[0] < 5:
        continue
    Td_S = sub["S"].to_numpy()
    Td_F = sub["F"].to_numpy()
    rho_SF, _ = spearmanr(Td_S, Td_F)
    r_SF, _ = pearsonr(Td_S, Td_F)
    sign_agree = int((np.sign(Td_S) == np.sign(Td_F)).sum())
    med_abs_S = float(np.median(np.abs(Td_S)))
    med_abs_F = float(np.median(np.abs(Td_F)))
    med_S = float(np.median(Td_S))
    med_F = float(np.median(Td_F))
    if med_abs_F > 2 * med_abs_S and med_abs_S < 0.05:
        flag = "amplitude-drift"
    elif sign_agree >= 7 and rho_SF > 0.5:
        flag = "structural-shift"
    else:
        flag = "ambiguous"
    contrast_rows.append({
        "band": band, "n_patients": int(sub.shape[0]),
        "median_T_d_S": round(med_S, 4),
        "median_T_d_F": round(med_F, 4),
        "median_abs_T_d_S": round(med_abs_S, 4),
        "median_abs_T_d_F": round(med_abs_F, 4),
        "sign_agree_S_F": f"{sign_agree}/{sub.shape[0]}",
        "spearman_rho_S_F": round(rho_SF, 3),
        "pearson_r_S_F": round(r_SF, 3),
        "flag": flag,
    })
contrast_df = pd.DataFrame(contrast_rows)
contrast_df.to_csv(OUT_BASE / "Td_dS_vs_dF_band_contrast.csv", index=False)
print(f"wrote {OUT_BASE / 'Td_dS_vs_dF_band_contrast.csv'}")


# ---------------------------------------------------------------------------
# Wilcoxon signed-rank significance test per (band, distance)
# ---------------------------------------------------------------------------
def _wilcoxon_p(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 5 or np.allclose(arr, 0):
        return float("nan")
    try:
        stat, p = wilcoxon(arr, alternative="greater")
        return float(p)
    except ValueError:
        return float("nan")


def _sign_p(values):
    """One-sided binomial sign test, alternative T_d > 0."""
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 5:
        return float("nan")
    n = arr.size
    n_pos = int((arr > 0).sum())
    # P(X >= n_pos | p=0.5) under one-sided "greater"
    from math import comb
    p = sum(comb(n, k) for k in range(n_pos, n + 1)) / 2 ** n
    return float(p)


sig_rows = []
for band in BRAIN_BANDS_NAMES:
    for d_lbl in DISTANCES:
        sub = td_df[td_df.band == band]
        vals = sub[d_lbl].dropna().values
        if vals.size < 5:
            continue
        p_w = _wilcoxon_p(vals)
        p_s = _sign_p(vals)
        sub_no_p06 = sub[sub.patient != "Pat_06"]
        vals_no = sub_no_p06[d_lbl].dropna().values
        p_w_no = _wilcoxon_p(vals_no)
        sig_rows.append({
            "band": band, "distance": d_lbl,
            "n": int(vals.size),
            "n_neg": int((vals < 0).sum()),
            "median": round(float(np.median(vals)), 4),
            "wilcoxon_p": round(p_w, 4) if np.isfinite(p_w) else np.nan,
            "sign_p": round(p_s, 4) if np.isfinite(p_s) else np.nan,
            "wilcoxon_p_no_pat06": round(p_w_no, 4) if np.isfinite(p_w_no) else np.nan,
        })
sig_df = pd.DataFrame(sig_rows)
# BH-FDR over the 18 (band x distance) cells
ps_w = sig_df["wilcoxon_p"].fillna(1.0).tolist()
ps_s = sig_df["sign_p"].fillna(1.0).tolist()
sig_df["wilcoxon_q_BH"] = np.round(bh_fdr(ps_w), 4)
sig_df["sign_q_BH"] = np.round(bh_fdr(ps_s), 4)
sig_df.to_csv(OUT_BASE / "Td_significance_tests.csv", index=False)
print(f"wrote {OUT_BASE / 'Td_significance_tests.csv'}")


# ---------------------------------------------------------------------------
# Final verdict table (cohort medians, sign counts, verdict per band)
# ---------------------------------------------------------------------------
def _med_pair(band, A, B):
    sub = distances_long[
        (distances_long.distance == "S") & (distances_long.band == band)
        & (((distances_long.phase_A == A) & (distances_long.phase_B == B))
           | ((distances_long.phase_A == B) & (distances_long.phase_B == A)))
    ]
    return float(np.median(sub.d_obs)) if not sub.empty else float("nan")


verdict_rows = []
for band in BRAIN_BANDS_NAMES:
    pre_tt = _med_pair(band, "rest_pre", "task_test")
    tt_post = _med_pair(band, "task_test", "rest_post")
    pre_post = _med_pair(band, "rest_pre", "rest_post")
    tl_tt = _med_pair(band, "task_learn", "task_test")
    sub_band = td_df[td_df.band == band]
    n_S = int((sub_band["S"].dropna() < 0).sum())
    n_P = int((sub_band["P"].dropna() < 0).sum())
    n_F = int((sub_band["F"].dropna() < 0).sum())
    n_pool = int(sub_band.shape[0])
    if n_S >= 6 and tt_post < pre_tt:
        verdict = "trace"
    elif n_S <= 3:
        verdict = "drift-only"
    else:
        verdict = "borderline"
    verdict_rows.append({
        "band": band,
        "d_S_RPre_TT": round(pre_tt, 3),
        "d_S_TT_RPost": round(tt_post, 3),
        "d_S_RPre_RPost": round(pre_post, 3),
        "d_S_TL_TT_within_task": round(tl_tt, 3),
        "T_d_median": round(tt_post - pre_tt, 3),
        "n_persist_d_S": f"{n_S}/{n_pool}",
        "n_persist_d_P": f"{n_P}/{n_pool}",
        "n_persist_d_F": f"{n_F}/{n_pool}",
        "verdict": verdict,
    })
verdict_df = pd.DataFrame(verdict_rows)
verdict_df.to_csv(OUT_BASE / "final_verdict_table.csv", index=False)
print(f"wrote {OUT_BASE / 'final_verdict_table.csv'}")

print()
print(f"=== {FC_METHOD} verdict (cohort median, n = {len(PATIENTS_4PHASE)}) ===")
print(verdict_df.to_string(index=False))
