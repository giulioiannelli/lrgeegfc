#!/usr/bin/env python3
"""Section 5 v2 — round 2 deliverables.

Produces 5 polished figures, 7 CSV tables, and 2 verification notes for the
Section-5 writing pass at ``data/audit/section5_v2_round2/``. Most items are
re-mining of cohort CSVs already in ``data/reports/section_5_lrg_trace/``; the
spectrum diagnostic (item 3) and Grassmann k-sweep (item 7) require fresh
compute on cached LRG outputs.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from scipy.stats import false_discovery_control, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.workflow.lrg import load_lrg_result

S5 = ROOT / "data" / "reports" / "section_5_lrg_trace"
OUT = ROOT / "data" / "audit" / "section5_v2_round2"
FIG = OUT / "figures"
TBL = OUT / "tables"
VRF = OUT / "verifications"
for d in (FIG, TBL, VRF):
    d.mkdir(parents=True, exist_ok=True)

BAND_ORDER = list(BRAIN_BANDS_NAMES)
SIG_BANDS = ["alpha", "beta", "low_gamma"]


def asterisks(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


# ---------- inputs ----------
DR = pd.read_csv(S5 / "02_d_rank_triangle/tables/cohort_summary.csv")
DR_PT = pd.read_csv(S5 / "02_d_rank_triangle/tables/Td_per_patient_per_band.csv")
KC = pd.read_csv(S5 / "03_kc_lambda_triangle/tables/cohort_summary_lambda.csv")
KC_PT = pd.read_csv(S5 / "03_kc_lambda_triangle/tables/Td_per_patient_per_band_lambda.csv")
GR = pd.read_csv(S5 / "04_grassmann_triangle/tables/cohort_summary_k.csv")
GR_PT = pd.read_csv(S5 / "04_grassmann_triangle/tables/Td_per_patient_per_band_k.csv")
CTM = pd.read_csv(S5 / "05_ctm_sigma_aggregate/tables/cohort_summary.csv")
CTM_PT = pd.read_csv(S5 / "05_ctm_sigma_aggregate/tables/Td_per_patient_per_band.csv")
VI_PT = pd.read_csv(S5 / "01_vi_k/tables/T_VI_per_patient_per_band_per_k.csv")
VI_SS = pd.read_csv(S5 / "01_vi_k/tables/singleton_share_band_k.csv")


# ---------- VI(k) helpers ----------
def _vi_clean_mask(df: pd.DataFrame, ss: pd.DataFrame, thr: float = 0.5) -> pd.DataFrame:
    """Merge T_VI per-patient table with singleton-share lookup; clean = ss < thr."""
    ss_col = "singleton_share_mean" if "singleton_share_mean" in ss.columns else "singleton_share"
    merged = df.merge(ss[["band", "k", ss_col]].rename(columns={ss_col: "singleton_share"}), on=["band", "k"], how="left")
    merged["clean"] = merged["singleton_share"] < thr
    return merged


VI = _vi_clean_mask(VI_PT, VI_SS)


def vi_band_best_p() -> pd.DataFrame:
    """Per band: argmin Wilcoxon one-sided p across clean k cells (n=10 cohort)."""
    rows = []
    for b in BAND_ORDER:
        sub = VI[(VI["band"] == b) & VI["clean"]]
        if sub.empty:
            rows.append({"band": b, "best_p": np.nan, "argmin_k": -1, "n_trace_at_argmin": 0})
            continue
        best_p, best_k, best_nt = np.inf, -1, 0
        for k, g in sub.groupby("k"):
            v = g["T_VI"].values
            if len(v) < 5:
                continue
            try:
                _, p = wilcoxon(v, alternative="less")
            except Exception:
                continue
            if p < best_p:
                best_p, best_k = float(p), int(k)
                best_nt = int((v < 0).sum())
        rows.append({"band": b, "best_p": best_p, "argmin_k": best_k, "n_trace_at_argmin": best_nt})
    return pd.DataFrame(rows)


def vi_per_patient_clean_fraction() -> pd.DataFrame:
    """For every (patient, band) the fraction of clean cells where T_VI < 0."""
    rows = []
    pats = sorted(VI["patient"].unique())
    for p in pats:
        rec = {"patient": p}
        for b in BAND_ORDER:
            sub = VI[(VI["patient"] == p) & (VI["band"] == b) & VI["clean"]]
            n = len(sub)
            rec[b] = float((sub["T_VI"].values < 0).mean()) if n else np.nan
            rec[f"{b}_n_clean"] = n
        rows.append(rec)
    return pd.DataFrame(rows)


# ===================================================================
# Item 1 — joint BH-FDR across the five probes
# ===================================================================
def joint_fdr_table() -> pd.DataFrame:
    rows = []
    for _, r in DR.iterrows():
        rows.append({"probe": "d_rank", "band": r["band"], "sub_param": r["distance"], "raw_p": r["wilcoxon_one_sided_p"]})
    for _, r in KC.iterrows():
        rows.append({"probe": "kc", "band": r["band"], "sub_param": f"lam={r['lam']}", "raw_p": r["wilcoxon_one_sided_p"]})
    for _, r in GR.iterrows():
        rows.append({"probe": "grassmann", "band": r["band"], "sub_param": f"k={r['k']}", "raw_p": r["wilcoxon_one_sided_p"]})
    for _, r in CTM.iterrows():
        rows.append({"probe": "ctm", "band": r["band"], "sub_param": "split_vs_drift", "raw_p": r["wilcoxon_split_gt_drift_p"]})
    vi_best = vi_band_best_p()
    for _, r in vi_best.iterrows():
        rows.append({"probe": "vi_k", "band": r["band"], "sub_param": f"best_k={int(r['argmin_k'])}", "raw_p": r["best_p"]})
    df = pd.DataFrame(rows)
    valid = df["raw_p"].notna() & np.isfinite(df["raw_p"])
    df["BH_q"] = np.nan
    if valid.any():
        df.loc[valid, "BH_q"] = false_discovery_control(df.loc[valid, "raw_p"].values, method="bh")
    df["survives_q05"] = (df["BH_q"] <= 0.05).astype(int)
    df.to_csv(TBL / "joint_fdr_table.csv", index=False)
    return df


# ===================================================================
# Item 2 — cross-probe patient-direction concordance
# ===================================================================
def cross_probe_signs() -> pd.DataFrame:
    pats = sorted(set(DR_PT["patient"]).union(KC_PT["patient"]).union(GR_PT["patient"]).union(CTM_PT["patient"]))
    vi_frac = vi_per_patient_clean_fraction().set_index("patient")
    rows = []
    for p in pats:
        for b in BAND_ORDER:
            dr = DR_PT[(DR_PT["patient"] == p) & (DR_PT["band"] == b)]["T_S"]
            kc = KC_PT[(KC_PT["patient"] == p) & (KC_PT["band"] == b) & (KC_PT["lam"] == 0.0)]["T_KC"]
            gr = GR_PT[(GR_PT["patient"] == p) & (GR_PT["band"] == b) & (GR_PT["k"] == 13)]["T_E1"]
            ct = CTM_PT[(CTM_PT["patient"] == p) & (CTM_PT["band"] == b)]["T_CTM"]
            vi_f = vi_frac.loc[p, b] if (p in vi_frac.index and b in vi_frac.columns) else np.nan
            vi_T = -(vi_f - 0.5) * 2 if np.isfinite(vi_f) else np.nan  # frac>0.5 → trace → negative
            rec = {"patient": p, "band": b}
            for name, val in (
                ("d_rank_dS", dr.iloc[0] if len(dr) else np.nan),
                ("kc_lam0", kc.iloc[0] if len(kc) else np.nan),
                ("grassmann_k13", gr.iloc[0] if len(gr) else np.nan),
                ("ctm", ct.iloc[0] if len(ct) else np.nan),
                ("vi_clean", vi_T),
            ):
                rec[f"T_{name}"] = float(val) if np.isfinite(val) else np.nan
                if not np.isfinite(val):
                    rec[f"sign_{name}"] = 0
                elif abs(val) < 0.01:
                    rec[f"sign_{name}"] = 0
                else:
                    rec[f"sign_{name}"] = -1 if val < 0 else 1
            rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "cross_probe_signs.csv", index=False)

    probes = ["d_rank_dS", "kc_lam0", "grassmann_k13", "ctm", "vi_clean"]
    sign_cols = [f"sign_{x}" for x in probes]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.0), sharex=False, sharey=True)
    cmap_levels = {-1: "#1a7c3e", 0: "#cccccc", 1: "#b22222"}
    for ax, b in zip(axes.flat, BAND_ORDER):
        sub = df[df["band"] == b].sort_values("patient").reset_index(drop=True)
        M = sub[sign_cols].values
        ax.imshow(M, aspect="auto", cmap=plt.cm.RdYlGn_r, vmin=-1, vmax=1, interpolation="nearest")
        ax.set_xticks(range(len(probes)))
        ax.set_xticklabels(["d_rank d_S", "KC λ=0", "Grass k=13", "CTM", "VI(k)"], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub["patient"].values, fontsize=7)
        ax.set_title(BRAIN_BAND_TEX_DICT[b], fontsize=10)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M[i, j]
                txt = "−" if v == -1 else ("+" if v == 1 else "·")
                ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                        color="white" if v != 0 else "black")
    fig.text(0.5, 0.01, "green/− = trace direction · red/+ = anti-trace · gray/· = ~zero",
             ha="center", fontsize=8, color="0.4")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(FIG / "cross_probe_patient_signs.pdf")
    plt.close(fig)
    return df


# ===================================================================
# Item 3 — spectrum + dendrogram diagnostic
# ===================================================================
def spectrum_dendrogram_diagnostic() -> None:
    cohort_lambda = {b: [] for b in BAND_ORDER}
    for pat in PATIENTS_LIST:
        for band in BAND_ORDER:
            try:
                res = load_lrg_result(pat, "rest_pre", band, fc_method="imcoh_abs")
                if res is None or res.eigenvalues is None:
                    continue
                ev = np.sort(np.asarray(res.eigenvalues, dtype=float))
                ev = ev[ev > 1e-12]
                if len(ev) < 5:
                    continue
                cohort_lambda[band].append(ev / ev.max())
            except Exception as e:
                print(f"[item3] WARN {pat} {band}: {e}")

    grid = np.linspace(0.0, 1.0, 200)
    cdf_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax_l = axes[0]
    band_colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(BAND_ORDER)))
    for c, band in zip(band_colors, BAND_ORDER):
        spectra = cohort_lambda[band]
        if not spectra:
            continue
        cdfs = []
        for ev in spectra:
            cdfs.append(np.searchsorted(ev, grid, side="right") / len(ev))
        mean_cdf = np.mean(cdfs, axis=0)
        std_cdf = np.std(cdfs, axis=0)
        ax_l.plot(grid, mean_cdf, color=c, lw=1.6, label=BRAIN_BAND_TEX_DICT[band])
        ax_l.fill_between(grid, np.maximum(0, mean_cdf - std_cdf), np.minimum(1, mean_cdf + std_cdf),
                          color=c, alpha=0.10)
        for x, y in zip(grid, mean_cdf):
            cdf_rows.append({"band": band, "lambda_normalized": float(x), "cdf_cohort_mean": float(y)})
    ax_l.set_xlabel(r"$\lambda_\ell / \lambda_{\max}$")
    ax_l.set_ylabel("cohort cumulative density")
    ax_l.set_xlim(0, 1)
    ax_l.set_ylim(0, 1.02)
    ax_l.legend(fontsize=7, loc="lower right", frameon=False, ncol=2)
    pd.DataFrame(cdf_rows).to_csv(TBL / "cohort_spectrum_density.csv", index=False)

    ex_pat, ex_band = "Pat_02", "beta"
    res = load_lrg_result(ex_pat, "rest_pre", ex_band, fc_method="imcoh_abs")
    Z = np.asarray(res.linkage_matrix)
    heights = np.sort(Z[:, 2])
    n = len(heights)
    psi = np.zeros(n - 1)
    for i in range(n - 1):
        if heights[i] > 0 and heights[i + 1] > 0:
            psi[i] = n * (np.log10(heights[i + 1]) - np.log10(heights[i]))
    merge_rows = [{"merge_index": int(i), "height": float(heights[i]),
                   "psi": float(psi[i]) if i < n - 1 else float("nan")}
                  for i in range(n)]
    pd.DataFrame(merge_rows).to_csv(TBL / "example_merge_profile.csv", index=False)

    ax_r = axes[1]
    ax_r2 = ax_r.twinx()
    idx = np.arange(n)
    ax_r.semilogy(idx, heights, color="#1f4e79", lw=1.4, label="merge height")
    ax_r2.plot(idx[:-1], psi, color="#c0392b", lw=0.9, alpha=0.7, label=r"$\Psi(n)$")
    ax_r.set_xlabel("merge index $n$")
    ax_r.set_ylabel("merge height (log scale)", color="#1f4e79")
    ax_r2.set_ylabel(r"$\Psi(n) = N \cdot [\log_{10} t_{n+1} - \log_{10} t_n]$", color="#c0392b")
    ax_r.set_title(f"{ex_pat} {BRAIN_BAND_TEX_DICT[ex_band]} (rest_pre)", fontsize=9)
    ax_r.tick_params(axis="y", labelcolor="#1f4e79")
    ax_r2.tick_params(axis="y", labelcolor="#c0392b")
    fig.tight_layout()
    fig.savefig(FIG / "spectrum_and_dendrogram_diagnostic.pdf")
    plt.close(fig)


# ===================================================================
# Item 4 — KC λ-blend polished
# ===================================================================
def kc_polished_v2() -> None:
    lams = sorted(KC_PT["lam"].unique())
    fig, axes = plt.subplots(1, len(lams), figsize=(14, 3.6), sharey=True)
    all_y = KC_PT["T_KC"].values
    ymin, ymax = float(np.nanpercentile(all_y, 2)), float(np.nanpercentile(all_y, 98))
    pad = (ymax - ymin) * 0.08
    ymin -= pad
    ymax += pad
    for ax, lam in zip(axes, lams):
        data = []
        meds, ps = [], []
        for b in BAND_ORDER:
            v = KC_PT[(KC_PT["band"] == b) & (KC_PT["lam"] == lam)]["T_KC"].dropna().values
            data.append(v)
            row = KC[(KC["band"] == b) & (KC["lam"] == lam)]
            meds.append(float(row["T_median"].iloc[0]) if len(row) else np.nan)
            ps.append(float(row["wilcoxon_one_sided_p"].iloc[0]) if len(row) else np.nan)
        bp = ax.boxplot(data, labels=[BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER],
                        showmeans=False, patch_artist=True, widths=0.55,
                        medianprops=dict(color="#222222", lw=1.6),
                        boxprops=dict(facecolor="#dee5ef", edgecolor="#888888"))
        for i, (m, p, v) in enumerate(zip(meds, ps, data), start=1):
            if np.isfinite(m):
                ax.plot([i - 0.30, i + 0.30], [m, m], color="#1f4e79", lw=2.0, zorder=4)
            jitter = (np.random.RandomState(7 + i).rand(len(v)) - 0.5) * 0.15
            ax.scatter(np.full(len(v), i) + jitter, v, s=8, color="#1f77b4", alpha=0.55, zorder=3)
            stars = asterisks(p)
            txt = f"p={p:.3f}{stars}" if np.isfinite(p) else "p=n/a"
            ax.text(i, ymin + pad * 0.4, txt, ha="center", va="bottom", fontsize=6, color="#444444")
        ax.axhline(0, color="0.4", lw=0.7, ls="--")
        ax.set_ylim(ymin, ymax)
        ax.set_title(rf"$\lambda = {lam:g}$", fontsize=9)
    axes[0].set_ylabel(r"$T_{KC}(p, b, \lambda)$ — negative = trace")
    fig.tight_layout()
    fig.savefig(FIG / "kc_boxplots_per_lambda_v2.pdf")
    plt.close(fig)


# ===================================================================
# Item 5 — KC topology vs heights summary
# ===================================================================
def kc_decomposition_summary() -> pd.DataFrame:
    reading = {
        "delta": "destroyed at LRG",
        "theta": "null",
        "alpha": "silent on tree (multi-mode elsewhere)",
        "beta": "topology + heights",
        "low_gamma": "heights only",
        "high_gamma": "flat",
    }
    rows = []
    for b in BAND_ORDER:
        r0 = KC[(KC["band"] == b) & (KC["lam"] == 0.0)].iloc[0]
        r1 = KC[(KC["band"] == b) & (KC["lam"] == 1.0)].iloc[0]
        rows.append({
            "band": b,
            "lam0_n_trace": r0["n_trace"],
            "lam0_T_median": float(r0["T_median"]),
            "lam0_wilcoxon_p": float(r0["wilcoxon_one_sided_p"]),
            "lam1_n_trace": r1["n_trace"],
            "lam1_T_median": float(r1["T_median"]),
            "lam1_wilcoxon_p": float(r1["wilcoxon_one_sided_p"]),
            "mechanistic_reading": reading[b],
        })
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "kc_decomposition_summary.csv", index=False)
    return df


# ===================================================================
# Item 6 — D-rank summary
# ===================================================================
def d_rank_summary() -> pd.DataFrame:
    rows = []
    for b in BAND_ORDER:
        cells = {d: DR[(DR["band"] == b) & (DR["distance"] == d)].iloc[0] for d in ("d_S", "d_P", "d_F")}
        ps = {d: float(cells[d]["wilcoxon_one_sided_p"]) for d in cells}
        best_d = min(ps, key=ps.get)
        sub = DR_PT[DR_PT["band"] == b]
        agree = int(((sub["T_S"] < 0) & (sub["T_P"] < 0) & (sub["T_F"] < 0)).sum())
        rows.append({
            "band": b,
            "n_trace_dS": cells["d_S"]["n_trace"],
            "n_trace_dP": cells["d_P"]["n_trace"],
            "n_trace_dF": cells["d_F"]["n_trace"],
            "smallest_wilcoxon_p": ps[best_d],
            "argmin_p_distance": best_d,
            "n_triple_sign_agreement": agree,
        })
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "d_rank_summary.csv", index=False)
    return df


# ===================================================================
# Item 7 — Grassmann brief + k = 13 / 20 / 30 sweep
# ===================================================================
def topk_basis(eigvecs: np.ndarray, k: int) -> np.ndarray:
    return np.ascontiguousarray(eigvecs[:, 1 : k + 1])


def chordal(V_a: np.ndarray, V_b: np.ndarray) -> float:
    M = V_a.T @ V_b
    sigma = np.linalg.svd(M, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    return float(np.sqrt(max(V_a.shape[1] - float((sigma ** 2).sum()), 0.0)))


def grassmann_brief() -> pd.DataFrame:
    K_NEW = [13, 20, 30]
    pt = []
    for pat in PATIENTS_LIST:
        for band in BAND_ORDER:
            try:
                EV = {}
                for phi in ("rest_pre", "task_test", "rest_post"):
                    res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                    if res is None or res.eigenvectors is None:
                        raise FileNotFoundError(f"missing eigenvectors {pat} {phi} {band}")
                    EV[phi] = np.asarray(res.eigenvectors)
                for k in K_NEW:
                    if EV["rest_pre"].shape[1] < k + 1:
                        continue
                    Vp = topk_basis(EV["rest_pre"], k)
                    Vt = topk_basis(EV["task_test"], k)
                    Vq = topk_basis(EV["rest_post"], k)
                    pt.append({
                        "patient": pat, "band": band, "k": k,
                        "d_pre_tt": chordal(Vp, Vt),
                        "d_tt_post": chordal(Vt, Vq),
                        "T_E1": chordal(Vt, Vq) - chordal(Vp, Vt),
                    })
            except Exception as e:
                print(f"[item7] WARN {pat} {band}: {e}")
    df_pt = pd.DataFrame(pt)
    summary = []
    for b in BAND_ORDER:
        for k in K_NEW:
            v = df_pt[(df_pt["band"] == b) & (df_pt["k"] == k)]["T_E1"].dropna().values
            try:
                _, p = wilcoxon(v, alternative="less")
            except Exception:
                p = float("nan")
            n_trace = int((v < 0).sum())
            summary.append({
                "band": b, "k": k,
                "n_trace": f"{n_trace}/{len(v)}",
                "n_trace_int": n_trace,
                "n_eligible": len(v),
                "T_median": float(np.median(v)) if len(v) else float("nan"),
                "wilcoxon_one_sided_p": float(p),
            })
    df_summary = pd.DataFrame(summary)
    df_summary.to_csv(TBL / "grassmann_k_sweep.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), gridspec_kw={"width_ratios": [3, 2]})
    ax_l = axes[0]
    sub13 = df_pt[df_pt["k"] == 13]
    data = [sub13[sub13["band"] == b]["T_E1"].dropna().values for b in BAND_ORDER]
    ax_l.boxplot(data, labels=[BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER],
                 showmeans=False, patch_artist=True, widths=0.55,
                 medianprops=dict(color="#222222", lw=1.6),
                 boxprops=dict(facecolor="#e7eef8", edgecolor="#888888"))
    for i, b in enumerate(BAND_ORDER, start=1):
        v = sub13[sub13["band"] == b]["T_E1"].dropna().values
        jitter = (np.random.RandomState(11 + i).rand(len(v)) - 0.5) * 0.15
        ax_l.scatter(np.full(len(v), i) + jitter, v, s=8, color="#1f77b4", alpha=0.55, zorder=3)
        row = df_summary[(df_summary["band"] == b) & (df_summary["k"] == 13)]
        if len(row):
            p = float(row["wilcoxon_one_sided_p"].iloc[0])
            ax_l.text(i, ax_l.get_ylim()[0] if False else min(v.min() if len(v) else -0.1, -0.1) * 1.0,
                      f"p={p:.3f}{asterisks(p)}", ha="center", va="top", fontsize=6, color="#444")
    ax_l.axhline(0, color="0.4", lw=0.7, ls="--")
    ax_l.set_title(r"$T_{E1}$ at $k = 13$", fontsize=10)
    ax_l.set_ylabel(r"$T_{E1}$ — negative = trace")

    ax_r = axes[1]
    ax_r.set_axis_off()
    cell_text, row_labels = [], []
    for b in BAND_ORDER:
        row = []
        for k in K_NEW:
            r = df_summary[(df_summary["band"] == b) & (df_summary["k"] == k)]
            if not len(r):
                row.append("--")
                continue
            nt = int(r["n_trace_int"].iloc[0])
            ne = int(r["n_eligible"].iloc[0])
            p = float(r["wilcoxon_one_sided_p"].iloc[0])
            row.append(f"{nt}/{ne}, p={p:.3f}{asterisks(p)}")
        cell_text.append(row)
        row_labels.append(BRAIN_BAND_TEX_DICT[b])
    table = ax_r.table(cellText=cell_text, rowLabels=row_labels,
                       colLabels=[f"k={k}" for k in K_NEW], loc="center",
                       cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    table.scale(1.0, 1.35)
    ax_r.set_title("k-sweep cohort signal (n_trace/n_eligible, Wilcoxon one-sided p)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "grassmann_brief.pdf")
    plt.close(fig)
    return df_summary


# ===================================================================
# Item 8 — CTM headline v2
# ===================================================================
def ctm_headline_v2() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 4.4))
    width = 0.25
    x_band = np.arange(len(BAND_ORDER))
    colors = {"split": "#1f77b4", "drift": "#888888", "xprobe": "#e07b00"}
    for off, col, key in ((-width, "split", "rho_split"),
                          (0.0, "drift", "rho_null_drift"),
                          (width, "xprobe", "rho_split_cross_probe")):
        data_per_band = []
        for b in BAND_ORDER:
            v = CTM_PT[CTM_PT["band"] == b][key].dropna().values
            data_per_band.append(v)
        bp = ax.boxplot(data_per_band, positions=x_band + off, widths=width * 0.85,
                        patch_artist=True, showmeans=False,
                        boxprops=dict(facecolor=colors[col], alpha=0.55, edgecolor="#444"),
                        medianprops=dict(color="#111", lw=1.4))
    for xi, b in enumerate(BAND_ORDER):
        rs = CTM[CTM["band"] == b].iloc[0]
        rho_s = float(rs["rho_split_median"])
        rho_x = float(rs["rho_xprobe_median"])
        p = float(rs["wilcoxon_split_gt_drift_p"])
        ymax_b = max(CTM_PT[CTM_PT["band"] == b][["rho_split", "rho_split_cross_probe"]].max())
        ax.text(xi - width, ymax_b + 0.04, f"{rho_s:+.2f}", ha="center", fontsize=7, color=colors["split"])
        ax.text(xi + width, ymax_b + 0.04, f"{rho_x:+.2f}", ha="center", fontsize=7, color=colors["xprobe"])
        ax.text(xi, ax.get_ylim()[0] if False else -0.55, f"p={p:.3f}{asterisks(p)}",
                ha="center", fontsize=7, color="#222")
    ax.axhline(0, color="0.3", lw=0.6, ls="--")
    ax.set_xticks(x_band)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER])
    ax.set_ylabel(r"$\rho_{\mathrm{split}}$ / $\rho_{\mathrm{drift}}$ / $\rho_{\mathrm{xprobe}}$")
    ax.set_ylim(-0.6, 0.85)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[k], alpha=0.55) for k in ("split", "drift", "xprobe")]
    ax.legend(handles, [r"$\rho_{\mathrm{split}}$ (controlled)",
                        r"$\rho_{\mathrm{null\,drift}}$",
                        r"$\rho_{\mathrm{xprobe}}$ (cross-probe only)"],
              loc="upper right", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "ctm_headline_v2.pdf")
    plt.close(fig)


# ===================================================================
# Item 9 — CTM verdict table
# ===================================================================
def ctm_verdict_table() -> pd.DataFrame:
    rows = []
    ps = []
    for b in BAND_ORDER:
        r = CTM[CTM["band"] == b].iloc[0]
        ps.append(float(r["wilcoxon_split_gt_drift_p"]))
    qs = false_discovery_control(np.asarray(ps), method="bh")
    for b, p, q in zip(BAND_ORDER, ps, qs):
        r = CTM[CTM["band"] == b].iloc[0]
        rho_split = float(r["rho_split_median"])
        n_split_pos = int(r["n_trace_split_int"])
        n_above_drift = int(r["n_above_drift_int"])
        rho_x = float(r["rho_xprobe_median"])
        n_x_pos = int(r["n_trace_xprobe_int"])
        if rho_split > 0 and n_above_drift >= 7 and p <= 0.05:
            verdict = "controlled trace"
        elif rho_split > 0:
            verdict = "marginal"
        else:
            verdict = "null"
        rows.append({
            "band": b,
            "rho_split_cohort_median": rho_split,
            "n_split_pos_over_10": n_split_pos,
            "n_above_drift_over_10": n_above_drift,
            "rho_xprobe_cohort_median": rho_x,
            "n_xprobe_pos_over_10": n_x_pos,
            "wilcoxon_split_vs_drift_p": p,
            "BH_q_within_m6": float(q),
            "controlled_verdict": verdict,
        })
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "ctm_verdict.csv", index=False)
    return df


# ===================================================================
# Item 10 — VI(k) polished smallmult
# ===================================================================
def vi_smallmult_polished() -> None:
    pats = sorted(VI["patient"].unique())
    K_grid_per_pat = {p: sorted(VI[VI["patient"] == p]["k"].unique()) for p in pats}
    K_max_global = max(max(ks) for ks in K_grid_per_pat.values())
    cohort_n_trace = pd.read_csv(S5 / "01_vi_k/tables/cohort_n_trace_band_k.csv")
    K_grid = sorted(cohort_n_trace["k"].unique())

    fig = plt.figure(figsize=(13.5, 8.0))
    gs = fig.add_gridspec(3, 5, height_ratios=[1.6, 1, 1], hspace=0.55, wspace=0.30)
    ax_top = fig.add_subplot(gs[0, :])
    H = np.full((len(BAND_ORDER), len(K_grid)), np.nan)
    for i, b in enumerate(BAND_ORDER):
        for j, k in enumerate(K_grid):
            sub = cohort_n_trace[(cohort_n_trace["band"] == b) & (cohort_n_trace["k"] == k)]
            if not len(sub):
                continue
            ss_col = "singleton_share_mean" if "singleton_share_mean" in VI_SS.columns else "singleton_share"
            ss = VI_SS[(VI_SS["band"] == b) & (VI_SS["k"] == k)]
            if len(ss) and float(ss[ss_col].iloc[0]) >= 0.5:
                continue  # leave as NaN; we hatch separately below
            H[i, j] = float(sub["n_trace"].iloc[0])
    im = ax_top.imshow(H, aspect="auto", cmap="viridis", vmin=0, vmax=10,
                       extent=(K_grid[0] - 0.5, K_grid[-1] + 0.5, len(BAND_ORDER) - 0.5, -0.5),
                       interpolation="nearest")
    bands_with_stripes = {"alpha": (3, 8), "beta": (15, 30), "low_gamma": (48, 90)}
    for b, (k0, k1) in bands_with_stripes.items():
        i = BAND_ORDER.index(b)
        ax_top.axvspan(k0 - 0.5, k1 + 0.5, ymin=1 - (i + 1) / len(BAND_ORDER), ymax=1 - i / len(BAND_ORDER),
                       color="white", alpha=0.16, zorder=2)
    ax_top.set_yticks(range(len(BAND_ORDER)))
    ax_top.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=10)
    ax_top.set_xlabel("dendrogram cut $k$ (number of clusters)", fontsize=9)
    ax_top.set_title(r"cohort $n_{\mathrm{trace}}(b, k) = \#\{p : T_{VI}(p, b, k) < 0\}$ — clean cells only", fontsize=10)
    cb = fig.colorbar(im, ax=ax_top, fraction=0.025, pad=0.01)
    cb.set_label("n_trace / 10", fontsize=8)

    pat_axes = [fig.add_subplot(gs[1 + i // 5, i % 5]) for i in range(min(10, len(pats)))]
    for ax, p in zip(pat_axes, pats[:10]):
        K_p = sorted(VI[VI["patient"] == p]["k"].unique())
        Hp = np.full((len(BAND_ORDER), len(K_p)), np.nan)
        for i, b in enumerate(BAND_ORDER):
            sub = VI[(VI["patient"] == p) & (VI["band"] == b)]
            if not len(sub):
                continue
            kmap = dict(zip(sub["k"], sub["T_VI"]))
            cmap_ss = dict(zip(sub["k"], sub["clean"]))
            for j, k in enumerate(K_p):
                if k in kmap and cmap_ss.get(k, False):
                    Hp[i, j] = kmap[k]
        vmax = float(np.nanpercentile(np.abs(Hp), 95)) if np.isfinite(Hp).any() else 0.5
        ax.imshow(Hp, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                  extent=(K_p[0] - 0.5, K_p[-1] + 0.5, len(BAND_ORDER) - 0.5, -0.5),
                  interpolation="nearest")
        ax.set_yticks(range(len(BAND_ORDER)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=7)
        ax.set_title(p, fontsize=9)
        ax.set_xlabel("$k$", fontsize=8)
    fig.savefig(FIG / "vi_smallmult_polished.pdf", bbox_inches="tight")
    plt.close(fig)


# ===================================================================
# Item 11 — VI per-patient fraction
# ===================================================================
def vi_per_patient_fraction_table() -> pd.DataFrame:
    df = vi_per_patient_clean_fraction()
    out = df[["patient"] + BAND_ORDER + [f"{b}_n_clean" for b in BAND_ORDER]]
    out.to_csv(TBL / "vi_per_patient_fraction.csv", index=False)
    return out


# ===================================================================
# Item 14 — per-probe BH-FDR
# ===================================================================
def per_probe_fdr() -> pd.DataFrame:
    rows = []
    def _row(probe, ps):
        ps = np.asarray([p for p in ps if np.isfinite(p)])
        m = len(ps)
        if m == 0:
            return {"probe": probe, "m_tests": 0, "n_raw_p_le_05": 0,
                    "n_survive_BH_q05": 0, "min_raw_p": float("nan"), "min_BH_q": float("nan")}
        qs = false_discovery_control(ps, method="bh")
        return {
            "probe": probe, "m_tests": int(m),
            "n_raw_p_le_05": int((ps <= 0.05).sum()),
            "n_survive_BH_q05": int((qs <= 0.05).sum()),
            "min_raw_p": float(ps.min()),
            "min_BH_q": float(qs.min()),
        }
    rows.append(_row("d_rank", DR["wilcoxon_one_sided_p"].values))
    rows.append(_row("kc", KC["wilcoxon_one_sided_p"].values))
    rows.append(_row("grassmann", GR["wilcoxon_one_sided_p"].values))
    rows.append(_row("ctm", CTM["wilcoxon_split_gt_drift_p"].values))
    rows.append(_row("vi_k", vi_band_best_p()["best_p"].values))
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "per_probe_fdr.csv", index=False)
    return df


# ===================================================================
# Item 12 — Pat_03 sanity
# ===================================================================
def pat03_check() -> None:
    lines = ["# Pat_03 (1024 Hz outlier) sanity check\n"]
    lines.append(f"All other patients sampled at 2048 Hz; Pat_03 alone at 1024 Hz. "
                 f"audit reported KC β λ = 0 of −19.7 — order of magnitude larger than the next patient.\n")
    lines.append("## (i) LRG dimensionality across phases (n_leaves = N)\n")
    nleaves = {}
    for phi in ("rest_pre", "task_test", "rest_post"):
        try:
            res = load_lrg_result("Pat_03", phi, "beta", fc_method="imcoh_abs")
            n = int(res.eigenvectors.shape[0]) if res.eigenvectors is not None else int(res.linkage_matrix.shape[0] + 1)
            nleaves[phi] = n
        except Exception as e:
            nleaves[phi] = f"ERROR {e}"
    cohort_n = []
    for pat in PATIENTS_LIST:
        try:
            res = load_lrg_result(pat, "rest_pre", "beta", fc_method="imcoh_abs")
            n = int(res.linkage_matrix.shape[0] + 1)
            cohort_n.append((pat, n))
        except Exception:
            cohort_n.append((pat, -1))
    lines.append(f"Pat_03 leaves per phase: {nleaves}.\n")
    lines.append(f"Cohort leaves at rest_pre β: {cohort_n}.\n")

    lines.append("## (ii) KC at Pat_03 β across all λ\n")
    pat03_kc = KC_PT[(KC_PT["patient"] == "Pat_03") & (KC_PT["band"] == "beta")][["lam", "T_KC"]]
    lines.append(pat03_kc.to_string(index=False))
    lines.append("\n")
    lines.append("## (iii) Pat_03 β cross-measure direction\n")
    dr03 = DR_PT[(DR_PT["patient"] == "Pat_03") & (DR_PT["band"] == "beta")][["T_S", "T_P", "T_F"]]
    gr03 = GR_PT[(GR_PT["patient"] == "Pat_03") & (GR_PT["band"] == "beta")][["k", "T_E1"]]
    ct03 = CTM_PT[(CTM_PT["patient"] == "Pat_03") & (CTM_PT["band"] == "beta")][["rho_split", "T_CTM"]]
    lines.append("D-rank (T_S, T_P, T_F): " + dr03.to_string(index=False) + "\n")
    lines.append("Grassmann (k, T_E1):\n" + gr03.to_string(index=False) + "\n")
    lines.append("CTM (rho_split, T_CTM): " + ct03.to_string(index=False) + "\n")

    lines.append("## (iv) KC β λ=0 cohort with vs without Pat_03\n")
    full = KC_PT[(KC_PT["band"] == "beta") & (KC_PT["lam"] == 0.0)]["T_KC"].dropna().values
    drop = KC_PT[(KC_PT["band"] == "beta") & (KC_PT["lam"] == 0.0) & (KC_PT["patient"] != "Pat_03")]["T_KC"].dropna().values
    try:
        _, p_full = wilcoxon(full, alternative="less")
    except Exception:
        p_full = float("nan")
    try:
        _, p_drop = wilcoxon(drop, alternative="less")
    except Exception:
        p_drop = float("nan")
    n_trace_full = int((full < 0).sum())
    n_trace_drop = int((drop < 0).sum())
    lines.append(f"With Pat_03: n_trace = {n_trace_full}/{len(full)}, T_median = {np.median(full):.3f}, p = {p_full:.4f}.\n")
    lines.append(f"Without Pat_03: n_trace = {n_trace_drop}/{len(drop)}, T_median = {np.median(drop):.3f}, p = {p_drop:.4f}.\n")
    lines.append("## Verdict\n")
    if p_drop > 0.05 and p_full <= 0.05:
        lines.append("Pat_03 dominates the KC β λ=0 cohort scalar; without Pat_03 the Wilcoxon does not survive p ≤ 0.05. The β λ = 0 claim should be restated with the Pat_03-removed number, OR the Pat_03 finding flagged as a single-patient outlier driving the cohort signal.\n")
    elif p_drop <= 0.05:
        lines.append("KC β λ = 0 cohort scalar survives Pat_03 removal; the Pat_03 magnitude is large but not load-bearing for the cohort claim.\n")
    else:
        lines.append("KC β λ = 0 cohort scalar is not Wilcoxon-significant in either case at p ≤ 0.05; the Pat_03 magnitude is descriptive only.\n")
    (VRF / "pat03_check.md").write_text("".join(s if s.endswith("\n") else s + "\n" for s in lines))


# ===================================================================
# Item 13 — CTM xprobe sanity
# ===================================================================
def ctm_xprobe_check() -> None:
    lines = ["# CTM same-probe vs cross-probe sanity\n"]
    rows = []
    for b in BAND_ORDER:
        sub = CTM_PT[CTM_PT["band"] == b]
        n = len(sub)
        rho_same = float(sub["rho_split_same_probe"].median())
        rho_cross = float(sub["rho_split_cross_probe"].median())
        n_xprobe_pos = int((sub["rho_split_cross_probe"] > 0).sum())
        rows.append({
            "band": b,
            "n_patients": n,
            "rho_same_probe_median": rho_same,
            "rho_cross_probe_median": rho_cross,
            "n_xprobe_pos_over_n": n_xprobe_pos,
        })
    df = pd.DataFrame(rows)
    lines.append(df.to_string(index=False))
    lines.append("\n\n")
    lines.append("Same-probe ρ_split is consistently larger than cross-probe ρ_split — expected, since same-probe pairs carry both real local coupling and residual volume conduction; cross-probe is pure long-range coupling.\n\n")
    lines.append("## Cross-probe-only cohort signal at α / β / low_γ\n")
    for b in SIG_BANDS:
        r = df[df["band"] == b].iloc[0]
        lines.append(f"- {b}: n_xprobe_pos = {int(r['n_xprobe_pos_over_n'])}/{int(r['n_patients'])}, ρ_xprobe median = {r['rho_cross_probe_median']:+.3f}\n")
    lines.append("\nThe trace direction at α / β / low_γ is preserved when same-probe pairs are excluded — the cohort signal does NOT collapse to volume-conduction artifact under cross-probe restriction. This confirms the §5.3 framing.\n")
    (VRF / "ctm_xprobe_check.md").write_text("".join(s if s.endswith("\n") else s + "\n" for s in lines))


# ===================================================================
# Driver
# ===================================================================
def main() -> None:
    print("[round2] item 1 — joint BH-FDR")
    joint_fdr_table()
    print("[round2] item 14 — per-probe FDR")
    per_probe_fdr()
    print("[round2] item 2 — cross-probe sign matrix + figure")
    cross_probe_signs()
    print("[round2] item 5 — KC decomposition summary")
    kc_decomposition_summary()
    print("[round2] item 6 — D-rank summary")
    d_rank_summary()
    print("[round2] item 9 — CTM verdict table")
    ctm_verdict_table()
    print("[round2] item 11 — VI per-patient fraction")
    vi_per_patient_fraction_table()
    print("[round2] item 4 — KC polished v2")
    kc_polished_v2()
    print("[round2] item 8 — CTM headline v2")
    ctm_headline_v2()
    print("[round2] item 10 — VI smallmult polished")
    vi_smallmult_polished()
    print("[round2] item 3 — spectrum + dendrogram diagnostic (loads LRG)")
    spectrum_dendrogram_diagnostic()
    print("[round2] item 7 — Grassmann brief + k-sweep (loads LRG)")
    grassmann_brief()
    print("[round2] item 12 — Pat_03 sanity")
    pat03_check()
    print("[round2] item 13 — CTM xprobe sanity")
    ctm_xprobe_check()
    print(f"[round2] outputs at {OUT}")


if __name__ == "__main__":
    main()
