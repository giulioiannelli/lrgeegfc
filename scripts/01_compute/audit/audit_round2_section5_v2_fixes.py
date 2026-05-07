#!/usr/bin/env python3
"""Round 2 figure fixes: VI smallmult, KC v2, Grassmann brief, CTM headline v2,
plus a per-(patient, band, phase) dendrogram + Ψ folder replacing the single-
example right panel of the spectrum diagnostic.

All figures re-render from the same cached CSVs / LRG cache. No new compute
beyond the per-(p, b, phi) dendrogram folder which loads existing LRG outputs.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.workflow.lrg import load_lrg_result

S5 = ROOT / "data" / "reports" / "section_5_lrg_trace"
OUT = ROOT / "data" / "audit" / "section5_v2_round2"
FIG = OUT / "figures"
TBL = OUT / "tables"
DEND_DIR = FIG / "dendrogram_psi"
DEND_DIR.mkdir(parents=True, exist_ok=True)

BAND_ORDER = list(BRAIN_BANDS_NAMES)


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
KC_PT = pd.read_csv(S5 / "03_kc_lambda_triangle/tables/Td_per_patient_per_band_lambda.csv")
KC = pd.read_csv(S5 / "03_kc_lambda_triangle/tables/cohort_summary_lambda.csv")
GR_SWEEP = pd.read_csv(TBL / "grassmann_k_sweep.csv")
CTM_PT = pd.read_csv(S5 / "05_ctm_sigma_aggregate/tables/Td_per_patient_per_band.csv")
CTM = pd.read_csv(S5 / "05_ctm_sigma_aggregate/tables/cohort_summary.csv")
VI_PT = pd.read_csv(S5 / "01_vi_k/tables/T_VI_per_patient_per_band_per_k.csv")
VI_SS = pd.read_csv(S5 / "01_vi_k/tables/singleton_share_band_k.csv")
VI_COHORT = pd.read_csv(S5 / "01_vi_k/tables/cohort_n_trace_band_k.csv")


# ===================================================================
# Item 1 — VI smallmult, trim K_grid to clean range
# ===================================================================
def vi_smallmult_polished_fix() -> None:
    ss_col = "singleton_share_mean" if "singleton_share_mean" in VI_SS.columns else "singleton_share"
    VI_SS_clean = VI_SS.merge(VI_COHORT[["band", "k"]].drop_duplicates(), on=["band", "k"], how="right")
    cohort_clean = VI_COHORT.merge(VI_SS[["band", "k", ss_col]], on=["band", "k"], how="left")
    cohort_clean["clean"] = cohort_clean[ss_col] < 0.5
    k_max_clean_per_band = cohort_clean[cohort_clean["clean"]].groupby("band")["k"].max()
    K_max_global_clean = int(k_max_clean_per_band.max())

    K_grid = sorted(VI_COHORT[VI_COHORT["k"] <= K_max_global_clean]["k"].unique())

    fig = plt.figure(figsize=(13.5, 8.0))
    gs = fig.add_gridspec(3, 5, height_ratios=[1.6, 1, 1], hspace=0.55, wspace=0.30)
    ax_top = fig.add_subplot(gs[0, :])
    H = np.full((len(BAND_ORDER), len(K_grid)), np.nan)
    for i, b in enumerate(BAND_ORDER):
        for j, k in enumerate(K_grid):
            sub = VI_COHORT[(VI_COHORT["band"] == b) & (VI_COHORT["k"] == k)]
            if not len(sub):
                continue
            ss = VI_SS[(VI_SS["band"] == b) & (VI_SS["k"] == k)]
            if len(ss) and float(ss[ss_col].iloc[0]) >= 0.5:
                continue
            H[i, j] = float(sub["n_trace"].iloc[0])
    im = ax_top.imshow(
        H, aspect="auto", cmap="inferno", vmin=0, vmax=10,
        extent=(K_grid[0] - 0.5, K_grid[-1] + 0.5, len(BAND_ORDER) - 0.5, -0.5),
        interpolation="nearest",
    )
    ax_top.set_xlim(K_grid[0] - 0.5, K_grid[-1] + 0.5)
    bands_with_stripes = {"alpha": (3, 8), "beta": (15, 30), "low_gamma": (48, 90)}
    for b, (k0, k1) in bands_with_stripes.items():
        if b not in BAND_ORDER:
            continue
        i = BAND_ORDER.index(b)
        ax_top.add_patch(mpatches.Rectangle(
            (k0 - 0.5, i - 0.5), (k1 - k0 + 1), 1,
            facecolor="white", alpha=0.18, edgecolor="none", zorder=2,
        ))
    ax_top.set_yticks(range(len(BAND_ORDER)))
    ax_top.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=10)
    ax_top.set_xlabel("dendrogram cut $k$ (number of clusters)", fontsize=9)
    ax_top.set_title(r"cohort $n_{\mathrm{trace}}(b, k) = \#\{p : T_{VI}(p, b, k) < 0\}$ — clean cells only",
                     fontsize=10)
    cb = fig.colorbar(im, ax=ax_top, fraction=0.025, pad=0.01)
    cb.set_label("n_trace / 10", fontsize=8)

    pats = sorted(VI_PT["patient"].unique())[:10]
    pat_axes = [fig.add_subplot(gs[1 + i // 5, i % 5]) for i in range(len(pats))]
    for ax, p in zip(pat_axes, pats):
        sub_all = VI_PT[VI_PT["patient"] == p].merge(VI_SS[["band", "k", ss_col]], on=["band", "k"], how="left")
        sub_all["clean"] = sub_all[ss_col] < 0.5
        sub_all = sub_all[sub_all["k"] <= K_max_global_clean]
        K_p = sorted(sub_all["k"].unique())
        if not K_p:
            ax.set_axis_off()
            continue
        Hp = np.full((len(BAND_ORDER), len(K_p)), np.nan)
        for i, b in enumerate(BAND_ORDER):
            band_sub = sub_all[sub_all["band"] == b]
            if not len(band_sub):
                continue
            kmap = dict(zip(band_sub["k"], band_sub["T_VI"]))
            cmap_ss = dict(zip(band_sub["k"], band_sub["clean"]))
            for j, k in enumerate(K_p):
                if k in kmap and bool(cmap_ss.get(k, False)):
                    Hp[i, j] = kmap[k]
        vmax = float(np.nanpercentile(np.abs(Hp), 95)) if np.isfinite(Hp).any() else 0.5
        ax.imshow(
            Hp, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
            extent=(K_p[0] - 0.5, K_p[-1] + 0.5, len(BAND_ORDER) - 0.5, -0.5),
            interpolation="nearest",
        )
        ax.set_xlim(K_p[0] - 0.5, K_p[-1] + 0.5)
        ax.set_yticks(range(len(BAND_ORDER)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=7)
        ax.set_title(p, fontsize=9)
        ax.set_xlabel("$k$", fontsize=8)
    fig.savefig(FIG / "vi_smallmult_polished.pdf", bbox_inches="tight")
    plt.close(fig)


# ===================================================================
# Item 2 — KC polished v2 (showfliers off, legend, fix annotation overlap)
# ===================================================================
def kc_polished_v2_fix() -> None:
    lams = sorted(KC_PT["lam"].unique())
    fig, axes = plt.subplots(1, len(lams), figsize=(14.5, 4.0), sharey=True)
    all_y = KC_PT["T_KC"].values
    p2 = float(np.nanpercentile(all_y, 2))
    p98 = float(np.nanpercentile(all_y, 98))
    pad = (p98 - p2) * 0.18
    ymin, ymax = p2 - pad * 1.4, p98 + pad * 0.6
    for ax, lam in zip(axes, lams):
        data, meds, ps = [], [], []
        for b in BAND_ORDER:
            v = KC_PT[(KC_PT["band"] == b) & (KC_PT["lam"] == lam)]["T_KC"].dropna().values
            data.append(v)
            row = KC[(KC["band"] == b) & (KC["lam"] == lam)]
            meds.append(float(row["T_median"].iloc[0]) if len(row) else np.nan)
            ps.append(float(row["wilcoxon_one_sided_p"].iloc[0]) if len(row) else np.nan)
        ax.boxplot(
            data, labels=[BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER],
            showmeans=False, patch_artist=True, widths=0.55, showfliers=False,
            medianprops=dict(color="#222222", lw=1.6),
            boxprops=dict(facecolor="#dee5ef", edgecolor="#888888"),
        )
        for i, (m, p, v) in enumerate(zip(meds, ps, data), start=1):
            if np.isfinite(m):
                ax.plot([i - 0.30, i + 0.30], [m, m], color="#1f4e79", lw=2.4, zorder=4)
            jitter = (np.random.RandomState(7 + i).rand(len(v)) - 0.5) * 0.16
            ax.scatter(np.full(len(v), i) + jitter, v, s=12, color="#1f77b4",
                       edgecolor="white", linewidth=0.4, alpha=0.75, zorder=3)
            stars = asterisks(p)
            txt = f"{p:.3f}{stars}" if np.isfinite(p) else "n/a"
            color = "#c0392b" if np.isfinite(p) and p <= 0.05 else "#555"
            ax.text(i, 0.985, txt, ha="center", va="top", fontsize=7, color=color,
                    rotation=45, transform=ax.get_xaxis_transform())
        ax.axhline(0, color="0.4", lw=0.7, ls="--")
        ax.set_ylim(ymin, ymax)
        ax.set_title(rf"$\lambda = {lam:g}$", fontsize=10, pad=14)
        ax.tick_params(axis="x", labelsize=9)
    axes[0].set_ylabel(r"$T_{KC}(p, b, \lambda)$ — negative = trace")
    handles = [
        mpatches.Patch(facecolor="#dee5ef", edgecolor="#888888", label="cohort IQR (box)"),
        mlines.Line2D([], [], color="#1f4e79", lw=2.4, label="cohort median"),
        mlines.Line2D([], [], color="#222222", lw=1.6, label="boxplot median"),
        mlines.Line2D([], [], color="#1f77b4", marker="o", linestyle="None", markersize=5,
                      markeredgecolor="white", label="patient (n=10)"),
    ]
    fig.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
        ncol=4, frameon=False, fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIG / "kc_boxplots_per_lambda_v2.pdf", bbox_inches="tight")
    plt.close(fig)


# ===================================================================
# Item 3 — Grassmann brief (showfliers off, legend, fix overlap)
# ===================================================================
def grassmann_brief_fix() -> None:
    K_NEW = [13, 20, 30]
    pt_path = TBL / "grassmann_k_sweep.csv"
    df_summary = pd.read_csv(pt_path)
    GR_E1 = pd.read_csv(S5 / "04_grassmann_triangle/tables/Td_per_patient_per_band_k.csv")
    sub13 = GR_E1[GR_E1["k"] == 13]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [3, 2]})
    ax_l = axes[0]
    data = [sub13[sub13["band"] == b]["T_E1"].dropna().values for b in BAND_ORDER]
    ax_l.boxplot(
        data, labels=[BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER],
        showmeans=False, patch_artist=True, widths=0.55, showfliers=False,
        medianprops=dict(color="#222222", lw=1.6),
        boxprops=dict(facecolor="#e7eef8", edgecolor="#888888"),
    )
    for i, b in enumerate(BAND_ORDER, start=1):
        v = sub13[sub13["band"] == b]["T_E1"].dropna().values
        jitter = (np.random.RandomState(11 + i).rand(len(v)) - 0.5) * 0.16
        ax_l.scatter(np.full(len(v), i) + jitter, v, s=12, color="#1f77b4",
                     edgecolor="white", linewidth=0.4, alpha=0.75, zorder=3)
        row = df_summary[(df_summary["band"] == b) & (df_summary["k"] == 13)]
        if len(row):
            p = float(row["wilcoxon_one_sided_p"].iloc[0])
            color = "#c0392b" if p <= 0.05 else "#222"
            ax_l.text(i, 0.06, f"p={p:.3f}{asterisks(p)}", ha="center", va="bottom",
                      fontsize=7.5, color=color, transform=ax_l.get_xaxis_transform())
    ax_l.axhline(0, color="0.4", lw=0.7, ls="--")
    ax_l.set_title(r"$T_{E1}$ at $k = 13$", fontsize=10)
    ax_l.set_ylabel(r"$T_{E1}$ — negative = trace")
    handles = [
        mpatches.Patch(facecolor="#e7eef8", edgecolor="#888888", label="cohort IQR"),
        mlines.Line2D([], [], color="#222222", lw=1.6, label="median"),
        mlines.Line2D([], [], color="#1f77b4", marker="o", linestyle="None", markersize=5,
                      markeredgecolor="white", label="patient (n=10)"),
    ]
    ax_l.legend(handles=handles, loc="upper right", frameon=False, fontsize=7)

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
    table = ax_r.table(
        cellText=cell_text, rowLabels=row_labels,
        colLabels=[f"k={k}" for k in K_NEW], loc="center", cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.45)
    ax_r.set_title("k-sweep cohort signal\n(n_trace/n_eligible, Wilcoxon one-sided p)",
                   fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "grassmann_brief.pdf", bbox_inches="tight")
    plt.close(fig)


# ===================================================================
# Item 4 — CTM headline v2 (showfliers off, all 3 medians, fix overlap)
# ===================================================================
def ctm_headline_v2_fix() -> None:
    fig, ax = plt.subplots(figsize=(12, 5.0))
    width = 0.25
    x_band = np.arange(len(BAND_ORDER))
    colors = {"split": "#1f77b4", "drift": "#888888", "xprobe": "#e07b00"}
    box_keys = (("split", "rho_split"), ("drift", "rho_null_drift"), ("xprobe", "rho_split_cross_probe"))
    for off, (col_name, key) in zip((-width, 0.0, width), box_keys):
        data_per_band = [CTM_PT[CTM_PT["band"] == b][key].dropna().values for b in BAND_ORDER]
        ax.boxplot(
            data_per_band, positions=x_band + off, widths=width * 0.85,
            patch_artist=True, showmeans=False, showfliers=False,
            boxprops=dict(facecolor=colors[col_name], alpha=0.55, edgecolor="#444"),
            medianprops=dict(color="#111", lw=1.4),
        )
    medians = {col_name: [] for col_name, _ in box_keys}
    for xi, b in enumerate(BAND_ORDER):
        rs = CTM[CTM["band"] == b].iloc[0]
        rho_s = float(rs["rho_split_median"])
        rho_x = float(rs["rho_xprobe_median"])
        rho_d = float(CTM_PT[CTM_PT["band"] == b]["rho_null_drift"].median())
        medians["split"].append(rho_s)
        medians["drift"].append(rho_d)
        medians["xprobe"].append(rho_x)
    ymax_top = 0.95
    for xi in range(len(BAND_ORDER)):
        ax.text(xi - width, ymax_top, f"{medians['split'][xi]:+.2f}",
                ha="center", va="bottom", fontsize=7, color=colors["split"],
                transform=ax.get_xaxis_transform())
        ax.text(xi, ymax_top, f"{medians['drift'][xi]:+.2f}",
                ha="center", va="bottom", fontsize=7, color=colors["drift"],
                transform=ax.get_xaxis_transform())
        ax.text(xi + width, ymax_top, f"{medians['xprobe'][xi]:+.2f}",
                ha="center", va="bottom", fontsize=7, color=colors["xprobe"],
                transform=ax.get_xaxis_transform())
    for xi, b in enumerate(BAND_ORDER):
        rs = CTM[CTM["band"] == b].iloc[0]
        p = float(rs["wilcoxon_split_gt_drift_p"])
        ax.text(xi, -0.10, f"p={p:.3f}{asterisks(p)}", ha="center", va="top",
                fontsize=8, color="#222", transform=ax.get_xaxis_transform())
    ax.axhline(0, color="0.3", lw=0.6, ls="--")
    ax.set_xticks(x_band)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=10)
    ax.set_ylabel(r"$\rho_{\mathrm{split}}$ / $\rho_{\mathrm{drift}}$ / $\rho_{\mathrm{xprobe}}$",
                  fontsize=10)
    all_v = []
    for _, key in box_keys:
        all_v.append(CTM_PT[key].dropna().values)
    all_v = np.concatenate(all_v)
    qmin, qmax = float(np.nanpercentile(all_v, 1)), float(np.nanpercentile(all_v, 99))
    ax.set_ylim(qmin - 0.10, qmax + 0.10)
    handles = [
        mpatches.Patch(facecolor=colors["split"], alpha=0.55, edgecolor="#444",
                       label=r"$\rho_{\mathrm{split}}$ (controlled)"),
        mpatches.Patch(facecolor=colors["drift"], alpha=0.55, edgecolor="#444",
                       label=r"$\rho_{\mathrm{null\,drift}}$"),
        mpatches.Patch(facecolor=colors["xprobe"], alpha=0.55, edgecolor="#444",
                       label=r"$\rho_{\mathrm{xprobe}}$ (cross-probe only)"),
        mlines.Line2D([], [], color="#111", lw=1.4, label="cohort median"),
    ]
    ax.text(-0.5, ymax_top + 0.05, "cohort medians (split / drift / xprobe):",
            ha="left", va="bottom", fontsize=7, color="#444",
            transform=ax.get_xaxis_transform())
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               ncol=4, frameon=False, fontsize=8.5)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIG / "ctm_headline_v2.pdf", bbox_inches="tight")
    plt.close(fig)


# ===================================================================
# Item 5 — per-(patient, band, phase) dendrogram + Ψ folder
# ===================================================================
def dendrogram_psi_per_pbphi() -> None:
    PHASES = ("rest_pre", "task_test", "rest_post")
    n_done = 0
    n_skip = 0
    for pat in PATIENTS_LIST:
        pat_dir = DEND_DIR / pat
        pat_dir.mkdir(exist_ok=True)
        for band in BAND_ORDER:
            for phi in PHASES:
                try:
                    res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                    if res is None or res.linkage_matrix is None:
                        raise FileNotFoundError("missing linkage")
                    Z = np.asarray(res.linkage_matrix)
                except Exception as e:
                    n_skip += 1
                    print(f"[dend_psi] WARN {pat} {band} {phi}: {e}")
                    continue
                heights = np.sort(Z[:, 2])
                n = len(heights)
                psi = np.zeros(n - 1)
                for i in range(n - 1):
                    if heights[i] > 0 and heights[i + 1] > 0:
                        psi[i] = n * (np.log10(heights[i + 1]) - np.log10(heights[i]))
                argmax_n = int(np.argmax(psi))
                cut_h = float(np.sqrt(heights[argmax_n] * heights[argmax_n + 1]))
                k_clusters = int(n + 1 - argmax_n - 1)

                fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.0),
                                         gridspec_kw={"width_ratios": [3, 2]})
                ax_d, ax_p = axes
                dendrogram(Z, ax=ax_d, no_labels=True,
                           color_threshold=cut_h, above_threshold_color="#888888")
                ax_d.axhline(cut_h, color="#c0392b", lw=1.1, ls="--",
                             label=f"Ψ-best cut: h={cut_h:.3f}, k={k_clusters} clusters")
                ax_d.set_yscale("log")
                tmin = max(heights[0] * 0.8, 1e-6)
                tmax = heights[-1] * 1.05
                ax_d.set_ylim(tmin, tmax)
                ax_d.set_title(f"{pat} {BRAIN_BAND_TEX_DICT[band]} ({phi})", fontsize=10)
                ax_d.set_ylabel("merge height (log)")
                ax_d.legend(fontsize=7, loc="upper left", frameon=False)

                idx = np.arange(n)
                ax_p.semilogy(idx, heights, color="#1f4e79", lw=1.4, label="merge height")
                ax_p.set_xlabel("merge index $n$")
                ax_p.set_ylabel("merge height (log)", color="#1f4e79")
                ax_p.tick_params(axis="y", labelcolor="#1f4e79")
                ax_p2 = ax_p.twinx()
                ax_p2.plot(idx[:-1], psi, color="#c0392b", lw=0.9, alpha=0.75)
                ax_p2.axvline(argmax_n, color="#c0392b", lw=0.6, ls=":")
                ax_p2.set_ylabel(r"$\Psi(n) = N \cdot [\log_{10} t_{n+1} - \log_{10} t_n]$",
                                 color="#c0392b", fontsize=8)
                ax_p2.tick_params(axis="y", labelcolor="#c0392b")
                psi_max = float(np.max(psi)) if len(psi) else 0.0
                ax_p2.text(argmax_n, psi_max * 1.05,
                           f"argmax n={argmax_n}\nΨ={psi_max:.2f}",
                           ha="center", va="bottom", fontsize=7, color="#c0392b")
                fig.tight_layout()
                fig.savefig(pat_dir / f"{band}_{phi}.pdf")
                plt.close(fig)
                n_done += 1
    print(f"[dend_psi] {n_done} figures emitted, {n_skip} skipped, output at {DEND_DIR}")


def main() -> None:
    print("[fixes] item 1 — VI smallmult trim")
    vi_smallmult_polished_fix()
    print("[fixes] item 2 — KC v2 polished")
    kc_polished_v2_fix()
    print("[fixes] item 3 — Grassmann brief polished")
    grassmann_brief_fix()
    print("[fixes] item 4 — CTM headline v2 polished")
    ctm_headline_v2_fix()
    print("[fixes] item 5 — per-(p, b, phi) dendrogram + Ψ folder")
    dendrogram_psi_per_pbphi()
    print(f"[fixes] outputs at {FIG}")


if __name__ == "__main__":
    main()
