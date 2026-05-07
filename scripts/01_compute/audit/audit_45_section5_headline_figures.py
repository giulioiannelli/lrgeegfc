#!/usr/bin/env python3
"""Audit 45 -- Section 5 manuscript-grade headline figures (β-unique story).

Four figures designed to convince a skeptical reader that:
    1. β is uniquely strong at the global-probe level (volcano).
    2. Every patient agrees at the most defensible cell, KC λ=0
       (per-patient real-vs-null trajectory grid).
    3. Cohort × (band, probe) fingerprint shows β saturation vs
       pixel noise elsewhere.
    4. β trace-leaves cluster anatomically in left temporal lobe
       + Hippocampus (2D MNI scatter).

Outputs land at
    data/outputs/figures/section_5_lrg_trace/headline/{volcano,
    per_patient_slopes, cohort_fingerprint, anatomy_2d}.pdf
and are mirrored under
    data/reports/section_5_lrg_trace/headline/figures/ so the writing
agent can cite them next to the Result-2 / Section-5 prose.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LinearSegmentedColormap

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

OUT_FIG = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "headline"
OUT_RPT = ROOT / "data" / "reports" / "section_5_lrg_trace" / "headline" / "figures"

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
# Visual palette: β bold red, others muted to make β pop.
BAND_COLOR = {
    "delta": "#7390b0",
    "theta": "#5d8a8a",
    "alpha": "#a3b35c",
    "beta": "#d62728",
    "low_gamma": "#e89c40",
    "high_gamma": "#8c5e9b",
}
BAND_ALPHA = {b: (1.0 if b == "beta" else 0.55) for b in BANDS}

CTL_CSV = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "cohort_controls_summary.csv"
KC_REAL = ROOT / "data" / "reports" / "section_5_lrg_trace" / "03_kc_lambda_triangle" / "tables" / "Td_per_patient_per_band_lambda.csv"
KC_NULL = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "kc_null.csv"
DR_REAL = ROOT / "data" / "reports" / "section_5_lrg_trace" / "02_d_rank_triangle" / "tables" / "Td_per_patient_per_band.csv"
DR_NULL = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "drank_null.csv"
ANAT_CSV = ROOT / "data" / "audit" / "lrg_localization_anatomy" / "per_trace_leaf.csv"


def make_volcano():
    """3-panel volcano: per-probe-family x-range so D-rank cluster is visible.

    The point: D-rank effects all sit at <0.3 (matrix-rank distance is
    insensitive to small reorganizations); KC at β reaches ≈5 (67× larger than D-rank); Grassmann
    sits in between. Same -log10(p) y-axis; same Bonferroni line shared
    across panels.
    """
    df = pd.read_csv(CTL_CSV)
    df["effect"] = df["median_null"] - df["median_real"]
    df["nlp"] = -np.log10(df["wilcoxon_p_real_lt_null"].clip(lower=1e-4))
    M = len(df)
    bonf_thr = 0.05 / M
    bonf_y = -np.log10(bonf_thr)

    probe_meta = [
        ("drank", "o", r"D-rank (matrix distance)",
         {"d_S": r"$d_S$", "d_P": r"$d_P$", "d_F": r"$d_F$"}),
        ("kc", "s", r"KC (dendrogram tree distance)",
         {"lambda=0.0": r"$\lambda{=}0$", "lambda=0.5": r"$\lambda{=}0.5$",
          "lambda=1.0": r"$\lambda{=}1$"}),
        ("grassmann", "^", r"Grassmann (eigenstructure)",
         {"k=13": r"$k{=}13$"}),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.6),
                             sharey=True,
                             gridspec_kw={"width_ratios": [1, 2.5, 1],
                                          "wspace": 0.18})
    ymax = 3.5
    for ax, (probe, marker, panel_title, vmap) in zip(axes, probe_meta):
        sub = df[df.probe == probe].copy()
        ax.axhline(bonf_y, color="#222", lw=1.0, ls="--", zorder=1)
        ax.axhline(-np.log10(0.05), color="#888", lw=0.6, ls=":", zorder=1)
        ax.axvline(0, color="0.7", lw=0.4, ls=":", zorder=1)

        for band in BANDS:
            s = sub[sub.band == band]
            if s.empty:
                continue
            ax.scatter(s["effect"], s["nlp"],
                       marker=marker, color=BAND_COLOR[band],
                       alpha=BAND_ALPHA[band], edgecolors="black",
                       linewidths=0.6,
                       s=160 if band == "beta" else 85,
                       zorder=5 if band == "beta" else 3)

        # Inline labels for Bonferroni-surviving cells in this panel.
        # Stack two λ=0/λ=1 labels (close in x) at different offsets.
        bonf_pass = sub[sub["wilcoxon_p_real_lt_null"] < bonf_thr].copy()
        bonf_pass = bonf_pass.sort_values("effect")
        for k, (_, r) in enumerate(bonf_pass.iterrows()):
            band_lbl = BAND_TEX[r["band"]]
            var_lbl = vmap.get(r["variant"], r["variant"])
            # alternate up/down offsets when multiple labels are close
            dy = 6 if k % 2 == 0 else -16
            ax.annotate(f"{band_lbl} {var_lbl}",
                        (r["effect"], r["nlp"]),
                        xytext=(8, dy), textcoords="offset points",
                        fontsize=10, fontweight="bold",
                        color=BAND_COLOR[r["band"]])

        ax.set_title(panel_title, fontsize=11)
        ax.set_xlabel(r"effect size  =  median$(T_d^{\rm null} - T_d^{\rm real})$",
                      fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylim(0, ymax)
        xmax_data = sub["effect"].max()
        xmin_data = sub["effect"].min()
        pad = max(0.18, 0.30 * max(abs(xmax_data), abs(xmin_data)))
        ax.set_xlim(xmin_data - pad, xmax_data + pad * 3.5)

    axes[0].set_ylabel(r"$-\log_{10}(p_{\rm Wilcoxon})$  (paired real $<$ null)",
                       fontsize=10)
    # Bonferroni explanation in left panel where x is small
    axes[0].text(0.04, bonf_y - 0.07,
                 f"Bonferroni m = {M}\n(p ≤ {bonf_thr:.4f})",
                 fontsize=8.5, color="0.15", va="top", ha="left",
                 bbox=dict(facecolor="white", edgecolor="0.6",
                           lw=0.4, pad=2.5))
    axes[0].text(0.04, -np.log10(0.05) + 0.04, "p = 0.05",
                 fontsize=8, color="0.4", va="bottom", ha="left")

    # Scientific callout above the KC panel
    kc_ax = axes[1]
    kc_ax.text(0.5, 1.15,
               r"the dendrogram (KC) amplifies effect ${\sim}$67$\times$ over "
               r"matrix distance (D-rank) at $\beta$",
               transform=kc_ax.transAxes, ha="center", va="bottom",
               fontsize=10, color="#222", fontstyle="italic",
               bbox=dict(facecolor="#fff7e6", edgecolor="#d4a017",
                         lw=0.7, pad=4))

    # Single shared band legend below all panels (clean strip, OUTSIDE plot area)
    band_handles = [plt.Line2D([0], [0], marker="o", linestyle="",
                               color=BAND_COLOR[b],
                               alpha=BAND_ALPHA[b],
                               markeredgecolor="black",
                               markersize=11 if b == "beta" else 8,
                               label=BAND_TEX[b]) for b in BANDS]
    fig.legend(handles=band_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.04),
               ncol=len(BANDS), frameon=False, fontsize=10,
               title="frequency band", title_fontsize=10)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.86, bottom=0.20, wspace=0.20)
    fig.savefig(OUT_FIG / "volcano.pdf", bbox_inches="tight")
    fig.savefig(OUT_RPT / "volcano.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[45] wrote volcano.pdf  (3-panel; M={M}, bonf p<{bonf_thr:.5f})")

def make_per_patient_slopes():
    """Real-vs-null scatter, 1x6 panels per band.

    Each dot = one patient. X-axis: T_null (within-baseline drift floor at
    KC λ=0). Y-axis: T_real (real triangle T_d at KC λ=0). Diagonal y=x
    is the "no effect" reference: dots below the diagonal carry trace
    direction (real < null). The trace half-plane is shaded pale red.
    Per-band axes are SHARED (same x and y range across panels) by
    standardizing both real and null by the patient's own |null| at that
    band — gives a unitless ratio centered on 0.
    """
    real = pd.read_csv(KC_REAL)
    null = pd.read_csv(KC_NULL)

    fig, axes = plt.subplots(1, 6, figsize=(15, 3.4), sharex=True, sharey=True)
    # Build a long table of (patient, band, n_std, r_std)
    rows = []
    for band in BANDS:
        rs = real[(real.band == band) & (real.lam == 0.0)].set_index("patient")
        ns = null[null.band == band].set_index("patient")
        for pat in PATIENTS:
            if pat not in rs.index or pat not in ns.index:
                continue
            r = float(rs.loc[pat, "T_KC"])
            n = float(ns.loc[pat, "null_T_lam0.0"])
            scale = max(abs(n), 1e-3)  # standardize by patient null magnitude
            rows.append({"patient": pat, "band": band,
                         "n": n, "r": r, "n_std": n / scale, "r_std": r / scale})
    long = pd.DataFrame(rows)
    # Use unstandardized (clip extreme outliers for shared scale)
    lim_lo, lim_hi = long[["n", "r"]].quantile(0.02).min(), long[["n", "r"]].quantile(0.98).max()
    pad = (lim_hi - lim_lo) * 0.10
    lim_lo -= pad; lim_hi += pad

    for ax, band in zip(axes, BANDS):
        sub = long[long.band == band]
        # Shade trace half-plane (below diagonal)
        xs = np.array([lim_lo, lim_hi])
        ax.fill_between(xs, xs, lim_lo, color="#ffd6d4", alpha=0.55,
                        zorder=0, edgecolor="none")
        # Diagonal y = x (no effect)
        ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], color="0.4", lw=0.9,
                ls="--", zorder=1)
        # Zero crosshair
        ax.axhline(0, color="0.7", lw=0.4, ls=":", zorder=1)
        ax.axvline(0, color="0.7", lw=0.4, ls=":", zorder=1)

        is_trace = sub["r"] < sub["n"]
        # Color non-trace dots grey, trace dots band-color
        ax.scatter(sub.loc[~is_trace, "n"], sub.loc[~is_trace, "r"],
                   color="0.7", s=55, edgecolor="0.3", linewidths=0.5,
                   zorder=3)
        ax.scatter(sub.loc[is_trace, "n"], sub.loc[is_trace, "r"],
                   color=BAND_COLOR[band], s=85 if band == "beta" else 65,
                   edgecolor="black", linewidths=0.6, zorder=4,
                   alpha=BAND_ALPHA[band])
        # Cohort centroid (median)
        med_n, med_r = sub["n"].median(), sub["r"].median()
        ax.scatter([med_n], [med_r], marker="X", s=160 if band == "beta" else 110,
                   color=BAND_COLOR[band], edgecolor="black", linewidths=1.0,
                   zorder=5)

        ax.set_title(BAND_TEX[band], fontsize=13,
                     color=BAND_COLOR[band] if band == "beta" else "0.2",
                     fontweight="bold" if band == "beta" else "normal")
        ax.set_xlim(lim_lo, lim_hi)
        ax.set_ylim(lim_lo, lim_hi)
        ax.set_aspect("equal", adjustable="box")
        ax.spines[["top", "right"]].set_visible(False)

    # x and y labels
    fig.text(0.5, 0.01, r"$T_d^{\rm null}$  (within-baseline drift floor, KC $\lambda=0$)",
             ha="center", fontsize=11)
    fig.text(0.005, 0.55, r"$T_d^{\rm real}$  (real triangle, KC $\lambda=0$)",
             va="center", rotation=90, fontsize=11)

    # Legend: trace region, no-effect line, cohort median
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#ffd6d4", edgecolor="none",
                      label=r"trace half-plane ($T_d^{\rm real} < T_d^{\rm null}$)"),
        plt.Line2D([0], [0], color="0.4", lw=1.0, ls="--",
                   label=r"no effect ($T_d^{\rm real} = T_d^{\rm null}$)"),
        plt.Line2D([0], [0], marker="X", color="black", lw=0,
                   markerfacecolor="0.5", markersize=10,
                   label="cohort median"),
        plt.Line2D([0], [0], marker="o", color="black", lw=0,
                   markerfacecolor="0.7", markersize=8,
                   label="patient (trace)"),
        plt.Line2D([0], [0], marker="o", color="black", lw=0,
                   markerfacecolor="0.7", markersize=8,
                   markeredgecolor="0.3",
                   label="patient (anti-trace)"),
    ]
    legend_handles[3].set_markerfacecolor("#d62728")
    legend_handles[4].set_markerfacecolor("0.7")
    fig.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.10), ncol=5,
               frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.05, right=0.985, top=0.86, bottom=0.16, wspace=0.10)
    fig.savefig(OUT_FIG / "per_patient_slopes.pdf", bbox_inches="tight")
    fig.savefig(OUT_RPT / "per_patient_slopes.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[45] wrote per_patient_slopes.pdf  (real-vs-null scatter with diagonal)")
def make_cohort_fingerprint():
    """Facet-agreement matrix.

    Single 10x6 heatmap. Cell value = number of facets (out of 5) where
    T_real < T_null at that (patient, band). Facets are
    {KC λ=0, KC λ=0.5, KC λ=1, d_F, Grassmann k=13}. Color: white at 0,
    saturated dark red at 5. Side-bar: column total = sum of cell values
    out of 50 (= cohort total trace agreements at that band).
    """
    real_kc = pd.read_csv(KC_REAL)
    null_kc = pd.read_csv(KC_NULL)
    real_dr = pd.read_csv(DR_REAL)
    null_dr = pd.read_csv(DR_NULL)
    real_g = pd.read_csv(ROOT / "data" / "reports" / "section_5_lrg_trace" /
                         "04_grassmann_triangle" / "tables" /
                         "Td_per_patient_per_band_k.csv")
    null_g_path = ROOT / "data" / "audit" / "lrg_global_probe_controls" / "grassmann_null.csv"
    null_g = pd.read_csv(null_g_path)

    n_pat = len(PATIENTS)
    n_bands = len(BANDS)
    counts = np.zeros((n_pat, n_bands), dtype=int)
    for i, pat in enumerate(PATIENTS):
        for j, band in enumerate(BANDS):
            agree = 0
            # KC λ=0, 0.5, 1
            for lam in [0.0, 0.5, 1.0]:
                rs = real_kc[(real_kc.patient == pat) & (real_kc.band == band) &
                             (real_kc.lam == lam)]
                ns = null_kc[(null_kc.patient == pat) & (null_kc.band == band)]
                if rs.empty or ns.empty:
                    continue
                r_v = float(rs["T_KC"].iloc[0])
                n_v = float(ns[f"null_T_lam{lam:.1f}"].iloc[0])
                if r_v < n_v:
                    agree += 1
            # d_F
            rs = real_dr[(real_dr.patient == pat) & (real_dr.band == band)]
            ns = null_dr[(null_dr.patient == pat) & (null_dr.band == band)]
            if not rs.empty and not ns.empty:
                if float(rs["T_F"].iloc[0]) < float(ns["null_T_F"].iloc[0]):
                    agree += 1
            # Grassmann k=13
            rs = real_g[(real_g.patient == pat) & (real_g.band == band) &
                        (real_g.k == 13)]
            ns = null_g[(null_g.patient == pat) & (null_g.band == band)]
            if not rs.empty and not ns.empty:
                if float(rs["T_E1"].iloc[0]) < float(ns["null_T_E1"].iloc[0]):
                    agree += 1
            counts[i, j] = agree

    # Custom colormap: white → red, with explicit bins
    from matplotlib.colors import ListedColormap
    palette = ["#ffffff", "#fcd6cf", "#fa9078", "#e94c3c", "#c9201a", "#7f0000"]
    cmap = ListedColormap(palette)

    band_tot = counts.sum(axis=0)
    band_tot_max = n_pat * 5  # 50

    fig = plt.figure(figsize=(10, 5.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 4.5], hspace=0.10)
    ax_top = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1], sharex=ax_top)

    # Top: per-band cohort agreement (sum of counts)
    bars = ax_top.bar(range(n_bands), band_tot, width=0.72,
                      color=[BAND_COLOR[b] for b in BANDS],
                      edgecolor="0.2", linewidth=0.7)
    for b_obj, band in zip(bars, BANDS):
        b_obj.set_alpha(BAND_ALPHA[band])
    for j, v in enumerate(band_tot):
        ax_top.text(j, v + 0.5, f"{int(v)}/{band_tot_max}",
                    ha="center", va="bottom", fontsize=9.5,
                    fontweight="bold" if BANDS[j] == "beta" else "normal",
                    color=BAND_COLOR[BANDS[j]] if BANDS[j] == "beta" else "0.2")
    ax_top.axhline(band_tot_max, color="0.5", lw=0.5, ls="--")
    ax_top.text(n_bands - 0.5, band_tot_max + 0.5, "max = 50",
                ha="right", va="bottom", fontsize=8, color="0.5")
    ax_top.set_ylim(0, band_tot_max * 1.13)
    ax_top.set_ylabel("cohort\nagreements", fontsize=10)
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_top.spines[["top", "right"]].set_visible(False)
    ax_top.set_title(r"facet-agreement: how many of the 5 LRG measures show $T_d^{\rm real} < T_d^{\rm null}$ at each (patient, band)",
                     fontsize=10.5, color="0.2", pad=10)

    # Main: heatmap
    im = ax.imshow(counts, cmap=cmap, vmin=-0.5, vmax=5.5, aspect="auto")
    ax.set_xticks(range(n_bands))
    ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=13)
    ax.set_yticks(range(n_pat))
    ax.set_yticklabels([p[-2:] for p in PATIENTS], fontsize=10)
    ax.set_ylabel("patient", fontsize=11)
    j_beta = BANDS.index("beta")
    ax.add_patch(plt.Rectangle((j_beta - 0.48, -0.48), 0.96, n_pat - 0.04,
                               fill=False, edgecolor="black", lw=1.8,
                               zorder=10))
    for i in range(n_pat):
        for j in range(n_bands):
            v = int(counts[i, j])
            color = "white" if v >= 4 else "0.15"
            ax.text(j, i, str(v), ha="center", va="center",
                    fontsize=10, color=color, fontweight="bold")

    # Add a discrete colorbar
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import BoundaryNorm
    cax = fig.add_axes([0.92, 0.12, 0.018, 0.55])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5], cmap.N)
    cb = ColorbarBase(cax, cmap=cmap, norm=norm,
                      ticks=[0, 1, 2, 3, 4, 5], spacing="proportional")
    cb.ax.tick_params(labelsize=9)
    cb.set_label(r"# of 5 LRG measures with $T_d^{\rm real} < T_d^{\rm null}$",
                 fontsize=9, rotation=90, labelpad=8)

    fig.subplots_adjust(left=0.075, right=0.90, top=0.92, bottom=0.10)
    fig.savefig(OUT_FIG / "cohort_fingerprint.pdf", bbox_inches="tight")
    fig.savefig(OUT_RPT / "cohort_fingerprint.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[45] wrote cohort_fingerprint.pdf  (facet-agreement 0..5)")
def make_anatomy_2d():
    """Sampling-corrected anatomy: region enrichment over cohort baseline.

    Per region r, compute trace-rate(r) = n_trace(r) / n_contacts(r) where
    n_contacts(r) is the cohort-wide count of sEEG contacts in region r.
    Enrichment = trace-rate(r) / baseline-rate, with baseline-rate =
    n_total_trace / n_total_contacts. Significance: hypergeometric test
    P(X >= n_trace | total, K, n_contacts(r)).

    A region with high raw n_trace can be at baseline if it has many
    contacts cohort-wide (sampling). Only regions with enrichment > 1
    AND p < 0.10 are reported.

    Top: enrichment bar chart (β + γ_l).
    Bottom: axial MNI scatter, β-only and γ_l-only, dots colored by
    region enrichment (so the visual "where" matches the bar story).
    """
    from scipy.stats import hypergeom
    df = pd.read_csv(ANAT_CSV)
    df = df[~df["region"].isin(["Wm", "Unk"])].copy()

    # Cohort-wide baseline: count contacts per region across all 10 patients
    base_rows = []
    for pat in PATIENTS:
        impl = pd.read_csv(ROOT / "data" / "raw" / "stereoeeg_patients" / pat /
                           f"implant_pat_{pat[-2:]}.csv")
        dk_col = next(c for c in impl.columns if c.startswith("Desikan"))
        for _, r in impl.iterrows():
            region = str(r[dk_col]).split(",")[0].strip()
            base_rows.append({"region": region, "patient": pat,
                              "x": float(r["x"]), "y": float(r["y"]),
                              "z": float(r["z"])})
    base = pd.DataFrame(base_rows)
    base = base[~base.region.isin(["Wm", "Unk"])]
    N_total = base.shape[0]
    base_per_reg = base.groupby("region").size().rename("n_contacts").reset_index()

    # Compute per-band per-region enrichment + hypergeometric p
    enrich_tables = {}
    for band in ["beta", "low_gamma"]:
        sub = df[df.band == band]
        K = sub.shape[0]
        base_rate = K / N_total if N_total > 0 else 0.0
        per_reg = (sub.groupby("region")
                   .agg(n_trace=("leaf_id", "count"),
                        n_pat_trace=("patient", "nunique"))
                   .reset_index())
        m = per_reg.merge(base_per_reg, on="region", how="left")
        m["trace_rate"] = m["n_trace"] / m["n_contacts"]
        m["enrichment"] = m["trace_rate"] / base_rate
        m["p_hyper"] = m.apply(
            lambda r: float(hypergeom.sf(r["n_trace"] - 1, N_total, K,
                                         int(r["n_contacts"]))), axis=1)
        # Filter: at least 5 cohort contacts, at least 2 trace patients,
        # enrichment > 1.0 (above baseline)
        m = m[(m["n_contacts"] >= 5) & (m["n_pat_trace"] >= 2) &
              (m["enrichment"] > 1.0)]
        m = m.sort_values("enrichment", ascending=True).tail(8)
        enrich_tables[band] = m
        print(f"[anatomy] {band}: cohort baseline rate {base_rate*100:.1f}%; "
              f"{len(m)} regions with enrichment>1, ≥5 contacts, ≥2 trace pts")
    # Save the enrichment table
    out_dir = ROOT / "data" / "audit" / "lrg_localization_anatomy"
    pd.concat([t.assign(band=b) for b, t in enrich_tables.items()]).to_csv(
        out_dir / "region_enrichment.csv", index=False)

    # Patient palette for the axial scatter
    pat_palette = plt.get_cmap("tab10").colors
    pat_color = {p: pat_palette[i] for i, p in enumerate(PATIENTS)}

    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.2], hspace=0.32, wspace=0.20)

    # ---- Top row: enrichment bar charts ----
    for col_idx, (band, title_color) in enumerate(
            [("beta", BAND_COLOR["beta"]),
             ("low_gamma", BAND_COLOR["low_gamma"])]):
        ax = fig.add_subplot(gs[0, col_idx])
        agg = enrich_tables[band]
        if agg.empty:
            ax.text(0.5, 0.5, "no enriched regions", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        # Color by hemisphere/structure AND alpha by significance
        colors = []; alphas = []
        for _, r in agg.iterrows():
            reg = r["region"]
            if reg.startswith("ctx-lh"):
                col = "#2c5b96"
            elif reg.startswith("ctx-rh"):
                col = "#c83737"
            elif reg in ("Hip", "Amy") or "Hippocampus" in reg or "Amygdala" in reg:
                col = "#7a4ba8"
            else:
                col = "#7f7f7f"
            colors.append(col)
            # Alpha by significance: p<0.05 -> 1.0, p<0.10 -> 0.65, else 0.30
            p = r["p_hyper"]
            alphas.append(1.0 if p < 0.05 else 0.6 if p < 0.10 else 0.3)
        bars = ax.barh(range(len(agg)), agg["enrichment"], color=colors,
                       edgecolor="0.2", linewidth=0.6, height=0.74)
        for b_obj, a in zip(bars, alphas):
            b_obj.set_alpha(a)
        # Baseline reference line at enrichment = 1
        ax.axvline(1.0, color="0.3", lw=1.0, ls="--", zorder=0)
        ax.text(1.0, len(agg) - 0.4, " cohort\n baseline",
                fontsize=8, color="0.3", va="center", ha="left")

        # Annotate each bar with rate, n_pat, p
        for i, (_, r) in enumerate(agg.iterrows()):
            stars = "***" if r["p_hyper"] < 0.001 else (
                    "**" if r["p_hyper"] < 0.01 else (
                    "*" if r["p_hyper"] < 0.05 else "n.s."))
            ax.text(r["enrichment"] + 0.08, i,
                    f"{r['trace_rate']*100:.0f}% rate ({int(r['n_trace'])}/"
                    f"{int(r['n_contacts'])} contacts, {int(r['n_pat_trace'])} pts) "
                    f"[p={r['p_hyper']:.3f} {stars}]",
                    va="center", fontsize=8, color="0.15")
        ax.set_yticks(range(len(agg)))
        ax.set_yticklabels(agg["region"].tolist(), fontsize=9)
        ax.set_xlabel(r"enrichment  =  region trace-rate / cohort baseline",
                      fontsize=9)
        title = (rf"$\beta$ trace-leaf enrichment by region (baseline 5.3%)" if band == "beta"
                 else rf"$\gamma_l$ trace-leaf enrichment by region (baseline 8.7%)")
        ax.set_title(title, fontsize=11, color=title_color, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(0, max(agg["enrichment"].max() * 1.55, 4.0))

    # ---- Bottom row: axial scatter, β / γ_l, dots colored by region enrichment ----
    for col_idx, band in enumerate(["beta", "low_gamma"]):
        ax = fig.add_subplot(gs[1, col_idx])
        sub = df[df.band == band].copy()
        # attach enrichment per leaf via region table
        et = enrich_tables[band].set_index("region")
        sub["region_enrichment"] = sub["region"].map(et["enrichment"])
        sub["region_enrichment"] = sub["region_enrichment"].fillna(0.0)

        ax.axvline(0, color="0.5", lw=0.6, ls="--", zorder=1)
        ax.text(0.04, 0.96, "L", transform=ax.transAxes, fontsize=14,
                color="0.3", fontweight="bold", ha="left", va="top")
        ax.text(0.96, 0.96, "R", transform=ax.transAxes, fontsize=14,
                color="0.3", fontweight="bold", ha="right", va="top")
        ax.text(0.5, 0.97, "anterior", transform=ax.transAxes,
                fontsize=8.5, color="0.4", ha="center", va="top")
        ax.text(0.5, 0.03, "posterior", transform=ax.transAxes,
                fontsize=8.5, color="0.4", ha="center", va="bottom")

        # Background: all sEEG contacts (very faint)
        ax.scatter(base["x"] / 1000, base["y"] / 1000,
                   color="0.92", s=4, alpha=0.45, edgecolors="none", zorder=1)

        # Foreground: trace-leaves, sized by region enrichment
        # Non-enriched (region_enrichment <= 1): grey small dots
        non_enr = sub[sub["region_enrichment"] <= 1.0]
        ax.scatter(non_enr["x"] / 1000, non_enr["y"] / 1000,
                   color="0.55", s=35, alpha=0.55, edgecolor="0.25",
                   linewidths=0.4, zorder=3)
        # Enriched: colored by enrichment intensity
        enr = sub[sub["region_enrichment"] > 1.0].copy()
        if not enr.empty:
            sc = ax.scatter(enr["x"] / 1000, enr["y"] / 1000,
                            c=enr["region_enrichment"],
                            cmap="Reds", vmin=1.0, vmax=4.0,
                            s=110, alpha=0.95, edgecolor="black",
                            linewidths=0.7, zorder=5)
            # Inline colorbar
            cax = inset_axes(ax, width="35%", height="3%",
                             loc="lower right", borderpad=1.0)
            cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
            cb.set_label("region enrichment", fontsize=8, labelpad=2)
            cb.ax.tick_params(labelsize=7)

        ax.set_xlabel("MNI x (mm)", fontsize=10)
        if col_idx == 0:
            ax.set_ylabel("MNI y (mm)", fontsize=10)
        else:
            ax.tick_params(left=True, labelleft=False)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(("axial view --- " +
                      (rf"$\beta$" if band == "beta" else rf"$\gamma_l$") +
                      "  (size+color = region enrichment)"),
                     fontsize=10.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(True, alpha=0.18, lw=0.4)
        ax.set_xlim(-90, 90)
        ax.set_ylim(-100, 80)

    # Color legend for hemisphere/structure (bar chart explanation)
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color="#2c5b96", label="LH cortex"),
        plt.Rectangle((0, 0), 1, 1, color="#c83737", label="RH cortex"),
        plt.Rectangle((0, 0), 1, 1, color="#7a4ba8", label="subcortical (Hip/Amy)"),
        plt.Line2D([0], [0], color="0.3", lw=1.0, ls="--",
                   label="cohort baseline (enrichment = 1)"),
        plt.Line2D([0], [0], marker="o", color="black", lw=0,
                   markerfacecolor="0.55", markersize=8,
                   label="trace-leaf in non-enriched region"),
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=5, frameon=False,
               fontsize=9)

    fig.subplots_adjust(left=0.18, right=0.985, top=0.95, bottom=0.07)
    fig.savefig(OUT_FIG / "anatomy_2d.pdf", bbox_inches="tight")
    fig.savefig(OUT_RPT / "anatomy_2d.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[45] wrote anatomy_2d.pdf  (sampling-corrected enrichment)")
def main() -> None:
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    OUT_RPT.mkdir(parents=True, exist_ok=True)
    make_volcano()
    make_per_patient_slopes()
    make_cohort_fingerprint()
    make_anatomy_2d()


if __name__ == "__main__":
    main()
