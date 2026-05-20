#!/usr/bin/env python3
"""Audit 64 — implant-geometry test of cohort anti-alignment (§5.3 / §6 settle).

Tests whether per-patient implant geometry (channel count, hemispheric and
lobar coverage, spatial dispersion, epi-zone fraction) correlates across the
n=10 cohort with per-patient trace direction at the three §5.3 bands
(α, β, γ_l). Pure correlational analysis on pre-existing outputs — no LRG,
no FC, no surrogate runs.

Inputs (existing):
    data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
    data/audit/ctm_triangle/Td_per_patient_per_band.csv
    data/audit/kc_triangle/Td_per_patient_per_band_lambda.csv
    data/audit/raw_fc_phase_distance/Td_per_patient_per_band.csv
    data/raw/stereoeeg_patients/Pat_NN/implant_pat_NN.csv  (cohort)

Outputs:
    data/audit/implant_geometry/per_patient_features.csv
    data/audit/implant_geometry/per_patient_trace_measures.csv
    data/audit/implant_geometry/correlation_matrix.csv
    data/audit/implant_geometry/group_comparison.csv
    data/audit/implant_geometry/anti_alignment_count.csv
    data/audit/implant_geometry/figures/correlation_heatmap.pdf
    data/audit/implant_geometry/figures/pro_vs_anti_panel.pdf
    data/audit/implant_geometry/README.md
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, mannwhitneyu

import matplotlib.pyplot as plt
from matplotlib import rc_context
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.utils.io import load_epileptic_nodes
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
OUT_DIR = ROOT / "data" / "audit" / "implant_geometry"
FIG_DIR = OUT_DIR / "figures"

LOBE_MAP = {
    "ctx-lh-superiorfrontal": "frontal", "ctx-rh-superiorfrontal": "frontal",
    "ctx-lh-rostralmiddlefrontal": "frontal",
    "ctx-rh-rostralmiddlefrontal": "frontal",
    "ctx-lh-caudalmiddlefrontal": "frontal",
    "ctx-rh-caudalmiddlefrontal": "frontal",
    "ctx-lh-parsopercularis": "frontal", "ctx-rh-parsopercularis": "frontal",
    "ctx-lh-parstriangularis": "frontal", "ctx-rh-parstriangularis": "frontal",
    "ctx-lh-parsorbitalis": "frontal", "ctx-rh-parsorbitalis": "frontal",
    "ctx-lh-lateralorbitofrontal": "frontal",
    "ctx-rh-lateralorbitofrontal": "frontal",
    "ctx-lh-medialorbitofrontal": "frontal",
    "ctx-rh-medialorbitofrontal": "frontal",
    "ctx-lh-precentral": "frontal", "ctx-rh-precentral": "frontal",
    "ctx-lh-paracentral": "frontal", "ctx-rh-paracentral": "frontal",
    "ctx-lh-frontalpole": "frontal", "ctx-rh-frontalpole": "frontal",
    "ctx-lh-superiorparietal": "parietal",
    "ctx-rh-superiorparietal": "parietal",
    "ctx-lh-inferiorparietal": "parietal",
    "ctx-rh-inferiorparietal": "parietal",
    "ctx-lh-supramarginal": "parietal", "ctx-rh-supramarginal": "parietal",
    "ctx-lh-postcentral": "parietal", "ctx-rh-postcentral": "parietal",
    "ctx-lh-precuneus": "parietal", "ctx-rh-precuneus": "parietal",
    "ctx-lh-superiortemporal": "temporal",
    "ctx-rh-superiortemporal": "temporal",
    "ctx-lh-middletemporal": "temporal", "ctx-rh-middletemporal": "temporal",
    "ctx-lh-inferiortemporal": "temporal",
    "ctx-rh-inferiortemporal": "temporal",
    "ctx-lh-fusiform": "temporal", "ctx-rh-fusiform": "temporal",
    "ctx-lh-transversetemporal": "temporal",
    "ctx-rh-transversetemporal": "temporal",
    "ctx-lh-bankssts": "temporal", "ctx-rh-bankssts": "temporal",
    "ctx-lh-entorhinal": "temporal", "ctx-rh-entorhinal": "temporal",
    "ctx-lh-parahippocampal": "temporal", "ctx-rh-parahippocampal": "temporal",
    "ctx-lh-temporalpole": "temporal", "ctx-rh-temporalpole": "temporal",
    "ctx-lh-lateraloccipital": "occipital",
    "ctx-rh-lateraloccipital": "occipital",
    "ctx-lh-lingual": "occipital", "ctx-rh-lingual": "occipital",
    "ctx-lh-cuneus": "occipital", "ctx-rh-cuneus": "occipital",
    "ctx-lh-pericalcarine": "occipital", "ctx-rh-pericalcarine": "occipital",
    "ctx-lh-isthmuscingulate": "limbic",
    "ctx-rh-isthmuscingulate": "limbic",
    "ctx-lh-caudalanteriorcingulate": "limbic",
    "ctx-rh-caudalanteriorcingulate": "limbic",
    "ctx-lh-posteriorcingulate": "limbic",
    "ctx-rh-posteriorcingulate": "limbic",
    "ctx-lh-rostralanteriorcingulate": "limbic",
    "ctx-rh-rostralanteriorcingulate": "limbic",
    "Hip": "limbic", "Amy": "limbic",
    "ctx-lh-insula": "insular", "ctx-rh-insula": "insular",
}
LOBES = ["frontal", "temporal", "parietal", "occipital",
         "limbic", "insular", "other"]


def hemi_of(region: str) -> str:
    if region.startswith("ctx-lh-"):
        return "left"
    if region.startswith("ctx-rh-"):
        return "right"
    if region in ("Hip", "Amy"):
        # all left-implanted in this cohort, per audit_56
        return "left"
    return "midline"


def lobe_of(region: str) -> str:
    return LOBE_MAP.get(region, "other")


# ---------------------------------------------------------------------------
# Step 1 — per-patient implant-geometry features
# ---------------------------------------------------------------------------
def build_features() -> pd.DataFrame:
    rows = []
    for pat in COHORT:
        impl_path = (ROOT / "data" / "raw" / "stereoeeg_patients" / pat
                     / f"implant_pat_{pat[-2:]}.csv")
        impl = pd.read_csv(impl_path)
        dk_col = next(c for c in impl.columns if c.startswith("Desikan"))
        regions = [str(r).split(",")[0].strip() for r in impl[dk_col]]
        hemis = [hemi_of(r) for r in regions]
        lobes = [lobe_of(r) for r in regions]
        # coordinates: implant CSVs are in micrometers (×1000 mm) per
        # observed magnitudes; normalize to mm for centroid + dispersion.
        x = np.asarray(impl["x"], dtype=float) / 1000.0
        y = np.asarray(impl["y"], dtype=float) / 1000.0
        z = np.asarray(impl["z"], dtype=float) / 1000.0
        coords = np.column_stack([x, y, z])

        N = len(regions)
        N_L = sum(h == "left" for h in hemis)
        N_R = sum(h == "right" for h in hemis)
        N_M = sum(h == "midline" for h in hemis)
        lobe_counts = {f"N_{lobe}": sum(lo == lobe for lo in lobes)
                       for lobe in LOBES}
        B_hemi = (N_L - N_R) / N if N else float("nan")
        # lobe_max ignores 'other' so the headline is interpretable
        # anatomically; if all 'other', falls back to 'other'.
        non_other = {k: v for k, v in lobe_counts.items()
                     if k != "N_other"}
        if any(non_other.values()):
            lobe_max_key = max(non_other, key=non_other.get)
        else:
            lobe_max_key = "N_other"
        lobe_max = lobe_max_key.replace("N_", "")
        frac_lobe_max = lobe_counts[lobe_max_key] / N if N else float("nan")

        centroid = coords.mean(axis=0)
        dists = np.linalg.norm(coords - centroid, axis=1)
        sigma_disp = float(dists.std(ddof=1))

        epi_set = set(load_epileptic_nodes(pat))
        labels_norm = [str(lbl).strip() for lbl in impl["label"]]
        N_epi = sum(lbl in epi_set for lbl in labels_norm)
        frac_epi = N_epi / N if N else float("nan")

        row = {
            "patient": pat,
            "N": N,
            "N_L": N_L,
            "N_R": N_R,
            "N_M": N_M,
            **lobe_counts,
            "B_hemi": B_hemi,
            "lobe_max": lobe_max,
            "frac_lobe_max": frac_lobe_max,
            "centroid_x": centroid[0],
            "centroid_y": centroid[1],
            "centroid_z": centroid[2],
            "sigma_disp": sigma_disp,
            "N_epi": N_epi,
            "frac_epi": frac_epi,
        }
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 2 — per-patient trace-direction scalars
# ---------------------------------------------------------------------------
def build_trace_measures() -> pd.DataFrame:
    # obs_rho, obs_z from matched-strength split-baseline
    mss = pd.read_csv(OUT_DIR.parent / "matched_strength_surrogate_split_baseline"
                      / "per_patient_per_band.csv")
    mss = mss[mss.band.isin(BANDS)][["patient", "band", "obs_rho", "obs_z"]]

    # KC λ=0 and λ=1: T_KC column. Sign convention: positive T_KC means
    # d_tt_post > d_pre_tt (= reset / drift away); negative T_KC means
    # d_tt_post < d_pre_tt (= trace direction). To match other scalars
    # (positive = trace direction), flip sign.
    kc = pd.read_csv(OUT_DIR.parent / "kc_triangle"
                     / "Td_per_patient_per_band_lambda.csv")
    kc = kc[kc.band.isin(BANDS) & kc.lam.isin([0.0, 1.0])].copy()
    kc["trace_score"] = -kc["T_KC"]
    kc_l0 = (kc[kc.lam == 0.0][["patient", "band", "trace_score"]]
             .rename(columns={"trace_score": "T_KC_l0"}))
    kc_l1 = (kc[kc.lam == 1.0][["patient", "band", "trace_score"]]
             .rename(columns={"trace_score": "T_KC_l1"}))

    # Substrate T_d^(d_S) from raw_fc_phase_distance. Column S in that
    # CSV is the raw triangle inequality residual T_d^(d_S); negative =
    # trace direction. Flip sign for consistency.
    rfc = pd.read_csv(OUT_DIR.parent / "raw_fc_phase_distance"
                      / "Td_per_patient_per_band.csv")
    rfc = rfc[rfc.band.isin(BANDS)][["patient", "band", "S"]].copy()
    rfc["T_d_dS"] = -rfc["S"]
    rfc = rfc[["patient", "band", "T_d_dS"]]

    out = (mss
           .merge(kc_l0, on=["patient", "band"])
           .merge(kc_l1, on=["patient", "band"])
           .merge(rfc, on=["patient", "band"]))
    return out


# ---------------------------------------------------------------------------
# Step 3 — correlation matrix
# ---------------------------------------------------------------------------
TRACE_MEASURES = ["obs_rho", "obs_z", "T_KC_l0", "T_KC_l1", "T_d_dS"]
FEATURES = ["N", "B_hemi", "frac_lobe_max", "sigma_disp", "frac_epi",
            "N_temporal", "N_frontal", "N_L", "N_R"]


def compute_correlations(feat: pd.DataFrame,
                         trace: pd.DataFrame) -> pd.DataFrame:
    merged = trace.merge(feat, on="patient")
    rows = []
    for band in BANDS:
        sub = merged[merged.band == band]
        for measure in TRACE_MEASURES:
            for feature in FEATURES:
                m_vals = sub[measure].to_numpy()
                f_vals = sub[feature].to_numpy()
                if np.isnan(m_vals).any() or np.isnan(f_vals).any():
                    continue
                rho_s, p_s = spearmanr(f_vals, m_vals)
                r_p, p_p = pearsonr(f_vals, m_vals)
                rows.append({
                    "band": band,
                    "trace_measure": measure,
                    "feature": feature,
                    "n_patients": len(sub),
                    "spearman_rho": rho_s,
                    "spearman_p": p_s,
                    "pearson_r": r_p,
                    "pearson_p": p_p,
                })
    df = pd.DataFrame(rows)
    df["bh_q"] = bh_fdr(df.spearman_p.tolist())
    df["bh_significant_q05"] = df.bh_q < 0.05
    return df


# ---------------------------------------------------------------------------
# Step 4 — pro-vs-anti group comparison + anti-alignment count
# ---------------------------------------------------------------------------
def group_comparison(feat: pd.DataFrame,
                     trace: pd.DataFrame) -> pd.DataFrame:
    merged = trace.merge(feat, on="patient")
    rows = []
    for band in BANDS:
        sub = merged[merged.band == band]
        pro = sub[sub.obs_rho > 0]
        anti = sub[sub.obs_rho <= 0]
        for feature in FEATURES:
            pro_v = pro[feature].to_numpy()
            anti_v = anti[feature].to_numpy()
            if len(pro_v) == 0 or len(anti_v) == 0:
                u_stat, p_val = float("nan"), float("nan")
            else:
                try:
                    u_stat, p_val = mannwhitneyu(
                        pro_v, anti_v, alternative="two-sided")
                except ValueError:
                    u_stat, p_val = float("nan"), float("nan")
            rows.append({
                "band": band,
                "feature": feature,
                "n_pro": int(len(pro_v)),
                "n_anti": int(len(anti_v)),
                "pro_median": float(np.median(pro_v)) if len(pro_v) else np.nan,
                "anti_median": float(np.median(anti_v)) if len(anti_v) else np.nan,
                "u_stat": float(u_stat) if u_stat == u_stat else np.nan,
                "p_two_sided": float(p_val) if p_val == p_val else np.nan,
            })
    return pd.DataFrame(rows)


def anti_alignment_count(trace: pd.DataFrame) -> pd.DataFrame:
    # 4 probes × 3 bands = 12 cells per patient. Probe is "anti" when
    # trace_score < 0 (= opposite of trace direction).
    rows = []
    for pat in COHORT:
        sub = trace[trace.patient == pat]
        cells = []
        for _, r in sub.iterrows():
            for probe in ("obs_rho", "T_KC_l0", "T_KC_l1", "T_d_dS"):
                cells.append(1 if r[probe] < 0 else 0)
        n_anti = sum(cells)
        n_total = len(cells)
        if n_anti >= 6:
            cls = "consistently_anti"
        elif n_anti <= 2:
            cls = "consistently_pro"
        else:
            cls = "mixed"
        rows.append({
            "patient": pat,
            "n_anti": n_anti,
            "n_total": n_total,
            "frac_anti": n_anti / n_total if n_total else np.nan,
            "class": cls,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 5 — visualizations
# ---------------------------------------------------------------------------
def plot_correlation_heatmap(corr: pd.DataFrame) -> None:
    n_measures = len(TRACE_MEASURES)
    fig, axes = plt.subplots(1, n_measures, figsize=(3.2 * n_measures, 4.0),
                             sharey=True)
    if n_measures == 1:
        axes = [axes]
    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    for ax, measure in zip(axes, TRACE_MEASURES):
        sub = corr[corr.trace_measure == measure]
        mat = (sub.pivot(index="feature", columns="band",
                         values="spearman_rho")
               .reindex(index=FEATURES, columns=BANDS))
        im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", norm=norm,
                       aspect="auto")
        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels(BANDS, rotation=30, ha="right")
        ax.set_title(measure, fontsize=9)
        if ax is axes[0]:
            ax.set_yticks(range(len(FEATURES)))
            ax.set_yticklabels(FEATURES)
        for i, feature in enumerate(FEATURES):
            for j, band in enumerate(BANDS):
                row = sub[(sub.feature == feature) & (sub.band == band)]
                if row.empty:
                    continue
                p = float(row.spearman_p.iloc[0])
                q = float(row.bh_q.iloc[0])
                marker = ""
                if q < 0.05:
                    marker = "★★"
                elif p < 0.05:
                    marker = "★"
                if marker:
                    ax.text(j, i, marker, ha="center", va="center",
                            color="black", fontsize=8)
    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label(r"Spearman $\rho$", fontsize=9)
    fig.savefig(FIG_DIR / "correlation_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_pro_vs_anti(feat: pd.DataFrame, trace: pd.DataFrame,
                     corr: pd.DataFrame,
                     anti_count: pd.DataFrame) -> None:
    beta_corr = corr[(corr.band == "beta")
                     & (corr.trace_measure == "obs_rho")]
    top3 = (beta_corr.assign(absrho=beta_corr.spearman_rho.abs())
            .sort_values("absrho", ascending=False).head(3))
    top_features = top3.feature.tolist()

    merged = trace.merge(feat, on="patient").merge(
        anti_count[["patient", "n_anti"]], on="patient")

    n_rows = len(top_features)
    n_cols = len(BANDS)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.6 * n_cols, 2.6 * n_rows),
                             sharex=False)
    if n_rows == 1:
        axes = np.array([axes])
    cmap = plt.get_cmap("coolwarm")
    norm = plt.Normalize(vmin=0, vmax=12)
    for i, feature in enumerate(top_features):
        for j, band in enumerate(BANDS):
            ax = axes[i, j]
            sub = merged[merged.band == band]
            x_jitter = np.random.default_rng(seed=42).normal(0, 0.02,
                                                             size=len(sub))
            for _, r in sub.iterrows():
                color = cmap(norm(r["n_anti"]))
                ax.scatter(r["obs_rho"], r[feature],
                           c=[color], s=60, edgecolor="k",
                           linewidth=0.5)
                ax.annotate(r["patient"].replace("Pat_", ""),
                            xy=(r["obs_rho"], r[feature]),
                            xytext=(3, 3), textcoords="offset points",
                            fontsize=6)
            ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
            ax.set_xlabel(r"$\rho_{\mathrm{split}}$" if i == n_rows - 1 else "",
                          fontsize=8)
            ax.set_ylabel(feature if j == 0 else "", fontsize=8)
            ax.set_title(f"{band}" if i == 0 else "", fontsize=9)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
    cbar.set_label("anti-alignment count (0–12)", fontsize=8)
    fig.savefig(FIG_DIR / "pro_vs_anti_panel.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Step 6+7 — verdict + README disclosure
# ---------------------------------------------------------------------------
def verdict_for_band(band: str, corr: pd.DataFrame) -> tuple[str, dict]:
    """Pick the verdict label for a band according to ticket criteria.

    strong: any feature × {obs_rho, obs_z} with |Spearman ρ| ≥ 0.7 AND
            bh_q < 0.05.
    moderate: any feature × {obs_rho, obs_z} with |Spearman ρ| ∈ [0.5, 0.7)
              AND uncorrected p < 0.05.
    null: otherwise.
    """
    sub = corr[(corr.band == band)
               & (corr.trace_measure.isin(["obs_rho", "obs_z"]))]
    sub = sub.assign(absrho=sub.spearman_rho.abs())
    strong = sub[(sub.absrho >= 0.7) & (sub.bh_q < 0.05)]
    if not strong.empty:
        top = strong.sort_values("absrho", ascending=False).iloc[0]
        return "strong", top.to_dict()
    moderate = sub[(sub.absrho >= 0.5) & (sub.absrho < 0.7)
                   & (sub.spearman_p < 0.05)]
    if not moderate.empty:
        top = moderate.sort_values("absrho", ascending=False).iloc[0]
        return "moderate", top.to_dict()
    top = sub.sort_values("absrho", ascending=False).iloc[0]
    return "null", top.to_dict()


def disclosure_paragraph(verdict: str, band: str, top: dict,
                         feat: pd.DataFrame, trace: pd.DataFrame) -> str:
    feature = top["feature"]
    rho = top["spearman_rho"]
    p = top["spearman_p"]
    q = top["bh_q"]
    merged = trace.merge(feat, on="patient")
    sub = merged[merged.band == band]
    pro = sub[sub.obs_rho > 0]
    anti = sub[sub.obs_rho <= 0]
    pro_med = float(pro[feature].median()) if not pro.empty else float("nan")
    anti_med = float(anti[feature].median()) if not anti.empty else float("nan")
    anti_pats = ", ".join(sorted(anti.patient.tolist()))

    if verdict == "strong":
        return (
            f"Cohort heterogeneity in the per-patient {band} trace direction "
            f"correlates with implant geometry: {feature} correlates with "
            f"the per-patient {band} obs_rho at Spearman ρ = {rho:+.3f} "
            f"(p = {p:.3g}, BH-FDR q = {q:.3g}). The "
            f"{len(anti)} anti-aligned patients ({anti_pats}) sit at the "
            f"{'negative' if rho > 0 else 'positive'} tail of {feature} "
            f"(median {anti_med:.3f} against pro-aligned median "
            f"{pro_med:.3f}). This is a result; the framing of the "
            f"anti-aligned sub-cohort as implant-geometry driven is "
            f"supported by the data, not by post-hoc speculation.")
    if verdict == "moderate":
        return (
            f"Cohort heterogeneity in the per-patient {band} trace direction "
            f"shows directional evidence of association with implant "
            f"geometry: {feature} correlates with per-patient {band} obs_rho "
            f"at Spearman ρ = {rho:+.3f} (p = {p:.3g} uncorrected, BH-FDR "
            f"q = {q:.3g}, not surviving q = 0.05 multiple-comparison "
            f"correction). The directional association is consistent with "
            f"the hypothesis that implant-geometry sampling contributes to "
            f"the per-patient distribution but is not formally established "
            f"at q = 0.05 in this 10-patient cohort.")
    return (
        f"Cohort heterogeneity in the per-patient {band} trace direction "
        f"does not correlate with any tested implant-geometry feature at "
        f"the 0.05 uncorrected level for the load-bearing trace scalars "
        f"(obs_rho, obs_z). The {len(anti)} anti-aligned patients "
        f"({anti_pats}) are not distinguishable from the {len(pro)} "
        f"pro-aligned patients on hemisphere balance, lobar coverage, "
        f"spatial dispersion, or epileptogenic-zone fraction. The cohort "
        f"heterogeneity is documented; the mechanism is not identified in "
        f"this 10-patient sample. The largest absolute Spearman ρ across "
        f"the load-bearing measures is "
        f"|ρ| = {abs(top['spearman_rho']):.3f} for {feature}.")


def write_readme(feat: pd.DataFrame, trace: pd.DataFrame,
                 corr: pd.DataFrame, group: pd.DataFrame,
                 anti_count: pd.DataFrame) -> None:
    band_verdicts = {b: verdict_for_band(b, corr) for b in BANDS}

    n_passing = sum(v in ("strong", "moderate")
                    for v, _ in band_verdicts.values())

    if n_passing == 3:
        global_verdict = (
            "All three trace bands (α, β, γ_l) admit at least moderate "
            "implant-geometry correlation. The §6.3 outlook can carry an "
            "implant-driven cohort-heterogeneity framing without bracket.")
    elif n_passing >= 1:
        bands_pass = ", ".join(
            b for b, (v, _) in band_verdicts.items()
            if v in ("strong", "moderate"))
        global_verdict = (
            f"{n_passing}/3 trace bands admit at least moderate "
            f"implant-geometry correlation ({bands_pass}). The §6.3 "
            "outlook framing is band-resolved: implant-driven for the "
            "passing bands, undetermined for the others.")
    else:
        global_verdict = (
            "No trace band admits implant-geometry correlation at the "
            "0.05 uncorrected level for the load-bearing trace scalars. "
            "The §6.3 outlook drops the implant-geometry framing and "
            "documents cohort heterogeneity without proposed mechanism.")

    lines = [
        "---",
        "name: implant_geometry_test",
        "scope: settle_section_5_section_6_implant_geometry_attribution",
        "date: 2026-05-11",
        "status: complete",
        "---",
        "",
        "# Implant-geometry test of cohort anti-alignment",
        "",
        "Head. Cohort heterogeneity in the per-patient trace direction "
        "at α, β, γ_l is tested against per-patient implant-geometry "
        "features (channel count, hemispheric balance, lobar coverage, "
        "spatial dispersion, epileptogenic-zone fraction). Pure "
        "correlational analysis on existing per-patient outputs: no LRG, "
        "no FC, no surrogate runs. Verdict per band picks one of "
        "`strong` / `moderate` / `null` against the load-bearing trace "
        "scalars (obs_rho, obs_z from the matched-strength split-baseline "
        "surrogate).",
        "",
        "## Per-band verdicts",
        "",
    ]
    for band in BANDS:
        v, top = band_verdicts[band]
        lines.append(
            f"- **{band}**: `{v}` "
            f"— top feature {top['feature']}, "
            f"Spearman ρ = {top['spearman_rho']:+.3f} "
            f"(p = {top['spearman_p']:.3g}, BH-FDR q = {top['bh_q']:.3g}) "
            f"vs {top['trace_measure']}."
        )
    lines.extend([
        "",
        "## Global verdict",
        "",
        global_verdict,
        "",
        "## Anti-alignment count per patient (across 12 = 4 probes × 3 bands)",
        "",
        "| patient | n_anti / 12 | class |",
        "|---|---|---|",
    ])
    for _, r in anti_count.iterrows():
        lines.append(
            f"| {r['patient']} | {r['n_anti']} / {r['n_total']} | {r['class']} |"
        )

    lines.extend([
        "",
        "## §5 / §6 / §6.3 disclosure paragraphs (band-resolved)",
        "",
    ])
    for band in BANDS:
        v, top = band_verdicts[band]
        para = disclosure_paragraph(v, band, top, feat, trace)
        lines.append(f"### {band} — `{v}` verdict")
        lines.append("")
        lines.append(para)
        lines.append("")

    lines.extend([
        "## Top three implant features at β (|Spearman ρ| against obs_rho)",
        "",
    ])
    beta_top = (corr[(corr.band == "beta")
                     & (corr.trace_measure == "obs_rho")]
                .assign(absrho=lambda d: d.spearman_rho.abs())
                .sort_values("absrho", ascending=False)
                .head(3))
    lines.append("| feature | Spearman ρ | p | BH-FDR q |")
    lines.append("|---|---|---|---|")
    for _, r in beta_top.iterrows():
        lines.append(
            f"| {r['feature']} | {r['spearman_rho']:+.3f} | "
            f"{r['spearman_p']:.3g} | {r['bh_q']:.3g} |"
        )

    lines.extend([
        "",
        "## Provenance",
        "",
        "- Cohort: " + ", ".join(COHORT),
        "- Trace bands: " + ", ".join(BANDS),
        "- Trace measures: " + ", ".join(TRACE_MEASURES),
        "- Implant features: " + ", ".join(FEATURES),
        "- Sources:",
        "  - `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` (obs_rho, obs_z)",
        "  - `data/audit/kc_triangle/Td_per_patient_per_band_lambda.csv` (T_KC λ=0, λ=1; sign-flipped)",
        "  - `data/audit/raw_fc_phase_distance/Td_per_patient_per_band.csv` (T_d^(d_S); sign-flipped)",
        "  - `data/raw/stereoeeg_patients/Pat_NN/implant_pat_NN.csv` (implant geometry)",
        "- BH-FDR over m = "
        + str(len(BANDS) * len(TRACE_MEASURES) * len(FEATURES))
        + " correlation cells at q = 0.05.",
        "- Build script: `scripts/01_compute/audit/audit_64_implant_geometry_test.py`",
        "",
    ])
    (OUT_DIR / "README.md").write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    feat = build_features()
    feat.to_csv(OUT_DIR / "per_patient_features.csv", index=False)
    print("[Step 1] per-patient features:")
    print(feat[["patient", "N", "N_L", "N_R", "B_hemi", "lobe_max",
                "frac_lobe_max", "sigma_disp", "frac_epi"]]
          .to_string(index=False))

    trace = build_trace_measures()
    trace.to_csv(OUT_DIR / "per_patient_trace_measures.csv", index=False)
    print("\n[Step 2] per-patient trace measures (head):")
    print(trace.head().to_string(index=False))

    corr = compute_correlations(feat, trace)
    corr.to_csv(OUT_DIR / "correlation_matrix.csv", index=False)
    print(f"\n[Step 3] correlation matrix rows: {len(corr)}")
    print("Top |Spearman ρ| cells:")
    print(corr.assign(absrho=corr.spearman_rho.abs())
          .sort_values("absrho", ascending=False).head(8)
          .to_string(index=False))

    group = group_comparison(feat, trace)
    group.to_csv(OUT_DIR / "group_comparison.csv", index=False)
    anti = anti_alignment_count(trace)
    anti.to_csv(OUT_DIR / "anti_alignment_count.csv", index=False)
    print("\n[Step 4] anti-alignment counts:")
    print(anti.to_string(index=False))

    print("\n[Step 5] rendering figures …")
    with rc_context({}):
        plot_correlation_heatmap(corr)
        plot_pro_vs_anti(feat, trace, corr, anti)
    print(f"  → {FIG_DIR / 'correlation_heatmap.pdf'}")
    print(f"  → {FIG_DIR / 'pro_vs_anti_panel.pdf'}")

    write_readme(feat, trace, corr, group, anti)
    print(f"\n[Step 6/7] {OUT_DIR / 'README.md'} written.")


if __name__ == "__main__":
    main()
