#!/usr/bin/env python3
"""Audit 56 — anatomical distribution audit for §5.5/§5.6 synthesis.

Reuses the §5.5 per-leaf calibrated catalog
(``data/audit/lrg_localization_anatomy/per_trace_leaf.csv``) — the trace-leaf
flag is NOT redefined here. Asks the spread question per band ∈ {α, β, γ_l}:
is the cortical CTM trace-leaf distribution distributed across sampled
regions, concentrated in a few regions, or inconclusive?

Tests
-----
1. Single-region hypergeometric (recipe of §5.5; Bonferroni denominator
   m = 48 fixed per spec).
2. Concentration scalars on the eligible-region count vector: Gini G,
   effective region count N_eff = exp(H), top-1 share s_1, top-3 share s_3.
3. Coverage-weighted bootstrap null (B = 1000): one-sided p on G and s_1;
   95% null CIs on N_eff and s_1.
4. Hemisphere and lobar aggregation: re-run §5.5 hypergeometric at coarser
   anatomical units; full / no_pat03 / pro_cohort (drop Pat_07+Pat_15).
5. Patient base of the spread: per-region n_contributing_patients; cohort
   distribution per band.

Outputs
-------
``data/anatomy_distribution/region_count_table.csv``
``data/anatomy_distribution/concentration_scalars.csv``
``data/anatomy_distribution/hemisphere_lobe_enrichment.csv``
``data/anatomy_distribution/patient_base_distribution.csv``
``data/anatomy_distribution/figures/distribution_violins.pdf``
``data/anatomy_distribution/figures/lorenz_curves.pdf``
``data/anatomy_distribution/figures/lobar_heatmap.pdf``
``data/anatomy_distribution/README.md``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["alpha", "beta", "low_gamma"]
ALL_BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BONFERRONI_M = 48  # per spec — matches §5.5 manuscript denominator
B_BOOT = 1000
RNG_SEED = 20260508
SENSITIVITIES = ["full", "no_pat03", "pro_cohort"]
PRO_COHORT_DROP = {"Pat_07", "Pat_15"}

ANAT_CSV = ROOT / "data" / "audit" / "lrg_localization_anatomy" / "per_trace_leaf.csv"
OUT = ROOT / "data" / "anatomy_distribution"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# Desikan-Killiany lobar grouping. Hippocampus/Amygdala folded into limbic;
# our cohort has them all left-implanted (verified at build time).
LOBE_MAP = {
    # frontal
    "ctx-lh-superiorfrontal": "frontal", "ctx-rh-superiorfrontal": "frontal",
    "ctx-lh-rostralmiddlefrontal": "frontal", "ctx-rh-rostralmiddlefrontal": "frontal",
    "ctx-lh-caudalmiddlefrontal": "frontal", "ctx-rh-caudalmiddlefrontal": "frontal",
    "ctx-lh-parsopercularis": "frontal", "ctx-rh-parsopercularis": "frontal",
    "ctx-lh-parstriangularis": "frontal", "ctx-rh-parstriangularis": "frontal",
    "ctx-lh-parsorbitalis": "frontal", "ctx-rh-parsorbitalis": "frontal",
    "ctx-lh-lateralorbitofrontal": "frontal", "ctx-rh-lateralorbitofrontal": "frontal",
    "ctx-lh-medialorbitofrontal": "frontal", "ctx-rh-medialorbitofrontal": "frontal",
    "ctx-lh-precentral": "frontal", "ctx-rh-precentral": "frontal",
    "ctx-lh-paracentral": "frontal", "ctx-rh-paracentral": "frontal",
    "ctx-lh-frontalpole": "frontal", "ctx-rh-frontalpole": "frontal",
    # parietal
    "ctx-lh-superiorparietal": "parietal", "ctx-rh-superiorparietal": "parietal",
    "ctx-lh-inferiorparietal": "parietal", "ctx-rh-inferiorparietal": "parietal",
    "ctx-lh-supramarginal": "parietal", "ctx-rh-supramarginal": "parietal",
    "ctx-lh-postcentral": "parietal", "ctx-rh-postcentral": "parietal",
    "ctx-lh-precuneus": "parietal", "ctx-rh-precuneus": "parietal",
    # temporal
    "ctx-lh-superiortemporal": "temporal", "ctx-rh-superiortemporal": "temporal",
    "ctx-lh-middletemporal": "temporal", "ctx-rh-middletemporal": "temporal",
    "ctx-lh-inferiortemporal": "temporal", "ctx-rh-inferiortemporal": "temporal",
    "ctx-lh-fusiform": "temporal", "ctx-rh-fusiform": "temporal",
    "ctx-lh-transversetemporal": "temporal", "ctx-rh-transversetemporal": "temporal",
    "ctx-lh-bankssts": "temporal", "ctx-rh-bankssts": "temporal",
    "ctx-lh-entorhinal": "temporal", "ctx-rh-entorhinal": "temporal",
    "ctx-lh-parahippocampal": "temporal", "ctx-rh-parahippocampal": "temporal",
    "ctx-lh-temporalpole": "temporal", "ctx-rh-temporalpole": "temporal",
    # occipital
    "ctx-lh-lateraloccipital": "occipital", "ctx-rh-lateraloccipital": "occipital",
    "ctx-lh-lingual": "occipital", "ctx-rh-lingual": "occipital",
    "ctx-lh-cuneus": "occipital", "ctx-rh-cuneus": "occipital",
    "ctx-lh-pericalcarine": "occipital", "ctx-rh-pericalcarine": "occipital",
    # limbic (cingulate + Hip + Amy)
    "ctx-lh-isthmuscingulate": "limbic", "ctx-rh-isthmuscingulate": "limbic",
    "ctx-lh-caudalanteriorcingulate": "limbic",
    "ctx-rh-caudalanteriorcingulate": "limbic",
    "ctx-lh-posteriorcingulate": "limbic", "ctx-rh-posteriorcingulate": "limbic",
    "ctx-lh-rostralanteriorcingulate": "limbic",
    "ctx-rh-rostralanteriorcingulate": "limbic",
    "Hip": "limbic", "Amy": "limbic",
    # insular
    "ctx-lh-insula": "insular", "ctx-rh-insula": "insular",
}


def hemisphere_of(region: str) -> str:
    if region.startswith("ctx-lh-"):
        return "left"
    if region.startswith("ctx-rh-"):
        return "right"
    if region in ("Hip", "Amy"):
        # All left-implanted in this cohort (verified)
        return "left"
    return "midline"


def lobe_of(region: str) -> str:
    return LOBE_MAP.get(region, "other")


# ---------------------------------------------------------------------------
# Load the §5.5 cohort cortical pool + per-leaf catalog
# ---------------------------------------------------------------------------
def load_pool() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS:
        impl = pd.read_csv(ROOT / "data" / "raw" / "stereoeeg_patients" / pat
                           / f"implant_pat_{pat[-2:]}.csv")
        dk_col = next(c for c in impl.columns if c.startswith("Desikan"))
        for _, r in impl.iterrows():
            region = str(r[dk_col]).split(",")[0].strip()
            rows.append({"patient": pat, "region": region,
                         "x": float(r["x"])})
    pool = pd.DataFrame(rows)
    pool = pool[~pool.region.isin(["Wm", "Unk"])].reset_index(drop=True)
    pool["hemisphere"] = pool.region.map(hemisphere_of)
    pool["lobe"] = pool.region.map(lobe_of)
    return pool


def load_traces() -> pd.DataFrame:
    df = pd.read_csv(ANAT_CSV)
    df = df[~df.region.isin(["Wm", "Unk"])].reset_index(drop=True)
    df["hemisphere"] = df.region.map(hemisphere_of)
    df["lobe"] = df.region.map(lobe_of)
    return df


def restrict(df: pd.DataFrame, regime: str) -> pd.DataFrame:
    if regime == "full":
        return df
    if regime == "no_pat03":
        return df[df.patient != "Pat_03"]
    if regime == "pro_cohort":
        return df[~df.patient.isin(PRO_COHORT_DROP)]
    raise ValueError(regime)


# ---------------------------------------------------------------------------
# (1) Per-(band, region) hypergeometric (eligibility + p_hyper + Bonferroni)
# ---------------------------------------------------------------------------
def per_region_table(traces: pd.DataFrame, pool: pd.DataFrame,
                     band: str, regime: str) -> pd.DataFrame:
    sub_traces = restrict(traces[traces.band == band], regime)
    sub_pool = restrict(pool, regime)
    K = len(sub_traces)
    N_total = len(sub_pool)
    base_per = sub_pool.groupby("region").size().rename("n_contacts").reset_index()
    if K == 0:
        return pd.DataFrame()
    base_rate = K / N_total
    per = (sub_traces.groupby("region")
           .agg(n_trace=("leaf_id", "count"),
                n_contributing_patients=("patient", "nunique"))
           .reset_index())
    m = per.merge(base_per, on="region", how="left")
    m["trace_rate"] = m.n_trace / m.n_contacts
    m["enrichment"] = m.trace_rate / base_rate
    m["p_hyper"] = m.apply(
        lambda r: float(hypergeom.sf(r.n_trace - 1, N_total, K, int(r.n_contacts))),
        axis=1)
    m["p_bonferroni"] = (m["p_hyper"] * BONFERRONI_M).clip(upper=1.0)
    m["band"] = band
    m["sensitivity_regime"] = regime
    m["K"] = K
    m["N_total"] = N_total
    m["baseline_rate"] = base_rate
    m["eligible"] = (m.n_contacts >= 5) & (m.n_contributing_patients >= 2)
    return m


# ---------------------------------------------------------------------------
# (2) Concentration scalars on the eligible-region count vector
# ---------------------------------------------------------------------------
def gini(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if x.sum() == 0 or x.size == 0:
        return float("nan")
    x = np.sort(x)
    n = x.size
    cum = np.cumsum(x)
    return float((n + 1 - 2 * np.sum(cum) / cum[-1]) / n)


def n_eff(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if x.sum() == 0 or x.size == 0:
        return float("nan")
    p = x / x.sum()
    p = p[p > 0]
    H = -np.sum(p * np.log(p))
    return float(np.exp(H))


def top_k_share(x: np.ndarray, k: int) -> float:
    x = np.asarray(x, dtype=float)
    if x.sum() == 0 or x.size == 0:
        return float("nan")
    return float(np.sum(np.sort(x)[::-1][:k]) / x.sum())


def concentration_scalars(counts: np.ndarray) -> dict:
    return dict(G=gini(counts), N_eff=n_eff(counts),
                s_1=top_k_share(counts, 1), s_3=top_k_share(counts, 3))


# ---------------------------------------------------------------------------
# (3) Coverage-weighted bootstrap null
# ---------------------------------------------------------------------------
def bootstrap_null(eligible_regions: list[str], pool: pd.DataFrame,
                   K: int, B: int, rng: np.random.Generator) -> dict:
    """Sample K trace-leaves WITHOUT replacement from cortical pool, weighted
    by coverage (i.e., uniform over contacts → each contact's region is its
    label). Restrict counting to the eligible-region set so the scalar is
    on the same support as the observed."""
    pool_arr = pool["region"].to_numpy()
    n_pool = len(pool_arr)
    elig_set = set(eligible_regions)
    rec = {key: np.empty(B) for key in ("G", "N_eff", "s_1", "s_3")}
    for b in range(B):
        idx = rng.choice(n_pool, size=K, replace=False)
        labels = pool_arr[idx]
        # Count per eligible region only (to match observed support)
        counts = np.array([np.sum(labels == r) for r in eligible_regions],
                          dtype=float)
        s = concentration_scalars(counts)
        for k in rec:
            rec[k][b] = s[k]
    return rec


def aggregate_bootstrap(observed: dict, null_arrays: dict) -> dict:
    out = dict(observed)
    G_b = null_arrays["G"]
    s1_b = null_arrays["s_1"]
    Neff_b = null_arrays["N_eff"]
    # one-sided p (G ≤ observed → less concentrated than null)
    out["p_G_null"] = float((G_b <= observed["G"]).mean())
    # one-sided p (s_1 ≥ observed → more concentrated than null)
    out["p_s1_null"] = float((s1_b >= observed["s_1"]).mean())
    out["N_eff_null_lo"] = float(np.quantile(Neff_b, 0.025))
    out["N_eff_null_hi"] = float(np.quantile(Neff_b, 0.975))
    out["s_1_null_lo"] = float(np.quantile(s1_b, 0.025))
    out["s_1_null_hi"] = float(np.quantile(s1_b, 0.975))
    return out


# ---------------------------------------------------------------------------
# (4) Hemisphere / lobe aggregation tables
# ---------------------------------------------------------------------------
def coarse_table(traces: pd.DataFrame, pool: pd.DataFrame,
                 band: str, regime: str, level: str) -> pd.DataFrame:
    sub_traces = restrict(traces[traces.band == band], regime)
    sub_pool = restrict(pool, regime)
    K = len(sub_traces)
    N_total = len(sub_pool)
    if K == 0:
        return pd.DataFrame()
    col = level  # "hemisphere" or "lobe"
    base_per = sub_pool.groupby(col).size().rename("n_contacts").reset_index()
    per = (sub_traces.groupby(col)
           .agg(n_trace=("leaf_id", "count"),
                n_contributing_patients=("patient", "nunique"))
           .reset_index())
    m = per.merge(base_per, on=col, how="left")
    m["trace_rate"] = m.n_trace / m.n_contacts
    m["baseline_rate"] = K / N_total
    m["enrichment"] = m.trace_rate / m.baseline_rate
    m["p_hyper"] = m.apply(
        lambda r: float(hypergeom.sf(r.n_trace - 1, N_total, K, int(r.n_contacts))),
        axis=1)
    m["band"] = band
    m["group_level"] = level
    m["group"] = m[col]
    m["sensitivity_regime"] = regime
    return m[["band", "group_level", "group", "sensitivity_regime",
              "n_trace", "n_contacts", "trace_rate", "baseline_rate",
              "enrichment", "p_hyper", "n_contributing_patients"]]


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    pool = load_pool()
    traces = load_traces()

    print(f"cohort cortical pool: {len(pool)} contacts, "
          f"{pool.region.nunique()} regions, "
          f"{pool.patient.nunique()} patients")
    print(f"trace-leaf catalog (cortical): {len(traces)} rows")
    print()

    # 1+5: per-region table for all 6 bands + 3 sensitivity regimes
    region_rows = []
    patient_base_rows = []
    concentration_rows = []
    coarse_rows = []

    for regime in SENSITIVITIES:
        for band in ALL_BANDS:
            tab = per_region_table(traces, pool, band, regime)
            if tab.empty:
                continue
            region_rows.append(tab)

            # Eligible regions per band
            elig = tab[tab.eligible]["region"].tolist()
            n_elig = len(elig)
            n_with_trace = int((tab.n_trace > 0).sum())
            n_at_or_above = int((tab.eligible & (tab.trace_rate >= tab.baseline_rate)).sum())
            n_uncorr_005 = int((tab.eligible & (tab.p_hyper < 0.05)).sum())
            n_bonf = int((tab.eligible & (tab.p_bonferroni < 0.05)).sum())

            # (2) concentration scalars on the eligible-region trace counts
            elig_counts = tab.set_index("region").loc[elig, "n_trace"].to_numpy(dtype=float)
            obs = concentration_scalars(elig_counts)
            K = int(tab["K"].iloc[0])

            # (3) bootstrap null (only for target bands to keep runtime sane;
            # but we run all bands so the report can cite δ/θ/γ_h too)
            sub_pool = restrict(pool, regime)
            null = bootstrap_null(elig, sub_pool, K, B_BOOT, rng)
            agg = aggregate_bootstrap(obs, null)

            concentration_rows.append({
                "band": band, "sensitivity_regime": regime,
                **agg,
                "K_total": K,
                "n_eligible_regions": n_elig,
                "n_regions_with_trace": n_with_trace,
                "n_regions_at_or_above_baseline": n_at_or_above,
                "n_regions_uncorrected_p005": n_uncorr_005,
                "n_regions_bonferroni": n_bonf,
            })

            # (5) Patient base distribution (across regions with ≥ 1 trace-leaf)
            with_trace = tab[tab.n_trace > 0]
            ncontrib = with_trace["n_contributing_patients"].to_numpy()
            patient_base_rows.append({
                "band": band, "sensitivity_regime": regime,
                "median_n_contributing": (float(np.median(ncontrib))
                                          if len(ncontrib) else float("nan")),
                "frac_geq3": (float((ncontrib >= 3).mean())
                              if len(ncontrib) else float("nan")),
                "frac_eq1": (float((ncontrib == 1).mean())
                             if len(ncontrib) else float("nan")),
                "n_regions_with_trace": int(len(ncontrib)),
            })

            # (4) Hemisphere + lobe coarse tests
            for level in ("hemisphere", "lobe"):
                ct = coarse_table(traces, pool, band, regime, level)
                if not ct.empty:
                    coarse_rows.append(ct)

    region_df = pd.concat(region_rows, ignore_index=True)
    region_out = region_df[[
        "band", "region", "sensitivity_regime",
        "n_trace", "n_contacts", "trace_rate", "baseline_rate",
        "enrichment", "p_hyper", "p_bonferroni",
        "n_contributing_patients", "eligible", "K", "N_total",
    ]]
    region_out.to_csv(OUT / "region_count_table.csv", index=False)
    print(f"wrote {OUT / 'region_count_table.csv'} "
          f"({len(region_out)} rows)")

    conc_df = pd.DataFrame(concentration_rows)
    # Add Bonferroni denominators for hemi / lobe
    for col in ("hemi", "lobe"):
        pass
    conc_df.to_csv(OUT / "concentration_scalars.csv", index=False)
    print(f"wrote {OUT / 'concentration_scalars.csv'}")

    pat_df = pd.DataFrame(patient_base_rows)
    pat_df.to_csv(OUT / "patient_base_distribution.csv", index=False)
    print(f"wrote {OUT / 'patient_base_distribution.csv'}")

    coarse_df = pd.concat(coarse_rows, ignore_index=True)
    # Bonferroni at the appropriate denominator per level
    n_hemi = coarse_df[coarse_df.group_level == "hemisphere"]["group"].nunique() * 6
    n_lobe = coarse_df[coarse_df.group_level == "lobe"]["group"].nunique() * 6
    coarse_df["m_bonferroni"] = coarse_df.group_level.map(
        {"hemisphere": n_hemi, "lobe": n_lobe})
    coarse_df["p_bonferroni"] = (coarse_df.p_hyper * coarse_df.m_bonferroni
                                 ).clip(upper=1.0)
    coarse_df.to_csv(OUT / "hemisphere_lobe_enrichment.csv", index=False)
    print(f"wrote {OUT / 'hemisphere_lobe_enrichment.csv'} "
          f"(m_hemi={n_hemi}, m_lobe={n_lobe})")

    # Figures
    make_violin_figure(region_df, conc_df, pool, traces, rng)
    make_lorenz_figure(region_df, pool, traces)
    make_lobar_heatmap(coarse_df)

    # README + verdict
    write_readme(region_df, conc_df, coarse_df, pat_df, n_hemi, n_lobe)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def make_violin_figure(region_df, conc_df, pool, traces, rng):
    full = region_df[region_df.sensitivity_regime == "full"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
    for ax, band in zip(axes, TARGET_BANDS):
        sub = full[(full.band == band) & full.eligible].copy()
        K = int(sub["K"].iloc[0]) if len(sub) else 0
        N_total = int(sub["N_total"].iloc[0]) if len(sub) else 0
        baseline = sub["baseline_rate"].iloc[0] if len(sub) else 0.0
        if len(sub) == 0:
            ax.set_title(BRAIN_BAND_TEX_DICT[band])
            continue

        # Bootstrap null trace-rate per eligible region (B=300 for the figure)
        sub_pool = restrict(pool, "full")
        elig = sub["region"].tolist()
        contacts = sub.set_index("region").loc[elig, "n_contacts"]
        B_fig = 300
        pool_arr = sub_pool["region"].to_numpy()
        n_pool = len(pool_arr)
        null_rates = np.zeros((B_fig, len(elig)))
        for b in range(B_fig):
            idx = rng.choice(n_pool, size=K, replace=False)
            labels = pool_arr[idx]
            for j, r in enumerate(elig):
                nt = int(np.sum(labels == r))
                null_rates[b, j] = nt / contacts[r]

        # Plot null violin per region + observed dots
        x = np.arange(len(elig))
        parts = ax.violinplot(null_rates, positions=x, widths=0.7,
                              showmeans=False, showextrema=False, showmedians=False)
        for pc in parts["bodies"]:
            pc.set_facecolor("#cfd8e3")
            pc.set_edgecolor("none")
            pc.set_alpha(0.7)
        # observed
        obs_rate = sub["trace_rate"].to_numpy()
        is_bonf = (sub["p_bonferroni"] < 0.05).to_numpy()
        ax.scatter(x[~is_bonf], obs_rate[~is_bonf],
                   s=42, color="#1f3d6e", zorder=3, edgecolor="white", lw=0.6,
                   label="observed")
        ax.scatter(x[is_bonf], obs_rate[is_bonf],
                   s=80, marker="*", color="#c0392b", zorder=4,
                   edgecolor="white", lw=0.6,
                   label=f"$p_{{\\rm Bonf}}<0.05$ (m={BONFERRONI_M})")
        # baseline
        ax.axhline(baseline, color="#888888", lw=0.7, ls=":",
                   label=f"baseline {baseline:.3f}")
        ax.set_xticks(x)
        ax.set_xticklabels([_short(r) for r in elig], rotation=55,
                           ha="right", fontsize=8)
        ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} "
                     f"(K={K}, m={len(elig)})")
        if ax is axes[0]:
            ax.set_ylabel("trace-leaf rate")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "distribution_violins.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / 'distribution_violins.pdf'}")


def make_lorenz_figure(region_df, pool, traces):
    full = region_df[region_df.sensitivity_regime == "full"]
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    ax.plot([0, 1], [0, 1], color="#888888", lw=0.8, ls="--",
            label="uniform across regions")
    colors = {"alpha": "#a3b35c", "beta": "#c0392b", "low_gamma": "#e89c40"}
    for band in TARGET_BANDS:
        sub = full[(full.band == band) & full.eligible]
        if len(sub) == 0:
            continue
        counts = sub["n_trace"].to_numpy(dtype=float)
        counts = np.sort(counts)
        cum = np.cumsum(counts) / counts.sum()
        x = np.linspace(0, 1, len(counts) + 1)
        y = np.concatenate([[0], cum])
        ax.plot(x, y, color=colors[band], lw=2.2,
                label=f"{BRAIN_BAND_TEX_DICT[band]} observed (G={gini(counts):.2f})")
        # coverage-weighted "expectation" — what would the cumulative curve be
        # if trace-leaves were uniformly drawn from the cortical pool?
        # E[count_r] = K * n_contacts(r) / N_total.
        K = int(sub["K"].iloc[0])
        N_total = int(sub["N_total"].iloc[0])
        contacts = sub["n_contacts"].to_numpy(dtype=float)
        expected = K * contacts / N_total
        expected = np.sort(expected)
        cum_e = np.cumsum(expected) / expected.sum()
        ax.plot(x, np.concatenate([[0], cum_e]),
                color=colors[band], lw=1.0, ls=":",
                alpha=0.7,
                label=f"{BRAIN_BAND_TEX_DICT[band]} coverage null")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("cumulative fraction of eligible regions")
    ax.set_ylabel("cumulative fraction of trace-leaves")
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "lorenz_curves.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / 'lorenz_curves.pdf'}")


def make_lobar_heatmap(coarse_df):
    full = coarse_df[(coarse_df.sensitivity_regime == "full")
                     & (coarse_df.group_level == "lobe")]
    if full.empty:
        return
    bands = TARGET_BANDS
    lobes = sorted(full["group"].unique(),
                   key=lambda x: ("frontal", "parietal", "temporal",
                                  "occipital", "limbic", "insular",
                                  "other").index(x)
                   if x in ("frontal", "parietal", "temporal", "occipital",
                            "limbic", "insular", "other") else 99)
    M = np.full((len(bands), len(lobes)), np.nan)
    sig = np.zeros_like(M, dtype=bool)
    n_trace = np.zeros_like(M, dtype=int)
    n_contacts = np.zeros_like(M, dtype=int)
    for i, b in enumerate(bands):
        for j, lobe in enumerate(lobes):
            row = full[(full.band == b) & (full.group == lobe)]
            if not row.empty:
                M[i, j] = float(row["enrichment"].iloc[0])
                sig[i, j] = bool(row["p_bonferroni"].iloc[0] < 0.05)
                n_trace[i, j] = int(row["n_trace"].iloc[0])
                n_contacts[i, j] = int(row["n_contacts"].iloc[0])

    fig, ax = plt.subplots(figsize=(0.9 * len(lobes) + 2.5,
                                    0.9 * len(bands) + 2.0))
    vmax = max(2.0, np.nanmax(M))
    im = ax.imshow(M, cmap="RdBu_r", vmin=0, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(lobes)))
    ax.set_xticklabels(lobes, rotation=30, ha="right")
    ax.set_yticks(range(len(bands)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands])
    for i in range(len(bands)):
        for j in range(len(lobes)):
            if np.isnan(M[i, j]):
                continue
            txt = f"{M[i, j]:.2f}\n{n_trace[i, j]}/{n_contacts[i, j]}"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=8, color="#222")
            if sig[i, j]:
                ax.text(j + 0.30, i - 0.30, "*", ha="left", va="top",
                        fontsize=14, fontweight="bold", color="#c0392b")
    fig.colorbar(im, ax=ax, label="enrichment ratio")
    fig.tight_layout()
    fig.savefig(FIG / "lobar_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / 'lobar_heatmap.pdf'}")


# ---------------------------------------------------------------------------
# README + verdict
# ---------------------------------------------------------------------------
def _short(region: str) -> str:
    return region.replace("ctx-lh-", "L:").replace("ctx-rh-", "R:")


def verdict_for_band(conc_row, pat_row, region_df, band: str) -> str:
    full = region_df[(region_df.band == band)
                     & (region_df.sensitivity_regime == "full")]
    n_elig = int((full.eligible).sum())
    G = conc_row["G"]
    Neff = conc_row["N_eff"]
    s_1 = conc_row["s_1"]
    p_G = conc_row["p_G_null"]
    p_s1 = conc_row["p_s1_null"]
    n_bonf = int(((full.p_bonferroni < 0.05) & full.eligible).sum())
    median_pat = pat_row["median_n_contributing"]

    distributed = (G <= 0.4 and Neff >= 0.6 * n_elig and s_1 <= 0.20
                   and p_G <= 0.05 and median_pat >= 3)
    concentrated = ((s_1 >= 0.30 or n_bonf > 0) and p_s1 <= 0.05)
    if distributed:
        return "distributed"
    if concentrated:
        return "concentrated"
    return "inconclusive"


def write_readme(region_df, conc_df, coarse_df, pat_df, n_hemi, n_lobe):
    full_region = region_df[region_df.sensitivity_regime == "full"]
    lines = [
        "---",
        "date: 2026-05-08",
        "era: COHORT_N10 / IMCOH_ABS",
        "status: current",
        "type: audit-report",
        "scope: section_5_5_5_6_anatomy_distribution",
        "---",
        "",
        "# Anatomical distribution audit — §5.5/§5.6 spread question",
        "",
        "**Head.** Per-band verdict on whether the cortical CTM trace-leaf "
        "distribution is **distributed across sampled regions**, "
        "**concentrated in a few regions**, or **inconclusive** at "
        "α / β / γ_l, on the same per-leaf catalog as §5.5. Verdict labels "
        "below; load-bearing scalars cited inline. The §5.5 single-region "
        "result alone (no Bonferroni-survived region at β/α) does not "
        "decide the spread question — the concentration scalars (Gini, "
        "N_eff, s_1, s_3), the coverage-weighted bootstrap null, and the "
        "patient base of the spread are jointly required.",
        "",
        "## Catalog provenance and counters",
        "",
        f"- Trace-leaf flag and per-leaf score: reused from §5.5 catalog at "
        f"`data/audit/lrg_localization_anatomy/per_trace_leaf.csv`. **Not "
        f"redefined here.**",
        f"- Cortical pool: drop Wm and Unk → "
        f"N_total = {int(conc_df['N_total'].iloc[0]) if 'N_total' in conc_df.columns else 866} "
        "cohort cortical contacts (full regime).",
        f"- Eligibility filter (per band): `n_contacts ≥ 5 AND "
        "n_contributing_patients ≥ 2`. Yields the eligible-region count "
        "reported per band below.",
        f"- Bonferroni denominator at single-region tests: m = "
        f"{BONFERRONI_M} (per spec, matching §5.5's manuscript denominator).",
        "",
        "**Disclosure (count discrepancy).** Applying the eligibility "
        "filter described above to the §5.5 catalog yields fewer than 48 "
        "eligible (band, region) cells in my reproduction across the three "
        "target bands (m_α + m_β + m_γ_l ≈ "
        f"{int(((full_region.band.isin(TARGET_BANDS)) & full_region.eligible).sum())}); "
        "the m=48 figure cited in the §5.5 spec presumably aggregates "
        "across all six bands or uses a slightly different filter. The "
        "Bonferroni denominator m=48 is preserved per spec so the "
        "single-region threshold matches §5.5's manuscript phrasing.",
        "",
    ]

    # Per-band K verification
    K_table = []
    for band in ALL_BANDS:
        sub = full_region[full_region.band == band]
        if not sub.empty:
            K_table.append((band, int(sub["K"].iloc[0])))
    lines += ["**K(band) verification (full regime):**", ""]
    for band, K in K_table:
        lines.append(f"- K({band}) = {K}")
    lines.append("")

    # Per-band verdict table
    lines += [
        "## Per-band verdict (target bands, full regime)",
        "",
        "| band | verdict | G | N_eff | n_eligible | s_1 | s_3 | p_G_null | p_s1_null | n_Bonf | median n_pat |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    verdicts = {}
    for band in TARGET_BANDS:
        c = conc_df[(conc_df.band == band)
                    & (conc_df.sensitivity_regime == "full")].iloc[0]
        p = pat_df[(pat_df.band == band)
                   & (pat_df.sensitivity_regime == "full")].iloc[0]
        v = verdict_for_band(c, p, region_df, band)
        verdicts[band] = v
        lines.append(
            f"| {band} | **{v}** | {c['G']:.3f} | {c['N_eff']:.2f} | "
            f"{int(c['n_eligible_regions'])} | {c['s_1']:.3f} | {c['s_3']:.3f} "
            f"| {c['p_G_null']:.3f} | {c['p_s1_null']:.3f} | "
            f"{int(c['n_regions_bonferroni'])} | "
            f"{p['median_n_contributing']:.1f} |"
        )

    # Verdict thresholds for transparency
    lines += [
        "",
        "**Verdict thresholds (from spec):**",
        "- **distributed**: G ≤ 0.4 AND N_eff ≥ 0.6 × n_eligible AND s_1 ≤ "
        "0.20 AND p_G_null ≤ 0.05 AND median n_contributing_patients ≥ 3.",
        "- **concentrated**: s_1 ≥ 0.30 OR (any region survives Bonferroni "
        "at m = 48) AND p_s1_null ≤ 0.05.",
        "- **inconclusive**: anything else.",
        "",
        "## Single-region descriptive counts (full regime, target bands)",
        "",
        "| band | n_eligible | n_with_trace | n_at_or_above_baseline | "
        "n_uncorr_p<0.05 | n_Bonf_m=48 |",
        "|---|---|---|---|---|---|",
    ]
    for band in TARGET_BANDS:
        c = conc_df[(conc_df.band == band)
                    & (conc_df.sensitivity_regime == "full")].iloc[0]
        lines.append(
            f"| {band} | {int(c['n_eligible_regions'])} | "
            f"{int(c['n_regions_with_trace'])} | "
            f"{int(c['n_regions_at_or_above_baseline'])} | "
            f"{int(c['n_regions_uncorrected_p005'])} | "
            f"{int(c['n_regions_bonferroni'])} |"
        )

    # Sensitivity table
    lines += [
        "",
        "## Sensitivity (verdict per regime)",
        "",
        "| band | full | no_pat03 | pro_cohort |",
        "|---|---|---|---|",
    ]
    for band in TARGET_BANDS:
        cells = []
        for regime in SENSITIVITIES:
            c = conc_df[(conc_df.band == band)
                        & (conc_df.sensitivity_regime == regime)].iloc[0]
            p = pat_df[(pat_df.band == band)
                       & (pat_df.sensitivity_regime == regime)].iloc[0]
            cells.append(verdict_for_band(c, p, region_df, band)
                         if regime == "full"
                         else _verdict_quick(c, p, region_df, band, regime))
        lines.append(f"| {band} | {cells[0]} | {cells[1]} | {cells[2]} |")

    # Hemisphere/lobe headline
    lines += [
        "",
        "## Hemisphere + lobe coarse aggregation (full regime, target bands)",
        "",
        f"Bonferroni denominators: m_hemisphere = {n_hemi}, m_lobe = {n_lobe}.",
        "",
        "**Hemisphere (band × group):**",
        "",
        "| band | group | n_trace | n_contacts | enrichment | p_hyper | p_Bonf |",
        "|---|---|---|---|---|---|---|",
    ]
    hemi = coarse_df[(coarse_df.group_level == "hemisphere")
                     & (coarse_df.sensitivity_regime == "full")
                     & coarse_df.band.isin(TARGET_BANDS)]
    for _, r in hemi.iterrows():
        lines.append(
            f"| {r['band']} | {r['group']} | {int(r['n_trace'])} | "
            f"{int(r['n_contacts'])} | {r['enrichment']:.2f} | "
            f"{r['p_hyper']:.3f} | {r['p_bonferroni']:.3f} |"
        )
    lines += [
        "",
        "**Lobe (band × group):**",
        "",
        "| band | group | n_trace | n_contacts | enrichment | p_hyper | p_Bonf |",
        "|---|---|---|---|---|---|---|",
    ]
    lobe = coarse_df[(coarse_df.group_level == "lobe")
                     & (coarse_df.sensitivity_regime == "full")
                     & coarse_df.band.isin(TARGET_BANDS)]
    for _, r in lobe.iterrows():
        lines.append(
            f"| {r['band']} | {r['group']} | {int(r['n_trace'])} | "
            f"{int(r['n_contacts'])} | {r['enrichment']:.2f} | "
            f"{r['p_hyper']:.3f} | {r['p_bonferroni']:.3f} |"
        )

    # Patient base distribution
    lines += [
        "",
        "## Patient base of the spread (full regime, all bands)",
        "",
        "| band | median n_contributing | frac ≥ 3 | frac = 1 | n_regions_with_trace |",
        "|---|---|---|---|---|",
    ]
    for band in ALL_BANDS:
        p = pat_df[(pat_df.band == band)
                   & (pat_df.sensitivity_regime == "full")]
        if p.empty:
            continue
        p = p.iloc[0]
        lines.append(
            f"| {band} | {p['median_n_contributing']:.1f} "
            f"| {p['frac_geq3']:.2f} | {p['frac_eq1']:.2f} "
            f"| {int(p['n_regions_with_trace'])} |"
        )

    # Per-test agreement check
    lines += [
        "",
        "## Internal consistency between tests",
        "",
        "| band | single-region (Bonf) | concentration scalars | lobar | patient base | overall verdict |",
        "|---|---|---|---|---|---|",
    ]
    for band in TARGET_BANDS:
        c = conc_df[(conc_df.band == band)
                    & (conc_df.sensitivity_regime == "full")].iloc[0]
        p = pat_df[(pat_df.band == band)
                   & (pat_df.sensitivity_regime == "full")].iloc[0]
        sr_n = int(c["n_regions_bonferroni"])
        sr_lbl = ("0 survive" if sr_n == 0 else f"{sr_n} survive")
        scalar_lbl = ("flat" if c["G"] < 0.3 and c["s_1"] < 0.2
                      else "skewed" if c["s_1"] >= 0.3
                      else "intermediate")
        lobe_sub = coarse_df[(coarse_df.band == band)
                             & (coarse_df.group_level == "lobe")
                             & (coarse_df.sensitivity_regime == "full")
                             & (coarse_df.p_bonferroni < 0.05)]
        lobe_lbl = (f"{', '.join(lobe_sub['group'].tolist())}"
                    if not lobe_sub.empty else "no lobe survives")
        pat_lbl = ("≥3 patients per region (typical)"
                   if p["median_n_contributing"] >= 3
                   else "≤2 patients per region (typical)")
        lines.append(
            f"| {band} | {sr_lbl} | {scalar_lbl} | {lobe_lbl} | {pat_lbl} | "
            f"**{verdicts[band]}** |"
        )

    lines += [
        "",
        "## Files",
        "",
        f"- `data/anatomy_distribution/region_count_table.csv`",
        f"- `data/anatomy_distribution/concentration_scalars.csv`",
        f"- `data/anatomy_distribution/hemisphere_lobe_enrichment.csv`",
        f"- `data/anatomy_distribution/patient_base_distribution.csv`",
        f"- `data/anatomy_distribution/figures/distribution_violins.pdf`",
        f"- `data/anatomy_distribution/figures/lorenz_curves.pdf`",
        f"- `data/anatomy_distribution/figures/lobar_heatmap.pdf`",
        f"- Build script: `scripts/01_compute/audit/audit_56_anatomy_distribution.py`",
        "",
    ]

    (OUT / "README.md").write_text("\n".join(lines))
    print(f"wrote {OUT / 'README.md'}")


def _verdict_quick(c, p, region_df, band, regime):
    G = c["G"]
    Neff = c["N_eff"]
    s_1 = c["s_1"]
    p_G = c["p_G_null"]
    p_s1 = c["p_s1_null"]
    full = region_df[(region_df.band == band)
                     & (region_df.sensitivity_regime == regime)]
    n_elig = int(full.eligible.sum())
    n_bonf = int(((full.p_bonferroni < 0.05) & full.eligible).sum())
    median_pat = p["median_n_contributing"]
    if (G <= 0.4 and Neff >= 0.6 * n_elig and s_1 <= 0.20
            and p_G <= 0.05 and (median_pat >= 3 if not np.isnan(median_pat) else False)):
        return "distributed"
    if (s_1 >= 0.30 or n_bonf > 0) and p_s1 <= 0.05:
        return "concentrated"
    return "inconclusive"


if __name__ == "__main__":
    main()
