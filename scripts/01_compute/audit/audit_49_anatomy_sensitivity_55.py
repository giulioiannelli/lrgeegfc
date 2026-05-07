#!/usr/bin/env python3
"""Audit 49 -- Section 5.5 anatomy sensitivity tests.

Three sensitivity tests + headline verification on top of measure 13's
calibrated trace-leaf anatomy:

A. Pro-cohort β enrichment (n=8, drop Pat_07 + Pat_15 anti-aligned)
B. Pat_02 dropout for low-γ ctx-lh-fusiform
C. Anti-aligned per-patient regional characterization (Pat_07, Pat_15)
+ verify 'with Pat_03' / 'without Pat_03' headline cells for β + low-γ

Methodology mirrors audit_45.make_anatomy_2d:
  - Cohort baseline rate K/N from per_trace_leaf rows + per-patient
    implant tables (Wm, Unk excluded).
  - Per region r: trace_rate = n_trace(r) / n_contacts(r),
    enrichment = trace_rate / baseline_rate.
  - Hypergeometric test P(X >= n_trace | N_total, K, n_contacts(r)).
  - Filter: n_contacts ≥ 5, n_pat_trace ≥ 2, enrichment > 1
    (relaxed to n_pat_trace ≥ 1 in single-patient Test C).

Outputs all enrichment tables under
    data/audit/lrg_localization_anatomy/sensitivity_55/
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
ANTI_ALIGNED = ["Pat_07", "Pat_15"]
PRO_COHORT = [p for p in PATIENTS if p not in ANTI_ALIGNED]

PER_LEAF_CSV = (ROOT / "data" / "audit" / "lrg_localization_anatomy"
                / "per_trace_leaf.csv")
RAW_ROOT = ROOT / "data" / "raw" / "stereoeeg_patients"
OUT_DIR = (ROOT / "data" / "audit" / "lrg_localization_anatomy"
           / "sensitivity_55")

# Bonferroni multiplicity for §5.5: 48 = 6 bands × 8 candidate regions.
BONF_M = 48
BONF_THR = 0.05 / BONF_M  # ≈ 0.00104


def load_cohort_contacts(patients: list[str]) -> pd.DataFrame:
    """Long-format contact table across the given patients (Wm/Unk dropped)."""
    rows = []
    for pat in patients:
        impl = pd.read_csv(RAW_ROOT / pat / f"implant_pat_{pat[-2:]}.csv")
        dk_col = next(c for c in impl.columns if c.startswith("Desikan"))
        for _, r in impl.iterrows():
            region = str(r[dk_col]).split(",")[0].strip()
            rows.append({"patient": pat, "region": region})
    df = pd.DataFrame(rows)
    return df[~df.region.isin(["Wm", "Unk"])].reset_index(drop=True)


def enrichment_table(
    leaf_df: pd.DataFrame,
    base_df: pd.DataFrame,
    band: str,
    *,
    min_contacts: int = 5,
    min_pat_trace: int = 2,
    min_enrichment: float = 1.0,
) -> tuple[pd.DataFrame, dict]:
    """Compute per-region enrichment for one band on the supplied subset.

    leaf_df: per_trace_leaf restricted to the patients we are testing.
    base_df: cohort contact table for the same patient subset.
    Returns (per_region_table, dict_of_summary_stats).
    """
    sub = leaf_df[leaf_df.band == band].copy()
    sub = sub[~sub.region.isin(["Wm", "Unk"])]
    K = len(sub)
    N_total = len(base_df)
    base_rate = K / N_total if N_total > 0 else 0.0

    base_per_reg = base_df.groupby("region").size().rename("n_contacts").reset_index()
    per_reg = (sub.groupby("region")
               .agg(n_trace=("leaf_id", "count"),
                    n_pat_trace=("patient", "nunique"))
               .reset_index())
    m = per_reg.merge(base_per_reg, on="region", how="left").fillna({"n_contacts": 0})
    m["n_contacts"] = m["n_contacts"].astype(int)
    m["trace_rate"] = m["n_trace"] / m["n_contacts"].replace(0, np.nan)
    m["enrichment"] = m["trace_rate"] / base_rate if base_rate > 0 else np.nan
    m["p_hyper"] = m.apply(
        lambda r: float(hypergeom.sf(r["n_trace"] - 1, N_total, K,
                                     int(r["n_contacts"]))) if r["n_contacts"] > 0 else np.nan,
        axis=1)
    m["bonf_pass"] = m["p_hyper"] < BONF_THR
    m_filt = m[(m["n_contacts"] >= min_contacts)
               & (m["n_pat_trace"] >= min_pat_trace)
               & (m["enrichment"].fillna(0) > min_enrichment)].copy()
    m_filt = m_filt.sort_values("enrichment", ascending=False)
    summary = {"band": band, "K": int(K), "N_total": int(N_total),
               "base_rate": float(base_rate)}
    return m_filt, summary


def run_subset(
    leaf_df_full: pd.DataFrame,
    label: str,
    patients: list[str],
    *,
    bands=("beta", "low_gamma"),
    min_pat_trace: int = 2,
    out_csv: Path | None = None,
) -> tuple[pd.DataFrame, list[dict]]:
    """Compute enrichment for the given patient subset across bands."""
    leaf = leaf_df_full[leaf_df_full.patient.isin(patients)].copy()
    base = load_cohort_contacts(patients)
    summaries = []
    rows = []
    for b in bands:
        tbl, summ = enrichment_table(leaf, base, b, min_pat_trace=min_pat_trace)
        summ["label"] = label
        summaries.append(summ)
        for _, r in tbl.iterrows():
            rows.append({"label": label, **r.to_dict()})
    out = pd.DataFrame(rows)
    if out_csv is not None:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_csv, index=False)
    return out, summaries


def per_region_lookup(
    leaf_df: pd.DataFrame,
    band: str,
    region: str,
    patients: list[str],
) -> dict:
    """Single-cell lookup with hypergeometric stat (no filter)."""
    base = load_cohort_contacts(patients)
    sub = leaf_df[(leaf_df.patient.isin(patients))
                  & (leaf_df.band == band)
                  & (~leaf_df.region.isin(["Wm", "Unk"]))]
    K = len(sub)
    N_total = len(base)
    base_rate = K / N_total if N_total > 0 else 0.0
    n_trace = int((sub.region == region).sum())
    n_contacts = int((base.region == region).sum())
    n_pat_trace = int(sub[sub.region == region].patient.nunique())
    if n_contacts == 0:
        return {"band": band, "region": region, "K": K, "N_total": N_total,
                "base_rate": base_rate, "n_trace": n_trace, "n_contacts": 0,
                "n_pat_trace": n_pat_trace,
                "trace_rate": np.nan, "enrichment": np.nan, "p_hyper": np.nan,
                "bonf_pass": False}
    trace_rate = n_trace / n_contacts
    enrich = trace_rate / base_rate if base_rate > 0 else np.nan
    p = float(hypergeom.sf(n_trace - 1, N_total, K, n_contacts))
    return {"band": band, "region": region, "K": K, "N_total": N_total,
            "base_rate": base_rate, "n_trace": n_trace, "n_contacts": n_contacts,
            "n_pat_trace": n_pat_trace,
            "trace_rate": trace_rate, "enrichment": enrich, "p_hyper": p,
            "bonf_pass": p < BONF_THR}


def fmt_p(p: float) -> str:
    if pd.isna(p):
        return "n/a"
    if p < 1e-5:
        return f"{p:.2e}"
    if p < 0.001:
        return f"{p:.4f}"
    return f"{p:.3f}"


def header(s: str) -> None:
    print(f"\n{'=' * len(s)}\n{s}\n{'=' * len(s)}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    leaf = pd.read_csv(PER_LEAF_CSV)

    header(f"Bonferroni m = {BONF_M},  α/m = {BONF_THR:.5f}")
    header("Per-patient β + low-γ trace-leaf totals (calibrated, p95)")
    print(leaf.groupby(["band", "patient"]).size()
          .unstack(fill_value=0)
          .loc[["beta", "low_gamma"]].T)

    # ------------------------------------------------------------------
    # 0. Verify locked headline cells
    # ------------------------------------------------------------------
    header("Headline verification")
    cells = [
        ("low_gamma", "ctx-lh-fusiform", PATIENTS, "low_γ fusiform | full n=10"),
        ("low_gamma", "ctx-lh-fusiform",
         [p for p in PATIENTS if p != "Pat_03"], "low_γ fusiform | drop Pat_03"),
        ("beta", "Hip", PATIENTS, "β Hippocampus | full n=10"),
        ("beta", "Hip", [p for p in PATIENTS if p != "Pat_03"],
         "β Hippocampus | drop Pat_03"),
        ("beta", "ctx-lh-fusiform", PATIENTS, "β fusiform | full n=10"),
        ("beta", "ctx-lh-fusiform", [p for p in PATIENTS if p != "Pat_03"],
         "β fusiform | drop Pat_03"),
        ("beta", "ctx-lh-superiortemporal", PATIENTS,
         "β L sup temp | full n=10"),
        ("beta", "ctx-lh-superiortemporal",
         [p for p in PATIENTS if p != "Pat_03"], "β L sup temp | drop Pat_03"),
    ]
    head_rows = []
    for band, region, pats, lab in cells:
        info = per_region_lookup(leaf, band, region, pats)
        info["scenario"] = lab
        info["n_patients_subset"] = len(pats)
        head_rows.append(info)
        print(f"  {lab:38s}  enrichment={info['enrichment']!s:>8.6}  "
              f"p={fmt_p(info['p_hyper'])}  bonf={info['bonf_pass']}  "
              f"({info['n_trace']}/{info['n_contacts']} contacts, "
              f"{info['n_pat_trace']} pts)")
    head_df = pd.DataFrame(head_rows)
    head_df.to_csv(OUT_DIR / "headline_verification.csv", index=False)

    # ------------------------------------------------------------------
    # Test A: pro-cohort β enrichment (n=8)
    # ------------------------------------------------------------------
    header("Test A — Pro-cohort β enrichment (n=8, drop Pat_07 + Pat_15)")
    print(f"  patients = {PRO_COHORT}")
    pro_full, pro_summ = run_subset(leaf, "pro_cohort_n8", PRO_COHORT,
                                    bands=("beta",),
                                    out_csv=OUT_DIR / "testA_pro_cohort_beta.csv")
    print(f"  K = {pro_summ[0]['K']} trace-leaves; "
          f"N_total = {pro_summ[0]['N_total']} cortical contacts; "
          f"baseline = {pro_summ[0]['base_rate']*100:.2f}%")
    if pro_full.empty:
        print("  no enriched regions pass filter")
    else:
        for _, r in pro_full.iterrows():
            print(f"  {r.region:32s}  enrichment={r.enrichment:.2f}×  "
                  f"p={fmt_p(r.p_hyper)}  bonf={bool(r.bonf_pass)}  "
                  f"({int(r.n_trace)}/{int(r.n_contacts)}, "
                  f"{int(r.n_pat_trace)} pts)")
    # Specifically also lookup the headline regions in the pro-cohort
    print()
    for region in ["Hip", "ctx-lh-fusiform", "ctx-lh-superiortemporal",
                   "ctx-lh-middletemporal", "ctx-lh-lateralorbitofrontal"]:
        info = per_region_lookup(leaf, "beta", region, PRO_COHORT)
        print(f"  [explicit] {region:32s}  "
              f"enrich={info['enrichment']!s:>6}  "
              f"p={fmt_p(info['p_hyper'])}  bonf={info['bonf_pass']}  "
              f"({info['n_trace']}/{info['n_contacts']}, "
              f"{info['n_pat_trace']} pts)")

    # ------------------------------------------------------------------
    # Test B: Pat_02 dropout for low-γ ctx-lh-fusiform
    # ------------------------------------------------------------------
    header("Test B — Pat_02 dropout for low-γ ctx-lh-fusiform")
    no_pat02 = [p for p in PATIENTS if p != "Pat_02"]
    print(f"  patients = {no_pat02}")
    bres = per_region_lookup(leaf, "low_gamma", "ctx-lh-fusiform", no_pat02)
    pd.DataFrame([bres]).to_csv(OUT_DIR / "testB_low_gamma_fusiform_no_pat02.csv",
                                index=False)
    print(f"  K = {bres['K']}, N_total = {bres['N_total']}, "
          f"baseline = {bres['base_rate']*100:.2f}%")
    print(f"  ctx-lh-fusiform: {bres['n_trace']} trace / "
          f"{bres['n_contacts']} contacts → "
          f"trace-rate = {bres['trace_rate']*100:.2f}%")
    print(f"  enrichment = {bres['enrichment']:.2f}×, "
          f"p_hyper = {fmt_p(bres['p_hyper'])},  "
          f"Bonferroni m=48 pass = {bres['bonf_pass']}, "
          f"{bres['n_pat_trace']} contributing patients")

    # Also do the full cohort enrichment for low-γ without Pat_02 to
    # see if any other region surfaces.
    nopat02_full, _ = run_subset(leaf, "no_pat02", no_pat02,
                                 bands=("low_gamma",),
                                 out_csv=OUT_DIR / "testB_full_low_gamma_no_pat02.csv")
    print("\n  All low-γ enriched regions w/o Pat_02:")
    for _, r in nopat02_full.iterrows():
        print(f"    {r.region:32s}  enrich={r.enrichment:.2f}×  "
              f"p={fmt_p(r.p_hyper)}  bonf={bool(r.bonf_pass)}  "
              f"({int(r.n_trace)}/{int(r.n_contacts)}, "
              f"{int(r.n_pat_trace)} pts)")

    # ------------------------------------------------------------------
    # Test C: anti-aligned per-patient regional breakdown
    # ------------------------------------------------------------------
    header("Test C — Anti-aligned single-patient regional characterization")
    cohort_base_per_reg = (load_cohort_contacts(PATIENTS)
                           .groupby("region").size()
                           .rename("n_contacts_cohort").reset_index())
    rows_c = []
    for pat in ANTI_ALIGNED:
        print(f"\n  {pat} — β trace-leaves")
        sub = leaf[(leaf.patient == pat) & (leaf.band == "beta")]
        sub_cort = sub[~sub.region.isin(["Wm", "Unk"])]
        n_total_leaves = len(sub)
        n_cortical = len(sub_cort)
        print(f"    {n_total_leaves} total β trace-leaves "
              f"({n_cortical} cortical, {n_total_leaves - n_cortical} Wm/Unk)")
        # per-region distribution
        reg_count = (sub_cort.groupby("region")
                     .size().rename("n_trace").reset_index()
                     .sort_values("n_trace", ascending=False))
        # patient's own contact distribution
        own_base = load_cohort_contacts([pat])
        own_per_reg = (own_base.groupby("region").size()
                       .rename("n_contacts_self").reset_index())
        K = n_cortical
        N_total = len(own_base)
        base_rate = K / N_total if N_total > 0 else 0.0
        merged = reg_count.merge(own_per_reg, on="region", how="left")\
                          .merge(cohort_base_per_reg, on="region", how="left")
        merged["trace_rate_self"] = merged.n_trace / merged.n_contacts_self
        merged["enrichment_self"] = merged.trace_rate_self / base_rate
        merged["p_hyper_self"] = merged.apply(
            lambda r: float(hypergeom.sf(r.n_trace - 1, N_total, K,
                                         int(r.n_contacts_self)))
            if r.n_contacts_self > 0 else np.nan, axis=1)
        merged["patient"] = pat
        merged["base_rate_self"] = base_rate
        merged["K_self"] = K
        merged["N_total_self"] = N_total
        merged["bonf_pass"] = merged.p_hyper_self < BONF_THR
        for _, r in merged.iterrows():
            print(f"    {r.region:36s}  "
                  f"{int(r.n_trace)}/{int(r.n_contacts_self)} self  "
                  f"({int(r.n_contacts_cohort) if not pd.isna(r.n_contacts_cohort) else 0} cohort)  "
                  f"rate={r.trace_rate_self*100:.0f}%  "
                  f"enrich={r.enrichment_self:.2f}×  "
                  f"p={fmt_p(r.p_hyper_self)}  bonf={bool(r.bonf_pass)}")
            rows_c.append(r.to_dict())
        print(f"    self baseline = {base_rate*100:.2f}% "
              f"(K={K}, N_total={N_total} cortical contacts)")
    pd.DataFrame(rows_c).to_csv(OUT_DIR / "testC_anti_aligned_per_patient.csv",
                                index=False)

    print(f"\n[49] all tables written under {OUT_DIR}")


if __name__ == "__main__":
    main()
