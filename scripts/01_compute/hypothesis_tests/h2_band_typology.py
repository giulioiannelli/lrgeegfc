#!/usr/bin/env python3
"""Band typology: classify each frequency band into three regimes.

Three quantities per (patient, band):
  1. M_task  = mean VI(rest_pre, task_test) across k
               → "how much did task reorganize the structure vs rest_pre?"
  2. rho     = Spearman ρ(Δ_rest, Δ_task) on ultrametric distances (H2c)
               → "did rest_post drift in the same direction as task?"
  3. Δρ      = mean Δρ across k (H2d)
               → "do specific task-induced pair co-activations persist?"

Band-level summary (mean across patients + rank). Classification:

    TRACE BAND   : M_task high, Δρ high   → task reorganizes AND rest_post inherits
    ERGODIC BAND : M_task high, Δρ low    → task reorganizes but rest_post erases
    DEAD BAND    : M_task low,  Δρ low    → task didn't reorganize in the first place
    (residual)   : M_task low,  Δρ high   → unlikely in practice; flag for inspection

"high"/"low" are relative to the inter-band median of that metric.

Outputs:
    data/reports/imcoh_vi/h2_band_typology.md
    data/reports/imcoh_vi/h2_band_typology.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT

OUT_DIR = REPORTS_ROOT / "imcoh_vi"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT


# ─────────────────────────── loaders ───────────────────────────

def load_m_task() -> pd.DataFrame:
    """k-averaged VI(rest_pre, task_test) per (patient, band)."""
    df = pd.read_csv(OUT_DIR / "vi_raw_profiles.csv")
    mask = (df["phase_a"] == "rest_pre") & (df["phase_b"] == "task_test")
    sub = df[mask].copy()
    wide = (
        sub.groupby(["patient", "band"])["vi"].mean()
        .unstack("band").reindex(columns=BANDS)
    )
    return wide


def load_rho() -> pd.DataFrame:
    """H2c Spearman ρ (target=task_test) per (patient, band)."""
    df = pd.read_csv(OUT_DIR / "h2c_ultrametric_drift_raw.csv")
    df = df[df["has_rpre"] & df["has_rpost"] & df["has_ttest"]].copy()
    wide = (
        df.pivot(index="patient", columns="band", values="rho_task")
        .reindex(columns=BANDS)
    )
    return wide


def load_delta_rho() -> pd.DataFrame:
    """H2d k-averaged Δρ per (patient, band)."""
    df = pd.read_csv(OUT_DIR / "h2d_persistence_raw.csv")
    wide = (
        df.groupby(["patient", "band"])["delta_rho"].mean()
        .unstack("band").reindex(columns=BANDS)
    )
    return wide


# ─────────────────────────── classification ───────────────────────────

def classify_band(m_task: float, delta_rho: float,
                  m_task_median: float, delta_rho_median: float) -> str:
    hi_mt = m_task >= m_task_median
    hi_dr = delta_rho >= delta_rho_median
    if hi_mt and hi_dr:
        return "TRACE"
    if hi_mt and not hi_dr:
        return "ERGODIC"
    if not hi_mt and not hi_dr:
        return "DEAD"
    return "atypical (low-M, high-Δρ)"


# ─────────────────────────── main ───────────────────────────

def main():
    print("Loading per-patient × band tables...")
    m_task_pat = load_m_task()
    rho_pat = load_rho()
    drho_pat = load_delta_rho()

    # Per-patient set common to all three metrics
    common = m_task_pat.index.intersection(rho_pat.index).intersection(drho_pat.index)
    print(f"Common patients (all 3 metrics defined): {len(common)} → {list(common)}")

    m_task_pat = m_task_pat.loc[common]
    rho_pat = rho_pat.loc[common]
    drho_pat = drho_pat.loc[common]

    # Band-level aggregation (mean across patients)
    summary = pd.DataFrame({
        "M_task (mean VI rpre-ttest)": m_task_pat.mean(axis=0),
        "ρ  (H2c direction)":           rho_pat.mean(axis=0),
        "Δρ (H2d persistence)":         drho_pat.mean(axis=0),
    }).reindex(BANDS)

    # Per-metric ranks (1 = lowest)
    ranks = pd.DataFrame({
        f"rank_{c}": summary[c].rank(method="min").astype(int)
        for c in summary.columns
    }).reindex(BANDS)

    # Classification relative to inter-band medians
    mt_med = summary["M_task (mean VI rpre-ttest)"].median()
    dr_med = summary["Δρ (H2d persistence)"].median()

    cls = {
        b: classify_band(summary.loc[b, "M_task (mean VI rpre-ttest)"],
                         summary.loc[b, "Δρ (H2d persistence)"],
                         mt_med, dr_med)
        for b in BANDS
    }

    # ─── write CSV ───────────────────────────────────────────
    out_csv = summary.join(ranks)
    out_csv["classification"] = pd.Series(cls)
    out_csv.index.name = "band"
    out_csv.to_csv(OUT_DIR / "h2_band_typology.csv")
    print(f"Wrote {OUT_DIR / 'h2_band_typology.csv'}")

    # ─── write MD report ─────────────────────────────────────
    md: list[str] = []
    ap = md.append

    ap("# H2 band typology — three-way classification")
    ap("")
    ap(f"n patients (common to M_task, ρ, Δρ) = **{len(common)}**: {', '.join(common)}")
    ap("")
    ap("Three quantities per (patient, band):")
    ap("- **M_task** = mean VI(rest_pre, task_test) across k. *How much did task reorganize vs rest_pre?*")
    ap("- **ρ** = Spearman ρ on ultrametric drift vectors (H2c). *Did rest_post drift same direction as task?*")
    ap("- **Δρ** = k-averaged block-persistence excess (H2d). *Do task-induced pair clusters persist at rest_post?*")
    ap("")
    ap("**Classification** (relative to inter-band median):")
    ap("- **TRACE** — M_task ≥ median *and* Δρ ≥ median → task reorganized, rest_post inherited")
    ap("- **ERGODIC** — M_task ≥ median *and* Δρ < median → task reorganized, rest_post erased")
    ap("- **DEAD** — M_task < median *and* Δρ < median → task didn't reorganize to begin with")
    ap("")
    ap(f"Inter-band medians: M_task = {mt_med:.3f}, Δρ = {dr_med:.3f}")
    ap("")

    # Main band table
    ap("## Band-level summary")
    ap("")
    ap("| band | M_task | rank(M) | ρ | rank(ρ) | Δρ | rank(Δρ) | class |")
    ap("|:-----|------:|-------:|--:|-------:|---:|--------:|:------|")
    for b in BANDS:
        ap(f"| {TEX[b]} | {summary.loc[b, 'M_task (mean VI rpre-ttest)']:.3f} "
           f"| {ranks.loc[b, 'rank_M_task (mean VI rpre-ttest)']} "
           f"| {summary.loc[b, 'ρ  (H2c direction)']:+.3f} "
           f"| {ranks.loc[b, 'rank_ρ  (H2c direction)']} "
           f"| {summary.loc[b, 'Δρ (H2d persistence)']:+.3f} "
           f"| {ranks.loc[b, 'rank_Δρ (H2d persistence)']} "
           f"| **{cls[b]}** |")
    ap("")

    # Per-patient tables
    def pp_table(label: str, df: pd.DataFrame, fmt: str) -> None:
        ap(f"### {label} (per patient)")
        ap("")
        ap("| patient | " + " | ".join(TEX[b] for b in BANDS) + " |")
        ap("|:--------|" + "-----:|" * len(BANDS))
        for pat in df.index:
            row = df.loc[pat]
            ap(f"| {pat} | " + " | ".join(
                fmt.format(row[b]) if np.isfinite(row[b]) else "—"
                for b in BANDS
            ) + " |")
        ap("")

    pp_table("M_task = mean VI(rest_pre, task_test) across k",
             m_task_pat, "{:.3f}")
    pp_table("ρ (H2c) per patient × band", rho_pat, "{:+.3f}")
    pp_table("Δρ (H2d) k-averaged per patient × band", drho_pat, "{:+.3f}")

    md_path = OUT_DIR / "h2_band_typology.md"
    md_path.write_text("\n".join(md))
    print(f"Wrote {md_path}")

    # ─── console summary ─────────────────────────────────────
    print("\n" + "=" * 60)
    print("BAND TYPOLOGY (mean across patients)")
    print("=" * 60)
    print(f"{'band':<12} {'M_task':>8} {'ρ':>8} {'Δρ':>8}  class")
    print("-" * 60)
    for b in BANDS:
        print(f"{b:<12} "
              f"{summary.loc[b, 'M_task (mean VI rpre-ttest)']:>8.3f} "
              f"{summary.loc[b, 'ρ  (H2c direction)']:>8.3f} "
              f"{summary.loc[b, 'Δρ (H2d persistence)']:>8.3f}  "
              f"{cls[b]}")
    print("=" * 60)


if __name__ == "__main__":
    main()
