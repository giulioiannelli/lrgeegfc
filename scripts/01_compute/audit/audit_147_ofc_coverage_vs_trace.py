#!/usr/bin/env python3
"""audit_147 — does per-patient OFC electrode coverage explain who-traces?

The Q1 follow-up (`.agents/reports/2026-06-25_multiphase-snr-reliability.md`):
the between-patient trace heterogeneity is NOT pure measurement reliability
(reliability explains ~6%; Pat_15 is the highest-reliability, lowest-trace
patient). Proposed explanation: the β trace LIVES in OFC (locked audit_83
β→OFC, LOO-robust), so a patient must IMPLANT OFC to show it; Pat_15 is the
right-hemisphere-only implant and would under-sample the bilateral OFC hub.
This script TESTS that, instead of asserting it.

5-point preamble
1. Claim: per-patient β `ρ_split` tracks OFC electrode coverage — patients who
   sample OFC trace, patients who don't (Pat_15) don't.
2. Null/baseline: no association (Spearman ρ≈0); coverage is irrelevant and the
   residual biology is something else.
3. Strongest alternative: the association is just node count / graph size (more
   contacts anywhere → bigger ρ). Controlled by also correlating ρ_split with
   TOTAL N and with non-OFC coverage, and by using OFC FRACTION (size-normalised).
4. Cannot: n=10, so this is a descriptive association with a wide CI, not a
   powered test; it cannot prove causation (coverage vs some correlated
   recording/biology factor); OFC system label granularity is the atlas's.
5. Falsify: if Spearman(ρ_split_β, OFC_frac) ≈ 0 or negative, or if total-N
   correlates as strongly as OFC, the coverage explanation is rejected and the
   residual biology stays unexplained.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
GATE = ROOT / "data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv"
REL = ROOT / "data/audit/multiphase_snr/per_patient_per_band_reliability.csv"
OUT = ROOT / "data/audit/ofc_coverage"


def coverage(pat, system):
    rdf = load_channel_regions(pat)
    sysv = rdf["system"].astype(str)
    n = len(rdf)
    cnt = int((sysv == system).sum())
    # left/right split if a hemisphere column exists
    hemi = None
    for col in ("hemisphere", "hemi"):
        if col in rdf.columns:
            hemi = rdf[col].astype(str)
            break
    lr = ""
    if hemi is not None:
        m = (sysv == system)
        ls = sorted(hemi[m].str[:1].str.upper().unique().tolist())
        lr = "/".join(ls)
    return n, cnt, cnt / n if n else np.nan, lr


def main():
    gate = pd.read_csv(GATE)
    rel = pd.read_csv(REL) if REL.exists() else None
    rows = []
    for pat in COHORT:
        n, ofc_n, ofc_frac, ofc_lr = coverage(pat, "OFC")
        _, pfc_n, pfc_frac, _ = coverage(pat, "PFC")
        b = gate[(gate.patient == pat) & (gate.band == "beta")]
        rho_b = float(b["obs_rho"].iloc[0]) if len(b) else np.nan
        relc = np.nan
        if rel is not None:
            rb = rel[(rel.patient == pat) & (rel.band == "beta")]
            if len(rb) and "rel_composite" in rb:
                relc = float(rb["rel_composite"].iloc[0])
        rows.append(dict(patient=pat, N=n, rho_beta=rho_b,
                         OFC_n=ofc_n, OFC_frac=ofc_frac, OFC_hemi=ofc_lr,
                         PFC_n=pfc_n, rel_composite_beta=relc))
    df = pd.DataFrame(rows).sort_values("rho_beta", ascending=False)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "ofc_coverage_vs_trace.csv", index=False)

    print(df.round(3).to_string(index=False))
    print("\n=== associations across the 10 patients ===")
    def sp(x, y):
        m = df[[x, y]].dropna()
        if len(m) < 3:
            return np.nan, np.nan
        r, p = spearmanr(m[x], m[y])
        return r, p
    for x in ["OFC_n", "OFC_frac", "PFC_n", "N"]:
        r, p = sp(x, "rho_beta")
        tag = "  <-- key" if x.startswith("OFC") else ("  (size-confound control)" if x == "N" else "")
        print(f"  Spearman(rho_beta, {x:9s}) = {r:+.3f}  (p={p:.3f}){tag}")
    # does OFC coverage explain beyond reliability? rank-partial via residuals
    m = df[["rho_beta", "OFC_frac", "rel_composite_beta"]].dropna()
    if len(m) >= 4:
        from scipy.stats import rankdata
        def resid(a, b):
            ra, rb = rankdata(a), rankdata(b)
            beta = np.polyfit(rb, ra, 1)
            return ra - np.polyval(beta, rb)
        rr = resid(m["rho_beta"], m["rel_composite_beta"])
        ro = resid(m["OFC_frac"], m["rel_composite_beta"])
        pr, pp = spearmanr(rr, ro)
        print(f"\n  partial Spearman(rho_beta, OFC_frac | reliability) = {pr:+.3f} (p={pp:.3f})")
    print("\n=== Pat_15 (the falsifier patient) ===")
    p15 = df[df.patient == "Pat_15"].iloc[0]
    print(f"  Pat_15: rho_beta={p15.rho_beta:+.3f}  OFC_n={int(p15.OFC_n)} "
          f"OFC_frac={p15.OFC_frac:.3f} hemi=[{p15.OFC_hemi}]  rel={p15.rel_composite_beta:.3f}")
    print(f"\nOutputs -> {OUT}/ofc_coverage_vs_trace.csv")


if __name__ == "__main__":
    main()
