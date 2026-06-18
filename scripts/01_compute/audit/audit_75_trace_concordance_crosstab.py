#!/usr/bin/env python3
"""Audit 75 — per-patient trace concordance cross-tab + blind-FC contrast.

SURFACING ONLY. No new statistic, no new null. This script aggregates
already-computed matched-strength surrogate outputs into a per-band,
per-patient concordance classification across the two LRG-layer probes
that carry the manuscript's universality claim, plus the blind raw-FC
edge baseline that is supposed to be unable to resolve the same richness.

Three sources, all matched-strength (R=200, 4-cycle +/- delta,
split-baseline Delta_task = X^tt - X^preA, Delta_rest = X^post - X^preB):

  rho_split^coph : data/audit/matched_strength_surrogate_split_baseline/
                       per_patient_per_band.csv                 (audit_63)
  rho_split^raw  : data/audit/raw_fc_matched_strength/
                       per_patient_per_band.csv                 (audit_67, BLIND)
  T_G*           : data/audit/grassmann_matched_strength_surrogate/
                       per_patient_per_band_per_k.csv           (audit_66 -> n_sig_k)
  T_G* gate/LOO  : data/audit/grassmann_cluster_extent/
                       cohort_summary.csv + loo_cluster_p_mass.csv

rho_split^coph and rho_split^raw are the SAME split-baseline Spearman
statistic under the SAME null; the only difference is the per-pair
representation (LRG cophenetic distance D_coph vs raw imcoh_abs edges
triu(A)). The contrast isolates what the hierarchy adds over flat FC.

Per-patient trace rules (verbatim from .agents/preprint/tables/
beta_per_patient.tex, the locked manuscript table):

  rho_split (coph & raw) : trace iff obs_z >= 1.96 vs own surrogate
  T_G*                   : trace iff n_sig_k >= 20, where n_sig_k counts
                           k-cells with obs_p_one_sided_upper < 0.05
                           (null expectation ~= 5.5 cells)

Integrity gate: the beta column recomputed here MUST reproduce
beta_per_patient.tex exactly (T_G* n_sig_k per patient + the 7/10
rho_split flags). HALT on mismatch.

Outputs (data/audit/trace_concordance_crosstab/):
  per_patient_flags.csv  - 60 rows, the three flags + LRG concordance class
  band_crosstab.csv      - 6 rows, both / only-one / neither tallies + blind
Prints the per-band single-patient breakdown to stdout. Provenance README.md
and the synthesis report are maintained alongside (not written by this script).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

AUDIT = ROOT / "data" / "audit"
COPH_CSV = AUDIT / "matched_strength_surrogate_split_baseline" / "per_patient_per_band.csv"
RAW_CSV = AUDIT / "raw_fc_matched_strength" / "per_patient_per_band.csv"
GRASS_PERK_CSV = AUDIT / "grassmann_matched_strength_surrogate" / "per_patient_per_band_per_k.csv"
GRASS_GATE_CSV = AUDIT / "grassmann_cluster_extent" / "cohort_summary.csv"
GRASS_LOO_CSV = AUDIT / "grassmann_cluster_extent" / "loo_cluster_p_mass.csv"

OUT = AUDIT / "trace_concordance_crosstab"

BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

Z_THRESH = 1.96          # rho_split per-patient significance (manuscript rule)
NSIGK_THRESH = 20        # T_G* per-patient significance (manuscript rule)
ALPHA_K = 0.05           # per-k cell threshold for n_sig_k

# beta column of beta_per_patient.tex (n_sig_k), for the integrity gate.
TEX_BETA_NSIGK = {
    "Pat_02": 108, "Pat_03": 96, "Pat_05": 87, "Pat_06": 96, "Pat_07": 62,
    "Pat_08": 62, "Pat_10": 71, "Pat_13": 61, "Pat_14": 78, "Pat_15": 1,
}
# beta rho_split^coph trace patients (z >= 1.96) from the same table = 7/10.
TEX_BETA_COPH_TRACE = {
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08", "Pat_13",
}


def load_rho_flags(csv_path: Path, label: str) -> pd.DataFrame:
    """Per-(patient, band) z and z>=1.96 trace flag for a rho_split table."""
    df = pd.read_csv(csv_path)
    out = df[["patient", "band", "obs_rho", "obs_z"]].copy()
    out = out.rename(columns={"obs_rho": f"{label}_rho", "obs_z": f"{label}_z"})
    out[f"{label}_trace"] = out[f"{label}_z"] >= Z_THRESH
    return out


def load_grassmann_nsigk() -> pd.DataFrame:
    """Aggregate per-(patient, band, k) Grassmann surrogate p into n_sig_k."""
    df = pd.read_csv(GRASS_PERK_CSV)
    df["sig"] = df["obs_p_one_sided_upper"] < ALPHA_K
    agg = (
        df.groupby(["patient", "band"])["sig"].sum().astype(int)
        .reset_index().rename(columns={"sig": "TG_nsigk"})
    )
    agg["TG_trace"] = agg["TG_nsigk"] >= NSIGK_THRESH
    return agg


def integrity_check(coph: pd.DataFrame, grass: pd.DataFrame) -> None:
    """Reproduce the locked beta_per_patient.tex column or HALT."""
    gb = grass[grass["band"] == "beta"].set_index("patient")["TG_nsigk"].to_dict()
    bad = {p: (gb.get(p), v) for p, v in TEX_BETA_NSIGK.items() if gb.get(p) != v}
    if bad:
        print("INTEGRITY FAILURE — beta T_G* n_sig_k does not match tex:", file=sys.stderr)
        for p, (got, want) in bad.items():
            print(f"  {p}: recomputed {got}, tex {want}", file=sys.stderr)
        sys.exit(1)

    cb = set(
        coph[(coph["band"] == "beta") & coph["coph_trace"]]["patient"].tolist()
    )
    if cb != TEX_BETA_COPH_TRACE:
        print("INTEGRITY FAILURE — beta rho_split^coph trace set != tex:", file=sys.stderr)
        print(f"  recomputed {sorted(cb)}", file=sys.stderr)
        print(f"  tex        {sorted(TEX_BETA_COPH_TRACE)}", file=sys.stderr)
        sys.exit(1)
    print("[integrity] beta column reproduces beta_per_patient.tex exactly "
          f"(T_G* n_sig_k 10/10 cells, rho_split^coph {len(cb)}/10).")


def classify(row: pd.Series) -> str:
    c, g = bool(row["coph_trace"]), bool(row["TG_trace"])
    if c and g:
        return "both"
    if c and not g:
        return "only_coph"
    if g and not c:
        return "only_TG"
    return "neither"


def loo_wilcoxon_p(obs: np.ndarray, surr_med: np.ndarray) -> float:
    """Max one-sided paired Wilcoxon p over leave-one-patient-out (descriptive)."""
    diffs = obs - surr_med
    worst = 0.0
    for i in range(len(diffs)):
        d = np.delete(diffs, i)
        try:
            p = wilcoxon(d, alternative="greater", zero_method="wilcox").pvalue
        except ValueError:
            p = np.nan
        if np.isfinite(p):
            worst = max(worst, p)
    return worst


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    coph = load_rho_flags(COPH_CSV, "coph")
    raw = load_rho_flags(RAW_CSV, "raw")
    grass = load_grassmann_nsigk()

    integrity_check(coph, grass)

    flags = (
        coph.merge(grass, on=["patient", "band"], how="inner")
        .merge(raw, on=["patient", "band"], how="inner")
    )
    flags["lrg_class"] = flags.apply(classify, axis=1)
    flags["lrg_any"] = flags["coph_trace"] | flags["TG_trace"]

    flags["band"] = pd.Categorical(flags["band"], categories=BAND_ORDER, ordered=True)
    flags = flags.sort_values(["band", "patient"]).reset_index(drop=True)
    flags.to_csv(OUT / "per_patient_flags.csv", index=False)

    # --- per-band cross-tab tallies + cohort gates -------------------------
    gate = pd.read_csv(GRASS_GATE_CSV).set_index("band")
    loo = pd.read_csv(GRASS_LOO_CSV)
    coph_full = pd.read_csv(COPH_CSV)
    raw_full = pd.read_csv(RAW_CSV)

    rows = []
    for band in BAND_ORDER:
        b = flags[flags["band"] == band]
        n_both = int((b["lrg_class"] == "both").sum())
        n_only_c = int((b["lrg_class"] == "only_coph").sum())
        n_only_g = int((b["lrg_class"] == "only_TG").sum())
        n_neither = int((b["lrg_class"] == "neither").sum())
        n_any = int(b["lrg_any"].sum())
        n_blind = int(b["raw_trace"].sum())

        cb = coph_full[coph_full["band"] == band]
        coph_loo = loo_wilcoxon_p(cb["obs_rho"].to_numpy(), cb["surr_p50"].to_numpy())
        rb = raw_full[raw_full["band"] == band]
        raw_loo = loo_wilcoxon_p(rb["obs_rho"].to_numpy(), rb["surr_p50"].to_numpy())

        g_p = float(gate.loc[band, "cluster_p_cluster_mass"])
        g_loo = float(gate.loc[band, "cluster_p_mass_loo_max"])
        g_loo_arg = gate.loc[band, "cluster_p_mass_loo_argmax_patient"]
        g_verdict = gate.loc[band, "verdict_cluster_extent"]

        rows.append(dict(
            band=band,
            n_both=n_both, n_only_coph=n_only_c, n_only_TG=n_only_g,
            n_neither=n_neither, n_lrg_any=n_any, n_blind=n_blind,
            coph_nsig=int(cb["obs_z"].ge(Z_THRESH).sum()),
            coph_loo_max_p=round(coph_loo, 4),
            TG_nsig=int(b["TG_trace"].sum()),
            TG_clusterp=g_p, TG_loo_max_p=g_loo, TG_loo_argmax=g_loo_arg,
            TG_verdict=g_verdict,
            blind_loo_max_p=round(raw_loo, 4),
        ))
    crosstab = pd.DataFrame(rows)
    crosstab.to_csv(OUT / "band_crosstab.csv", index=False)

    # --- stdout report -----------------------------------------------------
    pats = sorted(flags["patient"].unique())
    sym = {"both": "B", "only_coph": "C", "only_TG": "G", "neither": "."}
    print("\n================ PER-PATIENT CONCORDANCE (LRG layer) ================")
    print("  B=both probes  C=only rho_split^coph  G=only T_G*  .=neither")
    print("  (blind = raw-FC edge rho_split^raw, z>=1.96)\n")
    header = "band        " + " ".join(p.replace("Pat_", "") for p in pats)
    print(header)
    for band in BAND_ORDER:
        b = flags[flags["band"] == band].set_index("patient")
        cells = " ".join(f"{sym[b.loc[p,'lrg_class']]:>2}" for p in pats)
        print(f"{band:<11} {cells}")
    print("\n  blind raw-FC trace (z>=1.96):")
    for band in BAND_ORDER:
        b = flags[flags["band"] == band].set_index("patient")
        cells = " ".join(f"{'x' if b.loc[p,'raw_trace'] else '.':>2}" for p in pats)
        print(f"{band:<11} {cells}")

    print("\n================ BAND CROSS-TAB ================\n")
    cols = ["band", "n_both", "n_only_coph", "n_only_TG", "n_neither",
            "n_lrg_any", "n_blind", "coph_nsig", "TG_nsig",
            "TG_clusterp", "TG_loo_max_p", "coph_loo_max_p", "blind_loo_max_p"]
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(crosstab[cols].to_string(index=False))

    print("\nLegend: n_both/only/neither over the two LRG probes; "
          "n_lrg_any = >=1 LRG probe; n_blind = raw-FC edges.")
    print("coph_nsig = rho_split^coph patients z>=1.96; TG_nsig = T_G* patients n_sig_k>=20.")
    print(f"\nWrote:\n  {OUT/'per_patient_flags.csv'}\n  {OUT/'band_crosstab.csv'}")


if __name__ == "__main__":
    main()
