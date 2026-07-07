#!/usr/bin/env python3
"""Audit 146 — does the MULTISCALE (cophenetic) trace BEAT the RAW per-edge comparison?

Head: the LRG cophenetic ``rho_split`` and a raw per-edge ``rho`` are computed
from the IDENTICAL split-baseline construction (same patients, bands, rest_pre
halves); the only difference is cophenetic-distance-space vs raw-edge-space.
This script runs the head-to-head on three axes (magnitude, band-taxonomy story,
SNR-gating) so the multiscale tool's contribution can be claimed honestly or
withdrawn.

================================ 5-POINT CRITICAL PREAMBLE ====================
(1) CLAIM. The LRG cophenetic trace (rho_split) surfaces cross-phase
    persistence that a plain raw per-edge comparison (raw_rho) does NOT —
    i.e. the multiscale/Laplacian machinery adds value beyond a coordinate
    re-labelling of the same edges. Concretely, on the SAME patients we expect
    at least one of: (a) cophenetic rho_split clears its matched-strength null
    where raw does not; (b) the 3-tier consistency taxonomy
    [consistent beta/low_gamma | patient-specific delta/alpha/high_gamma |
    absent theta] is sharper/cleaner in cophenetic-space than in raw-space;
    (c) the SNR->trace gating (high-SNR patients carry the trace, low-SNR do
    not) is multiscale-specific or sharper.

(2) NULL / BASELINE. The baseline IS the raw measure. raw_rho =
    Spearman(triu(A_tt) - triu(A_preA), triu(A_post) - triu(A_preB)) on raw
    imcoh_abs edges. This is the COMPARISON BASELINE mandated by
    feedback_results_only_in_laplacian_framework: raw FC is never a result, it
    is the foil the multiscale tool must beat. Each measure ALSO has its own
    R=200 4-cycle +/-delta matched-strength surrogate null already cached
    (audit_63 for cophenetic, audit_67 for raw) so "clears a null" is
    comparable across the two spaces.

(3) STRONGEST ALTERNATIVE. The strongest alternative to the claim is:
    "rho_split is just raw_rho in a monotone-compressed coordinate system; the
    ultrametric/average-linkage step adds nothing — whatever the cophenetic
    trace says, the raw edges already said it, possibly with a LARGER effect
    size because cophenetic distances are quantized to N-1 merge heights and
    lose magnitude." If raw reproduces the band taxonomy AND the SNR-gating AND
    has >= the magnitude, the LRG machinery is a redundant re-labelling and the
    multiscale claim must be withdrawn.

(4) DOES THE COMPARISON CONTROL FOR IT (mechanism). Yes, by construction:
    every factor except the distance-space is held fixed. Same task_test FC,
    same rest_post FC, same two rest_pre half-FCs, same Spearman estimator,
    same N(N-1)/2 pair support feeding both deltas, same cohort. The ONLY
    difference is whether the two phase-difference vectors live on raw edges
    (raw_rho) or on canonical_cophenet(W) ultrametric distances (rho_split).
    Therefore any divergence in null-clearing, band pattern, or SNR-separation
    is attributable to the cophenetic transform alone -- the alternative in (3)
    is exactly what a paired raw-vs-cophenetic contrast on identical inputs
    tests. What it CANNOT do: it cannot show the cophenetic trace is
    *biologically* superior, only that it is *statistically* non-redundant on
    these three axes; and it cannot adjudicate the DOWNSTREAM multiscale claims
    (OFC localization, encoding/inference dissociation) which are not edge
    statistics at all -- those live or die on their own audits. A raw per-edge
    rho is the hardest like-for-like foil for the *trace magnitude* claim; it
    is not a foil for the *graph-structure* (community/path/hierarchy) claims,
    which raw edges cannot express by definition.

(5) WHAT WOULD FALSIFY "multiscale beats raw". (i) raw_rho clears its
    matched-strength null in every band where cophenetic does, AND with effect
    size >= cophenetic (magnitude axis lost); (ii) the per-band cohort pattern
    of raw_rho reproduces the 3-tier consistency taxonomy as cleanly or more
    cleanly than cophenetic (story axis lost); (iii) raw_rho separates the
    high-SNR-5 from the low-SNR-5 as sharply as rho_split (gating axis lost).
    If all three hold, the honest verdict is: the multiscale tool is redundant
    with raw on the trace axis and cannot be claimed as the contribution -- a
    serious finding, reported plainly, not hidden.
==============================================================================

Inputs (all cached, both produced by the SAME split-baseline pipeline):
- MULTISCALE rho_split + its matched-strength null:
    data/audit/matched_strength_surrogate_split_baseline/
        per_patient_per_band.csv   (col obs_rho, surr_p50/p95, obs_p_one_sided)
        cohort_summary.csv
- RAW per-edge rho + its matched-strength null (audit_67, identical construction):
    data/audit/raw_fc_matched_strength/
        per_patient_per_band.csv   (col obs_rho, surr_p50/p95, obs_p_one_sided)
        cohort_summary_all_bands.csv
- Half-FCs (for the independent recompute spot-check + SNR proxy):
    data/cache/imcoh_halves_fc/{pat}/{band}_rest_pre_{A|B}_imcoh_abs.npy

The two CSVs are node-matched (verified: 60 cells, 0 N_nodes mismatch); their
obs_rho columns reproduce an independent recompute bit-exactly (diffs ~1e-17).
This script reads them, recomputes a few cells to certify, runs the three axes,
and writes the deliverables. It does NOT re-run any surrogate (the matched-
strength nulls are cached and stable).

Output: data/audit/raw_vs_multiscale/
    per_patient_per_band.csv   (patient, band, rho_split_multiscale, raw_rho,
                                delta, snr_tier, + null flags)
    cohort_summary.csv         (per-band paired stats, full cohort + hiSNR-5)
    snr_gating.csv             (per-band high-vs-low separation, both measures)
    recompute_certificate.csv  (spot-check vs cached)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon, mannwhitneyu

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.metrics.node_localization import canonical_cophenet
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial
from lrg_eegfc.workflow.fc import load_fc_matrix

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
# Published per-patient SNR ranking (task-vs-baseline detectability).
HI_SNR = ["Pat_06", "Pat_05", "Pat_02", "Pat_03", "Pat_08"]
LO_SNR = ["Pat_15", "Pat_14", "Pat_07", "Pat_13", "Pat_10"]
BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
# 3-tier multiscale consistency taxonomy (the "story" we test raw against).
TAXONOMY = {
    "delta": "patient-specific", "theta": "absent",
    "alpha": "patient-specific", "beta": "consistent",
    "low_gamma": "consistent", "high_gamma": "patient-specific",
}

MS_CSV = (ROOT / "data" / "audit"
          / "matched_strength_surrogate_split_baseline"
          / "per_patient_per_band.csv")
RAW_CSV = (ROOT / "data" / "audit" / "raw_fc_matched_strength"
           / "per_patient_per_band.csv")
HALVES = CACHE_ROOT / "imcoh_halves_fc"

OUT = ROOT / "data" / "audit" / "raw_vs_multiscale"
OUT.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------- #
# Prep helper (match the canonical pipeline exactly)
# --------------------------------------------------------------------------- #
def _prep(W: np.ndarray) -> np.ndarray:
    W = np.asarray(W, dtype=np.float64)
    W = 0.5 * (W + W.T)
    np.fill_diagonal(W, 0.0)
    return np.clip(W, 0.0, 1.0)


def _recompute_cell(pat: str, band: str) -> tuple[float, float]:
    """Independent recompute of (raw_rho, rho_split) from the half-FC cache."""
    A_tt = _prep(load_fc_matrix(pat, "task_test", band, "imcoh_abs"))
    A_post = _prep(load_fc_matrix(pat, "rest_post", band, "imcoh_abs"))
    A_preA = _prep(np.load(HALVES / pat / f"{band}_rest_pre_A_imcoh_abs.npy"))
    A_preB = _prep(np.load(HALVES / pat / f"{band}_rest_pre_B_imcoh_abs.npy"))
    iu = np.triu_indices(A_tt.shape[0], k=1)
    raw_rho = spearmanr(A_tt[iu] - A_preA[iu],
                        A_post[iu] - A_preB[iu]).statistic
    d_tt, d_post = canonical_cophenet(A_tt), canonical_cophenet(A_post)
    d_preA, d_preB = canonical_cophenet(A_preA), canonical_cophenet(A_preB)
    rho_split = spearmanr(d_tt - d_preA, d_post - d_preB).statistic
    return float(raw_rho), float(rho_split)


# --------------------------------------------------------------------------- #
# Load + certify the cached pair
# --------------------------------------------------------------------------- #
def load_paired() -> pd.DataFrame:
    ms = pd.read_csv(MS_CSV)
    rw = pd.read_csv(RAW_CSV)
    keep = ["patient", "band", "N_nodes", "obs_rho",
            "surr_p50", "surr_p95", "obs_p_one_sided"]
    ms = ms[keep].rename(columns={
        "obs_rho": "rho_split_multiscale",
        "surr_p50": "ms_surr_p50", "surr_p95": "ms_surr_p95",
        "obs_p_one_sided": "ms_p_one_sided"})
    rw = rw[keep].rename(columns={
        "obs_rho": "raw_rho",
        "surr_p50": "raw_surr_p50", "surr_p95": "raw_surr_p95",
        "obs_p_one_sided": "raw_p_one_sided"})
    df = ms.merge(rw, on=["patient", "band", "N_nodes"], how="inner")
    assert len(df) == 60, f"expected 60 paired cells, got {len(df)}"
    df["delta"] = df["rho_split_multiscale"] - df["raw_rho"]
    df["snr_tier"] = np.where(df.patient.isin(HI_SNR), "high", "low")
    # per-measure: does the OBSERVED clear its OWN matched-strength p95?
    df["ms_clears_null"] = df.rho_split_multiscale > df.ms_surr_p95
    df["raw_clears_null"] = df.raw_rho > df.raw_surr_p95
    return df


def certify(df: pd.DataFrame) -> pd.DataFrame:
    """Recompute a few cells; assert bit-exactness vs cached obs_rho."""
    rows = []
    for pat, band in [("Pat_06", "beta"), ("Pat_02", "alpha"),
                      ("Pat_08", "delta"), ("Pat_05", "low_gamma"),
                      ("Pat_03", "high_gamma")]:
        raw_rho, rho_split = _recompute_cell(pat, band)
        c = df[(df.patient == pat) & (df.band == band)].iloc[0]
        rows.append(dict(
            patient=pat, band=band,
            raw_recompute=raw_rho, raw_cached=c.raw_rho,
            raw_absdiff=abs(raw_rho - c.raw_rho),
            ms_recompute=rho_split, ms_cached=c.rho_split_multiscale,
            ms_absdiff=abs(rho_split - c.rho_split_multiscale)))
    cert = pd.DataFrame(rows)
    max_diff = max(cert.raw_absdiff.max(), cert.ms_absdiff.max())
    assert max_diff < 1e-9, f"recompute drift {max_diff:.2e} > 1e-9 — STOP"
    print(f"[audit_146] recompute certificate OK (max |diff| = {max_diff:.2e})")
    return cert


# --------------------------------------------------------------------------- #
# AXIS 1 — magnitude: does cophenetic persist MORE than raw edges?
# --------------------------------------------------------------------------- #
def axis1_magnitude(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort_name, pats in [("full10", COHORT), ("hiSNR5", HI_SNR)]:
        sub = df[df.patient.isin(pats)]
        for band in BAND_ORDER:
            b = sub[sub.band == band]
            d = b["delta"].to_numpy()              # multiscale - raw, paired
            ms = b["rho_split_multiscale"].to_numpy()
            raw = b["raw_rho"].to_numpy()
            n = len(d)
            # Two-sided paired Wilcoxon on the difference (is there ANY gap?).
            try:
                w2 = wilcoxon(ms, raw, alternative="two-sided",
                              method="approx")
                p_two = float(w2.pvalue)
            except ValueError:
                p_two = np.nan
            # One-sided: is multiscale GREATER than raw? (our hoped direction)
            z_gt, p_gt = wilcoxon_z(d)            # H1: delta>0
            rb = rank_biserial(d)
            rows.append(dict(
                cohort=cohort_name, band=band, n=n,
                median_multiscale=float(np.median(ms)),
                median_raw=float(np.median(raw)),
                median_delta=float(np.median(d)),
                mean_delta=float(np.mean(d)),
                n_multiscale_gt_raw=int((d > 0).sum()),
                wilcoxon_z_delta_gt0=z_gt,
                p_delta_gt0_one_sided=p_gt,
                p_delta_two_sided=p_two,
                rank_biserial_delta=rb))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# AXIS 2 — story: is the band taxonomy multiscale-specific or in raw too?
# --------------------------------------------------------------------------- #
def axis2_story(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for band in BAND_ORDER:
        b = df[df.band == band]
        ms = b["rho_split_multiscale"].to_numpy()
        raw = b["raw_rho"].to_numpy()
        rows.append(dict(
            band=band,
            taxonomy_multiscale=TAXONOMY[band],
            ms_median=float(np.median(ms)),
            ms_iqr=float(np.subtract(*np.percentile(ms, [75, 25]))),
            ms_n_pos=int((ms > 0).sum()),
            ms_n_clears_null=int(b["ms_clears_null"].sum()),
            raw_median=float(np.median(raw)),
            raw_iqr=float(np.subtract(*np.percentile(raw, [75, 25]))),
            raw_n_pos=int((raw > 0).sum()),
            raw_n_clears_null=int(b["raw_clears_null"].sum())))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# AXIS 3 — SNR-gating: is high-vs-low separation a raw or multiscale property?
# --------------------------------------------------------------------------- #
def axis3_snr(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for band in BAND_ORDER:
        b = df[df.band == band]
        out = dict(band=band)
        for label, col in [("multiscale", "rho_split_multiscale"),
                           ("raw", "raw_rho")]:
            hi = b[b.snr_tier == "high"][col].to_numpy()
            lo = b[b.snr_tier == "low"][col].to_numpy()
            # Mann-Whitney one-sided: high > low (detectability gating).
            try:
                u = mannwhitneyu(hi, lo, alternative="greater")
                p_mw = float(u.pvalue)
                # rank-biserial effect size from U.
                rb = 1.0 - 2.0 * u.statistic / (len(hi) * len(lo))
                rb = -rb  # so + means high>low
            except ValueError:
                p_mw, rb = np.nan, np.nan
            out[f"{label}_hi_median"] = float(np.median(hi))
            out[f"{label}_lo_median"] = float(np.median(lo))
            out[f"{label}_hi_minus_lo"] = float(np.median(hi) - np.median(lo))
            out[f"{label}_mw_p_hi_gt_lo"] = p_mw
            out[f"{label}_mw_rankbiserial"] = rb
        rows.append(out)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    df = load_paired()
    cert = certify(df)

    a1 = axis1_magnitude(df)
    a2 = axis2_story(df)
    a3 = axis3_snr(df)

    # per-patient-per-band table (the headline deliverable)
    cols = ["patient", "band", "N_nodes", "rho_split_multiscale", "raw_rho",
            "delta", "snr_tier", "ms_clears_null", "raw_clears_null",
            "ms_p_one_sided", "raw_p_one_sided",
            "ms_surr_p95", "raw_surr_p95"]
    df_out = df[cols].copy()
    df_out["band"] = pd.Categorical(df_out["band"], BAND_ORDER, ordered=True)
    df_out = df_out.sort_values(["band", "patient"]).reset_index(drop=True)

    df_out.to_csv(OUT / "per_patient_per_band.csv", index=False)
    a1.to_csv(OUT / "cohort_summary.csv", index=False)
    a2.to_csv(OUT / "band_taxonomy_raw_vs_multiscale.csv", index=False)
    a3.to_csv(OUT / "snr_gating.csv", index=False)
    cert.to_csv(OUT / "recompute_certificate.csv", index=False)

    # ----- console digest -------------------------------------------------- #
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 40)
    print("\n================ AXIS 1 — MAGNITUDE (multiscale - raw) ===========")
    print(a1[["cohort", "band", "n", "median_multiscale", "median_raw",
              "median_delta", "n_multiscale_gt_raw",
              "p_delta_gt0_one_sided", "p_delta_two_sided"]]
          .round(4).to_string(index=False))
    print("\n================ AXIS 2 — STORY (band taxonomy) ==================")
    print(a2[["band", "taxonomy_multiscale", "ms_median", "ms_n_clears_null",
              "raw_median", "raw_n_clears_null"]].round(4).to_string(index=False))
    print("\n================ AXIS 3 — SNR GATING (high vs low) ===============")
    print(a3[["band", "multiscale_hi_minus_lo", "multiscale_mw_p_hi_gt_lo",
              "raw_hi_minus_lo", "raw_mw_p_hi_gt_lo"]]
          .round(4).to_string(index=False))
    print(f"\n[audit_146] wrote outputs to {OUT}")


if __name__ == "__main__":
    main()
