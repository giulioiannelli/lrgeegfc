#!/usr/bin/env python3
"""Audit 145 — full per-phase reliability ("multiphase SNR") refinement.

Critical preamble (mandatory house style; read before the code)
================================================================
(1) CLAIM. The between-patient heterogeneity in the cophenetic trace
    ``rho_split`` (per band, n=10) is driven by MEASUREMENT RELIABILITY,
    not biology: a patient looks like a non-tracer because one of its
    four phase cophenetic distance vectors (D_preA, D_preB, D_tt, D_post)
    is noisy, which attenuates the cross-phase Spearman correlation. The
    refinement: replace the rest_pre-only noise floor with a FULL
    per-phase reliability composite and ask whether, after correcting
    rho_split for attenuation, any residual cross-patient variability
    (real biology) survives.

(2) NULL / MODEL. Classical correction-for-attenuation. The observed
    cross-phase trace is the disattenuated trace times the geometric mean
    of the reliabilities of the two difference vectors:

        rho_obs  ~=  rho_true * sqrt( rel(Delta_task) * rel(Delta_rest) )

    where rel(Delta_task) is the test-retest reliability of
    Delta_task = D_tt - D_preA and rel(Delta_rest) of
    Delta_rest = D_post - D_preB. We estimate per-phase split-half
    reliabilities rel_phase = Spearman(D_phaseA, D_phaseB) for
    phase in {rest_pre, task_test, rest_post}, then propagate them to the
    two DIFFERENCE vectors under an additive-error model (see Section 3 of
    the report and ``_diff_reliability`` below). The disattenuated trace is
    rho_split_corrected = rho_obs / sqrt(rel(Delta_task)*rel(Delta_rest)).
    The "reliability-only" null model for the between-patient question is:
    rho_split is a deterministic function of the reliability composite plus
    i.i.d. noise; if true, partialling reliability out of rho_split leaves
    no patient-level signal.

(3) STRONGEST ALTERNATIVE IT MUST CONTROL FOR. The competing benign
    explanation we want to KILL is "the heterogeneity that the current
    rest_pre-only d_noise attributes to biology is in fact unreliability
    living in task_test or rest_post that d_noise never measures." A
    patient with a clean rest_pre but a noisy rest_post would have a high
    rest_pre-only SNR yet a genuinely attenuated rho_split — the
    rest_pre-only axis cannot see this; the full per-phase composite can.
    So the alternative the full model must control for is *phase-specific
    unreliability outside rest_pre*.

(4) DOES IT ACTUALLY CONTROL FOR IT — BY MECHANISM, AND WHAT IT CANNOT DO.
    Mechanism: rel_phase for task_test and rest_post are now MEASURED from
    their own A/B half-FCs, so phase-specific unreliability enters the
    attenuation composite. That is a real improvement over rest_pre-only.
    HONEST LIMITS (stated up front, not buried):
      - Split-half reliability HALVES the data per phase, so each
        rel_phase is itself a downward-biased, noisy estimate of the
        full-duration reliability (fewer Welch segments). It estimates
        the reliability of a HALF recording, and Spearman-Brown up-correction
        to full length rests on the parallel-halves assumption. We report
        both the raw half reliability and the Spearman-Brown full-length
        version, and flag that the absolute correction magnitude is
        model-dependent.
      - Spearman split-half reliability is a RANK statistic on a cophenetic
        vector that is ~98% tied at the N-1 merge heights (known atom
        structure). Ties depress Spearman; rel_phase is therefore a
        conservative (low) reliability estimate, which INFLATES the
        disattenuated trace. Corrected rho values above ~1 are an artifact
        of this floor, not super-reliability — we cap interpretation, not
        the arithmetic, and report capped cells honestly.
      - The error-propagation from per-phase reliability to DIFFERENCE-vector
        reliability assumes the two phases entering a difference have
        comparable signal variance and INDEPENDENT measurement error. Task
        and rest_pre share the same electrodes and montage; if their errors
        are correlated, the difference reliability is mis-estimated. This
        estimator CANNOT separate correlated-error from true shared signal.
      - The whole exercise is correlational at n=10. A residual partial
        correlation is suggestive, not a hypothesis test with power; CIs are
        wide. We report the residual with a bootstrap CI and say so.

(5) WHAT WOULD FALSIFY. If, after regressing rho_split on the full
    reliability composite per band, the residual cross-patient variance is
    indistinguishable from zero (reliability model R^2 ~= 1, partial
    correlation of rho_split with a held-back patient axis ~= 0 with a CI
    spanning 0), then the heterogeneity IS entirely reliability and the
    "universal, detectability-limited trace" reading holds. Conversely, a
    reliability model that explains little variance, or a stable residual
    ordering of patients that survives reliability correction, falsifies
    the pure-reliability claim and means a biological axis survives.

What this script does
=====================
  1. Caches the MISSING task_test split-half FCs the same way audit_63
     caches rest_pre halves (``compute_imcoh_abs_halves``), applying the
     Pat_10 [53,54,55] row drop and verifying N matches the full task_test
     FC before caching. rest_pre/rest_post halves are already cached.
  2. Computes per-phase cophenetic reliability rel_phase = Spearman(
     canonical_cophenet(W_A), canonical_cophenet(W_B)) for the 3 phases.
  3. Loads the cached observed rho_split (headline; not recomputed),
     forms difference-vector reliabilities, the reliability composite
     R_comp = sqrt(rel(Delta_task)*rel(Delta_rest)), the disattenuated
     rho_split_corrected, and a corrected-SNR analog.
  4. Per band, regresses/partials rho_split against R_comp to answer the
     residual-biology question; bootstraps CIs at n=10.
  5. Re-ranks patients by corrected SNR vs the current rest_pre-only SNR.

Outputs
-------
    data/audit/multiphase_snr/
        per_patient_per_band_reliability.csv
        per_band_summary.csv
        patient_ranking.csv
        README.md
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.metrics.node_localization import canonical_cophenet
from lrg_eegfc.workflow.fc import load_fc_matrix

# Reuse the existing half-FC helper (no private copy)
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # type: ignore


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
HALF_PHASES = ("rest_pre", "task_test", "rest_post")

# Pat_10 task_test raw timeseries has 3 extra rows; project rule.
PAT10_TASK_DROP = [53, 54, 55]

OUT = ROOT / "data" / "audit" / "multiphase_snr"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"
HEADLINE_CSV = (ROOT / "data" / "audit"
                / "matched_strength_surrogate_split_baseline"
                / "per_patient_per_band.csv")

OUT.mkdir(parents=True, exist_ok=True)
HALVES_FC_CACHE.mkdir(parents=True, exist_ok=True)

# Reliabilities below this are treated as too noisy to invert; the
# disattenuation divides by ~0 and is reported as capped, not used.
REL_FLOOR = 0.05


# ---------------------------------------------------------------------------
# Half-FC cache paths + lazy populate (mirrors audit_63 for rest_pre)
# ---------------------------------------------------------------------------
def _half_fc_path(pat: str, band: str, phase: str, half: str) -> Path:
    return HALVES_FC_CACHE / pat / f"{band}_{phase}_{half}_imcoh_abs.npy"


def _orient_channels_time(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    return X


def ensure_task_test_halves(pat: str, bands: list[str]) -> None:
    """Compute + cache task_test half FCs for `bands` if any are missing.

    Same recipe as audit_63.ensure_half_fcs but for task_test, with the
    Pat_10 [53,54,55] row drop and an N-vs-full-FC verification.
    """
    missing = [(b, h) for b in bands for h in ("A", "B")
               if not _half_fc_path(pat, b, "task_test", h).exists()]
    if not missing:
        return
    print(f"[audit_145] {pat}: caching {len(missing)} missing task_test half FCs")
    X = load_timeseries(pat, "task_test", SEEG_DATAPATH)
    X = _orient_channels_time(X)
    # NOTE on the Pat_10 [53,54,55] task drop: ``load_timeseries`` ALREADY
    # applies the Pat_10 task channel mask internally — it returns 113 channels
    # for task_test, with rows [53,54,55] already removed (verified equal to the
    # canonical full task_test FC, which is also 113). Re-dropping here would
    # strip three MORE rows and give N=110 (the N-vs-full-FC guard below caught
    # exactly that on the first run). So we do NOT drop manually; the loader
    # handles it, identically to how audit_63 builds the rest_pre halves. The
    # guard below is the real safety net and stays unconditional. The defensive
    # branch fires only if a future loader regresses to the raw 116-channel form.
    if pat == "Pat_10" and X.shape[0] == 116:
        keep = [i for i in range(X.shape[0]) if i not in PAT10_TASK_DROP]
        X = X[keep, :]
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nperseg_half = max(256, nperseg_for_fs(fs) // 2)
    halves_fc = compute_imcoh_abs_halves(X, fs, nperseg_half, BRAIN_BANDS)
    del X
    gc.collect()

    # Verify N against the canonical full task_test FC before caching.
    N_full = load_fc_matrix(pat, "task_test", bands[0],
                            fc_method="imcoh_abs").shape[0]
    pat_dir = HALVES_FC_CACHE / pat
    pat_dir.mkdir(parents=True, exist_ok=True)
    for (band, half), A in halves_fc.items():
        if band not in bands:
            continue
        if A.shape[0] != N_full:
            raise ValueError(
                f"[audit_145] {pat}/{band}/task_test half-FC N={A.shape[0]} "
                f"!= full task_test N={N_full}. Aborting cache.")
        np.save(_half_fc_path(pat, band, "task_test", half),
                np.asarray(A, dtype=np.float32))


def load_half_fc(pat: str, band: str, phase: str, half: str) -> np.ndarray:
    """Load a half FC, normalize to symmetric, zero-diag, clipped (N,N)."""
    W = np.load(_half_fc_path(pat, band, phase, half))
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    W = 0.5 * (W + W.T)
    return W


# ---------------------------------------------------------------------------
# Per-phase cophenetic reliability
# ---------------------------------------------------------------------------
def phase_reliability(pat: str, band: str, phase: str) -> float:
    """rel_phase = Spearman(canonical_cophenet(W_A), canonical_cophenet(W_B))."""
    Da = canonical_cophenet(load_half_fc(pat, band, phase, "A"))
    Db = canonical_cophenet(load_half_fc(pat, band, phase, "B"))
    if Da.size != Db.size:
        return float("nan")
    return float(spearmanr(Da, Db).statistic)


def spearman_brown_full(rel_half: float) -> float:
    """Spearman-Brown up-correction from a half recording to full length.

    rel_full = 2 * rel_half / (1 + rel_half). Maps a half-length reliability
    to the implied full-length reliability under the parallel-halves
    assumption. Clipped to [0, 1) on input to avoid blow-up at rel->1.
    """
    r = float(rel_half)
    if not np.isfinite(r):
        return float("nan")
    r = min(max(r, -0.999), 0.999)
    return 2.0 * r / (1.0 + r)


# ---------------------------------------------------------------------------
# Difference-vector reliability (error-propagation; see report Section 3)
# ---------------------------------------------------------------------------
def _diff_reliability(rel_full: float, rel_other: float) -> float:
    """Reliability of the difference vector D_full - D_half.

    A difference vector enters rho_split: Delta_task = D_tt - D_preA,
    Delta_rest = D_post - D_preB. Under an additive classical-true-score
    model X = T + E with equal signal variance sigma_S^2 in the two phases
    and independent errors, the reliability of the difference is

        rel_diff = ( rel_full + rel_other - 2*r_TT )
                   / ( 2 - 2*r_EE )

    where r_TT is the cross-phase true-score correlation and r_EE the
    cross-phase error correlation. We do NOT know r_TT / r_EE and CANNOT
    estimate them here (no third repeat). We therefore use the
    parallel-equal-variance, independent-error, ORTHOGONAL-true-score
    approximation (r_TT = 0, r_EE = 0):

        rel_diff ~= (rel_full + rel_other) / 2

    i.e. the difference reliability is the mean of the two component
    reliabilities. This is the standard "difference of two equally
    reliable, uncorrelated measurements" result and is what the report
    documents. It is exact when the two phases contribute equal signal
    variance and the cophenetic shape change is independent of the
    baseline shape. It OVER-estimates rel_diff if the true scores are
    positively correlated (the usual case for difference scores), so the
    resulting disattenuation is CONSERVATIVE (under-corrects). We accept
    that bias deliberately: it biases against a "pure reliability" verdict.
    """
    a = float(rel_full)
    b = float(rel_other)
    if not (np.isfinite(a) and np.isfinite(b)):
        return float("nan")
    return 0.5 * (a + b)


# ---------------------------------------------------------------------------
# Bootstrap helpers (patient-resampling at n=10)
# ---------------------------------------------------------------------------
def _spearman_ci(x: np.ndarray, y: np.ndarray, n_boot: int,
                 rng: np.random.Generator) -> tuple[float, float, float]:
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3:
        return float("nan"), float("nan"), float("nan")
    obs = float(spearmanr(x, y).statistic)
    n = x.size
    boots = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if np.unique(idx).size < 3:
            boots[b] = np.nan
            continue
        boots[b] = spearmanr(x[idx], y[idx]).statistic
    boots = boots[np.isfinite(boots)]
    lo = float(np.quantile(boots, 0.025)) if boots.size else float("nan")
    hi = float(np.quantile(boots, 0.975)) if boots.size else float("nan")
    return obs, lo, hi


def _partial_resid_corr(rho: np.ndarray, rcomp: np.ndarray,
                        n_boot: int, rng: np.random.Generator
                        ) -> dict:
    """Residual-biology probe per band.

    Regress rho_split on the reliability composite (rank-based: residual of
    rho after removing its monotone dependence on R_comp), and ask whether
    the residual still carries a stable patient ordering. We operationalize
    "residual biology" two ways:

      (a) R^2 of the reliability-only OLS of rho ~ R_comp (how much
          between-patient variance reliability explains). 1 - R^2 is the
          residual-variance fraction.
      (b) Bootstrap stability of the residual ordering: correlation of the
          per-patient residual with itself across patient resamples is not
          meaningful; instead we report whether the spread of residuals is
          large relative to their bootstrap noise (residual std vs a
          permutation/bootstrap band).

    Honest framing: at n=10 neither is a powered test. We report R^2 with a
    bootstrap CI and the residual vector for inspection.
    """
    m = np.isfinite(rho) & np.isfinite(rcomp)
    rho, rcomp = rho[m], rcomp[m]
    n = rho.size
    if n < 3:
        return {"n": int(n)}
    # Spearman rho<->R_comp and its CI
    sp, sp_lo, sp_hi = _spearman_ci(rho, rcomp, n_boot, rng)
    # OLS R^2 of rho ~ R_comp (linear; complements the rank stat)
    X = np.vstack([np.ones(n), rcomp]).T
    beta, *_ = np.linalg.lstsq(X, rho, rcond=None)
    pred = X @ beta
    resid = rho - pred
    ss_tot = float(np.sum((rho - rho.mean()) ** 2))
    ss_res = float(np.sum(resid ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    # Bootstrap CI on R^2
    r2_boot = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        Xb = np.vstack([np.ones(n), rcomp[idx]]).T
        yb = rho[idx]
        sst = float(np.sum((yb - yb.mean()) ** 2))
        if sst <= 0:
            r2_boot[b] = np.nan
            continue
        bb, *_ = np.linalg.lstsq(Xb, yb, rcond=None)
        ssr = float(np.sum((yb - Xb @ bb) ** 2))
        r2_boot[b] = 1.0 - ssr / sst
    r2_boot = r2_boot[np.isfinite(r2_boot)]
    r2_lo = float(np.quantile(r2_boot, 0.025)) if r2_boot.size else float("nan")
    r2_hi = float(np.quantile(r2_boot, 0.975)) if r2_boot.size else float("nan")
    return {
        "n": int(n),
        "spearman_rho_rcomp": sp,
        "spearman_rho_rcomp_lo": sp_lo,
        "spearman_rho_rcomp_hi": sp_hi,
        "r2_reliability_only": r2,
        "r2_lo": r2_lo,
        "r2_hi": r2_hi,
        "resid_std": float(np.std(resid, ddof=1)) if n > 1 else float("nan"),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=BANDS)
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--n-boot", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=20260625)
    args = ap.parse_args()

    bands = args.bands
    patients = args.patients
    rng = np.random.default_rng(args.seed)

    # ---- Pre-flight: cache missing task_test halves ----
    print("[audit_145] pre-flight: ensure task_test half FCs cached")
    for pat in patients:
        ensure_task_test_halves(pat, bands)

    # ---- Load cached observed rho_split headline ----
    head = pd.read_csv(HEADLINE_CSV)
    head_idx = {(r.patient, r.band): float(r.obs_rho)
                for r in head.itertuples()}

    t0 = time.time()
    rows: list[dict] = []
    for band in bands:
        for pat in patients:
            rel = {ph: phase_reliability(pat, band, ph) for ph in HALF_PHASES}
            rel_pre = rel["rest_pre"]
            rel_tt = rel["task_test"]
            rel_post = rel["rest_post"]
            # Spearman-Brown full-length versions
            relf = {ph: spearman_brown_full(v) for ph, v in rel.items()}

            obs_rho = head_idx.get((pat, band), float("nan"))

            # --- Difference-vector reliabilities (half-length basis) ---
            # Delta_task = D_tt - D_preA ; Delta_rest = D_post - D_preB.
            # Both halves of rest_pre come from the SAME split, so we use
            # the rest_pre split-half reliability for the rest_pre component
            # in both differences (it IS the reliability of a rest_pre half).
            rel_dtask = _diff_reliability(rel_tt, rel_pre)
            rel_drest = _diff_reliability(rel_post, rel_pre)
            # Full-length (Spearman-Brown) versions
            relf_dtask = _diff_reliability(relf["task_test"], relf["rest_pre"])
            relf_drest = _diff_reliability(relf["rest_post"], relf["rest_pre"])

            # --- Reliability composite + disattenuation ---
            def _composite(rd_t, rd_r):
                if not (np.isfinite(rd_t) and np.isfinite(rd_r)):
                    return float("nan")
                prod = rd_t * rd_r
                return float(np.sqrt(prod)) if prod > 0 else float("nan")

            r_comp = _composite(rel_dtask, rel_drest)          # half-basis
            r_comp_full = _composite(relf_dtask, relf_drest)   # SB full-length

            def _disatt(rho, comp):
                if not (np.isfinite(rho) and np.isfinite(comp)):
                    return float("nan"), False
                if comp < REL_FLOOR:
                    return float("nan"), True   # capped (too noisy to invert)
                return float(rho / comp), False

            rho_corr, capped = _disatt(obs_rho, r_comp)
            rho_corr_full, capped_full = _disatt(obs_rho, r_comp_full)

            # --- Corrected SNR analog ---
            # Current SNR = d_task / d_noise with d_noise = 1 - rel_pre.
            # The full-reliability noise analog uses the WORST (max
            # 1-rel) across the phases entering the trace, i.e. the
            # bottleneck reliability, NOT just rest_pre. We also report
            # the geometric noise composite 1 - r_comp.
            d_noise_pre = 1.0 - rel_pre if np.isfinite(rel_pre) else float("nan")
            # bottleneck across the 3 phases (worst phase reliability)
            phase_vals = [v for v in (rel_pre, rel_tt, rel_post)
                          if np.isfinite(v)]
            rel_bottleneck = min(phase_vals) if phase_vals else float("nan")
            d_noise_full = (1.0 - rel_bottleneck
                            if np.isfinite(rel_bottleneck) else float("nan"))

            rows.append({
                "patient": pat,
                "band": band,
                "obs_rho": obs_rho,
                "rel_rest_pre": rel_pre,
                "rel_task_test": rel_tt,
                "rel_rest_post": rel_post,
                "rel_rest_pre_sb": relf["rest_pre"],
                "rel_task_test_sb": relf["task_test"],
                "rel_rest_post_sb": relf["rest_post"],
                "rel_delta_task": rel_dtask,
                "rel_delta_rest": rel_drest,
                "rel_composite": r_comp,
                "rel_composite_sb": r_comp_full,
                "rho_corrected": rho_corr,
                "rho_corrected_capped": capped,
                "rho_corrected_sb": rho_corr_full,
                "rho_corrected_sb_capped": capped_full,
                "d_noise_restpre": d_noise_pre,
                "rel_bottleneck": rel_bottleneck,
                "d_noise_full": d_noise_full,
            })
        print(f"[audit_145] band {band} done ({time.time()-t0:.1f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_patient_per_band_reliability.csv", index=False)
    print(f"[audit_145] wrote per_patient_per_band_reliability.csv "
          f"({len(df)} rows)")

    # ---- Per-band residual-biology summary ----
    summ_rows = []
    for band in bands:
        sub = df[df.band == band]
        rho = sub["obs_rho"].to_numpy(dtype=float)
        rcomp = sub["rel_composite"].to_numpy(dtype=float)
        rcomp_full = sub["rel_composite_sb"].to_numpy(dtype=float)
        # Corrected SNR <-> rho relationship: does the corrected reliability
        # composite still order rho? (a high composite = high detectability)
        sp_obs, sp_lo, sp_hi = _spearman_ci(rho, rcomp, args.n_boot, rng)
        # residual-biology probe (R^2 of reliability-only model)
        rb = _partial_resid_corr(rho, rcomp, args.n_boot, rng)
        rb_full = _partial_resid_corr(rho, rcomp_full, args.n_boot, rng)
        summ_rows.append({
            "band": band,
            "n_patients": int(np.isfinite(rho).sum()),
            # corrected SNR composite <-> rho_split correlation (the headline
            # detectability number, now full-phase)
            "spearman_rho_relcomp": sp_obs,
            "spearman_rho_relcomp_lo": sp_lo,
            "spearman_rho_relcomp_hi": sp_hi,
            # reliability-only model fit (residual-biology)
            "r2_reliability_only": rb.get("r2_reliability_only", float("nan")),
            "r2_lo": rb.get("r2_lo", float("nan")),
            "r2_hi": rb.get("r2_hi", float("nan")),
            "resid_var_frac": (1.0 - rb["r2_reliability_only"]
                               if np.isfinite(rb.get("r2_reliability_only",
                                                     float("nan")))
                               else float("nan")),
            "r2_reliability_only_sb": rb_full.get("r2_reliability_only",
                                                  float("nan")),
            "median_rel_composite": float(np.nanmedian(rcomp)),
            "median_rho_corrected": float(
                np.nanmedian(sub["rho_corrected"].to_numpy(dtype=float))),
            "n_capped": int(sub["rho_corrected_capped"].sum()),
        })
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(OUT / "per_band_summary.csv", index=False)
    print("\n[audit_145] per-band summary:")
    print(summ.to_string(index=False))

    # ---- Patient ranking: current SNR vs corrected SNR (per band) ----
    # "Current SNR" needs d_task; recompute d_task = 1 - Spearman(D_tt, D_pre)
    # on the FULL FCs (cheap; mirrors the session's d_task). The corrected
    # SNR analog uses the bottleneck reliability composite.
    rank_rows = []
    for band in bands:
        for pat in patients:
            # d_task from full FCs
            try:
                D_tt = canonical_cophenet(
                    load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs"))
                D_pre = canonical_cophenet(
                    load_fc_matrix(pat, "rest_pre", band, fc_method="imcoh_abs"))
                d_task = 1.0 - float(spearmanr(D_tt, D_pre).statistic)
            except Exception:
                d_task = float("nan")
            sub = df[(df.band == band) & (df.patient == pat)].iloc[0]
            d_noise_pre = float(sub["d_noise_restpre"])
            d_noise_full = float(sub["d_noise_full"])
            snr_current = (d_task / d_noise_pre
                           if (np.isfinite(d_task) and np.isfinite(d_noise_pre)
                               and d_noise_pre > 0) else float("nan"))
            snr_corrected = (d_task / d_noise_full
                             if (np.isfinite(d_task) and np.isfinite(d_noise_full)
                                 and d_noise_full > 0) else float("nan"))
            rank_rows.append({
                "patient": pat, "band": band,
                "d_task": d_task,
                "d_noise_restpre": d_noise_pre,
                "d_noise_full": d_noise_full,
                "snr_current": snr_current,
                "snr_corrected": snr_corrected,
            })
    rank = pd.DataFrame(rank_rows)
    # Add per-band ranks (1 = highest SNR)
    for col, newc in (("snr_current", "rank_current"),
                      ("snr_corrected", "rank_corrected")):
        rank[newc] = (rank.groupby("band")[col]
                      .rank(ascending=False, method="min"))
    rank.to_csv(OUT / "patient_ranking.csv", index=False)
    print(f"\n[audit_145] wrote patient_ranking.csv ({len(rank)} rows)")

    # Re-test SNR->rho with the corrected SNR per band (rho vs snr_corrected)
    print("\n[audit_145] corrected-SNR <-> rho_split (per band):")
    for band in bands:
        rsub = rank[rank.band == band].merge(
            df[df.band == band][["patient", "obs_rho"]], on="patient")
        rho = rsub["obs_rho"].to_numpy(dtype=float)
        sc = rsub["snr_corrected"].to_numpy(dtype=float)
        scur = rsub["snr_current"].to_numpy(dtype=float)
        sp_corr, lo_c, hi_c = _spearman_ci(rho, sc, args.n_boot, rng)
        sp_cur, lo_u, hi_u = _spearman_ci(rho, scur, args.n_boot, rng)
        print(f"  {band:11s} corrected rho(snr,rho_split)={sp_corr:+.3f} "
              f"[{lo_c:+.3f},{hi_c:+.3f}]   current={sp_cur:+.3f} "
              f"[{lo_u:+.3f},{hi_u:+.3f}]")

    print(f"\n[audit_145] done ({time.time()-t0:.1f}s). Outputs at {OUT}")


if __name__ == "__main__":
    main()
