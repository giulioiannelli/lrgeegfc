#!/usr/bin/env python3
"""Audit 143 — distance-SWAP head-to-head on the team's OWN ρ_split statistic.

==========================================================================
MANDATORY 5-POINT CRITICAL PREAMBLE
==========================================================================

(1) CLAIM UNDER TEST.
    The paper's load-bearing methodological claim is that the FULL multiscale
    Laplacian read-out — the cophenetic per-pair distance fed into the team's
    split-baseline ρ_split — detects a cross-phase persistence trace that a
    GLOBAL-SPECTRAL view of the SAME adjacency (the leading-k Laplacian
    eigenvector embedding spectral clustering runs k-means on, BEFORE the cut)
    and a RAW-FC edge baseline CANNOT. Sharpest case = ALPHA (locked ρ^coph
    split-baseline cohort p≈0.00195; α is NOT visible to the global-eigenmode
    Grassmann probe). β is the positive control (trace on both probes); raw FC
    should be weakest everywhere (β is the weakest band at raw FC per the team's
    account). θ null everywhere.

    v1 (audit_141) was inconclusive because it compared methods on a distance
    TRIANGLE + spectral-clustering PARTITION similarity (ARI/NMI) — the WRONG
    common footing: the team's actual trace statistic is
        ρ_split = Spearman(Δ_task, Δ_rest)
    a correlation of per-pair cross-phase distance SHIFTS, materially more
    sensitive than the triangle. v1's triangle ρ^coph did not reproduce the
    locked numbers, so v1 pitted the baselines against a WEAKENED ρ^coph.

    v2 FIX — HOLD THE STATISTIC FIXED, SWAP ONLY THE DISTANCE. ρ_split is a
    recipe with three steps:
        step 1: per-pair distance D[i,j]                      ← SWAP THIS
        step 2: shifts Δ_task = D^tt − D^pre_A,
                       Δ_rest = D^post − D^pre_B               ← IDENTICAL
        step 3: ρ_split = Spearman(Δ_task, Δ_rest) over triu  ← IDENTICAL
    Only step 1 is "cophenetic". The fair head-to-head keeps steps 2–3 byte-
    identical across arms and swaps the distance:
        ARM 1 cophenetic   D_coph = cophenet(UPGMA(T(τ=1/λmax)))  [team's measure]
        ARM 2 spectral-k   D_spec^(k)[i,j] = ||Y_i − Y_j||_2, Y = leading-k
                           non-trivial normalized-Laplacian eigenvectors
                           (cols v_2..v_{k+1}) — the PRE-CUT spectral embedding;
                           sweep k ∈ {2,3,5,8,13}
        ARM 3 raw FC       D_raw[i,j] = 1 − |ImCoh|[i,j]  (no diffusion, no tree)
        ARM 4 diffusion-τ  D_diff[i,j] = (1/ρ)[i,j] at the SINGLE scale τ=1/λmax
                           (the heat-kernel communication distance BEFORE UPGMA)
                           — separates "diffusion at one scale" from "multiscale
                           hierarchical integration".

(2) THE NULL.
    The audit_63 split-baseline matched-strength family: R=200 strength-
    preserving 4-cycle ±δ rewired Laplacians per (patient, band, phase ∈
    {pre_A, pre_B, task_test, rest_post}), INDEPENDENT per phase, seed=20260511
    — REUSED from the cached eigendecompositions
    (data/cache/matched_strength_surrogate_lrg/.../*_seed20260511_*). Every arm
    recomputes its OWN ρ_split per surrogate from the SAME cached ensemble:
        - cophenetic   from cophenetic_condensed_from_eigs(eigs)
        - spectral-k   from the normalized-Laplacian eigvecs of the surrogate
                       adjacency reconstructed via adjacency_from_laplacian_eigs
        - raw FC       from 1 − adjacency_from_laplacian_eigs(eigs)
        - diffusion-τ  from the heat-kernel communication distance of eigs
    Per-patient upper-tail p = mean(ρ_surr >= ρ_obs); cohort verdict =
    one-sided Wilcoxon (wilcoxon_z) of (ρ_obs − surr_median) > 0; BH-FDR
    (bh_fdr) WITHIN each arm across the 6-band family.

(3) STRONGEST ALTERNATIVE THE NULL MUST CONTROL FOR.
    Per-node STRENGTH heterogeneity drift auto-correlated across phases (the
    KC-β collapse of 2026-05-11): a distance can look like it carries a trace
    purely because node strengths persist. The matched-strength null reproduces
    each phase's exact strength sequence while scrambling edge identity, so a
    surviving ρ_split cannot be attributed to strength alone. The team's
    mandatory minimum null for any cohort FC claim, applied IDENTICALLY to all
    four arms ⇒ APPLES-TO-APPLES.

(4) DOES THE NULL ACTUALLY CONTROL FOR IT — BY MECHANISM.
    Yes on the strength axis: 4-cycle ±δ preserves the strength sequence exactly
    (verify_strengths tol 1e-4 at cache build). WHAT IT CANNOT DO: it is
    INDEPENDENT per phase, so "clears the null" means "ρ_split exceeds an
    independent-strength-randomized baseline", NOT "specific edge-identity
    cross-phase structure carries the trace" — IDENTICAL caveat to audit_63, the
    same null, shared across arms, so the COMPARISON is fair even though no arm's
    null isolates edge identity. Arm-specific caveat: the spectral-k embedding
    distance ||Y_i − Y_j|| is sign/rotation-AMBIGUOUS in the eigenvectors, but
    the ambiguity cancels WITHIN a phase's distance matrix (Euclidean distances
    between embedded rows are invariant to a global orthogonal transform of the
    columns and to per-column sign), so the per-pair distance — and therefore
    the SHIFT and ρ_split — is well-defined without any cross-phase alignment.

(5) WHAT WOULD FALSIFY THE CLAIM / REMAINING LIMITATIONS.
    GATE (mandatory): the cophenetic arm MUST reproduce the locked audit_63
    split-baseline cohort p — α≈0.00195, β≈0.00488 — to reasonable tolerance
    (it reuses the SAME cached surrogate ensemble + the SAME ρ_split + the SAME
    half-FC inputs, so this should be ~exact). If it does not, STOP and debug;
    the whole point of v2 is to use the REAL statistic.
    The thesis is FALSIFIED / WEAKENED if the leading-k spectral embedding
    recovers the α trace (clears matched-strength) at ANY swept k — that would
    mean the full multiscale cophenetic distance is not adding detection power
    over the pre-cut global-spectral view in the sharpest case. We give the
    embedding its BEST k (report the smallest-p k) and the normalized-Laplacian
    (Ng-Jordan-Weiss) embedding it actually uses. Remaining limitations: the
    null is independent-per-phase (shared, fair, isolates strength not edge
    identity); the spectral embedding is row-normalized per NJW so near-zero-
    degree nodes get a stabilized unit row; raw FC and diffusion-τ are unhierar-
    chical single-view baselines by construction (that IS the contrast).

==========================================================================

Outputs
-------
    data/audit/spectral_distance_swap/
        cohort_summary.csv          # band × arm → ρ_split trace verdict (+ best k)
        per_patient_per_band.csv    # every (patient, band, arm, k) cell
        README.md

Library reuse (no private forks)
--------------------------------
- FC loading (full phases)   : audit_141 load_phase_fc idiom (workflow.fc.load_fc_matrix)
- half-FC (pre_A / pre_B)    : cached at data/cache/imcoh_halves_fc/ (audit_63 layout);
                               lazily (re)built via _fc_split_half.compute_imcoh_abs_halves
- matched-strength surrogate : utils.surrogate.matched_strength
                               (load_or_compute_eigs_at_path — reuses the cached
                                seed=20260511 split-baseline ensemble;
                                adjacency_from_laplacian_eigs;
                                cophenetic_condensed_from_eigs)
- cohort stats               : utils.metrics.hypothesis.wilcoxon_z, bh_fdr

Surrogate cache
---------------
This audit REUSES the canonical audit_63 split-baseline ensemble at
``data/cache/matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz``
for phase ∈ {rest_pre_A, rest_pre_B, task_test, rest_post}. The cophenetic arm
therefore reproduces the locked audit_63 number from the IDENTICAL surrogate
draws (the correctness gate). The seed/swap_factor default to 20260511 / 20 so a
cache hit is guaranteed; pass --seed to regenerate a fresh ensemble.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BAND_TEX_DICT, FS_OVERRIDES, nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, wilcoxon_z
from lrg_eegfc.utils.surrogate.matched_strength import (
    adjacency_from_laplacian_eigs,
    cophenetic_condensed_from_eigs,
    load_or_compute_eigs_at_path,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

# Reuse the existing half-FC helper (no copy) — same import audit_63 uses.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # type: ignore

SURROGATE_LRG_CACHE = CACHE_ROOT / "matched_strength_surrogate_lrg"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
# alpha = PRIMARY; beta = positive control (both probes); theta = null;
# delta / low_gamma / high_gamma = remaining 6-band family for BH-FDR.
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
# Split-baseline phases (audit_63 statistic of record).
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
K_GRID = (2, 3, 5, 8, 13)        # leading-k spectral embedding sweep
N_SURROGATES = 200
SWAP_FACTOR = 20
# Canonical audit_63 split-baseline seed: REUSE the cached ensemble so the
# cophenetic arm reproduces the locked number from the identical draws.
SEED = 20260511

# Locked audit_63 split-baseline ρ_split cohort p-values (the gate target).
# Source: data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
AUDIT63_LOCKED_P = {
    "delta": 0.2783203125, "theta": 0.7216796875, "alpha": 0.001953125,
    "beta": 0.0048828125, "low_gamma": 0.1162109375, "high_gamma": 0.24609375,
}
GATE_BANDS = ("alpha", "beta")   # bands checked against the locked p
GATE_TOL = 0.01                  # |p_coph − p_locked| tolerance for PASS

OUT = ROOT / "data" / "audit" / "spectral_distance_swap"
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Phase FC loading (full + half), symmetrized, clipped — audit_63 / audit_141 idiom
# ---------------------------------------------------------------------------
def _half_fc_path(pat: str, band: str, half: str) -> Path:
    return HALVES_FC_CACHE / pat / f"{band}_rest_pre_{half}_imcoh_abs.npy"


def ensure_half_fcs(pat: str, bands: list[str]) -> None:
    """Compute and cache rsPre half FCs for `bands` if any are missing
    (audit_63 pre-flight, verbatim helper reuse)."""
    missing = [(b, h) for b in bands for h in ("A", "B")
               if not _half_fc_path(pat, b, h).exists()]
    if not missing:
        return
    print(f"[audit_143] {pat}: caching {len(missing)} missing half FCs")
    X = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH),
                   dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nperseg_half = max(256, nperseg_for_fs(fs) // 2)
    halves_fc = compute_imcoh_abs_halves(X, fs, nperseg_half, BRAIN_BANDS)
    pat_dir = HALVES_FC_CACHE / pat
    pat_dir.mkdir(parents=True, exist_ok=True)
    for (band, half), A in halves_fc.items():
        if band not in bands:
            continue
        np.save(_half_fc_path(pat, band, half),
                np.asarray(A, dtype=np.float32))


def load_phase_fc(pat: str, phase: str, band: str) -> np.ndarray:
    """Phase FC (full or half) → (N,N) float64, zero-diag, clipped, symmetrized."""
    if phase in ("rest_pre_A", "rest_pre_B"):
        W = np.load(_half_fc_path(pat, band, phase[-1]))
    else:
        W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
        if W is None:
            raise FileNotFoundError(f"no cached FC for {pat}/{phase}/{band}")
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


# ---------------------------------------------------------------------------
# STEP 1 — the four per-pair condensed distances (the ONLY thing that varies).
# Each returns a condensed (upper-triangular) vector of length N(N−1)/2.
# ---------------------------------------------------------------------------
def _condensed_upper(D: np.ndarray) -> np.ndarray:
    iu = np.triu_indices(D.shape[0], k=1)
    return D[iu]


def dist_spectral_embedding(eigvals_L: np.ndarray, eigvecs_L: np.ndarray,
                            W: np.ndarray, k: int) -> np.ndarray:
    """ARM 2 — Euclidean distance in the leading-k NORMALIZED-Laplacian
    eigenvector embedding (the pre-cut spectral-clustering representation).

    L_sym = I − D^-1/2 W D^-1/2; take the k non-trivial smallest-eigenvalue
    eigenvectors (skip the trivial mode at index 0), row-normalize (Ng-Jordan-
    Weiss), then D_spec^(k)[i,j] = ||Y_i − Y_j||_2. Recomputed from W (the
    unnormalized-Laplacian eigs cached for the surrogate are L = D−W, NOT L_sym,
    so we form L_sym explicitly — identical to how observed and surrogate are
    treated). ``eigvals_L``/``eigvecs_L`` are accepted for signature symmetry
    but the normalized embedding is built from ``W`` directly.
    """
    deg = W.sum(axis=1)
    with np.errstate(divide="ignore"):
        d_inv_sqrt = np.where(deg > 0, 1.0 / np.sqrt(deg), 0.0)
    L_sym = np.eye(W.shape[0]) - (d_inv_sqrt[:, None] * W * d_inv_sqrt[None, :])
    L_sym = 0.5 * (L_sym + L_sym.T)
    _, evecs = np.linalg.eigh(L_sym)
    Y = evecs[:, 1:k + 1]                       # cols v_2 .. v_{k+1}
    norms = np.linalg.norm(Y, axis=1, keepdims=True)
    Y = Y / np.where(norms > 0, norms, 1.0)     # NJW row-normalization
    # Pairwise Euclidean distances (condensed upper triangle).
    sq = np.maximum(
        (Y * Y).sum(1)[:, None] + (Y * Y).sum(1)[None, :] - 2.0 * (Y @ Y.T),
        0.0)
    return _condensed_upper(np.sqrt(sq))


def dist_raw_fc(W: np.ndarray) -> np.ndarray:
    """ARM 3 — raw edge baseline D_raw[i,j] = 1 − |ImCoh|[i,j] (W ≥ 0 already)."""
    return _condensed_upper(1.0 - W)


def dist_diffusion_tau(eigvals_L: np.ndarray,
                       eigvecs_L: np.ndarray) -> np.ndarray:
    """ARM 4 — single-scale heat-kernel communication distance T = 1/ρ at
    τ=1/λmax, BEFORE the UPGMA hierarchy. This is the cophenetic arm's input
    distance with the multiscale UPGMA step REMOVED, isolating "diffusion at one
    scale" from "multiscale hierarchical integration" (the cophenetic arm = this
    + UPGMA). Same propagator construction as cophenetic_condensed_from_eigs."""
    lam_max = eigvals_L[-1]
    tau = 1.0 / lam_max
    diag_exp = np.exp(-tau * eigvals_L)
    rho = (eigvecs_L * diag_exp) @ eigvecs_L.T
    rho /= np.trace(rho)
    with np.errstate(divide="ignore"):
        Trho = 1.0 / rho
    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0.0)
    finite = np.isfinite(Trho)
    if not finite.all():
        cap = np.nanmax(Trho[finite]) if finite.any() else 1e6
        Trho = np.where(finite, Trho, cap)
    return _condensed_upper(Trho)


# ---------------------------------------------------------------------------
# Per-phase distance bundle: every arm's condensed distance vector for ONE
# adjacency, from its (cached or fresh) Laplacian eigendecomposition.
# ---------------------------------------------------------------------------
ARMS = ["cophenetic", "spectral", "raw_fc", "diffusion_tau"]


def distance_bundle(eigvals_L: np.ndarray, eigvecs_L: np.ndarray,
                    W: np.ndarray, k_grid: tuple[int, ...]) -> dict:
    """All arms' condensed per-pair distances for one phase adjacency."""
    bundle: dict = {
        "cophenetic": cophenetic_condensed_from_eigs(eigvals_L, eigvecs_L),
        "raw_fc": dist_raw_fc(W),
        "diffusion_tau": dist_diffusion_tau(eigvals_L, eigvecs_L),
        "spectral": {k: dist_spectral_embedding(eigvals_L, eigvecs_L, W, k)
                     for k in k_grid},
    }
    return bundle


def distance_bundle_from_W(W: np.ndarray, k_grid: tuple[int, ...]) -> dict:
    """Observed-side bundle: eigendecompose L = D − W, then `distance_bundle`."""
    L = np.diag(W.sum(axis=1)) - W
    eigvals, eigvecs = np.linalg.eigh(L)
    return distance_bundle(eigvals, eigvecs, W, k_grid)


# ---------------------------------------------------------------------------
# STEP 2 + 3 — IDENTICAL across arms: shifts then Spearman (audit_63 verbatim).
#   Δ_task = D^tt − D^pre_A ;  Δ_rest = D^post − D^pre_B ;
#   ρ_split = Spearman(Δ_task, Δ_rest)
# ---------------------------------------------------------------------------
def rho_split(d_pre_A: np.ndarray, d_pre_B: np.ndarray,
              d_tt: np.ndarray, d_post: np.ndarray) -> float:
    dD_task = d_tt - d_pre_A
    dD_rest = d_post - d_pre_B
    rho, _ = spearmanr(dD_task, dD_rest)
    return float(rho)


def _rho_for_arm(arm: str, k: int, B: dict[str, dict]) -> float:
    """ρ_split for one arm+k from the four-phase distance bundles."""
    def pull(phase: str):
        return B[phase][arm][k] if arm == "spectral" else B[phase][arm]
    return rho_split(pull("rest_pre_A"), pull("rest_pre_B"),
                     pull("task_test"), pull("rest_post"))


def arm_roster(k_grid: tuple[int, ...]) -> list[tuple[str, int]]:
    """(arm, k) pairs. cophenetic/raw_fc/diffusion_tau are k-free (k=0)."""
    roster: list[tuple[str, int]] = [
        ("cophenetic", 0), ("raw_fc", 0), ("diffusion_tau", 0)]
    roster += [("spectral", k) for k in k_grid]
    return roster


# ---------------------------------------------------------------------------
# Per-(patient, band): observed ρ_split + surrogate-tail p for every arm+k,
# reusing the cached audit_63 split-baseline ensemble for ALL arms.
# ---------------------------------------------------------------------------
def per_cell(pat: str, band: str, n_surr: int, swap_factor: int,
             seed: int, verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in PHASES}
    except Exception as e:
        if verbose:
            print(f"[audit_143] SKIP {pat}/{band}: {e}")
        return []
    N = Ws[PHASES[0]].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_143] SKIP {pat}/{band}: phase shape mismatch")
        return []

    roster = arm_roster(K_GRID)

    # Observed bundles + observed ρ_split per arm.
    obs_B = {ph: distance_bundle_from_W(Ws[ph], K_GRID) for ph in PHASES}
    Lc = obs_B[PHASES[0]]["cophenetic"].size
    if any(obs_B[ph]["cophenetic"].size != Lc for ph in PHASES):
        if verbose:
            print(f"[audit_143] SKIP {pat}/{band}: ultrametric size mismatch")
        return []
    obs_rho = {(arm, k): _rho_for_arm(arm, k, obs_B) for arm, k in roster}

    # Surrogate ensemble — load (cache HIT at seed=20260511) the audit_63
    # split-baseline matched-strength Laplacian eigs per phase, then derive
    # EVERY arm's surrogate ρ_split from the SAME draws.
    surr_eigs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    rngs = {ph: np.random.default_rng([seed, i])
            for i, ph in enumerate(PHASES)}
    for ph in PHASES:
        path = (SURROGATE_LRG_CACHE / pat
                / f"{band}_{ph}_R{n_surr}_swap{swap_factor}"
                  f"_seed{seed}_imcoh_abs.npz")
        evals, evecs = load_or_compute_eigs_at_path(
            path, Ws[ph], n_surr, swap_factor, rngs[ph], verbose=verbose)
        surr_eigs[ph] = (evals, evecs)

    surr_rho = {key: np.full(n_surr, np.nan) for key in obs_rho}
    for r in range(n_surr):
        B_r: dict[str, dict] = {}
        ok = True
        for ph in PHASES:
            evals, evecs = surr_eigs[ph]
            er, vr = evals[r], evecs[r]
            if not (np.isfinite(er).all() and np.isfinite(vr).all()):
                ok = False
                break
            W_r = adjacency_from_laplacian_eigs(er, vr)
            B_r[ph] = distance_bundle(er, vr, W_r, K_GRID)
        if not ok:
            continue
        for arm, k in roster:
            surr_rho[(arm, k)][r] = _rho_for_arm(arm, k, B_r)

    rows = []
    for arm, k in roster:
        s = surr_rho[(arm, k)]
        s_finite = s[np.isfinite(s)]
        if s_finite.size == 0:
            continue
        o = obs_rho[(arm, k)]
        p_upper = float(np.mean(s_finite >= o))   # team convention
        rows.append({
            "patient": pat, "band": band, "arm": arm, "K": int(k),
            "N_nodes": int(N), "n_pairs": int(Lc),
            "n_surrogates": int(s_finite.size),
            "obs_rho": float(o),
            "surr_rho_mean": float(np.mean(s_finite)),
            "surr_rho_p50": float(np.median(s_finite)),
            "surr_rho_p95": float(np.quantile(s_finite, 0.95)),
            "obs_p_upper": p_upper,
        })
    if verbose:
        rc = next((r for r in rows if r["arm"] == "cophenetic"), None)
        if rc:
            print(f"[audit_143] {pat}/{band}: coph ρ={rc['obs_rho']:+.3f} "
                  f"p={rc['obs_p_upper']:.3f}")
    return rows


# ---------------------------------------------------------------------------
# Cohort aggregation: per (band, arm) cohort Wilcoxon of (ρ_obs − surr_p50) > 0.
# For the swept spectral arm, report BEST-FOR-BASELINE k (smallest cohort p) so
# the baseline gets its strongest shot at recovering the trace.
# ---------------------------------------------------------------------------
def _cohort_for_subset(sub: pd.DataFrame) -> dict:
    contrast = (sub["obs_rho"] - sub["surr_rho_p50"]).to_numpy()
    z, p = wilcoxon_z(contrast)
    return {
        "n_patients": len(sub),
        "n_pos_obs_rho": int((sub["obs_rho"] > 0).sum()),       # sign /10
        "n_above_surrogate": int((sub["obs_p_upper"] < 0.05).sum()),
        "obs_median_rho": float(np.median(sub["obs_rho"])),
        "surr_median_rho": float(np.median(sub["surr_rho_p50"])),
        "wilcoxon_z": float(z),
        "wilcoxon_p": float(p),
    }


def cohort_summary(per_pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in BANDS:
        for arm in ARMS:
            fam = per_pat[(per_pat.band == band) & (per_pat.arm == arm)]
            if fam.empty:
                continue
            if arm == "spectral":
                best, best_K = None, None
                for k in sorted(fam.K.unique()):
                    s = fam[fam.K == k]
                    if len(s) < 3:
                        continue
                    res = _cohort_for_subset(s)
                    if best is None or res["wilcoxon_p"] < best["wilcoxon_p"]:
                        best, best_K = res, int(k)
                if best is None:
                    continue
            else:
                best, best_K = _cohort_for_subset(fam), 0
            row = {"band": band, "arm": arm, "best_K": best_K, **best}
            if arm == "cophenetic":
                row["audit63_locked_p"] = AUDIT63_LOCKED_P.get(band, np.nan)
            out.append(row)
    df = pd.DataFrame(out)
    # BH-FDR WITHIN each arm across the 6-band family.
    df["wilcoxon_q_bh"] = np.nan
    for arm in df["arm"].unique():
        mask = df["arm"] == arm
        df.loc[mask, "wilcoxon_q_bh"] = bh_fdr(
            df.loc[mask, "wilcoxon_p"].tolist())
    df["detects_trace"] = df["wilcoxon_p"] < 0.05
    df["detects_trace_q"] = df["wilcoxon_q_bh"] < 0.05
    return df


# ---------------------------------------------------------------------------
# Correctness gate: cophenetic arm vs locked audit_63 split-baseline p (α, β).
# ---------------------------------------------------------------------------
def gate_report(cohort: pd.DataFrame) -> tuple[bool, list[str]]:
    lines, ok_all = [], True
    coph = cohort[cohort.arm == "cophenetic"].set_index("band")
    for band in GATE_BANDS:
        locked = AUDIT63_LOCKED_P[band]
        if band not in coph.index:
            lines.append(f"  {band:<6} MISSING from cophenetic arm — GATE FAIL")
            ok_all = False
            continue
        got = float(coph.loc[band, "wilcoxon_p"])
        ok = abs(got - locked) <= GATE_TOL
        ok_all &= ok
        lines.append(
            f"  {band:<6} cophenetic p={got:.6f}  locked={locked:.6f}  "
            f"|Δ|={abs(got - locked):.2e}  {'PASS' if ok else 'FAIL'}")
    return ok_all, lines


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(cohort: pd.DataFrame, gate_ok: bool, gate_lines: list[str],
                 n_surr: int, swap_factor: int, seed: int,
                 runtime_s: float) -> None:
    arm_label = {"cophenetic": "cophenetic (multiscale UPGMA)",
                 "spectral": "leading-k spectral embedding",
                 "raw_fc": "raw FC (1−|ImCoh|)",
                 "diffusion_tau": "single-τ diffusion distance"}
    lines = [
        "---",
        "name: spectral_distance_swap",
        "scope: headtohead_rho_split_distance_swap_cophenetic_vs_spectral_embedding_vs_rawfc",
        "date: 2026-06-23",
        "status: first_pass",
        "---",
        "",
        "# Distance-swap head-to-head on the team's OWN ρ_split statistic (v2)",
        "",
        "**Head.** Holding the team's split-baseline ρ_split = "
        "Spearman(Δ_task, Δ_rest) FIXED and swapping ONLY the per-pair "
        "distance, does the FULL multiscale **cophenetic** distance recover "
        "the α/β persistence that the pre-cut **leading-k spectral embedding** "
        "(the representation spectral clustering runs k-means on, at any "
        f"k∈{{{', '.join(map(str, K_GRID))}}}) and the **raw-FC** edge baseline "
        "cannot, under the audit_63 matched-strength null? v1 (audit_141) was "
        "inconclusive because it used a distance TRIANGLE + ARI/NMI partition "
        "similarity — the wrong common footing; this v2 uses the REAL "
        "ρ_split shift-correlation, reproduces the locked numbers on the "
        "cophenetic arm as a gate, then swaps the distance.",
        "",
        "## Correctness gate (cophenetic arm vs locked audit_63 split-baseline)",
        "",
        f"**GATE: {'PASS' if gate_ok else 'FAIL'}** "
        f"(tolerance |Δp| ≤ {GATE_TOL}; same cached surrogate ensemble, "
        f"seed={seed}).",
        "",
        "```",
        *gate_lines,
        "```",
        "",
        "## Band × distance-arm cohort table",
        "",
        "ρ_split per-patient → cohort one-sided Wilcoxon of "
        "(ρ_obs − surr_median) > 0; matched-strength upper-tail p; BH-FDR "
        "within each arm across the 6-band family. Spectral arm reports the "
        "BEST-FOR-BASELINE k.",
        "",
        "| band | arm | best k | sign (#ρ>0) | n>surr | obs median ρ | "
        "Wilcoxon p | q (BH) | detects? |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in cohort.iterrows():
        lines.append(
            f"| {r['band']} | {arm_label[r['arm']]} | "
            f"{r['best_K'] if r['arm'] == 'spectral' else '—'} "
            f"| {r['n_pos_obs_rho']}/{r['n_patients']} "
            f"| {r['n_above_surrogate']}/{r['n_patients']} "
            f"| {r['obs_median_rho']:+.3f} "
            f"| {r['wilcoxon_p']:.4f} "
            f"| {r['wilcoxon_q_bh']:.4f} "
            f"| {'YES' if r['detects_trace'] else 'no'} |")
    lines.extend([
        "",
        "## Distance-arm definitions (STEP 1 — the only thing that varies)",
        "",
        "Steps 2–3 are byte-identical across arms: "
        "Δ_task = D^tt − D^pre_A, Δ_rest = D^post − D^pre_B, "
        "ρ_split = Spearman(Δ_task, Δ_rest) over the upper-triangle pairs "
        "(disjoint rest_pre halves A/B — the audit_63 construction).",
        "",
        "- **cophenetic** — D_coph = cophenet(UPGMA(T(τ=1/λmax))); the team's "
        "multiscale read-out. The number of record.",
        "- **spectral (leading-k)** — D_spec^(k)[i,j] = ||Y_i − Y_j||₂, "
        "Y = leading-k non-trivial normalized-Laplacian eigenvectors "
        "(cols v₂..v_{k+1}), Ng-Jordan-Weiss row-normalized — exactly the "
        f"pre-cut spectral-clustering representation. Swept k∈{{{', '.join(map(str, K_GRID))}}}, "
        "best-for-baseline k reported.",
        "- **raw FC** — D_raw[i,j] = 1 − |ImCoh|[i,j]; edge-wise, no diffusion, "
        "no hierarchy.",
        "- **single-τ diffusion** — D_diff[i,j] = (1/ρ)[i,j] at τ=1/λmax (heat-"
        "kernel communication distance BEFORE UPGMA); separates one-scale "
        "diffusion from multiscale hierarchical integration.",
        "",
        "## Reading the verdict",
        "",
        "Thesis SURVIVES if **cophenetic** detects the α trace (clears matched-"
        "strength, ideally BH-q<0.05) while the **leading-k spectral embedding** "
        "MISSES α at every k AND **raw FC** misses it. β is the positive control "
        "(cophenetic expected strongest; raw FC weakest). The demonstration is "
        "WEAKENED if the spectral embedding catches α at any k — reported "
        "plainly; the baseline is given its best k and the embedding it "
        "actually uses.",
        "",
        "## Provenance",
        "",
        f"- N_surrogates = {n_surr}; n_swaps = {swap_factor}·N(N−1)/2; "
        f"seed = {seed} (REUSES the cached audit_63 split-baseline ensemble)",
        "- Cohort: " + ", ".join(COHORT),
        "- Bands: " + ", ".join(BANDS),
        "- Phases: rest_pre_A, rest_pre_B, task_test, rest_post (split-baseline)",
        "- FC method: imcoh_abs",
        "- k grid: " + ", ".join(str(k) for k in K_GRID),
        f"- Wall-clock runtime: {runtime_s:.1f} s",
        "- Build: scripts/01_compute/audit/audit_143_spectral_distance_swap_headtohead.py",
        "",
        "## Caveats (referee-facing)",
        "",
        "- The matched-strength null is INDEPENDENT per phase (the audit_63 "
        "family, shared identically across all four arms) — it controls per-node "
        "strength heterogeneity but does NOT isolate edge-identity cross-phase "
        "structure for ANY arm. 'Clears the null' = exceeds an independent-"
        "strength-randomized baseline. Fair because every arm gets the same null.",
        "- The spectral-embedding distance ||Y_i − Y_j|| is sign/rotation-"
        "ambiguous in the eigenvectors, but the ambiguity cancels WITHIN a "
        "phase's distance matrix, so the per-pair SHIFT and ρ_split are "
        "well-defined without cross-phase alignment.",
        "- raw FC and single-τ diffusion are single-view (unhierarchical) "
        "baselines by construction — that IS the contrast against the "
        "multiscale cophenetic arm.",
        "",
    ])
    (OUT / "README.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=BANDS)
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--n-surrogates", type=int, default=N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=SWAP_FACTOR)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--gate-only", action="store_true",
                    help="run α+β only, print the correctness gate, exit.")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    bands = list(GATE_BANDS) if args.gate_only else args.bands

    # Pre-flight: ensure half FCs cached (audit_63 helper).
    for pat in args.patients:
        ensure_half_fcs(pat, bands)

    t0 = time.time()
    all_rows: list[dict] = []
    for band in bands:
        for pat in args.patients:
            t_cell = time.time()
            rows = per_cell(pat, band, args.n_surrogates, args.swap_factor,
                            args.seed, args.verbose)
            all_rows.extend(rows)
            if rows:
                print(f"[audit_143] {pat}/{band}: {len(rows)} arm-cells "
                      f"({time.time() - t_cell:.1f}s)")

    runtime = time.time() - t0
    per_pat = pd.DataFrame(all_rows)
    cohort = cohort_summary(per_pat)
    gate_ok, gate_lines = gate_report(cohort)

    print("\n=== Correctness gate (cophenetic vs locked audit_63) ===")
    print("\n".join(gate_lines))
    print(f"GATE: {'PASS' if gate_ok else 'FAIL'}\n")

    if args.gate_only:
        print(cohort[cohort.arm == "cophenetic"].to_string(index=False))
        print(f"\n[audit_143] gate-only {runtime:.1f}s")
        return

    per_pat.to_csv(OUT / "per_patient_per_band.csv", index=False)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)
    print(cohort.to_string(index=False))
    write_readme(cohort, gate_ok, gate_lines, args.n_surrogates,
                 args.swap_factor, args.seed, runtime)
    print(f"\n[audit_143] total {runtime:.1f}s; outputs at {OUT}")


if __name__ == "__main__":
    main()
