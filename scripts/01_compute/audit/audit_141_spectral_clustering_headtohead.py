#!/usr/bin/env python3
"""Audit 141 — ρ^coph vs off-the-shelf spectral/PCA: head-to-head trace detection.

==========================================================================
MANDATORY 5-POINT CRITICAL PREAMBLE
==========================================================================

(1) CLAIM UNDER TEST.
    The paper's load-bearing methodological claim is that the team's
    multiscale Laplacian read-out — the cophenetic per-pair distance
    ρ^coph (UPGMA on the heat-kernel ultrametric at τ = 1/λ_max) —
    detects a cross-phase "persistence trace" that a NAIVE spectral
    analysis of the SAME adjacency MISSES. The sharpest case is the
    ALPHA band: the locked verdict is α carries the trace on ρ^coph
    (audit_63 split-baseline cohort p = 0.00195) but NOT on the
    global-eigenmode Grassmann probe. Prediction: off-the-shelf
    spectral clustering (k-means on Laplacian eigenvectors → ARI/NMI
    partition distance) and a PCA/leading-eigenmode subspace baseline
    should FAIL to recover the α trace that ρ^coph recovers. β (trace
    on BOTH probes) should be partly visible to the spectral baseline;
    θ should be null everywhere; δ / γ_low (Grassmann-flavoured) should
    be visible to the PCA/subspace baseline since that ≈ the Grassmann
    probe.

    Trace statistic (team convention, sign LOCKED):
        T = d(rest_pre, task_test) − d(task_test, rest_post),  T>0 = trace.
    Per patient T per method; cohort = one-sided Wilcoxon
    (alternative='greater'); per-method cohort p vs a matched-strength
    surrogate upper-tail null.

(2) THE NULL.
    Matched-strength surrogate adjacencies (4-cycle ±δ strength-
    preserving rewiring, n_swaps = 20·N(N−1)/2, R=200) per phase,
    INDEPENDENT per phase — exactly the audit_63 / audit_66 family.
    For each method we recompute its per-patient T on the surrogate
    ensemble and take the per-patient upper-tail p = mean(s ≥ T_obs);
    cohort verdict = one-sided Wilcoxon of (T_obs − surr_median) > 0.

(3) STRONGEST ALTERNATIVE THE NULL MUST CONTROL FOR.
    The dominant confound for ANY cross-phase FC statistic is per-node
    STRENGTH heterogeneity drift (the KC-β collapse of 2026-05-11):
    a method can look like it detects a "trace" purely because node
    strengths are auto-correlated across phases. The matched-strength
    null reproduces each phase's exact strength sequence while
    scrambling edge identity, so a surviving T cannot be attributed to
    strength alone. This is the team's mandatory minimum null for any
    cohort FC claim.

(4) DOES THE NULL ACTUALLY CONTROL FOR IT — BY MECHANISM.
    Yes for the strength axis: 4-cycle ±δ preserves the strength
    sequence exactly (verify_strengths tol 1e-4), so per-method T on
    surrogates isolates strength-independent structure. WHAT IT CANNOT
    DO: it is INDEPENDENT per phase, so it destroys cross-phase backbone
    coherence in addition to per-phase shape. Therefore a "clears the
    null" verdict here means "T exceeds an independent-strength-
    randomized baseline", NOT "specific edge-identity structure carries
    the cross-phase trace". This is identical to the caveat the team
    already documented for audit_63; it is the same null, so the
    head-to-head is APPLES-TO-APPLES across the three methods.
    SECOND caveat specific to this audit: ARI / NMI partition distances
    are normally FORBIDDEN for the team's OWN results
    (feedback_no_partition_metrics_use_rho_coph) because at any fixed K
    one cluster dominates and the partition distance is dominated by
    boundary churn. We use them HERE ONLY to instantiate the strawman
    spectral-clustering baseline we aim to beat — they characterize the
    naive method, they are not a team result. We sweep K and report the
    BEST-CASE-for-the-baseline K so the strawman gets its strongest shot.

(5) WHAT WOULD FALSIFY THE CLAIM / REMAINING LIMITATIONS.
    The claim is FALSIFIED if a fair spectral-clustering or PCA baseline
    recovers the α trace at cohort p < 0.05 (matched-strength) at any
    swept K — that would mean ρ^coph is not adding detection power over
    off-the-shelf spectral methods in the sharpest test case, weakening
    the paper's thesis. We give the baseline every reasonable advantage
    (normalized AND unnormalized Laplacian; K∈{2,3,4,5,6,8,10}; ARI AND
    NMI; PCA chordal-subspace sweep) and report the most favourable
    cell. Remaining limitations: (a) the surrogate is independent-per-
    phase (shared with the reference, so fair, but neither isolates
    edge-identity structure); (b) k-means has random init — we fix the
    seed and use n_init=4 — validated to give the same cross-phase ARI/NMI
    as n_init=10 to ~3rd decimal on Pat_02/α, with k-means the sole runtime
    bottleneck; (c) we compare on the FULL-rest_pre triangle
    for clean common ground, and separately cite the locked audit_63
    split-baseline ρ^coph p so the reference column is the team's actual
    locked number; (d) ARI/NMI are discrete/quantized so the surrogate
    tail can be coarse — reported honestly.

==========================================================================

Outputs
-------
    data/audit/spectral_headtohead/
        cohort_summary.csv          # band × method → trace verdict
        per_patient_per_band.csv    # every (patient, band, method, K) cell
        README.md

Library reuse (no private forks)
--------------------------------
- FC loading                : workflow.fc.load_fc_matrix
- matched-strength surrogate : utils.surrogate.matched_strength
                               (load_or_compute_eigs_at_path — caches the
                                expensive 4-cycle ±δ shuffle + eigh to disk;
                                adjacency_from_laplacian_eigs;
                                cophenetic_condensed_from_eigs /
                                cophenetic_condensed_from_adjacency)
- chordal Grassmann distance : utils.metrics.spectral.chordal_distance
- cohort stats               : utils.metrics.hypothesis.wilcoxon_z, bh_fdr

Surrogate cache
---------------
The strength-preserving shuffle (~290 ms × 3 phases × R) dominates runtime,
so the surrogate Laplacian eigendecompositions are cached once per
(patient, band, phase) at
``data/cache/matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R{R}_swap{SF}_seed{S}_imcoh_abs.npz``
via the library ``load_or_compute_eigs_at_path``. Every method's surrogate
statistic is then derived from the cached eigs (adjacency reconstructed for
the spectral/PCA methods; cophenetic distance straight from the eigs), so
re-runs are I/O-bound. NOTE: this audit uses the FULL ``rest_pre`` phase, so
its ``rest_pre`` cache file is DISTINCT from the audit_66 split-baseline
``rest_pre_A`` ensemble even at the same (R, SF, seed); ``task_test`` /
``rest_post`` could in principle collide with audit_66 only if the same seed
AND rng-draw order were used — they are not (seed 20260622 here vs 20260511).
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr, wilcoxon_z
from lrg_eegfc.utils.metrics.spectral import chordal_distance
from lrg_eegfc.utils.surrogate.matched_strength import (
    adjacency_from_laplacian_eigs,
    cophenetic_condensed_from_adjacency,
    cophenetic_condensed_from_eigs,
    load_or_compute_eigs_at_path,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

SURROGATE_LRG_CACHE = CACHE_ROOT / "matched_strength_surrogate_lrg"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
# delta, gamma_low = Grassmann-flavoured controls; alpha = PRIMARY;
# beta = positive control (both probes); theta = negative control.
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PHASES = ("rest_pre", "task_test", "rest_post")
K_GRID = (2, 3, 4, 5, 6, 8, 10)
KMEANS_N_INIT = 4   # k-means restarts (verdict stable from 4↑; see phase_bundle)
N_SURROGATES = 200
SWAP_FACTOR = 20
SEED = 20260622

# Locked audit_63 split-baseline ρ^coph cohort p-values (the team's
# actual reference numbers; read from
# data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv).
# Cited as the reference column; NOT recomputed here (the team's number
# of record). low_gamma / theta etc. included for completeness.
RHO_COPH_REF_AUDIT63 = {
    "delta": 0.2783, "theta": 0.7217, "alpha": 0.001953,
    "beta": 0.004883, "low_gamma": 0.1162, "high_gamma": 0.2461,
}
RHO_COPH_REF_NABOVE = {
    "delta": "4/10", "theta": "2/10", "alpha": "5/10",
    "beta": "7/10", "low_gamma": "5/10", "high_gamma": "4/10",
}

OUT = ROOT / "data" / "audit" / "spectral_headtohead"
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Phase FC loading (full phases, symmetrized, clipped) — same idiom as audit_63
# ---------------------------------------------------------------------------
def load_phase_fc(pat: str, phase: str, band: str) -> np.ndarray:
    W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
    if W is None:
        raise FileNotFoundError(f"no cached FC for {pat}/{phase}/{band}")
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


# ---------------------------------------------------------------------------
# Method 1 — cophenetic ρ^coph distance (recomputed on the FULL triangle so it
# shares inputs with the baselines). d = 1 − Spearman(D_coph_A, D_coph_B);
# D_coph via `cophenetic_condensed_from_adjacency` (UPGMA on the heat-kernel
# ultrametric at τ=1/λ_max). Assembled in `phase_bundle` / `_T_from_bundles`.
# The LOCKED reference p is the audit_63 split-baseline ρ_split.

# ---------------------------------------------------------------------------
# Method 2 — spectral-clustering partition distance (the strawman baseline).
# k-means on the leading-K eigenvectors of L_sym (normalized) or L (unnorm).
# Cross-phase distance = 1 − ARI or 1 − NMI between partitions.
# ARI/NMI are FORBIDDEN for team results — used here only to characterize
# the off-the-shelf baseline we aim to beat (see preamble point 4).
# ---------------------------------------------------------------------------
def _spectral_embedding(W: np.ndarray, k: int, normalized: bool) -> np.ndarray:
    """Leading-k non-trivial eigenvectors of the (normalized) Laplacian."""
    deg = W.sum(axis=1)
    if normalized:
        with np.errstate(divide="ignore"):
            d_inv_sqrt = np.where(deg > 0, 1.0 / np.sqrt(deg), 0.0)
        L = np.eye(W.shape[0]) - (d_inv_sqrt[:, None] * W * d_inv_sqrt[None, :])
        L = 0.5 * (L + L.T)
    else:
        L = np.diag(deg) - W
    eigvals, eigvecs = np.linalg.eigh(L)
    # Skip the trivial smallest mode (index 0), take next k.
    return eigvecs[:, 1:k + 1]


# The spectral-clustering embedding is row-normalized (Ng-Jordan-Weiss) for
# the normalized Laplacian inside `phase_bundle`; the unnormalized variant
# uses the raw Fiedler block. Method 3 (PCA) uses the top-K adjacency
# eigenvectors → chordal Grassmann distance (≈ the Grassmann probe), also
# assembled in `phase_bundle`.


# ---------------------------------------------------------------------------
# Per-phase spectral bundle.
# For one adjacency W, compute ONCE all the per-method/per-K ingredients:
#   - cophenetic condensed distance (rho_coph)
#   - top-Kmax adjacency eigenvectors (pca_chordal)
#   - spectral-clustering partitions at every K for L_sym and L_unnorm
# so every method+K downstream reads precomputed objects (no redundant eigh /
# KMeans across the ~29 method-cells × R surrogates).
# ---------------------------------------------------------------------------
def phase_bundle(W: np.ndarray, k_grid: tuple[int, ...], seed: int,
                 coph: np.ndarray | None = None) -> dict:
    """All per-method/per-K ingredients for one adjacency.

    ``coph`` may be passed precomputed (from cached Laplacian eigs via
    `cophenetic_condensed_from_eigs`) to skip the redundant eigh; otherwise it
    is computed from ``W``.
    """
    kmax = max(k_grid)
    if coph is None:
        coph = cophenetic_condensed_from_adjacency(W)
    # PCA: top-kmax adjacency eigenvectors (largest-magnitude → last columns)
    _, w_evecs = np.linalg.eigh(W)
    pca_V = {k: np.ascontiguousarray(w_evecs[:, -k:]) for k in k_grid}
    # Spectral-clustering partitions per (norm, K). n_init=KMEANS_N_INIT
    # (4): the cross-phase ARI/NMI verdict is stable from n_init=4 upward
    # (validated on Pat_02/α — n_init 10 vs 4 agree to ~3rd decimal),
    # and KMeans is the sole runtime bottleneck (everything else <0.1s/cell),
    # so 4 restarts keeps the strawman fair while halving cost. Applied to
    # BOTH observed and surrogate sides for consistency.
    parts: dict[tuple[str, int], np.ndarray] = {}
    for norm_tag, normd in (("sym", True), ("unnorm", False)):
        emb_full = _spectral_embedding(W, kmax, normd)  # leading kmax modes
        for k in k_grid:
            emb = emb_full[:, :k]
            if normd:
                norms = np.linalg.norm(emb, axis=1, keepdims=True)
                emb = emb / np.where(norms > 0, norms, 1.0)
            km = KMeans(n_clusters=k, n_init=KMEANS_N_INIT, random_state=seed)
            parts[(norm_tag, k)] = km.fit_predict(emb)
    return {"coph": coph, "pca_V": pca_V, "parts": parts}


def bundle_from_laplacian_eigs(eigvals: np.ndarray, eigvecs: np.ndarray,
                               k_grid: tuple[int, ...], seed: int) -> dict:
    """Phase bundle from a cached unnormalized-Laplacian eigendecomposition.

    Reconstructs the adjacency ``W`` from ``L = D − W`` eigs
    (`adjacency_from_laplacian_eigs`) for the spectral / PCA methods, and
    computes the cophenetic distance straight from the eigs
    (`cophenetic_condensed_from_eigs`, no re-eigh). Used for the cached
    surrogate ensemble so the expensive 4-cycle ±δ shuffle never re-runs.
    """
    coph = cophenetic_condensed_from_eigs(eigvals, eigvecs)
    W = adjacency_from_laplacian_eigs(eigvals, eigvecs)
    return phase_bundle(W, k_grid, seed, coph=coph)


def _T_from_bundles(method: str, k: int,
                    Bp: dict, Bt: dict, Bq: dict) -> float:
    """Triangle T = d(pre, tt) − d(tt, post) for one method+K from bundles."""
    if method == "rho_coph":
        d_pre_tt = 1.0 - float(spearmanr(Bp["coph"], Bt["coph"])[0])
        d_tt_post = 1.0 - float(spearmanr(Bt["coph"], Bq["coph"])[0])
    elif method == "pca_chordal":
        d_pre_tt = chordal_distance(Bp["pca_V"][k], Bt["pca_V"][k])
        d_tt_post = chordal_distance(Bt["pca_V"][k], Bq["pca_V"][k])
    elif method.startswith("spec_"):
        # spec_{sym,unnorm}_{ari,nmi}
        _, norm_tag, metric = method.split("_")
        sim_fn = (adjusted_rand_score if metric == "ari"
                  else normalized_mutual_info_score)
        d_pre_tt = 1.0 - float(sim_fn(Bp["parts"][(norm_tag, k)],
                                      Bt["parts"][(norm_tag, k)]))
        d_tt_post = 1.0 - float(sim_fn(Bt["parts"][(norm_tag, k)],
                                       Bq["parts"][(norm_tag, k)]))
    else:
        raise ValueError(method)
    return float(d_pre_tt - d_tt_post)


def method_roster(k_grid: tuple[int, ...]) -> list[tuple[str, int]]:
    """(method_name, K) pairs. rho_coph is K-free (K=0)."""
    roster: list[tuple[str, int]] = [("rho_coph", 0)]
    for k in k_grid:
        for norm_tag in ("sym", "unnorm"):
            for metric in ("ari", "nmi"):
                roster.append((f"spec_{norm_tag}_{metric}", k))
    for k in k_grid:
        roster.append(("pca_chordal", k))
    return roster


# ---------------------------------------------------------------------------
# Per-(patient, band) — compute obs T and surrogate-tail p for every method.
# Surrogates: independent per-phase matched-strength rewiring (audit_63
# family), generated once and shared across all methods.
# ---------------------------------------------------------------------------
def per_cell(pat: str, band: str, n_surr: int, swap_factor: int,
             seed: int, verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in PHASES}
    except Exception as e:
        if verbose:
            print(f"[audit_141] SKIP {pat}/{band}: {e}")
        return []
    N = Ws[PHASES[0]].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_141] SKIP {pat}/{band}: shape mismatch")
        return []

    roster = method_roster(K_GRID)

    # Observed bundles + observed T per method
    obs_B = {ph: phase_bundle(Ws[ph], K_GRID, seed) for ph in PHASES}
    obs_T = {
        (name, k): _T_from_bundles(name, k, obs_B["rest_pre"],
                                   obs_B["task_test"], obs_B["rest_post"])
        for name, k in roster
    }

    # Surrogate ensemble — load (or generate-and-cache) the matched-strength
    # Laplacian eigendecompositions per phase, then derive every method's
    # surrogate T from the cached eigs. The expensive 4-cycle ±δ shuffle runs
    # at most once per (patient, band, phase) and is reused on re-runs.
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

    surr_T = {key: np.full(n_surr, np.nan) for key in obs_T}
    for r in range(n_surr):
        B_r: dict[str, dict] = {}
        ok = True
        for ph in PHASES:
            evals, evecs = surr_eigs[ph]
            if not np.isfinite(evals[r]).all() or not np.isfinite(evecs[r]).all():
                ok = False
                break
            B_r[ph] = bundle_from_laplacian_eigs(evals[r], evecs[r],
                                                 K_GRID, seed)
        if not ok:
            continue
        for name, k in roster:
            surr_T[(name, k)][r] = _T_from_bundles(
                name, k, B_r["rest_pre"], B_r["task_test"], B_r["rest_post"])

    rows = []
    for name, k in roster:
        s = surr_T[(name, k)]
        s_finite = s[np.isfinite(s)]
        if s_finite.size == 0:
            continue
        o = obs_T[(name, k)]
        # Team convention: upper-tail p = mean(s_finite >= obs_T).
        p_upper = float(np.mean(s_finite >= o))
        rows.append({
            "patient": pat, "band": band, "method": name, "K": int(k),
            "N_nodes": int(N), "n_surrogates": int(s_finite.size),
            "obs_T": float(o),
            "surr_T_mean": float(np.mean(s_finite)),
            "surr_T_p50": float(np.median(s_finite)),
            "surr_T_p95": float(np.quantile(s_finite, 0.95)),
            "obs_p_upper": p_upper,
        })
    if verbose:
        rc = next((r for r in rows if r["method"] == "rho_coph"), None)
        if rc:
            print(f"[audit_141] {pat}/{band}: rho_coph T={rc['obs_T']:+.3f} "
                  f"p={rc['obs_p_upper']:.3f}")
    return rows


# ---------------------------------------------------------------------------
# Cohort aggregation: per (band, method-family) cohort Wilcoxon of
# (T_obs − surr_p50) > 0. For the swept methods, report BEST-FOR-BASELINE K
# (smallest cohort p) so the strawman gets its strongest shot.
# ---------------------------------------------------------------------------
METHOD_FAMILIES = [
    ("rho_coph", "rho_coph"),                 # reference (recomputed triangle)
    ("spec_sym_ari", "spec_sym_ari"),
    ("spec_sym_nmi", "spec_sym_nmi"),
    ("spec_unnorm_ari", "spec_unnorm_ari"),
    ("spec_unnorm_nmi", "spec_unnorm_nmi"),
    ("pca_chordal", "pca_chordal"),
]


def _cohort_for_subset(sub: pd.DataFrame) -> dict:
    """Cohort Wilcoxon of (obs_T − surr_p50) across patients for one cell."""
    contrast = (sub["obs_T"] - sub["surr_T_p50"]).to_numpy()
    z, p = wilcoxon_z(contrast)
    n_above = int((sub["obs_p_upper"] < 0.05).sum())
    n_pos = int((sub["obs_T"] > 0).sum())
    return {
        "n_patients": len(sub),
        "n_pos_obs_T": n_pos,                 # sign count /10
        "n_above_surrogate": n_above,
        "obs_median_T": float(np.median(sub["obs_T"])),
        "surr_median_T": float(np.median(sub["surr_T_p50"])),
        "wilcoxon_z": float(z),
        "wilcoxon_p": float(p),
    }


def cohort_summary(per_pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in BANDS:
        for fam, fam_label in METHOD_FAMILIES:
            fam_rows = per_pat[(per_pat.band == band)
                               & (per_pat.method == fam)]
            if fam_rows.empty:
                continue
            if fam == "rho_coph":
                best = _cohort_for_subset(fam_rows)
                best_K = 0
            else:
                # Sweep K; pick the K with the SMALLEST cohort p
                # (best-case for the baseline).
                best = None
                best_K = None
                for k in sorted(fam_rows.K.unique()):
                    sub = fam_rows[fam_rows.K == k]
                    if len(sub) < 3:
                        continue
                    res = _cohort_for_subset(sub)
                    if (best is None
                            or res["wilcoxon_p"] < best["wilcoxon_p"]):
                        best = res
                        best_K = int(k)
                if best is None:
                    continue
            row = {
                "band": band,
                "method": fam_label,
                "best_K": best_K,
                **best,
            }
            # Attach the locked audit_63 ρ^coph reference for context.
            if fam == "rho_coph":
                row["audit63_split_p"] = RHO_COPH_REF_AUDIT63.get(band, np.nan)
                row["audit63_split_n_above"] = RHO_COPH_REF_NABOVE.get(band, "")
            out.append(row)
    df = pd.DataFrame(out)
    # BH-FDR WITHIN each method across the 6 bands (the multiplicity family:
    # "does method X detect a trace in band Y, correcting for testing 6
    # bands"). FDR-ing reference + strawman together would dilute, so each
    # method gets its own 6-band family.
    df["wilcoxon_q_bh"] = np.nan
    for fam_label in df["method"].unique():
        mask = df["method"] == fam_label
        df.loc[mask, "wilcoxon_q_bh"] = bh_fdr(
            df.loc[mask, "wilcoxon_p"].tolist())
    df["detects_trace"] = (df["wilcoxon_p"] < 0.05)
    df["detects_trace_q"] = (df["wilcoxon_q_bh"] < 0.05)
    return df


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(cohort: pd.DataFrame, n_surr: int, swap_factor: int,
                 runtime_s: float) -> None:
    lines = [
        "---",
        "name: spectral_headtohead",
        "scope: methodological_headtohead_rhocoph_vs_spectral_pca_trace_detection",
        "date: 2026-06-22",
        "status: first_pass",
        "---",
        "",
        "# ρ^coph vs off-the-shelf spectral / PCA — trace detection head-to-head",
        "",
        "**Head.** Does the team's cophenetic ρ^coph read-out detect a "
        "cross-phase persistence trace that a NAIVE spectral-clustering "
        "(ARI/NMI on Laplacian-eigenvector k-means) or PCA-subspace "
        "(chordal Grassmann on adjacency eigenvectors) baseline MISSES? "
        "Sharpest case = the ALPHA band (ρ^coph-only per the locked "
        "verdict). Each method's per-patient triangle "
        "`T = d(rest_pre, task_test) − d(task_test, rest_post)` (T>0=trace) "
        f"is tested against an R={n_surr} matched-strength surrogate null "
        "(independent per-phase 4-cycle ±δ rewiring, the audit_63 family); "
        "cohort verdict = one-sided Wilcoxon of (T_obs − surr_median) > 0. "
        "Swept baselines report the BEST-CASE-for-the-baseline K. "
        "**ARI/NMI are used ONLY to characterize the strawman baseline** "
        "(team results never use partition-cut metrics; see "
        "`feedback_no_partition_metrics_use_rho_coph`).",
        "",
        "## Band × method cohort table",
        "",
        "| band | method | best K | sign count (#T>0) | n above surrogate | "
        "obs median T | cohort Wilcoxon p | q (BH) | detects trace? |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in cohort.iterrows():
        lines.append(
            f"| {r['band']} | {r['method']} | {r['best_K']} "
            f"| {r['n_pos_obs_T']}/{r['n_patients']} "
            f"| {r['n_above_surrogate']}/{r['n_patients']} "
            f"| {r['obs_median_T']:+.3f} "
            f"| {r['wilcoxon_p']:.4f} "
            f"| {r['wilcoxon_q_bh']:.4f} "
            f"| {'YES' if r['detects_trace'] else 'no'} |"
        )
    lines.extend([
        "",
        "## Locked ρ^coph reference (audit_63 split-baseline)",
        "",
        "The `rho_coph` rows above are the FULL-rest_pre triangle recomputed "
        "here so all three methods share identical inputs. The team's "
        "number of record is the split-baseline ρ_split from "
        "`matched_strength_surrogate_split_baseline` (audit_63):",
        "",
        "| band | audit_63 split-baseline p | n above own surrogate |",
        "|---|---|---|",
    ])
    for band in BANDS:
        lines.append(
            f"| {band} | {RHO_COPH_REF_AUDIT63.get(band, float('nan')):.4f} "
            f"| {RHO_COPH_REF_NABOVE.get(band, '')} |")
    lines.extend([
        "",
        "## Method definitions",
        "",
        "- **rho_coph** — d = 1 − Spearman(D_coph_A[triu], D_coph_B[triu]); "
        "D_coph from UPGMA on the heat-kernel ultrametric at τ=1/λ_max "
        "(`cophenetic_condensed_from_adjacency`). The team's read-out.",
        "- **spec_{sym,unnorm}_{ari,nmi}** — k-means (n_init=4) on the "
        "leading-K eigenvectors of the symmetric-normalized (Ng-Jordan-Weiss, "
        "row-normalized embedding) or unnormalized Laplacian; cross-phase "
        "d = 1 − ARI (or 1 − NMI). Swept K∈{2,3,4,5,6,8,10}, best-for-baseline "
        "K reported. **Strawman baseline only.**",
        "- **pca_chordal** — top-K eigenvectors of the adjacency W; "
        "cross-phase d = chordal Grassmann distance "
        "(`utils.metrics.spectral.chordal_distance`). ≈ the Grassmann probe.",
        "",
        "## Reading the verdict",
        "",
        "The paper's thesis SURVIVES if `rho_coph` detects the α trace "
        "(p<0.05) while every spectral-clustering and PCA baseline MISSES it "
        "(p≥0.05 at every K). The thesis is WEAKENED if any fair baseline "
        "recovers the α trace. β is expected to be partly visible to the "
        "spectral baseline (trace on both probes); θ null everywhere; "
        "δ / γ_low may be visible to pca_chordal (≈ Grassmann).",
        "",
        "## Provenance",
        "",
        f"- N_surrogates = {n_surr}; n_swaps = {swap_factor}·N(N−1)/2; "
        f"seed = {SEED}",
        "- Cohort: " + ", ".join(COHORT),
        "- Bands: " + ", ".join(BANDS),
        "- Phases: rest_pre, task_test, rest_post (FULL phases)",
        "- FC method: imcoh_abs",
        "- K grid: " + ", ".join(str(k) for k in K_GRID),
        f"- Wall-clock runtime: {runtime_s:.1f} s",
        "- Build: scripts/01_compute/audit/audit_141_spectral_clustering_headtohead.py",
        "",
        "## Caveats (referee-facing)",
        "",
        "- The surrogate is INDEPENDENT per phase (shared with audit_63, so "
        "the comparison is fair) — it controls strength heterogeneity but "
        "does NOT isolate edge-identity cross-phase structure for ANY method.",
        "- ARI/NMI are quantized; surrogate upper-tail p can be coarse.",
        "- k-means random init fixed (seed, n_init=4; verdict-stable vs "
        "n_init=10); the baseline is given "
        "L_sym AND L_unnorm, ARI AND NMI, the full K sweep, and the "
        "best-case-for-baseline K — a genuinely strong shot, not a strawman "
        "knee-capped by one bad K.",
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
    ap.add_argument("--out-suffix", default="",
                    help="suffix for the per-patient CSV so parallel "
                         "band-streams don't clobber each other; merged "
                         "afterwards with --merge.")
    ap.add_argument("--merge", action="store_true",
                    help="skip compute; merge all per_patient_per_band*.csv "
                         "shards into the final CSVs + README.")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if args.merge:
        shards = sorted(OUT.glob("per_patient_per_band_*.csv"))
        if not shards:
            print("[audit_141] no shards to merge")
            return
        per_pat = pd.concat([pd.read_csv(s) for s in shards],
                            ignore_index=True)
        per_pat.to_csv(OUT / "per_patient_per_band.csv", index=False)
        cohort = cohort_summary(per_pat)
        cohort.to_csv(OUT / "cohort_summary.csv", index=False)
        print("\n" + cohort.to_string(index=False))
        write_readme(cohort, args.n_surrogates, args.swap_factor, 0.0)
        print(f"\n[audit_141] merged {len(shards)} shards; outputs at {OUT}")
        return

    t0 = time.time()
    all_rows: list[dict] = []
    for band in args.bands:
        for pat in args.patients:
            t_cell = time.time()
            rows = per_cell(pat, band, args.n_surrogates, args.swap_factor,
                            args.seed, args.verbose)
            all_rows.extend(rows)
            if rows:
                print(f"[audit_141] {pat}/{band}: {len(rows)} method-cells "
                      f"({time.time() - t_cell:.1f}s)")

    runtime = time.time() - t0
    per_pat = pd.DataFrame(all_rows)
    suffix = f"_{args.out_suffix}" if args.out_suffix else ""
    per_pat.to_csv(OUT / f"per_patient_per_band{suffix}.csv", index=False)

    # When sharding (out-suffix set), defer cohort summary + README to the
    # --merge pass so all bands are present.
    if args.out_suffix:
        print(f"\n[audit_141] shard '{args.out_suffix}' done in {runtime:.1f}s "
              f"→ per_patient_per_band{suffix}.csv")
        return

    cohort = cohort_summary(per_pat)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)
    print("\n" + cohort.to_string(index=False))

    write_readme(cohort, args.n_surrogates, args.swap_factor, runtime)
    print(f"\n[audit_141] total {runtime:.1f}s; outputs at {OUT}")


if __name__ == "__main__":
    main()
