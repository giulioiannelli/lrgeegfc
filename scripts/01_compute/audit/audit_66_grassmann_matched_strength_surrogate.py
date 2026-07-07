#!/usr/bin/env python3
"""Audit 66 — Grassmann matched-strength surrogate null (§5.4 spectral settle).

Companion to audit_63 (ρ_split) and audit_65 (KC tree distance). The §5.4
Grassmann chordal-distance result currently rests on a within-session null.
audit_63 returned cohort-significant β / α at the per-pair ρ_split layer.
audit_65 returned no separated cell at the KC tree-distance layer (β λ=0
contaminated; α λ=0 also_negative). This script puts the §5.4 Grassmann
subspace under the same R=200 strength-preserving 4-cycle ±δ surrogate
across the full §5.4 k-grid (k = 2..80).

For each (patient, band ∈ {α, β, γ_l}) cell, generate R=200 surrogate
Laplacians per phase Φ ∈ {rsPre_A, taskTest, rsPost} via the same
algorithm as audit_63 / audit_65 (n_swaps = SWAP_FACTOR · N(N−1)/2,
SWAP_FACTOR = 20). For each surrogate r and each k ∈ K_GRID:

    T_G_surr_r(k) = d_chord(k; U^tt_surr, U^post_surr)
                  - d_chord(k; U^pre_A_surr, U^tt_surr)

with d_chord(k; A, B) = sqrt(k − sum_i σ_i²) where σ_i are the singular
values of A^T B (A, B ∈ R^{N×k}). Eigenvectors are the slowest k non-trivial
modes (skipping the zero-mode at index 0), matching audit_37 / audit_46.

Outputs
-------
    data/audit/grassmann_matched_strength_surrogate/
        cohort_summary.csv
        per_patient_per_band_per_k.csv
        joint_signature.csv
        figures/grassmann_cohort_distribution.pdf
        figures/grassmann_per_patient_panel.pdf
        figures/grassmann_kspan_verdict.pdf
        README.md

Surrogate cache (introduced here for reuse)
-------------------------------------------
This is the first script that caches the surrogate Laplacian
eigendecompositions to disk so any future probe can compute its statistic
on the SAME ensemble without regenerating the surrogates. Cache layout:

    data/cache/matched_strength_surrogate_lrg/Pat_NN/
        {band}_{phase}_R{R}_swap{SF}_seed{S}_imcoh_abs.npz
        # contains: eigvals (R, N) float64, eigvecs (R, N, N) float64

New audits SHOULD import the helpers (now extracted to a library module):

    from lrg_eegfc.utils.surrogate.matched_strength import (
        load_or_compute_surrogate_eigs,
        strength_preserving_shuffle,
    )

The local copies of `strength_preserving_shuffle`, `verify_strengths`,
`load_or_compute_surrogate_eigs`, `surrogate_cache_path` in this script
are kept (rather than swapped for library imports) only to preserve
exact byte-for-byte reproducibility of the 2026-05-11 first-pass run.
A subsequent refactor will swap them for library imports.

Audits 62 / 63 / 65 predate this cache; their docstrings now point here.
The canonical R=200, SWAP_FACTOR=20, seed=20260511 ensemble for trace
bands {α, β, γ_l} × phases {rest_pre_A, task_test, rest_post} is the
cache produced by THIS run.
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BAND_TEX_DICT, FS_OVERRIDES, nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.workflow.fc import load_fc_matrix

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # type: ignore


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["alpha", "beta", "low_gamma"]
# Canonical EEG band order for figure / report iteration (delta..high_gamma).
CANONICAL_BAND_ORDER = ["delta", "theta", "alpha", "beta",
                         "low_gamma", "high_gamma"]
PHASES = ("rest_pre_A", "task_test", "rest_post")


def _band_order_present(per_pat: pd.DataFrame) -> list[str]:
    """Bands appearing in `per_pat` ordered by canonical EEG order."""
    present = set(per_pat.band.unique())
    return [b for b in CANONICAL_BAND_ORDER if b in present]
# Full spectrum: min(N) across cohort = 113 (Pat_10) → max k = 112.
K_GRID = list(range(2, 113))
K_REFERENCE = (3, 20, 60, 100)
N_SURROGATES = 200
SWAP_FACTOR = 20

OUT = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate"
FIG = OUT / "figures"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"
SURROGATE_LRG_CACHE = CACHE_ROOT / "matched_strength_surrogate_lrg"
AUDIT_63_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
AUDIT_65_DIR = ROOT / "data" / "audit" / "kc_matched_strength_surrogate"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
HALVES_FC_CACHE.mkdir(parents=True, exist_ok=True)
SURROGATE_LRG_CACHE.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Half-FC cache (lazy populate, identical to audit_63 / audit_65)
# ---------------------------------------------------------------------------
def _half_fc_path(pat: str, band: str, half: str) -> Path:
    return (HALVES_FC_CACHE / pat
            / f"{band}_rest_pre_{half}_imcoh_abs.npy")


def ensure_half_fcs(pat: str, bands: list[str]) -> None:
    missing = [(b, h) for b in bands for h in ("A", "B")
               if not _half_fc_path(pat, b, h).exists()]
    if not missing:
        return
    print(f"[audit_66] {pat}: caching {len(missing)} missing half FCs")
    X = load_timeseries(pat, "rest_pre", SEEG_DATAPATH)
    X = np.asarray(X, dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nperseg_half = max(256, nperseg_for_fs(fs) // 2)
    halves_fc = compute_imcoh_abs_halves(X, fs, nperseg_half, BRAIN_BANDS)
    del X
    gc.collect()
    pat_dir = HALVES_FC_CACHE / pat
    pat_dir.mkdir(parents=True, exist_ok=True)
    for (band, half), A in halves_fc.items():
        if band not in bands:
            continue
        path = _half_fc_path(pat, band, half)
        np.save(path, np.asarray(A, dtype=np.float32))


def load_phase_fc(pat: str, phase: str, band: str) -> np.ndarray:
    if phase in ("rest_pre_A", "rest_pre_B"):
        half = phase[-1]
        W = np.load(_half_fc_path(pat, band, half))
    else:
        W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    W = 0.5 * (W + W.T)
    return W


# ---------------------------------------------------------------------------
# Laplacian eigendecomposition (full; we cache full eigvals + eigvecs so
# downstream probes — Grassmann, KC, ρ-propagator, etc. — can reuse the
# same surrogate ensemble without re-running the shuffle + eigendecomp.)
# ---------------------------------------------------------------------------
def laplacian_eig(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    eigvals, eigvecs = np.linalg.eigh(L)
    return eigvals, eigvecs


def topk_basis(eigvecs: np.ndarray, k_max: int) -> np.ndarray:
    """Return eigvecs[:, 1:k_max+1] (slowest k_max non-trivial modes)."""
    return np.ascontiguousarray(eigvecs[:, 1:k_max + 1])


# ---------------------------------------------------------------------------
# Strength-preserving 4-cycle ±δ rewiring (identical to audit_63 / audit_65)
# ---------------------------------------------------------------------------
def strength_preserving_shuffle(W: np.ndarray, n_swaps: int,
                                rng: np.random.Generator,
                                w_max: float = 1.0) -> np.ndarray:
    W = W.copy()
    N = W.shape[0]
    samples = rng.integers(0, N, size=(n_swaps, 4))
    fracs = rng.uniform(0.0, 1.0, size=n_swaps)
    for i in range(n_swaps):
        a, b, c, d = samples[i]
        if a == b or a == c or a == d or b == c or b == d or c == d:
            continue
        w1 = W[a, b]; w2 = W[c, d]; w3 = W[a, d]; w4 = W[c, b]
        lo = -w1 if -w1 > -w2 else -w2
        if w3 - w_max > lo:
            lo = w3 - w_max
        if w4 - w_max > lo:
            lo = w4 - w_max
        hi = w_max - w1 if w_max - w1 < w_max - w2 else w_max - w2
        if w3 < hi:
            hi = w3
        if w4 < hi:
            hi = w4
        if lo >= hi:
            continue
        delta = lo + fracs[i] * (hi - lo)
        n1 = w1 + delta; n2 = w2 + delta
        n3 = w3 - delta; n4 = w4 - delta
        W[a, b] = n1; W[b, a] = n1
        W[c, d] = n2; W[d, c] = n2
        W[a, d] = n3; W[d, a] = n3
        W[c, b] = n4; W[b, c] = n4
    return W


def verify_strengths(W_obs: np.ndarray, W_surr: np.ndarray,
                      tol: float = 1e-4) -> bool:
    s_obs = W_obs.sum(axis=1)
    s_surr = W_surr.sum(axis=1)
    return bool(np.max(np.abs(s_obs - s_surr)) < tol)


# ---------------------------------------------------------------------------
# Grassmann chordal distance d(A_k, B_k) = sqrt(k - sum sigma_i^2)
# ---------------------------------------------------------------------------
def chordal(V_a: np.ndarray, V_b: np.ndarray, k: int) -> float:
    A = V_a[:, :k]
    B = V_b[:, :k]
    M = A.T @ B
    sigma = np.linalg.svd(M, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    return float(np.sqrt(max(k - float((sigma ** 2).sum()), 0.0)))


def t_g_at_all_k(U_pre: np.ndarray, U_tt: np.ndarray,
                  U_post: np.ndarray, k_grid: list[int]) -> np.ndarray:
    """Return T_G[i] for each k in k_grid.

    T_G = d(rest_pre, task) - d(task, rest_post). Project-wide convention:
    T_G > 0 = trace (rsPost closer to task than rsPre is to task).
    """
    out = np.empty(len(k_grid), dtype=float)
    for i, k in enumerate(k_grid):
        d_pre_tt = chordal(U_pre, U_tt, k)
        d_tt_post = chordal(U_tt, U_post, k)
        out[i] = d_pre_tt - d_tt_post
    return out


# ---------------------------------------------------------------------------
# Surrogate eigendecomposition cache
# ---------------------------------------------------------------------------
def _surr_cache_path(pat: str, band: str, phase: str,
                      n_surr: int, swap_factor: int, seed: int) -> Path:
    pat_dir = SURROGATE_LRG_CACHE / pat
    return (pat_dir / f"{band}_{phase}_R{n_surr}_swap{swap_factor}"
                       f"_seed{seed}_imcoh_abs.npz")


def load_or_compute_surrogate_eigs(pat: str, band: str, phase: str,
                                    W: np.ndarray, n_surr: int,
                                    swap_factor: int, seed: int,
                                    rng: np.random.Generator,
                                    verbose: bool = False
                                    ) -> tuple[np.ndarray, np.ndarray]:
    """Load cached (eigvals, eigvecs) of shape (R, N), (R, N, N) if present.

    Otherwise compute R surrogates of W via the strength-preserving shuffle,
    eigendecompose each, save to cache, and return.
    """
    path = _surr_cache_path(pat, band, phase, n_surr, swap_factor, seed)
    if path.exists():
        with np.load(path) as data:
            evals = data["eigvals"]
            evecs = data["eigvecs"]
        if (evals.shape[0] == n_surr and evecs.shape[0] == n_surr
                and evals.shape[1] == W.shape[0]
                and evecs.shape[1] == W.shape[0]
                and evecs.shape[2] == W.shape[0]):
            return evals.astype(np.float64), evecs.astype(np.float64)
        if verbose:
            print(f"[audit_66] cache shape mismatch at {path}; recomputing")

    N = W.shape[0]
    n_swaps = swap_factor * (N * (N - 1)) // 2
    evals = np.empty((n_surr, N), dtype=np.float64)
    evecs = np.empty((n_surr, N, N), dtype=np.float64)
    ok_count = 0
    for r in range(n_surr):
        W_s = strength_preserving_shuffle(W, n_swaps, rng)
        if not verify_strengths(W, W_s, tol=1e-4):
            evals[r, :] = np.nan
            evecs[r, :, :] = np.nan
            continue
        ev, ec = laplacian_eig(W_s)
        evals[r, :] = ev
        evecs[r, :, :] = ec
        ok_count += 1

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, eigvals=evals, eigvecs=evecs)
    if verbose:
        print(f"[audit_66] cached {ok_count}/{n_surr} surrogates → {path.name}")
    return evals, evecs


# ---------------------------------------------------------------------------
# Per-cell pipeline
# ---------------------------------------------------------------------------
def per_cell(pat: str, band: str, n_surr: int, swap_factor: int,
             k_grid: list[int], seed: int, rng: np.random.Generator,
             verbose: bool = False) -> dict | None:
    Ws: dict[str, np.ndarray] = {}
    for phase in PHASES:
        try:
            Ws[phase] = load_phase_fc(pat, phase, band)
        except Exception as e:
            if verbose:
                print(f"[audit_66] SKIP {pat}/{band}/{phase}: {e}")
            return None
    N = Ws[PHASES[0]].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        return None
    k_max = max(k_grid)
    if k_max + 1 > N:
        if verbose:
            print(f"[audit_66] SKIP {pat}/{band}: N={N} < k_max+1={k_max+1}")
        return None

    # Observed (split-baseline; pre_A as rsPre leg)
    U_obs = {phase: topk_basis(laplacian_eig(Ws[phase])[1], k_max)
              for phase in PHASES}
    obs_T_G = t_g_at_all_k(U_obs["rest_pre_A"], U_obs["task_test"],
                            U_obs["rest_post"], k_grid)

    # Surrogates — load from cache or compute + cache
    n_swaps = swap_factor * (N * (N - 1)) // 2
    surr_eigvecs: dict[str, np.ndarray] = {}
    for phase in PHASES:
        _, evecs = load_or_compute_surrogate_eigs(
            pat, band, phase, Ws[phase],
            n_surr, swap_factor, seed, rng, verbose=verbose)
        surr_eigvecs[phase] = evecs

    surr_T_G = np.empty((n_surr, len(k_grid)), dtype=float)
    for r in range(n_surr):
        U_surr = {phase: topk_basis(surr_eigvecs[phase][r], k_max)
                  for phase in PHASES}
        if any(np.isnan(U_surr[phase]).any() for phase in PHASES):
            surr_T_G[r, :] = np.nan
            continue
        surr_T_G[r, :] = t_g_at_all_k(
            U_surr["rest_pre_A"], U_surr["task_test"],
            U_surr["rest_post"], k_grid)

    # Per-k stats
    rows: list[dict] = []
    for i, k in enumerate(k_grid):
        s = surr_T_G[:, i]
        s_finite = s[np.isfinite(s)]
        if s_finite.size == 0:
            continue
        s_mean = float(np.mean(s_finite))
        s_std = float(np.std(s_finite, ddof=1))
        obs_T = float(obs_T_G[i])
        z = (obs_T - s_mean) / s_std if s_std > 0 else float("nan")
        # Trace direction: T_G > 0; one-sided upper-tail.
        p_upper = float(np.mean(s_finite >= obs_T))
        rows.append({
            "patient": pat,
            "band": band,
            "k": int(k),
            "N_nodes": int(N),
            "n_swaps_per_surrogate": int(n_swaps),
            "n_surrogates": int(s_finite.size),
            "obs_T_G": obs_T,
            "surr_T_G_mean": s_mean,
            "surr_T_G_std": s_std,
            "surr_T_G_p5": float(np.quantile(s_finite, 0.05)),
            "surr_T_G_p25": float(np.quantile(s_finite, 0.25)),
            "surr_T_G_p50": float(np.quantile(s_finite, 0.50)),
            "surr_T_G_p75": float(np.quantile(s_finite, 0.75)),
            "surr_T_G_p95": float(np.quantile(s_finite, 0.95)),
            "obs_z": z,
            "obs_p_one_sided_upper": p_upper,
        })
    return {
        "rows": rows,
        "obs_T_G": obs_T_G,
        "surr_T_G": surr_T_G,
        "patient": pat,
        "band": band,
    }


# ---------------------------------------------------------------------------
# Cohort summary per (band, k)
# ---------------------------------------------------------------------------
def cohort_summary(per_pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in _band_order_present(per_pat):
        for k in K_GRID:
            sub = per_pat[(per_pat.band == band) & (per_pat.k == k)]
            if sub.empty:
                continue
            obs_T = sub.obs_T_G.values
            surr_med = sub.surr_T_G_p50.values
            # obs_p_one_sided_upper was the lower-tail p; under the flipped
            # T_G>0=trace convention, the "trace direction" tail is upper.
            n_above = int((sub.obs_p_one_sided_upper < 0.05).sum())
            try:
                wz, wp = wilcoxon(obs_T - surr_med, alternative="greater")
                cohort_z = float(wz)
                cohort_p = float(wp)
            except Exception:
                cohort_z = float("nan")
                cohort_p = float("nan")
            med_obs = float(np.median(obs_T))
            med_surr = float(np.median(surr_med))

            # Verdict — same labels as audit_65 (trace direction = positive)
            if (cohort_p < 0.05 and abs(med_surr) < 0.05 * max(1.0, abs(med_obs))
                    and n_above >= 8):
                verdict = "separated"
            elif med_obs > 0 and med_surr > 0.5 * med_obs:
                verdict = "also_positive"
            else:
                verdict = "intermediate"
            out.append({
                "band": band,
                "k": k,
                "n_patients": len(sub),
                "obs_median_T_G": med_obs,
                "surr_median_T_G_per_patient_median": med_surr,
                "n_patients_above_own_surrogate": n_above,
                "n_patients_above_own_surrogate_str": f"{n_above}/{len(sub)}",
                "paired_wilcoxon_z": cohort_z,
                "paired_wilcoxon_p": cohort_p,
                "verdict": verdict,
            })
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Joint signature with audit_63 + audit_65
# ---------------------------------------------------------------------------
def joint_signature(per_pat: pd.DataFrame) -> pd.DataFrame:
    a63 = (pd.read_csv(AUDIT_63_DIR / "per_patient_per_band.csv")
           if (AUDIT_63_DIR / "per_patient_per_band.csv").exists()
           else pd.DataFrame())
    a65 = (pd.read_csv(AUDIT_65_DIR / "per_patient_per_band.csv")
           if (AUDIT_65_DIR / "per_patient_per_band.csv").exists()
           else pd.DataFrame())

    rows = []
    for pat in COHORT:
        for band in _band_order_present(per_pat):
            r63 = a63[(a63.patient == pat) & (a63.band == band)] if not a63.empty else pd.DataFrame()
            l0 = (a65[(a65.patient == pat) & (a65.band == band)
                       & (a65["lambda"] == 0.0)]
                  if not a65.empty else pd.DataFrame())
            l1 = (a65[(a65.patient == pat) & (a65.band == band)
                       & (a65["lambda"] == 1.0)]
                  if not a65.empty else pd.DataFrame())

            rho_z = float(r63.iloc[0]["obs_z"]) if not r63.empty else float("nan")
            rho_p = float(r63.iloc[0]["obs_p_one_sided"]) if (
                not r63.empty and "obs_p_one_sided" in r63.columns
            ) else float("nan")
            kc0_z = float(l0.iloc[0]["obs_z"]) if not l0.empty else float("nan")
            kc0_p = float(l0.iloc[0]["obs_p_one_sided_upper"]) if not l0.empty else float("nan")
            kc1_z = float(l1.iloc[0]["obs_z"]) if not l1.empty else float("nan")
            kc1_p = float(l1.iloc[0]["obs_p_one_sided_upper"]) if not l1.empty else float("nan")

            grass_zs = {}
            grass_ps = {}
            for k_ref in K_REFERENCE:
                cell = per_pat[(per_pat.patient == pat)
                               & (per_pat.band == band)
                               & (per_pat.k == k_ref)]
                if cell.empty:
                    grass_zs[k_ref] = float("nan")
                    grass_ps[k_ref] = float("nan")
                else:
                    grass_zs[k_ref] = float(cell.iloc[0]["obs_z"])
                    grass_ps[k_ref] = float(cell.iloc[0]["obs_p_one_sided_upper"])

            n_sig = 0
            if np.isfinite(rho_p) and rho_p < 0.05:
                n_sig += 1
            if np.isfinite(kc0_p) and kc0_p < 0.05:
                n_sig += 1
            if np.isfinite(kc1_p) and kc1_p < 0.05:
                n_sig += 1
            for k_ref in K_REFERENCE:
                if np.isfinite(grass_ps[k_ref]) and grass_ps[k_ref] < 0.05:
                    n_sig += 1

            rows.append({
                "patient": pat,
                "band": band,
                "rho_split_z": rho_z,
                "rho_split_p_one_sided_upper": rho_p,
                "T_KC_l0_z": kc0_z,
                "T_KC_l0_p_one_sided_lower": kc0_p,
                "T_KC_l1_z": kc1_z,
                "T_KC_l1_p_one_sided_lower": kc1_p,
                "T_G_k3_z": grass_zs[3],
                "T_G_k3_p_one_sided_lower": grass_ps[3],
                "T_G_k20_z": grass_zs[20],
                "T_G_k20_p_one_sided_lower": grass_ps[20],
                "T_G_k60_z": grass_zs[60],
                "T_G_k60_p_one_sided_lower": grass_ps[60],
                "n_probes_significant_at_p05": n_sig,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def make_cohort_figure(per_pat: pd.DataFrame, cohort: pd.DataFrame,
                        out_path: Path) -> None:
    bands = _band_order_present(per_pat)
    n_b = len(bands)
    fig, axes = plt.subplots(n_b, 1, figsize=(7.5, 2.6 * n_b), sharex=True)
    if n_b == 1:
        axes = [axes]
    for ax, band in zip(axes, bands):
        sub_co = cohort[cohort.band == band].sort_values("k")
        if sub_co.empty:
            ax.set_axis_off()
            continue
        ks = sub_co.k.values
        obs_med = sub_co.obs_median_T_G.values
        surr_med = sub_co.surr_median_T_G_per_patient_median.values
        # Per-band per-k cohort 5-95% across patients
        p5 = []
        p95 = []
        for k in ks:
            cells = per_pat[(per_pat.band == band) & (per_pat.k == k)]
            obs_vals = cells.obs_T_G.values
            p5.append(np.quantile(obs_vals, 0.05) if obs_vals.size else np.nan)
            p95.append(np.quantile(obs_vals, 0.95) if obs_vals.size else np.nan)
        ax.fill_between(ks, p5, p95, color="#1f3d6e", alpha=0.18,
                        label="obs 5–95%")
        ax.plot(ks, obs_med, color="#1f3d6e", lw=1.4, label="obs median")
        ax.plot(ks, surr_med, color="#aa4444", lw=1.2,
                label="surr median (per-pat med)")
        ax.axhline(0, color="0.6", lw=0.6, ls="--")
        ax.set_ylabel(rf"$T_G$ — {BRAIN_BAND_TEX_DICT[band]}")
        ax.spines[["top", "right"]].set_visible(False)
        if band == bands[0]:
            ax.legend(frameon=False, fontsize=8, loc="lower right")
    axes[-1].set_xlabel(r"spectral cutoff $k$")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_per_patient_panel(per_pat: pd.DataFrame, out_path: Path) -> None:
    bands = _band_order_present(per_pat)
    n_pat = len(COHORT)
    n_b = len(bands)
    n_k = len(K_REFERENCE)
    fig, axes = plt.subplots(n_pat, n_b * n_k,
                             figsize=(2.6 * n_b * n_k, 1.55 * n_pat),
                             sharex=False)
    for i, pat in enumerate(COHORT):
        for j, band in enumerate(bands):
            for kk, k_ref in enumerate(K_REFERENCE):
                col = j * n_k + kk
                ax = axes[i, col]
                cell = per_pat[(per_pat.patient == pat)
                               & (per_pat.band == band)
                               & (per_pat.k == k_ref)]
                if cell.empty:
                    ax.set_axis_off()
                    continue
                r = cell.iloc[0]
                # Show per-cell p5 / p50 / p95 plus obs
                ax.bar(["p5", "p50", "p95"],
                       [r.surr_T_G_p5, r.surr_T_G_p50, r.surr_T_G_p95],
                       color="#cccccc", edgecolor="#888888")
                ax.axhline(r.obs_T_G, color="#1f3d6e", lw=1.6,
                           label=f"obs={r.obs_T_G:+.2f}")
                ax.axhline(0, color="#888888", lw=0.6, ls="--")
                if i == 0:
                    ax.set_title(
                        f"{BRAIN_BAND_TEX_DICT[band]} k={k_ref}",
                        fontsize=9)
                if col == 0:
                    ax.set_ylabel(pat.replace("Pat_", "P"), fontsize=8)
                ax.tick_params(labelsize=7)
                ax.legend(loc="upper left", frameon=False, fontsize=6)
                ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_kspan_verdict_figure(cohort: pd.DataFrame, out_path: Path) -> None:
    bands = [b for b in CANONICAL_BAND_ORDER if b in set(cohort.band.unique())]
    n_b = len(bands)
    fig, axes = plt.subplots(n_b, 1, figsize=(8.5, 1.8 * n_b), sharex=True)
    if n_b == 1:
        axes = [axes]
    color_map = {"separated": "#1f7a1f", "intermediate": "#888888",
                 "also_positive": "#aa4444"}
    for ax, band in zip(axes, bands):
        sub = cohort[cohort.band == band].sort_values("k")
        if sub.empty:
            ax.set_axis_off()
            continue
        for _, r in sub.iterrows():
            ax.bar(int(r.k), 1, width=1.0,
                   color=color_map.get(str(r.verdict), "#cccccc"),
                   edgecolor="none")
        # Annotate cohort p curve
        ax2 = ax.twinx()
        ks = sub.k.values
        ps = sub.paired_wilcoxon_p.values
        ax2.plot(ks, ps, color="#1f3d6e", lw=1.2)
        ax2.axhline(0.05, color="#aa4444", lw=0.5, ls="--")
        ax2.set_ylabel(r"cohort Wilcoxon $p$", fontsize=8)
        ax2.set_ylim(0, 1.0)
        ax.set_yticks([])
        ax.set_ylabel(BRAIN_BAND_TEX_DICT[band], fontsize=10)
        ax.spines[["top", "left", "right"]].set_visible(False)
    axes[-1].set_xlabel(r"spectral cutoff $k$")
    # Legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, label=l)
               for l, c in color_map.items()]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(cohort: pd.DataFrame, per_pat: pd.DataFrame,
                 n_surr: int, swap_factor: int,
                 runtime_s: float) -> None:
    lines = [
        "---",
        "name: grassmann_matched_strength_surrogate",
        "scope: section_5_4_grassmann_robustness_against_matched_strength_null",
        "date: 2026-05-11",
        "status: first_pass",
        "---",
        "",
        "# Grassmann subspace under matched-strength surrogacy",
        "",
        f"**Head.** R={n_surr} strength-preserving (4-cycle ±δ, "
        f"SWAP_FACTOR={swap_factor}) surrogate Laplacians per (patient, "
        f"band ∈ {{{', '.join(_band_order_present(per_pat))}}}, "
        f"phase ∈ {{rsPre_A, taskTest, rsPost}}). "
        f"For each surrogate, the Grassmann chordal distance T_G(k) is "
        f"computed at every k ∈ K_GRID = "
        f"{K_GRID[0]}..{K_GRID[-1]} across the §5.4 spectral continuum. "
        f"Observed counterpart uses U^pre_A as the rsPre leg to match "
        f"the surrogate baseline.",
        "",
        "## Per-band cohort surface",
        "",
        "| band | k_min | k_max | k_separated | k_int_p_lt_05 | k_also_pos | n_k_sep | longest_run_sep |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for band in _band_order_present(per_pat):
        sub = cohort[cohort.band == band]
        if sub.empty:
            continue
        ks = sub.k.values
        verdicts = sub.verdict.values
        ps = sub.paired_wilcoxon_p.values
        sep_mask = (verdicts == "separated")
        int_p_mask = (verdicts == "intermediate") & (ps < 0.05)
        pos_mask = (verdicts == "also_positive")
        n_sep = int(sep_mask.sum())
        # Longest contiguous separated run
        longest = 0
        cur = 0
        for v in sep_mask:
            cur = cur + 1 if v else 0
            if cur > longest:
                longest = cur
        k_sep_str = ",".join(str(int(k)) for k, v in zip(ks, sep_mask) if v)
        k_int_p_str = ",".join(str(int(k)) for k, v in zip(ks, int_p_mask) if v)
        k_pos_str = ",".join(str(int(k)) for k, v in zip(ks, pos_mask) if v)
        lines.append(
            f"| {band} | {int(ks.min())} | {int(ks.max())} "
            f"| {k_sep_str or '—'} "
            f"| {k_int_p_str or '—'} "
            f"| {k_pos_str or '—'} "
            f"| {n_sep}/{len(ks)} | {longest} |"
        )

    lines.extend([
        "",
        "## Reference cutoffs (k=3, 20, 60) per band",
        "",
        "| band | k | obs median T_G | surr median (per-pat med) | n above | Wilcoxon p | verdict |",
        "|---|---|---|---|---|---|---|",
    ])
    for band in _band_order_present(per_pat):
        for k_ref in K_REFERENCE:
            row = cohort[(cohort.band == band) & (cohort.k == k_ref)]
            if row.empty:
                continue
            r = row.iloc[0]
            lines.append(
                f"| {band} | {int(r.k)} "
                f"| {r['obs_median_T_G']:+.3f} "
                f"| {r['surr_median_T_G_per_patient_median']:+.3f} "
                f"| {r['n_patients_above_own_surrogate_str']} "
                f"| {r['paired_wilcoxon_p']:.4f} "
                f"| **{r['verdict']}** |"
            )

    lines.extend([
        "",
        "## Verdict labels",
        "",
        "- **separated**: cohort Wilcoxon p < 0.05 (one-sided, T_obs > T_surr; "
        "trace direction T_G > 0 per feedback_td_sign_convention.md) "
        "AND |median surr T_G across patients| < 0.05 × |median obs T_G| AND "
        "≥ 8/10 patients individually above their own surrogate at p<0.05.",
        "- **also_positive**: median surr T_G > 0.5 × median obs T_G AND "
        "median obs T_G > 0 (surrogate substantially recovers trace direction).",
        "- **intermediate**: anything else.",
        "",
        "## Test choice",
        "",
        "Algorithm matches audit_63 (ρ_split) and audit_65 (KC) exactly. "
        "Independent per-phase 4-cycle ±δ rewiring; SWAP_FACTOR = 20; "
        "n_swaps = SWAP_FACTOR × N(N−1)/2 per surrogate per phase. "
        "Observed split-baseline uses U^pre_A as the rsPre leg.",
        "",
        "Eigenvectors of the symmetric Laplacian L = D − W (zero-mode "
        "skipped at index 0; slowest k_max = 80 non-trivial modes "
        "retained). Chordal distance d(A_k, B_k) = sqrt(k − Σ σ_i²) "
        "where σ_i = svdvals(A^T B). Equivalent to sqrt(Σ sin²θ_i) for "
        "principal angles θ_i; standard Grassmann chordal distance.",
        "",
        "## Provenance",
        "",
        f"- N_surrogates = {n_surr} per cell",
        f"- n_swaps_per_surrogate = SWAP_FACTOR ({swap_factor}) × N(N−1)/2",
        "- Cohort: " + ", ".join(COHORT),
        "- Bands: " + ", ".join(_band_order_present(per_pat)),
        "- Phases: " + ", ".join(PHASES),
        f"- K_GRID = {K_GRID[0]}..{K_GRID[-1]} (n={len(K_GRID)})",
        "- FC method: imcoh_abs",
        "- Half-FC cache reused: `data/cache/imcoh_halves_fc/Pat_NN/`",
        "- audit_63 / audit_65 z-scores joined into `joint_signature.csv`",
        f"- Wall-clock runtime: {runtime_s:.1f} s",
        "- Build script: `scripts/01_compute/audit/audit_66_grassmann_matched_strength_surrogate.py`",
        "",
        "## Files",
        "",
        "- `per_patient_per_band_per_k.csv` — observed + per-cell surrogate stats (one row per patient × band × k)",
        "- `cohort_summary.csv` — per-(band, k) cohort verdict",
        "- `joint_signature.csv` — extended four-probe table (ρ_split + KC λ=0/1 + Grass k=3/20/60)",
        "- `figures/grassmann_cohort_distribution.pdf` — per-band cohort surface across k",
        "- `figures/grassmann_per_patient_panel.pdf` — per-(patient, band, k_ref) cell summary",
        "- `figures/grassmann_kspan_verdict.pdf` — per-band verdict color strip across k",
        "",
    ])
    (OUT / "README.md").write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=TARGET_BANDS)
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--n-surrogates", type=int, default=N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=SWAP_FACTOR)
    ap.add_argument("--seed", type=int, default=20260511)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    bands = args.bands
    patients = args.patients
    n_surr = args.n_surrogates
    swap_factor = args.swap_factor

    print("[audit_66] pre-flight: ensure half FCs cached")
    for pat in patients:
        ensure_half_fcs(pat, bands)

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    all_rows: list[dict] = []
    for band in bands:
        for pat in patients:
            t_cell = time.time()
            cell = per_cell(pat, band, n_surr, swap_factor, K_GRID,
                             args.seed, rng, verbose=args.verbose)
            if cell is None:
                continue
            all_rows.extend(cell["rows"])
            dt = time.time() - t_cell
            # Quick per-cell summary at three reference k
            summary_bits = []
            for k_ref in K_REFERENCE:
                row = next((r for r in cell["rows"] if r["k"] == k_ref), None)
                if row is not None:
                    summary_bits.append(
                        f"k{k_ref}:T={row['obs_T_G']:+.2f} z={row['obs_z']:+.2f}"
                    )
            print(f"[audit_66] {pat}/{band}: {' | '.join(summary_bits)} ({dt:.1f}s)")

    runtime = time.time() - t0
    print(f"[audit_66] total: {runtime:.1f}s for {len({(r['patient'], r['band']) for r in all_rows})} cells")

    per_pat = pd.DataFrame(all_rows)
    per_pat.to_csv(OUT / "per_patient_per_band_per_k.csv", index=False)

    cohort = cohort_summary(per_pat)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)

    joint = joint_signature(per_pat)
    joint.to_csv(OUT / "joint_signature.csv", index=False)

    make_cohort_figure(per_pat, cohort, FIG / "grassmann_cohort_distribution.pdf")
    make_per_patient_panel(per_pat, FIG / "grassmann_per_patient_panel.pdf")
    make_kspan_verdict_figure(cohort, FIG / "grassmann_kspan_verdict.pdf")

    write_readme(cohort, per_pat, n_surr, swap_factor, runtime)
    print(f"[audit_66] outputs at {OUT}")


if __name__ == "__main__":
    main()
