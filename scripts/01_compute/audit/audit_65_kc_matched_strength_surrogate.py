#!/usr/bin/env python3
"""Audit 64 — KC tree-distance matched-strength Laplacian surrogate null.

Companion to audit_63. The §5.2 β KC 10/10 cohort claim (q=0.006 at m=12,
ten patients individually below their own within-baseline split-half null
at both λ=0 and λ=1) currently rests on the *narrow* within-baseline
null. audit_63 returned an *intermediate* verdict for the §5.3 ρ_split
test against a wider matched-strength null (β cohort z=+0.48, n_above
7/10). This script puts the §5.2 KC tree distance under the same wider
null so the two probes are evaluated against the same surrogate family.

For each (patient, band ∈ {α, β, γ_l}) cell, generate R=200
surrogate Laplacians per phase Φ ∈ {rsPre_A, taskTest, rsPost} via
the **same 4-cycle ±δ strength-preserving rewiring as audit_63**
(n_swaps = SWAP_FACTOR * N(N-1)/2, SWAP_FACTOR=20). For each surrogate
r and each λ ∈ {0, 1}:

    T_KC_surr_r(λ) = d_KC(λ; D^tt_surr_r, D^post_surr_r)
                    - d_KC(λ; D^pre_A_surr_r, D^tt_surr_r)

The observed counterpart uses D^pre_A as the rsPre leg to match the
surrogate baseline (NOT the §5.2 full rsPre leg — see README for the
implication).

Cache reuse
-----------
Half FC matrices (rsPre_A, rsPre_B) are reused from audit_63's lazy
cache at `data/cache/imcoh_halves_fc/Pat_NN/`. Surrogate adjacency
matrices themselves are NOT cached by audit_63 (they were streamed),
so this script regenerates them by the same algorithm with a fresh
deterministic seed. Exact rng-state matching with audit_63 is not
possible; surrogate-ensemble distributions are matched in the
statistical sense only.

Outputs
-------
    data/audit/kc_matched_strength_surrogate/
        cohort_summary.csv
        per_patient_per_band.csv
        joint_signature.csv
        figures/kc_cohort_distribution.pdf
        figures/kc_per_patient_panel.pdf
        README.md

Surrogate cache reuse (post-2026-05-11)
---------------------------------------
This script regenerates surrogate adjacency matrices on every invocation
because it predates the disk cache introduced in audit_66. New audits
that need the same matched-strength ensemble MUST instead import from
``lrg_eegfc.utils.surrogate.matched_strength`` and call
``load_or_compute_surrogate_eigs(pat, band, phase, W, R, SF, seed, rng)``,
which returns ``(eigvals[R, N], eigvecs[R, N, N])`` cached at
``data/cache/matched_strength_surrogate_lrg/Pat_NN/``
``{band}_{phase}_R{R}_swap{SF}_seed{S}_imcoh_abs.npz``.
Any downstream statistic (KC, Grassmann, ρ propagator, eigenmode embedding)
on the SAME (R, swap_factor, seed) ensemble runs in minutes from cache
instead of hours from regeneration. The canonical R=200, SWAP_FACTOR=20
ensemble for trace bands {α, β, γ_l} is keyed on seed=20260511. This
script's R=200, SWAP_FACTOR=20 ensemble (seed=20260511, three phases
{pre_A, tt, post}) was NOT cached when produced; rerun against the
audit_66 cache to read the same eigendecompositions and recompute KC
T_KC(λ) at λ ∈ {0, 1} in minutes.
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
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BAND_TEX_DICT, FS_OVERRIDES, nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.metrics.tree_distance import kc_distance
from lrg_eegfc.workflow.fc import load_fc_matrix

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # type: ignore


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["alpha", "beta", "low_gamma"]
PHASES = ("rest_pre_A", "task_test", "rest_post")
LAMBDAS = (0.0, 1.0)
N_SURROGATES = 200
SWAP_FACTOR = 20

OUT = ROOT / "data" / "audit" / "kc_matched_strength_surrogate"
FIG = OUT / "figures"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"
AUDIT_63_DIR = ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
HALVES_FC_CACHE.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Half-FC cache (lazy populate, identical to audit_63)
# ---------------------------------------------------------------------------
def _half_fc_path(pat: str, band: str, half: str) -> Path:
    return (HALVES_FC_CACHE / pat
            / f"{band}_rest_pre_{half}_imcoh_abs.npy")


def ensure_half_fcs(pat: str, bands: list[str]) -> None:
    """Compute and cache rsPre half FCs for `bands` if any are missing."""
    missing = [(b, h) for b in bands for h in ("A", "B")
               if not _half_fc_path(pat, b, h).exists()]
    if not missing:
        return
    print(f"[audit_64] {pat}: caching {len(missing)} missing half FCs")
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
    """Load a phase FC matrix (full or half), normalize to (N,N) float64."""
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
# LRG ultrametric → linkage matrix Z (KC needs the full linkage, not cophenet)
# ---------------------------------------------------------------------------
def lrg_linkage(W: np.ndarray) -> np.ndarray:
    """Eigendecomp(L) → ρ(τ=1/λ_max) → Trho=1/ρ → average linkage Z."""
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    eigvals, eigvecs = np.linalg.eigh(L)
    lam_max = eigvals[-1]
    tau = 1.0 / lam_max
    diag_exp = np.exp(-tau * eigvals)
    rho = (eigvecs * diag_exp) @ eigvecs.T
    rho /= np.trace(rho)
    with np.errstate(divide="ignore"):
        Trho = 1.0 / rho
    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0.0)
    finite = np.isfinite(Trho)
    if not finite.all():
        cap = np.nanmax(Trho[finite]) if finite.any() else 1e6
        Trho = np.where(finite, Trho, cap)
    Trho_condensed = squareform(Trho, checks=False)
    Z = linkage(Trho_condensed, method="average")
    return Z


# ---------------------------------------------------------------------------
# Strength-preserving 4-cycle ±δ rewiring (identical to audit_63)
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
# T_KC triangle from three phase linkages
# ---------------------------------------------------------------------------
def t_kc_from_linkages(Z_pre: np.ndarray, Z_tt: np.ndarray,
                        Z_post: np.ndarray, lam: float) -> tuple[float, float, float]:
    """Returns (d_pre_tt, d_tt_post, T_KC).

    T_KC = d(rest_pre, task) - d(task, rest_post). Project-wide convention:
    T_KC > 0 = trace (rsPost closer to task than rsPre is to task).
    """
    d_pre_tt = kc_distance(Z_pre, Z_tt, lam=lam, normalize=True)
    d_tt_post = kc_distance(Z_tt, Z_post, lam=lam, normalize=True)
    return float(d_pre_tt), float(d_tt_post), float(d_pre_tt - d_tt_post)


# ---------------------------------------------------------------------------
# Per-cell pipeline
# ---------------------------------------------------------------------------
def per_cell(pat: str, band: str, n_surr: int, swap_factor: int,
             rng: np.random.Generator, verbose: bool = False) -> dict | None:
    Ws: dict[str, np.ndarray] = {}
    for phase in PHASES:
        try:
            Ws[phase] = load_phase_fc(pat, phase, band)
        except Exception as e:
            if verbose:
                print(f"[audit_64] SKIP {pat}/{band}/{phase}: {e}")
            return None

    N = Ws[PHASES[0]].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_64] SKIP {pat}/{band}: phase shape mismatch")
        return None

    # Observed
    Z_obs = {phase: lrg_linkage(Ws[phase]) for phase in PHASES}
    obs = {}
    for lam in LAMBDAS:
        d_pre_tt, d_tt_post, t_kc = t_kc_from_linkages(
            Z_obs["rest_pre_A"], Z_obs["task_test"],
            Z_obs["rest_post"], lam=lam)
        obs[lam] = {
            "d_pre_tt": d_pre_tt,
            "d_tt_post": d_tt_post,
            "T_KC": t_kc,
        }

    # Surrogates
    n_swaps = swap_factor * (N * (N - 1)) // 2
    surr_T_KC = {lam: np.empty(n_surr, dtype=float) for lam in LAMBDAS}
    surr_d_pre_tt = {lam: np.empty(n_surr, dtype=float) for lam in LAMBDAS}
    surr_d_tt_post = {lam: np.empty(n_surr, dtype=float) for lam in LAMBDAS}
    for r in range(n_surr):
        Z_surr: dict[str, np.ndarray] = {}
        ok = True
        for phase in PHASES:
            W_s = strength_preserving_shuffle(Ws[phase], n_swaps, rng)
            if not verify_strengths(Ws[phase], W_s, tol=1e-4):
                ok = False
                break
            Z_surr[phase] = lrg_linkage(W_s)
        if not ok:
            for lam in LAMBDAS:
                surr_T_KC[lam][r] = np.nan
                surr_d_pre_tt[lam][r] = np.nan
                surr_d_tt_post[lam][r] = np.nan
            continue
        for lam in LAMBDAS:
            dp, dt, tk = t_kc_from_linkages(
                Z_surr["rest_pre_A"], Z_surr["task_test"],
                Z_surr["rest_post"], lam=lam)
            surr_T_KC[lam][r] = tk
            surr_d_pre_tt[lam][r] = dp
            surr_d_tt_post[lam][r] = dt

    # Per-lambda stats
    out_rows: list[dict] = []
    surr_arrays: dict[float, np.ndarray] = {}
    for lam in LAMBDAS:
        s = surr_T_KC[lam]
        s_finite = s[np.isfinite(s)]
        if s_finite.size == 0:
            continue
        surr_arrays[lam] = s_finite
        s_mean = float(np.mean(s_finite))
        s_std = float(np.std(s_finite, ddof=1))
        obs_T = obs[lam]["T_KC"]
        z = (obs_T - s_mean) / s_std if s_std > 0 else float("nan")
        # Trace direction: T_KC > 0; one-sided upper-tail
        p_upper = float(np.mean(s_finite >= obs_T))
        out_rows.append({
            "patient": pat,
            "band": band,
            "lambda": lam,
            "N_nodes": int(N),
            "n_swaps_per_surrogate": int(n_swaps),
            "n_surrogates": int(s_finite.size),
            "obs_T_KC": obs_T,
            "obs_d_pre_tt": obs[lam]["d_pre_tt"],
            "obs_d_tt_post": obs[lam]["d_tt_post"],
            "surr_T_KC_mean": s_mean,
            "surr_T_KC_std": s_std,
            "surr_T_KC_p5": float(np.quantile(s_finite, 0.05)),
            "surr_T_KC_p25": float(np.quantile(s_finite, 0.25)),
            "surr_T_KC_p50": float(np.quantile(s_finite, 0.50)),
            "surr_T_KC_p75": float(np.quantile(s_finite, 0.75)),
            "surr_T_KC_p95": float(np.quantile(s_finite, 0.95)),
            "obs_z": z,
            "obs_p_one_sided_upper": p_upper,
        })
    return {"rows": out_rows, "surr_arrays": surr_arrays,
            "patient": pat, "band": band, "obs": obs}


# ---------------------------------------------------------------------------
# Cohort summary
# ---------------------------------------------------------------------------
def cohort_summary(per_pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    for band in TARGET_BANDS:
        for lam in LAMBDAS:
            sub = per_pat[(per_pat.band == band)
                          & (per_pat["lambda"] == lam)]
            if sub.empty:
                continue
            obs_T = sub.obs_T_KC.values
            surr_med_per_pat = sub.surr_T_KC_p50.values
            n_above = int((sub.obs_p_one_sided_upper < 0.05).sum())
            try:
                wz, wp = wilcoxon(obs_T - surr_med_per_pat,
                                  alternative="greater")
                cohort_z = float(wz)
                cohort_p = float(wp)
            except Exception:
                cohort_z = float("nan")
                cohort_p = float("nan")
            med_obs = float(np.median(obs_T))
            med_surr = float(np.median(surr_med_per_pat))
            p95_surr = float(np.quantile(surr_med_per_pat, 0.95))

            # Verdict per spec (trace direction = positive under project convention):
            # - separated: cohort_p<0.05 AND |median surr T_KC| < 0.05 AND
            #   ≥ 8/10 patients individually above their surrogate at p<0.05.
            # - also_positive: surr median T_KC > 0.5 * obs median T_KC
            #   AND obs median T_KC > 0.
            # - intermediate: anything else.
            if (cohort_p < 0.05 and abs(med_surr) < 0.05 and n_above >= 8):
                verdict = "separated"
            elif med_obs > 0 and med_surr > 0.5 * med_obs:
                verdict = "also_positive"
            else:
                verdict = "intermediate"
            out.append({
                "band": band,
                "lambda": lam,
                "n_patients": len(sub),
                "obs_median_T_KC": med_obs,
                "surr_median_T_KC_per_patient_median": med_surr,
                "surr_median_T_KC_per_patient_p95": p95_surr,
                "n_patients_above_own_surrogate": f"{n_above}/{len(sub)}",
                "paired_wilcoxon_z": cohort_z,
                "paired_wilcoxon_p": cohort_p,
                "verdict": verdict,
            })
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Joint signature with audit_63 ρ_split
# ---------------------------------------------------------------------------
def joint_signature(per_pat: pd.DataFrame) -> pd.DataFrame:
    audit_63_csv = AUDIT_63_DIR / "per_patient_per_band.csv"
    if not audit_63_csv.exists():
        print(f"[audit_64] WARN: audit_63 per_patient CSV not found at {audit_63_csv}")
        a63 = pd.DataFrame(columns=["patient", "band", "obs_z"])
    else:
        a63 = pd.read_csv(audit_63_csv)
    rows = []
    for pat in COHORT:
        for band in TARGET_BANDS:
            l0 = per_pat[(per_pat.patient == pat) & (per_pat.band == band)
                         & (per_pat["lambda"] == 0.0)]
            l1 = per_pat[(per_pat.patient == pat) & (per_pat.band == band)
                         & (per_pat["lambda"] == 1.0)]
            r63 = a63[(a63.patient == pat) & (a63.band == band)]
            rho_z = float(r63.iloc[0]["obs_z"]) if not r63.empty else float("nan")
            rho_p = float(r63.iloc[0]["obs_p_one_sided"]) if (
                not r63.empty and "obs_p_one_sided" in r63.columns
            ) else float("nan")
            kc0_z = float(l0.iloc[0]["obs_z"]) if not l0.empty else float("nan")
            kc0_p = float(l0.iloc[0]["obs_p_one_sided_upper"]) if not l0.empty else float("nan")
            kc1_z = float(l1.iloc[0]["obs_z"]) if not l1.empty else float("nan")
            kc1_p = float(l1.iloc[0]["obs_p_one_sided_upper"]) if not l1.empty else float("nan")

            # n_probes_significant_at_p05
            # - rho_split: one-sided upper-tail at p<0.05 (audit_63 obs_p)
            # - T_KC λ=0/1: one-sided upper-tail at p<0.05 (positive = trace)
            n_sig = 0
            if np.isfinite(rho_p) and rho_p < 0.05:
                n_sig += 1
            if np.isfinite(kc0_p) and kc0_p < 0.05:
                n_sig += 1
            if np.isfinite(kc1_p) and kc1_p < 0.05:
                n_sig += 1
            rows.append({
                "patient": pat,
                "band": band,
                "rho_split_z": rho_z,
                "rho_split_p_one_sided_upper": rho_p,
                "T_KC_l0_z": kc0_z,
                "T_KC_l0_p_one_sided_upper": kc0_p,
                "T_KC_l1_z": kc1_z,
                "T_KC_l1_p_one_sided_upper": kc1_p,
                "n_probes_significant_at_p05": n_sig,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def make_cohort_figure(per_pat: pd.DataFrame, cohort: pd.DataFrame,
                        out_path: Path) -> None:
    n_lam = len(LAMBDAS)
    n_b = len(TARGET_BANDS)
    fig, axes = plt.subplots(n_lam, n_b, figsize=(4.5 * n_b, 4.4 * n_lam),
                             sharey=False)
    if n_lam == 1:
        axes = np.array([axes])
    if n_b == 1:
        axes = axes.reshape(n_lam, 1)
    for i, lam in enumerate(LAMBDAS):
        for j, band in enumerate(TARGET_BANDS):
            ax = axes[i, j]
            sub = per_pat[(per_pat.band == band)
                          & (per_pat["lambda"] == lam)].copy()
            if sub.empty:
                ax.set_axis_off()
                continue
            sub = sub.sort_values("obs_T_KC", ascending=True).reset_index(drop=True)
            ys = np.arange(len(sub))
            for y, (_, r) in zip(ys, sub.iterrows()):
                ax.hlines(y, r.surr_T_KC_p5, r.surr_T_KC_p95,
                          color="#aaaaaa", lw=4, alpha=0.85)
                ax.plot(r.surr_T_KC_p50, y, marker="|", color="#444444",
                        markersize=10, mew=1.2)
                ax.plot(r.obs_T_KC, y, marker="o", color="#1f3d6e",
                        markersize=7, mec="white", mew=0.8)
            ax.set_yticks(ys)
            ax.set_yticklabels([p.replace("Pat_", "P")
                                for p in sub.patient], fontsize=8)
            ax.axvline(0, color="0.6", lw=0.7, ls="--", zorder=0)
            cohort_row = cohort[(cohort.band == band)
                                & (cohort["lambda"] == lam)].iloc[0]
            verdict = str(cohort_row["verdict"])
            wp = float(cohort_row["paired_wilcoxon_p"])
            ax.set_xlabel(rf"$T_{{KC}}(\lambda={lam:g})$")
            ax.set_title(
                f"{BRAIN_BAND_TEX_DICT[band]} — verdict: {verdict}\n"
                rf"Wilcoxon $p={wp:.4f}$,  $n_{{<}} = "
                rf"{cohort_row['n_patients_above_own_surrogate']}$",
                fontsize=10)
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_per_patient_panel(cell_results: list[dict],
                            per_pat: pd.DataFrame,
                            out_path: Path) -> None:
    n_lam = len(LAMBDAS)
    n_pat = len(COHORT)
    n_b = len(TARGET_BANDS)
    fig, axes = plt.subplots(n_pat, n_b * n_lam,
                             figsize=(2.8 * n_b * n_lam, 1.55 * n_pat),
                             sharex=False)
    by_pat_band: dict[tuple[str, str], dict] = {}
    for cell in cell_results:
        by_pat_band[(cell["patient"], cell["band"])] = cell
    for i, pat in enumerate(COHORT):
        for j, band in enumerate(TARGET_BANDS):
            for k, lam in enumerate(LAMBDAS):
                col = j * n_lam + k
                ax = axes[i, col]
                cell = by_pat_band.get((pat, band))
                if cell is None or lam not in cell["surr_arrays"]:
                    ax.set_axis_off()
                    continue
                surr = cell["surr_arrays"][lam]
                obs_T = cell["obs"][lam]["T_KC"]
                ax.hist(surr, bins=30, color="#cccccc",
                        edgecolor="#888888", alpha=0.85)
                ax.axvline(obs_T, color="#1f3d6e", lw=1.6,
                           label=f"obs={obs_T:+.3f}")
                ax.axvline(0, color="#888888", lw=0.6, ls="--")
                if i == n_pat - 1:
                    ax.set_xlabel(rf"$T_{{KC}}(\lambda={lam:g})$")
                if col == 0:
                    ax.set_ylabel(pat.replace("Pat_", "P"), fontsize=8)
                if i == 0:
                    ax.set_title(
                        f"{BRAIN_BAND_TEX_DICT[band]} $\\lambda={lam:g}$",
                        fontsize=9)
                ax.tick_params(labelsize=7)
                ax.legend(loc="upper left", frameon=False, fontsize=7)
                ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(cohort: pd.DataFrame, joint: pd.DataFrame,
                 n_surr: int, swap_factor: int, runtime_s: float) -> None:
    lines = [
        "---",
        "name: kc_matched_strength_surrogate",
        "scope: section_5_2_kc_robustness_against_matched_strength_null",
        "date: 2026-05-11",
        "status: first_pass",
        "---",
        "",
        "# KC tree-distance under matched-strength surrogacy",
        "",
        f"**Head.** R={n_surr} strength-preserving (4-cycle ±δ, "
        f"SWAP_FACTOR={swap_factor}) surrogate Laplacians per (patient, "
        f"band ∈ {{α, β, γ_l}}, phase ∈ {{rsPre_A, taskTest, rsPost}}). "
        f"For each surrogate, KC tree distance T_KC(λ) is computed at "
        f"λ=0 (pure topology) and λ=1 (pure heights). Observed counterpart "
        f"uses D^pre_A (first half of rsPre) as the rsPre leg to match "
        f"the surrogate baseline. The §5.2 published 10/10 cohort claim "
        f"used D^pre (full rsPre) — the split-baseline observed numbers "
        f"may differ slightly; the comparison here is split-baseline "
        f"observed vs split-baseline surrogate, both with the same baseline.",
        "",
        "## Cohort verdict",
        "",
        "| band | λ | obs median T_KC | surr median (per-pat med) | n below own surrogate | Wilcoxon p | verdict |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, r in cohort.iterrows():
        lines.append(
            f"| {r['band']} | {r['lambda']:g} "
            f"| {r['obs_median_T_KC']:+.3f} "
            f"| {r['surr_median_T_KC_per_patient_median']:+.3f} "
            f"| {r['n_patients_above_own_surrogate']} "
            f"| {r['paired_wilcoxon_p']:.4f} "
            f"| **{r['verdict']}** |"
        )
    lines.extend([
        "",
        "## Verdict labels (spec)",
        "",
        "- **separated**: Wilcoxon p < 0.05 (one-sided, T_obs < T_surr) AND "
        "|median surr T_KC across patients| < 0.05 AND ≥ 8/10 patients "
        "individually below their own surrogate at p<0.05 lower-tail.",
        "- **also_negative**: median surr T_KC across patients < 0.5 × "
        "median obs T_KC AND median obs T_KC < 0 (the surrogate "
        "reproduces or substantially recovers the observed trace direction).",
        "- **intermediate**: anything else — both per-patient counts and "
        "cohort-paired statistic must be inspected verbally.",
        "",
        "Note on the ticket spec: ticket asks for 'cohort z > 2' as part of "
        "the separated criterion. The paired Wilcoxon test returns a rank "
        "sum statistic, not a z-score. We use cohort p < 0.05 as the "
        "equivalent tail criterion. The Wilcoxon z column is still saved "
        "for cross-reference but is the raw scipy `statistic` (rank sum).",
        "",
        "## Joint signature across three probes",
        "",
        "Per (patient, band), we collect z-scores from:",
        "",
        "- ρ_split (audit_63): one-sided upper-tail (large positive ρ_split "
        "= trace).",
        "- T_KC(λ=0): one-sided upper-tail (large positive T_KC = trace).",
        "- T_KC(λ=1): one-sided upper-tail.",
        "",
        "`joint_signature.csv` columns:",
        "`patient, band, rho_split_z, rho_split_p_one_sided_upper, "
        "T_KC_l0_z, T_KC_l0_p_one_sided_upper, T_KC_l1_z, "
        "T_KC_l1_p_one_sided_upper, n_probes_significant_at_p05`",
        "",
        "## Algorithm choice (vs §5.2 within-baseline null)",
        "",
        "The §5.2 within-baseline null is built from a single rsPre split-half: "
        "T_KC^null = d_KC(D^pre_A, D^pre_B) − d_KC(D^pre_B, D^pre_A) ≈ 0 by "
        "construction. The matched-strength surrogate here is wider: each "
        "phase's adjacency is independently rewired to preserve the per-node "
        "strength sequence, so T_KC_surr distributes around what a typical "
        "matched-strength random Laplacian gives across the triangle. A "
        "'separated' verdict means the §5.2 KC trace direction is not "
        "explained by per-phase strength heterogeneity propagating through "
        "the average-linkage construction.",
        "",
        "Caveat (independent per-phase rewiring vs cross-phase coordination): "
        "matching audit_63, the surrogate is constructed independently per "
        "phase. There is no preserved cross-phase identity at the edge level. "
        "The surrogate therefore decorrelates phases at the edge level by "
        "construction, which can pull T_KC_surr toward zero relative to "
        "observed (which has cross-phase coherence). This is intrinsic to "
        "the strength-preserving family and is the same caveat as in audit_63.",
        "",
        "## Provenance",
        "",
        f"- N_surrogates = {n_surr} per cell",
        f"- n_swaps_per_surrogate = SWAP_FACTOR ({swap_factor}) × N(N-1)/2",
        "- Cohort: " + ", ".join(COHORT),
        "- Bands: " + ", ".join(TARGET_BANDS),
        "- Phases: " + ", ".join(PHASES),
        "- λ values: " + ", ".join(f"{x:g}" for x in LAMBDAS),
        "- KC normalization: True (Kendall & Colijn 2016 default)",
        "- FC method: imcoh_abs",
        "- LRG: τ = 1/λ_max, ultrametric via average linkage on Trho.",
        "- Half-FC cache reused: `data/cache/imcoh_halves_fc/Pat_NN/`.",
        "- audit_63 ρ_split z-scores joined from "
        "`data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv`",
        f"- Wall-clock runtime: {runtime_s:.1f} s",
        "- Build script: `scripts/01_compute/audit/audit_64_kc_matched_strength_surrogate.py`",
        "",
        "## Files",
        "",
        "- `per_patient_per_band.csv` — observed + per-cell surrogate stats (one row per patient × band × λ)",
        "- `cohort_summary.csv` — per-(band, λ) cohort verdict",
        "- `joint_signature.csv` — three-probe z-score table per (patient, band)",
        "- `figures/kc_cohort_distribution.pdf` — per-(band, λ) cohort overlay",
        "- `figures/kc_per_patient_panel.pdf` — per-(patient, band, λ) histograms",
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

    print("[audit_64] pre-flight: ensure half FCs cached")
    for pat in patients:
        ensure_half_fcs(pat, bands)

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    cell_results: list[dict] = []
    all_rows: list[dict] = []
    for band in bands:
        for pat in patients:
            t_cell = time.time()
            cell = per_cell(pat, band, n_surr, swap_factor, rng,
                             verbose=args.verbose)
            if cell is None:
                continue
            cell_results.append(cell)
            all_rows.extend(cell["rows"])
            dt = time.time() - t_cell
            l0 = next((r for r in cell["rows"]
                        if r["lambda"] == 0.0), None)
            l1 = next((r for r in cell["rows"]
                        if r["lambda"] == 1.0), None)
            l0s = (f"l0: T={l0['obs_T_KC']:+.3f} z={l0['obs_z']:+.2f} "
                   f"p={l0['obs_p_one_sided_upper']:.3f}"
                   if l0 else "l0: NA")
            l1s = (f"l1: T={l1['obs_T_KC']:+.3f} z={l1['obs_z']:+.2f} "
                   f"p={l1['obs_p_one_sided_upper']:.3f}"
                   if l1 else "l1: NA")
            print(f"[audit_64] {pat}/{band}: {l0s} | {l1s} ({dt:.1f}s)")

    runtime = time.time() - t0
    print(f"[audit_64] total: {runtime:.1f}s for {len(cell_results)} cells")

    per_pat = pd.DataFrame(all_rows)
    per_pat.to_csv(OUT / "per_patient_per_band.csv", index=False)

    cohort = cohort_summary(per_pat)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)
    print(cohort.to_string(index=False))

    joint = joint_signature(per_pat)
    joint.to_csv(OUT / "joint_signature.csv", index=False)

    make_cohort_figure(per_pat, cohort, FIG / "kc_cohort_distribution.pdf")
    make_per_patient_panel(cell_results, per_pat,
                            FIG / "kc_per_patient_panel.pdf")

    write_readme(cohort, joint, n_surr, swap_factor, runtime)
    print(f"[audit_64] outputs at {OUT}")


if __name__ == "__main__":
    main()
