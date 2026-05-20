#!/usr/bin/env python3
"""Audit 67 — Grassmann matched-strength surrogate AFTER epileptic-zone
exclusion (§5.4 γ_h sensitivity to HFO contamination, generalized to all
bands).

Critical preamble (per CLAUDE.md rule, before any code)
=======================================================

(1) **Claim under test.** The audit_66 Grassmann trace T_G(k) =
    d_chord(U^tt, U^post; k) - d_chord(U^pre_A, U^tt; k) survives R=200
    4-cycle ±δ matched-strength surrogacy at γ_h on 9 strict-separated
    k cells (k=19..27, audit_66 verdict). γ_h spans 80-300 Hz in this
    project (physiological broadband gamma + hippocampal ripple range
    80-200 Hz overlapping memory consolidation literature, plus fast-
    gamma / lower HFO 200-300 Hz overlapping epileptogenic HFO band).
    The band-averaged |ImCoh| statistic cannot separate these regimes.
    This audit asks: does the γ_h trace survive after dropping
    epileptic-zone contacts from the FC matrix entirely? If yes → the
    trace is carried by non-epi contacts → physiological-consolidation
    interpretation strengthened. If no → the trace requires
    epileptic-zone participation → HFO-contamination plausible.
    Extended to all 6 bands for full sensitivity surface.

(2) **The null.** Matched-strength 4-cycle ±δ surrogate on the REDUCED
    FC W_red = W[non_epi, non_epi] for each phase, R=200 per cell.
    Same algorithm as audit_63/65/66. Per-node strengths in the
    non-epi subgraph are preserved.

(3) **Strongest plausible alternative.** The γ_h trace could be (a)
    purely physiological signal robust to dropping epi nodes; (b)
    pure HFO contamination carried by epi contacts → dies under
    exclusion; (c) interface phenomenon (epi nodes couple to non-epi
    during memory consolidation, both physiological and pathological,
    edges removed under exclusion) → also dies under exclusion even
    if physiological; (d) genuinely physiological signal that
    happens to be most strongly localized to epi-adjacent tissue.

(4) **Null's mechanical reach.**
    - The matched-strength null on W_red controls for amplitude / per-
      node strength effects within the non-epi subset. If the trace
      direction in W_red survives the null, it is a property of the
      shape of the non-epi subnetwork, not its strength distribution.
    - The null does NOT control for (c): interface effects die under
      exclusion regardless of pathological origin. A negative result
      here is therefore weak evidence against (a); it does not
      cleanly separate (b) from (c).
    - The k axis is NOT directly comparable to audit_66. audit_66 ran
      k=2..112 on full FC; audit_67 runs k=2..88 (min N_reduced=89,
      driven by Pat_13 with 30 epi contacts). Relative k ranges
      (e.g., γ_h "intermediate-k 19-27") map approximately but not
      exactly to the same scale.
    - Pat_15 has 0 epi contacts → epi-exclusion is identity. Pat_15's
      audit_67 row should match audit_66 row up to surrogate
      regeneration noise.

(5) **Falsification + limitations.**
    - Falsification of "γ_h trace is physiological": if γ_h cohort
      Wilcoxon p > 0.05 across all of k=19..27 after exclusion AND
      n_below << audit_66 baseline. Then γ_h trace requires epi-
      contact participation; physiological-consolidation claim is
      weakened.
    - Strengthening: if γ_h survives in non-epi subset with cohort
      p<0.05 at ≥half the audit_66-significant k cells, physiological
      interpretation is robust to HFO ablation.
    - Limitations: (i) Pat_13 30-epi heavy load may dominate cohort
      Wilcoxon; (ii) varying N_reduced per patient subtly biases
      cross-patient k-comparability; (iii) interface effects (c)
      indistinguishable from (b) under this design; (iv) we can't
      rule out residual HFO leakage from non-epi contacts spatially
      adjacent to epileptic foci.

Outputs
-------
``data/audit/grassmann_epi_exclusion/``
    per_patient_per_band_per_k.csv     observed + surrogate stats
    cohort_summary.csv                 per-(band, k) cohort verdict
    epi_counts.csv                     N, N_epi, N_reduced per patient
    sensitivity.csv                    audit_66 vs audit_67 verdict join
    figures/grassmann_cohort_distribution_epiX.pdf
    figures/grassmann_kspan_verdict_epiX.pdf
    figures/grassmann_kspan_sensitivity_vs_audit66.pdf
    README.md

Surrogate cache (separate from audit_66 to avoid clobber)
---------------------------------------------------------
``data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/
    {band}_{phase}_epiX_R{R}_swap{SF}_seed{S}_imcoh_abs.npz``
contains ``eigvals (R, N_reduced)``, ``eigvecs (R, N_reduced, N_reduced)``.

For Pat_15 (n_epi = 0), W_red = W and the surrogate ensemble is
mathematically equivalent to audit_66's Pat_15 cache up to RNG state.
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
from lrg_eegfc.utils.io import load_channel_labels, load_epileptic_nodes
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.surrogate.matched_strength import (
    strength_preserving_shuffle, verify_strengths,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # type: ignore


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
ALL_BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
CANONICAL_BAND_ORDER = ALL_BANDS  # identical
PHASES = ("rest_pre_A", "task_test", "rest_post")

# Preflight (2026-05-14): min N_reduced = 89 (Pat_13 has 30 epi contacts).
# K_GRID upper bound set conservatively to 88. γ_h finding (k=19..27 in
# audit_66) and β manuscript range (k=27..55) both fully covered.
K_GRID = list(range(2, 89))
K_REFERENCE = (3, 20, 25, 40, 60, 80)

N_SURROGATES = 200
SWAP_FACTOR = 20

OUT = ROOT / "data" / "audit" / "grassmann_epi_exclusion"
FIG = OUT / "figures"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"
SURROGATE_LRG_CACHE_EPIX = (CACHE_ROOT
                             / "matched_strength_surrogate_epi_excluded_lrg")
AUDIT_66_DIR = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
SURROGATE_LRG_CACHE_EPIX.mkdir(parents=True, exist_ok=True)


def _band_order_present(per_pat: pd.DataFrame) -> list[str]:
    present = set(per_pat.band.unique())
    return [b for b in CANONICAL_BAND_ORDER if b in present]


# ---------------------------------------------------------------------------
# Epi-exclusion mask per patient
# ---------------------------------------------------------------------------
def epi_keep_mask(pat: str, n_full: int) -> np.ndarray:
    """Return bool mask of length n_full, True for non-epi contacts."""
    labels = load_channel_labels(pat)
    if len(labels) != n_full:
        raise ValueError(
            f"{pat}: channel_labels length {len(labels)} != FC N {n_full}"
        )
    epi_set = set(load_epileptic_nodes(pat))
    keep = np.array([lbl not in epi_set for lbl in labels], dtype=bool)
    return keep


# ---------------------------------------------------------------------------
# Half-FC cache (reused from audit_66; computed on FULL recording, then
# masked at load time)
# ---------------------------------------------------------------------------
def _half_fc_path(pat: str, band: str, half: str) -> Path:
    return (HALVES_FC_CACHE / pat
            / f"{band}_rest_pre_{half}_imcoh_abs.npy")


def ensure_half_fcs(pat: str, bands: list[str]) -> None:
    missing = [(b, h) for b in bands for h in ("A", "B")
               if not _half_fc_path(pat, b, h).exists()]
    if not missing:
        return
    print(f"[audit_67] {pat}: caching {len(missing)} missing half FCs")
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


def load_phase_fc_full(pat: str, phase: str, band: str) -> np.ndarray:
    """Load FULL FC (before epi-exclusion), symmetric, [0, 1], zero-diag."""
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
# Laplacian + Grassmann (identical to audit_66)
# ---------------------------------------------------------------------------
def laplacian_eig(W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    return np.linalg.eigh(L)


def topk_basis(eigvecs: np.ndarray, k_max: int) -> np.ndarray:
    return np.ascontiguousarray(eigvecs[:, 1:k_max + 1])


def chordal(V_a: np.ndarray, V_b: np.ndarray, k: int) -> float:
    A = V_a[:, :k]
    B = V_b[:, :k]
    M = A.T @ B
    sigma = np.linalg.svd(M, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    return float(np.sqrt(max(k - float((sigma ** 2).sum()), 0.0)))


def t_g_at_all_k(U_pre: np.ndarray, U_tt: np.ndarray,
                  U_post: np.ndarray, k_grid: list[int]) -> np.ndarray:
    out = np.empty(len(k_grid), dtype=float)
    for i, k in enumerate(k_grid):
        d_pre_tt = chordal(U_pre, U_tt, k)
        d_tt_post = chordal(U_tt, U_post, k)
        out[i] = d_tt_post - d_pre_tt
    return out


# ---------------------------------------------------------------------------
# Surrogate eigendecomposition cache (epi-excluded variant)
# ---------------------------------------------------------------------------
def _surr_cache_path_epiX(pat: str, band: str, phase: str,
                           n_surr: int, swap_factor: int, seed: int) -> Path:
    pat_dir = SURROGATE_LRG_CACHE_EPIX / pat
    return (pat_dir / f"{band}_{phase}_epiX_R{n_surr}_swap{swap_factor}"
                       f"_seed{seed}_imcoh_abs.npz")


def load_or_compute_surrogate_eigs_epiX(pat: str, band: str, phase: str,
                                         W_red: np.ndarray, n_surr: int,
                                         swap_factor: int, seed: int,
                                         rng: np.random.Generator,
                                         verbose: bool = False
                                         ) -> tuple[np.ndarray, np.ndarray]:
    path = _surr_cache_path_epiX(pat, band, phase, n_surr, swap_factor, seed)
    N = W_red.shape[0]
    if path.exists():
        with np.load(path) as data:
            evals = data["eigvals"]
            evecs = data["eigvecs"]
        if (evals.shape == (n_surr, N) and evecs.shape == (n_surr, N, N)):
            return evals.astype(np.float64), evecs.astype(np.float64)
        if verbose:
            print(f"[audit_67] cache shape mismatch at {path}; recomputing")

    n_swaps = swap_factor * (N * (N - 1)) // 2
    evals = np.empty((n_surr, N), dtype=np.float64)
    evecs = np.empty((n_surr, N, N), dtype=np.float64)
    ok_count = 0
    for r in range(n_surr):
        W_s = strength_preserving_shuffle(W_red, n_swaps, rng)
        if not verify_strengths(W_red, W_s, tol=1e-4):
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
        print(f"[audit_67] cached {ok_count}/{n_surr} surrogates -> {path.name}")
    return evals, evecs


# ---------------------------------------------------------------------------
# Per-cell pipeline (one (patient, band))
# ---------------------------------------------------------------------------
def per_cell(pat: str, band: str, n_surr: int, swap_factor: int,
             k_grid: list[int], seed: int, rng: np.random.Generator,
             verbose: bool = False) -> dict | None:
    Ws_full: dict[str, np.ndarray] = {}
    for phase in PHASES:
        try:
            Ws_full[phase] = load_phase_fc_full(pat, phase, band)
        except Exception as e:
            if verbose:
                print(f"[audit_67] SKIP {pat}/{band}/{phase}: {e}")
            return None

    N_full = Ws_full[PHASES[0]].shape[0]
    if any(W.shape != (N_full, N_full) for W in Ws_full.values()):
        return None

    keep = epi_keep_mask(pat, N_full)
    N_red = int(keep.sum())
    n_epi = N_full - N_red
    Ws_red = {phase: Ws_full[phase][np.ix_(keep, keep)] for phase in PHASES}

    k_grid_eff = [k for k in k_grid if k + 1 <= N_red]
    if not k_grid_eff:
        if verbose:
            print(f"[audit_67] SKIP {pat}/{band}: N_red={N_red} < min k+1")
        return None
    k_max = max(k_grid_eff)

    # Observed (split-baseline; pre_A as rsPre leg)
    U_obs = {phase: topk_basis(laplacian_eig(Ws_red[phase])[1], k_max)
              for phase in PHASES}
    obs_T_G = t_g_at_all_k(U_obs["rest_pre_A"], U_obs["task_test"],
                            U_obs["rest_post"], k_grid_eff)

    # Surrogates — load from epi-excluded cache or compute + cache
    n_swaps = swap_factor * (N_red * (N_red - 1)) // 2
    surr_eigvecs: dict[str, np.ndarray] = {}
    for phase in PHASES:
        _, evecs = load_or_compute_surrogate_eigs_epiX(
            pat, band, phase, Ws_red[phase],
            n_surr, swap_factor, seed, rng, verbose=verbose)
        surr_eigvecs[phase] = evecs

    surr_T_G = np.empty((n_surr, len(k_grid_eff)), dtype=float)
    for r in range(n_surr):
        U_surr = {phase: topk_basis(surr_eigvecs[phase][r], k_max)
                  for phase in PHASES}
        if any(np.isnan(U_surr[phase]).any() for phase in PHASES):
            surr_T_G[r, :] = np.nan
            continue
        surr_T_G[r, :] = t_g_at_all_k(
            U_surr["rest_pre_A"], U_surr["task_test"],
            U_surr["rest_post"], k_grid_eff)

    rows: list[dict] = []
    for i, k in enumerate(k_grid_eff):
        s = surr_T_G[:, i]
        s_finite = s[np.isfinite(s)]
        if s_finite.size == 0:
            continue
        s_mean = float(np.mean(s_finite))
        s_std = float(np.std(s_finite, ddof=1))
        obs_T = float(obs_T_G[i])
        z = (obs_T - s_mean) / s_std if s_std > 0 else float("nan")
        p_low = float(np.mean(s_finite <= obs_T))
        rows.append({
            "patient": pat,
            "band": band,
            "k": int(k),
            "N_full": int(N_full),
            "N_epi": int(n_epi),
            "N_reduced": int(N_red),
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
            "obs_p_one_sided_lower": p_low,
        })
    return {"rows": rows, "patient": pat, "band": band}


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
            n_below = int((sub.obs_p_one_sided_lower < 0.05).sum())
            try:
                wz, wp = wilcoxon(obs_T - surr_med, alternative="less")
                cohort_z = float(wz)
                cohort_p = float(wp)
            except Exception:
                cohort_z = float("nan")
                cohort_p = float("nan")
            med_obs = float(np.median(obs_T))
            med_surr = float(np.median(surr_med))
            if (cohort_p < 0.05 and abs(med_surr) < 0.05 * max(1.0, abs(med_obs))
                    and n_below >= 8):
                verdict = "separated"
            elif med_obs < 0 and med_surr < 0.5 * med_obs:
                verdict = "also_negative"
            else:
                verdict = "intermediate"
            out.append({
                "band": band,
                "k": k,
                "n_patients": len(sub),
                "obs_median_T_G": med_obs,
                "surr_median_T_G_per_patient_median": med_surr,
                "n_patients_below_own_surrogate": n_below,
                "n_patients_below_own_surrogate_str": f"{n_below}/{len(sub)}",
                "paired_wilcoxon_z": cohort_z,
                "paired_wilcoxon_p": cohort_p,
                "verdict": verdict,
            })
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Sensitivity table: audit_66 cohort verdict vs audit_67 cohort verdict
# at each (band, k) cell that exists in both.
# ---------------------------------------------------------------------------
def sensitivity_table(cohort_67: pd.DataFrame) -> pd.DataFrame:
    a66_path = AUDIT_66_DIR / "cohort_summary.csv"
    if not a66_path.exists():
        print(f"[audit_67] WARN: audit_66 cohort_summary missing at {a66_path}"
              "; sensitivity table will only have audit_67 columns.")
        sens = cohort_67.copy()
        sens.columns = [f"{c}_epiX" for c in sens.columns]
        return sens
    a66 = pd.read_csv(a66_path)
    rows = []
    for _, r67 in cohort_67.iterrows():
        match = a66[(a66.band == r67.band) & (a66.k == r67.k)]
        if match.empty:
            r66 = {}
        else:
            r66 = match.iloc[0].to_dict()
        rows.append({
            "band": r67["band"],
            "k": int(r67["k"]),
            # audit_66 (full FC)
            "full_obs_median_T_G": r66.get("obs_median_T_G", np.nan),
            "full_surr_median_T_G": r66.get(
                "surr_median_T_G_per_patient_median", np.nan),
            "full_n_below": r66.get("n_patients_below_own_surrogate", np.nan),
            "full_paired_wilcoxon_p": r66.get("paired_wilcoxon_p", np.nan),
            "full_verdict": r66.get("verdict", "missing"),
            # audit_67 (epi-excluded)
            "epiX_obs_median_T_G": r67["obs_median_T_G"],
            "epiX_surr_median_T_G": r67[
                "surr_median_T_G_per_patient_median"],
            "epiX_n_below": r67["n_patients_below_own_surrogate"],
            "epiX_paired_wilcoxon_p": r67["paired_wilcoxon_p"],
            "epiX_verdict": r67["verdict"],
        })
    sens = pd.DataFrame(rows)
    # Sensitivity flags: did the trace strengthen / weaken / persist / die?
    def _flag(row):
        f = str(row["full_verdict"])
        e = str(row["epiX_verdict"])
        fp = row["full_paired_wilcoxon_p"]
        ep = row["epiX_paired_wilcoxon_p"]
        try:
            fp = float(fp); ep = float(ep)
        except Exception:
            return "unknown"
        if not (np.isfinite(fp) and np.isfinite(ep)):
            return "unknown"
        if fp < 0.05 and ep < 0.05:
            return "persist"
        if fp < 0.05 and ep >= 0.05:
            return "weaken"
        if fp >= 0.05 and ep < 0.05:
            return "emerge"
        return "absent"
    sens["sensitivity_flag"] = sens.apply(_flag, axis=1)
    return sens


# ---------------------------------------------------------------------------
# Epi-count summary (one row per patient)
# ---------------------------------------------------------------------------
def write_epi_counts() -> pd.DataFrame:
    rows = []
    for pat in COHORT:
        labels = load_channel_labels(pat)
        epi_set = set(load_epileptic_nodes(pat))
        n_full = len(labels)
        n_epi = sum(1 for l in labels if l in epi_set)
        rows.append({
            "patient": pat,
            "N_full": n_full,
            "N_epi": n_epi,
            "N_reduced": n_full - n_epi,
            "epi_fraction": n_epi / n_full,
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "epi_counts.csv", index=False)
    return df


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
        p5, p95 = [], []
        for k in ks:
            cells = per_pat[(per_pat.band == band) & (per_pat.k == k)]
            obs_vals = cells.obs_T_G.values
            p5.append(np.quantile(obs_vals, 0.05) if obs_vals.size else np.nan)
            p95.append(np.quantile(obs_vals, 0.95) if obs_vals.size else np.nan)
        ax.fill_between(ks, p5, p95, color="#1f3d6e", alpha=0.18,
                        label="obs 5-95%")
        ax.plot(ks, obs_med, color="#1f3d6e", lw=1.4, label="obs median")
        ax.plot(ks, surr_med, color="#aa4444", lw=1.2,
                label="surr median (per-pat med)")
        ax.axhline(0, color="0.6", lw=0.6, ls="--")
        ax.set_ylabel(rf"$T_G$ - {BRAIN_BAND_TEX_DICT[band]}")
        ax.spines[["top", "right"]].set_visible(False)
        if band == bands[0]:
            ax.legend(frameon=False, fontsize=8, loc="lower right")
    axes[-1].set_xlabel(r"spectral cutoff $k$ (epi-excluded)")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_kspan_verdict_figure(cohort: pd.DataFrame, out_path: Path,
                               title_suffix: str = "epi-excluded") -> None:
    bands = [b for b in CANONICAL_BAND_ORDER if b in set(cohort.band.unique())]
    n_b = len(bands)
    fig, axes = plt.subplots(n_b, 1, figsize=(8.5, 1.8 * n_b), sharex=True)
    if n_b == 1:
        axes = [axes]
    color_map = {"separated": "#1f7a1f", "intermediate": "#888888",
                 "also_negative": "#aa4444"}
    for ax, band in zip(axes, bands):
        sub = cohort[cohort.band == band].sort_values("k")
        if sub.empty:
            ax.set_axis_off()
            continue
        for _, r in sub.iterrows():
            ax.bar(int(r.k), 1, width=1.0,
                   color=color_map.get(str(r.verdict), "#cccccc"),
                   edgecolor="none")
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
    axes[-1].set_xlabel(rf"spectral cutoff $k$ ({title_suffix})")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, label=l)
               for l, c in color_map.items()]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_sensitivity_figure(sens: pd.DataFrame, out_path: Path) -> None:
    """Side-by-side audit_66 vs audit_67 p-curves, per band, with flag strip."""
    bands = [b for b in CANONICAL_BAND_ORDER if b in set(sens.band.unique())]
    n_b = len(bands)
    fig, axes = plt.subplots(n_b, 1, figsize=(8.5, 2.4 * n_b), sharex=True)
    if n_b == 1:
        axes = [axes]
    for ax, band in zip(axes, bands):
        sub = sens[sens.band == band].sort_values("k")
        if sub.empty:
            ax.set_axis_off()
            continue
        ks = sub.k.values
        fp = sub.full_paired_wilcoxon_p.values
        ep = sub.epiX_paired_wilcoxon_p.values
        ax.plot(ks, fp, color="#1f3d6e", lw=1.2, label="full FC (audit_66)")
        ax.plot(ks, ep, color="#cc6633", lw=1.2, label="epi-excluded")
        ax.axhline(0.05, color="#aa4444", lw=0.5, ls="--")
        ax.set_ylim(0, 1.0)
        ax.set_ylabel(rf"cohort $p$ - {BRAIN_BAND_TEX_DICT[band]}")
        ax.spines[["top", "right"]].set_visible(False)
        # Flag strip across bottom
        flag_color = {"persist": "#1f7a1f", "weaken": "#aa4444",
                       "emerge": "#cc6633", "absent": "#cccccc",
                       "unknown": "#999999"}
        for k, f in zip(ks, sub.sensitivity_flag.values):
            ax.bar(int(k), 0.03, bottom=-0.04, width=1.0,
                   color=flag_color.get(str(f), "#999999"), edgecolor="none")
        if band == bands[0]:
            ax.legend(frameon=False, fontsize=8, loc="upper right")
    axes[-1].set_xlabel(r"spectral cutoff $k$")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, label=l)
               for l, c in {"persist": "#1f7a1f", "weaken": "#aa4444",
                              "emerge": "#cc6633", "absent": "#cccccc"}.items()]
    fig.legend(handles=handles, loc="upper center", ncol=4,
               frameon=False, bbox_to_anchor=(0.5, 1.005))
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(cohort: pd.DataFrame, per_pat: pd.DataFrame,
                 sens: pd.DataFrame, epi_counts: pd.DataFrame,
                 n_surr: int, swap_factor: int,
                 runtime_s: float) -> None:
    lines = [
        "---",
        "name: grassmann_epi_exclusion",
        "scope: section_5_4_grassmann_robustness_to_epileptic_zone_exclusion",
        "date: 2026-05-14",
        "status: first_pass",
        "---",
        "",
        "# Grassmann subspace under matched-strength surrogacy AFTER "
        "epileptic-zone exclusion",
        "",
        f"**Head.** R={n_surr} strength-preserving (4-cycle ±δ, "
        f"SWAP_FACTOR={swap_factor}) surrogate Laplacians on the "
        f"epi-EXCLUDED FC W_red = W[non_epi, non_epi] per (patient, "
        f"band, phase). γ_h-priority: did the audit_66 trace at k=19..27 "
        f"survive after removing epileptic-zone contacts? Extended to all "
        f"6 bands for full sensitivity surface.",
        "",
        "## Epi counts per patient",
        "",
        "| patient | N_full | N_epi | N_reduced | epi_fraction |",
        "|---|---|---|---|---|",
    ]
    for _, r in epi_counts.iterrows():
        lines.append(
            f"| {r['patient']} | {int(r['N_full'])} | {int(r['N_epi'])} "
            f"| {int(r['N_reduced'])} | {r['epi_fraction']:.3f} |"
        )
    lines.extend([
        "",
        "## Per-band cohort surface (epi-excluded)",
        "",
        "| band | k_min | k_max | k_separated | k_int_p_lt_05 | k_also_neg | "
        "n_k_sep | longest_run_sep |",
        "|---|---|---|---|---|---|---|---|",
    ])
    for band in _band_order_present(per_pat):
        sub = cohort[cohort.band == band]
        if sub.empty:
            continue
        ks = sub.k.values
        verdicts = sub.verdict.values
        ps = sub.paired_wilcoxon_p.values
        sep_mask = (verdicts == "separated")
        int_p_mask = (verdicts == "intermediate") & (ps < 0.05)
        neg_mask = (verdicts == "also_negative")
        n_sep = int(sep_mask.sum())
        longest, cur = 0, 0
        for v in sep_mask:
            cur = cur + 1 if v else 0
            if cur > longest:
                longest = cur
        k_sep_str = ",".join(str(int(k)) for k, v in zip(ks, sep_mask) if v)
        k_int_p_str = ",".join(str(int(k)) for k, v in zip(ks, int_p_mask) if v)
        k_neg_str = ",".join(str(int(k)) for k, v in zip(ks, neg_mask) if v)
        lines.append(
            f"| {band} | {int(ks.min())} | {int(ks.max())} "
            f"| {k_sep_str or '—'} "
            f"| {k_int_p_str or '—'} "
            f"| {k_neg_str or '—'} "
            f"| {n_sep}/{len(ks)} | {longest} |"
        )
    lines.extend([
        "",
        "## Sensitivity vs audit_66 (full FC)",
        "",
        "| band | n_k | persist | weaken | emerge | absent | longest_persist |",
        "|---|---|---|---|---|---|---|",
    ])
    for band in _band_order_present(per_pat):
        sub = sens[sens.band == band].sort_values("k")
        if sub.empty:
            continue
        flags = sub.sensitivity_flag.values
        n_per = int((flags == "persist").sum())
        n_wk = int((flags == "weaken").sum())
        n_em = int((flags == "emerge").sum())
        n_ab = int((flags == "absent").sum())
        longest, cur = 0, 0
        for f in flags:
            cur = cur + 1 if f == "persist" else 0
            if cur > longest:
                longest = cur
        lines.append(
            f"| {band} | {len(sub)} | {n_per} | {n_wk} | {n_em} | {n_ab} "
            f"| {longest} |"
        )
    lines.extend([
        "",
        "## γ_h zoom — k=19..27 sensitivity (audit_66 strict-separated window)",
        "",
        "| k | full_p | full_n_below | epiX_p | epiX_n_below | flag |",
        "|---|---|---|---|---|---|",
    ])
    gh = sens[(sens.band == "high_gamma") & (sens.k.between(19, 27))]
    for _, r in gh.iterrows():
        try:
            fp_str = f"{float(r['full_paired_wilcoxon_p']):.4f}"
            ep_str = f"{float(r['epiX_paired_wilcoxon_p']):.4f}"
        except Exception:
            fp_str = "—"; ep_str = "—"
        lines.append(
            f"| {int(r['k'])} | {fp_str} | {r['full_n_below']} "
            f"| {ep_str} | {r['epiX_n_below']} "
            f"| **{r['sensitivity_flag']}** |"
        )

    lines.extend([
        "",
        "## β zoom — k=27..55 sensitivity (audit_66 manuscript longest-run)",
        "",
        "| k | full_p | full_n_below | epiX_p | epiX_n_below | flag |",
        "|---|---|---|---|---|---|",
    ])
    bb = sens[(sens.band == "beta") & (sens.k.between(27, 55))]
    for _, r in bb.iterrows():
        try:
            fp_str = f"{float(r['full_paired_wilcoxon_p']):.4f}"
            ep_str = f"{float(r['epiX_paired_wilcoxon_p']):.4f}"
        except Exception:
            fp_str = "—"; ep_str = "—"
        lines.append(
            f"| {int(r['k'])} | {fp_str} | {r['full_n_below']} "
            f"| {ep_str} | {r['epiX_n_below']} "
            f"| **{r['sensitivity_flag']}** |"
        )

    lines.extend([
        "",
        "## Verdict labels",
        "",
        "- **persist**: full_p < 0.05 AND epiX_p < 0.05 — trace survives epi-exclusion.",
        "- **weaken**: full_p < 0.05 AND epiX_p ≥ 0.05 — trace requires epi-contact participation.",
        "- **emerge**: full_p ≥ 0.05 AND epiX_p < 0.05 — non-epi subset reveals trace masked by epi contacts (rare).",
        "- **absent**: full_p ≥ 0.05 AND epiX_p ≥ 0.05 — no trace in either.",
        "",
        "## Test choice",
        "",
        "Identical to audit_66 except W is replaced by W_red = W[non_epi, non_epi] "
        "before LRG; surrogate ensemble is generated on W_red. Epi mask from "
        "`load_epileptic_nodes` (red-font label in implant XLSX) joined with "
        "`load_channel_labels` to FC index order. K_GRID = "
        f"{K_GRID[0]}..{K_GRID[-1]} (capped by min N_reduced=89, Pat_13). "
        "Pat_15 has 0 epi → identity transform; included as anchor.",
        "",
        "## Provenance",
        "",
        f"- N_surrogates = {n_surr} per cell",
        f"- n_swaps_per_surrogate = SWAP_FACTOR ({swap_factor}) × N_red(N_red−1)/2",
        "- Cohort: " + ", ".join(COHORT),
        "- Bands: " + ", ".join(_band_order_present(per_pat)),
        "- Phases: " + ", ".join(PHASES),
        f"- K_GRID = {K_GRID[0]}..{K_GRID[-1]} (n={len(K_GRID)})",
        "- FC method: imcoh_abs",
        "- Cache (epi-excluded): "
        "`data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/`",
        f"- Wall-clock runtime: {runtime_s:.1f} s",
        "- Build script: `scripts/01_compute/audit/audit_67_grassmann_epi_exclusion.py`",
        "",
        "## Files",
        "",
        "- `per_patient_per_band_per_k.csv` — observed + surrogate stats",
        "- `cohort_summary.csv` — per-(band, k) cohort verdict",
        "- `sensitivity.csv` — audit_66 vs audit_67 join + flag",
        "- `epi_counts.csv` — N_full / N_epi / N_reduced per patient",
        "- `figures/grassmann_cohort_distribution_epiX.pdf`",
        "- `figures/grassmann_kspan_verdict_epiX.pdf`",
        "- `figures/grassmann_kspan_sensitivity_vs_audit66.pdf`",
        "",
    ])
    (OUT / "README.md").write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=ALL_BANDS)
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--n-surrogates", type=int, default=N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=SWAP_FACTOR)
    ap.add_argument("--seed", type=int, default=20260514)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    bands = args.bands
    patients = args.patients
    n_surr = args.n_surrogates
    swap_factor = args.swap_factor

    print("[audit_67] pre-flight: ensure half FCs cached")
    for pat in patients:
        ensure_half_fcs(pat, bands)

    print("[audit_67] writing epi counts")
    epi_counts = write_epi_counts()

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
            summary_bits = []
            for k_ref in K_REFERENCE:
                row = next((r for r in cell["rows"] if r["k"] == k_ref), None)
                if row is not None:
                    summary_bits.append(
                        f"k{k_ref}:T={row['obs_T_G']:+.2f} z={row['obs_z']:+.2f}"
                    )
            print(f"[audit_67] {pat}/{band}: {' | '.join(summary_bits)} ({dt:.1f}s)")

    runtime = time.time() - t0
    n_cells = len({(r['patient'], r['band']) for r in all_rows})
    print(f"[audit_67] total: {runtime:.1f}s for {n_cells} cells")

    per_pat = pd.DataFrame(all_rows)
    per_pat.to_csv(OUT / "per_patient_per_band_per_k.csv", index=False)

    cohort = cohort_summary(per_pat)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)

    sens = sensitivity_table(cohort)
    sens.to_csv(OUT / "sensitivity.csv", index=False)

    make_cohort_figure(per_pat, cohort,
                        FIG / "grassmann_cohort_distribution_epiX.pdf")
    make_kspan_verdict_figure(cohort,
                               FIG / "grassmann_kspan_verdict_epiX.pdf")
    make_sensitivity_figure(sens,
                             FIG / "grassmann_kspan_sensitivity_vs_audit66.pdf")

    write_readme(cohort, per_pat, sens, epi_counts,
                 n_surr, swap_factor, runtime)
    print(f"[audit_67] outputs at {OUT}")


if __name__ == "__main__":
    main()
