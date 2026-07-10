#!/usr/bin/env python3
"""Audit 165 — Grassmann subspace analog of the inference-specific β trace.

Second, structurally-independent witness for the flagship inference-specific
consolidation currently carried by the SINGLE cophenetic probe
``T_infspec·e = ρ_sym(f, p·e)`` (β Wilcoxon p = 0.0098, β-only; audit_152).
Extends the whole-task Grassmann probe (audit_66/70) from the collapsed task
to the four-phase arc ``rest_pre → task_learn → task_test → rest_post``.

Scope report (statistic, caveats, gate, pseudocode):
    .agents/guides/task-persistence-investigation/2026-07-09_grassmann-inference-arc.md

5-point critical preamble
-------------------------
1. CLAIM. The leading-k Laplacian eigenmode subspace of ``rest_post`` sits
   preferentially closer to the inference subspace (``task_test``) than to the
   encoding subspace (``task_learn``), beyond the rest baseline, in β and not in
   the null bands — corroborating the cophenetic dissociation (α consolidates
   encoding; β consolidates encoding + inference-specific).
2. NULL. Matched-strength 4-cycle ±δ (R=200, swap 20, seed 20260511), the exact
   cached ensemble that backs audit_66 / the cophenetic arc; independent
   per-phase rewiring, paired by realization index.
3. STRONGEST ALTERNATIVES. (i) encoding leakage TL↔TT; (ii) distance-difference
   noise (scalar-per-patient, low power); (iii) task_test 1.4–2.5× longer than
   task_learn biasing d_G(·,TT) vs d_G(·,TL); (iv) strength drift; (v) eigen
   near-degeneracy at the k-cut.
4. NULL'S REACH. Matched-strength kills (iv) only. (i) cancels by construction —
   Δ differences closeness-to-TT against closeness-to-TL; (ii) within-rest floor
   d_G(preA,preB) + cluster-mass over 111 k; (iii) the MANDATORY length-matched
   recompute of U_k^TT from head-truncated task_test; (v) k-grid sweep.
5. FALSIFICATION. Requires T_G^infspec(k) > 0 over a contiguous k-band with
   cluster_p_mass < 0.05, AND β-selectivity, AND survival of length-matching.
   Falsified if RP is equidistant from TL/TT beyond baseline, or the β effect
   collapses under length-matching, or it fires in θ (the absent band).
   Non-corroboration is AMBIGUOUS (probe blind-spot / power deficit), never a
   refutation of the cophenetic result; corroboration is strong.

The measure (per patient, per band, per k ∈ {2..112})
-----------------------------------------------------
Phases preA=rest_pre_A, preB=rest_pre_B, TL=task_learn, TT=task_test,
RP=rest_post. Chordal d_G(x,y;k) = sqrt(k − ‖U_x^{(k)T} U_y^{(k)}‖_F²), U from
combinatorial L=D−W, eigh ascending, drop trivial column 0, columns [1:k+1].

    Δ(r;k)          = d_G(r,TL;k) − d_G(r,TT;k)
    T_G^infspec(k)  = Δ(RP;k) − 0.5·[Δ(preA;k)+Δ(preB;k)]      HEADLINE
    T_G^enc(k)      = d_G(preA,TL;k) − d_G(TL,RP;k)            encoding echo
    T_G^onl(k)      = d_G(preA,TT;k) − d_G(TT,RP;k)            = whole-task T_G
    floor(k)        = d_G(preA,preB;k)                         within-rest floor
Sign convention: positive = trace (one-sided Wilcoxon 'greater').

Length-matched variant: recompute U_k^TT from ``_tt_fc_lenmatched`` (head-
truncated to n_TL, unchanged nperseg) and use cached ``task_test_lm_head``
surrogate eigs; everything else identical → T_G^infspec_lenmatched.

Cross-check gate (asserted before any new number is trusted): T_G^onl(k) here
MUST equal audit_66's whole-task obs_T_G per (patient, band, k) to < 1e-6.

Cohort aggregation (reused verbatim from audit_70): per-k paired one-sided
Wilcoxon (obs > surr-mean) → p_k; cluster_mass = Σ_{p_k<0.05} −log10(p_k);
empirical cluster_p_mass via the phantom-surrogate null (each surrogate as obs,
mean of the rest as reference); normalized T_G* = mass/(111·log10(R+1));
Decision-8 gate strong<0.01 / weak / none≥0.05; LOO max p (strong-tier
precondition; computed only where the gate clears — its sole role).

Outputs: data/audit/grassmann_inference_arc/{per_patient.csv,
cohort_summary.csv, README.md}. No figures.
"""
from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (  # noqa: E402
    BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH  # noqa: E402
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch  # noqa: E402
from lrg_eegfc.utils.io.patient import load_timeseries  # noqa: E402
from lrg_eegfc.utils.metrics.spectral import (  # noqa: E402
    chordal_distances_for_k_grid, laplacian_eig,
)
from lrg_eegfc.utils.surrogate.matched_strength import (  # noqa: E402
    load_or_compute_surrogate_eigs, surrogate_cache_path,
)

# Reuse the audit_70 cohort-aggregation helpers VERBATIM (no fork).
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_70_grassmann_cluster_extent import (  # type: ignore  # noqa: E402
    cluster_mass, longest_run_below, wilcoxon_per_k_greater,
)
# Reuse the audit_66 FC loaders (same construction that produced the CSV the
# cross-check validates against), the audit_113 windowing, and the audit_113b
# length-matched task_test FC (spot-checked below against our band-hoisted core).
from audit_66_grassmann_matched_strength_surrogate import (  # type: ignore  # noqa: E402
    ensure_half_fcs, load_phase_fc,
)
from audit_113_inference_length_control import _window  # type: ignore  # noqa: E402
from audit_113b_inference_lenmatched_null import (  # type: ignore  # noqa: E402
    PLACEMENT, _tt_fc_lenmatched,
)

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import _band_abs  # type: ignore  # noqa: E402

# --- configuration (frozen, matches audit_66 / audit_70) -------------------
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
K_GRID = list(range(2, 113))          # 111 values
R = 200
SWAP = 20
SEED = 20260511
ALPHA_K = 0.05
FC = "imcoh_abs"

# phase → surrogate cache phase-string
OBS_PHASES = ["rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post"]
LM_PHASE = "task_test_lm_head"        # length-matched TT surrogate ensemble

# gated functionals (full matched-strength cluster-mass test)
GATED = ["enc", "onl", "infspec", "infspec_lm"]
# LOO (strong-tier precondition) only for the headline + its length-matched twin
LOO_FUNCS = ["infspec", "infspec_lm"]

OBS_CSV = (ROOT / "data/audit/grassmann_matched_strength_surrogate"
           / "per_patient_per_band_per_k.csv")
OUT = ROOT / "data/audit/grassmann_inference_arc"
# One-time cache of the observed length-matched task_test eigvecs (the Welch
# imcoh spectrum is band-independent, so it is computed ONCE per patient and
# band-sliced — hoisting the band loop out of the 39-s Welch).
LM_EIG_CACHE = CACHE_ROOT / "lenmatched_tt_obs_eig"


def _clean_fc(W: np.ndarray) -> np.ndarray:
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def _tt_coh_lenmatched(pat: str):
    """Band-independent freq-resolved imcoh of the head-truncated task_test.

    The band-hoisted core of ``audit_113b._tt_fc_lenmatched``: identical
    windowing (`_window`, PLACEMENT='head', n = n_task_learn) and identical
    `compute_msc_welch(metric='imcoh')`, but returns the full ``(freqs, Coh)``
    so all six bands slice one Welch. Spot-checked bit-for-bit against
    ``_tt_fc_lenmatched`` in the preflight."""
    X = np.asarray(load_timeseries(pat, "task_test", SEEG_DATAPATH), float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    n_TL = int(max(load_timeseries(pat, "task_learn", SEEG_DATAPATH).shape))
    fs = float(FS_OVERRIDES.get(pat, 2048.0))
    Xw = np.ascontiguousarray(_window(X, n_TL, PLACEMENT))
    freqs, Coh = compute_msc_welch(Xw, fs, nperseg=nperseg_for_fs(fs),
                                   metric="imcoh")
    del X, Xw
    return freqs, Coh


def _lm_eig_path(pat: str, band: str) -> Path:
    return LM_EIG_CACHE / pat / f"{band}_{LM_PHASE}_obs_eigvecs.npz"


def ensure_lm_eigvecs(pat: str, bands: list[str],
                      spot_check_band: str | None = None) -> None:
    """Compute + cache observed length-matched task_test eigvecs for all bands.

    ``spot_check_band`` (once, at startup) recomputes ``_tt_fc_lenmatched`` for
    that single band and asserts our band-hoisted FC is bit-identical — one band
    proves the hoisting (only the ``_band_abs`` slice differs per band)."""
    paths = {b: _lm_eig_path(pat, b) for b in bands}
    if all(p.exists() for p in paths.values()):
        return
    freqs, Coh = _tt_coh_lenmatched(pat)
    for b in bands:
        flo, fhi = BRAIN_BANDS[b]
        W = _clean_fc(_band_abs(Coh, freqs, flo, fhi))
        if b == spot_check_band:
            W_ref = _tt_fc_lenmatched(pat, b)
            dev = float(np.max(np.abs(W - W_ref)))
            if not (dev < 1e-9):
                raise SystemExit(
                    f"lm FC spot-check FAILED {pat}/{b}: max|hoisted−ref|={dev:.2e}")
            print(f"    lm FC spot-check {pat}/{b}: max dev {dev:.2e} (OK)",
                  flush=True)
        evecs = laplacian_eig(W)[1]
        p = _lm_eig_path(pat, b)
        p.parent.mkdir(parents=True, exist_ok=True)
        np.savez(p, eigvecs=evecs.astype(np.float64))
    del Coh


def load_lm_eigvecs(pat: str, band: str) -> np.ndarray:
    with np.load(_lm_eig_path(pat, band)) as d:
        return d["eigvecs"].astype(np.float64)


# ---------------------------------------------------------------------------
# functionals from a dict of chordal k-grids per ordered phase-pair
# ---------------------------------------------------------------------------
def _functionals(d: dict[tuple[str, str], np.ndarray],
                 tt: str = "task_test") -> dict[str, np.ndarray]:
    """Build the four arc functionals (+floor) from chordal distances ``d``.

    ``d[(x, y)]`` is the chordal k-grid for the phase pair (symmetric, so only
    one ordering is stored). ``tt`` selects which key plays the inference
    (task_test) leg — ``"task_test"`` for the standard run, ``"task_test_lm"``
    for the length-matched run.
    """
    def g(x, y):
        return d[(x, y)] if (x, y) in d else d[(y, x)]

    d_RP_TL = g("rest_post", "task_learn")
    d_RP_TT = g("rest_post", tt)
    d_pA_TL = g("rest_pre_A", "task_learn")
    d_pA_TT = g("rest_pre_A", tt)
    d_pB_TL = g("rest_pre_B", "task_learn")
    d_pB_TT = g("rest_pre_B", tt)
    delta_RP = d_RP_TL - d_RP_TT
    delta_pA = d_pA_TL - d_pA_TT
    delta_pB = d_pB_TL - d_pB_TT
    return {
        "enc": d_pA_TL - d_RP_TL,                       # d(preA,TL) − d(TL,RP)
        "onl": d_pA_TT - d_RP_TT,                       # d(preA,TT) − d(TT,RP)
        "infspec": delta_RP - 0.5 * (delta_pA + delta_pB),
        "floor": g("rest_pre_A", "rest_pre_B"),
    }


# pairs required to build every functional (TT leg keyed generically)
def _pairs(tt: str) -> list[tuple[str, str]]:
    return [("rest_post", "task_learn"), ("rest_post", tt),
            ("rest_pre_A", "task_learn"), ("rest_pre_A", tt),
            ("rest_pre_B", "task_learn"), ("rest_pre_B", tt),
            ("rest_pre_A", "rest_pre_B")]


# ---------------------------------------------------------------------------
# cross-check gate (anti-hallucination) — do FIRST
# ---------------------------------------------------------------------------
def cross_check() -> None:
    print("=== CROSS-CHECK: T_G^onl(k) vs audit_66 whole-task obs_T_G ===",
          flush=True)
    csv = pd.read_csv(OBS_CSV)
    max_dev = 0.0
    worst = None
    for pat in PATIENTS:                    # full β cohort (cheap, thorough)
        evA = laplacian_eig(load_phase_fc(pat, "rest_pre_A", "beta"))[1]
        evT = laplacian_eig(load_phase_fc(pat, "task_test", "beta"))[1]
        evP = laplacian_eig(load_phase_fc(pat, "rest_post", "beta"))[1]
        d_preA_TT = chordal_distances_for_k_grid(evA, evT, K_GRID)
        d_TT_RP = chordal_distances_for_k_grid(evT, evP, K_GRID)
        my_onl = d_preA_TT - d_TT_RP
        sub = csv[(csv.patient == pat) & (csv.band == "beta")].sort_values("k")
        csv_obs = sub.set_index("k")["obs_T_G"].reindex(K_GRID).values
        dev = float(np.nanmax(np.abs(my_onl - csv_obs)))
        if dev > max_dev:
            max_dev, worst = dev, pat
    print(f"    max |T_G^onl − audit_66 obs_T_G| over β cohort = {max_dev:.3e} "
          f"(worst {worst})", flush=True)
    if not (max_dev < 1e-6):
        raise SystemExit(
            f"CROSS-CHECK FAILED: deviation {max_dev:.3e} >= 1e-6 — the chordal "
            f"/ topk_basis / laplacian_eig pipeline does not reproduce the "
            f"whole-task Grassmann probe. Refusing to report new numbers.")
    print("    PASS (< 1e-6) — arc functionals share the audit_66 geometry.\n",
          flush=True)


# ---------------------------------------------------------------------------
# per-cell: observed + surrogate functionals for one (patient, band)
# ---------------------------------------------------------------------------
def _preflight_caches(pat: str, band: str) -> None:
    need = [surrogate_cache_path(pat, band, ph, R, SWAP, SEED, FC)
            for ph in OBS_PHASES + [LM_PHASE]]
    missing = [p for p in need if not p.exists()]
    if missing:
        raise SystemExit(
            "Missing surrogate cache (refusing to recompute silently):\n  "
            + "\n  ".join(str(m) for m in missing))


def compute_cell(pat: str, band: str, rng: np.random.Generator
                 ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    _preflight_caches(pat, band)

    # --- observed eigvecs (fresh eigh on observed W; lm eigvecs from cache) ---
    W = {ph: load_phase_fc(pat, ph, band) for ph in OBS_PHASES}
    ev = {ph: laplacian_eig(W[ph])[1] for ph in OBS_PHASES}
    ev["task_test_lm"] = load_lm_eigvecs(pat, band)   # band-hoisted Welch cache
    if ev["task_test_lm"].shape[0] != W["task_test"].shape[0]:
        raise SystemExit(f"{pat}/{band}: lm task_test N mismatch")

    d_obs = {p: chordal_distances_for_k_grid(ev[p[0]], ev[p[1]], K_GRID)
             for p in _pairs("task_test")}
    d_obs_lm = {p: chordal_distances_for_k_grid(ev[p[0]], ev[p[1]], K_GRID)
                for p in _pairs("task_test_lm")}
    f_std = _functionals(d_obs, tt="task_test")
    f_lm = _functionals(d_obs_lm, tt="task_test_lm")
    obs = {"enc": f_std["enc"], "onl": f_std["onl"],
           "infspec": f_std["infspec"], "infspec_lm": f_lm["infspec"],
           "floor": f_std["floor"]}

    # --- surrogate eigvecs (cached matched-strength ensemble; lm uses its own
    #     cached ensemble — full-TT W passed only for the N shape check on hit) ---
    se: dict[str, np.ndarray] = {}
    for ph in OBS_PHASES:
        _, se[ph] = load_or_compute_surrogate_eigs(
            pat, band, ph, W[ph], R, SWAP, SEED, rng)
    _, se["task_test_lm"] = load_or_compute_surrogate_eigs(
        pat, band, LM_PHASE, W["task_test"], R, SWAP, SEED, rng)

    surr = {f: np.full((R, len(K_GRID)), np.nan) for f in GATED}
    pairs_std = _pairs("task_test")
    pairs_lm = _pairs("task_test_lm")
    for r in range(R):
        # nan eigvecs (rare failed-strength surrogate) propagate to nan
        # functionals; nanmean/nanmedian drop them downstream, per-functional.
        d_r = {p: chordal_distances_for_k_grid(se[p[0]][r], se[p[1]][r], K_GRID)
               for p in pairs_std}
        d_r_lm = {p: chordal_distances_for_k_grid(se[p[0]][r], se[p[1]][r], K_GRID)
                  for p in pairs_lm}
        fr = _functionals(d_r, tt="task_test")
        fr_lm = _functionals(d_r_lm, tt="task_test_lm")
        surr["enc"][r] = fr["enc"]
        surr["onl"][r] = fr["onl"]
        surr["infspec"][r] = fr["infspec"]
        surr["infspec_lm"][r] = fr_lm["infspec"]

    del W, ev, se
    gc.collect()
    return obs, surr


# ---------------------------------------------------------------------------
# aggregation (audit_70 verbatim) for one (functional, band)
# ---------------------------------------------------------------------------
def phantom_null(surr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Phantom-surrogate null (audit_70): each surrogate r is the phantom
    observation, reference = mean of the remaining R−1. Returns
    ``(null_mass[R], null_LR[R])``."""
    Rn = surr.shape[1]
    null_mass = np.zeros(Rn)
    null_LR = np.zeros(Rn, dtype=int)
    for r in range(Rn):
        phantom = surr[:, r, :]
        mask = np.ones(Rn, dtype=bool); mask[r] = False
        ref = np.nanmean(surr[:, mask, :], axis=1)
        pk = wilcoxon_per_k_greater(phantom, ref)
        null_mass[r] = cluster_mass(pk, ALPHA_K)
        null_LR[r] = longest_run_below(pk, ALPHA_K)
    return null_mass, null_LR


def aggregate(func: str, obs: np.ndarray, surr: np.ndarray) -> dict:
    surr_mean = np.nanmean(surr, axis=1)                # (P, K)
    p_k = wilcoxon_per_k_greater(obs, surr_mean)
    obs_mass = cluster_mass(p_k, ALPHA_K)
    obs_LR = longest_run_below(p_k, ALPHA_K)
    null_mass, null_LR = phantom_null(surr)
    cluster_p_mass = (1 + int(np.sum(null_mass >= obs_mass))) / (R + 1)
    cluster_p_LR = (1 + int(np.sum(null_LR >= obs_LR))) / (R + 1)
    t_g_star = obs_mass / (len(K_GRID) * np.log10(R + 1))

    verdict = ("strong" if cluster_p_mass < 0.01
               else "weak" if cluster_p_mass < 0.05 else "no_trace")

    # per-patient sign counts (obs functional median over k; and vs own surr)
    med_obs_pp = np.nanmedian(obs, axis=1)              # (P,)
    med_rel_pp = np.nanmedian(obs - surr_mean, axis=1)  # obs above own surr-mean
    n_pos_median = int(np.sum(med_obs_pp > 0))
    n_pos_vs_surr = int(np.sum(med_rel_pp > 0))

    return {
        "obs_cluster_mass": obs_mass, "obs_longest_run": obs_LR,
        "T_G_star_normalized": t_g_star,
        "cluster_p_mass": cluster_p_mass, "cluster_p_LR": cluster_p_LR,
        "verdict": verdict,
        "cohort_median_obs": float(np.median(med_obs_pp)),
        "cohort_median_surr": float(np.median(np.nanmedian(surr_mean, axis=1))),
        "n_pos_median": n_pos_median, "n_pos_vs_surr": n_pos_vs_surr,
        "_med_obs_pp": med_obs_pp, "_med_rel_pp": med_rel_pp,
    }


def loo_max_p(obs: np.ndarray, surr: np.ndarray) -> tuple[float, str]:
    """Leave-one-patient-out max cluster_p_mass (strong-tier precondition)."""
    P = obs.shape[0]
    loo_p = np.zeros(P)
    for drop in range(P):
        keep = np.ones(P, dtype=bool); keep[drop] = False
        o = obs[keep]; s = surr[keep]
        sm = np.nanmean(s, axis=1)
        pk = wilcoxon_per_k_greater(o, sm)
        m = cluster_mass(pk, ALPHA_K)
        nm = np.zeros(R)
        for r in range(R):
            ph = s[:, r, :]
            mask = np.ones(R, dtype=bool); mask[r] = False
            ref = np.nanmean(s[:, mask, :], axis=1)
            nm[r] = cluster_mass(wilcoxon_per_k_greater(ph, ref), ALPHA_K)
        loo_p[drop] = (1 + int(np.sum(nm >= m))) / (R + 1)
    j = int(np.argmax(loo_p))
    return float(loo_p[j]), PATIENTS[j]


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[audit_165] output → {OUT}", flush=True)
    print(f"[audit_165] pre-flight: ensure rest_pre half FCs cached", flush=True)
    for pat in PATIENTS:
        ensure_half_fcs(pat, BANDS)

    print(f"[audit_165] pre-flight: cache observed length-matched task_test "
          f"eigvecs (band-hoisted Welch, one-time)", flush=True)
    t_lm = time.time()
    for i, pat in enumerate(PATIENTS):
        tc = time.time()
        ensure_lm_eigvecs(pat, BANDS, spot_check_band=("beta" if i == 0 else None))
        el = time.time() - t_lm
        print(f"  lm-eig [{i+1}/{len(PATIENTS)}] {pat} {time.time()-tc:.1f}s "
              f"| elapsed {el:.0f}s", flush=True)

    cross_check()

    rng = np.random.default_rng(SEED)
    summary_rows: list[dict] = []
    per_patient_rows: list[dict] = []
    t_start = time.time()
    n_cells = len(BANDS) * len(PATIENTS)
    cell_i = 0

    for band in BANDS:
        tb = time.time()
        P, K = len(PATIENTS), len(K_GRID)
        obs_f = {f: np.full((P, K), np.nan) for f in GATED + ["floor"]}
        surr_f = {f: np.full((P, R, K), np.nan) for f in GATED}

        for pi, pat in enumerate(PATIENTS):
            tc = time.time()
            obs, surr = compute_cell(pat, band, rng)
            for f in GATED:
                obs_f[f][pi] = obs[f]
                surr_f[f][pi] = surr[f]
            obs_f["floor"][pi] = obs["floor"]
            cell_i += 1
            el = time.time() - t_start
            eta = el / cell_i * (n_cells - cell_i)
            print(f"  [{cell_i}/{n_cells}] {pat}/{band} "
                  f"cell {time.time()-tc:.2f}s | elapsed {el:.0f}s ETA {eta:.0f}s",
                  flush=True)

        # within-rest floor (descriptive)
        floor_med_pp = np.nanmedian(obs_f["floor"], axis=1)
        cohort_floor = float(np.median(floor_med_pp))

        print(f"  --- aggregating {band} ---", flush=True)
        for f in GATED:
            ta = time.time()
            agg = aggregate(f, obs_f[f], surr_f[f])
            lp, lpat = (np.nan, "")
            if f in LOO_FUNCS and agg["cluster_p_mass"] < 0.05:
                lp, lpat = loo_max_p(obs_f[f], surr_f[f])
            # cohort median |infspec| vs floor (abstain diagnostic)
            cohort_abs = float(np.median(np.abs(agg["_med_obs_pp"])))
            summary_rows.append({
                "band": band, "functional": f,
                "obs_cluster_mass": agg["obs_cluster_mass"],
                "obs_longest_run": agg["obs_longest_run"],
                "T_G_star_normalized": agg["T_G_star_normalized"],
                "cluster_p_mass": agg["cluster_p_mass"],
                "cluster_p_LR": agg["cluster_p_LR"],
                "verdict": agg["verdict"],
                "cluster_p_mass_loo_max": lp,
                "cluster_p_mass_loo_argmax_patient": lpat,
                "cohort_median_obs": agg["cohort_median_obs"],
                "cohort_median_surr": agg["cohort_median_surr"],
                "n_pos_median_of10": agg["n_pos_median"],
                "n_pos_vs_surr_of10": agg["n_pos_vs_surr"],
                "cohort_median_abs_obs": cohort_abs,
                "within_rest_floor_median": cohort_floor,
            })
            for pi, pat in enumerate(PATIENTS):
                per_patient_rows.append({
                    "patient": pat, "band": band, "functional": f,
                    "median_obs_over_k": float(agg["_med_obs_pp"][pi]),
                    "median_obs_minus_surr_over_k": float(agg["_med_rel_pp"][pi]),
                    "floor_median_over_k": float(floor_med_pp[pi]),
                })
            loo_s = f" LOOmax={lp:.4f}({lpat})" if f in LOO_FUNCS and lp == lp else ""
            print(f"    {f:11s} mass={agg['obs_cluster_mass']:7.2f} "
                  f"T_G*={agg['T_G_star_normalized']:.4f} "
                  f"p_mass={agg['cluster_p_mass']:.4f} "
                  f"[{agg['verdict']}] "
                  f"pos_med={agg['n_pos_median']}/10 "
                  f"pos_vs_surr={agg['n_pos_vs_surr']}/10{loo_s} "
                  f"({time.time()-ta:.1f}s)", flush=True)

        print(f"  === {band} done ({time.time()-tb:.0f}s, "
              f"floor med {cohort_floor:+.3f}) ===\n", flush=True)

    pd.DataFrame(summary_rows).to_csv(OUT / "cohort_summary.csv", index=False)
    pd.DataFrame(per_patient_rows).to_csv(OUT / "per_patient.csv", index=False)
    write_readme(pd.DataFrame(summary_rows), time.time() - t_start)
    print(f"[audit_165] total {(time.time()-t_start)/60:.1f} min", flush=True)
    print(f"[audit_165] wrote {OUT/'cohort_summary.csv'}", flush=True)
    print(f"[audit_165] wrote {OUT/'per_patient.csv'}", flush=True)


def write_readme(summ: pd.DataFrame, runtime_s: float) -> None:
    def row(band, f):
        r = summ[(summ.band == band) & (summ.functional == f)]
        return r.iloc[0] if len(r) else None
    lines = [
        "---",
        "name: grassmann-inference-arc",
        "era: IMCOH_ABS_COHORT_N10",
        "status: complete",
        "date: 2026-07-09",
        "scope: .agents/guides/task-persistence-investigation/2026-07-09_grassmann-inference-arc.md",
        "companion: ../grassmann_cluster_extent/  ../grassmann_matched_strength_surrogate/",
        "---",
        "",
        "# Grassmann subspace arc — second witness for the inference-specific β trace",
        "",
        "Subspace-geometry sibling of the cophenetic `T_infspec·e = ρ_sym(f, p·e)` "
        "(β Wilcoxon p=0.0098, β-only; audit_152). Extends the whole-task Grassmann "
        "probe (audit_66/70) to the four-phase arc. Matched-strength R=200 null "
        "(seed 20260511), reusing the cached ensemble byte-for-byte. Cross-check: "
        "`T_G^onl` reproduces audit_66 whole-task `obs_T_G` to < 1e-6.",
        "",
        "## Headline — T_G^infspec = Δ(RP) − ½[Δ(preA)+Δ(preB)], Δ(r)=d_G(r,TL)−d_G(r,TT)",
        "",
        "| band | cluster_p_mass | T_G* | verdict | LOO max p | pos_median/10 | pos_vs_surr/10 |",
        "|---|---|---|---|---|---|---|",
    ]
    for b in BANDS:
        r = row(b, "infspec")
        if r is None:
            continue
        lp = f"{r.cluster_p_mass_loo_max:.4f}" if r.cluster_p_mass_loo_max == r.cluster_p_mass_loo_max else "—"
        lines.append(
            f"| {b} | {r.cluster_p_mass:.4f} | {r.T_G_star_normalized:.4f} "
            f"| {r.verdict} | {lp} | {int(r.n_pos_median_of10)} "
            f"| {int(r.n_pos_vs_surr_of10)} |")
    lines += ["", "## Length-matched headline — T_G^infspec_lenmatched (head-truncated TT)",
              "",
              "| band | cluster_p_mass | T_G* | verdict | LOO max p |",
              "|---|---|---|---|---|"]
    for b in BANDS:
        r = row(b, "infspec_lm")
        if r is None:
            continue
        lp = f"{r.cluster_p_mass_loo_max:.4f}" if r.cluster_p_mass_loo_max == r.cluster_p_mass_loo_max else "—"
        lines.append(f"| {b} | {r.cluster_p_mass:.4f} | {r.T_G_star_normalized:.4f} "
                     f"| {r.verdict} | {lp} |")
    lines += ["", "## Companions — encoding echo & inference-online (whole-task T_G)",
              "",
              "| band | enc p_mass | enc verdict | onl p_mass | onl verdict |",
              "|---|---|---|---|---|"]
    for b in BANDS:
        re_, ro = row(b, "enc"), row(b, "onl")
        if re_ is None:
            continue
        lines.append(f"| {b} | {re_.cluster_p_mass:.4f} | {re_.verdict} "
                     f"| {ro.cluster_p_mass:.4f} | {ro.verdict} |")
    lines += [
        "",
        "## Gate (Decision-8, mass-only)",
        "- strong ⇔ cluster_p_mass < 0.01; weak ⇔ 0.01 ≤ p < 0.05; no_trace ⇔ p ≥ 0.05.",
        "- LOO max p (strong-tier precondition) computed only where the gate clears "
        "(p_mass < 0.05); its sole role is guarding a strong verdict against "
        "single-patient dependence.",
        "- `floor` = within-rest d_G(preA,preB), descriptive noise band "
        "(`within_rest_floor_median` column).",
        "",
        "## Parameters",
        f"- k-grid 2..112 ({len(K_GRID)} values); R={R}; SWAP={SWAP}; seed={SEED}; FC={FC}.",
        f"- length-matched TT: head-truncation to n(task_learn), PLACEMENT='{PLACEMENT}'.",
        f"- cohort: {', '.join(PATIENTS)} (n={len(PATIENTS)}).",
        f"- runtime {runtime_s/60:.1f} min.",
        "",
        "## Files",
        "- `cohort_summary.csv` — per (band, functional) mass / T_G* / cluster_p_mass "
        "/ cluster_p_LR / verdict / LOO / per-patient sign counts / floor.",
        "- `per_patient.csv` — per (patient, band, functional) median-over-k of obs, "
        "obs−surr, and floor.",
        "",
    ]
    (OUT / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
