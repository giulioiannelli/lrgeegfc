#!/usr/bin/env python3
"""audit_156 — epi/WM pair-class tissue stratification under the rho_sym estimator.

Migrates the two ρ_split-era pair-class results that R1.4/R1.6 cite to the
arbitrary-half-free rho_sym estimator, so no ρ_split number lands in the paper:

  * matched-strength pair-class verdict  (ports audit_77 epi + audit_83 WM)
  * random same-size pair-subset null p_pair  (ports audit_85 --mode pairclass)

for the pair classes  epi: {nonepi_nonepi, cross, epi_epi}  and
                      wm : {gray_gray,     cross, wm_wm}.

The ONLY change vs the ρ_split originals is the per-CLASS statistic: instead of a
single-arm Spearman it becomes the SYMMETRIC average of both split-half arm
assignments (mirrors rho_sym = 1/2[rho_AB + rho_BA]):

    rho_sym(pk) = 1/2 [ spearman((D_task-D_preA)[pk], (D_post-D_preB)[pk])
                      + spearman((D_task-D_preB)[pk], (D_post-D_preA)[pk]) ]

applied IDENTICALLY to the observed trace, every cached matched-strength surrogate
realization (canonical seed-20260511 full-graph ensemble, REUSED — no
regeneration), and every random-pair draw. Old audits 77/83/85 are NOT edited;
new script, new output dir. arm1 (the A→task,B→rest assignment) is emitted too as
a bit-level cross-check that must reproduce the locked ρ_split CSVs.

5-point preamble
1. Claim: the tissue verdicts hold under rho_sym — gray↔gray β CONCENTRATED, the
   seizure core does NOT carry β (epi_epi weaken/depleted), epi_epi α RECRUITS the
   core (concentrated + strengthens), wm_wm β below-baseline-but-strength-clean.
2. Null: the verdicts were an artifact of the arbitrary A/B split-half arm choice.
3. Strongest alternative: (i) node strength (matched-strength surrogate controls it,
   estimator-orthogonal); (ii) pair count for small classes (pair-decimation null
   controls it). The estimator swap is orthogonal to both.
4. Cannot: reuses R=200 cached surrogate (same as locked audit_77/83). Pair-class
   null is global-rewiring, not within-class (inherited caveat). RAW substrate is
   BASELINE ONLY (framework purity) — cophenetic is the result.
5. Falsify: if arm1 fails to reproduce the ρ_split CSVs (pipeline bug), or if a
   cited verdict flips sign of conclusion under rho_sym.

Output: data/audit/epi_wm_stratified_rhosym/
    pairclass_rhosym_per_patient.csv
    pairclass_rhosym_cohort.csv
    README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from numba import njit, prange
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
    lrg_ultrametric_condensed,
)
from audit_77_epi_stratified_cophenetic import _surr_cophenetic_stack  # type: ignore
import _epi_stratify as es  # type: ignore
import _wm_stratify as ws  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_wm_stratified_rhosym"
SEED_PAIR = 20260612          # match audit_85 _pair_rng for arm1 cross-check
R_PAIR_DEFAULT = 1000
PHASES_4 = es.PHASES_4        # (rest_pre_A, rest_pre_B, task_test, rest_post)

# stratify registry: (helper module, pairclass configs, node-mask fn, strat index
# used by audit_85's _pair_rng salt so arm1 draws match bit-for-bit).
STRAT = {
    "wm": dict(mod=ws, configs=ws.PAIRCLASS_CONFIGS, mask=ws.wm_mask_for,
               pair=ws.pair_keep_for_config, strat_i=0),
    "epi": dict(mod=es, configs=es.PAIRCLASS_CONFIGS, mask=es.epi_mask_for,
                pair=es.pair_keep_for_config, strat_i=1),
}


# ---------------------------------------------------------------------------
# Spearman over rows (average-tie ranks, BIT-IDENTICAL to scipy.spearmanr).
# numba-parallel over rows: ~8x faster than scipy.rankdata(axis=1) (single-
# threaded C) on the (R, M) surrogate/pair-decimation batches; verified
# max|nb - scipy| = 0.0 on tied cophenetic-difference data. RNG-free → the
# JIT introduces no nondeterminism (feedback_numba_for_surrogates).
# ---------------------------------------------------------------------------
@njit(cache=False, fastmath=False)
def _rank_avg_1d(x, order):
    n = x.size; r = np.empty(n, np.float64); i = 0
    while i < n:
        j = i
        while j + 1 < n and x[order[j + 1]] == x[order[i]]:
            j += 1
        avg = 0.5 * (i + j) + 1.0            # average rank of the tie block (1-based)
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


@njit(parallel=True, cache=False, fastmath=False)
def _spearman_rows_nb(A, B):
    R = A.shape[0]; out = np.empty(R, np.float64)
    for t in prange(R):
        ra = _rank_avg_1d(A[t], np.argsort(A[t], kind="mergesort"))
        rb = _rank_avg_1d(B[t], np.argsort(B[t], kind="mergesort"))
        ra -= ra.mean(); rb -= rb.mean()
        num = (ra * rb).sum()
        den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
        out[t] = num / den if den > 0 else np.nan
    return out


def _spearman_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Spearman ρ per row of (R, M) arrays; average ties (== scipy)."""
    return _spearman_rows_nb(np.ascontiguousarray(A, np.float64),
                             np.ascontiguousarray(B, np.float64))


def _rho_sym_rows(xA, yB, xB, yA, sel) -> np.ndarray:
    """rho_sym per row over the pair subset `sel` (bool or index)."""
    s1 = _spearman_rows(xA[:, sel], yB[:, sel])
    s2 = _spearman_rows(xB[:, sel], yA[:, sel])
    return 0.5 * (s1 + s2)


def _rho_sym_one(xA, yB, xB, yA, sel) -> tuple[float, float]:
    """(rho_sym, arm1) for a single observed vector over subset `sel`."""
    a1, _ = spearmanr(xA[sel], yB[sel])
    a2, _ = spearmanr(xB[sel], yA[sel])
    return float(0.5 * (a1 + a2)), float(a1)


def _pair_rng(pat, band, strat_i, cfg_i):
    """Replicates audit_85._pair_rng so arm1 draws match the ρ_split CSV."""
    pat_n = int(pat.split("_")[-1])
    band_i = es.ALL_BANDS.index(band)
    return np.random.default_rng(
        np.random.SeedSequence([SEED_PAIR, strat_i, pat_n, band_i, cfg_i]))


# ---------------------------------------------------------------------------
# Per (patient, band): every pair class, both nulls
# ---------------------------------------------------------------------------
def per_patient_band(pat, band, r_pair, verbose):
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in PHASES_4}
    except Exception as e:
        if verbose:
            print(f"  SKIP {pat}/{band}: load failed: {e}")
        return [], {}
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        return [], {}

    # observed full-graph cophenetic + the four difference-arm vectors
    D = {ph: lrg_ultrametric_condensed(Ws[ph]) for ph in PHASES_4}
    xA = D["task_test"] - D["rest_pre_A"]
    yB = D["rest_post"] - D["rest_pre_B"]
    xB = D["task_test"] - D["rest_pre_B"]
    yA = D["rest_post"] - D["rest_pre_A"]
    n_pairs_full = xA.size

    # cached full-graph matched-strength surrogate cophenetic (R, n_pairs) x4 phases
    surr = {}
    for ph in PHASES_4:
        surr[ph] = _surr_cophenetic_stack(
            Ws[ph], es.surr_eig_path("full", pat, band, ph),
            es.cell_rng(pat, band, ph, "full"),
            es.N_SURROGATES, es.SWAP_FACTOR, verbose)
    sxA = surr["task_test"] - surr["rest_pre_A"]
    syB = surr["rest_post"] - surr["rest_pre_B"]
    sxB = surr["task_test"] - surr["rest_pre_B"]
    syA = surr["rest_post"] - surr["rest_pre_A"]
    good = np.all(np.isfinite(sxA), axis=1) & np.all(np.isfinite(syB), axis=1)

    rows, store = [], {}
    for stratify, cfg in STRAT.items():
        mask = cfg["mask"](pat, n_expected=N)
        for cfg_i, config in enumerate(cfg["configs"]):
            pk = cfg["pair"](mask, config)
            if pk is None:
                continue
            M = int(pk.sum())
            # --- observed rho_sym over the class ---
            obs_sym, obs_arm1 = _rho_sym_one(xA, yB, xB, yA, pk)
            # --- matched-strength null (surrogate rho_sym over the class) ---
            ss = _rho_sym_rows(sxA[good], syB[good], sxB[good], syA[good], pk)
            ss = ss[np.isfinite(ss)]
            surr_med = float(np.median(ss)) if ss.size else np.nan
            p_ms = float(np.mean(ss >= obs_sym)) if ss.size else np.nan
            # --- pair-decimation null (random M-pair subsets, matched rng) ---
            rng = _pair_rng(pat, band, cfg["strat_i"], cfg_i)
            idx = np.array([rng.choice(n_pairs_full, size=M, replace=False)
                            for _ in range(r_pair)])
            xAn, yBn = xA[idx], yB[idx]              # (R_pair, M)
            xBn, yAn = xB[idx], yA[idx]
            null = 0.5 * (_spearman_rows(xAn, yBn) + _spearman_rows(xBn, yAn))
            null = null[np.isfinite(null)]
            p_pair = float(np.mean(null >= obs_sym)) if null.size else np.nan
            rows.append({
                "patient": pat, "band": band, "stratify": stratify,
                "config": config, "M_class_pairs": M, "n_pairs_full": n_pairs_full,
                "obs_rho_sym": obs_sym, "obs_arm1_rho_split": obs_arm1,
                "ms_surr_median": surr_med, "ms_p_upper": p_ms,
                "pair_null_median": float(np.median(null)) if null.size else np.nan,
                "p_pair": p_pair, "n_surr": int(ss.size), "n_pair_draws": int(null.size),
            })
            store[(stratify, config)] = {
                "obs": obs_sym,
                "null_full": np.pad(null, (0, r_pair - null.size),
                                    constant_values=np.nan) if null.size else np.full(r_pair, np.nan),
            }
    return rows, store


# ---------------------------------------------------------------------------
# Cohort summary
# ---------------------------------------------------------------------------
def cohort_summary(per_pat, store):
    out = []
    for stratify, cfg in STRAT.items():
        for band in es.ALL_BANDS:
            for config in cfg["configs"]:
                d = per_pat[(per_pat.stratify == stratify)
                            & (per_pat.band == band)
                            & (per_pat.config == config)]
                if d.empty:
                    continue
                obs = d.obs_rho_sym.to_numpy(float)
                surr_med = d.ms_surr_median.to_numpy(float)
                n_above = int((d.ms_p_upper < 0.05).sum())
                v = es.cohort_verdict(obs, surr_med, n_above)
                # cohort pair-decimation: median-over-patients per draw
                cell = store.get((band, stratify, config), [])
                cohort_obs = float(np.median([c["obs"] for c in cell])) if cell else np.nan
                cohort_p_pair = cohort_null_med = np.nan
                if cell:
                    nm = np.vstack([c["null_full"] for c in cell])
                    cmn = np.nanmedian(nm, axis=0)
                    cmn = cmn[np.isfinite(cmn)]
                    if cmn.size and np.isfinite(cohort_obs):
                        cohort_p_pair = float(np.mean(cmn >= cohort_obs))
                        cohort_null_med = float(np.median(cmn))
                if np.isfinite(cohort_p_pair) and cohort_p_pair < 0.05:
                    pc_verdict = "concentrated"
                elif np.isfinite(cohort_p_pair) and cohort_p_pair > 0.95:
                    pc_verdict = "depleted"
                else:
                    pc_verdict = "generic_paircount"
                out.append({
                    "stratify": stratify, "band": band, "config": config,
                    "n_defined": v["n_defined"],
                    "cohort_obs_rho_sym": v["med_obs"],
                    "ms_surr_median": v["med_surr"],
                    "n_above_own_surrogate": f"{n_above}/{v['n_defined']}",
                    "ms_wilcoxon_p": v["wilcoxon_p"],
                    "ms_verdict": v["verdict"],
                    "cohort_obs_class": cohort_obs,
                    "cohort_pair_null_median": cohort_null_med,
                    "cohort_p_pair": cohort_p_pair,
                    "pairclass_verdict": pc_verdict,
                })
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--r-pair", type=int, default=R_PAIR_DEFAULT)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    _spearman_rows(np.random.default_rng(0).standard_normal((4, 16)),
                   np.random.default_rng(1).standard_normal((4, 16)))  # warm JIT
    print("[audit_156] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    all_rows, store = [], {}
    n_cells = len(args.bands) * len(args.patients)
    k = 0
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            rows, pat_store = per_patient_band(pat, band, args.r_pair, args.verbose)
            all_rows.extend(rows)
            for (stratify, config), s in pat_store.items():
                store.setdefault((band, stratify, config), []).append(s)
            k += 1
            el = time.time() - t0
            eta = el / k * (n_cells - k)
            bits = " ".join(
                f"{r['stratify'][0]}:{r['config'][:4]}={r['obs_rho_sym']:+.2f}"
                f"(ms{r['ms_p_upper']:.2f},pp{r['p_pair']:.2f})" for r in rows)
            print(f"[audit_156 {k}/{n_cells}] {pat}/{band}: {bits} "
                  f"({time.time()-tc:.1f}s, elapsed {el:.0f}s ETA {eta:.0f}s)",
                  flush=True)

    per_pat = pd.DataFrame(all_rows)
    per_pat.to_csv(OUT / "pairclass_rhosym_per_patient.csv", index=False)
    cohort = cohort_summary(per_pat, store)
    cohort.to_csv(OUT / "pairclass_rhosym_cohort.csv", index=False)
    runtime = time.time() - t0
    print(f"\n[audit_156] {len(per_pat)} rows, {runtime:.1f}s -> {OUT}")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(cohort[["stratify", "band", "config", "cohort_obs_rho_sym",
                      "ms_wilcoxon_p", "ms_verdict", "cohort_p_pair",
                      "pairclass_verdict"]].to_string(index=False))


if __name__ == "__main__":
    main()
