#!/usr/bin/env python3
"""Audit 85 — node-count + pair-count confound controls for tissue stratification.

Two complementary controls behind one CLI flag (`--mode`):

* ``--mode node`` (default) — **size-matched random-node-decimation control for
  SUBSET exclusion**, documented in full below. Resolves the node-count honesty
  flag on the subgraph configs (``exclude_wm`` / ``exclude_epi``).
* ``--mode pairclass`` — **random same-size PAIR-subset null for the PAIR-CLASS
  configs** (``gray_gray`` / ``cross`` / ``wm_wm``; ``nonepi_nonepi`` / ``cross``
  / ``epi_epi``). Pair-class restricts the *pair set* on the FULL-graph distances
  (no node removal → no node-count confound), but a small class (epi_epi ≈ 45
  pairs) may just be a noisy M-pair estimate of the whole-graph trace. This mode
  draws M random pairs from the full pair set R times and asks whether the
  observed class ρ_split is genuinely CONCENTRATED (upper tail), DEPLETED (lower
  tail), or indistinguishable from a random M-pair subset (``generic_paircount``
  → underpowered, not a localization). Orthogonal to matched-strength. Writes
  ``pairclass_decimation_{cohort,per_patient}.csv`` + ``README_pairclass.md``.
  See the ``MODE = pairclass`` section below.

Resolves the node-count honesty flag on BOTH suspect-tissue exclusion controls:
**C6 white-matter exclusion** (`--stratify wm`) and **C5 epileptic-zone exclusion**
(`--stratify epi`). Both `exclude_wm` and `exclude_epi` drop a chunk of nodes, so
the full-vs-exclude *change* in ρ_split mixes any tissue biology with a generic
node-reduction effect (a smaller graph gives a denser / cleaner LRG
renormalization, or a different-variance cohort Spearman, regardless of *which*
nodes go). A `strengthen` or `emerge` under exclusion is therefore ambiguous
between "the removed tissue was diluting the trace" and "removing any K nodes does
this". This audit isolates the node-count effect.

For each (patient, band, stratification): let K = #excluded-tissue nodes. Remove K
nodes **at random** R_dec times (the removed set need NOT be the tissue),
recompute the split-baseline ρ_split (cophenetic + raw) on each N−K-node subgraph,
and compare the observed `exclude_<tissue>` ρ_split to that random-decimation
distribution — per patient (`p_dec`) and at the cohort level (`cohort_p_dec` =
P(random cohort-median ρ_split ≥ the exclude cohort-median), the statistic that
matches the verdict gate).

Critical preamble (per CLAUDE.md rule, before any code)
=======================================================
(1) **Claim.** The `exclude_<tissue>` sharpening / emergence of the ρ_split trace
    is SPECIFIC to removing that tissue — not a generic consequence of shrinking
    the montage.
(2) **Null.** Size-matched random-node decimation: drop K = #tissue nodes chosen
    uniformly at random, recompute ρ_split, R_dec=200 draws → an empirical
    distribution of the trace statistic under random K-node loss. NOT a
    matched-strength null (that is audit_77/83's job, regenerated per submatrix);
    this is the ORTHOGONAL control they lack — it holds *node count* fixed and
    randomizes *which* nodes are removed.
(3) **Strongest plausible alternative the null should control for.** That the
    observed `exclude_<tissue>` ρ_split is whatever you would get by removing
    *any* K nodes — i.e. the effect is the graph getting smaller, with no role
    for tissue type.
(4) **Whether the null controls for it — by mechanism.** The random-decimation
    subgraphs have EXACTLY the same node count N−K as `exclude_<tissue>` and are
    run through the IDENTICAL pipeline (submatrix → LRG ultrametric / raw
    squareform → split-baseline Spearman). The ONLY difference between obs and
    null is *which* nodes were removed. Clean apples-to-apples isolation of the
    tissue-identity effect with node count held fixed. It does NOT control for
    per-node strength (matched-strength does that) — deliberately complementary:
    matched-strength says the trace is not a strength artifact; this says it is
    (or is not) a node-count artifact. What it CANNOT reject: a confound between
    tissue membership and node strength, so a "tissue_specific" verdict means
    "specific to the set of nodes that happen to be that tissue".
(5) **Falsification + limitations.** "exclude_<tissue> sharpening is
    tissue-specific" requires the observed cohort-median in the UPPER tail of the
    random distribution (`cohort_p_dec < 0.05`). "Generic node-count" → obs ≈
    random median (`p_dec` ≈ 0.5). "Tissue was carrying the trace" → obs in the
    LOWER tail (`p_dec` ≈ 1). Limitations: (i) random decimation can sever a whole
    probe / hemisphere — intentional (it is the generic-node-loss reference), so
    the null is broad; (ii) single tissue definition per stratification (WM =
    atlas-dominant `Wm`; epi = canonical red-font SOZ via `build_epi_masks`);
    (iii) patients with K=0 (e.g. Pat_15 has 0 epi nodes) are skipped — exclusion
    is identity there, so the epi cohort is effectively n=9.

Outputs (per stratification)
----------------------------
``data/audit/{wm,epi}_stratified/``
    decimation_control_per_patient.csv   one row per pat×band×substrate
    decimation_control_cohort.csv        per (band, substrate) cohort verdict
    README_decimation.md
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

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

# Reuse the proven pipeline (no forks): audit_63 cophenetic + half-FC + ρ_split,
# audit_83 raw condensed, _wm_stratify / _epi_stratify masks + vocabulary.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
    lrg_ultrametric_condensed,
    rho_split_from_phases,
)
from audit_83_wm_stratified_cophenetic import raw_condensed  # type: ignore
import _wm_stratify as ws  # type: ignore
import _epi_stratify as es  # type: ignore


SUBSTRATES = ("cophenetic", "raw")
SEED_DEC = 20260610  # dedicated node-decimation seed (distinct from audit_77/83 ensembles)
SEED_PAIR = 20260612  # dedicated PAIR-decimation seed (mode=pairclass)
R_DEC_DEFAULT = 200
R_PAIR_DEFAULT = 1000  # pair draws are cheap (re-index + Spearman; no eig, no surrogate)

# Stratification registry — the ONLY difference between WM-X and epi-X is the
# K-defining tissue mask, the pair-class membership rule, the output directory,
# and labels. Everything statistical is shared. COHORT / ALL_BANDS / PHASES_4 /
# MIN_NODES are identical in both _wm_stratify and _epi_stratify, so we read them
# from `ws`.
STRAT = {
    "wm": dict(
        mask_for=ws.wm_mask_for, exclude_cfg="exclude_wm",
        pair_for=ws.pair_keep_for_config, pairclass_configs=ws.PAIRCLASS_CONFIGS,
        out=ROOT / "data" / "audit" / "wm_stratified",
        tissue="white matter (dominant Desikan-Killiany tissue == 'Wm', atlas argmax)",
        companions="audit_83 (cophenetic+raw), audit_84 (Grassmann), audit_86 (gate)",
    ),
    "epi": dict(
        mask_for=es.epi_mask_for, exclude_cfg="exclude_epi",
        pair_for=es.pair_keep_for_config, pairclass_configs=es.PAIRCLASS_CONFIGS,
        out=ROOT / "data" / "audit" / "epi_stratified",
        tissue="epileptic zone (canonical red-font SOZ via build_epi_masks)",
        companions="audit_77 (cophenetic+raw), audit_78 (Grassmann), audit_67 (Grassmann gate)",
    ),
}


def _dec_rng(pat: str, band: str, stratify: str) -> np.random.Generator:
    """Deterministic per-cell generator for the random-node draws.

    WM preserves its ORIGINAL seed entropy (no stratify salt) so the
    already-cascaded WM numbers stay bit-identical; epi gets a distinct salt."""
    pat_n = int(pat.split("_")[-1])
    band_i = ws.ALL_BANDS.index(band)
    if stratify == "wm":
        entropy = [SEED_DEC, pat_n, band_i]          # original WM stream — frozen
    else:
        entropy = [SEED_DEC, 991, pat_n, band_i]     # epi salt = 991
    return np.random.default_rng(np.random.SeedSequence(entropy))


def _pair_rng(pat: str, band: str, stratify: str, config: str
              ) -> np.random.Generator:
    """Deterministic per-(cell, config) generator for the random-PAIR-subset
    draws (mode=pairclass). Salted distinctly from the node-decimation stream."""
    pat_n = int(pat.split("_")[-1])
    band_i = ws.ALL_BANDS.index(band)
    strat_i = 0 if stratify == "wm" else 1
    cfg_i = list(STRAT[stratify]["pairclass_configs"]).index(config)
    return np.random.default_rng(
        np.random.SeedSequence([SEED_PAIR, strat_i, pat_n, band_i, cfg_i]))


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman ρ — same call audit_77/83 use to produce the observed pair-class
    ρ_split, so obs recomputed here matches the cohort CSVs bit-for-bit."""
    rho, _ = spearmanr(a, b)
    return float(rho)


def _rho_splits_for_keep(Ws: dict[str, np.ndarray], keep: np.ndarray
                         ) -> tuple[float, float]:
    """(cophenetic ρ_split, raw ρ_split) on the node submatrix W[keep, keep].

    Both substrates come off the one submatrix extraction; raw needs no eigh."""
    idx = np.ix_(keep, keep)
    Wk = {ph: np.ascontiguousarray(Ws[ph][idx]) for ph in ws.PHASES_4}
    Dc = {ph: lrg_ultrametric_condensed(Wk[ph]) for ph in ws.PHASES_4}
    Dr = {ph: raw_condensed(Wk[ph]) for ph in ws.PHASES_4}
    rc = rho_split_from_phases(Dc["rest_pre_A"], Dc["rest_pre_B"],
                               Dc["task_test"], Dc["rest_post"])
    rr = rho_split_from_phases(Dr["rest_pre_A"], Dr["rest_pre_B"],
                               Dr["task_test"], Dr["rest_post"])
    return rc, rr


def per_patient_band(pat: str, band: str, r_dec: int, stratify: str,
                     verbose: bool) -> tuple[list[dict], dict]:
    """Return (per-patient rows, store) where store[substrate] =
    {"obs": float, "dec_full": np.ndarray(r_dec) with NaN preserved} so the
    cohort-level decimation null can be assembled across patients."""
    mask_for = STRAT[stratify]["mask_for"]
    exclude_cfg = STRAT[stratify]["exclude_cfg"]
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in ws.PHASES_4}
    except Exception as e:
        if verbose:
            print(f"[audit_85:{stratify}] SKIP {pat}/{band}: load failed: {e}")
        return [], {}
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_85:{stratify}] SKIP {pat}/{band}: phase shape mismatch")
        return [], {}

    mask = mask_for(pat, n_expected=N)   # True = excluded-tissue node
    K = int(mask.sum())
    Nk = N - K
    if K == 0 or Nk < ws.MIN_NODES:
        if verbose:
            print(f"[audit_85:{stratify}] SKIP {pat}/{band}: K={K}, Nk={Nk} "
                  "(exclusion is identity or degenerate)")
        return [], {}

    # Observed exclude_<tissue> statistic (drop the K tissue nodes).
    obs_c, obs_r = _rho_splits_for_keep(Ws, ~mask)

    # Random-decimation null: drop K *random* nodes, R_dec times.
    rng = _dec_rng(pat, band, stratify)
    dec_c = np.full(r_dec, np.nan)
    dec_r = np.full(r_dec, np.nan)
    for i in range(r_dec):
        drop = rng.permutation(N)[:K]
        keep = np.ones(N, dtype=bool)
        keep[drop] = False
        rc, rr = _rho_splits_for_keep(Ws, keep)
        dec_c[i] = rc
        dec_r[i] = rr
    gc.collect()

    rows = []
    store = {}
    for substrate, obs, dec in (("cophenetic", obs_c, dec_c),
                                ("raw", obs_r, dec_r)):
        d = dec[np.isfinite(dec)]
        if d.size == 0:
            continue
        p_dec = float(np.mean(d >= obs))
        qs = np.quantile(d, [0.05, 0.25, 0.50, 0.75, 0.95])
        rows.append({
            "patient": pat, "band": band, "substrate": substrate,
            "stratify": stratify, "exclude_config": exclude_cfg,
            "N_full": int(N), "K_dropped": int(K), "n_nodes_kept": int(Nk),
            "frac_dropped": float(K / N),
            "obs_exclude": float(obs),
            "dec_mean": float(np.mean(d)), "dec_std": float(np.std(d, ddof=1)),
            "dec_p5": float(qs[0]), "dec_p25": float(qs[1]),
            "dec_p50": float(qs[2]), "dec_p75": float(qs[3]),
            "dec_p95": float(qs[4]),
            "obs_minus_dec_median": float(obs - qs[2]),
            "p_dec": p_dec,
            "obs_above_dec_median": bool(obs > qs[2]),
            "obs_above_dec_p95": bool(obs > qs[4]),
            "n_decimations": int(d.size),
        })
        store[substrate] = {"obs": float(obs), "dec_full": dec.copy()}
    return rows, store


# ---------------------------------------------------------------------------
# Cohort summary (per band, substrate)
# ---------------------------------------------------------------------------
def _band_order(df: pd.DataFrame) -> list[str]:
    present = set(df.band.unique())
    return [b for b in ws.ALL_BANDS if b in present]


def cohort_summary(per_pat: pd.DataFrame, store: dict, stratify: str
                   ) -> pd.DataFrame:
    """Per (band, substrate) cohort verdict.

    PRIMARY statistic = cohort-level ``cohort_p_dec`` (cohort-median
    ``exclude_<tissue>`` ρ_split vs the null of cohort-median ρ_split under
    independent random K-node decimation per patient — one draw r → one median
    over patients → R_dec cohort samples). This matches the C5/C6 flag (about
    the cohort verdict). ``median_p_dec_per_patient`` is the secondary read."""
    out = []
    for substrate in SUBSTRATES:
        sp = per_pat[per_pat.substrate == substrate]
        for band in _band_order(per_pat):
            d = sp[sp.band == band]
            if d.empty:
                continue
            n = len(d)
            med_p_dec = float(np.median(d.p_dec))
            n_above_med = int(d.obs_above_dec_median.sum())
            n_above_p95 = int(d.obs_above_dec_p95.sum())

            cell = store.get((band, substrate), [])
            obs_vec = np.array([c["obs"] for c in cell], dtype=float)
            cohort_obs = float(np.median(obs_vec)) if obs_vec.size else np.nan
            cohort_p_dec = cohort_dec_med = float("nan")
            cohort_dec_p5 = cohort_dec_p95 = float("nan")
            if cell:
                dec_mat = np.vstack([c["dec_full"] for c in cell])  # (n_pat, R)
                cohort_med_dec = np.nanmedian(dec_mat, axis=0)       # (R,)
                cmd = cohort_med_dec[np.isfinite(cohort_med_dec)]
                if cmd.size and np.isfinite(cohort_obs):
                    cohort_p_dec = float(np.mean(cmd >= cohort_obs))
                    cohort_dec_med = float(np.median(cmd))
                    cohort_dec_p5 = float(np.quantile(cmd, 0.05))
                    cohort_dec_p95 = float(np.quantile(cmd, 0.95))

            tag = "WM" if stratify == "wm" else "epi"
            if np.isfinite(cohort_p_dec) and cohort_p_dec < 0.05:
                verdict = f"{tag}_specific"
            elif np.isfinite(cohort_p_dec) and cohort_p_dec > 0.95:
                verdict = f"{tag}_carried"
            else:
                verdict = "generic_nodecount"
            out.append({
                "stratify": stratify, "substrate": substrate, "band": band,
                "n_patients": n,
                "median_frac_dropped": float(np.median(d.frac_dropped)),
                "cohort_obs_exclude": cohort_obs,
                "cohort_dec_median": cohort_dec_med,
                "cohort_dec_p5": cohort_dec_p5, "cohort_dec_p95": cohort_dec_p95,
                "cohort_obs_minus_dec_median": (cohort_obs - cohort_dec_med
                    if np.isfinite(cohort_dec_med) else float("nan")),
                "cohort_p_dec": cohort_p_dec,
                "median_p_dec_per_patient": med_p_dec,
                "n_above_dec_median": f"{n_above_med}/{n}",
                "n_above_dec_p95": f"{n_above_p95}/{n}",
                "decimation_verdict": verdict,
            })
    return pd.DataFrame(out)


# ===========================================================================
# MODE = pairclass — random-same-size-PAIR-subset null (NO node removal)
# ===========================================================================
# Pair-class configs (epi: nonepi_nonepi / cross / epi_epi; WM: gray_gray /
# cross / wm_wm) restrict the PAIR SET on the FULL-graph cophenetic / raw
# distances — so they have NO node-count confound (the node-decimation control
# above does not apply). Their open hole is a DIFFERENT one: a class with M
# pairs (epi_epi is only ~45) might just be a noisy M-pair estimate of the
# whole-graph trace, not a genuine concentration. This control isolates that:
# draw M random pairs from the FULL pair set R_pair times, recompute ρ_split,
# and ask where the observed class ρ_split sits.
#   p_pair < 0.05  -> CONCENTRATED  (class carries the trace MORE than a random
#                     M-pair subset of the whole graph)
#   p_pair > 0.95  -> DEPLETED      (class carries it LESS; e.g. a genuinely
#                     "spared" tissue-pair class)
#   else           -> generic_paircount (indistinguishable from a random M-pair
#                     subset -> the class value is underpowered, NOT a
#                     localization). This is ORTHOGONAL to the matched-strength
#                     verdict: 'separated' means "not a strength artifact";
#                     'concentrated' means "more trace-bearing than random pairs".
def per_patient_band_pairclass(pat: str, band: str, r_pair: int, stratify: str,
                               verbose: bool) -> tuple[list[dict], dict]:
    """Return (rows, store) for the pair-decimation null on every pair-class
    config. store[(config, substrate)] = {"obs", "null_full"(R,) NaN-padded}."""
    pair_for = STRAT[stratify]["pair_for"]
    configs = STRAT[stratify]["pairclass_configs"]
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in ws.PHASES_4}
    except Exception as e:
        if verbose:
            print(f"[audit_85:pairclass:{stratify}] SKIP {pat}/{band}: "
                  f"load failed: {e}")
        return [], {}
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_85:pairclass:{stratify}] SKIP {pat}/{band}: "
                  "phase shape mismatch")
        return [], {}

    mask = STRAT[stratify]["mask_for"](pat, n_expected=N)  # True = tissue node

    # Full-graph condensed distances (cophenetic + raw), computed ONCE; the
    # pair-class restriction and the random subsets both index into these — so
    # there is never a submatrix rebuild (the whole point: no node-count effect).
    Dc = {ph: lrg_ultrametric_condensed(Ws[ph]) for ph in ws.PHASES_4}
    Dr = {ph: raw_condensed(Ws[ph]) for ph in ws.PHASES_4}
    dT = {
        "cophenetic": Dc["task_test"] - Dc["rest_pre_A"],
        "raw": Dr["task_test"] - Dr["rest_pre_A"],
    }
    dR = {
        "cophenetic": Dc["rest_post"] - Dc["rest_pre_B"],
        "raw": Dr["rest_post"] - Dr["rest_pre_B"],
    }
    n_pairs_full = int(dT["cophenetic"].size)

    rows: list[dict] = []
    store: dict = {}
    for config in configs:
        pk = pair_for(mask, config)
        if pk is None:                      # < MIN_PAIRS for this patient
            if verbose:
                print(f"[audit_85:pairclass:{stratify}] {pat}/{band} {config}: "
                      "undefined (too few class pairs)")
            continue
        M = int(pk.sum())
        # One set of random M-pair index draws, shared across substrates so the
        # cophenetic vs raw comparison sits on identical pair subsets.
        rng = _pair_rng(pat, band, stratify, config)
        draws = [rng.choice(n_pairs_full, size=M, replace=False)
                 for _ in range(r_pair)]
        for substrate in SUBSTRATES:
            a, b = dT[substrate], dR[substrate]
            obs = _rho(a[pk], b[pk])
            null = np.array([_rho(a[idx], b[idx]) for idx in draws], dtype=float)
            d = null[np.isfinite(null)]
            if d.size == 0 or not np.isfinite(obs):
                continue
            p_pair = float(np.mean(d >= obs))
            qs = np.quantile(d, [0.05, 0.25, 0.50, 0.75, 0.95])
            rows.append({
                "patient": pat, "band": band, "config": config,
                "substrate": substrate, "stratify": stratify,
                "N_full": int(N), "n_pairs_full": n_pairs_full,
                "M_class_pairs": M, "frac_pairs": float(M / n_pairs_full),
                "obs_class": float(obs),
                "null_mean": float(np.mean(d)),
                "null_std": float(np.std(d, ddof=1)) if d.size > 1 else float("nan"),
                "null_p5": float(qs[0]), "null_p25": float(qs[1]),
                "null_p50": float(qs[2]), "null_p75": float(qs[3]),
                "null_p95": float(qs[4]),
                "obs_minus_null_median": float(obs - qs[2]),
                "p_pair": p_pair,
                "obs_above_null_median": bool(obs > qs[2]),
                "obs_above_null_p95": bool(obs > qs[4]),
                "obs_below_null_p5": bool(obs < qs[0]),
                "n_draws": int(d.size),
            })
            store[(config, substrate)] = {"obs": float(obs),
                                          "null_full": null.copy()}
    gc.collect()
    return rows, store


def cohort_summary_pairclass(per_pat: pd.DataFrame, store: dict, stratify: str
                             ) -> pd.DataFrame:
    """Per (config, substrate, band) cohort verdict for the pair-decimation null.

    Cohort statistic mirrors the node control: one draw r -> median over patients
    of the random-M-pair ρ_split -> R cohort samples; cohort_p_pair =
    P(cohort-median random >= cohort-median observed)."""
    configs = STRAT[stratify]["pairclass_configs"]
    out = []
    for substrate in SUBSTRATES:
        sp = per_pat[per_pat.substrate == substrate]
        for band in _band_order(per_pat):
            for config in configs:
                d = sp[(sp.band == band) & (sp.config == config)]
                if d.empty:
                    continue
                n = len(d)
                med_p_pair = float(np.median(d.p_pair))
                n_conc = int(d.obs_above_null_p95.sum())
                n_depl = int(d.obs_below_null_p5.sum())

                cell = store.get((band, config, substrate), [])
                obs_vec = np.array([c["obs"] for c in cell], dtype=float)
                cohort_obs = float(np.median(obs_vec)) if obs_vec.size else np.nan
                cohort_p_pair = cohort_null_med = float("nan")
                cohort_null_p5 = cohort_null_p95 = float("nan")
                if cell:
                    null_mat = np.vstack([c["null_full"] for c in cell])  # (n,R)
                    cohort_med_null = np.nanmedian(null_mat, axis=0)       # (R,)
                    cmn = cohort_med_null[np.isfinite(cohort_med_null)]
                    if cmn.size and np.isfinite(cohort_obs):
                        cohort_p_pair = float(np.mean(cmn >= cohort_obs))
                        cohort_null_med = float(np.median(cmn))
                        cohort_null_p5 = float(np.quantile(cmn, 0.05))
                        cohort_null_p95 = float(np.quantile(cmn, 0.95))

                if np.isfinite(cohort_p_pair) and cohort_p_pair < 0.05:
                    verdict = "concentrated"
                elif np.isfinite(cohort_p_pair) and cohort_p_pair > 0.95:
                    verdict = "depleted"
                else:
                    verdict = "generic_paircount"
                out.append({
                    "stratify": stratify, "substrate": substrate, "band": band,
                    "config": config, "n_patients": n,
                    "median_M_class_pairs": float(np.median(d.M_class_pairs)),
                    "median_frac_pairs": float(np.median(d.frac_pairs)),
                    "cohort_obs_class": cohort_obs,
                    "cohort_null_median": cohort_null_med,
                    "cohort_null_p5": cohort_null_p5,
                    "cohort_null_p95": cohort_null_p95,
                    "cohort_obs_minus_null_median": (cohort_obs - cohort_null_med
                        if np.isfinite(cohort_null_med) else float("nan")),
                    "cohort_p_pair": cohort_p_pair,
                    "median_p_pair_per_patient": med_p_pair,
                    "n_concentrated": f"{n_conc}/{n}",
                    "n_depleted": f"{n_depl}/{n}",
                    "pairclass_verdict": verdict,
                })
    return pd.DataFrame(out)


def write_readme_pairclass(per_pat: pd.DataFrame, cohort: pd.DataFrame,
                           stratify: str, r_pair: int, runtime_s: float) -> None:
    cfg = STRAT[stratify]
    out = cfg["out"]
    configs = list(cfg["pairclass_configs"])
    L: list[str] = []
    a = L.append
    a("---")
    a(f"name: {stratify}_pairclass_decimation_control")
    a("scope: random_same_size_pair_subset_null_for_pairclass_localization")
    a(f"stratify: {stratify}")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: first_pass")
    a("build_script: scripts/01_compute/audit/audit_85_wm_decimation_control.py --mode pairclass")
    a(f"companions: {cfg['companions']}")
    a("---")
    a("")
    a(f"# Pair-class localization control — random {r_pair}-draw same-size "
      "pair-subset null (n≤10)")
    a("")
    a("**Head.** Pair-class configs (" + ", ".join(f"`{c}`" for c in configs)
      + ") restrict the *pair set* on the FULL-graph cophenetic / raw distances, "
      "so unlike the subgraph exclusions they have **no node-count confound**. "
      "Their open hole is power: a class with M pairs may just be a noisy M-pair "
      "estimate of the whole-graph trace. For each (patient, band, class) we draw "
      f"M random pairs from the full pair set R={r_pair} times, recompute the "
      "split-baseline ρ_split, and locate the observed class ρ_split in that "
      "random-subset distribution.")
    a("")
    a("- **concentrated** (cohort_p_pair < 0.05): the class carries the trace "
      "MORE than a random M-pair subset of the whole graph ⇒ genuine "
      "localization, not a pair-count artifact.")
    a("- **depleted** (cohort_p_pair > 0.95): the class carries it LESS ⇒ a "
      "genuinely spared / anti pair-class.")
    a("- **generic_paircount** (0.05 ≤ cohort_p_pair ≤ 0.95): indistinguishable "
      "from a random M-pair subset ⇒ the class value is **underpowered**, NOT a "
      "localization. Do not read it as sparing or concentration.")
    a("")
    a("**Orthogonal to matched-strength.** A `separated` matched-strength "
      "verdict means 'not a strength artifact'; `concentrated` here means 'more "
      "trace-bearing than random pairs of the same count'. A class can be "
      "`separated` yet `generic_paircount` (survives the strength null but is "
      "not more localized than the graph average).")
    a("")
    a(f"Tissue: {cfg['tissue']}.")
    a("")
    a("## Cohort verdict (substrate × band × class)")
    a("")
    a("| substrate | band | class | n | M pairs | cohort obs | null med | "
      "obs−null | cohort p_pair | verdict |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for substrate in SUBSTRATES:
        cs = cohort[cohort.substrate == substrate]
        for band in _band_order(per_pat):
            for config in configs:
                r = cs[(cs.band == band) & (cs.config == config)]
                if r.empty:
                    continue
                r = r.iloc[0]
                a(f"| {substrate} | {BRAIN_BAND_TEX_DICT.get(band, band)} "
                  f"| {config} | {r['n_patients']} "
                  f"| {r['median_M_class_pairs']:.0f} "
                  f"| {r['cohort_obs_class']:+.3f} "
                  f"| {r['cohort_null_median']:+.3f} "
                  f"| {r['cohort_obs_minus_null_median']:+.3f} "
                  f"| {r['cohort_p_pair']:.3f} "
                  f"| **{r['pairclass_verdict']}** |")
    a("")
    a("## Caveats")
    a("- The null draws from the FULL pair set, so its mean ≈ the `full`-config "
      "ρ_split; this control therefore asks 'is the class more trace-bearing "
      "than the graph average at matched pair count', the relevant localization "
      "question. It does NOT replace matched-strength (orthogonal axis).")
    a("- Small-M classes (epi_epi ~45 pairs) have a wide null by construction; a "
      "`generic_paircount` verdict there is the honest 'cannot tell' answer.")
    if stratify == "epi":
        a("- Pat_15 has 0 epi nodes → `epi_epi` / `cross` undefined → that class "
          "cohort is n=9.")
    a("")
    a("## Provenance")
    a(f"- R = {r_pair} random pair-subset draws per (patient, band, class); seed "
      f"{SEED_PAIR} (per-cell SeedSequence, salted by stratify + config).")
    a("- Cohort: " + ", ".join(ws.COHORT))
    a("- ρ_split = Spearman(D_task − D_preA, D_post − D_preB) restricted to the "
      "pair subset; positive = TRACE. FC = imcoh_abs; τ = 1/λ_max.")
    a("- Substrates: LRG cophenetic ρ_split, raw |ImCoh| ρ_split (full-graph "
      "distances, pairs restricted — no submatrix rebuild).")
    a(f"- Wall-clock: {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `pairclass_decimation_per_patient.csv` — one row per pat×band×class×substrate")
    a("- `pairclass_decimation_cohort.csv` — per (substrate, band, class) verdict")
    (out / "README_pairclass.md").write_text("\n".join(L) + "\n",
                                             encoding="utf-8")


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(per_pat: pd.DataFrame, cohort: pd.DataFrame, stratify: str,
                 r_dec: int, runtime_s: float) -> None:
    cfg = STRAT[stratify]
    out = cfg["out"]
    L: list[str] = []
    a = L.append
    a("---")
    a(f"name: {stratify}_decimation_control")
    a("scope: size_matched_random_node_decimation_control_for_subset_exclusion")
    a(f"stratify: {stratify} ({cfg['exclude_cfg']})")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: first_pass")
    a("build_script: scripts/01_compute/audit/audit_85_wm_decimation_control.py")
    a(f"companions: {cfg['companions']}")
    a("---")
    a("")
    a(f"# {cfg['exclude_cfg']} node-count control — random-node decimation (n≤10)")
    a("")
    a("**Head.** Isolates the node-count confound on the "
      f"`{cfg['exclude_cfg']}` exclusion control. For each (patient, band) we "
      f"drop K = #{stratify} nodes *at random* R_dec={r_dec} times, recompute "
      "the split-baseline ρ_split (cophenetic + raw) on each N−K-node subgraph, "
      f"and ask where the observed `{cfg['exclude_cfg']}` ρ_split sits in that "
      "random-decimation distribution. obs in the UPPER tail (cohort_p_dec "
      f"small) ⇒ removing {stratify} helps MORE than removing random nodes ⇒ "
      "tissue-specific. obs ≈ random median (cohort_p_dec ≈ 0.5) ⇒ a generic "
      "node-count effect, not tissue. obs in the LOWER tail (≈ 1) ⇒ the tissue "
      "was carrying the trace.")
    a("")
    a(f"Tissue: {cfg['tissue']}.")
    a("")
    a("## Cohort verdict (band × substrate)")
    a("")
    a("Primary statistic = cohort-level `cohort_p_dec`.")
    a("")
    a("| substrate | band | n | frac drop | cohort obs | cohort dec med | "
      "obs−dec | cohort p_dec | med p_dec (per-pat) | verdict |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for substrate in SUBSTRATES:
        cs = cohort[cohort.substrate == substrate]
        for band in _band_order(per_pat):
            r = cs[cs.band == band]
            if r.empty:
                continue
            r = r.iloc[0]
            a(f"| {substrate} | {BRAIN_BAND_TEX_DICT.get(band, band)} "
              f"| {r['n_patients']} | {r['median_frac_dropped']:.2f} "
              f"| {r['cohort_obs_exclude']:+.3f} | {r['cohort_dec_median']:+.3f} "
              f"| {r['cohort_obs_minus_dec_median']:+.3f} "
              f"| {r['cohort_p_dec']:.3f} | {r['median_p_dec_per_patient']:.3f} "
              f"| **{r['decimation_verdict']}** |")
    a("")
    tag = "WM" if stratify == "wm" else "epi"
    a("## Verdict labels (gate on cohort_p_dec)")
    a("")
    a(f"- **{tag}_specific**: cohort_p_dec < 0.05 — the cohort-median "
      f"`{cfg['exclude_cfg']}` ρ_split is in the upper tail of the random "
      f"K-node-loss distribution. Removing {stratify} raises the cohort trace "
      "more than removing the same count of random nodes ⇒ tissue-specific.")
    a("- **generic_nodecount**: 0.05 ≤ cohort_p_dec ≤ 0.95 — typical of removing "
      "any K nodes ⇒ the change is the graph getting smaller, not tissue "
      "membership. Demote any `emerge`/`strengthen` to 'not tissue-specific'.")
    a(f"- **{tag}_carried**: cohort_p_dec > 0.95 — lower tail ⇒ removing "
      f"{stratify} *hurts* relative to random removal ⇒ tissue carried the trace.")
    a("")
    a("## Caveats")
    a("- Random decimation can sever a whole probe / hemisphere; the null "
      "deliberately includes such subgraphs (the generic-node-loss reference), "
      "so it is broad.")
    a("- This control holds NODE COUNT fixed and randomizes tissue identity; it "
      "is the complement of the matched-strength null (which holds per-node "
      "strength fixed). A `tissue_specific` verdict is 'specific to the tissue "
      "node set', which the matched-strength control further constrains.")
    if stratify == "epi":
        a("- Pat_15 has 0 epi nodes → exclude_epi is identity → skipped; the epi "
          "decimation cohort is effectively n=9.")
        a("- This control tests the **subgraph** `exclude_epi` claim only. The "
          "C5 **pair-class** decomposition (CROSS / NONEPI / EPI↔EPI) is a "
          "*same-graph* analysis with NO node-count confound and is unaffected.")
    a("")
    a("## Provenance")
    a(f"- R_dec = {r_dec} random decimations per (patient, band); seed "
      f"{SEED_DEC} (per-cell SeedSequence; epi salt 991, WM unsalted/frozen).")
    a("- Cohort: " + ", ".join(ws.COHORT))
    a("- Phases: " + ", ".join(ws.PHASES_4))
    a("- Substrates: LRG cophenetic ρ_split, raw |ImCoh| ρ_split.")
    a("- ρ_split = Spearman(D_task − D_preA, D_post − D_preB); positive = TRACE.")
    a("- FC method: imcoh_abs; τ = 1/λ_max; average-linkage ultrametric.")
    a(f"- Tissue: {cfg['tissue']}.")
    a(f"- Wall-clock: {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `decimation_control_per_patient.csv` — one row per pat×band×substrate")
    a("- `decimation_control_cohort.csv` — per (band, substrate) verdict")
    (out / "README_decimation.md").write_text("\n".join(L) + "\n",
                                              encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _run_stratify(stratify: str, bands: list[str], patients: list[str],
                  r_dec: int, verbose: bool) -> None:
    out = STRAT[stratify]["out"]
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    all_rows: list[dict] = []
    store: dict = {}
    for band in bands:
        for pat in patients:
            tc = time.time()
            rows, pat_store = per_patient_band(pat, band, r_dec, stratify,
                                               verbose)
            all_rows.extend(rows)
            for substrate, s in pat_store.items():
                store.setdefault((band, substrate), []).append(s)
            bits = " ".join(
                f"{r['substrate'][:4]}:obs{r['obs_exclude']:+.2f}"
                f"(dec{r['dec_p50']:+.2f},p{r['p_dec']:.2f})"
                for r in rows)
            if rows:
                print(f"[audit_85:{stratify}] {pat}/{band}: {bits} "
                      f"({time.time()-tc:.1f}s)")
    runtime = time.time() - t0
    per_pat = pd.DataFrame(all_rows)
    per_pat.to_csv(out / "decimation_control_per_patient.csv", index=False)
    cohort = cohort_summary(per_pat, store, stratify)
    cohort.to_csv(out / "decimation_control_cohort.csv", index=False)
    write_readme(per_pat, cohort, stratify, r_dec, runtime)
    print(f"[audit_85:{stratify}] {len(per_pat)} rows, {runtime:.1f}s -> {out}")
    if not cohort.empty:
        print(cohort[["substrate", "band", "cohort_obs_exclude",
                      "cohort_dec_median", "cohort_p_dec",
                      "decimation_verdict"]].to_string(index=False))


def _run_stratify_pairclass(stratify: str, bands: list[str],
                            patients: list[str], r_pair: int, verbose: bool
                            ) -> None:
    out = STRAT[stratify]["out"]
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    all_rows: list[dict] = []
    store: dict = {}
    for band in bands:
        for pat in patients:
            tc = time.time()
            rows, pat_store = per_patient_band_pairclass(
                pat, band, r_pair, stratify, verbose)
            all_rows.extend(rows)
            for (config, substrate), s in pat_store.items():
                store.setdefault((band, config, substrate), []).append(s)
            bits = " ".join(
                f"{r['config']}/{r['substrate'][:4]}:obs{r['obs_class']:+.2f}"
                f"(null{r['null_p50']:+.2f},p{r['p_pair']:.2f})"
                for r in rows if r["substrate"] == "cophenetic")
            if rows:
                print(f"[audit_85:pairclass:{stratify}] {pat}/{band}: {bits} "
                      f"({time.time()-tc:.1f}s)")
    runtime = time.time() - t0
    per_pat = pd.DataFrame(all_rows)
    per_pat.to_csv(out / "pairclass_decimation_per_patient.csv", index=False)
    cohort = cohort_summary_pairclass(per_pat, store, stratify)
    cohort.to_csv(out / "pairclass_decimation_cohort.csv", index=False)
    write_readme_pairclass(per_pat, cohort, stratify, r_pair, runtime)
    print(f"[audit_85:pairclass:{stratify}] {len(per_pat)} rows, "
          f"{runtime:.1f}s -> {out}")
    if not cohort.empty:
        print(cohort[["substrate", "band", "config", "cohort_obs_class",
                      "cohort_null_median", "cohort_p_pair",
                      "pairclass_verdict"]].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="node", choices=["node", "pairclass"],
                    help="node = size-matched random-NODE decimation (subgraph "
                         "exclusion confound); pairclass = random same-size "
                         "PAIR-subset null (pair-class localization).")
    ap.add_argument("--stratify", nargs="+", default=["wm", "epi"],
                    choices=["wm", "epi"])
    ap.add_argument("--bands", nargs="+", default=list(ws.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(ws.COHORT))
    ap.add_argument("--r-dec", type=int, default=R_DEC_DEFAULT)
    ap.add_argument("--r-pair", type=int, default=R_PAIR_DEFAULT)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("[audit_85] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    for stratify in args.stratify:
        if args.mode == "node":
            print(f"\n========== mode=node  stratify = {stratify} "
                  f"({STRAT[stratify]['exclude_cfg']}) ==========")
            _run_stratify(stratify, args.bands, args.patients, args.r_dec,
                          args.verbose)
        else:
            print(f"\n========== mode=pairclass  stratify = {stratify} "
                  f"({', '.join(STRAT[stratify]['pairclass_configs'])}) "
                  "==========")
            _run_stratify_pairclass(stratify, args.bands, args.patients,
                                    args.r_pair, args.verbose)


if __name__ == "__main__":
    main()
