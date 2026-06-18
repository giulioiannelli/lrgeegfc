#!/usr/bin/env python3
"""Audit 77 — epileptic-node-stratified cophenetic ρ_split trace.

Unifies the epi-stratified cophenetic LRG trace across six configurations,
all six bands, the full n=10 cohort, into ONE tidy cache feeding the preprint
figures. Subsumes/extends audit_68_{alpha,beta}_epi_exclusion (exclude-epi,
α/β only) to all bands, and adds the epi-only subgraph and the epi↔non-epi
interface-pair decomposition.

Critical preamble (per CLAUDE.md rule, before any code)
=======================================================
(1) **Claim.** The established cophenetic ρ_split trace (positive = TRACE;
    ρ_split = Spearman(D_task − D_preA, D_post − D_preB), split-baseline) is
    (a) robust to removing epileptic-zone contacts (exclude_epi), (b) present
    or absent within the epileptic subgraph (epi_only), and (c) carried
    disproportionately by a particular epi-stratified pair class — non-epi↔
    non-epi, epi↔non-epi interface (cross), or epi↔epi.

(2) **Null.** Matched-strength 4-cycle ±δ surrogate, R=200. For subgraph
    configs the ensemble is regenerated on the submatrix W[keep, keep]
    (own cache dir, never clobbering the canonical full-graph cache). For
    pair-class configs the null is the FULL-graph matched-strength ensemble
    with ρ_split recomputed on the restricted pair set (the audit_71
    pair_keep mechanism) — i.e. it holds the global per-node strength
    sequence fixed and asks whether the class-restricted co-movement survives.

(3) **Strongest plausible alternatives.** (i) Class size: non-epi↔non-epi
    always has the most pairs, so a more stable Spearman, independent of any
    trace localization. (ii) The global strength sequence rather than
    class-specific topology. (iii) For epi_only, small-N Spearman noise
    (6–30 nodes; 15–435 pairs) masquerading as signal, with a near-degenerate
    matched-strength null (few feasible 4-cycles on a tiny graph).

(4) **Null's mechanical reach.** The matched-strength null fixes per-node
    strength exactly, so a surviving statistic is NOT explained by (ii). It
    does NOT control for (i) — class-size differences are reported as n_pairs
    per class and must be read alongside the verdict — nor does the pair-class
    null randomize *within* the class (it is a global-rewiring null). For
    epi_only it controls strength but the surrogate is near-degenerate, so
    epi_only "separated" verdicts are NOT trustworthy and are flagged
    exploratory.

(5) **Falsification + limitations.** "Trace is interface-carried" requires
    cross = separated AND nonepi_nonepi weakening vs full. "Trace survives in
    healthy tissue" requires exclude_epi = persist. Limitations: pair-class
    null is global-rewiring not within-class; epi_only underpowered; class
    sizes unequal; Pat_13 (30 epi, 25%) leverages every epi config → a
    leave-Pat_13-out cohort p is reported.

Outputs
-------
``data/audit/epi_stratified/``
    cophenetic_per_patient.csv   tidy long-format (one row per pat×band×config)
    cophenetic_cohort.csv        per (band, config) cohort verdict + sensitivity
    README.md

Surrogate caches (per config, no clobber)
-----------------------------------------
- full / pair-class → ``data/cache/matched_strength_surrogate_lrg/`` (seed 20260511)
- exclude_epi       → ``data/cache/matched_strength_surrogate_epi_excluded_lrg/`` (seed 20260514, audit_67)
- epi_only          → ``data/cache/matched_strength_surrogate_epi_only_lrg/`` (seed 20260511)
Pat_15 has 0 epi contacts → epi_only / cross / epi_epi are undefined (n=9 by
construction; emitted with defined=False).
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
from lrg_eegfc.utils.surrogate import (
    cophenetic_condensed_from_eigs,
    load_or_compute_eigs_at_path,
)

# Reuse audit_63 cophenetic + half-FC pipeline (no copy).
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
    lrg_ultrametric_condensed,
)
import _epi_stratify as es  # type: ignore


SUBSTRATE = "cophenetic"
OUT = ROOT / "data" / "audit" / "epi_stratified"
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Per-surrogate cophenetic reconstruction from a cached eigendecomposition
# ---------------------------------------------------------------------------
def _surr_cophenetic_stack(W: np.ndarray, cache_path: Path,
                           rng: np.random.Generator, n_surr: int,
                           swap_factor: int, verbose: bool) -> np.ndarray:
    """Return (R, n_pairs) condensed cophenetic vectors for the matched-strength
    ensemble of `W`, loading/generating its eigendecomposition at `cache_path`.
    Rows for strength-violating surrogates are NaN.
    """
    evals, evecs = load_or_compute_eigs_at_path(
        cache_path, W, n_surr, swap_factor, rng, verbose=verbose)
    R = evals.shape[0]
    N = W.shape[0]
    n_pairs = N * (N - 1) // 2
    S = np.full((R, n_pairs), np.nan, dtype=np.float64)
    for r in range(R):
        if not np.isfinite(evals[r]).all():
            continue
        S[r] = cophenetic_condensed_from_eigs(evals[r], evecs[r])
    return S


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    rho, _ = spearmanr(a, b)
    return float(rho)


def _stats_row(pat: str, band: str, config: str, dT_obs: np.ndarray,
               dR_obs: np.ndarray, dT_surr: np.ndarray, dR_surr: np.ndarray,
               pair_keep: np.ndarray | None, n_nodes: int,
               surrogate_source: str, n_swaps: int, seed: int) -> dict:
    """ρ_split observed + surrogate stats for one config cell."""
    if pair_keep is None:
        sel = slice(None)
        n_pairs = int(dT_obs.size)
    else:
        sel = pair_keep
        n_pairs = int(pair_keep.sum())
    obs = _rho(dT_obs[sel], dR_obs[sel])
    R = dT_surr.shape[0]
    srho = np.full(R, np.nan)
    for r in range(R):
        a, b = dT_surr[r][sel], dR_surr[r][sel]
        if np.isfinite(a).all() and np.isfinite(b).all():
            srho[r] = _rho(a, b)
    surr = srho[np.isfinite(srho)]
    return _row(pat, band, config, obs, surr, n_nodes, n_pairs,
                surrogate_source, n_swaps, seed, defined=True)


def _row(pat: str, band: str, config: str, obs: float, surr: np.ndarray,
         n_nodes: int, n_pairs: int, surrogate_source: str, n_swaps: int,
         seed: int, defined: bool) -> dict:
    if surr is not None and surr.size > 0:
        smean = float(np.mean(surr)); sstd = float(np.std(surr, ddof=1))
        qs = np.quantile(surr, [0.05, 0.25, 0.50, 0.75, 0.95])
        z = (obs - smean) / sstd if sstd > 0 else float("nan")
        p_up = float(np.mean(surr >= obs))
        n_s = int(surr.size)
    else:
        smean = sstd = float("nan")
        qs = [float("nan")] * 5
        z = float("nan"); p_up = float("nan"); n_s = 0
    return {
        "patient": pat, "band": band, "config": config,
        "config_class": es.CONFIG_CLASS[config], "substrate": SUBSTRATE,
        "k": np.nan, "n_nodes": int(n_nodes), "n_pairs": int(n_pairs),
        "obs_stat": float(obs) if obs == obs else float("nan"),
        "surr_mean": smean, "surr_std": sstd,
        "surr_p5": float(qs[0]), "surr_p25": float(qs[1]),
        "surr_p50": float(qs[2]), "surr_p75": float(qs[3]),
        "surr_p95": float(qs[4]),
        "obs_z": z, "obs_p_one_sided": p_up,
        "n_surrogates": n_s, "n_swaps_per_surrogate": int(n_swaps),
        "surrogate_source": surrogate_source, "seed": int(seed),
        "defined": bool(defined),
    }


def _undefined_row(pat: str, band: str, config: str, n_nodes: int) -> dict:
    return _row(pat, band, config, float("nan"), None, n_nodes, 0,
                "undefined", 0, es.config_seed(config), defined=False)


# ---------------------------------------------------------------------------
# Per (patient, band): all requested configs, sharing loaded FCs / full eigs
# ---------------------------------------------------------------------------
def per_patient_band(pat: str, band: str, configs: list[str], n_surr: int,
                     swap_factor: int, verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in es.PHASES_4}
    except Exception as e:
        if verbose:
            print(f"[audit_77] SKIP {pat}/{band}: load failed: {e}")
        return []
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_77] SKIP {pat}/{band}: phase shape mismatch")
        return []
    epi = es.epi_mask_for(pat, n_expected=N)
    rows: list[dict] = []

    # ---- Full-graph cophenetic: serves 'full' + pair-class configs ----
    pc = [c for c in configs if c in es.PAIRCLASS_CONFIGS]
    if "full" in configs or pc:
        D = {ph: lrg_ultrametric_condensed(Ws[ph]) for ph in es.PHASES_4}
        dT_o = D["task_test"] - D["rest_pre_A"]
        dR_o = D["rest_post"] - D["rest_pre_B"]
        surr = {}
        for ph in es.PHASES_4:
            surr[ph] = _surr_cophenetic_stack(
                Ws[ph], es.surr_eig_path("full", pat, band, ph, n_surr, swap_factor),
                es.cell_rng(pat, band, ph, "full"), n_surr, swap_factor, verbose)
        dT_s = surr["task_test"] - surr["rest_pre_A"]
        dR_s = surr["rest_post"] - surr["rest_pre_B"]
        n_swaps_full = swap_factor * (N * (N - 1)) // 2
        if "full" in configs:
            rows.append(_stats_row(pat, band, "full", dT_o, dR_o, dT_s, dR_s,
                                   None, N, "fullgraph", n_swaps_full,
                                   es.config_seed("full")))
        for c in pc:
            pk = es.pair_keep_for_config(epi, c)
            if pk is None:
                rows.append(_undefined_row(pat, band, c, N))
            else:
                rows.append(_stats_row(pat, band, c, dT_o, dR_o, dT_s, dR_s,
                                       pk, N, "fullgraph_restricted",
                                       n_swaps_full, es.config_seed(c)))
        del surr, dT_s, dR_s
        gc.collect()

    # ---- Subgraph configs: exclude_epi, epi_only ----
    for c in [c for c in configs if c in ("exclude_epi", "epi_only")]:
        keep = es.node_mask_for_config(epi, c)
        if keep is None:
            rows.append(_undefined_row(pat, band, c, 0))
            continue
        idx = np.ix_(keep, keep)
        Wk = {ph: np.ascontiguousarray(Ws[ph][idx]) for ph in es.PHASES_4}
        Nk = int(keep.sum())
        Dk = {ph: lrg_ultrametric_condensed(Wk[ph]) for ph in es.PHASES_4}
        dT_o = Dk["task_test"] - Dk["rest_pre_A"]
        dR_o = Dk["rest_post"] - Dk["rest_pre_B"]
        surr = {}
        for ph in es.PHASES_4:
            surr[ph] = _surr_cophenetic_stack(
                Wk[ph], es.surr_eig_path(c, pat, band, ph, n_surr, swap_factor),
                es.cell_rng(pat, band, ph, c), n_surr, swap_factor, verbose)
        dT_s = surr["task_test"] - surr["rest_pre_A"]
        dR_s = surr["rest_post"] - surr["rest_pre_B"]
        n_swaps = swap_factor * (Nk * (Nk - 1)) // 2
        rows.append(_stats_row(pat, band, c, dT_o, dR_o, dT_s, dR_s, None, Nk,
                               "submatrix", n_swaps, es.config_seed(c)))
        del surr, dT_s, dR_s
        gc.collect()

    return rows


# ---------------------------------------------------------------------------
# Cohort summary
# ---------------------------------------------------------------------------
def _band_order(df: pd.DataFrame) -> list[str]:
    present = set(df.band.unique())
    return [b for b in es.ALL_BANDS if b in present]


def cohort_summary(per_pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    bands = _band_order(per_pat)
    # full-config cohort p per band (the sensitivity baseline)
    full_p: dict[str, float] = {}
    for band in bands:
        d = per_pat[(per_pat.band == band) & (per_pat.config == "full")
                    & per_pat.defined]
        n_above = int((d.obs_p_one_sided < 0.05).sum())
        full_p[band] = es.cohort_verdict(
            d.obs_stat.values, d.surr_p50.values, n_above)["wilcoxon_p"]
    for band in bands:
        for config in [c for c in es.ALL_CONFIGS
                       if c in set(per_pat[per_pat.band == band].config)]:
            sub = per_pat[(per_pat.band == band) & (per_pat.config == config)]
            d = sub[sub.defined]
            n_total = len(sub)
            n_above = int((d.obs_p_one_sided < 0.05).sum())
            v = es.cohort_verdict(d.obs_stat.values, d.surr_p50.values, n_above)
            d2 = d[d.patient != "Pat_13"]
            n_above2 = int((d2.obs_p_one_sided < 0.05).sum())
            v2 = es.cohort_verdict(d2.obs_stat.values, d2.surr_p50.values,
                                   n_above2)
            flag = ("baseline" if config == "full"
                    else es.sensitivity_flag(
                        full_p.get(band, np.nan), v["wilcoxon_p"],
                        cfg_defined=(v["verdict"] != "undefined")))
            out.append({
                "band": band, "config": config,
                "config_class": es.CONFIG_CLASS[config], "substrate": SUBSTRATE,
                "n_patients": n_total, "n_defined": v["n_defined"],
                "obs_median": v["med_obs"],
                "surr_median_median": v["med_surr"],
                "median_n_pairs": float(np.median(d.n_pairs)) if len(d) else np.nan,
                "n_above_own_surrogate": f"{n_above}/{v['n_defined']}",
                "paired_wilcoxon_z": v["wilcoxon_z"],
                "paired_wilcoxon_p": v["wilcoxon_p"],
                "lopat13_wilcoxon_p": v2["wilcoxon_p"],
                "verdict": v["verdict"], "sensitivity_flag": flag,
            })
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(per_pat: pd.DataFrame, cohort: pd.DataFrame,
                 n_surr: int, swap_factor: int, runtime_s: float) -> None:
    L: list[str] = []
    a = L.append
    a("---")
    a("name: epi_stratified_cophenetic")
    a("scope: epi_node_stratified_cophenetic_rho_split_trace_all_bands")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: first_pass")
    a("build_script: scripts/01_compute/audit/audit_77_epi_stratified_cophenetic.py")
    a("---")
    a("")
    a("# Epi-node-stratified cophenetic ρ_split trace (all bands, n=10)")
    a("")
    a("**Head.** The split-baseline cophenetic ρ_split trace (positive = TRACE) "
      "re-run across six epi views — full graph, non-epi↔non-epi pairs, "
      "epi-excluded subgraph, epi↔non-epi interface pairs, epi↔epi pairs, "
      "epi-only subgraph — for all six bands under R="
      f"{n_surr} matched-strength surrogacy. Subgraph configs regenerate the "
      "null on the submatrix; pair-class configs restrict the full-graph null "
      "(holding the global strength sequence fixed). epi_only is exploratory "
      "(tiny graphs, near-degenerate null).")
    a("")
    a("## Cohort verdict — band × config")
    a("")
    a("| band | config | n_def | obs med | surr med | n>surr | Wilcoxon p | "
      "LO-P13 p | verdict | vs full |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for band in _band_order(per_pat):
        for _, r in cohort[cohort.band == band].iterrows():
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {r['config']} "
              f"| {r['n_defined']} | {r['obs_median']:+.3f} "
              f"| {r['surr_median_median']:+.3f} | {r['n_above_own_surrogate']} "
              f"| {r['paired_wilcoxon_p']:.4f} | {r['lopat13_wilcoxon_p']:.4f} "
              f"| {r['verdict']} | **{r['sensitivity_flag']}** |")
    a("")
    a("## Sensitivity labels (config vs full)")
    a("")
    a("- **persist**: full p<0.05 AND config p<0.05 — trace present in this view.")
    a("- **weaken**: full p<0.05 AND config p≥0.05 — trace needs the removed pairs/nodes.")
    a("- **emerge**: full p≥0.05 AND config p<0.05 — view reveals a trace masked at full graph.")
    a("- **absent**: neither significant. **undefined**: <8 defined patients (e.g. Pat_15 epi configs).")
    a("")
    a("## Caveats")
    a("")
    a("- **epi_only is exploratory**: 6–30 nodes (Pat_15=0), 15–435 pairs; "
      "matched-strength on tiny graphs is near-degenerate. Do not trust "
      "epi_only 'separated' verdicts.")
    a("- **Pair-class null is global-rewiring**, not within-class; it tests "
      "class localization against a strength-preserving rewiring of the whole "
      "graph. Class sizes are unequal — read n_pairs (median per config below).")
    a("- **Pat_13** (30 epi, 25%) leverages every epi config → LO-P13 p column "
      "is the leave-Pat_13-out cohort robustness line.")
    a("")
    a("## Median pairs per config (class-size context)")
    a("")
    a("| band | " + " | ".join(es.ALL_CONFIGS) + " |")
    a("|---|" + "|".join(["---"] * len(es.ALL_CONFIGS)) + "|")
    for band in _band_order(per_pat):
        cells = []
        for c in es.ALL_CONFIGS:
            row = cohort[(cohort.band == band) & (cohort.config == c)]
            cells.append("—" if row.empty or not np.isfinite(row.iloc[0]["median_n_pairs"])
                         else f"{int(row.iloc[0]['median_n_pairs'])}")
        a(f"| {band} | " + " | ".join(cells) + " |")
    a("")
    a("## Provenance")
    a(f"- N_surrogates = {n_surr}, SWAP_FACTOR = {swap_factor}")
    a("- Cohort: " + ", ".join(es.COHORT))
    a("- Phases: " + ", ".join(es.PHASES_4))
    a("- FC method: imcoh_abs; τ = 1/λ_max; average-linkage ultrametric.")
    a("- Seeds: full/pair-class 20260511, exclude_epi 20260514 (audit_67 reuse), "
      "epi_only 20260511.")
    a(f"- Wall-clock: {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `cophenetic_per_patient.csv` — tidy long-format per pat×band×config")
    a("- `cophenetic_cohort.csv` — per (band, config) cohort verdict + sensitivity")
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--configs", nargs="+", default=list(es.ALL_CONFIGS))
    ap.add_argument("--n-surrogates", type=int, default=es.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=es.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("[audit_77] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    all_rows: list[dict] = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            rows = per_patient_band(pat, band, args.configs,
                                    args.n_surrogates, args.swap_factor,
                                    args.verbose)
            all_rows.extend(rows)
            bits = " ".join(
                f"{r['config'][:5]}:{r['obs_stat']:+.2f}"
                f"(p{r['obs_p_one_sided']:.2f})"
                for r in rows if r["defined"])
            print(f"[audit_77] {pat}/{band}: {bits} ({time.time()-tc:.1f}s)")

    runtime = time.time() - t0
    per_pat = pd.DataFrame(all_rows)
    per_pat.to_csv(OUT / "cophenetic_per_patient.csv", index=False)
    cohort = cohort_summary(per_pat)
    cohort.to_csv(OUT / "cophenetic_cohort.csv", index=False)
    write_readme(per_pat, cohort, args.n_surrogates, args.swap_factor, runtime)
    print(f"[audit_77] {len(per_pat)} rows, {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
