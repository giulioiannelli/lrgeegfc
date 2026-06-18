#!/usr/bin/env python3
"""Audit 78 — epileptic-node-stratified Grassmann T_G(k) trace.

Grassmann companion to audit_77. The chordal subspace statistic
T_G(k) = d_chord(rest_pre, task; k) − d_chord(task, rest_post; k)
(positive = TRACE) lives on the whole k-dimensional slow-mode subspace of the
Laplacian, so it admits only the SUBGRAPH views (full / exclude_epi /
epi_only) — there is no pair-restriction analogue, so the pair-class
decomposition is cophenetic-only (audit_77).

To avoid recomputing what is already validated:
- ``full``        is re-emitted from audit_66 (grassmann_matched_strength_surrogate).
- ``exclude_epi`` is re-emitted from audit_67 (grassmann_epi_exclusion).
- ``epi_only``    is computed fresh here (reusing the epi_only eigendecomposition
  cache that audit_77 builds at seed 20260511).

Critical preamble (per CLAUDE.md rule)
======================================
(1) **Claim.** The Grassmann subspace trace (β, γ_l load-bearing in audit_66)
    survives epi-exclusion (audit_67 already shows the all-band surface) and is
    present/absent within the epileptic subgraph alone (epi_only).
(2) **Null.** R=200 matched-strength 4-cycle ±δ on the relevant graph (full;
    W[non_epi]; W[epi]).
(3) **Strongest alternative.** epi_only T_G is dominated by small-N subspace
    instability: with n_epi∈{6..30}, k must satisfy k+1≤n_epi, so the cohort
    has ≥8 patients only at k≲5 (Pat_15 has 0 epi and is structurally absent),
    and the matched-strength null on a 6–30 node graph is near-degenerate.
(4) **Mechanical reach.** The null fixes per-node strength, so a surviving T_G
    is not a strength artifact; it does NOT rescue epi_only from small-N
    subspace noise. epi_only cohort verdicts past the ≥8-patient frontier are
    NOT trustworthy.
(5) **Falsification + limits.** "Trace lives in healthy tissue" → exclude_epi
    persists (audit_67 surface). epi_only is exploratory only; its k-coverage
    collapse is reported explicitly. Pat_13 (30 epi) is the only patient with
    deep epi_only k-coverage → it dominates high-k epi_only cells.

Outputs
-------
``data/audit/epi_stratified/``
    grassmann_per_patient_per_k.csv   tidy long-format (one row per pat×band×config×k)
    grassmann_cohort.csv              per (band, config, k) cohort verdict + sensitivity
    README_grassmann.md
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_67_grassmann_epi_exclusion import (  # type: ignore
    K_GRID,
    chordal,            # noqa: F401  (used transitively via t_g_at_all_k)
    ensure_half_fcs,
    laplacian_eig,
    load_phase_fc_full,
    t_g_at_all_k,
    topk_basis,
)
import _epi_stratify as es  # type: ignore


SUBSTRATE = "grassmann"
OUT = ROOT / "data" / "audit" / "epi_stratified"
OUT.mkdir(parents=True, exist_ok=True)
A66 = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate" \
    / "per_patient_per_band_per_k.csv"
A67 = ROOT / "data" / "audit" / "grassmann_epi_exclusion" \
    / "per_patient_per_band_per_k.csv"

_SCHEMA = ["patient", "band", "config", "config_class", "substrate", "k",
           "n_nodes", "n_pairs", "obs_stat", "surr_mean", "surr_std",
           "surr_p5", "surr_p25", "surr_p50", "surr_p75", "surr_p95",
           "obs_z", "obs_p_one_sided", "n_surrogates",
           "n_swaps_per_surrogate", "surrogate_source", "seed", "defined"]


# ---------------------------------------------------------------------------
# Re-emit prior Grassmann CSVs into the unified schema
# ---------------------------------------------------------------------------
def _map_prior(df: pd.DataFrame, config: str, n_nodes_col: str,
               surrogate_source: str, seed: int) -> pd.DataFrame:
    out = pd.DataFrame({
        "patient": df.patient, "band": df.band, "config": config,
        "config_class": "subgraph", "substrate": SUBSTRATE, "k": df.k,
        "n_nodes": df[n_nodes_col], "n_pairs": np.nan,
        "obs_stat": df.obs_T_G, "surr_mean": df.surr_T_G_mean,
        "surr_std": df.surr_T_G_std, "surr_p5": df.surr_T_G_p5,
        "surr_p25": df.surr_T_G_p25, "surr_p50": df.surr_T_G_p50,
        "surr_p75": df.surr_T_G_p75, "surr_p95": df.surr_T_G_p95,
        "obs_z": df.obs_z, "obs_p_one_sided": df.obs_p_one_sided_upper,
        "n_surrogates": df.n_surrogates,
        "n_swaps_per_surrogate": df.n_swaps_per_surrogate,
        "surrogate_source": surrogate_source, "seed": seed, "defined": True,
    })
    return out[_SCHEMA]


# ---------------------------------------------------------------------------
# epi_only Grassmann (fresh; reuses the epi_only eigendecomp cache)
# ---------------------------------------------------------------------------
def epi_only_rows(pat: str, band: str, n_surr: int, swap_factor: int,
                  verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc_full(pat, ph, band) for ph in es.PHASES_3}
    except Exception as e:
        if verbose:
            print(f"[audit_78] SKIP {pat}/{band} epi_only: {e}")
        return []
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        return []
    epi = es.epi_mask_for(pat, n_expected=N)
    keep = es.node_mask_for_config(epi, "epi_only")
    if keep is None:                      # Pat_15 or n_epi < MIN_NODES
        return []
    idx = np.ix_(keep, keep)
    Wk = {ph: np.ascontiguousarray(Ws[ph][idx]) for ph in es.PHASES_3}
    Nk = int(keep.sum())
    kgrid = [k for k in K_GRID if k + 1 <= Nk]
    if not kgrid:
        return []
    kmax = max(kgrid)

    U_obs = {ph: topk_basis(laplacian_eig(Wk[ph])[1], kmax) for ph in es.PHASES_3}
    obs_TG = t_g_at_all_k(U_obs["rest_pre_A"], U_obs["task_test"],
                          U_obs["rest_post"], kgrid)

    surr_evecs = {}
    for ph in es.PHASES_3:
        _, evecs = load_or_compute_eigs_at_path(
            es.surr_eig_path("epi_only", pat, band, ph, n_surr, swap_factor),
            Wk[ph], n_surr, swap_factor,
            es.cell_rng(pat, band, ph, "epi_only"), verbose=verbose)
        surr_evecs[ph] = evecs
    surr_TG = np.full((n_surr, len(kgrid)), np.nan)
    for r in range(n_surr):
        U = {ph: topk_basis(surr_evecs[ph][r], kmax) for ph in es.PHASES_3}
        if any(np.isnan(U[ph]).any() for ph in es.PHASES_3):
            continue
        surr_TG[r] = t_g_at_all_k(U["rest_pre_A"], U["task_test"],
                                  U["rest_post"], kgrid)

    n_swaps = swap_factor * (Nk * (Nk - 1)) // 2
    rows = []
    for i, k in enumerate(kgrid):
        s = surr_TG[:, i]
        s = s[np.isfinite(s)]
        obs = float(obs_TG[i])
        if s.size > 0:
            smean = float(np.mean(s)); sstd = float(np.std(s, ddof=1))
            qs = np.quantile(s, [0.05, 0.25, 0.50, 0.75, 0.95])
            z = (obs - smean) / sstd if sstd > 0 else float("nan")
            p_up = float(np.mean(s >= obs))
        else:
            smean = sstd = float("nan"); qs = [float("nan")] * 5
            z = p_up = float("nan")
        rows.append({
            "patient": pat, "band": band, "config": "epi_only",
            "config_class": "subgraph", "substrate": SUBSTRATE, "k": int(k),
            "n_nodes": Nk, "n_pairs": np.nan, "obs_stat": obs,
            "surr_mean": smean, "surr_std": sstd, "surr_p5": float(qs[0]),
            "surr_p25": float(qs[1]), "surr_p50": float(qs[2]),
            "surr_p75": float(qs[3]), "surr_p95": float(qs[4]),
            "obs_z": z, "obs_p_one_sided": p_up, "n_surrogates": int(s.size),
            "n_swaps_per_surrogate": int(n_swaps),
            "surrogate_source": "submatrix", "seed": es.config_seed("epi_only"),
            "defined": True,
        })
    return rows


# ---------------------------------------------------------------------------
# Cohort summary (per band, config, k)
# ---------------------------------------------------------------------------
def _band_order(df: pd.DataFrame) -> list[str]:
    present = set(df.band.unique())
    return [b for b in es.ALL_BANDS if b in present]


def cohort_summary(per_pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    full = per_pat[per_pat.config == "full"]
    full_p = {}
    for (band, k), d in full.groupby(["band", "k"]):
        d = d[d.defined]
        n_above = int((d.obs_p_one_sided < 0.05).sum())
        full_p[(band, k)] = es.cohort_verdict(
            d.obs_stat.values, d.surr_p50.values, n_above)["wilcoxon_p"]
    for band in _band_order(per_pat):
        for config in [c for c in es.SUBGRAPH_CONFIGS
                       if c in set(per_pat[per_pat.band == band].config)]:
            sub = per_pat[(per_pat.band == band) & (per_pat.config == config)]
            for k in sorted(sub.k.unique()):
                cell = sub[sub.k == k]
                d = cell[cell.defined]
                n_above = int((d.obs_p_one_sided < 0.05).sum())
                v = es.cohort_verdict(d.obs_stat.values, d.surr_p50.values,
                                      n_above)
                d2 = d[d.patient != "Pat_13"]
                n_above2 = int((d2.obs_p_one_sided < 0.05).sum())
                v2 = es.cohort_verdict(d2.obs_stat.values, d2.surr_p50.values,
                                       n_above2)
                flag = ("baseline" if config == "full"
                        else es.sensitivity_flag(
                            full_p.get((band, int(k)), np.nan),
                            v["wilcoxon_p"],
                            cfg_defined=(v["verdict"] != "undefined")))
                out.append({
                    "band": band, "config": config, "substrate": SUBSTRATE,
                    "k": int(k), "n_defined": v["n_defined"],
                    "obs_median": v["med_obs"],
                    "surr_median_median": v["med_surr"],
                    "n_above_own_surrogate": n_above,
                    "paired_wilcoxon_z": v["wilcoxon_z"],
                    "paired_wilcoxon_p": v["wilcoxon_p"],
                    "lopat13_wilcoxon_p": v2["wilcoxon_p"],
                    "verdict": v["verdict"], "sensitivity_flag": flag,
                })
    return pd.DataFrame(out)


def _kspan(cohort_band_cfg: pd.DataFrame) -> tuple[int, int, int]:
    """(#k separated, longest separated run, max k with ≥8 defined)."""
    sub = cohort_band_cfg.sort_values("k")
    sep = (sub.verdict == "separated").values
    n_sep = int(sep.sum())
    longest = cur = 0
    for v in sep:
        cur = cur + 1 if v else 0
        longest = max(longest, cur)
    pw = sub[sub.n_defined >= 8]
    kmax8 = int(pw.k.max()) if len(pw) else 0
    return n_sep, longest, kmax8


def write_readme(per_pat: pd.DataFrame, cohort: pd.DataFrame,
                 n_surr: int, runtime_s: float) -> None:
    L: list[str] = []
    a = L.append
    a("---")
    a("name: epi_stratified_grassmann")
    a("scope: epi_node_stratified_grassmann_T_G_trace_all_bands")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: first_pass")
    a("build_script: scripts/01_compute/audit/audit_78_epi_stratified_grassmann.py")
    a("---")
    a("")
    a("# Epi-node-stratified Grassmann T_G(k) trace (all bands, n=10)")
    a("")
    a("**Head.** Chordal subspace trace T_G(k) (positive = TRACE) across the "
      "three subgraph views. `full` re-emitted from audit_66, `exclude_epi` "
      "from audit_67, `epi_only` computed fresh. Pair-class is cophenetic-only "
      "(audit_77): the Grassmann statistic has no pair-restriction analogue. "
      "`epi_only` is exploratory — k-coverage collapses past the ≥8-patient "
      "frontier (n_epi∈{6..30}, Pat_15=0).")
    a("")
    a("## k-span summary — band × config")
    a("")
    a("| band | config | #k separated | longest run | max k (≥8 pat) |")
    a("|---|---|---|---|---|")
    for band in _band_order(per_pat):
        for config in es.SUBGRAPH_CONFIGS:
            cc = cohort[(cohort.band == band) & (cohort.config == config)]
            if cc.empty:
                continue
            n_sep, longest, kmax8 = _kspan(cc)
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {config} "
              f"| {n_sep} | {longest} | {kmax8} |")
    a("")
    a("## Caveats")
    a("- **epi_only k-coverage**: ≥8 patients only at low k (need n_epi>k; "
      "Pat_15 absent). High-k epi_only cells are Pat_13-dominated and "
      "underpowered — not trustworthy.")
    a("- exclude_epi k≤88 (min N_reduced=89, Pat_13); full k≤112. k axes are "
      "approximately but not exactly comparable across configs (varying N).")
    a("")
    a("## Provenance")
    a(f"- N_surrogates = {n_surr}; phases " + ", ".join(es.PHASES_3))
    a("- full: audit_66 (seed 20260511); exclude_epi: audit_67 (seed 20260514); "
      "epi_only: fresh (seed 20260511, reuses audit_77 epi_only eig cache).")
    a(f"- Wall-clock (epi_only compute + merge): {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `grassmann_per_patient_per_k.csv` — tidy long-format")
    a("- `grassmann_cohort.csv` — per (band, config, k) cohort verdict + sensitivity")
    (OUT / "README_grassmann.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(es.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--n-surrogates", type=int, default=es.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=es.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if not A66.exists() or not A67.exists():
        raise SystemExit(f"[audit_78] need prior CSVs:\n  {A66}\n  {A67}")

    frames = [
        _map_prior(pd.read_csv(A66), "full", "N_nodes", "fullgraph", 20260511),
        _map_prior(pd.read_csv(A67), "exclude_epi", "N_reduced", "submatrix",
                   20260514),
    ]
    frames = [f[f.band.isin(args.bands) & f.patient.isin(args.patients)]
              for f in frames]

    print("[audit_78] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    epi_rows: list[dict] = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            rows = epi_only_rows(pat, band, args.n_surrogates,
                                 args.swap_factor, args.verbose)
            epi_rows.extend(rows)
            if rows:
                ks = [r["k"] for r in rows]
                print(f"[audit_78] {pat}/{band} epi_only: Nk={rows[0]['n_nodes']} "
                      f"k={min(ks)}..{max(ks)} ({time.time()-tc:.1f}s)")
    runtime = time.time() - t0

    if epi_rows:
        frames.append(pd.DataFrame(epi_rows)[_SCHEMA])
    per_pat = pd.concat(frames, ignore_index=True)
    per_pat.to_csv(OUT / "grassmann_per_patient_per_k.csv", index=False)

    cohort = cohort_summary(per_pat)
    cohort.to_csv(OUT / "grassmann_cohort.csv", index=False)
    write_readme(per_pat, cohort, args.n_surrogates, runtime)
    print(f"[audit_78] {len(per_pat)} rows ({len(epi_rows)} epi_only), "
          f"{runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
