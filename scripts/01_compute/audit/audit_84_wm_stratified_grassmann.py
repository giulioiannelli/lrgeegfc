#!/usr/bin/env python3
"""Audit 84 — white-matter-node-stratified Grassmann T_G(k) trace.

Grassmann companion to audit_83 (white-matter analogue of audit_78). The
chordal subspace statistic
T_G(k) = d_chord(rest_pre, task; k) − d_chord(task, rest_post; k)
(positive = TRACE) lives on the whole k-dimensional slow-mode subspace of the
Laplacian, so it admits only the SUBGRAPH views (full / exclude_wm / wm_only) —
there is no pair-restriction analogue, so the pair-class decomposition is
cophenetic/raw only (audit_83).

To avoid recomputing what is already validated:
- ``full``       is re-emitted from audit_66 (grassmann_matched_strength_surrogate).
- ``exclude_wm`` and ``wm_only`` are computed fresh here, REUSING the
  matched-strength eigendecomposition caches that audit_83 builds on the WM
  submatrices (``wmX`` / ``wmONLY``). Run audit_83 first for guaranteed cache
  hits; otherwise this script generates them (deterministic via the shared
  ``_wm_stratify.cell_rng``).

Critical preamble (per CLAUDE.md rule)
======================================
(1) **Claim.** The Grassmann subspace trace (β, γ_l load-bearing in audit_66)
    survives white-matter exclusion (exclude_wm) and is present/absent within
    the white-matter subgraph alone (wm_only).
(2) **Null.** R=200 matched-strength 4-cycle ±δ on the relevant graph (full;
    W[gray]; W[wm]).
(3) **Strongest alternative.** exclude_wm/wm_only change N (30–57 % of nodes are
    WM), so the slow-mode subspace dimension available shrinks and the k-axis is
    not exactly comparable across configs; a T_G change could be subspace-size,
    not tissue. Reported via per-config N and the ≥8-patient k-frontier.
(4) **Mechanical reach.** The null fixes per-node strength, so a surviving T_G
    is not a strength artifact; it does NOT make k-axes comparable across N.
(5) **Falsification + limits.** "Trace lives in gray matter" → exclude_wm
    persists across the manuscript k-range. wm_only is the secondary arm; unlike
    epi_only it keeps 36–67 nodes so its k-coverage reaches the cohort frontier
    for moderate k. Pat_15 (most WM, β-anti) → leave-Pat_15-out p reported.

Outputs
-------
``data/audit/wm_stratified/``
    grassmann_per_patient_per_k.csv   tidy long (one row per pat×band×config×k)
    grassmann_cohort.csv              per (band, config, k) cohort verdict
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
    ensure_half_fcs,
    laplacian_eig,
    load_phase_fc_full,
    t_g_at_all_k,
    topk_basis,
)
import _wm_stratify as ws  # type: ignore


SUBSTRATE = "grassmann"
OUT = ROOT / "data" / "audit" / "wm_stratified"
OUT.mkdir(parents=True, exist_ok=True)
A66 = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate" \
    / "per_patient_per_band_per_k.csv"

LO_PATIENT = "Pat_15"

_SCHEMA = ["patient", "band", "config", "config_class", "substrate", "k",
           "n_nodes", "n_pairs", "obs_stat", "surr_mean", "surr_std",
           "surr_p5", "surr_p25", "surr_p50", "surr_p75", "surr_p95",
           "obs_z", "obs_p_one_sided", "n_surrogates",
           "n_swaps_per_surrogate", "surrogate_source", "seed", "defined"]


# ---------------------------------------------------------------------------
# Re-emit the audit_66 full-graph Grassmann CSV into the unified schema
# ---------------------------------------------------------------------------
def _map_prior(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame({
        "patient": df.patient, "band": df.band, "config": "full",
        "config_class": "subgraph", "substrate": SUBSTRATE, "k": df.k,
        "n_nodes": df.N_nodes, "n_pairs": np.nan,
        "obs_stat": df.obs_T_G, "surr_mean": df.surr_T_G_mean,
        "surr_std": df.surr_T_G_std, "surr_p5": df.surr_T_G_p5,
        "surr_p25": df.surr_T_G_p25, "surr_p50": df.surr_T_G_p50,
        "surr_p75": df.surr_T_G_p75, "surr_p95": df.surr_T_G_p95,
        "obs_z": df.obs_z, "obs_p_one_sided": df.obs_p_one_sided_upper,
        "n_surrogates": df.n_surrogates,
        "n_swaps_per_surrogate": df.n_swaps_per_surrogate,
        "surrogate_source": "fullgraph", "seed": 20260511, "defined": True,
    })
    return out[_SCHEMA]


# ---------------------------------------------------------------------------
# Fresh subgraph Grassmann (exclude_wm / wm_only); reuses audit_83 eig caches
# ---------------------------------------------------------------------------
def subgraph_rows(pat: str, band: str, config: str, n_surr: int,
                  swap_factor: int, verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc_full(pat, ph, band) for ph in ws.PHASES_3}
    except Exception as e:
        if verbose:
            print(f"[audit_84] SKIP {pat}/{band} {config}: {e}")
        return []
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        return []
    wm = ws.wm_mask_for(pat, n_expected=N)
    keep = ws.node_mask_for_config(wm, config)
    if keep is None:
        return []
    idx = np.ix_(keep, keep)
    Wk = {ph: np.ascontiguousarray(Ws[ph][idx]) for ph in ws.PHASES_3}
    Nk = int(keep.sum())
    kgrid = [k for k in K_GRID if k + 1 <= Nk]
    if not kgrid:
        return []
    kmax = max(kgrid)

    U_obs = {ph: topk_basis(laplacian_eig(Wk[ph])[1], kmax) for ph in ws.PHASES_3}
    obs_TG = t_g_at_all_k(U_obs["rest_pre_A"], U_obs["task_test"],
                          U_obs["rest_post"], kgrid)

    surr_evecs = {}
    for ph in ws.PHASES_3:
        _, evecs = load_or_compute_eigs_at_path(
            ws.surr_eig_path(config, pat, band, ph, n_surr, swap_factor),
            Wk[ph], n_surr, swap_factor,
            ws.cell_rng(pat, band, ph, config), verbose=verbose)
        surr_evecs[ph] = evecs
    surr_TG = np.full((n_surr, len(kgrid)), np.nan)
    for r in range(n_surr):
        U = {ph: topk_basis(surr_evecs[ph][r], kmax) for ph in ws.PHASES_3}
        if any(np.isnan(U[ph]).any() for ph in ws.PHASES_3):
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
            "patient": pat, "band": band, "config": config,
            "config_class": "subgraph", "substrate": SUBSTRATE, "k": int(k),
            "n_nodes": Nk, "n_pairs": np.nan, "obs_stat": obs,
            "surr_mean": smean, "surr_std": sstd, "surr_p5": float(qs[0]),
            "surr_p25": float(qs[1]), "surr_p50": float(qs[2]),
            "surr_p75": float(qs[3]), "surr_p95": float(qs[4]),
            "obs_z": z, "obs_p_one_sided": p_up, "n_surrogates": int(s.size),
            "n_swaps_per_surrogate": int(n_swaps),
            "surrogate_source": "submatrix", "seed": ws.config_seed(config),
            "defined": True,
        })
    return rows


# ---------------------------------------------------------------------------
# Cohort summary (per band, config, k)
# ---------------------------------------------------------------------------
def _band_order(df: pd.DataFrame) -> list[str]:
    present = set(df.band.unique())
    return [b for b in ws.ALL_BANDS if b in present]


def cohort_summary(per_pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    full = per_pat[per_pat.config == "full"]
    full_p = {}
    for (band, k), d in full.groupby(["band", "k"]):
        d = d[d.defined]
        n_above = int((d.obs_p_one_sided < 0.05).sum())
        full_p[(band, k)] = ws.cohort_verdict(
            d.obs_stat.values, d.surr_p50.values, n_above)["wilcoxon_p"]
    for band in _band_order(per_pat):
        for config in [c for c in ws.SUBGRAPH_CONFIGS
                       if c in set(per_pat[per_pat.band == band].config)]:
            sub = per_pat[(per_pat.band == band) & (per_pat.config == config)]
            for k in sorted(sub.k.unique()):
                cell = sub[sub.k == k]
                d = cell[cell.defined]
                n_above = int((d.obs_p_one_sided < 0.05).sum())
                v = ws.cohort_verdict(d.obs_stat.values, d.surr_p50.values, n_above)
                d2 = d[d.patient != LO_PATIENT]
                n_above2 = int((d2.obs_p_one_sided < 0.05).sum())
                v2 = ws.cohort_verdict(d2.obs_stat.values, d2.surr_p50.values, n_above2)
                flag = ("baseline" if config == "full"
                        else ws.sensitivity_flag(
                            full_p.get((band, int(k)), np.nan), v["wilcoxon_p"],
                            cfg_defined=(v["verdict"] != "undefined")))
                out.append({
                    "band": band, "config": config, "substrate": SUBSTRATE,
                    "k": int(k), "n_defined": v["n_defined"],
                    "obs_median": v["med_obs"],
                    "surr_median_median": v["med_surr"],
                    "n_above_own_surrogate": n_above,
                    "paired_wilcoxon_z": v["wilcoxon_z"],
                    "paired_wilcoxon_p": v["wilcoxon_p"],
                    "lopat15_wilcoxon_p": v2["wilcoxon_p"],
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
    a("name: wm_stratified_grassmann")
    a("scope: white_matter_node_stratified_grassmann_T_G_trace_all_bands")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: first_pass")
    a("build_script: scripts/01_compute/audit/audit_84_wm_stratified_grassmann.py")
    a("---")
    a("")
    a("# White-matter-node-stratified Grassmann T_G(k) trace (all bands, n=10)")
    a("")
    a("**Head.** Chordal subspace trace T_G(k) (positive = TRACE) across the "
      "three subgraph views. `full` re-emitted from audit_66; `exclude_wm` and "
      "`wm_only` computed fresh (reusing the audit_83 WM-submatrix eig caches). "
      "Pair-class is cophenetic/raw only (audit_83): the Grassmann statistic has "
      "no pair-restriction analogue. WM is 30–57 % of nodes, so exclude_wm/"
      "wm_only change the subspace dimension and the k-axis is only approximately "
      "comparable across configs.")
    a("")
    a("## k-span summary — band × config")
    a("")
    a("| band | config | #k separated | longest run | max k (≥8 pat) |")
    a("|---|---|---|---|---|")
    for band in _band_order(per_pat):
        for config in ws.SUBGRAPH_CONFIGS:
            cc = cohort[(cohort.band == band) & (cohort.config == config)]
            if cc.empty:
                continue
            n_sep, longest, kmax8 = _kspan(cc)
            a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {config} "
              f"| {n_sep} | {longest} | {kmax8} |")
    a("")
    a("## Caveats")
    a(f"- **{LO_PATIENT}** (most WM, β-LRG anti) → LO-P15 column is the "
      "leave-one-out robustness line.")
    a("- exclude_wm/wm_only change N → k-axes only approximately comparable; "
      "read max k (≥8 pat) per config.")
    a("")
    a("## Provenance")
    a(f"- N_surrogates = {n_surr}; phases " + ", ".join(ws.PHASES_3))
    a("- full: audit_66 (seed 20260511); exclude_wm/wm_only: fresh "
      "(seed 20260608, reuse audit_83 WM eig caches).")
    a("- WM = dominant Desikan-Killiany tissue == 'Wm' (atlas argmax).")
    a(f"- Wall-clock (subgraph compute + merge): {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `grassmann_per_patient_per_k.csv` — tidy long-format")
    a("- `grassmann_cohort.csv` — per (band, config, k) cohort verdict")
    (OUT / "README_grassmann.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(ws.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(ws.COHORT))
    ap.add_argument("--configs", nargs="+",
                    default=["exclude_wm", "wm_only"])
    ap.add_argument("--n-surrogates", type=int, default=ws.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=ws.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if not A66.exists():
        raise SystemExit(f"[audit_84] need prior full-graph CSV: {A66}")

    full = _map_prior(pd.read_csv(A66))
    full = full[full.band.isin(args.bands) & full.patient.isin(args.patients)]
    frames = [full]

    print("[audit_84] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    fresh: list[dict] = []
    for config in [c for c in ("exclude_wm", "wm_only") if c in args.configs]:
        for band in args.bands:
            for pat in args.patients:
                tc = time.time()
                rows = subgraph_rows(pat, band, config, args.n_surrogates,
                                     args.swap_factor, args.verbose)
                fresh.extend(rows)
                if rows:
                    ks = [r["k"] for r in rows]
                    print(f"[audit_84] {pat}/{band} {config}: Nk={rows[0]['n_nodes']} "
                          f"k={min(ks)}..{max(ks)} ({time.time()-tc:.1f}s)")
    runtime = time.time() - t0

    if fresh:
        frames.append(pd.DataFrame(fresh)[_SCHEMA])
    per_pat = pd.concat(frames, ignore_index=True)
    per_pat.to_csv(OUT / "grassmann_per_patient_per_k.csv", index=False)

    cohort = cohort_summary(per_pat)
    cohort.to_csv(OUT / "grassmann_cohort.csv", index=False)
    write_readme(per_pat, cohort, args.n_surrogates, runtime)
    print(f"[audit_84] {len(per_pat)} rows ({len(fresh)} fresh subgraph), "
          f"{runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
