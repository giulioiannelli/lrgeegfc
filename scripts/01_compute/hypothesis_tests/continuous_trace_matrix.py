#!/usr/bin/env python3
"""Continuous-trace matrix — per-cell distance-shift summary.

Per (patient, band) under ``imcoh_abs`` computes:
    D_φ           = ultrametric distance matrices (φ ∈ {pre, learn, test, post})
    Δ_task(i,j)   = D_test(i,j) − D_pre(i,j)
    Δ_rest(i,j)   = D_post(i,j) − D_pre(i,j)
    σ(i,j)        = sgn(Δ_task) · sgn(Δ_rest)            ∈ {−1, 0, +1}
    ρ             = Spearman(Δ_task, Δ_rest) on upper-triangle pairs

Run modes
---------
* default (--mode shared): both Δ vectors share the same D_pre baseline.
  This is the descriptive companion to H2c. ρ contains a contribution
  from the shared baseline noise term (any noise in D_pre appears in
  both Δ_task and Δ_rest with the same sign).

* --mode split-baseline (Run A control): split rest_pre time series in
  halves, compute D_pre_A and D_pre_B from the halves cache (already
  populated by h2e_split_half.py). Then::

      Δ_task = D_test − D_pre_A
      Δ_rest = D_post − D_pre_B

  D_pre_A ⊥ D_pre_B by construction (independent estimators on disjoint
  samples), so the shared-baseline-noise contribution to ρ is removed.
  Population baseline still cancels in expectation. The remaining ρ
  reflects only "task-induced shift in test that persists in post".

Outputs (suffix encodes mode):
    data/reports/imcoh_continuous_trace/per_cell_summary[_split].csv
    data/reports/imcoh_continuous_trace/per_pair[_split]/{patient}_{band}.npz
        (per-cell Δ_task, Δ_rest, σ vectors; loaded by the figure script)

Scope: ``.agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md``
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE
from lrg_eegfc.config.paths import CACHE_ROOT, IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result


HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"


def load_D_square(pat: str, phase: str, band: str,
                  cache_root=IMCOH_LRG_CACHE) -> np.ndarray | None:
    """Load LRG ultrametric matrix, return as N×N square form."""
    try:
        r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs",
                            cache_root=cache_root)
    except Exception:
        return None
    if r is None:
        return None
    um = np.asarray(r.ultrametric_matrix)
    if um.ndim == 1:
        return squareform(um)
    return um


def _per_cell_shared(pat: str, band: str) -> dict | None:
    """Original mode: shared D_pre baseline."""
    Ds = {phase: load_D_square(pat, phase, band)
          for phase in ["rest_pre", "task_test", "rest_post"]}
    if any(D is None for D in Ds.values()):
        return None
    N = Ds["rest_pre"].shape[0]
    if not all(D.shape == (N, N) for D in Ds.values()):
        return None

    iu = np.triu_indices(N, k=1)
    dD_task = (Ds["task_test"] - Ds["rest_pre"])[iu]
    dD_rest = (Ds["rest_post"] - Ds["rest_pre"])[iu]
    dmax_pre = float(Ds["rest_pre"][iu].max())
    return _finalize(pat, band, N, iu, dD_task, dD_rest, dmax_pre)


def _per_cell_split(pat: str, band: str) -> dict | None:
    """Run A: independent half-baselines.

    Δ_task = D_test − D_pre_A,  Δ_rest = D_post − D_pre_B,
    where D_pre_A and D_pre_B are LRG ultrametric matrices computed on
    disjoint halves of the rest_pre time series (cache populated by
    h2e_split_half.py at ``data/cache/imcoh_lrg_halves/``).
    """
    D_test = load_D_square(pat, "task_test", band, cache_root=IMCOH_LRG_CACHE)
    D_post = load_D_square(pat, "rest_post", band, cache_root=IMCOH_LRG_CACHE)
    D_pre_A = load_D_square(pat, "rest_pre_A", band, cache_root=HALVES_CACHE)
    D_pre_B = load_D_square(pat, "rest_pre_B", band, cache_root=HALVES_CACHE)
    if any(D is None for D in (D_test, D_post, D_pre_A, D_pre_B)):
        return None
    N = D_test.shape[0]
    if not all(D.shape == (N, N) for D in (D_post, D_pre_A, D_pre_B)):
        return None

    iu = np.triu_indices(N, k=1)
    dD_task = (D_test - D_pre_A)[iu]
    dD_rest = (D_post - D_pre_B)[iu]
    # dmax_pre uses the average of A and B for cohort-level normalisation
    dmax_pre = float(0.5 * (D_pre_A[iu].max() + D_pre_B[iu].max()))
    return _finalize(pat, band, N, iu, dD_task, dD_rest, dmax_pre)


def _finalize(pat: str, band: str, N: int, iu, dD_task, dD_rest,
              dmax_pre: float) -> dict:

    rho = float(spearmanr(dD_task, dD_rest).statistic)
    sigma = np.sign(dD_task) * np.sign(dD_rest)
    frac_pos = float((sigma > 0).mean())
    frac_neg = float((sigma < 0).mean())
    frac_zero = float((sigma == 0).mean())
    m_p = len(dD_task)

    return {
        "summary": {
            "patient": pat, "band": band,
            "N_p": N, "m_p": m_p,
            "rho": rho,
            "frac_pos_sigma": frac_pos,
            "frac_neg_sigma": frac_neg,
            "frac_zero_sigma": frac_zero,
            "dmax_pre": dmax_pre,
        },
        "vectors": {
            "dD_task": dD_task,
            "dD_rest": dD_rest,
            "sigma": sigma.astype(np.int8),
            "iu_i": iu[0].astype(np.int32),
            "iu_j": iu[1].astype(np.int32),
        },
    }


def per_cell(pat: str, band: str, mode: str = "shared") -> dict | None:
    """Dispatch to mode-specific per-cell computation."""
    if mode == "shared":
        return _per_cell_shared(pat, band)
    if mode == "split-baseline":
        return _per_cell_split(pat, band)
    raise ValueError(f"Unknown mode {mode!r}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", default=None)
    ap.add_argument("--bands", default=None)
    ap.add_argument("--mode", choices=["shared", "split-baseline"],
                    default="shared",
                    help="shared: original Δ_task/Δ_rest sharing D_pre. "
                         "split-baseline: Run A control with independent "
                         "D_pre_A / D_pre_B from rest_pre halves.")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    patients = (args.patients.split(",") if args.patients
                else list(PATIENTS_4PHASE))
    bands = args.bands.split(",") if args.bands else list(BRAIN_BANDS_NAMES)

    suffix = "" if args.mode == "shared" else "_split"
    out_dir = REPORTS_ROOT / "imcoh_continuous_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    pair_dir = out_dir / f"per_pair{suffix}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    summaries: list[dict] = []
    skipped = []
    for pat in patients:
        for band in bands:
            r = per_cell(pat, band, mode=args.mode)
            if r is None:
                skipped.append((pat, band))
                if args.verbose:
                    print(f"[ctm] SKIP {pat:7s} {band:10s} — incomplete cache")
                continue
            summaries.append(r["summary"])
            vecs = r["vectors"]
            np.savez_compressed(
                pair_dir / f"{pat}_{band}.npz",
                **vecs,
            )
            if args.verbose:
                s = r["summary"]
                print(f"[ctm] {pat:7s} {band:10s}  N={s['N_p']:3d}  "
                      f"m={s['m_p']:5d}  ρ={s['rho']:+.3f}  "
                      f"σ_pos={s['frac_pos_sigma']:.2f}  "
                      f"σ_neg={s['frac_neg_sigma']:.2f}")

    df = pd.DataFrame(summaries)
    csv_path = out_dir / f"per_cell_summary{suffix}.csv"
    df.to_csv(csv_path, index=False)
    print()
    print(f"[ctm] mode={args.mode}  wrote {csv_path}  ({len(df)} cells)")
    if skipped:
        print(f"[ctm] skipped: {skipped}")
    print()
    print("[ctm] per-band cohort summary:")
    print(f"{'band':<10s}  {'n':>3s}  {'med_ρ':>7s}  {'mean_ρ':>7s}  "
          f"{'n_ρ>0':>5s}  {'med_σ_pos':>10s}")
    for band in bands:
        sub = df[df["band"] == band]
        if sub.empty:
            print(f"{band:<10s}  0")
            continue
        n = len(sub)
        med_r = sub["rho"].median()
        mn_r  = sub["rho"].mean()
        n_pos = int((sub["rho"] > 0).sum())
        med_sp = sub["frac_pos_sigma"].median()
        print(f"{band:<10s}  {n:>3d}  {med_r:>+.3f}  {mn_r:>+.3f}  "
              f"{n_pos}/{n}    {med_sp:>+.3f}")


if __name__ == "__main__":
    main()
