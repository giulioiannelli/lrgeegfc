#!/usr/bin/env python3
"""What is a cohort-median scale slope actually worth?

A scale-dependence claim is routinely read off a *cohort-median profile*: take
the median margin across patients at each scale, correlate it with log s, quote
the Spearman rho. Lane E's per-scale grid, and the lead this lane was reopened
on, are both stated that way -- rho(log s) = +0.87, -0.83, +0.64.

Two things make that number unreadable as evidence, and this script measures
both instead of asserting them.

1. **The median suppresses the per-patient noise that the axis is made of.**
   The cohort-median curve is smoother than any individual patient's curve, so
   its Spearman against log s is systematically larger in absolute value than
   the typical patient's. The inflation factor is measured directly, per band
   and per object, as ``|rho(median profile)| / mean_k |rho(patient k)|``.
2. **It has almost no degrees of freedom.** Part A measured this scale axis as
   worth roughly one independent test, so a 16-point Spearman is not a
   16-point Spearman. Rather than argue about the right df, the null
   distribution of the cohort-median rho is constructed directly: promote a
   held-out matched-strength realization to the observed slot, form the same
   cohort-median profile, take the same Spearman. The width of that
   distribution is the answer.

The output is one table per (object, band): the observed cohort-median rho, its
null median and its two-sided tail probability, the mean per-patient rho with
its own null, and the inflation factor between them. Nothing is tested against
zero; every reference is the same statistic on the object's own surrogates.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from w0s_07_contrast_scale_structure import (                  # noqa: E402
    BANDS, ENC, FunctionalCells, _per_patient_slopes,
)

OUT = Path(os.environ.get(
    "W0S_OUT_CONTRAST",
    ROOT / "data" / "paper_final" / "lane_s_scale" / "contrast"))
N_DRAWS = int(os.environ.get("W0S_NDRAWS", "200"))


def _median_slope(M: np.ndarray, ls: np.ndarray) -> float:
    """Spearman between the cohort-median profile and ``log s`` -- the quoted number."""
    med = np.nanmedian(M, axis=0)
    ok = np.isfinite(med)
    if ok.sum() < 4:
        return np.nan
    return float(spearmanr(med[ok], ls[ok])[0])


def _two_sided(obs: float, null: np.ndarray) -> float:
    n = np.asarray(null, float)
    n = n[np.isfinite(n)]
    if not n.size or not np.isfinite(obs):
        return np.nan
    return float((1 + int((np.abs(n) >= abs(obs)).sum())) / (n.size + 1))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cs = FunctionalCells(ENC / "grid", PATIENTS_4PHASE, BANDS)
    s = np.asarray(cs.s, float)
    ls = np.log(s)
    print(f"[w0s-medslope] {len(cs.readouts)} objects x {len(cs.bands)} bands "
          f"| {cs.n_scales} scales | {N_DRAWS} held-out draws", flush=True)

    rows = []
    for m in cs.readouts:
        for b in cs.bands:
            if not cs.have(b):
                continue
            M, _ = cs.margins(b, m)
            obs_med = _median_slope(M, ls)
            pp = _per_patient_slopes(M, s)
            obs_pat = float(np.nanmean(np.abs(pp)))
            n = min(N_DRAWS, cs.R)
            nm = np.full(n, np.nan)
            npat = np.full(n, np.nan)
            for i in range(n):
                Mn = cs.heldout_margins(b, m, i)
                nm[i] = _median_slope(Mn, ls)
                npat[i] = float(np.nanmean(np.abs(_per_patient_slopes(Mn, s))))
            rows.append(dict(
                object=m, band=b,
                rho_cohort_median=obs_med,
                rho_median_null_med=float(np.nanmedian(nm)),
                rho_median_null_absmed=float(np.nanmedian(np.abs(nm))),
                rho_median_null_abs_p95=float(np.nanpercentile(np.abs(nm), 95)),
                p_rho_median=_two_sided(obs_med, nm),
                mean_abs_rho_patient=obs_pat,
                mean_abs_rho_patient_null=float(np.nanmedian(npat)),
                inflation=(abs(obs_med) / obs_pat if obs_pat else np.nan),
                inflation_null=float(np.nanmedian(np.abs(nm) / npat)),
                n_draws=int(np.isfinite(nm).sum())))
            print(f"  {m:22s} {b:11s} rho_med={obs_med:+.3f} "
                  f"null|rho|med={rows[-1]['rho_median_null_absmed']:.3f} "
                  f"p95={rows[-1]['rho_median_null_abs_p95']:.3f} "
                  f"p={rows[-1]['p_rho_median']:.3f} "
                  f"inflation={rows[-1]['inflation']:.2f}", flush=True)
    d = pd.DataFrame(rows)
    d["q_rho_median"] = np.nan
    for m in d.object.unique():
        sel = (d.object == m) & d.p_rho_median.notna()
        if sel.any():
            d.loc[sel, "q_rho_median"] = np.asarray(
                bh_fdr(d.loc[sel, "p_rho_median"].to_numpy()), float)
    d.to_csv(OUT / "cohort_median_slope_null.csv", index=False)
    print(f"\n  inflation factor: median {d.inflation.median():.2f}x observed, "
          f"{d.inflation_null.median():.2f}x under the null "
          f"(the median smooths noise whether or not there is signal)",
          flush=True)
    print(f"  cells with q < 0.05: {int((d.q_rho_median < 0.05).sum())}/{len(d)}",
          flush=True)
    print(f"[w0s-medslope] wrote {OUT/'cohort_median_slope_null.csv'}", flush=True)


if __name__ == "__main__":
    main()
