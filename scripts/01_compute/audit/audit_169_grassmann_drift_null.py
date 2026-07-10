#!/usr/bin/env python3
"""audit_169 — within-recording DRIFT null for the GRASSMANN subspace trace.

Fills the gap flagged 2026-07-10: the drift null (audit_167 / C2) was run on raw
edges + cophenetic, but NEVER on the Grassmann subspace probe. The talk's
dual-null figure ("a band clears only if it overcomes BOTH matched-strength AND
drift") needs the Grassmann drift arm. Mirrors audit_167's construction exactly,
substituting the Grassmann whole-task T_G for the cophenetic/raw arc.

5-point preamble
1. Claim: the Grassmann subspace trace T_G (whole-task, beta/low_gamma/delta clear
   matched-strength) is task-driven, not session drift.
2. Null (drift): window a SINGLE pre-task rest_pre into 5 ordered segments, map
   pre=w0, task_test=w3, rest_post=w4, and compute the identical T_G through them.
   rest_pre is BEFORE any task, so a positive T_G here is drift by construction.
3. Strongest alt this addresses: slow within-session nonstationarity rotates the
   leading subspace monotonically, so a later-vs-earlier subspace contrast is
   positive with no task involved.
4. Reach / limits: whole-task T_G is a SIMPLE (non-conditional) contrast, so a
   windowed drift null is valid for it (unlike the retracted inference-specific
   partial correlation, audit_168). SNR caveat: 5 sub-windows are ~1/5 duration ->
   noisier FC -> T_drift biased LOW -> conservative for a positive drift finding
   (a STRICTER bar than the full-duration C2). Real baseline = rest_pre_A (1/2),
   drift baseline = w0 (1/5): baseline-duration mismatch, same caveat direction.
5. Falsify (of the task claim): if grassmann_drift reproduces grassmann_real
   (paired real>drift Wilcoxon fails), the Grassmann trace in that band is
   drift-driven. Expectation from coph (audit_167): beta clears, others fragile.

Real T_G per patient loaded from the canonical grassmann_inference_arc (functional
'onl', median-over-k); drift T_G recomputed here through the identical
laplacian_eig -> topk_basis -> t_g_at_all_k -> median-over-k pipeline.

Reads : grassmann_inference_arc/per_patient.csv (real, functional 'onl')
        raw timeseries (rest_pre) via load_timeseries
Writes: data/audit/grassmann_drift_null/{per_patient_per_band,cohort_summary}.csv
"""
from __future__ import annotations
import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import BRAIN_BANDS, nperseg_for_fs, FS_OVERRIDES
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_66_grassmann_matched_strength_surrogate import (  # type: ignore
    laplacian_eig, topk_basis, t_g_at_all_k,
)

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08",
          "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = list(BRAIN_BANDS)
NWIN = 5
KCAP = 113  # k grid 2..min(N, KCAP) -- matches the real grassmann pipeline cap
REAL_SRC = ROOT / "data/audit/grassmann_inference_arc/per_patient.csv"
OUTDIR = ROOT / "data/audit/grassmann_drift_null"


def band_abs(Coh, freqs, flo, fhi):
    m = (freqs >= flo) & (freqs <= fhi)
    A = np.abs(Coh[:, :, m]).mean(axis=-1)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return np.clip(A, 0.0, 1.0)


def tg_median(Wpre, Wtt, Wpost):
    """median-over-k whole-task T_G (T_G>0 = trace), identical collapse to real."""
    N = Wpre.shape[0]
    kgrid = list(range(2, min(N, KCAP)))
    kmax = kgrid[-1]
    U = lambda W: topk_basis(laplacian_eig(W)[1], kmax)
    tg = t_g_at_all_k(U(Wpre), U(Wtt), U(Wpost), kgrid)
    return float(np.median(tg))


def windows(pat):
    """5 imcoh cross-spectra from ordered windows of rest_pre (Welch once/window)."""
    X = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nper = nperseg_for_fs(fs)
    T = X.shape[1]
    out = []
    for i in range(NWIN):
        Xi = X[:, i * T // NWIN:(i + 1) * T // NWIN]
        fr, C = compute_msc_welch(Xi, fs, nperseg=min(nper, Xi.shape[1]),
                                  metric="imcoh")
        out.append((fr, C))
    return out


def load_real():
    df = pd.read_csv(REAL_SRC)
    df = df[df.functional == "onl"]
    return {(r.patient, r.band): float(r.median_obs_over_k) for r in df.itertuples()}


def main(limit=None):
    real = load_real()
    pats = COHORT[:limit] if limit else COHORT
    print(f"[audit_169] Grassmann within-rest_pre {NWIN}-window DRIFT null "
          f"({len(pats)} patients)\n")
    rows = []
    t0 = time.time()
    for pi, pat in enumerate(pats, 1):
        tp = time.time()
        wins = windows(pat)
        for band in BANDS:
            flo, fhi = BRAIN_BANDS[band]
            Wd = [band_abs(C, fr, flo, fhi) for (fr, C) in wins]
            drift = tg_median(Wd[0], Wd[3], Wd[4])  # pre=w0, task_test=w3, post=w4
            rows.append(dict(patient=pat, band=band,
                             real_TG=real.get((pat, band), np.nan),
                             drift_TG=drift))
        el = time.time() - tp
        eta = (time.time() - t0) / pi * (len(pats) - pi)
        print(f"  [{pi}/{len(pats)}] {pat}  {el:5.1f}s   ETA {eta:5.1f}s", flush=True)

    df = pd.DataFrame(rows)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTDIR / "per_patient_per_band.csv", index=False)

    # cohort drift gate: paired one-sided Wilcoxon real > drift
    crows = []
    print(f"\n{'band':<11}{'real_med':>9}{'drift_med':>10}{'drift_p':>9}  verdict")
    for band in BANDS:
        g = df[df.band == band].dropna(subset=["real_TG", "drift_TG"])
        rr, dd = g.real_TG.values, g.drift_TG.values
        try:
            p = float(wilcoxon(rr - dd, alternative="greater").pvalue)
        except ValueError:
            p = np.nan
        verdict = "clears_drift" if p < 0.05 else "drift"
        crows.append(dict(band=band, n=len(g), real_median=float(np.median(rr)),
                          drift_median=float(np.median(dd)), drift_gate_p=p,
                          verdict=verdict))
        star = "  <== clears" if p < 0.05 else ""
        print(f"  {band:<11}{np.median(rr):>9.3f}{np.median(dd):>10.3f}"
              f"{p:>9.3f}  {verdict}{star}")
    pd.DataFrame(crows).to_csv(OUTDIR / "cohort_summary.csv", index=False)
    print(f"\n[audit_169] total {time.time() - t0:.1f}s -> {OUTDIR}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=None,
                    help="run only the first N patients (timing/extrapolation)")
    args = ap.parse_args()
    main(limit=args.limit)
