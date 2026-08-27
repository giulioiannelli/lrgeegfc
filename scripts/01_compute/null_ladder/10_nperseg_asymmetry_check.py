#!/usr/bin/env python3
"""W0-B side-check: is the canonical rho_sym inflated by the A/B nperseg asymmetry?

THE ISSUE. In the production pipeline the two rest_pre halves A and B are
estimated at ``nperseg_for_fs(fs) // 2`` (audit_63 pre-flight, to keep the
segment count comparable on a half-length recording), while ``task_test`` and
``rest_post`` use the full ``nperseg_for_fs(fs)``. A different Welch segment
length means a different frequency resolution and therefore a different
estimation bias in the FC matrix.

Now look at where those matrices land:

    rho_sym = 1/2[ rho(D_task - D_A, D_post - D_B) + rho(D_task - D_B, D_post - D_A) ]

BOTH differences inside each Spearman have the form (full-nperseg estimate) minus
(half-nperseg estimate). Any systematic component of that estimator difference is
COMMON to the two arguments of the correlation, and a common additive component
inflates a correlation. If the effect is material, part of the canonical trace is
an artifact of the segmentation choice rather than a property of the phases.

This script measures it directly: recompute the OBSERVED rho_sym(s) with A and B
estimated at the SAME nperseg as task_test / rest_post, and compare against the
canonical mixed-nperseg construction, per patient and per scale. No surrogates —
this is about the observed statistic only.

5-point critical preamble
(1) CLAIM: the canonical rho_sym reflects cross-phase structure, not the A/B
    segmentation choice.
(2) NULL: matching nperseg across all four phases leaves rho_sym unchanged.
(3) STRONGEST ALTERNATIVE: a shared estimator-bias vector between the two
    difference terms manufactures correlation.
(4) CONTROLS FOR IT BY MECHANISM: matching nperseg removes the shared
    full-minus-half component while leaving the phases and their durations
    untouched, so any drop is attributable to the segmentation.
    CANNOT REJECT: matching nperseg also HALVES the number of segments in A and B
    relative to the canonical setting (a 4096-sample window on a ~300 s half gives
    about half as many windows as a 2048-sample one), so a drop mixes the removed
    bias with added estimator noise. A drop is therefore an upper bound on the
    artifact, not a point estimate of it; an absence of a drop is the clean
    result.
(5) FALSIFICATION: a large, consistent drop would mean the incumbent gate's
    observed statistic is partly a segmentation artifact and the A/B construction
    must be respecified before any lane uses it.

Outputs: data/paper_final/w0b_nulls/nperseg_check/per_patient_scale.csv
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))

from lrg_eegfc.config.const import (
    BRAIN_BANDS, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales
from lrg_eegfc.utils.surrogate.timeseries_nulls import (
    coherency_from_csd, csd_from_segment_subset, imcoh_abs_from_coherency, segment_fft,
)
from audit_150_rho_sym_gate import COHORT

PHASES = ("A", "B", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)
OUT = ROOT / "data" / "paper_final" / "w0b_nulls" / "nperseg_check"


def eig(W):
    return laplacian_eig(select_backbone(W, "mst020", frac=0.20))


def W_of(X, fs, nps, bnd):
    _, F, sc = segment_fft(X, fs, nps, band=bnd)
    W = imcoh_abs_from_coherency(coherency_from_csd(csd_from_segment_subset(F, None, sc)))
    del F
    return W


def per_cell(job):
    pat, band = job
    t0 = time.time()
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs); nps_h = max(256, nps // 2)
    bnd = BRAIN_BANDS[band]
    try:
        Xr = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), dtype=np.float32)
        if Xr.shape[0] > Xr.shape[1]:
            Xr = Xr.T
        T = Xr.shape[1]
        halves = {"A": Xr[:, :T // 2], "B": Xr[:, T // 2:]}
        W_mixed = {t: W_of(halves[t], fs, nps_h, bnd) for t in ("A", "B")}   # canonical
        W_match = {t: W_of(halves[t], fs, nps, bnd) for t in ("A", "B")}     # matched
        del Xr, halves
        for ph in ("task_test", "rest_post"):
            X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), dtype=np.float32)
            if X.shape[0] > X.shape[1]:
                X = X.T
            w = W_of(X, fs, nps, bnd)
            W_mixed[ph] = w; W_match[ph] = w
            del X
    except Exception as e:
        return None, f"{pat}/{band}: FAILED ({type(e).__name__}: {e})"

    r_mix = rho_sym_over_scales({p: eig(W_mixed[p]) for p in PHASES}, SGRID)
    r_mat = rho_sym_over_scales({p: eig(W_match[p]) for p in PHASES}, SGRID)
    rows = [dict(patient=pat, band=band, s=float(s),
                 rho_canonical_mixed_nperseg=float(r_mix[j]),
                 rho_matched_nperseg=float(r_mat[j]),
                 delta=float(r_mix[j] - r_mat[j]))
            for j, s in enumerate(SGRID)]
    return rows, (f"{pat}/{band} canonical[s=1]={r_mix[0]:+.3f} "
                  f"matched[s=1]={r_mat[0]:+.3f} d={r_mix[0]-r_mat[0]:+.3f} "
                  f"({time.time()-t0:.0f}s)")


def main():
    bands = os.environ.get("NL_BANDS", "alpha,beta").split(",")
    jobs = [(p, b) for b in bands for p in COHORT]
    workers = int(os.environ.get("NL_WORKERS", 4))
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[nperseg] {len(jobs)} cells, bands={bands}, workers={workers}", flush=True)
    t0, rows = time.time(), []
    with Pool(workers) as pool:
        for i, (rl, msg) in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                rows.extend(rl)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {msg} | {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        print("[nperseg] no rows"); return
    df.to_csv(OUT / "per_patient_scale.csv", index=False)

    print(f"\n{'='*88}\nA/B nperseg asymmetry — canonical (mixed) vs matched, per band\n{'='*88}")
    print(f"  {'band':8s} {'scale':>7s} {'canon_med':>10s} {'match_med':>10s} "
          f"{'d_med':>8s} {'p(d>0)':>8s}")
    for band in bands:
        for s in np.sort(df.s.unique()):
            x = df[(df.band == band) & np.isclose(df.s, s)].dropna()
            if len(x) < 5:
                continue
            try:
                p = float(wilcoxon(x.delta.values, alternative="greater")[1])
            except Exception:
                p = np.nan
            print(f"  {band:8s} {s:7.1f} {x.rho_canonical_mixed_nperseg.median():+10.3f} "
                  f"{x.rho_matched_nperseg.median():+10.3f} {x.delta.median():+8.3f} {p:8.4f}")
    print(f"\nwrote {OUT/'per_patient_scale.csv'}  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
