#!/usr/bin/env python3
"""W0-B: the lag-destroying null (N1) applied to the FIVE-PHASE cross-phase
functionals -- T_test, T_learn, T_infspec, T_infspec_pe.

Same rung as script 03's N1 (per-channel circular shift realized as its exact
frequency-domain phase ramp: ``|C_ij(f)|`` and all per-channel spectra preserved
exactly, only the lag randomized), now over the five-phase arc
{A, B, task_learn, task_test, rest_post} and scoring all four functionals per
scale. The observed values are recomputed from timeseries through the identical
code path, never read from a prior table.

READ THE CALIBRATION FIRST. ``T_infspec_pe`` is a PARTIAL correlation, and script
07 measures what it returns on a no-signal input. If it is not zero-centred there,
the p-values this script produces for it are not interpretable as "evidence of an
inference-specific trace", and the report must say so rather than quoting them.
This script therefore emits, for every functional, BOTH the surrogate-referenced
p and the raw observed value, so a reader can apply the calibration offset.

================== 5-point critical preamble (pre-registered) ==================
(1) CLAIM. The inference-specific persistence ``T_infspec_pe`` (and the encoding
    echo ``T_learn``) reflect task-induced reorganization of TIME-LAGGED
    interaction, not of coherence magnitude.
(2) NULL. N1: per-channel circular shift as an exact frequency-domain phase ramp,
    drawn independently per phase, preserving every ``|C_ij(f)|`` and every
    per-channel spectrum, destroying only lag.
(3) STRONGEST ALTERNATIVE. That these functionals are driven by the coherence
    MAGNITUDE geometry -- largely volume conduction and electrode placement --
    which is shared across all five phases and can generate cross-phase
    correlations with no lagged-interaction change at all.
(4) DOES IT CONTROL FOR IT, BY MECHANISM -- AND WHAT IT CANNOT REJECT. Yes for the
    magnitude-vs-lag question: the surrogate keeps each phase's own ``|C|`` and
    randomizes the phase, so a functional that survives requires lag structure.
    Licensed downstream because the readout is exactly invariant to the global
    rescale N1 induces (verified residual 0.0e+00).
    CANNOT REJECT: (i) an account in which the effect lives in cross-phase changes
    of coherence MAGNITUDE -- N1 holds those fixed by construction; (ii) any
    session-level confound (drift, nonstationarity, boundary placement), since all
    five phases keep their true boundaries -- that is N3's and the calibration's
    job; (iii) for ``T_infspec_pe`` specifically, N1 cannot repair a functional
    that is not zero-centred to begin with. A surrogate-referenced p-value tests
    "observed above surrogate", which is well-posed even for a biased estimator
    ONLY IF the bias is common to observed and surrogate; the calibration
    measures whether that holds.
(5) FALSIFICATION. If ``T_infspec_pe`` and ``T_learn`` fail the cohort gate under
    N1 at every scale, the encoding/inference decomposition is not a lag-carried
    phenomenon and cannot be presented as one.

Outputs: data/paper_final/w0b_nulls/func_n1/{per_patient_scale.csv, cohort_gate.csv}
"""
from __future__ import annotations
import argparse, json, os, sys, time
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
from lrg_eegfc.utils.fc.heat_multiscale import (
    CROSS_PHASE_FUNCTIONALS, cross_phase_functionals_over_scales, laplacian_eig,
)
from lrg_eegfc.utils.surrogate.timeseries_nulls import (
    coherency_from_csd, csd_from_segment_subset, imcoh_abs_from_coherency,
    lag_randomized_coherency, phase_randomized_coherency, segment_fft,
)
from audit_150_rho_sym_gate import COHORT, BANDS

PHASES5 = ("A", "B", "task_learn", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)
FRAC, BACKBONE = 0.20, "mst020"
R = 200
BASE_SEED = 20260826
OUT_ROOT = ROOT / "data" / "paper_final" / "w0b_nulls"


def eig(W):
    return laplacian_eig(select_backbone(W, BACKBONE, frac=FRAC))


def load_segs5(pat, band):
    """Per-segment band FFTs for the five canonical phases (production nperseg)."""
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs); nps_h = max(256, nps // 2)
    bnd = BRAIN_BANDS[band]
    out = {}
    Xr = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    if Xr.shape[0] > Xr.shape[1]:
        Xr = Xr.T
    T = Xr.shape[1]
    for tag, sl in (("A", slice(0, T // 2)), ("B", slice(T // 2, T))):
        out[tag] = segment_fft(Xr[:, sl], fs, nps_h, band=bnd) + (fs, nps_h)
    del Xr
    for ph in ("task_learn", "task_test", "rest_post"):
        X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
        if X.shape[0] > X.shape[1]:
            X = X.T
        out[ph] = segment_fft(X, fs, nps, band=bnd) + (fs, nps)
        del X
    return out


def per_cell(job):
    idx, pat, band, rung = job
    t0 = time.time()
    try:
        segs = load_segs5(pat, band)
    except Exception as e:
        return None, f"{pat}/{band}/{rung}: load failed ({type(e).__name__}: {e})"

    Cobs, Wobs = {}, {}
    for ph, (freqs, F, scale, fs, nps) in segs.items():
        Cobs[ph] = coherency_from_csd(csd_from_segment_subset(F, None, scale))
        Wobs[ph] = imcoh_abs_from_coherency(Cobs[ph])
    obs = cross_phase_functionals_over_scales({p: eig(Wobs[p]) for p in PHASES5}, SGRID)

    rng = np.random.default_rng([BASE_SEED, idx])
    surr = {k: np.full((R, SGRID.size), np.nan) for k in CROSS_PHASE_FUNCTIONALS}
    for r in range(R):
        Ws = {}
        for ph, (freqs, F, scale, fs, nps) in segs.items():
            if rung == "n1":
                Cs = lag_randomized_coherency(Cobs[ph], freqs, fs, rng, nperseg=nps)
            elif rung == "n2":
                Cs = phase_randomized_coherency(Cobs[ph], rng)
            else:
                raise ValueError(rung)
            Ws[ph] = imcoh_abs_from_coherency(Cs)
        so = cross_phase_functionals_over_scales({p: eig(Ws[p]) for p in PHASES5}, SGRID)
        for k in CROSS_PHASE_FUNCTIONALS:
            surr[k][r] = so[k]

    rows = []
    for k in CROSS_PHASE_FUNCTIONALS:
        for j, s in enumerate(SGRID):
            col = surr[k][:, j][np.isfinite(surr[k][:, j])]
            o = obs[k][j]
            rows.append(dict(
                patient=pat, band=band, rung=rung, func=k, s=float(s),
                obs=float(o) if np.isfinite(o) else np.nan,
                surr_p50=float(np.percentile(col, 50)) if col.size else np.nan,
                surr_p95=float(np.percentile(col, 95)) if col.size else np.nan,
                p=float(np.mean(col >= o)) if col.size and np.isfinite(o) else np.nan,
                n_surr=int(col.size)))
    return rows, (f"{pat}/{band}/{rung} T_test={obs['T_test'][0]:+.3f} "
                  f"T_learn={obs['T_learn'][0]:+.3f} "
                  f"T_isp={obs['T_infspec'][0]:+.3f} "
                  f"T_ispe={obs['T_infspec_pe'][0]:+.3f} ({time.time()-t0:.0f}s)")


def cohort_gate(df):
    out = []
    for (rung, fn, band), g in df.groupby(["rung", "func", "band"]):
        for s in np.sort(g.s.unique()):
            x = g[np.isclose(g.s, s)].dropna(subset=["obs", "surr_p50"])
            if len(x) < 5:
                continue
            d = x.obs.values - x.surr_p50.values
            try:
                p = float(wilcoxon(d, alternative="greater")[1]) if np.any(d != 0) else np.nan
            except Exception:
                p = np.nan
            out.append(dict(rung=rung, func=fn, band=band, s=float(s), gate_p=p,
                            n_pat=len(x), n_above=int((x.p < 0.05).sum()),
                            obs_med=float(np.median(x.obs)),
                            surr_med=float(np.median(x.surr_p50)),
                            margin_med=float(np.median(d))))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rungs", default="n1")
    ap.add_argument("--bands", default="delta,theta,alpha,beta,low_gamma")
    ap.add_argument("--patients", default="")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("NL_WORKERS", 6)))
    a = ap.parse_args()
    rungs = a.rungs.split(","); bands = a.bands.split(",")
    pats = a.patients.split(",") if a.patients else COHORT
    jobs = [(i, p, b, rg) for i, (rg, b, p) in
            enumerate((rg, b, p) for rg in rungs for b in bands for p in pats)]
    print(f"[func-null] {len(jobs)} cells | rungs={rungs} bands={bands} R={R} "
          f"funcs={CROSS_PHASE_FUNCTIONALS} workers={a.workers}", flush=True)

    t0, rows = time.time(), []
    with Pool(a.workers) as pool:
        for i, (rl, msg) in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                rows.extend(rl)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {msg} | {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)

    df = pd.DataFrame(rows)
    if df.empty:
        print("[func-null] no rows"); return
    for rung in rungs:
        d = df[df.rung == rung]
        if d.empty:
            continue
        od = OUT_ROOT / f"func_{rung}"
        od.mkdir(parents=True, exist_ok=True)
        pp = od / "per_patient_scale.csv"
        if pp.exists():
            d = (pd.concat([pd.read_csv(pp), d], ignore_index=True)
                   .drop_duplicates(subset=["patient", "band", "rung", "func", "s"],
                                    keep="last"))
        d.to_csv(pp, index=False)
        g = cohort_gate(d)
        g.to_csv(od / "cohort_gate.csv", index=False)
        (od / "config.json").write_text(json.dumps(dict(
            rung=rung, phases=list(PHASES5), functionals=list(CROSS_PHASE_FUNCTIONALS),
            backbone=BACKBONE, frac=FRAC, R=R, s_grid=[float(x) for x in SGRID],
            cohort=pats, bands=bands, base_seed=BASE_SEED,
            caveat="T_infspec_pe is a partial correlation; read calibration/ first"),
            indent=2))
        if g.empty or "func" not in g.columns:
            print("  (gate needs >=5 patients)"); continue
        print(f"\n=== {rung}: cohort gate by functional x band (scales cleared / 16) ===",
              flush=True)
        print(f"  {'func':14s} {'band':11s} {'clears':>8s} {'min p':>8s} "
              f"{'obs_med':>9s} {'surr_med':>9s}", flush=True)
        for fn in CROSS_PHASE_FUNCTIONALS:
            for band in bands:
                gb = g[(g.func == fn) & (g.band == band)].dropna(subset=["gate_p"])
                if gb.empty:
                    continue
                b = gb.loc[gb.gate_p.idxmin()]
                print(f"  {fn:14s} {band:11s} {int((gb.gate_p<0.05).sum()):>4d}/{len(gb):<3d} "
                      f"{b.gate_p:8.4f} {b.obs_med:+9.3f} {b.surr_med:+9.3f}", flush=True)
    print(f"\n[func-null] {len(df)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
