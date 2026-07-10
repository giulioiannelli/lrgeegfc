#!/usr/bin/env python3
"""audit_167 — within-recording DRIFT null for the pairwise-descriptor ladder.

Tests the user's hypothesis (2026-07-09): the raw-edge inference/whole-task trace
is session DRIFT (raw FC "always shows trace"), which matched-strength does NOT
control for. Builds a pure-drift arc by windowing a SINGLE pre-task resting
recording (rest_pre) into 5 temporally ordered segments and running the IDENTICAL
rho_sym arc through it — no task, no consolidation possible, so any positive T is
drift / within-recording nonstationarity. Run through BOTH raw edges and
cophenetic.

Reading:
  - if raw_drift ~ raw_real  -> the raw-edge trace is drift.
  - if coph_drift << coph_real -> cophenetic FILTERS drift (multiscale necessary).
  - phase map: preA=w1, preB=w2, task_learn=w3, task_test=w4, rest_post=w5 (all
    inside rest_pre). SNR caveat: 5 sub-windows are ~1/5 duration -> noisier FC ->
    T_drift biased LOW (conservative for a positive drift finding).

5-point preamble
1. Claim: the ladder traces (raw + cophenetic) are task consolidation, not drift.
2. Null (drift): a same-length arc built from pure pre-task rest reproduces T.
3. Strongest alt this addresses: monotonic session nonstationarity makes any
   later-vs-earlier contrast correlate (raw FC always shows trace).
4. Reach: rest_pre is BEFORE any task -> a positive T here CANNOT be consolidation;
   it is drift by construction. Does NOT match per-phase duration exactly (windows
   shorter) -> SNR caveat noted, conservative direction.
5. Falsify (of the task claim): if raw_drift reproduces raw_real, the raw trace is
   drift; if coph_drift also reproduces coph_real, even cophenetic is drift-driven.
"""
from __future__ import annotations
import sys, time
import numpy as np
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import BRAIN_BANDS, nperseg_for_fs, FS_OVERRIDES
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.surrogate.matched_strength import cophenetic_condensed_from_adjacency as coph
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore

COHORT = ["Pat_02","Pat_03","Pat_05","Pat_06","Pat_07","Pat_08","Pat_10","Pat_13","Pat_14","Pat_15"]
BANDS = list(BRAIN_BANDS)
NWIN = 5
PH = ("rest_pre_A","rest_pre_B","task_learn","task_test","rest_post")


def _band_abs(Coh, freqs, flo, fhi):
    m = (freqs >= flo) & (freqs <= fhi)
    A = np.abs(Coh[:, :, m]).mean(axis=-1)
    A = 0.5 * (A + A.T); np.fill_diagonal(A, 0.0)
    return np.clip(A, 0.0, 1.0)


def _iu(W):
    return np.asarray(W, float)[np.triu_indices(W.shape[0], 1)]


def _rho(a, b):
    return float(spearmanr(a, b)[0])


def _prho(a, b, c):
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    d = np.sqrt(max(0.0, (1 - rac**2) * (1 - rbc**2)))
    return (rab - rac * rbc) / d if d > 0 else np.nan


def _arc(D):  # sym over arms: returns (T_test, T_infspec_pe)
    def arm(bt, br):
        e = D["task_learn"] - D[bt]; f = D["task_test"] - D["task_learn"]
        g = D["task_test"] - D[bt]; p = D["rest_post"] - D[br]
        return _rho(g, p), _prho(f, p, e)
    a1 = arm("rest_pre_A", "rest_pre_B"); a2 = arm("rest_pre_B", "rest_pre_A")
    return 0.5*(a1[0]+a2[0]), 0.5*(a1[1]+a2[1])


def windowed_fc(pat, band):
    """5 imcoh_abs FCs from ordered windows of rest_pre (pure pre-task drift)."""
    X = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nper = nperseg_for_fs(fs)
    T = X.shape[1]
    flo, fhi = BRAIN_BANDS[band]
    # validation: full-recording FC via this pipeline vs cached loader
    freqs, Coh = compute_msc_welch(X, fs, nperseg=nper, metric="imcoh")
    full = _band_abs(Coh, freqs, flo, fhi)
    cached = load_phase_fc(pat, "rest_pre", band)
    val = _rho(_iu(full), _iu(cached)) if cached.shape == full.shape else np.nan
    Ws = []
    for i in range(NWIN):
        Xi = X[:, i * T // NWIN:(i + 1) * T // NWIN]
        fr, C = compute_msc_welch(Xi, fs, nperseg=min(nper, Xi.shape[1]), metric="imcoh")
        Ws.append(_band_abs(C, fr, flo, fhi))
    return Ws, val


def main():
    print(f"[audit_167] within-rest_pre {NWIN}-window DRIFT null, raw + cophenetic\n")
    rows = []
    t0 = time.time()
    for pi, pat in enumerate(COHORT, 1):
        for band in BANDS:
            try:
                Ws, val = windowed_fc(pat, band)
                # drift arc: preA=w0 preB=w1 TL=w2 TT=w3 RP=w4
                Draw = {ph: _iu(Ws[i]) for i, ph in enumerate(PH)}
                Dcop = {ph: coph(Ws[i]) for i, ph in enumerate(PH)}
                tt_r, ip_r = _arc(Draw)
                tt_c, ip_c = _arc(Dcop)
                # real (cached phases)
                Rraw = {ph: _iu(load_phase_fc(pat, ph, band)) for ph in PH}
                Rcop = {ph: coph(load_phase_fc(pat, ph, band)) for ph in PH}
                TT_r, IP_r = _arc(Rraw); TT_c, IP_c = _arc(Rcop)
                rows.append(dict(patient=pat, band=band, val=val,
                                 raw_drift_Ttest=tt_r, raw_real_Ttest=TT_r,
                                 raw_drift_ispe=ip_r, raw_real_ispe=IP_r,
                                 coph_drift_Ttest=tt_c, coph_real_Ttest=TT_c,
                                 coph_drift_ispe=ip_c, coph_real_ispe=IP_c))
            except Exception as exc:  # noqa
                print(f"  FAIL {pat}/{band}: {type(exc).__name__}: {exc}", flush=True)
        el = time.time() - t0
        print(f"[{pi}/10] {pat} done  {el:.0f}s ETA {el/pi*(10-pi):.0f}s "
              f"(val Spearman full-vs-cached ~{np.nanmedian([r['val'] for r in rows if r['patient']==pat]):.3f})",
              flush=True)

    import pandas as pd
    df = pd.DataFrame(rows)
    OUT = ROOT / "data" / "audit" / "drift_null_ladder"
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)

    def wp(x):
        x = x[np.isfinite(x)]
        try: return wilcoxon(x, alternative="greater")[1]
        except Exception: return float("nan")
    print("\n=== DRIFT vs REAL (cohort median; Wilcoxon p that DRIFT>0) ===")
    print(f"{'band':<11}{'raw_realTt':>11}{'raw_drTt':>10}{'p':>7}   "
          f"{'coph_realTt':>12}{'coph_drTt':>11}{'p':>7}   "
          f"{'raw_dr_ispe':>12}{'coph_dr_ispe':>13}")
    for band in BANDS:
        s = df[df.band == band]
        print(f"{band:<11}{s.raw_real_Ttest.median():>+11.3f}{s.raw_drift_Ttest.median():>+10.3f}"
              f"{wp(s.raw_drift_Ttest.values):>7.3f}   "
              f"{s.coph_real_Ttest.median():>+12.3f}{s.coph_drift_Ttest.median():>+11.3f}"
              f"{wp(s.coph_drift_Ttest.values):>7.3f}   "
              f"{s.raw_drift_ispe.median():>+12.3f}{s.coph_drift_ispe.median():>+13.3f}")
    print(f"\nfull-vs-cached validation Spearman: median {df.val.median():.3f} "
          f"(min {df.val.min():.3f}) -- >0.99 = channels/pipeline aligned")
    print("\nREAD: raw_drTt ~ raw_realTt => raw whole-task trace is DRIFT.")
    print("      coph_drTt << coph_realTt => cophenetic filters drift.")
    print("      *_dr_ispe = inference-specific under pure drift (compare to real raw .124 / coph .091).")
    print(f"[audit_167] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
