#!/usr/bin/env python3
"""audit_168 — DURATION-MATCHED drift null for the inference-specific decomposition.

Refines audit_167 (which used 5 EQUAL windows) to adjudicate the one open drift
gap: the inference-specific functional T_infspec_pe = partial_rho(f, p | e) with
f = D_TT - D_TL (task_test minus task_learn). audit_167's equal-window drift arc
gave task_learn / task_test / the baseline halves the SAME length, so it did NOT
reproduce the real arc's per-phase SNR asymmetry -- exactly the asymmetry the
inference-specific difference f is built from. This audit windows a single pre-task
resting recording (rest_pre, pure drift, no task -> no consolidation possible) into
5 temporally-ordered, non-overlapping slots whose lengths are PROPORTIONAL to the
real arc-slot durations [L_pre/2, L_pre/2, L_TL, L_TT, L_RP]. So each sham slot's
FC noise level matches the real phase it stands in for (duration-RATIO-matched),
and preA/preB get the two largest slots (~L_pre/2 each) as in the real half-split.

Construction caveat (honest, conservative): the 5 real durations sum to the WHOLE
session (> L_pre), so all slots are shrunk by s = L_pre / sum(d) < 1 to tile a
pure pre-task recording. Absolute lengths are ~s x shorter than the real phases ->
FC uniformly noisier -> drift T biased LOW. This is conservative for a positive
drift finding (if drift already reproduces the real trace at reduced SNR, the true
drift is >= that). Absolute duration-matching is impossible with clean (pre-task)
data: no single task-free recording is as long as the full session.

Money test: paired, per (patient, band), real_ispe - drift_ispe, cohort Wilcoxon
one-sided greater. If real inference-specific trace does NOT significantly exceed
the duration-matched drift arc, the N2 inference-specific claim is drift-unverified.

5-point critical preamble
1. Claim: the inference-specific beta refinement (T_infspec_pe, D_TT vs D_TL
   correlating with the persistent rest change) is task consolidation, not drift.
2. Null (duration-matched drift): a 5-slot arc built from pure pre-task rest, slots
   sized in proportion to the real phase durations, reproduces T_infspec_pe.
3. Strongest alt: monotonic within-session nonstationarity makes any later-vs-earlier
   contrast correlate; the equal-window audit_167 could not match the TL/TT/baseline
   SNR asymmetry that f = D_TT - D_TL depends on -> this fixes that.
4. Reach: rest_pre precedes any task -> positive T here is drift by construction.
   Duration-RATIO matched (not absolute); absolute shrink is conservative (noisier
   -> drift underestimated). Does NOT prove residual real trace is drift-free; a
   non-significant real-minus-drift = "not shown to beat drift", not "is drift".
5. Falsify (of the task claim): if drift_ispe ~ real_ispe (paired real-drift not > 0),
   the inference-specific trace is not shown to exceed drift -> drift-unverified.
   If real_ispe >> drift_ispe (paired p < 0.05), the N2 claim survives duration-matched
   drift and graduates from "drift-unverified" to drift-controlled.
"""
from __future__ import annotations
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
from lrg_eegfc.utils.surrogate.matched_strength import (
    cophenetic_condensed_from_adjacency as coph,
)

# Reuse audit_167's arc helpers verbatim (no fork) + audit_63's phase FC loader.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_167_drift_null_ladder import _band_abs, _iu, _rho, _arc  # type: ignore
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = list(BRAIN_BANDS)
PH = ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")
# real durations that set the slot proportions (rest_pre is split in two halves)
DUR_PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")


def phase_lengths(pat):
    """Sample count (T) of each real phase; sets duration-ratio of the sham slots."""
    out = {}
    for ph in DUR_PHASES:
        X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
        out[ph] = X.shape[1]
    return out


def proportional_bounds(L_pre, dur):
    """5 ordered, non-overlapping slot boundaries in [0, L_pre], lengths proportional
    to the real arc-slot durations [L_pre/2, L_pre/2, L_TL, L_TT, L_RP]."""
    d = np.array([dur["rest_pre"] / 2.0, dur["rest_pre"] / 2.0,
                  dur["task_learn"], dur["task_test"], dur["rest_post"]], float)
    cum = np.concatenate([[0.0], np.cumsum(d / d.sum())])
    b = np.round(cum * L_pre).astype(int)
    b[-1] = L_pre
    return [(int(b[i]), int(b[i + 1])) for i in range(5)], d


def windowed_specs(X, fs, bounds):
    """One Welch cross-spectrum per window (band-agnostic); slice bands later.
    Returns list of (freqs, Coh) for the 5 slots. ~6x fewer Welch than audit_167."""
    nper = nperseg_for_fs(fs)
    specs = []
    for a, b in bounds:
        Xi = X[:, a:b]
        fr, C = compute_msc_welch(Xi, fs, nperseg=min(nper, Xi.shape[1]), metric="imcoh")
        specs.append((fr, C))
    return specs


def per_patient(pat):
    """All 6 bands for one patient. Drift arc from proportional rest_pre windows +
    real arc from cached phase FCs. Returns list of row dicts."""
    X = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    fs = FS_OVERRIDES.get(pat, 2048.0)
    dur = phase_lengths(pat)
    bounds, d = proportional_bounds(X.shape[1], dur)
    shrink = X.shape[1] / d.sum()                 # absolute-length shrink factor s
    specs = windowed_specs(X, fs, bounds)         # 5 Welch, reused across bands
    win_len = [b - a for a, b in bounds]

    rows = []
    for band in BANDS:
        flo, fhi = BRAIN_BANDS[band]
        Ws = [_band_abs(C, fr, flo, fhi) for (fr, C) in specs]     # drift-slot FCs
        Draw = {ph: _iu(Ws[i]) for i, ph in enumerate(PH)}
        Dcop = {ph: coph(Ws[i]) for i, ph in enumerate(PH)}
        tt_r, ip_r = _arc(Draw)                                    # drift
        tt_c, ip_c = _arc(Dcop)
        # real arc from cached phase FCs (task genuinely happened)
        Rraw = {ph: _iu(load_phase_fc(pat, ph, band)) for ph in PH}
        Rcop = {ph: coph(load_phase_fc(pat, ph, band)) for ph in PH}
        TT_r, IP_r = _arc(Rraw)
        TT_c, IP_c = _arc(Rcop)
        rows.append(dict(
            patient=pat, band=band, shrink=shrink,
            win_min=int(min(win_len)), win_max=int(max(win_len)),
            raw_real_Ttest=TT_r, raw_drift_Ttest=tt_r,
            coph_real_Ttest=TT_c, coph_drift_Ttest=tt_c,
            raw_real_ispe=IP_r, raw_drift_ispe=ip_r,
            coph_real_ispe=IP_c, coph_drift_ispe=ip_c,
        ))
    return rows


def _wilcox_gt(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return float("nan")
    try:
        return float(wilcoxon(x, alternative="greater")[1])
    except Exception:
        return float("nan")


def main():
    print("[audit_168] DURATION-MATCHED drift null (inference-specific), raw + cophenetic")
    print("  slots proportional to real [L_pre/2, L_pre/2, L_TL, L_TT, L_RP]; "
          "Welch once/window, sliced over 6 bands\n", flush=True)
    t0 = time.time()
    rows = []
    for pi, pat in enumerate(COHORT, 1):
        try:
            rows.extend(per_patient(pat))
        except Exception as exc:  # noqa: BLE001
            print(f"  FAIL {pat}: {type(exc).__name__}: {exc}", flush=True)
        el = time.time() - t0
        sh = np.median([r["shrink"] for r in rows if r["patient"] == pat] or [np.nan])
        print(f"[{pi}/10] {pat} done  {el:.0f}s  ETA {el / pi * (10 - pi):.0f}s  "
              f"(shrink s~{sh:.2f})", flush=True)

    df = pd.DataFrame(rows)
    OUT = ROOT / "data" / "audit" / "duration_matched_drift_infspec"
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)

    # ---- cohort gates ----
    print("\n=== DURATION-MATCHED DRIFT vs REAL "
          "(cohort median; paired real-drift Wilcoxon>0 is the money test) ===")
    hdr = (f"{'band':<11}"
           f"{'coph_realIP':>12}{'coph_drIP':>11}{'r-d p':>7}{'dr>0 p':>8}   "
           f"{'raw_realIP':>11}{'raw_drIP':>10}{'r-d p':>7}{'dr>0 p':>8}")
    print(hdr)
    summ = []
    for band in BANDS:
        s = df[df.band == band]
        cr, cd = s.coph_real_ispe, s.coph_drift_ispe
        rr, rd = s.raw_real_ispe, s.raw_drift_ispe
        row = dict(
            band=band,
            coph_real_ispe=float(cr.median()), coph_drift_ispe=float(cd.median()),
            coph_real_minus_drift_p=_wilcox_gt((cr - cd).values),
            coph_drift_gt0_p=_wilcox_gt(cd.values),
            raw_real_ispe=float(rr.median()), raw_drift_ispe=float(rd.median()),
            raw_real_minus_drift_p=_wilcox_gt((rr - rd).values),
            raw_drift_gt0_p=_wilcox_gt(rd.values),
        )
        summ.append(row)
        print(f"{band:<11}"
              f"{row['coph_real_ispe']:>+12.3f}{row['coph_drift_ispe']:>+11.3f}"
              f"{row['coph_real_minus_drift_p']:>7.3f}{row['coph_drift_gt0_p']:>8.3f}   "
              f"{row['raw_real_ispe']:>+11.3f}{row['raw_drift_ispe']:>+10.3f}"
              f"{row['raw_real_minus_drift_p']:>7.3f}{row['raw_drift_gt0_p']:>8.3f}")
    pd.DataFrame(summ).to_csv(OUT / "cohort_gates.csv", index=False)

    print("\n=== whole-task T_test (continuity vs audit_167 / locked C2) ===")
    print(f"{'band':<11}{'coph_realTt':>12}{'coph_drTt':>11}{'r-d p':>7}   "
          f"{'raw_realTt':>11}{'raw_drTt':>10}{'r-d p':>7}")
    for band in BANDS:
        s = df[df.band == band]
        cr, cd = s.coph_real_Ttest, s.coph_drift_Ttest
        rr, rd = s.raw_real_Ttest, s.raw_drift_Ttest
        print(f"{band:<11}{cr.median():>+12.3f}{cd.median():>+11.3f}"
              f"{_wilcox_gt((cr - cd).values):>7.3f}   "
              f"{rr.median():>+11.3f}{rd.median():>+10.3f}"
              f"{_wilcox_gt((rr - rd).values):>7.3f}")

    sh = df.shrink
    print(f"\nabsolute-shrink s: median {sh.median():.2f} (min {sh.min():.2f}) "
          f"-- slots are ~s x real phase length (conservative: noisier -> drift low)")
    print("READ: coph_real_minus_drift_p (beta, inference-specific) < 0.05 => "
          "N2 inference-specific trace BEATS duration-matched drift (drift-controlled).")
    print("      >= 0.05 => not shown to beat drift => stays drift-unverified.")
    print(f"[audit_168] {len(df)} rows in {time.time() - t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
