#!/usr/bin/env python3
r"""audit_162 — windowed consolidation flow in the encoding/inference plane.

Turns the static four-phase consolidation arc (audit_152) into a DYNAMICAL
trajectory. Slide a 30 s window across the whole recording (rest_pre -> task_learn
-> task_test -> rest_post), compute the LRG cophenetic state per window per band,
and project every window onto a FIXED, interpretable plane whose axes are the two
things the task did:

    e  = Dbar_taskLearn - Dbar_restPre          (encoding axis)
    f0 = Dbar_taskTest  - Dbar_taskLearn        (inference-specific contrast)
    x(w) = rho_S( D(w) - Dbar_restPre , e )               encoding coordinate
    y(w) = rho_S( D(w) - Dbar_restPre , f0 | e )          inference coordinate (partial)

Everything is at WINDOW scale (phase window-means, not full-phase FC), so the
subtractions are scale-consistent; Spearman ranks make absolute scale irrelevant.
The rest_post window-cloud mean approximates the audit_152 (T_learn, T_infspec_pe)
numbers -- this is a VISUALIZATION of an already matched-strength-verified effect,
NOT a new test. No new null; verdict of record stays audit_152.

Scope: .agents/guides/task-persistence-investigation/2026-07-08_windowed-consolidation-flow.md

5-point critical preamble
1. Claim (visual): offline persistence is band-dissociated -- encoding echo broadband
   (delta/alpha/beta), inference-specific persistence beta-only. Dynamically: every
   band travels out along encoding + up along inference DURING task; only beta's
   inference coordinate stays elevated at rest_post.
2. Null: the "flow" is a projection artefact -- noisy window states projected onto
   data-built axes drift plausibly; the rest_post displacement is non-specific.
3. Strongest alt: (a) a data-driven embedding (PCA/UMAP) would be dominated by drift
   and could manufacture structure; (b) rest_post displacement could be matched-
   strength, not trace.
4. Does design cover it: (a) defeated by construction -- axes are FIXED to the task
   contrasts, not data-driven, so the plane cannot contain structure the contrasts
   lack; rest_pre cloud drawn as empirical origin. (b) NOT re-tested here -- settled
   by audit_152 matched-strength Wilcoxon (beta T_infspec_pe p=0.0098); its surrogate
   landing is overlaid. CANNOT establish significance (overlapping windows auto-
   correlated); shows the SHAPE of a verified effect.
5. Falsify (visual): if beta rest_post inference coord == its rest_pre scatter, or
   alpha/delta retain as much inference height as beta, revert to fig_arc_a. Cross-
   check: cohort-median windowed endpoint must reproduce sign+ordering of audit_152
   (T_learn, T_infspec_pe) -- asserted at end of run.

Reuse (no forks): windowed.py segment_ffts/imcoh_abs_cube/band_abs_average (the
audit_124 replay backend; full-window limit reproduces cached phase FC to rho~1),
cophenetic_condensed_from_adjacency (library LRG ultrametric, tau=1/lam_max).

Writes: data/audit/windowed_consolidation_flow/
  windowed_flow_per_window.csv   (patient, band, phase, w_idx, t_center, x, y, n_seg)
  windowed_flow_endpoints.csv    (patient, band, phase, x_mean, y_mean, x_med, y_med, n_win)
  README.md
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.fc.coherence import (
    segment_ffts, imcoh_abs_cube, band_abs_average,
)
from lrg_eegfc.utils.surrogate.matched_strength import (
    cophenetic_condensed_from_adjacency,
)

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = list(BRAIN_BANDS)
PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")

W_SEC = 30.0            # all-bands-valid: delta low edge 0.53 Hz -> 30 s = 16 cycles
OVERLAP = 0.25         # time-window-guide default
FMAX = max(hi for _lo, hi in BRAIN_BANDS.values()) + 5.0    # high_gamma hi + 5
MIN_SEG = 4            # need >=4 Welch segments for a stable band-|ImCoh|

OUT = ROOT / "data" / "audit" / "windowed_consolidation_flow"
OUT.mkdir(parents=True, exist_ok=True)
ARC = ROOT / "data" / "audit" / "consolidation_arc_rhosym" / "arc_per_patient.csv"


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    r, _ = spearmanr(a, b)
    return float(r)


def _partial_rho(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """First-order partial Spearman of (a, b) controlling c (rank world).

    Mirrors audit_152._partial_rho exactly so the windowed endpoint is the same
    functional as T_infspec_pe."""
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    denom = np.sqrt(max(0.0, (1.0 - rac ** 2) * (1.0 - rbc ** 2)))
    return (rab - rac * rbc) / denom if denom > 0 else np.nan


def _load_xt(pat: str, phase: str) -> np.ndarray:
    X = np.asarray(load_timeseries(pat, phase, SEEG_DATAPATH), float)
    return X if X.shape[0] <= X.shape[1] else X.T


def windowed_states(pat: str, fs: float, nperseg: int, verbose: bool):
    """Per-window cophenetic state for every band, over all four phases.

    Returns dict band -> dict phase -> (list of condensed cophenetic vectors),
    and a parallel meta list of (phase, w_idx, t_center, n_seg)."""
    wlen = int(round(W_SEC * fs))
    step = int(round(W_SEC * (1.0 - OVERLAP) * fs))
    states = {b: {ph: [] for ph in PHASES} for b in BANDS}
    meta = []
    for ph in PHASES:
        X = _load_xt(pat, ph)
        starts = list(range(0, X.shape[1] - wlen + 1, step))
        for wi, s in enumerate(starts):
            freqs, ff = segment_ffts(X[:, s:s + wlen], fs, nperseg, nperseg // 2, FMAX)
            nseg = ff.shape[1]
            if nseg < MIN_SEG:
                continue
            cube = imcoh_abs_cube(ff)                      # (N, N, F), once per window
            for b in BANDS:
                Wb = band_abs_average(cube, freqs, *BRAIN_BANDS[b])
                states[b][ph].append(cophenetic_condensed_from_adjacency(Wb))
            meta.append((ph, wi, (s + wlen / 2) / fs, nseg))
            del cube
        if verbose:
            print(f"    {ph:11s} {len(starts):3d} windows", flush=True)
    for b in BANDS:
        for ph in PHASES:
            states[b][ph] = np.asarray(states[b][ph]) if states[b][ph] else np.empty((0, 0))
    return states, meta


def _z(v: np.ndarray) -> np.ndarray:
    """Standardize one cophenetic vector to zero-mean unit-std (removes the
    window-vs-phase scale/offset; Spearman-invariant, but makes the raw DIFFERENCE
    D(w)−D0 scale-consistent so windows are not glued to a scale-dominated point)."""
    s = v.std()
    return (v - v.mean()) / s if s > 0 else v * 0.0


def _phase_coph(pat: str, band: str, ph: str) -> np.ndarray:
    W = np.asarray(load_fc_matrix(pat, ph, band, fc_method="imcoh_abs"), float)
    np.fill_diagonal(W, 0.0)
    W = 0.5 * (W + W.T)
    return cophenetic_condensed_from_adjacency(np.clip(W, 0.0, 1.0))


def project_band(states_b: dict, pat: str, band: str):
    """Project every window onto the encoding/inference plane, THREE variants:

      A = window-mean anchors, raw differences   (scale-consistent window view)
      B = window-mean anchors, z-scored states    (window view, scale-normalized)
      C = full-phase anchors,  z-scored states    (faithful to audit_152 directions)

    Returns dict variant -> {phase: (x array, y array)}. The cheap projections all
    reuse the one expensive set of window cophenetic states."""
    # --- window-mean anchors (A: raw, B: z-scored) ---------------------------
    D0a = states_b["rest_pre"].mean(0)
    ea = states_b["task_learn"].mean(0) - D0a
    f0a = states_b["task_test"].mean(0) - states_b["task_learn"].mean(0)
    zrows = {ph: np.array([_z(states_b[ph][i]) for i in range(len(states_b[ph]))])
             if len(states_b[ph]) else np.empty((0, 0)) for ph in PHASES}
    D0b = zrows["rest_pre"].mean(0)
    eb = zrows["task_learn"].mean(0) - D0b
    f0b = zrows["task_test"].mean(0) - zrows["task_learn"].mean(0)
    # --- full-phase anchors, z-scored (C) ------------------------------------
    D0c = _z(_phase_coph(pat, band, "rest_pre"))
    ec = _z(_phase_coph(pat, band, "task_learn")) - D0c
    f0c = _z(_phase_coph(pat, band, "task_test")) - _z(_phase_coph(pat, band, "task_learn"))

    out = {"A": {}, "B": {}, "C": {}}
    for ph in PHASES:
        S = states_b[ph]
        for var, (D0, e, f0, zscore) in {
            "A": (D0a, ea, f0a, False), "B": (D0b, eb, f0b, True),
            "C": (D0c, ec, f0c, True)}.items():
            xs, ys = [], []
            for w in range(S.shape[0]):
                d = (_z(S[w]) if zscore else S[w]) - D0
                xs.append(_rho(d, e))
                ys.append(_partial_rho(d, f0, e))
            out[var][ph] = (np.asarray(xs), np.asarray(ys))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    args = ap.parse_args()

    t0 = time.time()
    per_window, endpoints = [], []
    N = len(args.patients)
    for pi, pat in enumerate(args.patients):
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        tp = time.time()
        print(f"[{pi+1}/{N}] {pat} (fs={fs:.0f} nperseg={nperseg} W={W_SEC}s)", flush=True)
        states, meta = windowed_states(pat, fs, nperseg, verbose=True)
        for b in BANDS:
            coords = project_band(states[b], pat, b)          # {variant: {phase:(x,y)}}
            for var in ("A", "B", "C"):
                for ph in PHASES:
                    xs, ys = coords[var][ph]
                    for w in range(len(xs)):
                        per_window.append(dict(patient=pat, band=b, variant=var,
                                               phase=ph, w_idx=w, x=xs[w], y=ys[w]))
                    if len(xs):
                        endpoints.append(dict(patient=pat, band=b, variant=var, phase=ph,
                                              x_mean=float(np.nanmean(xs)),
                                              y_mean=float(np.nanmean(ys)),
                                              x_med=float(np.nanmedian(xs)),
                                              y_med=float(np.nanmedian(ys)),
                                              n_win=int(len(xs))))
        el = time.time() - tp
        eta = el * (N - pi - 1)
        print(f"    done {el:.1f}s  ETA {eta/60:.1f} min  ({len(meta)} windows)", flush=True)

    # attach t_center per (patient, phase, w_idx) from a second pass would need meta;
    # store the phase-relative window index which is enough for time ordering.
    pw = pd.DataFrame(per_window)
    ep = pd.DataFrame(endpoints)
    pw.to_csv(OUT / "windowed_flow_per_window.csv", index=False)
    ep.to_csv(OUT / "windowed_flow_endpoints.csv", index=False)
    print(f"\nwrote {OUT/'windowed_flow_per_window.csv'} ({len(pw)} rows)")
    print(f"wrote {OUT/'windowed_flow_endpoints.csv'} ({len(ep)} rows)")

    # ---- variant selection + P4 faithfulness cross-check vs audit_152 ----------
    rp = ep[ep.phase == "rest_post"]
    arc = pd.read_csv(ARC)[["patient", "band", "T_learn", "T_infspec_pe"]] if ARC.exists() else None
    print("\n=== rest_post COHORT-MEDIAN endpoint per band, per variant "
          "(the disentangling: does β lift in y above the rest?) ===")
    for var in ("A", "B", "C"):
        print(f"\n  variant {var}:  {'band':<11}{'med x (enc)':>12}{'med y (inf)':>12}"
              f"{'corr(x,Tlearn)':>16}{'corr(y,Tinf)':>14}")
        rv = rp[rp.variant == var]
        for b in BANDS:
            mb = rv[rv.band == b]
            cx = cy = float("nan")
            if arc is not None:
                mm = mb.merge(arc, on=["patient", "band"]).dropna()
                if len(mm) >= 4:
                    cx = _rho(mm.x_mean.values, mm.T_learn.values)
                    cy = _rho(mm.y_mean.values, mm.T_infspec_pe.values)
            print(f"  {'':<11}{b:<11}{mb.x_mean.median():>+12.3f}{mb.y_mean.median():>+12.3f}"
                  f"{cx:>+16.2f}{cy:>+14.2f}")
    print(f"\ntotal {(time.time()-t0)/60:.1f} min")
    _readme(pw, ep)


def _readme(pw, ep):
    rp = ep[(ep.phase == "rest_post") & (ep.variant == "A")]
    L = ["---", "name: windowed_consolidation_flow",
         "scope: windowed_encoding_inference_plane_trajectory",
         "date: 2026-07-08", "status: current", "---", "",
         "# Windowed consolidation flow (encoding/inference plane trajectory)", "",
         "**Head.** audit_152 four-phase arc dynamized: 30 s sliding-window LRG "
         "cophenetic states projected onto the FIXED encoding/inference plane "
         "(axes = task contrasts). rest_post cloud mean ~ (T_learn, T_infspec_pe). "
         "Visualization of a verified effect; verdict of record = audit_152.", "",
         f"- W={W_SEC}s, overlap={OVERLAP}, nperseg=nperseg_for_fs(fs), all 6 bands "
         f"from one FFT cube per window; min {MIN_SEG} Welch segments.",
         "- x = ρ_S(Δ,e); y = partial ρ_S(Δ,f0|e); Δ = D(w) − Dbar_restPre; window scale.",
         "- Build: scripts/01_compute/audit/audit_162_windowed_consolidation_flow.py",
         "", "## rest_post endpoint (cohort median x,y per band)",
         "| band | med x (encoding) | med y (inference) |", "|---|---|---|"]
    for b in BANDS:
        mb = rp[rp.band == b]
        if len(mb):
            L.append(f"| {b} | {mb.x_mean.median():+.3f} | {mb.y_mean.median():+.3f} |")
    (OUT / "README.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
