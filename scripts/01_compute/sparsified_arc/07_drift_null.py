#!/usr/bin/env python3
"""D2-drift -- the DECISIVE band-selectivity test: does the tau-trace exceed a
temporal-drift floor, per band, per scale?

The matched-strength gate on the percolation backbone LOST band-selectivity
(delta/theta/beta/gamma_high all clear) -- and "all bands trace" is a NEGATIVE
result. The drift null is the stronger control that should restore it: window a
single rest_pre recording into 5 ordered non-overlapping segments and map
A=w0, B=w1, task_test=w3, rest_post=w4 (audit_167 construction). Since every
window precedes any task, ANY positive rho_sym(s) on this arc is PURE DRIFT.
A genuine task trace must exceed this floor.

VALID only for the whole-task standard trace T_test = rho_sym (a simple
later-vs-earlier correlation). NOT applied to the conditional T_infspec_pe
(that estimator entangles drift removal -- use matched-strength; locked caveat
S2_drift_controls.md).

Reuses the cached D2 observed obs(s) (trace_arc/{band}/{pat}.npz) so the real
arm is bit-identical to D2; computes the drift arm fresh on the SAME per-cell
scale grid via the same kernel (percolation_backbone -> rho_sym_over_scales).
Cohort gate per scale = Wilcoxon(obs - drift, greater); band passes drift if it
exceeds drift where it also clears matched-strength.

Outputs:
    data/sparsified_arc/trace_arc/drift/{band}/{patient}.npz  (s, obs, drift)
    data/sparsified_arc/trace_arc/drift/drift_gate_vs_s.csv
    data/sparsified_arc/trace_arc/drift/drift_verdict.csv
"""
from __future__ import annotations
import json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.fc.backbone import percolation_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales

ARC = ROOT / "data" / "sparsified_arc" / "trace_arc"
OUT = ARC / "drift"
NWIN = 5   # audit_167 layout: preA=0 preB=1 task_learn=2 task_test=3 rest_post=4


def _band_abs(Coh, freqs, flo, fhi):        # verbatim audit_167
    m = (freqs >= flo) & (freqs <= fhi)
    A = np.abs(Coh[:, :, m]).mean(axis=-1)
    A = 0.5 * (A + A.T); np.fill_diagonal(A, 0.0)
    return np.clip(A, 0.0, 1.0)


def _eig(W):
    return laplacian_eig(percolation_backbone(W)[0])


def per_cell(job):
    pat, band = job
    cache = ARC / band / f"{pat}.npz"
    if not cache.exists():
        return None
    z = np.load(cache)
    s_grid = z["s"]; obs = z["obs"]
    try:
        X = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    except Exception:
        return None
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nper = nperseg_for_fs(fs); T = X.shape[1]; flo, fhi = BRAIN_BANDS[band]
    Ws = []
    for i in range(NWIN):
        Xi = X[:, i * T // NWIN:(i + 1) * T // NWIN]
        fr, C = compute_msc_welch(Xi, fs, nperseg=min(nper, Xi.shape[1]), metric="imcoh")
        Ws.append(_band_abs(C, fr, flo, fhi))
    drift_eig = {"A": _eig(Ws[0]), "B": _eig(Ws[1]),
                 "task_test": _eig(Ws[3]), "rest_post": _eig(Ws[4])}
    drift = rho_sym_over_scales(drift_eig, s_grid)
    cell = OUT / band; cell.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell / f"{pat}.npz", s=s_grid, obs=obs, drift=drift)
    return dict(patient=pat, band=band,
                s=";".join(f"{v:.3f}" for v in s_grid),
                obs=";".join(f"{v:.4f}" for v in obs),
                drift=";".join(f"{v:.4f}" for v in drift))


def cohort(df):
    vs_s, verdict = [], []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        S = np.array([[float(v) for v in r.split(";")] for r in x.s])
        OBS = np.array([[float(v) for v in r.split(";")] for r in x.obs])
        DFT = np.array([[float(v) for v in r.split(";")] for r in x.drift])
        s_axis = np.nanmedian(S, 0); n_s = s_axis.size
        gate_p = np.full(n_s, np.nan)
        for j in range(n_s):
            d = OBS[:, j] - DFT[:, j]; d = d[np.isfinite(d)]
            if d.size >= 5 and np.any(d != 0):
                try:
                    _, gate_p[j] = wilcoxon(d, alternative="greater")
                except Exception:
                    pass
            vs_s.append(dict(band=band, s=float(s_axis[j]), drift_gate_p=float(gate_p[j]),
                             obs_median=float(np.nanmedian(OBS[:, j])),
                             drift_median=float(np.nanmedian(DFT[:, j]))))
        obs_max = np.nanmax(OBS, 1); dft_max = np.nanmax(DFT, 1)
        try:
            _, p_sm = wilcoxon(obs_max - dft_max, alternative="greater")
        except Exception:
            p_sm = np.nan
        fire = np.isfinite(gate_p) & (gate_p < 0.05)
        verdict.append(dict(band=band, n=len(x),
                            drift_gate_scalemax=float(p_sm),
                            n_scales_beat_drift=int(fire.sum()),
                            frac_beat_drift=round(float(fire.mean()), 3),
                            obs_max_med=float(np.median(obs_max)),
                            drift_max_med=float(np.median(dft_max)),
                            n_pat_obs_gt_drift=int((obs_max > dft_max).sum()),
                            passes_drift=bool(np.isfinite(p_sm) and p_sm < 0.05)))
    return pd.DataFrame(vs_s), pd.DataFrame(verdict)


def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, b) for b in BANDS for p in COHORT]
    if a.limit:
        jobs = jobs[:a.limit]
    ncpu = int(os.environ.get("SA_WORKERS", min(14, (os.cpu_count() or 4) - 2)))
    print(f"[D2-drift] {len(jobs)} cells, windowed rest_pre drift arc, {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
            if i % 10 == 0 or i == len(jobs):
                print(f"[{i}/{len(jobs)}] {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    vs_s, verdict = cohort(df)
    vs_s.to_csv(OUT / "drift_gate_vs_s.csv", index=False)
    verdict.to_csv(OUT / "drift_verdict.csv", index=False)
    print(f"\n[D2-drift] {len(rows)} cells in {time.time()-t0:.0f}s", flush=True)
    print("\n=== DRIFT NULL: does the whole-task trace exceed the drift floor? ===", flush=True)
    print(verdict[["band", "drift_gate_scalemax", "n_scales_beat_drift", "frac_beat_drift",
                   "obs_max_med", "drift_max_med", "n_pat_obs_gt_drift", "passes_drift"]].to_string(index=False),
          flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="D2_drift_null", construction="audit_167 5-window rest_pre drift arc",
        statistic="rho_sym whole-task (VALID); NOT conditional T_infspec_pe",
        gate="cohort Wilcoxon(obs - drift, greater) per scale + scale-max", cohort=COHORT, bands=BANDS,
        note="obs reused from D2 trace_arc cache (bit-identical real arm).",
    ), indent=2))
    print(f"[D2-drift] outputs -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
