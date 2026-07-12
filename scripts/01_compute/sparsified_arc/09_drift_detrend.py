#!/usr/bin/env python3
"""D2-drift-DETREND -- the SNR-clean drift control: does the trace survive removing
the session-long LINEAR FC drift, at full coverage?

The within-rest window drift null (07/08) is SNR-confounded: full-length obs vs
short drift windows says "trace", SNR-matched says "drift", and it can't decide.
This control removes drift a different way -- WITHOUT shortening the signal:

  * tile each phase (rest_pre, task_test, rest_post) with K windows, band-average
    -> a time-resolved edge series over the ordered session (global time t);
  * per edge, OLS slope b_e vs t  = the linear session drift of that edge;
  * base phase FC = mean of that phase's K windows (uses ALL data);
  * detrended phase FC = base - b_e * (t_bar_phase - t_bar_global)  (clip [0,1]) --
    the linear-in-time drift removed, everything else (task-specific reorganisation)
    kept, at full coverage;
  * rho_sym(s) recomputed with the SAME estimator, base vs detrended.

5-point preamble
1. CLAIM: the whole-task trace rho_sym(s) is NOT merely the session-long monotonic
   FC drift that makes task_test & rest_post both "later" than rest_pre.
2. NULL/CONTROL: remove the per-edge linear-in-session-time component; re-measure.
3. STRONGEST ALTERNATIVE: rho_sym>0 because FC drifts monotonically pre->task->post
   (drift aligns D_task-D_pre with D_post-D_pre), no task needed.
4. DOES IT CONTROL IT: yes for the LINEAR drift component, at full coverage (base FC
   uses all data; only the slope is window-estimated) -> neither SNR-matched-conservative
   nor short-window-noisy. CANNOT reject: nonlinear/step nonstationarity coincident with
   the task; does not test node strength (that's the matched-strength rung).
5. FALSIFY: if rho_sym_detrend collapses to ~0 (retained<<1, not sig>0) the trace WAS
   linear drift. If rho_sym_detrend stays ~rho_sym_base and sig>0, it is NOT linear drift.

Reuses the D2 s-grid (trace_arc/{band}/{pat}.npz). Same 1-Welch-bank-per-window /
all-bands optimisation as 08. Deterministic (no resampling) -> fast.

Outputs:
    data/sparsified_arc/trace_arc/drift_detrend/{band}/{patient}.npz
        (s, rho_base, rho_detrend, obs_full)
    data/sparsified_arc/trace_arc/drift_detrend/gate_vs_s.csv
    data/sparsified_arc/trace_arc/drift_detrend/verdict.csv
    data/sparsified_arc/trace_arc/drift_detrend/config.json
"""
from __future__ import annotations
import argparse, json, os, sys, time
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
OUT = ARC / "drift_detrend"
PHASES_ORDERED = ("rest_pre", "task_test", "rest_post")
K = 8   # windows per phase (tile -> time-resolved edge series for the slope)


def _band_abs(Coh, freqs, flo, fhi):
    m = (freqs >= flo) & (freqs <= fhi)
    A = np.abs(Coh[:, :, m]).mean(axis=-1)
    A = 0.5 * (A + A.T); np.fill_diagonal(A, 0.0)
    return np.clip(A, 0.0, 1.0)


def _clean(W):
    W = 0.5 * (W + W.T); np.fill_diagonal(W, 0.0)
    return np.clip(W, 0.0, 1.0)


def _eig_bb(W):
    return laplacian_eig(percolation_backbone(_clean(W))[0])


def per_patient(pat):
    sgrids = {}
    for band in BANDS:
        cache = ARC / band / f"{pat}.npz"
        if not cache.exists():
            return []
        sgrids[band] = np.load(cache)["s"]
    fs = FS_OVERRIDES.get(pat, 2048.0); nper = nperseg_for_fs(fs)

    # windowed band-averaged FCs for every phase (1 Welch/window, all bands derived)
    seqs = {b: {ph: [] for ph in PHASES_ORDERED} for b in BANDS}
    for ph in PHASES_ORDERED:
        X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
        T = X.shape[1]; L = T // K
        for w in range(K):
            Xw = X[:, w * L:(w + 1) * L]
            fr, C = compute_msc_welch(Xw, fs, nperseg=min(nper, Xw.shape[1]), metric="imcoh")
            for b in BANDS:
                flo, fhi = BRAIN_BANDS[b]
                seqs[b][ph].append(_band_abs(C, fr, flo, fhi))
            del C
        del X

    # global time: rest_pre=0..K-1, task=K..2K-1, post=2K..3K-1
    t = np.arange(3 * K, dtype=float); tbar = t.mean()
    tden = ((t - tbar) ** 2).sum()
    tbar_pre = t[:K].mean(); tbar_task = t[K:2 * K].mean(); tbar_post = t[2 * K:].mean()

    rows = []
    for band in BANDS:
        Wpre = np.stack(seqs[band]["rest_pre"])
        Wtask = np.stack(seqs[band]["task_test"])
        Wpost = np.stack(seqs[band]["rest_post"])
        allW = np.concatenate([Wpre, Wtask, Wpost], axis=0)          # (3K,N,N)
        b_e = ((t - tbar)[:, None, None] * (allW - allW.mean(0))).sum(0) / tden  # (N,N) slope
        # base phase FCs (use all windows); A,B = interleaved split-half of rest_pre
        A = Wpre[0::2].mean(0); B = Wpre[1::2].mean(0)
        task = Wtask.mean(0); post = Wpost.mean(0)
        # detrended: remove linear-in-time drift at each phase's mean time
        Ad = A - b_e * (tbar_pre - tbar); Bd = B - b_e * (tbar_pre - tbar)
        td = task - b_e * (tbar_task - tbar); pd_ = post - b_e * (tbar_post - tbar)
        s_grid = sgrids[band]
        try:
            base_eig = {"A": _eig_bb(A), "B": _eig_bb(B),
                        "task_test": _eig_bb(task), "rest_post": _eig_bb(post)}
            detr_eig = {"A": _eig_bb(Ad), "B": _eig_bb(Bd),
                        "task_test": _eig_bb(td), "rest_post": _eig_bb(pd_)}
            rho_base = rho_sym_over_scales(base_eig, s_grid)
            rho_detr = rho_sym_over_scales(detr_eig, s_grid)
        except Exception:
            continue
        obs_full = np.load(ARC / band / f"{pat}.npz")["obs"]
        cell = OUT / band; cell.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cell / f"{pat}.npz", s=s_grid, rho_base=rho_base,
                            rho_detrend=rho_detr, obs_full=obs_full)
        rows.append(dict(patient=pat, band=band,
                         s=";".join(f"{v:.4f}" for v in s_grid),
                         rho_base=";".join(f"{v:.4f}" for v in rho_base),
                         rho_detrend=";".join(f"{v:.4f}" for v in rho_detr)))
    return rows


def cohort(df):
    vs_s, verdict = [], []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        def M(col):
            return np.array([[float(v) for v in r.split(";")] for r in x[col]])
        S = M("s"); BASE = M("rho_base"); DETR = M("rho_detrend")
        s_axis = np.nanmedian(S, 0); n_s = s_axis.size
        p_pos = np.full(n_s, np.nan)   # detrended still > 0 ?
        p_drop = np.full(n_s, np.nan)  # base > detrended ? (did detrend reduce it)
        for j in range(n_s):
            d = DETR[:, j]; d = d[np.isfinite(d)]
            if d.size >= 5 and np.any(d != 0):
                try:
                    p_pos[j] = wilcoxon(d, alternative="greater")[1]
                except Exception:
                    pass
            dd = BASE[:, j] - DETR[:, j]; dd = dd[np.isfinite(dd)]
            if dd.size >= 5 and np.any(dd != 0):
                try:
                    p_drop[j] = wilcoxon(dd, alternative="greater")[1]
                except Exception:
                    pass
            vs_s.append(dict(band=band, s=float(s_axis[j]),
                             rho_base_med=float(np.nanmedian(BASE[:, j])),
                             rho_detrend_med=float(np.nanmedian(DETR[:, j])),
                             p_detrend_pos=float(p_pos[j]), p_base_gt_detrend=float(p_drop[j])))
        js = int(np.nanargmax(np.nanmedian(BASE, 0)))   # emergence scale (base)
        base_med = np.nanmedian(BASE, 0); detr_med = np.nanmedian(DETR, 0)
        retained = np.where(np.abs(base_med) > 1e-6, detr_med / base_med, np.nan)
        surv = np.isfinite(p_pos) & (p_pos < 0.05)
        verdict.append(dict(band=band, n=len(x),
                            s_emerge=float(s_axis[js]),
                            rho_base_emerge=float(base_med[js]),
                            rho_detrend_emerge=float(detr_med[js]),
                            retained_emerge=float(retained[js]),
                            p_detrend_pos_emerge=float(p_pos[js]),
                            n_scales_detrend_sig=int(surv.sum()),
                            survives_detrend=bool(surv.any())))
    return pd.DataFrame(vs_s), pd.DataFrame(verdict)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    pats = COHORT[:a.limit] if a.limit else COHORT
    ncpu = int(os.environ.get("SA_WORKERS", 4))
    print(f"[D2-detrend] {len(pats)} patients x {len(BANDS)} bands | K={K} windows/phase, "
          f"linear session-drift removed at full coverage | {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, rlist in enumerate(pool.imap_unordered(per_patient, pats), 1):
            rows.extend(rlist)
            el = time.time() - t0
            print(f"[{i}/{len(pats)} pt] {el:.0f}s  ETA {el/i*(len(pats)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        print("[D2-detrend] no cells", flush=True); return
    vs_s, verdict = cohort(df)
    vs_s.to_csv(OUT / "gate_vs_s.csv", index=False)
    verdict.to_csv(OUT / "verdict.csv", index=False)
    print(f"\n[D2-detrend] {len(rows)} cells in {time.time()-t0:.0f}s\n", flush=True)
    print("=== DOES THE TRACE SURVIVE LINEAR-DRIFT REMOVAL (full coverage)? ===", flush=True)
    print(verdict[["band", "s_emerge", "rho_base_emerge", "rho_detrend_emerge",
                   "retained_emerge", "p_detrend_pos_emerge", "n_scales_detrend_sig",
                   "survives_detrend"]].to_string(index=False), flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="D2_drift_detrend", control="per-edge linear session-drift removal, full coverage",
        K_windows_per_phase=K, test="rho_sym_detrend still > 0 (survives) vs collapses (was drift)",
        cohort=COHORT, bands=BANDS,
        note="base FC = mean of K windows (all data); only slope is window-estimated.",
    ), indent=2))
    print(f"[D2-detrend] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
