#!/usr/bin/env python3
"""Audit 125 — N4 replay states: per-window task-likeness + BURSTINESS make-or-break.

Computes, per non-overlapping window of a phase, the LRG cophenetic config
c[x,w,b] from windowed |ImCoh| (recipe == N1, anchored in audit_124), and derives
THREE task-likeness detectors of the window vs the task-state config (T=task_test),
referenced to rest_pre:

    d_w = c_w - cbar_pre ,  d_T = cbar_T - cbar_pre
    (1) GLOBAL SPEARMAN  s = spearman(d_w, d_T)              # the scoped measure
    (2) GLOBAL COSINE    cos(d_w, d_T)
    (3) LOCALIZED COSINE cos over the top-K |d_T| pairs       # task-signature subspace

Per phase the multiband score S[x,w] = mean_b zscore_within_phase(detector). A
REPLAY STATE is a transient burst in S. The make-or-break test (scope N4.2): does
S show temporal structure beyond a STATIONARY process? Tested by:

  - lag-1 autocorrelation AC1(S) vs a TIME-SHUFFLE null (random permutation;
    valid because windows are NON-overlapping) -> p_time
  - longest run above the within-phase median (dwell) vs the same null

Time-shuffle preserves the marginal of S exactly (same mean, same event RATE), so
only genuine temporal CLUSTERING can beat it. Elevated mean alone is conceded to
N1 and is NOT the N4 claim.

Specificity (scope N_altconfig, alternative B = arousal/drift toward ANY config):
AC1 for the real task target vs an ensemble of SHUFFLED-target directions. Real
must exceed shuffled.

Consolidation (scope N4.4): cross-phase-comparable RAW detector mean + event rate
at a rest_pre-derived threshold; rest_post vs rest_pre.

Caches per-window cophenetic vectors (the expensive part) so cohort scale-up and
matched-strength nulls reuse them.

Outputs
-------
    data/audit/replay_states/
        per_window.csv              # every window's detectors, all phases
        burstiness_summary.csv      # AC1 + p_time + specificity, per patient/phase/detector
        cache/<pat>_<band>_<phase>_coph.npz
        figures/timecourse_<pat>.pdf
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import segment_ffts, imcoh_abs_cube, band_abs_average
from lrg_eegfc.utils.surrogate.matched_strength import cophenetic_condensed_from_adjacency
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()

OUT = ROOT / "data" / "audit" / "replay_states"
FIG = OUT / "figures"
CACHE = OUT / "cache"
for d in (OUT, FIG, CACHE):
    d.mkdir(parents=True, exist_ok=True)

PHASES = ["rest_pre", "task_test", "rest_post"]


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def cosine(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    if na == 0 or nb == 0:
        return np.nan
    return float(np.dot(a, b) / (na * nb))


def ac1(x):
    """Lag-1 autocorrelation (mean-removed, normalized)."""
    x = np.asarray(x, float)
    x = x - x.mean()
    denom = np.dot(x, x)
    if denom == 0:
        return 0.0
    return float(np.dot(x[:-1], x[1:]) / denom)


def longest_run_above_median(x):
    b = np.asarray(x) > np.median(x)
    best = cur = 0
    for v in b:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return int(best)


def perm_pvalue(x, stat_fn, n_perm, rng):
    obs = stat_fn(x)
    null = np.array([stat_fn(rng.permutation(x)) for _ in range(n_perm)])
    p = (1 + np.sum(null >= obs)) / (n_perm + 1)
    return float(obs), float(p), null


def compute_window_cophenets(pat, phase, bands, L_sec, fs, nperseg, max_windows):
    """Return {band: (coph_array (n_win, P), times (n_win,))} for non-overlapping
    windows; cache to disk."""
    X = load_timeseries(pat, phase, SEEG_DATAPATH)
    X = np.asarray(X, float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    N, T = X.shape
    wlen = int(round(L_sec * fs))
    starts = list(range(0, T - wlen + 1, wlen))
    if max_windows and len(starts) > max_windows:
        starts = list(np.linspace(0, T - wlen, max_windows, dtype=int))
    fmax = max(BRAIN_BANDS[b][1] for b in bands) + 5.0
    out = {b: [] for b in bands}
    times = []
    for s in starts:
        win = X[:, s:s + wlen]
        freqs, ff = segment_ffts(win, fs, nperseg, fmax_keep=fmax)
        if ff.shape[1] < 4:
            continue
        cube = imcoh_abs_cube(ff)
        for b in bands:
            lo, hi = BRAIN_BANDS[b]
            W = band_abs_average(cube, freqs, lo, hi)
            out[b].append(coph(W))
        times.append((s + wlen / 2) / fs)
    res = {}
    for b in bands:
        arr = np.array(out[b])
        np.savez_compressed(CACHE / f"{pat}_{b}_{phase}_coph.npz",
                            coph=arr.astype(np.float32), times=np.array(times),
                            L_sec=L_sec, N=N)
        res[b] = (arr, np.array(times))
    return res


def detectors_for_phase(coph_by_band, refs, topk_frac):
    """Return dict detector -> {band: array(n_win)} for one phase."""
    bands = list(coph_by_band.keys())
    det = {"spear_global": {}, "cos_global": {}, "cos_local": {}}
    for b in bands:
        arr, _ = coph_by_band[b]          # (n_win, P)
        cbar_pre, cbar_T = refs[b]
        dT = cbar_T - cbar_pre
        k = max(50, int(topk_frac * dT.size))
        idx = np.argsort(np.abs(dT))[-k:]
        dT_loc = dT[idx]
        sg, cg, cl = [], [], []
        for c in arr:
            dw = c - cbar_pre
            sg.append(spearmanr(dw, dT)[0])
            cg.append(cosine(dw, dT))
            cl.append(cosine(dw[idx], dT_loc))
        det["spear_global"][b] = np.array(sg)
        det["cos_global"][b] = np.array(cg)
        det["cos_local"][b] = np.array(cl)
    return det


def zmean_multiband(det_band_dict):
    """mean over bands of within-phase z-scored detector -> (n_win,)."""
    mats = []
    for b, v in det_band_dict.items():
        v = np.asarray(v, float)
        sd = v.std()
        mats.append((v - v.mean()) / sd if sd > 0 else v * 0.0)
    return np.mean(mats, axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=["Pat_05", "Pat_02", "Pat_13"])
    ap.add_argument("--bands", nargs="+", default=["alpha", "beta"])
    ap.add_argument("--L-sec", type=float, default=20.0)
    ap.add_argument("--max-windows", type=int, default=0, help="0 = all non-overlapping")
    ap.add_argument("--topk-frac", type=float, default=0.10)
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--n-shuf-target", type=int, default=200)
    ap.add_argument("--target-phase", default="task_test")
    ap.add_argument("--seed", type=int, default=20260622)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    DET_NAMES = ["spear_global", "cos_global", "cos_local"]
    t0 = time.time()
    win_rows = []
    burst_rows = []

    for pat in args.patients:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        print(f"\n=== {pat} (fs={fs:.0f}, L={args.L_sec:.0f}s) ===")
        # references
        refs = {}
        ok = True
        for b in args.bands:
            Wpre = load_fc_matrix(pat, "rest_pre", b, fc_method="imcoh_abs")
            WT = load_fc_matrix(pat, args.target_phase, b, fc_method="imcoh_abs")
            if Wpre is None or WT is None:
                ok = False; break
            refs[b] = (coph(Wpre), coph(WT))
        if not ok:
            print(f"  [skip] missing refs"); continue

        det_by_phase = {}
        coph_cache = {}
        for ph in PHASES:
            cb = compute_window_cophenets(pat, ph, args.bands, args.L_sec, fs, nperseg, args.max_windows)
            coph_cache[ph] = cb
            det = detectors_for_phase(cb, refs, args.topk_frac)
            det_by_phase[ph] = det
            nwin = len(next(iter(cb.values()))[1])
            print(f"  {ph}: {nwin} windows")
            # record per-window detectors (multiband mean, raw)
            for dn in DET_NAMES:
                raw = np.mean([det[dn][b] for b in args.bands], axis=0)
                times = next(iter(cb.values()))[1]
                for w in range(len(raw)):
                    win_rows.append(dict(patient=pat, phase=ph, detector=dn,
                                         w=w, t_sec=float(times[w]), value=float(raw[w])))

        # ---- burstiness make-or-break per detector, per phase ----
        for dn in DET_NAMES:
            # specificity: build shuffled-target detector AC1 ensemble on rest_post
            for ph in PHASES:
                S = zmean_multiband(det_by_phase[ph][dn])
                if len(S) < 8:
                    continue
                ac, p_ac, _ = perm_pvalue(S, ac1, args.n_perm, rng)
                lr, p_lr, _ = perm_pvalue(S, longest_run_above_median, args.n_perm, rng)
                # specificity only meaningful in rest phases
                spec_pct = np.nan
                if ph == "rest_post":
                    shuf_acs = []
                    for _ in range(args.n_shuf_target):
                        # shuffle target direction pairs per band -> recompute detector
                        det_shuf = {}
                        for b in args.bands:
                            arr, _ = coph_cache[ph][b]
                            cbar_pre, cbar_T = refs[b]
                            dT = cbar_T - cbar_pre
                            perm = rng.permutation(dT.size)
                            dT_s = dT[perm]
                            if dn == "spear_global":
                                vals = [spearmanr(c - cbar_pre, dT_s)[0] for c in arr]
                            elif dn == "cos_global":
                                vals = [cosine(c - cbar_pre, dT_s) for c in arr]
                            else:
                                k = max(50, int(args.topk_frac * dT.size))
                                idx = np.argsort(np.abs(dT_s))[-k:]
                                vals = [cosine((c - cbar_pre)[idx], dT_s[idx]) for c in arr]
                            det_shuf[b] = np.array(vals)
                        shuf_acs.append(ac1(zmean_multiband(det_shuf)))
                    shuf_acs = np.array(shuf_acs)
                    spec_pct = float(np.mean(shuf_acs < ac))  # frac of shuffled below real
                burst_rows.append(dict(patient=pat, phase=ph, detector=dn, n_win=len(S),
                                       AC1=ac, p_time_AC1=p_ac,
                                       longest_run=lr, p_time_run=p_lr,
                                       specificity_pct=spec_pct))
                tag = "<<<" if (p_ac < 0.05 and ph == "rest_post") else ""
                print(f"    [{dn:12s} {ph:9s}] AC1={ac:+.3f} p_time={p_ac:.3f}  "
                      f"run={lr} p={p_lr:.3f}  spec={spec_pct if np.isnan(spec_pct) else f'{spec_pct:.2f}'} {tag}")

        # ---- consolidation contrast: raw detector mean + event rate ----
        for dn in DET_NAMES:
            raw = {ph: np.mean([det_by_phase[ph][dn][b] for b in args.bands], axis=0) for ph in PHASES}
            thr = np.quantile(raw["rest_pre"], 0.75)
            R_pre = float(np.mean(raw["rest_pre"] > thr))
            R_post = float(np.mean(raw["rest_post"] > thr))
            print(f"    [{dn:12s} CONSOLID] mean_pre={raw['rest_pre'].mean():+.3f} "
                  f"mean_post={raw['rest_post'].mean():+.3f} mean_task={raw['task_test'].mean():+.3f}  "
                  f"R_pre={R_pre:.2f} R_post={R_post:.2f}")
            burst_rows.append(dict(patient=pat, phase="consolidation", detector=dn, n_win=0,
                                   AC1=np.nan, p_time_AC1=np.nan, longest_run=0, p_time_run=np.nan,
                                   specificity_pct=np.nan,
                                   mean_pre=float(raw["rest_pre"].mean()),
                                   mean_post=float(raw["rest_post"].mean()),
                                   mean_task=float(raw["task_test"].mean()),
                                   R_pre=R_pre, R_post=R_post))

        # ---- figure: timecourses (z-scored S) for the scoped detector ----
        fig, axes = plt.subplots(len(DET_NAMES), 1, figsize=(8.5, 7.2), sharex=False)
        for ax, dn in zip(axes, DET_NAMES):
            for ph, col in zip(["rest_pre", "rest_post"], ["#888888", "#1f3d6e"]):
                S = zmean_multiband(det_by_phase[ph][dn])
                ax.plot(np.arange(len(S)), S, "-o", ms=3, color=col, label=ph)
                thr = np.median(S)
                ax.fill_between(np.arange(len(S)), thr, S, where=S > thr, color=col, alpha=0.15)
            ax.axhline(0, color="0.8", lw=0.6)
            ax.set_ylabel(dn, fontsize=8)
            ax.legend(frameon=False, fontsize=7, ncol=2)
            ax.spines[["top", "right"]].set_visible(False)
        axes[-1].set_xlabel("window index")
        fig.tight_layout()
        fig.savefig(FIG / f"timecourse_{pat}.pdf", transparent=True)
        plt.close(fig)

    pd.DataFrame(win_rows).to_csv(OUT / "per_window.csv", index=False)
    bdf = pd.DataFrame(burst_rows)
    bdf.to_csv(OUT / "burstiness_summary.csv", index=False)

    print("\n=== BURSTINESS (rest_post) — AC1 vs time-shuffle null, by detector ===")
    sub = bdf[(bdf.phase == "rest_post")][["patient", "detector", "AC1", "p_time_AC1", "specificity_pct"]]
    print(sub.to_string(index=False))

    # ---- COHORT-LEVEL verdict (the actual verification) ----
    from scipy.stats import wilcoxon
    print("\n=== COHORT VERDICT (n={}) ===".format(bdf.patient.nunique()))
    coh_rows = []
    for dn in DET_NAMES:
        post = bdf[(bdf.phase == "rest_post") & (bdf.detector == dn)].set_index("patient")
        pre = bdf[(bdf.phase == "rest_pre") & (bdf.detector == dn)].set_index("patient")
        pats = [p for p in post.index if p in pre.index]
        if len(pats) < 5:
            continue
        ac_post = post.loc[pats, "AC1"].to_numpy()
        ac_pre = pre.loc[pats, "AC1"].to_numpy()
        n_sig = int((post.loc[pats, "p_time_AC1"] < 0.05).sum())
        # is rest_post AC1 > 0 across cohort?
        try:
            _, p_pos = wilcoxon(ac_post, alternative="greater")
        except Exception:
            p_pos = np.nan
        # is rest_post AC1 > rest_pre AC1 (paired, consolidation of burstiness)?
        try:
            _, p_paired = wilcoxon(ac_post - ac_pre, alternative="greater")
        except Exception:
            p_paired = np.nan
        # LOO-max of the paired statistic (drop-one robustness)
        loo_p = []
        for i in range(len(pats)):
            d = np.delete(ac_post - ac_pre, i)
            try:
                loo_p.append(wilcoxon(d, alternative="greater")[1])
            except Exception:
                loo_p.append(np.nan)
        loo_max = float(np.nanmax(loo_p)) if loo_p else np.nan
        # consolidation means
        cons = bdf[(bdf.phase == "consolidation") & (bdf.detector == dn)].set_index("patient")
        cpats = [p for p in cons.index if p in pats]
        dmean = float(np.median((cons.loc[cpats, "mean_post"] - cons.loc[cpats, "mean_pre"]).to_numpy())) if cpats else np.nan
        coh_rows.append(dict(detector=dn, n=len(pats),
                             median_AC1_post=float(np.median(ac_post)),
                             median_AC1_pre=float(np.median(ac_pre)),
                             n_patients_sig=f"{n_sig}/{len(pats)}",
                             wilcoxon_p_AC1_pos=p_pos,
                             wilcoxon_p_post_gt_pre=p_paired,
                             LOO_max_p_post_gt_pre=loo_max,
                             median_dmean_post_minus_pre=dmean))
    cohdf = pd.DataFrame(coh_rows)
    cohdf.to_csv(OUT / "cohort_verdict.csv", index=False)
    print(cohdf.to_string(index=False))
    print(f"\n[audit_125] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
