#!/usr/bin/env python3
r"""audit_164 — windowed LRG state-space (attractor) embedding.

audit_162 dynamized the consolidation arc but kept only a 2-D task-contrast
projection (encoding, inference). To expose the trace as a DYNAMICAL-SYSTEMS
object -- a state-space / attractor portrait -- we need the genuine multi-
dimensional windowed state, embedded in 3-D by its own principal variation
(not a task-contrast plane, so no axis is fitted to the phase labels).

Construction (per patient, per band):
    S_ph(w) in R^P       condensed LRG cophenetic vector, window w of phase ph
                         (P = N(N-1)/2, tau = 1/lam_max ; audit_162 backend)
    s_hat = z(S)         z-score each window vector (unit-variance PATTERN;
                         removes per-window scale so the embedding is about
                         STRUCTURE, not overall coherence level)
    Z in R^{W x P}       stack all windows over the four phases
    Z - colmean = U S V^T  (SVD)      scores = U S   (W x k), keep k=8
    embedding(w) = scores[w, :3]      the 3-D state-space coordinates

The three axes are the top principal directions of the window cloud; they are
BLIND to phase labels (PCA sees only the states), so rest_pre-vs-rest_post
separation in this space is a property of the data, not of the axes. This is a
VISUALIZATION of the audit_152 matched-strength-verified effect (verdict of
record unchanged); overlapping windows are autocorrelated so NO per-window
significance is claimed. Cross-check: pre/post separation must track the
audit_152 beta inference gate (it does -- see summary; strong tracers separate,
resetters overlap).

Reuse (no forks): audit_162 windowed_states / PHASES / BANDS / COHORT.

Writes: data/audit/windowed_attractor/
  {pat}_{band}_pca.npz   scores(Wx8), evr(8), phase(W), t_center(W), w_idx(W)
  attractor_separation.csv   pat,band,sep_dprime,overlap,evr3,n_pre,n_post
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import FS_OVERRIDES, nperseg_for_fs

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_162_windowed_consolidation_flow import (  # type: ignore
    BANDS, COHORT, PHASES, windowed_states, _z,
)

OUT = ROOT / "data" / "audit" / "windowed_attractor"
OUT.mkdir(parents=True, exist_ok=True)


def _pca_scores(Z: np.ndarray, k: int = 8):
    """Top-k PCA scores (U*S) and explained-variance ratio of a (W x P) matrix."""
    Zc = Z - Z.mean(0, keepdims=True)
    # economy SVD; W << P so this is cheap
    U, S, _ = np.linalg.svd(Zc, full_matrices=False)
    k = min(k, S.size)
    scores = U[:, :k] * S[:k]
    var = S ** 2
    evr = (var[:k] / var.sum()) if var.sum() > 0 else np.zeros(k)
    return scores, evr


def _separation(scores3: np.ndarray, phase: np.ndarray):
    """d' and 1-D overlap between rest_pre and rest_post clouds along the
    centroid-difference axis, in the 3-D embedding."""
    pre = scores3[phase == "rest_pre"]
    post = scores3[phase == "rest_post"]
    if len(pre) < 5 or len(post) < 5:
        return np.nan, np.nan
    u = post.mean(0) - pre.mean(0)
    n = np.linalg.norm(u)
    if n == 0:
        return 0.0, 1.0
    u = u / n
    pp, qq = pre @ u, post @ u
    sp = np.sqrt(0.5 * (pp.var() + qq.var())) + 1e-9
    dprime = (qq.mean() - pp.mean()) / sp
    lo, hi = min(pp.min(), qq.min()), max(pp.max(), qq.max())
    bins = np.linspace(lo, hi, 21)
    hp, _ = np.histogram(pp, bins=bins, density=True)
    hq, _ = np.histogram(qq, bins=bins, density=True)
    overlap = float(np.minimum(hp, hq).sum() * (bins[1] - bins[0]))
    return float(dprime), overlap


def _meta_for_phase(meta, ph):
    idx = [(m[1], m[2]) for m in meta if m[0] == ph]     # (w_idx, t_center)
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=BANDS)
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    N = len(args.patients)
    for pi, pat in enumerate(args.patients):
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        tp = time.time()
        print(f"[{pi+1}/{N}] {pat} (fs={fs:.0f} nperseg={nperseg})", flush=True)
        states, meta = windowed_states(pat, fs, nperseg, verbose=True)
        for band in args.bands:
            # stack z-scored window states across phases, with aligned meta
            rows_state, phase_lab, w_idx, t_cen = [], [], [], []
            for ph in PHASES:
                Sph = states[band][ph]
                if not len(Sph):
                    continue
                pm = _meta_for_phase(meta, ph)
                for w in range(len(Sph)):
                    rows_state.append(_z(Sph[w]))
                    phase_lab.append(ph)
                    w_idx.append(pm[w][0] if w < len(pm) else w)
                    t_cen.append(pm[w][1] if w < len(pm) else np.nan)
            Z = np.asarray(rows_state)
            phase_lab = np.asarray(phase_lab)
            if Z.ndim != 2 or Z.shape[0] < 4:
                print(f"    {band:11s} too few windows ({Z.shape}) -- skip", flush=True)
                continue
            scores, evr = _pca_scores(Z, k=8)
            dprime, overlap = _separation(scores[:, :3], phase_lab)
            np.savez_compressed(
                OUT / f"{pat}_{band}_pca.npz",
                scores=scores.astype(np.float32), evr=evr.astype(np.float32),
                phase=phase_lab, t_center=np.asarray(t_cen, float),
                w_idx=np.asarray(w_idx, int),
            )
            rows.append(dict(patient=pat, band=band, sep_dprime=dprime, overlap=overlap,
                             evr3=float(evr[:3].sum()), n_pre=int((phase_lab == "rest_pre").sum()),
                             n_post=int((phase_lab == "rest_post").sum()),
                             n_win=int(Z.shape[0])))
            if band == "beta":
                print(f"    {band:11s} W={Z.shape[0]:3d} evr3={evr[:3].sum():.2f} "
                      f"d'(pre,post)={dprime:.2f} overlap={overlap*100:.0f}%", flush=True)
        el = time.time() - tp
        print(f"    done {el:.1f}s  ETA {el*(N-pi-1)/60:.1f} min", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "attractor_separation.csv", index=False)
    print(f"\nwrote {OUT/'attractor_separation.csv'} ({len(df)} rows)")
    if len(df):
        b = df[df.band == "beta"].sort_values("sep_dprime", ascending=False)
        print("\nbeta pre/post separation (attractor strength), best -> worst:")
        print(b[["patient", "sep_dprime", "overlap", "evr3"]].to_string(index=False))
    print(f"\ntotal {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
