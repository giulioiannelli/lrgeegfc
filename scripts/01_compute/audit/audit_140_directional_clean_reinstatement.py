#!/usr/bin/env python3
"""audit_140 — N4 v3: is the DIRECTIONAL sustained reinstatement real BEYOND magnitude?

Scope: .agents/guides/task-persistence-investigation/2026-06-22_replay-states-v3-scale-target-pair.md

P1 (audit_136) found the magnetic-Laplacian reinstatement LEVEL positive (rest_post holds
the task's directional-flow geometry). Caveat flagged: |e^{−τL_H}| mixes magnitude (D̄=Σ|A|)
AND direction (iA), so the shift could merely inherit the already-known (very strong)
MAGNITUDE reinstatement. This isolates direction with a magnitude-matched control:

  signed band-averaged ImCoh A (antisymmetric).
  DIRECTIONAL operator : L_H  = D̄ − iA      (magnetic; carries lead/lag)
  MAGNITUDE-matched    : L_mag = D̄ − |A|     (SAME D̄; sign/direction removed)
  Only the off-diagonal (iA vs |A|) differs ⇒ the gap isolates DIRECTION.

Per patient, band, τ=α/λmax: cophenetic of each operator's propagator; task-likeness
  s_X(w) = ρcoph(c^X_w, c^X_T) − ρcoph(c^X_w, c^X_{R-}),  X ∈ {dir, mag}
  shift_X = mean_post s_X − mean_pre s_X
  directional_excess = shift_dir − shift_mag   (>0 ⇒ direction reinstated BEYOND magnitude)

Cohort one-sided Wilcoxon + LOO-max: (a) shift_dir>0 (directional reinstatement exists),
(b) directional_excess>0 (THE clean claim), (c) shift_mag>0 (magnitude sanity). A clean
positive on (b) = a genuine NEW result: rest_post reinstates the task's information-flow
DIRECTION, not just which contacts couple. (Sustained, not a transient replay state.)

Output: data/audit/replay_states/directional_clean_{per_patient,cohort}.csv
"""
from __future__ import annotations

import argparse
import time
import warnings

import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import segment_ffts, imcoh_signed_cube, band_signed_average
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrgsglib.core import compute_normalized_linkage, extract_ultrametric_matrix

OUT = ROOT / "data" / "audit" / "replay_states"
OUT.mkdir(parents=True, exist_ok=True)
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
ALPHAS = [1.0, 2.0, 4.0]
VAR_FLOOR = 1e-10
_EMPTY_G: dict[int, nx.Graph] = {}


def _giant(n):
    if n not in _EMPTY_G:
        _EMPTY_G[n] = nx.empty_graph(n)
    return _EMPTY_G[n]


def mag_eig(A):
    """Magnitude-matched Laplacian L_mag = D̄ − |A| (D̄=Σ|A|): real eig."""
    B = np.abs(0.5 * (np.asarray(A, float) - np.asarray(A, float).T))
    np.fill_diagonal(B, 0.0)
    L = np.diag(B.sum(1)) - B
    return np.linalg.eigh(L)


def dir_eig(A):
    """Magnetic Laplacian L_H = D̄ − iA (same D̄): real eigvals, complex eigvecs."""
    A = 0.5 * (np.asarray(A, float) - np.asarray(A, float).T)
    np.fill_diagonal(A, 0.0)
    dbar = np.abs(A).sum(1)
    return np.linalg.eigh(np.diag(dbar).astype(complex) - 1j * A)


def D_raw(ev, U, tau, complex_U=False):
    diag = np.exp(-tau * ev)
    P = (U * diag) @ (U.conj().T if complex_U else U.T)
    rho = np.abs(P) if complex_U else P
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    fin = np.isfinite(D)
    if not fin.all():
        D = np.where(fin, D, np.nanmax(D[fin]) if fin.any() else 1e12)
    return D


def D_coph(D, n):
    Zlink, _, _ = compute_normalized_linkage(squareform(D), _giant(n))
    return squareform(extract_ultrametric_matrix(Zlink, n))


def rho_coph(a, b):
    if np.nanvar(a) < VAR_FLOOR or np.nanvar(b) < VAR_FLOOR:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spearmanr(a, b).statistic)


def windowed_signed_A(X, fs, L_sec, nperseg, band):
    wlen = int(round(L_sec * fs))
    fmax = BRAIN_BANDS[band][1] + 5.0
    out = []
    for s in range(0, X.shape[1] - wlen + 1, wlen):
        freqs, ff = segment_ffts(X[:, s:s + wlen], fs, nperseg, nperseg // 2, fmax)
        if ff.shape[1] < 4:
            continue
        out.append(band_signed_average(imcoh_signed_cube(ff), freqs, *BRAIN_BANDS[band]))
    return out


def coph_legs(A, n, a):
    """(dir, mag) condensed cophenetic at τ=a/λmax for one signed matrix A."""
    evd, Ud = dir_eig(A); evm, Um = mag_eig(A)
    cd = D_coph(D_raw(evd, Ud, a / evd[-1], complex_U=True), n)
    cm = D_coph(D_raw(evm, Um, a / evm[-1], complex_U=False), n)
    return cd, cm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=["alpha", "beta"])
    ap.add_argument("--L-sec", type=float, default=20.0)
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    for pat in args.patients:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        try:
            def _load(ph):
                Xx = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
                return Xx if Xx.shape[0] <= Xx.shape[1] else Xx.T
            Xpre, Xpost = _load("rest_pre"), _load("rest_post")
        except Exception as e:
            print(f"[skip] {pat}: {e}"); continue
        for band in args.bands:
            At = load_fc_matrix(pat, "task_test", band, fc_method="imcoh")
            Ap = load_fc_matrix(pat, "rest_pre", band, fc_method="imcoh")
            if At is None or Ap is None:
                print(f"  [skip] {pat} {band}: signed refs missing"); continue
            At = np.asarray(At, float); Ap = np.asarray(Ap, float); n = At.shape[0]
            Aw_pre = windowed_signed_A(Xpre, fs, args.L_sec, nperseg, band)
            Aw_post = windowed_signed_A(Xpost, fs, args.L_sec, nperseg, band)
            for a in ALPHAS:
                cdT, cmT = coph_legs(At, n, a)
                cdP, cmP = coph_legs(Ap, n, a)
                wd_pre = [coph_legs(A, n, a) for A in Aw_pre]
                wd_post = [coph_legs(A, n, a) for A in Aw_post]
                s_dir_pre = np.array([rho_coph(cd, cdT) - rho_coph(cd, cdP) for cd, cm in wd_pre])
                s_dir_post = np.array([rho_coph(cd, cdT) - rho_coph(cd, cdP) for cd, cm in wd_post])
                s_mag_pre = np.array([rho_coph(cm, cmT) - rho_coph(cm, cmP) for cd, cm in wd_pre])
                s_mag_post = np.array([rho_coph(cm, cmT) - rho_coph(cm, cmP) for cd, cm in wd_post])
                shift_dir = float(np.nanmean(s_dir_post) - np.nanmean(s_dir_pre))
                shift_mag = float(np.nanmean(s_mag_post) - np.nanmean(s_mag_pre))
                rows.append(dict(patient=pat, band=band, alpha=a,
                                 shift_dir=shift_dir, shift_mag=shift_mag,
                                 directional_excess=shift_dir - shift_mag))
            print(f"  {pat} {band}: {len(Aw_post)} post-win ({time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "directional_clean_per_patient.csv", index=False)

    def wt(d):
        d = np.asarray(d, float); d = d[np.isfinite(d)]
        if d.size < 5 or np.allclose(d, 0):
            return np.nan, 0, len(d), np.nan
        p = float(wilcoxon(d, alternative="greater").pvalue)
        loo = float(np.nanmax([wilcoxon(np.delete(d, i), alternative="greater").pvalue
                               for i in range(len(d))]))
        return p, int((d > 0).sum()), len(d), loo

    coh = []
    for band in args.bands:
        for a in ALPHAS:
            sub = df[(df.band == band) & (df.alpha == a)]
            if sub.empty:
                continue
            for label, col in [("dir_reinstated", "shift_dir"),
                               ("DIRECTIONAL_excess(beyond_mag)", "directional_excess"),
                               ("mag_reinstated(sanity)", "shift_mag")]:
                p, npos, K, loo = wt(sub[col])
                coh.append(dict(band=band, alpha=a, test=label,
                                median=float(sub[col].median()), n_pos=f"{npos}/{K}",
                                wilcoxon_p=p, loo_max_p=loo))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "directional_clean_cohort.csv", index=False)
    print("\n=== directional reinstatement BEYOND magnitude — cohort ===")
    for band in args.bands:
        for a in ALPHAS:
            blk = cohdf[(cohdf.band == band) & (cohdf.alpha == a)]
            if blk.empty:
                continue
            print(f"\n  [{band} | τ=α/λmax, α={a}]")
            for _, r in blk.iterrows():
                star = "*" if (r.wilcoxon_p < 0.05 and r.loo_max_p < 0.05) else " "
                print(f"   {star}{r.test:32s} med={r['median']:+.4f} {r.n_pos:>5} "
                      f"p={r.wilcoxon_p:.3f} LOO={r.loo_max_p:.3f}")
    print(f"\n[audit_140] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
