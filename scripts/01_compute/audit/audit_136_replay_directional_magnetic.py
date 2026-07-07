#!/usr/bin/env python3
"""audit_136 — N4 v3 / P1: DIRECTIONAL replay via the MAGNETIC Laplacian propagator.

Scope: .agents/guides/task-persistence-investigation/2026-06-22_replay-states-v3-scale-target-pair.md

All prior replay tests used |ImCoh| — UNDIRECTED coupling magnitude — and found no
transient state at any scale (audit_133). The SIGNED ImCoh keeps phase lead/lag = WHO
DRIVES WHOM. P1 asks whether rest_post reinstates the task's DIRECTIONAL FLOW geometry,
the closest in-framework analog to the SEQUENCE that replay fundamentally is — staying
strictly Laplacian-propagator based.

Operator: signed band-averaged ImCoh A (real, antisymmetric) → magnetic Laplacian
  L_H = D̄ − iA,  D̄ = diag(Σ_j|A_ij|)   (Hermitian, PSD; gauge-clean across phases since
  same node order/charge is used for window AND references → convention cancels in the
  window-vs-task similarity). Propagator P_H(τ) = V e^{−τλ} V†; |P_H| is symmetric ≥0 →
  communication distance Δ = 1/|P_H| → cophenetic Dᶜ_H(τ). Everything else == audit_133:
  per τ=α/λmax, task-likeness contrast sᵂ(τ)=ρcoph(Dᶜ_w,Dᶜ_T)−ρcoph(Dᶜ_w,Dᶜ_{R-});
  shift (directional N1?), burstiness vs time-shuffle (directional replay STATE?),
  node-permuted-target PLACEBO (coarse-τ collapse control).

Two possible wins: a directional SUSTAINED reinstatement (new positive even w/o transient)
OR a directional transient state. Caveat: band-averaging signed ImCoh nets sub-band sign
flips (lossy by design; same as the signed-flow audits).

Output: data/audit/replay_states/directional_{per_patient,cohort}.csv
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
ALPHA = np.unique(np.r_[0.5, 1.0, np.geomspace(1.0, 40.0, 11)]).round(4)
VAR_FLOOR = 1e-10
_EMPTY_G: dict[int, nx.Graph] = {}


def _giant(n):
    if n not in _EMPTY_G:
        _EMPTY_G[n] = nx.empty_graph(n)
    return _EMPTY_G[n]


def signed_eig(A):
    """Magnetic Laplacian L_H = D̄ − iA of antisymmetric A → (real eigvals, complex eigvecs)."""
    A = 0.5 * (np.asarray(A, float) - np.asarray(A, float).T)
    np.fill_diagonal(A, 0.0)
    dbar = np.abs(A).sum(1)
    L_H = np.diag(dbar).astype(complex) - 1j * A
    return np.linalg.eigh(L_H)            # Hermitian → real eigvals, complex eigvecs


def D_raw_mag(ev, U, tau):
    """Magnetic propagator communication distance Δ = 1/|e^{−τL_H}| (symmetric, real)."""
    diag = np.exp(-tau * ev)
    P = (U * diag) @ U.conj().T           # Hermitian complex
    rho = np.abs(P)                       # symmetric, non-negative
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
    return extract_ultrametric_matrix(Zlink, n)


def rho_coph(a, b):
    if np.nanvar(a) < VAR_FLOOR or np.nanvar(b) < VAR_FLOOR:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spearmanr(a, b).statistic)


def ac1(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; x = x - x.mean()
    d = np.dot(x, x)
    return float(np.dot(x[:-1], x[1:]) / d) if d > 0 and x.size > 3 else 0.0


def perm_p(x, fn, n, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    obs = fn(x)
    null = np.array([fn(rng.permutation(x)) for _ in range(n)])
    return obs, (1 + np.sum(null >= obs)) / (n + 1)


def windowed_signed_eigs(X, fs, L_sec, nperseg, band):
    wlen = int(round(L_sec * fs))
    fmax = BRAIN_BANDS[band][1] + 5.0
    out = []
    for s in range(0, X.shape[1] - wlen + 1, wlen):
        freqs, ff = segment_ffts(X[:, s:s + wlen], fs, nperseg, nperseg // 2, fmax)
        if ff.shape[1] < 4:
            continue
        cube = imcoh_signed_cube(ff)
        A = band_signed_average(cube, freqs, *BRAIN_BANDS[band])
        ev, U = signed_eig(A)
        out.append((ev, U, A.shape[0]))
    return out


def coph_by_alpha(eigs_list, alphas):
    out = {a: [] for a in alphas}
    for ev, U, n in eigs_list:
        lam = ev[-1]
        for a in alphas:
            out[a].append(squareform(D_coph(D_raw_mag(ev, U, a / lam), n)))
    return {a: np.asarray(v) for a, v in out.items()}


def ref_coph_by_alpha(ev, U, n, alphas):
    return {a: squareform(D_coph(D_raw_mag(ev, U, a / ev[-1]), n)) for a in alphas}


def zscore_pair(a, b):
    pooled = np.concatenate([a, b])
    mu, sd = np.nanmean(pooled), np.nanstd(pooled)
    return (a - mu) / sd if sd > 0 else a * 0, (b - mu) / sd if sd > 0 else b * 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=["alpha", "beta"])
    ap.add_argument("--L-sec", type=float, default=20.0)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260622)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

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
            At = np.asarray(At, float); Ap = np.asarray(Ap, float)
            n = At.shape[0]
            ev_t, U_t = signed_eig(At)
            ev_p, U_p = signed_eig(Ap)
            perm = rng.permutation(n)
            ev_tp, U_tp = signed_eig(At[np.ix_(perm, perm)])
            cT = ref_coph_by_alpha(ev_t, U_t, n, ALPHA)
            cP = ref_coph_by_alpha(ev_p, U_p, n, ALPHA)
            cTp = ref_coph_by_alpha(ev_tp, U_tp, n, ALPHA)
            cw_pre = coph_by_alpha(windowed_signed_eigs(Xpre, fs, args.L_sec, nperseg, band), ALPHA)
            cw_post = coph_by_alpha(windowed_signed_eigs(Xpost, fs, args.L_sec, nperseg, band), ALPHA)
            for a in ALPHA:
                def contrast(arr, tgt):
                    return np.array([rho_coph(c, tgt) - rho_coph(c, cP[a]) for c in arr])
                s_pre, s_post = contrast(cw_pre[a], cT[a]), contrast(cw_post[a], cT[a])
                sp_post = contrast(cw_post[a], cTp[a])
                if np.sum(np.isfinite(s_pre)) < 5 or np.sum(np.isfinite(s_post)) < 5:
                    continue
                zpre, zpost = zscore_pair(s_pre, s_post)
                shift = float(np.nanmean(zpost) - np.nanmean(zpre))
                ac_post, p_post = perm_p(zpost, ac1, args.n_perm, rng)
                ac_pre, _ = perm_p(zpre, ac1, args.n_perm, rng)
                spp = sp_post[np.isfinite(sp_post)]
                spp = (spp - spp.mean()) / spp.std() if spp.std() > 0 else spp * 0
                ac_plac = ac1(spp)
                gp = zpost[np.isfinite(zpost)]; q95 = np.nanquantile(zpre, 0.95)
                excess_flash = float(np.mean(gp > q95) - np.mean((zpre[np.isfinite(zpre)] + shift) > q95))
                frac_coll = float(np.mean([np.nanvar(c) < VAR_FLOOR for c in cw_post[a]]))
                rows.append(dict(patient=pat, band=band, alpha=float(a),
                                 shift=shift, AC1_post=ac_post, AC1_pre=ac_pre,
                                 AC1_placebo=ac_plac, p_time_post=p_post,
                                 excess_flash=excess_flash,
                                 spread_excess=float(np.nanstd(zpost) - np.nanstd(zpre)),
                                 frac_collapsed=frac_coll,
                                 n_pre=int(np.isfinite(s_pre).sum()),
                                 n_post=int(np.isfinite(s_post).sum())))
            print(f"  {pat} {band}: {len(cw_post[ALPHA[0]])} post-win × {len(ALPHA)} τ "
                  f"(directional, {time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "directional_per_patient.csv", index=False)

    def wtest(d):
        d = np.asarray(d, float); d = d[np.isfinite(d)]
        if d.size < 5 or np.allclose(d, 0):
            return np.nan, 0, len(d)
        try:
            return float(wilcoxon(d, alternative="greater").pvalue), int((d > 0).sum()), len(d)
        except ValueError:
            return np.nan, int((d > 0).sum()), len(d)

    coh = []
    for band in args.bands:
        for a in ALPHA:
            sub = df[(df.band == band) & (np.isclose(df.alpha, a))]
            if sub.empty:
                continue
            p_sh, n_sh, K = wtest(sub["shift"])
            dac = (sub.AC1_post - sub.AC1_pre).to_numpy()
            p_b, n_b, _ = wtest(dac)
            dspec = (sub.AC1_post - sub.AC1_placebo).to_numpy()
            p_sp, n_sp, _ = wtest(dspec)
            p_ef, n_ef, _ = wtest(sub.excess_flash)
            coh.append(dict(band=band, alpha=float(a), K=K,
                            shift_med=float(sub["shift"].median()), shift_npos=f"{n_sh}/{K}", shift_p=p_sh,
                            burst_dAC1_med=float(np.nanmedian(dac)), burst_npos=f"{n_b}/{K}", burst_p=p_b,
                            specific_p=p_sp, flash_npos=f"{n_ef}/{K}", flash_p=p_ef,
                            n_indiv_bursty=int((sub.p_time_post < 0.05).sum()),
                            frac_collapsed=float(sub.frac_collapsed.median())))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "directional_cohort.csv", index=False)

    print("\n=== P1 DIRECTIONAL (magnetic-Laplacian) replay — cohort τ-profile ===")
    print("  shift = directional reinstatement LEVEL (N1-analog); burst = transient STATE;")
    print("  specific = real>placebo; * = burst_p<0.05 AND specific_p<0.05")
    for band in args.bands:
        blk = cohdf[cohdf.band == band].sort_values("alpha")
        if blk.empty:
            continue
        print(f"\n  [{band}]")
        print("   alpha  shift(npos,p)        burst dAC1(npos,p)    specific  flash_p  indiv  collapse")
        for _, r in blk.iterrows():
            hit = "*" if (r.burst_p < 0.05 and r.specific_p < 0.05) else " "
            print(f" {hit}{r.alpha:6.2f}  {r.shift_med:+.3f} {r.shift_npos:>5} p={r.shift_p:.3f}  "
                  f"{r.burst_dAC1_med:+.3f} {r.burst_npos:>5} p={r.burst_p:.3f}  "
                  f"p={r.specific_p:.3f}   {r.flash_p:.3f}   {r.n_indiv_bursty:2d}/{r.K}  {r.frac_collapsed:.2f}")
    print(f"\n[audit_136] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
