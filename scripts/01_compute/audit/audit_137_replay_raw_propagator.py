#!/usr/bin/env python3
"""audit_137 — N4 v3 / P5: RAW-PROPAGATOR transient (skip the cophenetic quantization).

Scope: .agents/guides/task-persistence-investigation/2026-06-22_replay-states-v3-scale-target-pair.md

Every prior replay test went through the LRG COPHENETIC/dendrogram — a heavily QUANTIZED
ultrametric (per-pair cophenetic is ~98% tied to N−1 merge heights; perpair_cophenetic
memory). A transient reinstatement could be SMOOTHED AWAY by that quantization yet survive
in the RAW, continuous diffusion propagator. P5 repeats the τ-resolved burstiness battery
(audit_133) but on the RAW propagator communication distance Δ(τ)=1/ρ̂(τ) — NO dendrogram.

Per patient, band, α (τ = α/λmax_window):
  Δ_w(α) = triu(1/ρ̂_w(α)), ρ̂=e^{−τL}/Z, of windowed |ImCoh| (20 s); refs Δ_T(α), Δ_{R-}(α).
  s_w(α) = ρ(Δ_w(α), Δ_T(α)) − ρ(Δ_w(α), Δ_{R-}(α))   (Spearman of raw triu, z-scored)
  shift(α) = mean_post − mean_pre  (raw-leg sustained level);
  AC1/excess/spread vs TIME-SHUFFLE = transient STATES; node-perm-target PLACEBO control.

Verdict: a scale τ* clears only if burstiness > time-shuffle AND target-specific (real >
placebo). This is the raw-FC leg (audit_121 T_d_raw / rho_split_raw) at the WINDOW level,
which the cophenetic cache never tested. Raw FC is a COMPARISON baseline not a result
([[feedback_results_only_in_laplacian_framework]]) — but for a NEGATIVE it rules out the
"quantization hid the transient" alternative; a POSITIVE would demand a cophenetic re-test.

Output: data/audit/replay_states/raw_propagator_{per_patient,cohort}.csv
"""
from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import segment_ffts, imcoh_abs_cube, band_abs_average
from lrg_eegfc.workflow.fc import load_fc_matrix

OUT = ROOT / "data" / "audit" / "replay_states"
OUT.mkdir(parents=True, exist_ok=True)
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
# fine → past-Fiedler (collapse). α=1 = locked finest scale.
ALPHA = np.unique(np.r_[0.5, 1.0, np.geomspace(1.0, 40.0, 11)]).round(4)
VAR_FLOOR = 1e-10
_TRIU_CACHE: dict[int, tuple] = {}


def _triu(n):
    if n not in _TRIU_CACHE:
        _TRIU_CACHE[n] = np.triu_indices(n, 1)
    return _TRIU_CACHE[n]


def D_raw(ev, U, tau):                       # audit_121 verbatim
    diag = np.exp(-tau * ev)
    Z = float(diag.sum())
    rho = (U * diag) @ U.T / Z
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 1.0 / rho
    D = np.maximum(D, D.T)
    np.fill_diagonal(D, 0.0)
    fin = np.isfinite(D)
    if not fin.all():
        D = np.where(fin, D, np.nanmax(D[fin]) if fin.any() else 1e12)
    return D


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


def eig_of_W(W):
    W = 0.5 * (W + W.T); np.fill_diagonal(W, 0.0)
    L = np.diag(W.sum(1)) - W
    return np.linalg.eigh(L)


def windowed_eigs(X, fs, L_sec, nperseg, band):
    wlen = int(round(L_sec * fs))
    fmax = BRAIN_BANDS[band][1] + 5.0
    out = []
    for s in range(0, X.shape[1] - wlen + 1, wlen):
        freqs, ff = segment_ffts(X[:, s:s + wlen], fs, nperseg, nperseg // 2, fmax)
        if ff.shape[1] < 4:
            continue
        W = band_abs_average(imcoh_abs_cube(ff), freqs, *BRAIN_BANDS[band])
        ev, U = eig_of_W(W)
        out.append((ev, U, W.shape[0]))
    return out


def raw_by_alpha(eigs_list, alphas):
    """dict α -> array(n_win, n_pairs) RAW propagator distance triu at τ=α/λmax_window."""
    out = {a: [] for a in alphas}
    for ev, U, n in eigs_list:
        lam = ev[-1]; iu = _triu(n)
        for a in alphas:
            out[a].append(D_raw(ev, U, a / lam)[iu])
    return {a: np.asarray(v) for a, v in out.items()}


def ref_raw_by_alpha(ev, U, n, alphas):
    iu = _triu(n)
    return {a: D_raw(ev, U, a / ev[-1])[iu] for a in alphas}


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
            Wt = load_fc_matrix(pat, "task_test", band, fc_method="imcoh_abs")
            Wp = load_fc_matrix(pat, "rest_pre", band, fc_method="imcoh_abs")
            if Wt is None or Wp is None:
                print(f"  [skip] {pat} {band}: refs missing"); continue
            n = Wt.shape[0]
            ev_t, U_t = eig_of_W(np.asarray(Wt, float))
            ev_p, U_p = eig_of_W(np.asarray(Wp, float))
            perm = rng.permutation(n)                     # placebo: scramble task anatomy
            ev_tp, U_tp = eig_of_W(np.asarray(Wt, float)[np.ix_(perm, perm)])
            cT = ref_raw_by_alpha(ev_t, U_t, n, ALPHA)
            cP = ref_raw_by_alpha(ev_p, U_p, n, ALPHA)
            cTp = ref_raw_by_alpha(ev_tp, U_tp, n, ALPHA)
            cw_pre = raw_by_alpha(windowed_eigs(Xpre, fs, args.L_sec, nperseg, band), ALPHA)
            cw_post = raw_by_alpha(windowed_eigs(Xpost, fs, args.L_sec, nperseg, band), ALPHA)
            for a in ALPHA:
                def contrast(arr, tgt):
                    return np.array([rho_coph(c, tgt) - rho_coph(c, cP[a]) for c in arr])
                s_pre, s_post = contrast(cw_pre[a], cT[a]), contrast(cw_post[a], cT[a])
                sp_post = contrast(cw_post[a], cTp[a])     # placebo target
                if np.sum(np.isfinite(s_pre)) < 5 or np.sum(np.isfinite(s_post)) < 5:
                    continue
                zpre, zpost = zscore_pair(s_pre, s_post)
                shift = float(np.nanmean(zpost) - np.nanmean(zpre))
                ac_post, p_post = perm_p(zpost, ac1, args.n_perm, rng)
                ac_pre, _ = perm_p(zpre, ac1, args.n_perm, rng)
                # placebo burstiness (z-scored on its own pool)
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
                  f"({time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "raw_propagator_per_patient.csv", index=False)

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
            p_sh, n_sh, _ = wtest(sub["shift"])
            dac = (sub.AC1_post - sub.AC1_pre).to_numpy()
            p_b, n_b, K = wtest(dac)
            dspec = (sub.AC1_post - sub.AC1_placebo).to_numpy()
            p_sp, n_sp, _ = wtest(dspec)
            p_ef, n_ef, _ = wtest(sub.excess_flash)
            coh.append(dict(band=band, alpha=float(a), K=K,
                            shift_med=float(sub["shift"].median()), shift_npos=f"{n_sh}/{K}", shift_p=p_sh,
                            burst_dAC1_med=float(np.nanmedian(dac)), burst_npos=f"{n_b}/{K}", burst_p=p_b,
                            specific_med=float(np.nanmedian(dspec)), specific_npos=f"{n_sp}/{K}", specific_p=p_sp,
                            flash_npos=f"{n_ef}/{K}", flash_p=p_ef,
                            n_indiv_bursty=int((sub.p_time_post < 0.05).sum()),
                            frac_collapsed=float(sub.frac_collapsed.median())))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "raw_propagator_cohort.csv", index=False)

    print("\n=== P5 RAW-PROPAGATOR transient — cohort τ-profile ===")
    print("  (shift = N1 level @scale; burst = AC1 post>pre vs shuffle; "
          "specific = real>placebo; * burst_p<0.05 AND specific_p<0.05)")
    for band in args.bands:
        blk = cohdf[cohdf.band == band].sort_values("alpha")
        if blk.empty:
            continue
        print(f"\n  [{band}]  α=1 is the locked finest scale (anchor = sustained 10/10)")
        print("   alpha  shift(npos,p)        burst dAC1(npos,p)    specific(p)   flash_p  indiv  collapse")
        for _, r in blk.iterrows():
            hit = "*" if (r.burst_p < 0.05 and r.specific_p < 0.05) else " "
            print(f" {hit}{r.alpha:6.2f}  {r.shift_med:+.3f} {r.shift_npos:>5} p={r.shift_p:.3f}  "
                  f"{r.burst_dAC1_med:+.3f} {r.burst_npos:>5} p={r.burst_p:.3f}  "
                  f"p={r.specific_p:.3f}   {r.flash_p:.3f}   {r.n_indiv_bursty:2d}/{r.K}  {r.frac_collapsed:.2f}")
    print(f"\n[audit_137] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
