#!/usr/bin/env python
"""
diag_diffusion_trace_ms_tau_sweep.py — EXPLORATION, pass 4 (MS margin vs tau).

Closes the gap in pass 3: matched-strength was only tested at alpha=1 (tau=1/lmax),
the fine-scale/strength-dominated endpoint. This sweeps alpha = tau*lmax and, at
EACH alpha, computes the diffusion rho_sym against BOTH nulls:

  * matched-strength margin:  median_obs - median(strength-surrogate)    <- the test
  * collapse placebo:         rho_indep = Spearman(D_preA-D_preB, D_X-D_post)  (obs)

Hypothesis under test (user, 2026-07-07): the fine-scale diffusion distance is
strength (D_tau^2 ~ strength+adjacency, my scope section 4), but as tau grows the
low-lambda COMMUNITY modes dominate - and matched-strength (which preserves node
strength but scrambles community structure) should NOT reproduce those. So the
obs-vs-strength margin may OPEN at some intermediate alpha even though it is ~0 at
alpha=1. Trustworthy signal = MS margin > 0 AND placebo ~ 0 AND alpha < alpha_Fiedler.

Vectorized over the R=200 cached surrogate eigenpairs (batched Gram + row-wise
rank correlation) so the whole sweep runs in minutes. Everything cached; no
surrogate generation; does NOT touch core artifacts.
Scope: .agents/guides/task-persistence-investigation/2026-07-06_diffusion-distance-trace.md
"""
from __future__ import annotations

import sys
import time
import warnings
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import CACHE_ROOT, REPORTS_ROOT, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.utils.lrg.diffusion import diffusion_distance_condensed
from lrg_eegfc.utils.surrogate.matched_strength import surrogate_cache_path
from lrg_eegfc.visuals.styles import use_lrg_style

FC_METHOD = "imcoh_abs"
BANDS = list(BRAIN_BANDS_NAMES)
PATIENTS = list(PATIENTS_4PHASE)
ARCS = (("learn", "task_learn"), ("test", "task_test"))
PHASES = ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")
HALF_PHASES = ("rest_pre_A", "rest_pre_B")
R, SWAP, SEED = 200, 20, 20260511
ALPHA = np.geomspace(0.4, 6.0, 16)
NA = len(ALPHA)
VAR_FLOOR = 1e-10
HALVES_LRG_CACHE = Path(CACHE_ROOT) / "imcoh_lrg_halves"
OUT_DATA = Path(REPORTS_ROOT) / "diffusion_distance_trace"
OUT_FIG = Path(FIGURES_ROOT) / "diffusion_trace"
IU = None  # upper-tri index, set per-N


def _tri(n):
    return np.triu_indices(n, k=1)


def _ranks(A):
    """Row-wise ranks (axis=1); ties broken by order (fine for continuous D)."""
    order = A.argsort(axis=1)
    rk = np.empty_like(A, dtype=float)
    cols = np.arange(A.shape[1], dtype=float)
    np.put_along_axis(rk, order, np.broadcast_to(cols, A.shape), axis=1)
    return rk


def _rowcorr(A, B):
    """Row-wise Pearson (A,B ranked -> Spearman). (R,M),(R,M)->(R,)."""
    Az = A - A.mean(1, keepdims=True)
    Bz = B - B.mean(1, keepdims=True)
    num = (Az * Bz).sum(1)
    den = np.sqrt((Az ** 2).sum(1) * (Bz ** 2).sum(1))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / den, np.nan)


def _diff_condensed_batch(eigvals, eigvecs, alpha, iu):
    """Vectorized diffusion distance for an ensemble. eigvals (R,N), eigvecs
    (R,N,N) -> condensed (R,M) at tau_r = alpha / lmax_r."""
    lam = eigvals[:, 1:]                                   # (R, N-1)
    tau = alpha / eigvals[:, -1]                           # (R,)
    w = np.exp(-2.0 * tau[:, None] * lam)                  # (R, N-1)
    psi = eigvecs[:, :, 1:] * np.sqrt(w)[:, None, :]       # (R, N, N-1)
    gram = np.matmul(psi, psi.transpose(0, 2, 1))          # (R, N, N)
    g = np.diagonal(gram, axis1=1, axis2=2)                # (R, N)
    d2 = g[:, :, None] + g[:, None, :] - 2.0 * gram
    d2 = np.clip(d2, 0.0, None)
    return np.sqrt(d2[:, iu[0], iu[1]])                    # (R, M)


def safe_spearman(a, b):
    if np.nanvar(a) < VAR_FLOOR or np.nanvar(b) < VAR_FLOOR:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spearmanr(a, b).statistic)


def rho_sym_null(Dn, ph_task):
    """Vectorized symmetric rho_sym for the surrogate ensemble.
    Dn: dict phase -> (R,M). Returns (R,)."""
    xA = _ranks(Dn[ph_task] - Dn["rest_pre_A"]); yB = _ranks(Dn["rest_post"] - Dn["rest_pre_B"])
    xB = _ranks(Dn[ph_task] - Dn["rest_pre_B"]); yA = _ranks(Dn["rest_post"] - Dn["rest_pre_A"])
    return 0.5 * (_rowcorr(xA, yB) + _rowcorr(xB, yA))


def load_obs_eig(patient, band, phase):
    if phase in HALF_PHASES:
        p = HALVES_LRG_CACHE / patient / f"{band}_{phase}_lrg_imcoh-abs.npz"
        if not p.exists():
            return None
        z = np.load(p, allow_pickle=True)
        return np.asarray(z["eigenvalues"], float), np.asarray(z["eigenvectors"], float)
    res = load_lrg_result(patient, phase, band, FC_METHOD)
    if res is None or res.eigenvalues is None:
        return None
    return np.asarray(res.eigenvalues, float), np.asarray(res.eigenvectors, float)


def load_surr_eig(patient, band, phase):
    p = surrogate_cache_path(patient, band, phase, R, SWAP, SEED, FC_METHOD)
    if not p.exists():
        return None
    with np.load(p) as z:
        return np.asarray(z["eigvals"], float), np.asarray(z["eigvecs"], float)


def main() -> int:
    OUT_DATA.mkdir(parents=True, exist_ok=True); OUT_FIG.mkdir(parents=True, exist_ok=True)
    npat = len(PATIENTS)
    obs = {(b, a): np.full((npat, NA), np.nan) for b in BANDS for a, _ in ARCS}
    nul = {(b, a): np.full((npat, NA, R), np.nan) for b in BANDS for a, _ in ARCS}
    pla = {(b, a): np.full((npat, NA), np.nan) for b in BANDS for a, _ in ARCS}
    aF = {b: np.full(npat, np.nan) for b in BANDS}

    cells = [(b, pi, p) for b in BANDS for pi, p in enumerate(PATIENTS)]
    t0 = time.time(); skipped = []
    for idx, (band, pi, patient) in enumerate(cells):
        tc = time.time()
        oe = {ph: load_obs_eig(patient, band, ph) for ph in PHASES}
        se = {ph: load_surr_eig(patient, band, ph) for ph in PHASES}
        if any(oe[ph] is None or se[ph] is None for ph in PHASES):
            skipped.append((patient, band)); print(f"[{idx+1}/{len(cells)}] SKIP {band} {patient}", flush=True); continue
        N = oe["rest_pre_A"][0].shape[0]; iu = _tri(N)
        lgap = {ph: oe[ph][0][oe[ph][0] > 1e-10].min() for ph in PHASES}
        lmax = {ph: oe[ph][0][-1] for ph in PHASES}
        aF[band][pi] = float(np.exp(np.mean([np.log(lmax[ph] / lgap[ph])
                                             for ph in ("rest_pre_A", "task_test", "rest_post")])))
        for ai, al in enumerate(ALPHA):
            Do = {ph: diffusion_distance_condensed(*oe[ph], al / lmax[ph]) for ph in PHASES}
            Dn = {ph: _diff_condensed_batch(se[ph][0], se[ph][1], al, iu) for ph in PHASES}
            for arc, pt in ARCS:
                r1 = safe_spearman(Do[pt] - Do["rest_pre_A"], Do["rest_post"] - Do["rest_pre_B"])
                r2 = safe_spearman(Do[pt] - Do["rest_pre_B"], Do["rest_post"] - Do["rest_pre_A"])
                obs[(band, arc)][pi, ai] = np.nanmean([r1, r2])
                nul[(band, arc)][pi, ai] = rho_sym_null(Dn, pt)
                pla[(band, arc)][pi, ai] = safe_spearman(
                    Do["rest_pre_A"] - Do["rest_pre_B"], Do[pt] - Do["rest_post"])
        dt = time.time() - tc
        if idx == 0:
            print(f"  [timing] first cell {dt:.1f}s -> ETA ~{dt*len(cells)/60:.1f} min", flush=True)
        el = time.time() - t0
        print(f"[{idx+1}/{len(cells)}] {band:11s} {patient} {dt:.1f}s "
              f"elapsed {el/60:.1f}m ETA {el/(idx+1)*(len(cells)-idx-1)/60:.1f}m", flush=True)

    np.savez_compressed(OUT_DATA / "rho_diff_ms_tau_sweep.npz", alpha=ALPHA,
                        bands=np.array(BANDS), patients=np.array(PATIENTS),
                        **{f"obs_{b}_{a}": obs[(b, a)] for b in BANDS for a, _ in ARCS},
                        **{f"nul_{b}_{a}": nul[(b, a)] for b in BANDS for a, _ in ARCS},
                        **{f"pla_{b}_{a}": pla[(b, a)] for b in BANDS for a, _ in ARCS},
                        **{f"aF_{b}": aF[b] for b in BANDS})
    if skipped:
        print(f"skipped: {skipped}", flush=True)
    _summary_plot(obs, nul, pla, aF)
    print(f"\nwall clock {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


def _summary_plot(obs, nul, pla, aF):
    print("\n" + "=" * 96)
    print("PASS-4  MATCHED-STRENGTH MARGIN vs TAU  (does obs-null open at intermediate alpha?)")
    print("=" * 96)
    print(f"{'band':11s} {'arc':5s} {'aFied':>6s} | {'margin@a=1':>10s} "
          f"{'MAXmargin':>9s} {'@alpha':>7s} {'ratio':>6s} {'placebo':>8s} {'p_MS':>6s} {'clean?':>7s}")
    for band in BANDS:
        aFm = np.nanmedian(aF[band])
        a1 = int(np.argmin(np.abs(ALPHA - 1.0)))
        clean = ALPHA < aFm
        for arc, _ in ARCS:
            mo = np.nanmedian(obs[(band, arc)], axis=0)                 # (NA,)
            mn = np.nanmedian(np.nanmedian(nul[(band, arc)], axis=0), axis=-1)  # (NA,)
            margin = mo - mn
            mp = np.nanmedian(pla[(band, arc)], axis=0)
            mwin = np.where(clean, margin, -np.inf)
            kmx = int(np.nanargmax(mwin)) if np.isfinite(mwin).any() else a1
            # MS p at the max-margin alpha (cohort-median null distribution)
            nd = np.nanmedian(nul[(band, arc)][:, kmx, :], axis=0)
            nd = nd[np.isfinite(nd)]
            p = (1 + np.sum(nd >= mo[kmx])) / (nd.size + 1) if nd.size else np.nan
            ratio = mo[kmx] / mn[kmx] if mn[kmx] > 1e-6 else np.nan
            print(f"{band:11s} {arc:5s} {aFm:6.1f} | {margin[a1]:10.3f} "
                  f"{margin[kmx]:9.3f} {ALPHA[kmx]:7.2f} {ratio:6.2f} {mp[kmx]:8.3f} "
                  f"{p:6.3f} {'yes' if clean[kmx] else 'COLL':>7s}")
    print("=" * 96)
    print("margin = median(obs) - median(strength-surrogate). Signal we want: a clean-window")
    print("(alpha<aFied) alpha where margin>>0, ratio>>1, placebo~0, p_MS<.05. margin~0 at all")
    print("clean alpha = diffusion is strength at every scale (multiscale hope dead).")

    import matplotlib.pyplot as plt
    use_lrg_style()
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.4), sharex=True)
    for ax, band in zip(axes.ravel(), BANDS):
        aFm = np.nanmedian(aF[band])
        if np.isfinite(aFm):
            ax.axvspan(aFm, ALPHA.max(), color="0.87", zorder=0)
        ax.axvline(1.0, color="0.6", lw=0.7, ls=":")
        ax.axhline(0.0, color="0.5", lw=0.8)
        for arc, _ in ARCS:
            mo = np.nanmedian(obs[(band, arc)], axis=0)
            mn = np.nanmedian(np.nanmedian(nul[(band, arc)], axis=0), axis=-1)
            mp = np.nanmedian(pla[(band, arc)], axis=0)
            ls = "-" if arc == "test" else "--"
            (ln,) = ax.plot(ALPHA, mo - mn, ls=ls, lw=2.0, label=f"{arc} MS margin")
            ax.plot(ALPHA, mp, ls=ls, lw=1.0, color=ln.get_color(), alpha=0.5)
        ax.set_xscale("log")
        ax.text(0.04, 0.93, BRAIN_BAND_TEX_DICT[band], transform=ax.transAxes,
                fontsize=13, va="top")
    for ax in axes[-1]:
        ax.set_xlabel(r"$\alpha=\tau\lambda_{\max}$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"obs $-$ strength-null  (thin=placebo)")
    h, l = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = OUT_FIG / "rho_diff_ms_tau_sweep.pdf"
    fig.savefig(out, transparent=True); plt.close(fig)
    print(f"\nsaved figure -> {out}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
