#!/usr/bin/env python
"""
diag_diffusion_trace_matched_strength.py — EXPLORATION, pass 3 (matched-strength).

The decisive control. Fine-scale (alpha=1, tau=1/lambda_max per phase) symmetric
split-half rho_sym on the DIFFUSION distance, gated against the mandatory
matched-strength null — and computed head-to-head with the COPHENETIC rho_sym on
the *identical* cached surrogate ensemble, to test the pass-2 observation that the
diffusion distance separates bands cleanly (beta/alpha genuine, gammas fake,
delta/theta small) where cophenetic rho_sym is a muddy near-tie.

Why fine-scale only (user, 2026-07-06): a fully connected weighted graph has ~one
characteristic diffusion scale (a single specific-heat peak), so the propagator's
value here is mode-selective band discrimination at that one scale, not multi-scale
resolution (pass 2 established there is no clean coarser-scale trace).

Everything is cached — NO surrogate generation:
  * observed eigenpairs: rest_pre_{A,B} (imcoh_lrg_halves), task_learn/task_test/
    rest_post (load_lrg_result).
  * surrogate eigenpairs: matched_strength_surrogate_lrg R=200 swap=20 seed=20260511
    for all five phases (rest_pre halves included).
Same ensemble drives diffusion + cophenetic (adjacency-from-eigs not needed; both
read the eigenpairs). rho_sym is the symmetric split-baseline form (audit_156);
cohort MS p = (1 + #{median_surr >= median_obs}) / (R+1), one-sided upper tail.

Observed-vs-null only; no new null construction. Does NOT touch core artifacts.
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
from lrg_eegfc.utils.surrogate.matched_strength import (
    surrogate_cache_path, cophenetic_condensed_from_eigs,
)
from lrg_eegfc.visuals.styles import use_lrg_style

FC_METHOD = "imcoh_abs"
BANDS = list(BRAIN_BANDS_NAMES)
PATIENTS = list(PATIENTS_4PHASE)
ARCS = (("learn", "task_learn"), ("test", "task_test"))
PHASES = ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")
HALF_PHASES = ("rest_pre_A", "rest_pre_B")
R, SWAP, SEED = 200, 20, 20260511
VAR_FLOOR = 1e-10
HALVES_LRG_CACHE = Path(CACHE_ROOT) / "imcoh_lrg_halves"
OUT_DATA = Path(REPORTS_ROOT) / "diffusion_distance_trace"
OUT_FIG = Path(FIGURES_ROOT) / "diffusion_trace"


def safe_spearman(a, b):
    if np.nanvar(a) < VAR_FLOOR or np.nanvar(b) < VAR_FLOOR:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spearmanr(a, b).statistic)


def rho_sym(D, ph_task):
    """Symmetric split-baseline rho_sym on a dict of condensed distance vectors."""
    r1 = safe_spearman(D[ph_task] - D["rest_pre_A"], D["rest_post"] - D["rest_pre_B"])
    r2 = safe_spearman(D[ph_task] - D["rest_pre_B"], D["rest_post"] - D["rest_pre_A"])
    return float(np.nanmean([r1, r2]))


def diff_condensed(ev, U):
    """Fine-scale (alpha=1) diffusion distance: tau = 1 / lambda_max."""
    return diffusion_distance_condensed(ev, U, 1.0 / ev[-1])


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
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    npat = len(PATIENTS)

    # per (band, arc): observed[pat] and null[pat, r] for both measures
    obs = {(b, a, m): np.full(npat, np.nan)
           for b in BANDS for a, _ in ARCS for m in ("diff", "coph")}
    null = {(b, a, m): np.full((npat, R), np.nan)
            for b in BANDS for a, _ in ARCS for m in ("diff", "coph")}

    cells = [(b, pi, p) for b in BANDS for pi, p in enumerate(PATIENTS)]
    t0 = time.time()
    skipped = []
    for idx, (band, pi, patient) in enumerate(cells):
        tc = time.time()
        oe = {ph: load_obs_eig(patient, band, ph) for ph in PHASES}
        se = {ph: load_surr_eig(patient, band, ph) for ph in PHASES}
        if any(oe[ph] is None or se[ph] is None for ph in PHASES):
            skipped.append((patient, band))
            print(f"[{idx+1}/{len(cells)}] {band:11s} {patient} SKIP", flush=True)
            continue

        # observed condensed vectors
        Dd = {ph: diff_condensed(*oe[ph]) for ph in PHASES}
        Dc = {ph: cophenetic_condensed_from_eigs(*oe[ph]) for ph in PHASES}
        for arc, ph_task in ARCS:
            obs[(band, arc, "diff")][pi] = rho_sym(Dd, ph_task)
            obs[(band, arc, "coph")][pi] = rho_sym(Dc, ph_task)

        # surrogate ensemble (r-th surrogate paired across phases)
        for r in range(R):
            evs = {ph: (se[ph][0][r], se[ph][1][r]) for ph in PHASES}
            if any(not np.all(np.isfinite(evs[ph][0])) for ph in PHASES):
                continue
            Dd_r = {ph: diff_condensed(*evs[ph]) for ph in PHASES}
            Dc_r = {ph: cophenetic_condensed_from_eigs(*evs[ph]) for ph in PHASES}
            for arc, ph_task in ARCS:
                null[(band, arc, "diff")][pi, r] = rho_sym(Dd_r, ph_task)
                null[(band, arc, "coph")][pi, r] = rho_sym(Dc_r, ph_task)

        dt = time.time() - tc
        if idx == 0:
            print(f"  [timing] first cell {dt:.1f}s -> ETA ~{dt*len(cells)/60:.1f} min",
                  flush=True)
        el = time.time() - t0
        print(f"[{idx+1}/{len(cells)}] {band:11s} {patient} {dt:.1f}s "
              f"elapsed {el/60:.1f}m ETA {el/(idx+1)*(len(cells)-idx-1)/60:.1f}m", flush=True)

    np.savez_compressed(
        OUT_DATA / "rho_diff_matched_strength.npz",
        bands=np.array(BANDS), patients=np.array(PATIENTS), R=R,
        **{f"obs_{b}_{a}_{m}": obs[(b, a, m)]
           for b in BANDS for a, _ in ARCS for m in ("diff", "coph")},
        **{f"null_{b}_{a}_{m}": null[(b, a, m)]
           for b in BANDS for a, _ in ARCS for m in ("diff", "coph")},
    )
    if skipped:
        print(f"skipped: {skipped}", flush=True)

    _summary_and_plot(obs, null)
    print(f"\nwall clock {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


def _cohort_gate(obs_v, null_v):
    """Cohort median obs, MS one-sided upper-tail p on the cohort-median null."""
    med_obs = np.nanmedian(obs_v)
    med_null = np.nanmedian(null_v, axis=0)          # (R,) cohort-median per surrogate
    med_null = med_null[np.isfinite(med_null)]
    if med_null.size == 0 or not np.isfinite(med_obs):
        return med_obs, np.nan, np.nan
    p = (1 + np.sum(med_null >= med_obs)) / (med_null.size + 1)
    return med_obs, float(np.nanmedian(med_null)), float(p)


def _summary_and_plot(obs, null):
    print("\n" + "=" * 92)
    print("PASS-3 MATCHED-STRENGTH  fine-scale rho_sym  (DIFFUSION vs COPHENETIC, same surrogates)")
    print("=" * 92)
    print(f"{'band':11s} {'arc':5s} | {'DIFF obs':>8s} {'null':>7s} {'p_MS':>7s} {'sig':>4s}"
          f" | {'COPH obs':>8s} {'null':>7s} {'p_MS':>7s} {'sig':>4s}")
    res = {}
    for band in BANDS:
        for arc, _ in ARCS:
            row = []
            for m in ("diff", "coph"):
                o, n, p = _cohort_gate(obs[(band, arc, m)], null[(band, arc, m)])
                res[(band, arc, m)] = (o, n, p)
                row += [o, n, p]
            sd = "***" if row[2] < 0.05 else ("." if row[2] < 0.10 else "")
            sc = "***" if row[5] < 0.05 else ("." if row[5] < 0.10 else "")
            print(f"{band:11s} {arc:5s} | {row[0]:8.3f} {row[1]:7.3f} {row[2]:7.3f} {sd:>4s}"
                  f" | {row[3]:8.3f} {row[4]:7.3f} {row[5]:7.3f} {sc:>4s}")
    print("=" * 92)
    print("obs = cohort-median rho_sym; null = median surrogate; p_MS one-sided (>=). ")
    print("*** p<.05, . p<.10. The question: does DIFF give a clean high/low band split")
    print("that COPH blurs? (p_MS is descriptive here — no BH across the 6-band family.)")

    import matplotlib.pyplot as plt
    use_lrg_style()
    x = np.arange(len(BANDS))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, (arc, _) in zip(axes, ARCS):
        for j, (m, off, mk) in enumerate([("diff", -0.16, "o"), ("coph", 0.16, "s")]):
            obv = [res[(b, arc, m)][0] for b in BANDS]
            nuv = [res[(b, arc, m)][1] for b in BANDS]
            pv = [res[(b, arc, m)][2] for b in BANDS]
            # null band per band: 5-95 pct of cohort-median null
            lo = [np.nanpercentile(np.nanmedian(null[(b, arc, m)], 0), 5) for b in BANDS]
            hi = [np.nanpercentile(np.nanmedian(null[(b, arc, m)], 0), 95) for b in BANDS]
            ax.vlines(x + off, lo, hi, color="0.6" if m == "diff" else "0.8", lw=3,
                      alpha=0.8, zorder=1)
            filled = [(m if p < 0.05 else "none") for p in pv]
            for xi, ov, fc in zip(x + off, obv, filled):
                ax.plot(xi, ov, mk, ms=8, mfc=("C0" if m == "diff" else "none"),
                        mec=("C0" if m == "diff" else "C3"), mew=1.6, zorder=3,
                        fillstyle=("full" if fc != "none" else "none"))
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
        ax.text(0.02, 0.96, {"learn": "learn", "test": "test"}[arc],
                transform=ax.transAxes, va="top", fontsize=11)
    axes[0].set_ylabel(r"$\rho_{\mathrm{sym}}$ (fine scale)")
    fig.tight_layout()
    out = OUT_FIG / "rho_diff_matched_strength_bandmap.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print(f"\nsaved figure -> {out}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
