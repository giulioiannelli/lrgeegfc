#!/usr/bin/env python
"""
diag_diffusion_distance_trace_desplit.py — EXPLORATION, pass 2 (de-confounded).

Follows pass-1 (diag_diffusion_distance_trace.py) with the two controls that
pass 1 lacked, ported faithfully from audit_121 (the cophenetic tau-sweep):

  1. SPLIT-HALF baseline + SYMMETRIC rho_sym (de-inflates the shared-baseline
     floor AND cancels the first-half/second-half non-equivalence). On the
     diffusion distance D_tau (lrg_eegfc.utils.lrg.diffusion), the object under
     test, instead of audit_121's cophenetic D_coph:

        rho_sym^diff(alpha, arc X) = 0.5 * [
            Spearman(D_X - D_preA,  D_post - D_preB)
          + Spearman(D_X - D_preB,  D_post - D_preA) ]
        X in {task_learn, task_test};  preA/preB = the two rest_pre halves.

  2. TRACE-FREE COLLAPSE PLACEBO (audit_121): manifestly trace-free difference
     vectors; must stay ~0. If it climbs to match rho_sym as alpha coarsens, that
     alpha is in the Fiedler-collapse regime and is NOT interpretable.

        rho_indep^diff(alpha, arc X) = Spearman(D_preA - D_preB, D_X - D_post)

Scale axis: alpha = tau * lambda_max, per-phase tau = alpha / lambda_max_phase
(audit_121 convention; alpha=1 = canonical 1/lambda_max anchor). The clean window
(audit_121) is alpha <~ 10-15; read rho_sym ONLY where the placebo ~ 0.

This is the OBSERVED de-confounded statistic. NO matched-strength null yet (that is
pass 3, mandatory before any cohort claim). Reuses cached full-phase eigenpairs
(load_lrg_result) and cached rest_pre-half eigenpairs (imcoh_lrg_halves). Does NOT
touch rho_sym, the preprint, or any core artifact.

Scope + 5-point preamble + wedge:
  .agents/guides/task-persistence-investigation/2026-07-06_diffusion-distance-trace.md
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
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.lrg.diffusion import diffusion_distance_condensed, effective_modes
from lrg_eegfc.visuals.styles import use_lrg_style

FC_METHOD = "imcoh_abs"
BANDS = list(BRAIN_BANDS_NAMES)
PATIENTS = list(PATIENTS_4PHASE)
ARCS = (("learn", "task_learn"), ("test", "task_test"))
HALF_PHASES = ("rest_pre_A", "rest_pre_B")

# audit_121 alpha grid: anchor alpha=1, finer below, into collapse (alpha=40).
ALPHA = np.unique(np.concatenate([[0.5, 0.7, 0.85, 1.0],
                                  np.geomspace(1.0, 40.0, 26)]))
NA = len(ALPHA)
VAR_FLOOR = 1e-10
HALVES_LRG_CACHE = Path(CACHE_ROOT) / "imcoh_lrg_halves"
OUT_DATA = Path(REPORTS_ROOT) / "diffusion_distance_trace"
OUT_FIG = Path(FIGURES_ROOT) / "diffusion_trace"


def safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman with zero-variance guard (collapsed geometry -> nan)."""
    if np.nanvar(a) < VAR_FLOOR or np.nanvar(b) < VAR_FLOOR:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(spearmanr(a, b).statistic)


def load_eigpair(patient: str, band: str, phase: str):
    """(eigvals ascending, eigvecs) for combinatorial L_hat, or None."""
    if phase in HALF_PHASES:
        p = HALVES_LRG_CACHE / patient / f"{band}_{phase}_lrg_imcoh-abs.npz"
        if not p.exists():
            return None
        z = np.load(p, allow_pickle=True)
        return np.asarray(z["eigenvalues"], float), np.asarray(z["eigenvectors"], float)
    res = load_lrg_result(patient, phase, band, FC_METHOD)
    if res is not None and res.eigenvalues is not None:
        return np.asarray(res.eigenvalues, float), np.asarray(res.eigenvectors, float)
    W = load_fc_matrix(patient, phase, band, FC_METHOD)
    if W is None:
        return None
    W = np.asarray(W, float)
    np.fill_diagonal(W, 0.0)
    lam, U = np.linalg.eigh(np.diag(W.sum(1)) - W)
    return lam, U


def main() -> int:
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    all_phases = list(HALF_PHASES) + ["task_learn", "task_test", "rest_post"]
    np_ = len(PATIENTS)

    rho_sym = {(b, a): np.full((np_, NA), np.nan) for b in BANDS for a, _ in ARCS}
    rho_ind = {(b, a): np.full((np_, NA), np.nan) for b in BANDS for a, _ in ARCS}
    gself = {b: np.full((np_, NA), np.nan) for b in BANDS}
    neff = {b: np.full((np_, NA), np.nan) for b in BANDS}
    alpha_F = {b: np.full(np_, np.nan) for b in BANDS}

    cells = [(b, pi, p) for b in BANDS for pi, p in enumerate(PATIENTS)]
    t_start = time.time()
    skipped = []

    for idx, (band, pi, patient) in enumerate(cells):
        t0 = time.time()
        eig = {ph: load_eigpair(patient, band, ph) for ph in all_phases}
        if any(eig[ph] is None for ph in all_phases):
            skipped.append((patient, band, "missing phase"))
            print(f"[{idx+1}/{len(cells)}] {band:11s} {patient} SKIP missing", flush=True)
            continue
        n_nodes = {ph: eig[ph][0].shape[0] for ph in all_phases}
        if len(set(n_nodes.values())) != 1:
            skipped.append((patient, band, f"N mismatch {n_nodes}"))
            print(f"[{idx+1}/{len(cells)}] {band:11s} {patient} SKIP Nmismatch", flush=True)
            continue

        lam_max = {ph: eig[ph][0].max() for ph in all_phases}
        lam_gap = {ph: eig[ph][0][eig[ph][0] > 1e-10].min() for ph in all_phases}
        # geometric-mean Fiedler collapse scale across the 3 full phases
        full3 = ["rest_pre_A", "task_test", "rest_post"]
        alpha_F[band][pi] = float(np.exp(np.mean([
            np.log(lam_max[ph] / lam_gap[ph]) for ph in full3])))

        # reference (alpha=1) diffusion vectors for self-similarity
        D_ref = {ph: diffusion_distance_condensed(*eig[ph], 1.0 / lam_max[ph])
                 for ph in all_phases}

        for ai, al in enumerate(ALPHA):
            D = {ph: diffusion_distance_condensed(*eig[ph], al / lam_max[ph])
                 for ph in all_phases}
            preA, preB, post = D["rest_pre_A"], D["rest_pre_B"], D["rest_post"]
            for arc, ph_task in ARCS:
                task = D[ph_task]
                # symmetric split-baseline rho_sym (audit_156 form on D_tau)
                r1 = safe_spearman(task - preA, post - preB)
                r2 = safe_spearman(task - preB, post - preA)
                rho_sym[(band, arc)][pi, ai] = np.nanmean([r1, r2])
                # trace-free collapse placebo (audit_121)
                rho_ind[(band, arc)][pi, ai] = safe_spearman(preA - preB, task - post)
            gself[band][pi, ai] = safe_spearman(D["rest_pre_A"], D_ref["rest_pre_A"])
            neff[band][pi, ai] = effective_modes(eig["rest_pre_A"][0], al / lam_max["rest_pre_A"])

        dt = time.time() - t0
        if idx == 0:
            print(f"  [timing] first cell {dt:.2f}s -> ETA ~{dt*len(cells):.0f}s", flush=True)
        el = time.time() - t_start
        print(f"[{idx+1}/{len(cells)}] {band:11s} {patient} {dt:.2f}s "
              f"elapsed {el:.0f}s ETA {el/(idx+1)*(len(cells)-idx-1):.0f}s", flush=True)

    np.savez_compressed(
        OUT_DATA / "rho_diff_desplit_pass2.npz",
        alpha=ALPHA, bands=np.array(BANDS), patients=np.array(PATIENTS),
        **{f"rhosym_{b}_{a}": rho_sym[(b, a)] for b in BANDS for a, _ in ARCS},
        **{f"rhoind_{b}_{a}": rho_ind[(b, a)] for b in BANDS for a, _ in ARCS},
        **{f"gself_{b}": gself[b] for b in BANDS},
        **{f"neff_{b}": neff[b] for b in BANDS},
        **{f"alphaF_{b}": alpha_F[b] for b in BANDS},
    )
    print(f"\nsaved -> {OUT_DATA/'rho_diff_desplit_pass2.npz'}", flush=True)
    if skipped:
        print(f"skipped {len(skipped)}: {skipped}", flush=True)

    _plot(rho_sym, rho_ind, alpha_F)
    _summary(rho_sym, rho_ind, alpha_F)
    print(f"\nwall clock {time.time()-t_start:.0f}s", flush=True)
    return 0


def _clean_alpha_idx(band, alpha_F):
    """Largest alpha index still clean: alpha < median Fiedler-collapse scale."""
    aF = np.nanmedian(alpha_F[band])
    return np.where(ALPHA < aF, np.arange(NA), -1).max(), aF


def _plot(rho_sym, rho_ind, alpha_F):
    import matplotlib.pyplot as plt
    use_lrg_style()
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.4), sharex=True)
    for ax, band in zip(axes.ravel(), BANDS):
        aF = np.nanmedian(alpha_F[band])
        if np.isfinite(aF):
            ax.axvspan(aF, ALPHA.max(), color="0.87", zorder=0)  # collapse regime
        ax.axvline(1.0, color="0.6", lw=0.7, ls=":", zorder=1)   # canonical anchor
        for arc, _ in ARCS:
            med = np.nanmedian(rho_sym[(band, arc)], axis=0)
            ind = np.nanmedian(rho_ind[(band, arc)], axis=0)
            ls = "-" if arc == "test" else "--"
            (ln,) = ax.plot(ALPHA, med, ls=ls, lw=1.9, label=f"{arc} ρ_sym")
            ax.plot(ALPHA, ind, ls=ls, lw=1.0, color=ln.get_color(), alpha=0.55)
        ax.axhline(0.0, color="0.5", lw=0.8)
        ax.set_xscale("log")
        ax.text(0.04, 0.93, BRAIN_BAND_TEX_DICT[band], transform=ax.transAxes,
                fontsize=13, va="top")
    for ax in axes[-1]:
        ax.set_xlabel(r"$\alpha = \tau\,\lambda_{\max}$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\rho^{\mathrm{diff}}_{\mathrm{sym}}$  (thin = placebo)")
    h, l = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=2,
               frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = OUT_FIG / "rho_diff_desplit_pass2.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print(f"saved figure -> {out}", flush=True)


def _summary(rho_sym, rho_ind, alpha_F):
    print("\n" + "=" * 82)
    print("PASS-2 DE-CONFOUNDED SUMMARY  (split-half rho_sym + collapse placebo; NO MS)")
    print("=" * 82)
    print(f"{'band':11s} {'arc':5s} {'sym@a=1':>8s} {'ind@a=1':>8s} "
          f"{'sym@peak':>9s} {'a@peak':>7s} {'ind@peak':>9s} {'aFied':>6s} {'clean?':>7s}")
    for band in BANDS:
        kclean, aF = _clean_alpha_idx(band, alpha_F)
        a1 = int(np.argmin(np.abs(ALPHA - 1.0)))
        for arc, _ in ARCS:
            ms = np.nanmedian(rho_sym[(band, arc)], axis=0)
            mi = np.nanmedian(rho_ind[(band, arc)], axis=0)
            win = ms.copy()
            win[kclean + 1:] = -np.inf          # restrict peak to clean window
            kpk = int(np.nanargmax(win)) if np.isfinite(win).any() else a1
            clean = "yes" if kpk <= kclean else "COLLAPSE"
            print(f"{band:11s} {arc:5s} {ms[a1]:8.3f} {mi[a1]:8.3f} "
                  f"{ms[kpk]:9.3f} {ALPHA[kpk]:7.2f} {mi[kpk]:9.3f} {aF:6.1f} {clean:>7s}")
    print("=" * 82)
    print("Read: sym@a=1 = de-inflated fine-scale trace (compare bands). A clean-")
    print("window peak with sym>>ind and a@peak away from 1 = genuine coarser-scale")
    print("lead. sym@peak~ind@peak, or peak flagged COLLAPSE = artifact. Placebo")
    print("(ind) must stay ~0 in the clean window. Pass 3 = matched-strength.")


if __name__ == "__main__":
    sys.exit(main())
