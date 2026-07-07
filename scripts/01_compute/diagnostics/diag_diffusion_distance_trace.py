#!/usr/bin/env python
"""
diag_diffusion_distance_trace.py — EXPLORATION, pass 1 (shape reconnaissance).

rho^diff(tau): the tau-resolved diffusion-distance version of rho_split.

    D_tau(i,j)^2 = sum_{a>=2} exp(-2*tau*lambda_a) * (u_a(i) - u_a(j))^2

built from the eigenpairs of the combinatorial Laplacian L_hat = D_hat - W of the
imcoh_abs FC matrix (the same operator the pipeline / Grassmann diagonalise). The
cross-phase functional is the drop-in rho_split generalization:

    dtask(tau) = D_tau[X]        - D_tau[rest_pre]     (X in {task_learn, task_test})
    drest(tau) = D_tau[rest_post] - D_tau[rest_pre]
    rho^diff_{p,b,X}(tau) = Spearman( triu dtask, triu drest )

PASS 1 SCOPE — this is shape reconnaissance, NOT a claim:
  * shared full rest_pre baseline (cheap; reuses cached eigenpairs). This inflates
    absolute rho by a common additive term, so we read CONTRASTS ONLY:
    per-band heterogeneity, learn-vs-test separation, and the location of any
    intermediate-tau hump. Absolute rho is not interpreted.
  * collapse guard: n_eff(tau) is recorded; the region where median_p n_eff < 2 is
    the Fiedler-collapse tail (audit_121) and is shaded "do not interpret".
  * NO matched-strength null and NO split-baseline placebo here — those are pass 2,
    mandatory before any cohort claim.

Does NOT touch the core paper: no rho_sym, no preprint, no existing results.
Outputs land under data/reports/diffusion_distance_trace/ and
data/outputs/figures/diffusion_trace/ only.

Full derivation, 5-point critical preamble, wedge argument and guardrails:
  .agents/guides/task-persistence-investigation/2026-07-06_diffusion-distance-trace.md
"""
from __future__ import annotations

import sys
import time
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.visuals.styles import use_lrg_style

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
FC_METHOD = "imcoh_abs"
PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
ARCS = (("learn", "task_learn"), ("test", "task_test"))
BANDS = list(BRAIN_BANDS_NAMES)
PATIENTS = list(PATIENTS_4PHASE)

# Dimensionless diffusion-scale ladder s = tau * lambda2(rest_pre), common to all
# patients by construction (tau_grid = S_GRID / lambda2_ref per cell). Spans the
# fine-scale/strength regime (small s) through the Fiedler-collapse tail (large s).
S_GRID = np.geomspace(0.02, 6.0, 32)
K = len(S_GRID)

OUT_DATA = Path(REPORTS_ROOT) / "diffusion_distance_trace"
OUT_FIG = Path(FIGURES_ROOT) / "diffusion_trace"


# ---------------------------------------------------------------------------
# Diffusion-distance primitive (LOCAL for pass 1; promote to
# src/lrg_eegfc/utils/lrg/diffusion.py on its pass-2 second caller per the
# library-first rule — general graph concept, general name).
# ---------------------------------------------------------------------------
def diffusion_distance_condensed(eigvals: np.ndarray, eigvecs: np.ndarray,
                                 tau: float) -> np.ndarray:
    """Upper-triangle vector of D_tau(i,j) = sqrt(sum_{a>=2} e^{-2 tau lam_a} (u_a(i)-u_a(j))^2).

    eigvals : (N,) ascending; eigvecs : (N,N) columns = unit-norm modes, col 0 trivial.
    Drops the trivial constant mode (a=0), which contributes 0 to every distance.
    """
    lam = eigvals[1:]                       # (N-1,) drop trivial mode
    U = eigvecs[:, 1:]                       # (N, N-1)
    w = np.exp(-2.0 * tau * lam)             # (N-1,) diffusion weights
    psi = U * np.sqrt(w)[None, :]            # (N, N-1) diffusion coordinates
    gram = psi @ psi.T                       # (N, N)
    g = np.diag(gram)
    d2 = g[:, None] + g[None, :] - 2.0 * gram
    iu = np.triu_indices(d2.shape[0], k=1)
    return np.sqrt(np.clip(d2[iu], 0.0, None))


def n_eff(eigvals: np.ndarray, tau: float) -> float:
    """Participation ratio of the diffusion weights = effective number of active
    modes. Runs from N-1 (tau->0, all modes) down to ~1 (Fiedler-only collapse)."""
    w = np.exp(-2.0 * tau * eigvals[1:])
    return float(w.sum() ** 2 / np.sum(w ** 2))


# ---------------------------------------------------------------------------
# Eigenpair loading (cached; eigh fallback for legacy caches)
# ---------------------------------------------------------------------------
def load_eigs(patient: str, phase: str, band: str):
    """Return (eigvals ascending, eigvecs) for the combinatorial L_hat, or None."""
    res = load_lrg_result(patient, phase, band, FC_METHOD)
    if res is not None and res.eigenvalues is not None and res.eigenvectors is not None:
        return np.asarray(res.eigenvalues), np.asarray(res.eigenvectors)
    W = load_fc_matrix(patient, phase, band, FC_METHOD)
    if W is None:
        return None
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    lap = np.diag(W.sum(axis=1)) - W
    lam, U = np.linalg.eigh(lap)             # ascending, orthonormal columns
    return lam, U


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
def main() -> int:
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)

    n_pat = len(PATIENTS)
    # rho[(band, arc)] -> (n_pat, K) ; nan where a cell is unavailable
    rho = {(b, a): np.full((n_pat, K), np.nan) for b in BANDS for a, _ in ARCS}
    neff = {b: np.full((n_pat, K), np.nan) for b in BANDS}

    cells = [(b, pi, p) for b in BANDS for pi, p in enumerate(PATIENTS)]
    n_cells = len(cells)
    t_start = time.time()
    skipped = []

    for idx, (band, pi, patient) in enumerate(cells):
        t0 = time.time()

        eig = {ph: load_eigs(patient, ph, band) for ph in PHASES}
        if any(eig[ph] is None for ph in PHASES):
            skipped.append((patient, band, "missing phase"))
            print(f"[{idx+1}/{n_cells}] {band:11s} {patient}  SKIP (missing phase)",
                  flush=True)
            continue

        # N must match across the 4 phases to compare the same pair index.
        n_nodes = {ph: eig[ph][0].shape[0] for ph in PHASES}
        if len(set(n_nodes.values())) != 1:
            skipped.append((patient, band, f"N mismatch {n_nodes}"))
            print(f"[{idx+1}/{n_cells}] {band:11s} {patient}  SKIP (N mismatch {n_nodes})",
                  flush=True)
            continue

        lam2_ref = eig["rest_pre"][0][1]     # Fiedler value of rest_pre
        if not np.isfinite(lam2_ref) or lam2_ref <= 0:
            skipped.append((patient, band, "bad lambda2"))
            print(f"[{idx+1}/{n_cells}] {band:11s} {patient}  SKIP (bad lambda2)",
                  flush=True)
            continue
        tau_grid = S_GRID / lam2_ref

        # Diffusion distances per phase per tau -> (K, M)
        D = {ph: np.stack([diffusion_distance_condensed(*eig[ph], t) for t in tau_grid])
             for ph in PHASES}
        neff[band][pi] = [n_eff(eig["rest_pre"][0], t) for t in tau_grid]

        d_pre, d_post = D["rest_pre"], D["rest_post"]
        for arc, ph_task in ARCS:
            d_task = D[ph_task]
            for k in range(K):
                dtask = d_task[k] - d_pre[k]
                drest = d_post[k] - d_pre[k]
                r = spearmanr(dtask, drest).statistic
                rho[(band, arc)][pi, k] = r

        dt = time.time() - t0
        if idx == 0:
            eta = dt * n_cells
            print(f"  [timing] first cell {dt:.2f}s  ->  full-sweep ETA ~{eta:.0f}s "
                  f"({n_cells} cells)", flush=True)
        elapsed = time.time() - t_start
        eta = elapsed / (idx + 1) * (n_cells - idx - 1)
        print(f"[{idx+1}/{n_cells}] {band:11s} {patient}  {dt:.2f}s  "
              f"elapsed {elapsed:.0f}s  ETA {eta:.0f}s", flush=True)

    # ---- save arrays ----
    np.savez_compressed(
        OUT_DATA / "rho_diff_pass1.npz",
        s_grid=S_GRID,
        bands=np.array(BANDS),
        patients=np.array(PATIENTS),
        **{f"rho_{b}_{a}": rho[(b, a)] for b in BANDS for a, _ in ARCS},
        **{f"neff_{b}": neff[b] for b in BANDS},
    )
    print(f"\nsaved arrays -> {OUT_DATA/'rho_diff_pass1.npz'}", flush=True)
    if skipped:
        print(f"skipped {len(skipped)} cells: {skipped}", flush=True)

    _plot(rho, neff)
    _print_contrast_summary(rho, neff)
    print(f"\nwall clock {time.time()-t_start:.0f}s", flush=True)
    return 0


def _collapse_mask(neff_band: np.ndarray) -> np.ndarray:
    """True where median_p n_eff >= 2 (interpretable, not collapsed)."""
    med = np.nanmedian(neff_band, axis=0)
    return med >= 2.0


def _plot(rho, neff):
    import matplotlib.pyplot as plt
    use_lrg_style()

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.2), sharex=True)
    for ax, band in zip(axes.ravel(), BANDS):
        ok = _collapse_mask(neff[band])
        # shade the collapse tail (median n_eff < 2)
        if (~ok).any():
            s_collapse = S_GRID[~ok]
            ax.axvspan(s_collapse.min(), S_GRID.max(), color="0.85", zorder=0)
        for arc, _ in ARCS:
            med = np.nanmedian(rho[(band, arc)], axis=0)
            q1 = np.nanpercentile(rho[(band, arc)], 25, axis=0)
            q3 = np.nanpercentile(rho[(band, arc)], 75, axis=0)
            ls = "-" if arc == "test" else "--"
            (line,) = ax.plot(S_GRID, med, ls=ls, lw=1.8, label=arc)
            ax.fill_between(S_GRID, q1, q3, alpha=0.15, color=line.get_color())
        ax.axhline(0.0, color="0.5", lw=0.8, zorder=1)
        ax.set_xscale("log")
        ax.text(0.05, 0.90, BRAIN_BAND_TEX_DICT[band], transform=ax.transAxes,
                fontsize=13, va="top")
    for ax in axes[-1]:
        ax.set_xlabel(r"$s = \tau\,\lambda_2$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\rho^{\mathrm{diff}}(\tau)$")

    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = OUT_FIG / "rho_diff_tau_profiles_pass1.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print(f"saved figure -> {out}", flush=True)


def _print_contrast_summary(rho, neff):
    """Text read-out of the contrasts we are actually allowed to interpret in
    pass 1: per-band spread of the peak (non-collapsed) rho, and the learn-vs-test
    gap. Absolute rho is inflated (shared baseline) — reported for shape only."""
    print("\n" + "=" * 74)
    print("PASS-1 CONTRAST SUMMARY  (shared-baseline; read contrasts, not levels)")
    print("=" * 74)
    print(f"{'band':11s} {'arc':5s} {'peak_rho':>9s} {'s@peak':>7s} "
          f"{'rho@s=0.02':>11s} {'n_ok':>5s}")
    for band in BANDS:
        ok = _collapse_mask(neff[band])
        for arc, _ in ARCS:
            med = np.nanmedian(rho[(band, arc)], axis=0)
            med_ok = np.where(ok, med, -np.inf)
            kpk = int(np.nanargmax(med_ok)) if np.isfinite(med_ok).any() else 0
            n_ok = int(np.sum(np.isfinite(rho[(band, arc)][:, kpk])))
            print(f"{band:11s} {arc:5s} {med[kpk]:9.3f} {S_GRID[kpk]:7.3f} "
                  f"{med[0]:11.3f} {n_ok:5d}")
    print("=" * 74)
    print("Reading: a band-specific INTERMEDIATE-s peak (s@peak away from the")
    print("small-s edge and outside the shaded collapse tail) is the candidate")
    print("higher-order signal. A peak pinned at s=0.02 = fine-scale/strength")
    print("regime = no gain over raw FC (the audit_121 falsifier). Learn-vs-test")
    print("separation at intermediate s = the second core. Pass 2 gates with")
    print("split-baseline rho_sym + collapse placebo + matched-strength.")


if __name__ == "__main__":
    sys.exit(main())
