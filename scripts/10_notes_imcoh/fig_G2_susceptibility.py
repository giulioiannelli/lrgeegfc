"""fig_G2 — entropy/susceptibility twinx plot, MSC vs ImCoh.

Reuses the recipe from `scripts/07_figures/gen_mslcd_figures.py`
(`_make_entropy_susceptibility_fig`): τ-grid spans [1/λ_max·0.3, 1/λ_2·50],
Ŝ = S/ln N, C̃ = -dŜ/d log10 τ, twin y-axis with (1-Ŝ) on the left.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from _shared import METHOD_COLORS, SECTION2_ROOT, apply_pub_style, save_fig
from lrg_eegfc.workflow.fc import load_fc_matrix


def entropy_at_tau(eigenvalues, tau):
    log_boltz = -tau * eigenvalues
    log_Z = np.logaddexp.reduce(log_boltz)
    log_p = log_boltz - log_Z
    p = np.exp(log_p)
    return -np.sum(p * np.where(p > 1e-300, log_p, 0))


def laplacian_eigvals(A):
    A = np.abs(A).copy()
    np.fill_diagonal(A, 0.0)
    L = np.diag(A.sum(axis=1)) - A
    ev = np.sort(np.linalg.eigvalsh(L))
    return np.maximum(ev, 0.0)


def compute_curves(eig_vals):
    N_b = len(eig_vals)
    lnN_b = np.log(N_b)
    tau_min_b = 1.0 / eig_vals[-1]
    tau_max_b = 1.0 / eig_vals[1]
    tau_grid = np.logspace(np.log10(tau_min_b * 0.3),
                           np.log10(tau_max_b * 50), 700)
    S = np.array([entropy_at_tau(eig_vals, t) for t in tau_grid]) / lnN_b
    log10_t = np.log10(tau_grid)
    C = -np.gradient(S, log10_t)
    return tau_grid, 1.0 - S, C


def main():
    apply_pub_style()
    patient, phase, band = "Pat_08", "task_learn", "alpha"

    COLOR_ENT = "#4363d8"
    COLOR_C = "#e6194b"
    LS = {"msc": "-", "imcoh_abs": "--"}
    DISPLAY = {"msc": "MSC", "imcoh_abs": r"$|\mathrm{ImCoh}|$"}

    fig, ax_ent = plt.subplots(figsize=(10, 5))
    ax_C = ax_ent.twinx()

    for method in ("msc", "imcoh_abs"):
        A = load_fc_matrix(patient, phase, band, method)
        ev = laplacian_eigvals(A)
        tau, one_minus_S, C = compute_curves(ev)

        ax_ent.plot(tau, one_minus_S, ls=LS[method], color=COLOR_ENT, lw=2.5,
                    label=fr"$1-\tilde S$ {DISPLAY[method]}")
        ax_C.plot(tau, C, ls=LS[method], color=COLOR_C, lw=2.0,
                  label=fr"$\tilde C$ {DISPLAY[method]}")

    ax_ent.set_xscale("log")
    ax_ent.set_xlabel(r"Diffusion time $\tau$", fontsize=15)
    ax_ent.set_ylabel(r"$1 - \tilde S(\tau)$", fontsize=15, color=COLOR_ENT)
    ax_C.set_ylabel(r"$\tilde C(\tau) = -d\tilde S/d\log_{10}\tau$",
                    fontsize=15, color=COLOR_C)
    ax_ent.tick_params(axis="y", labelcolor=COLOR_ENT, labelsize=12)
    ax_C.tick_params(axis="y", labelcolor=COLOR_C, labelsize=12)
    ax_ent.tick_params(axis="x", labelsize=12)

    h1, l1 = ax_ent.get_legend_handles_labels()
    h2, l2 = ax_C.get_legend_handles_labels()
    ax_ent.legend(h1 + h2, l1 + l2, fontsize=11)

    out = SECTION2_ROOT / "fig_G" / f"fig_G2_susceptibility_{patient}_{band}_{phase}_MSC_vs_ImCoh.pdf"
    save_fig(fig, out)
    print(f"OK wrote {out}")


if __name__ == "__main__":
    main()
