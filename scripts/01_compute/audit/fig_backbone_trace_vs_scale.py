#!/usr/bin/env python3
"""Cross-phase trace rho_sym(tau) vs the specific-heat scale, on a cycle-rich backbone.

The user's multiscale goal, precisely: on a moderate-density backbone that KEEPS
non-trivial cycles (NOT a tree), sweep the ONE parameter we like -- tau -- and ask
whether the cross-phase trace rho_sym appears at a CHARACTERISTIC SCALE, and whether
that scale coincides with the mesoscale peak of the specific heat C(tau).

Eigendecompose each phase's backbone Laplacian ONCE; every scale s=tau*lambda_max is
then just K=V exp(-(s/lmax)Lambda) V^T -> cophenetic -> rho_sym. Cohort median +/- IQR.
Overlaid: cohort-median C(tau). d in {0.10, 0.20} (cycle-rich, plateau of the gate).

Output: data/audit/sparse_backbone_propagator/trace_vs_scale.pdf
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase, rho_sym_split, COHORT
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction
from lrg_eegfc.workflow.diagnostics import compute_entropy_curve
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.visuals.styles import band_color

OUT = ROOT / "data" / "audit" / "sparse_backbone_propagator"
PHASES = ("A", "B", "task_test", "rest_post")
S_GRID = np.logspace(0.0, np.log10(200.0), 28)      # s = tau*lambda_max
RHO_FLOOR = 1e-30


def eig_phase(W):
    deg = W.sum(1); ev, V = np.linalg.eigh(np.diag(deg) - W)
    return np.maximum(ev, 0.0), V


def coph_at_s(ev, V, s):
    tau = s / ev[-1]
    rho = (V * np.exp(-tau * ev)) @ V.T; rho /= np.trace(rho)
    rho = np.where(rho > RHO_FLOOR, rho, RHO_FLOOR)
    T = 1.0 / rho; np.fill_diagonal(T, 0.0); T = np.maximum(T, T.T)
    return cophenet(linkage(squareform(T, checks=False), method="average"))


def sweep_cell(pat, band, d):
    A = {ph: mst_union_top_fraction(load_phase(pat, ph, band), d) for ph in PHASES}
    eig = {ph: eig_phase(A[ph]) for ph in PHASES}
    rho = np.full(len(S_GRID), np.nan)
    for i, s in enumerate(S_GRID):
        D = {ph: coph_at_s(*eig[ph], s) for ph in PHASES}
        rho[i], _ = rho_sym_split(D["A"], D["B"], D["task_test"], D["rest_post"])
    ev = eig["A"][0]
    ec = compute_entropy_curve(ev)
    s_C = ec["tau"] * ev[-1]; C = ec["C"]
    C_on_grid = np.interp(S_GRID, s_C, C)
    return rho, C_on_grid


def main():
    use_lrg_style()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, band in zip(axes, ("beta", "alpha")):
        for d, ls in [(0.10, "-"), (0.20, "--")]:
            R = []; Cc = []
            for pat in COHORT:
                try:
                    r, c = sweep_cell(pat, band, d)
                except Exception:
                    continue
                R.append(r); Cc.append(c)
            R = np.array(R); Cc = np.array(Cc)
            med = np.nanmedian(R, 0)
            col = band_color(band)
            ax.plot(S_GRID, med, ls, color=col, lw=2.0,
                    label=fr"$\rho_{{sym}}$ d={d:.2f}")
            if d == 0.10:
                q1, q3 = np.nanpercentile(R, [25, 75], axis=0)
                ax.fill_between(S_GRID, q1, q3, color=col, alpha=0.15, lw=0)
                # C(tau) cohort median on twin axis
                ax2 = ax.twinx()
                ax2.plot(S_GRID, np.nanmedian(Cc, 0), ":", color="0.45", lw=1.4)
                ax2.set_ylabel(r"$C(\tau)$ (grey)", color="0.45", fontsize=9)
                ax2.tick_params(axis="y", labelcolor="0.45")
        ax.axhline(0, color="0.7", lw=0.8, zorder=0)
        ax.set_xscale("log")
        ax.set_xlabel(r"scale $s=\tau\,\lambda_{\max}$")
        ax.set_ylabel(r"cohort median $\rho_{sym}$")
        ax.set_title(f"{band}  —  trace vs scale (d=0.10/0.20 backbone)",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    out = OUT / "trace_vs_scale.pdf"
    fig.savefig(out, transparent=True); plt.close(fig)
    print(f"[fig] trace vs scale -> {out}")


if __name__ == "__main__":
    main()
