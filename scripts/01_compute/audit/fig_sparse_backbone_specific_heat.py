#!/usr/bin/env python3
"""Specific heat C(tau) of sparsified FC backbones — does the multiscale ladder emerge?

Villegas 2025 App. A/F: C(tau) = -dS/dlog10(tau). Fully-connected high-<k> FC ->
Wigner-semicircle spectrum -> ONE collapse peak. Retaining non-trivial CYCLES
(moderate-density backbone, NOT a tree) can split it into a multi-peak mesoscale
ladder -- the telescopic-scanner regime the LRG is built for.

Shows C(tau) vs s = tau*lambda_max for a grid of densities (MST-union top-fraction,
which keeps cycles and stays connected) for representative (patient, band) cells.
The resolution window [s=1 (tau_min), s=lambda_max/lambda_2 (tau_max)] is shaded.

Output: data/audit/sparse_backbone_propagator/specific_heat_curves.pdf
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction, backbone_density
from lrg_eegfc.workflow.diagnostics import compute_entropy_curve
from lrg_eegfc.visuals.styles import use_lrg_style

OUT = ROOT / "data" / "audit" / "sparse_backbone_propagator"
CELLS = [("Pat_05", "beta"), ("Pat_08", "beta"), ("Pat_13", "beta"),
         ("Pat_06", "alpha"), ("Pat_02", "alpha"), ("Pat_05", "alpha")]
DENS = [1.0, 0.35, 0.20, 0.10, 0.05]     # full -> cycle-rich moderate; no tree
CMAP = plt.get_cmap("viridis")


def c_of_s(W):
    """C(tau) vs s=tau*lambda_max on the extended grid + resolution window + peak count."""
    deg = W.sum(1); ev = np.maximum(np.linalg.eigvalsh(np.diag(deg) - W), 0.0)
    ec = compute_entropy_curve(ev)
    s = ec["tau"] * ev[-1]                       # s = tau * lambda_max (s=1 == tau_min)
    C = ec["C"]
    s_lo, s_hi = 1.0, ev[-1] / ev[1]             # resolution window in s
    win = (ec["tau"] >= ec["tau_min"]) & (ec["tau"] <= ec["tau_max"])
    Cw = C.copy(); Cw[~win] = 0.0
    npk = len(find_peaks(Cw, prominence=0.03 * Cw.max())[0]) if Cw.max() > 0 else 0
    return s, C, s_lo, s_hi, npk


def main():
    use_lrg_style()
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2))
    for ax, (pat, band) in zip(axes.flat, CELLS):
        A = load_phase(pat, "rest_pre", band)
        for k, d in enumerate(DENS):
            B = A if d >= 1.0 else mst_union_top_fraction(A, d)
            s, C, s_lo, s_hi, npk = c_of_s(B)
            col = CMAP(0.12 + 0.8 * (1 - k / (len(DENS) - 1)))
            lab = f"full ({npk}pk)" if d >= 1.0 else f"d={backbone_density(B):.2f} ({npk}pk)"
            ax.plot(s, C / max(C.max(), 1e-9), color=col, lw=1.6, label=lab)
        ax.axvspan(1.0, s_hi, color="0.85", alpha=0.35, zorder=0)   # resolution window
        ax.set_xscale("log")
        ax.set_xlabel(r"$s=\tau\,\lambda_{\max}$")
        ax.set_ylabel(r"$C(\tau)$ (norm.)")
        ax.set_title(f"{pat}  {band}", fontsize=10, fontweight="bold")
        ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    fig.tight_layout()
    out = OUT / "specific_heat_curves.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print(f"[fig] specific heat curves -> {out}")


if __name__ == "__main__":
    main()
