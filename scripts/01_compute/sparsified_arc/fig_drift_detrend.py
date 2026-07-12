#!/usr/bin/env python3
"""DETREND figure -- does the trace survive removing the session-long linear FC drift?

Per band (2x3): cohort median rho_sym BEFORE detrend (rho_base, solid) vs AFTER
removing the per-edge linear session drift (rho_detrend, dashed), at full coverage.
A genuine (non-drift) trace stays positive after detrend; a drift-driven one collapses
to <=0. Scales where detrended rho_sym is still cohort-significant (>0, p<0.05) are
ticked. Reads drift_detrend/gate_vs_s.csv.
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import BANDS
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

DD = ROOT / "data" / "sparsified_arc" / "trace_arc" / "drift_detrend"
FIGS = ROOT / "data" / "sparsified_arc" / "figures" / "drift_detrend"


def main():
    use_lrg_style()
    gv = pd.read_csv(DD / "gate_vs_s.csv")
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, band in zip(axes.ravel(), BANDS):
        x = gv[gv.band == band].sort_values("s")
        if x.empty:
            ax.set_axis_off(); continue
        s = x.s.values; col = band_color(band)
        ax.plot(s, x.rho_base_med, "-", color=col, lw=2.4, label=r"$\rho_{sym}$ base")
        ax.plot(s, x.rho_detrend_med, "--", color=col, lw=1.8, alpha=0.8,
                label=r"$\rho_{sym}$ detrended")
        ax.fill_between(s, x.rho_detrend_med, x.rho_base_med, color=col, alpha=0.12, lw=0)
        beat = x.p_detrend_pos.values < 0.05
        if beat.any():
            ax.plot(s[beat], np.full(beat.sum(), np.nanmax(x.rho_base_med) * 1.05),
                    "v", color=col, ms=4)
        ax.axhline(0, color="0.5", lw=0.8); ax.axvline(1.0, color="0.6", lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"scale $s=\tau\,\lambda_{\max}$"); ax.set_ylabel(r"cohort $\rho_{sym}$")
        bp = float(x.p_detrend_pos.min())
        verdict = "survives" if bp < 0.05 else "= DRIFT (dies)"
        ax.set_title(f"{band}  —  detrend best p={bp:.3f}  [{verdict}]",
                     fontweight="bold", fontsize=9.5, color=col)
        ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "drift_detrend_summary.pdf"
    fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


if __name__ == "__main__":
    main()
