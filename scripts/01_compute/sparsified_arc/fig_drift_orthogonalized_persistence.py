#!/usr/bin/env python3
"""Test 2 figure -- drift-orthogonalised persistence: obs(post) vs held-out-rest null vs r.

One figure per space (fc, coph). Per band (2x3): cohort median T_obs(r) (solid, post
retains the task-specific direction) vs T_null(r) (dashed, a task-naive held-out rest
segment), across the repertoire dimension r (drift-immunity knob). r where the cohort
gate (Wilcoxon obs-null, p<0.05) fires is ticked; novelty (fraction of task change
outside the resting repertoire) annotated. A drift-immune trace = obs stays above null
as r grows. Reads orth_persist/gate_vs_r.csv.
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

OP = ROOT / "data" / "sparsified_arc" / "trace_arc" / "orth_persist"
FIGS = ROOT / "data" / "sparsified_arc" / "figures" / "orth_persist"
SPACE_LABEL = {"fc": "FC-edge (drift-clean)", "coph": "cophenetic@τ_min (rank)"}


def main():
    use_lrg_style()
    gv = pd.read_csv(OP / "gate_vs_r.csv")
    FIGS.mkdir(parents=True, exist_ok=True)
    for sp in ("fc", "coph"):
        g = gv[gv.space == sp]
        if g.empty:
            continue
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        for ax, band in zip(axes.ravel(), BANDS):
            x = g[g.band == band].sort_values("r")
            if x.empty:
                ax.set_axis_off(); continue
            r = x.r.values; col = band_color(band)
            ax.plot(r, x.T_obs_med, "-o", color=col, lw=2.2, ms=4,
                    label=r"$T_\perp$ obs (post)")
            ax.plot(r, x.T_null_med, "--s", color="0.45", lw=1.6, ms=3,
                    label="held-out rest null")
            ax.fill_between(r, x.T_null_med, x.T_obs_med,
                            where=(x.T_obs_med.values > x.T_null_med.values),
                            color=col, alpha=0.12, lw=0)
            beat = x.gate_p.values < 0.05
            if beat.any():
                yb = np.nanmax(x.T_obs_med.values) * 1.05
                ax.plot(r[beat], np.full(beat.sum(), yb), "v", color=col, ms=6)
            ax.axhline(0, color="0.7", lw=0.6)
            ax.set_xlabel(r"repertoire dim $r$ (drift removed)")
            ax.set_ylabel(r"$T_\perp$ persistence")
            ax.set_xticks(r)
            bp = float(x.gate_p.min()) if x.gate_p.notna().any() else np.nan
            nov = float(x.novelty_med.iloc[min(1, len(x) - 1)])
            verdict = "trace" if (np.isfinite(bp) and bp < 0.05) else "n.s."
            ax.set_title(f"{band} — best p={bp:.3f} [{verdict}], novelty≈{nov:.2f}",
                         fontweight="bold", fontsize=9.5, color=col)
            ax.legend(frameon=False, fontsize=7.5, loc="upper right")
        fig.tight_layout()
        out = FIGS / f"orth_persist_{sp}.pdf"
        fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


if __name__ == "__main__":
    main()
