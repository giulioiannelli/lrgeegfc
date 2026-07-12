#!/usr/bin/env python3
"""FAIR drift-null figure -- per-scale, distribution-based, SNR-matched.

Per band (2x3): cohort median obs_matched rho_sym(s) (solid) vs the FAIR drift
distribution p50-p95 band (grey), with the full-length D2 obs as a faint dashed
reference. Scales where the cohort obs_matched beats the drift p50 (per-scale
Wilcoxon, p<0.05) are ticked; the emergence scale (argmax obs_matched) is marked.
Mirrors trace_vs_tau_summary so the FAIR drift verdict is read the same way the
matched-strength arc was -- and, unlike 07's figure, uses NO scale-max statistic.

Reads data/sparsified_arc/trace_arc/drift_fair/gate_vs_s.csv (cohort curves+gate)
and the per-cell npz (full-length obs reference).
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

DF = ROOT / "data" / "sparsified_arc" / "trace_arc" / "drift_fair"
FIGS = ROOT / "data" / "sparsified_arc" / "figures" / "drift_fair"


def _obs_full_med(band):
    s_ref, vals = None, []
    for pat in COHORT:
        f = DF / band / f"{pat}.npz"
        if not f.exists():
            continue
        z = np.load(f)
        if s_ref is None:
            s_ref = z["s"]
        vals.append(np.interp(s_ref, z["s"], z["obs_full"]))
    return s_ref, (np.nanmedian(vals, 0) if vals else None)


def main():
    use_lrg_style()
    gv = pd.read_csv(DF / "gate_vs_s.csv")
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, band in zip(axes.ravel(), BANDS):
        x = gv[gv.band == band].sort_values("s")
        if x.empty:
            ax.set_axis_off(); continue
        s = x.s.values
        col = band_color(band)
        ax.fill_between(s, x.drift_p50_med, x.drift_p95_med, color="0.6", alpha=0.35,
                        lw=0, label="fair drift p50–p95")
        ax.plot(s, x.drift_p50_med, "--", color="0.35", lw=1.3, label="drift median")
        ax.plot(s, x.obs_matched_med, "-", color=col, lw=2.4,
                label=r"obs$_{\rm matched}\,\rho_{sym}(\tau)$")
        sf, of = _obs_full_med(band)
        if of is not None:
            ax.plot(sf, of, ":", color=col, lw=1.2, alpha=0.6,
                    label=r"obs$_{\rm full}$ (ref)")
        beat = x.gate_p.values < 0.05
        if beat.any():
            yb = np.nanmax(x.obs_matched_med.values) * 1.06
            ax.plot(s[beat], np.full(beat.sum(), yb), "v", color=col, ms=4)
        j = int(np.nanargmax(x.obs_matched_med.values))
        ax.axvline(s[j], color=col, lw=0.8, ls="-", alpha=0.4)
        ax.axhline(0, color="0.7", lw=0.6); ax.axvline(1.0, color="0.6", lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"scale $s=\tau\,\lambda_{\max}$")
        ax.set_ylabel(r"cohort $\rho_{sym}$")
        n_beat = int(beat.sum())
        ge = float(x.gate_p.values[j])
        verdict = "beats drift" if n_beat else "= drift"
        ax.set_title(f"{band}  —  {n_beat}/{len(s)} scales beat drift,  "
                     f"gate@emerge(s={s[j]:.0f})={ge:.3f}  [{verdict}]",
                     fontweight="bold", fontsize=9, color=col)
        ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "drift_fair_vs_tau_summary.pdf"
    fig.savefig(out, transparent=True); plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
