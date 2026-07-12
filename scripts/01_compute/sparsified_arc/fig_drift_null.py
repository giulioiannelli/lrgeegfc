#!/usr/bin/env python3
"""D2-drift figure -- the decisive control: whole-task trace vs temporal-drift floor.

Per band (2x3), cohort obs rho_sym(s) (solid) vs the within-rest_pre drift floor
(dashed) with the per-patient drift spread (p50-p95 shaded). A genuine task trace
must sit ABOVE the drift band. Scales where the cohort obs beats drift
(Wilcoxon obs-drift, greater, p<0.05) are ticked; the per-band scale-max drift
gate p is in the title. Mirrors trace_vs_tau_summary so the two are comparable.

Reads data/sparsified_arc/trace_arc/drift/{band}/{pat}.npz (s, obs, drift).
"""
from __future__ import annotations
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

DR = ROOT / "data" / "sparsified_arc" / "trace_arc" / "drift"
FIGS = ROOT / "data" / "sparsified_arc" / "figures" / "drift"


def _load_band(band):
    s_ref, obs, dft = None, [], []
    for pat in COHORT:
        f = DR / band / f"{pat}.npz"
        if not f.exists():
            continue
        z = np.load(f)
        if s_ref is None:
            s_ref = z["s"]
        obs.append(np.interp(s_ref, z["s"], z["obs"]))
        dft.append(np.interp(s_ref, z["s"], z["drift"]))
    return s_ref, np.array(obs), np.array(dft)


def main():
    use_lrg_style()
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, band in zip(axes.ravel(), BANDS):
        s, OBS, DFT = _load_band(band)
        if s is None:
            ax.set_axis_off(); continue
        col = band_color(band)
        obs_med = np.nanmedian(OBS, 0)
        d_p50 = np.nanpercentile(DFT, 50, 0)
        d_p95 = np.nanpercentile(DFT, 95, 0)
        # per-scale gate: obs beats drift?
        gate = np.full(s.size, np.nan)
        for j in range(s.size):
            d = OBS[:, j] - DFT[:, j]; d = d[np.isfinite(d)]
            if d.size >= 5 and np.any(d != 0):
                try:
                    gate[j] = wilcoxon(d, alternative="greater")[1]
                except Exception:
                    pass
        # scale-max drift gate
        try:
            psmax = wilcoxon(np.nanmax(OBS, 1) - np.nanmax(DFT, 1),
                             alternative="greater")[1]
        except Exception:
            psmax = np.nan
        n_pat = int((np.nanmax(OBS, 1) > np.nanmax(DFT, 1)).sum())

        ax.fill_between(s, d_p50, d_p95, color="0.6", alpha=0.35, lw=0,
                        label="drift floor p50–p95")
        ax.plot(s, d_p50, "--", color="0.35", lw=1.4, label="drift median")
        ax.plot(s, obs_med, "-", color=col, lw=2.4, label=r"obs $\rho_{sym}(\tau)$")
        beat = np.isfinite(gate) & (gate < 0.05)
        if beat.any():
            ytick = ax.get_ylim()[1]
            ax.plot(s[beat], np.full(beat.sum(), obs_med.max() * 1.05), "v",
                    color=col, ms=4)
        ax.axhline(0, color="0.7", lw=0.6)
        ax.axvline(1.0, color="0.6", lw=0.7, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(r"scale $s=\tau\,\lambda_{\max}$")
        ax.set_ylabel(r"cohort $\rho_{sym}$")
        verdict = "beats drift" if (np.isfinite(psmax) and psmax < 0.05) else "= DRIFT"
        ax.set_title(f"{band}  —  {n_pat}/10 pt > drift,  $p_{{smax}}$={psmax:.3f}  [{verdict}]",
                     fontweight="bold", fontsize=9.5, color=col)
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "drift_vs_tau_summary.pdf"
    fig.savefig(out, transparent=True); plt.close(fig)
    print(f"[fig] {out}")


if __name__ == "__main__":
    main()
