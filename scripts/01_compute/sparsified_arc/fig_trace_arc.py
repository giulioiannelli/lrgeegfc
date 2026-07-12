#!/usr/bin/env python3
"""D2 figures -- the tau-resolved trace, per band (the "role of tau" money figure).

  1. trace_vs_tau_summary.pdf -- one panel per band: cohort-median observed
     rho_sym(s) (colour) with IQR, matched-strength surrogate p50-p95 band (grey),
     firing scales (cohort gate p<0.05) marked. Title = classification +
     s=1 (tau_min) gate p + scale-max gate p. Instantly shows single- vs
     multi- vs continuous-multiscale.
  2. trace_gate_heatmap.pdf -- band x scale grid coloured by -log10 cohort gate p;
     the whole tau-trace landscape in one image.
  3. trace_curves_{band}.pdf -- per-patient control grid: obs(s) vs surrogate band.
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

ARC = ROOT / "data" / "sparsified_arc" / "trace_arc"
FIGS = ROOT / "data" / "sparsified_arc" / "figures" / "trace_arc"


def _cell(pat, band):
    f = ARC / band / f"{pat}.npz"
    return np.load(f) if f.exists() else None


def _stack(band, key, s_common):
    out = []
    for pat in COHORT:
        z = _cell(pat, band)
        if z is None:
            continue
        out.append(np.interp(s_common, z["s"], z[key], left=np.nan, right=np.nan))
    return np.array(out) if out else np.full((1, s_common.size), np.nan)


def _effect_size_class(obs_med, n_above95, n_pat):
    """Honest label from EFFECT SIZE, not the per-scale Wilcoxon (which rewards a
    small consistent offset). A scale 'traces' if a majority of patients exceed
    their OWN null p95 there. frac = fraction of scales that do."""
    fire = n_above95 >= np.ceil(n_pat / 2)
    frac = float(fire.mean())
    # contiguous firing bands
    nb = 0; prev = False
    for f in fire:
        if f and not prev:
            nb += 1
        prev = f
    if frac < 0.08:
        cls = "marginal ≈ null" if np.nanmax(n_above95) >= 2 else "none"
    elif frac >= 0.6:
        cls = "broad multiscale"
    elif nb >= 2:
        cls = "multi-band"
    else:
        cls = "mesoscale"
    return cls, frac, fire


def summary_fig(nature, gate):
    use_lrg_style()
    s_common = np.logspace(0, np.log10(200), 200)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    corrected = []
    for ax, band in zip(axes.ravel(), BANDS):
        col = band_color(band)
        OBS = _stack(band, "obs", s_common)
        P50 = _stack(band, "surr_p50", s_common)
        P95 = _stack(band, "surr_p95", s_common)
        obs_med = np.nanmedian(OBS, 0)
        q1, q3 = np.nanpercentile(OBS, [25, 75], axis=0)
        ax.fill_between(s_common, q1, q3, color=col, alpha=0.18, lw=0)
        ax.plot(s_common, obs_med, "-", color=col, lw=2.4, label=r"obs $\rho_{sym}(\tau)$")
        ax.fill_between(s_common, np.nanmedian(P50, 0), np.nanmedian(P95, 0),
                        color="0.6", alpha=0.35, lw=0, label="null p50–p95")
        ax.plot(s_common, np.nanmedian(P50, 0), ":", color="0.4", lw=1.2)
        ax.axhline(0, color="0.8", lw=0.7); ax.axvline(1.0, color="0.6", lw=0.8, ls=":")

        # EFFECT-SIZE firing: per-patient obs > own null p95, majority of cohort
        above95 = (OBS > P95).astype(float)
        n_above95 = np.nansum(above95, axis=0)
        n_pat = np.sum(np.isfinite(OBS).any(axis=1))
        cls, frac95, fire = _effect_size_class(obs_med, n_above95, n_pat)
        ytop = ax.get_ylim()[1]
        ax.plot(s_common[fire], np.full(fire.sum(), ytop * 0.96), marker="v",
                color=col, ms=4, ls="none")
        row = nature[nature.band == band].iloc[0]
        corrected.append(dict(band=band, effect_class=cls, frac_scales_obs_gt_p95=round(frac95, 3),
                              max_n_pat_above_p95=int(np.nanmax(n_above95)),
                              gate_p_scalemax=float(row["gate_p_scalemax"])))
        ax.set_title(f"{band}  [{cls}]   "
                     f"{int(np.nanmax(n_above95))}/{int(n_pat)} pt > own null   "
                     f"$p_{{smax}}$={row['gate_p_scalemax']:.3f}",
                     fontweight="bold", fontsize=10, color=col)
        ax.set_xscale("log"); ax.set_xlabel(r"scale $s=\tau\,\lambda_{\max}$")
        ax.set_ylabel(r"cohort $\rho_{sym}$")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "trace_vs_tau_summary.pdf"
    fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")
    pd.DataFrame(corrected).to_csv(ARC / "band_effectsize_class.csv", index=False)
    print("[corrected classification]\n" + pd.DataFrame(corrected).to_string(index=False))


def gate_heatmap(gate):
    use_lrg_style()
    s_vals = np.sort(gate.s.unique())
    M = np.full((len(BANDS), len(s_vals)), np.nan)
    for i, b in enumerate(BANDS):
        gb = gate[gate.band == b]
        for _, r in gb.iterrows():
            j = int(np.argmin(np.abs(s_vals - r["s"])))
            M[i, j] = -np.log10(max(r["gate_p"], 1e-4)) if np.isfinite(r["gate_p"]) else np.nan
    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(M, aspect="auto", cmap="magma", norm=Normalize(0, 3),
                   extent=[np.log10(s_vals.min()), np.log10(s_vals.max()), len(BANDS) - 0.5, -0.5])
    ax.set_yticks(range(len(BANDS))); ax.set_yticklabels(BANDS)
    ax.set_xlabel(r"$\log_{10}\,s=\log_{10}(\tau\lambda_{\max})$")
    ax.axvline(0.0, color="w", lw=1.0, ls=":")   # s=1
    cb = fig.colorbar(im, ax=ax, label=r"$-\log_{10}$ cohort gate $p$")
    cb.ax.axhline(-np.log10(0.05), color="c", lw=1.5)
    ax.set_title(r"$\tau$-trace landscape (matched-strength; white dotted $=\tau_{\min}$)",
                 fontweight="bold", fontsize=10)
    fig.tight_layout()
    out = FIGS / "trace_gate_heatmap.pdf"
    fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


def per_band_curves():
    use_lrg_style()
    for band in BANDS:
        col = band_color(band)
        fig, axes = plt.subplots(2, 5, figsize=(18, 6.5), sharex=True)
        for ax, pat in zip(axes.ravel(), COHORT):
            z = _cell(pat, band)
            if z is None:
                ax.set_axis_off(); continue
            ax.fill_between(z["s"], z["surr_p50"], z["surr_p95"], color="0.6", alpha=0.35, lw=0)
            ax.plot(z["s"], z["obs"], "-", color=col, lw=1.8)
            ax.plot(z["s"], z["surr_p50"], ":", color="0.4", lw=1.0)
            ax.axhline(0, color="0.85", lw=0.7); ax.axvline(1.0, color="0.6", lw=0.7, ls=":")
            ax.set_xscale("log"); ax.set_title(pat.replace("Pat_", "P"), fontsize=9)
        fig.supxlabel(r"scale $s=\tau\lambda_{\max}$   —   obs $\rho_{sym}$ colour, null p50–p95 grey", fontsize=9)
        fig.supylabel(r"$\rho_{sym}$", fontsize=9)
        fig.tight_layout()
        out = FIGS / f"trace_curves_{band}.pdf"
        fig.savefig(out, transparent=True); plt.close(fig); print(f"[fig] {out}")


def main():
    nature = pd.read_csv(ARC / "band_nature.csv")
    gate = pd.read_csv(ARC / "cohort_gate_vs_s.csv")
    summary_fig(nature, gate)
    gate_heatmap(gate)
    per_band_curves()
    print(f"[D2-figs] -> {FIGS}")


if __name__ == "__main__":
    main()
