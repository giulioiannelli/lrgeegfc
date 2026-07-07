"""audit_121b — figure for the τ-sensitivity-of-the-trace diagnostic (audit_121).

2×2: (a) β headline — ρ_split, T_d, placebo vs α with per-patient spaghetti;
(b) ρ_split(α) all bands; (c) T_d(α) all bands; (d) self-similarity G_self(α) +
collapsed fraction. α=1 (current convention) marked in every panel; the
collapse / placebo-contaminated zone is shaded out.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.notebook import move_to_rootf
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()
move_to_rootf(pathname="lrgeegfc")

OUT = Path("data/audit/tau_sweep_trace")
pp = pd.read_csv(OUT / "per_patient.csv")
co = pd.read_csv(OUT / "cohort.csv")

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
LAB = {"delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
       "beta": r"$\beta$", "low_gamma": r"$\gamma_{\mathrm{lo}}$",
       "high_gamma": r"$\gamma_{\mathrm{hi}}$"}
COL = {"delta": "#1f3a93", "theta": "#159a8c", "alpha": "#2e8b3d",
       "beta": "#c0392b", "low_gamma": "#e07b00", "high_gamma": "#7d3c98"}


def band_curve(band, col):
    s = co[co.band == band].sort_values("alpha")
    return s["alpha"].values, s[col].values, s


def collapse_alpha(band):
    """First α where the geometry stops being trustworthy."""
    s = co[co.band == band].sort_values("alpha")
    bad = s[(s.frac_collapsed > 0.05) | (s.median_G_self < 0.78)]
    return float(bad.alpha.min()) if len(bad) else float(s.alpha.max())


fig, axes = plt.subplots(2, 2, figsize=(11, 8))
ax_a, ax_b, ax_c, ax_d = axes.ravel()

# ---- (a) β headline ----------------------------------------------------------
b = "beta"
aF = float(co[co.band == b].median_alpha_fiedler.iloc[0])
ca = collapse_alpha(b)
for pat, g in pp[pp.band == b].groupby("patient"):
    g = g.sort_values("alpha")
    ax_a.plot(g.alpha, g.rho_split_coph, color=COL[b], lw=0.6, alpha=0.22)
al, rs, _ = band_curve(b, "median_rho_split_coph")
_, td, _ = band_curve(b, "median_T_d_coph")
_, ri, _ = band_curve(b, "median_rho_indep_coph")
ax_a.plot(al, rs, color=COL[b], lw=2.6, label=r"$\rho_{\mathrm{split}}$ (trace)")
ax_a.plot(al, td, color="#222222", lw=2.0, label=r"$T_d$ (trace)")
ax_a.plot(al, ri, color="#888888", lw=1.8, ls="--",
          label=r"$\rho_{\mathrm{indep}}$ (placebo)")
ax_a.axhline(0, color="k", lw=0.6, alpha=0.4)

# ---- (b) ρ_split all bands ; (c) T_d all bands -------------------------------
for ax, col, ttl in [(ax_b, "median_rho_split_coph", "rho_split"),
                     (ax_c, "median_T_d_coph", "T_d")]:
    for band in BANDS:
        al, y, _ = band_curve(band, col)
        lw = 2.6 if band == "beta" else 1.6
        ax.plot(al, y, color=COL[band], lw=lw, label=LAB[band])
    ax.axhline(0, color="k", lw=0.6, alpha=0.4)

# ---- (d) self-similarity + collapsed fraction --------------------------------
for band in BANDS:
    al, g, _ = band_curve(band, "median_G_self")
    ax_d.plot(al, g, color=COL[band], lw=2.6 if band == "beta" else 1.6)
ax_d.set_ylabel(r"self-similarity $\rho^{\mathrm{coph}}(\tau,\,1/\lambda_{\max})$")
ax_dr = ax_d.twinx()
fc = co.groupby("alpha").frac_collapsed.mean().sort_index()
ax_dr.fill_between(fc.index, 0, fc.values, color="#bbbbbb", alpha=0.35, lw=0)
ax_dr.set_ylabel("fraction collapsed", color="#777777")
ax_dr.set_ylim(0, 1)
ax_dr.tick_params(axis="y", colors="#777777")

# ---- shared cosmetics --------------------------------------------------------
panel_tags = ["a", "b", "c", "d"]
ylabels = [r"trace statistic", r"$\rho_{\mathrm{split}}$ (cophenetic)",
           r"$T_d$ (cophenetic)", None]
for k, ax in enumerate([ax_a, ax_b, ax_c, ax_d]):
    ax.set_xscale("log")
    ax.set_xlabel(r"$\alpha=\tau\,\lambda_{\max}$   (current convention: $\alpha=1$)")
    if ylabels[k]:
        ax.set_ylabel(ylabels[k])
    ax.axvline(1.0, color="k", lw=1.2, ls=":")            # current τ=1/λmax
    # collapse / contamination shade (cohort-wide onset ~ where >10% collapsed)
    a_bad = float(co.groupby("alpha").frac_collapsed.mean().pipe(
        lambda s: s[s > 0.1].index.min() if (s > 0.1).any() else co.alpha.max()))
    ax.axvspan(a_bad, co.alpha.max(), color="#dddddd", alpha=0.5, lw=0)
    ax.text(0.02, 0.95, panel_tags[k], transform=ax.transAxes,
            fontweight="bold", va="top", ha="left")

# Fiedler marker on the all-band panels (cohort-median per band is similar ~3)
for ax in (ax_b, ax_c, ax_d):
    ax.axvline(co.median_alpha_fiedler.median(), color="k", lw=0.8,
               ls="--", alpha=0.4)

h_a, l_a = ax_a.get_legend_handles_labels()
fig.legend(h_a, l_a, loc="upper center", bbox_to_anchor=(0.27, 1.0),
           ncol=3, frameon=False, fontsize=9)
h_b = [plt.Line2D([0], [0], color=COL[bd], lw=2.4) for bd in BANDS]
fig.legend(h_b, [LAB[bd] for bd in BANDS], loc="upper center",
           bbox_to_anchor=(0.74, 1.0), ncol=6, frameon=False, fontsize=9)

fig.tight_layout(rect=[0, 0, 1, 0.95])
out = OUT / "audit_121_tau_sweep_trace.pdf"
fig.savefig(out, transparent=True, bbox_inches="tight")
plt.close(fig)
print(f"Wrote {out}")
print(f"  β Fiedler α≈{aF:.1f}, collapse-onset α≈{ca:.1f}")
