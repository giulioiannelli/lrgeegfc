#!/usr/bin/env python3
"""Pillar 6 -- does implant geometry explain the inter-patient trace spread?

Joins the mst@0.20 per-patient trace (ms_mst020, at the mesoscale s=5.6 where alpha/
beta are jointly strongest) with the implant-geometry features (audit implant_geometry:
N_L/N_R hemisphere counts, centroid, lobe coverage) and asks whether the spread in
who-traces tracks electrode laterality -- reproducing audit_148 (laterality = the
driver of beta heterogeneity) on the recovered trace. No FC/LRG recompute.

Outputs (data/sparsified_arc/patient_variation/):
  master.csv, correlations.csv, fig_patient_variation.pdf
"""
from __future__ import annotations
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

MS = ROOT / "data" / "sparsified_arc" / "ms_mst020" / "per_patient_scale.csv"
GEO = ROOT / "data" / "audit" / "implant_geometry" / "per_patient_features.csv"
OUT = ROOT / "data" / "sparsified_arc" / "patient_variation"
S_REPORT = 5.6
BANDS = ["alpha", "beta"]


def main():
    use_lrg_style()
    OUT.mkdir(parents=True, exist_ok=True)
    ms = pd.read_csv(MS)
    geo = pd.read_csv(GEO).set_index("patient")
    s = np.sort(ms.s.unique()); sN = float(s[np.argmin(np.abs(s - S_REPORT))])

    M = pd.DataFrame(index=geo.index)
    for b in BANDS:
        x = ms[(ms.band == b) & np.isclose(ms.s, sN)].set_index("patient")
        M[f"{b}_rho"] = x["obs_rho"]; M[f"{b}_p"] = x["p"]
    for c in ["N", "N_L", "N_R", "N_M", "B_hemi", "centroid_x", "N_frontal", "N_epi"]:
        if c in geo.columns:
            M[c] = geo[c]
    M["left_frac"] = M["N_L"] / (M["N_L"] + M["N_R"]).replace(0, np.nan)
    M["right_frac"] = M["N_R"] / (M["N_L"] + M["N_R"]).replace(0, np.nan)
    M.to_csv(OUT / "master.csv")

    covs = ["left_frac", "right_frac", "N_R", "N_L", "B_hemi", "centroid_x", "N_frontal", "N"]
    rows = []
    for b in BANDS:
        y = M[f"{b}_rho"].astype(float)
        for c in covs:
            x = M[c].astype(float); ok = x.notna() & y.notna()
            r, p = spearmanr(x[ok], y[ok])
            rows.append(dict(band=b, covariate=c, rho=float(r), p=float(p), n=int(ok.sum())))
    corr = pd.DataFrame(rows)
    corr.to_csv(OUT / "correlations.csv", index=False)

    # figure: trace vs left_frac, per band
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, b in zip(axes, BANDS):
        col = band_color(b)
        for pat, r in M.iterrows():
            sig = r[f"{b}_p"] < 0.05
            ax.scatter(r.left_frac, r[f"{b}_rho"], s=90,
                       c=col if sig else "white", edgecolor=col, lw=1.8, zorder=3)
            ax.annotate(pat.replace("Pat_", ""), (r.left_frac, r[f"{b}_rho"]),
                        fontsize=8, xytext=(5, 3), textcoords="offset points")
        rr, pp = spearmanr(M.left_frac, M[f"{b}_rho"])
        ax.axhline(0, color="0.7", lw=0.7)
        ax.set_xlabel("left-hemisphere contact fraction  $N_L/(N_L{+}N_R)$")
        ax.set_ylabel(fr"{b} trace $\rho_{{\rm sym}}$ (s={sN:.1f})")
        ax.set_title(fr"$\bf{{{b[0]}}}$  {b}: $\rho$={rr:+.2f}, p={pp:.2f}",
                     color=col, fontweight="bold", fontsize=12, loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "fig_patient_variation.pdf", transparent=True); plt.close(fig)

    print(f"=== pillar 6: implant geometry vs mst@0.20 trace (s={sN:.1f}) ===")
    print(M[["beta_rho", "beta_p", "alpha_rho", "left_frac", "N_R", "N_L", "centroid_x"]]
          .sort_values("beta_rho", ascending=False).round(3).to_string())
    print("\nTop covariates by |rho| (per band):")
    print(corr.reindex(corr.rho.abs().sort_values(ascending=False).index)
          .groupby("band").head(4).round(3).to_string(index=False))
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
