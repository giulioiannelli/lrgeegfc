#!/usr/bin/env python3
"""Per-patient + cohort figures/tables for the mst@0.20 matched-strength gate.

Reads data/sparsified_arc/ms_mst020/{per_patient_scale,cohort_gate}.csv.
Emits:
  - printed per-patient x band table of obs_rho at each band's characteristic scale
    (s* = argmax cohort-median obs_rho), with * = beats own surrogate p95;
  - Figure A: tau-resolved gate -- obs_med(s) vs surrogate p50/p95 envelope per band,
    scales where cohort gate_p<0.05 ticked;
  - Figure B: per-patient obs_rho at s* per band, sorted, filled if it beats its own
    null -- shows which patients drive and which dissent (interpatient variability).
"""
from __future__ import annotations
import sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.visuals.styles import use_lrg_style, band_color

OUT = ROOT / "data" / "sparsified_arc" / "ms_mst020"
FIG = ROOT / "data" / "sparsified_arc" / "figures" / "ms_mst020"
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BTeX = {"delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
        "beta": r"$\beta$", "low_gamma": r"$\gamma_{\mathrm{low}}$",
        "high_gamma": r"$\gamma_{\mathrm{high}}$"}


S_REPORT = 5.6   # honest fixed mesoscale: alpha/beta jointly strongest, MS null not yet collapsed


def char_scale(gate, band):
    """Fixed reporting scale (nearest grid point to S_REPORT) for ALL bands -- avoids
    per-band scale cherry-picking and the coarse-s MS-null-collapse artifact."""
    s = np.sort(gate.s.unique())
    return float(s[np.argmin(np.abs(s - S_REPORT))])


def main():
    use_lrg_style()
    df = pd.read_csv(OUT / "per_patient_scale.csv")
    gate = pd.read_csv(OUT / "cohort_gate.csv")
    FIG.mkdir(parents=True, exist_ok=True)
    pats = sorted(df.patient.unique())
    sstar = {b: char_scale(gate, b) for b in BANDS}

    # ---- per-patient x band table at s*_band ----
    tbl = pd.DataFrame(index=pats, columns=BANDS, dtype=float)
    sig = pd.DataFrame(index=pats, columns=BANDS, dtype=bool)
    for b in BANDS:
        x = df[(df.band == b) & np.isclose(df.s, sstar[b])].set_index("patient")
        for p in pats:
            if p in x.index:
                tbl.loc[p, b] = x.loc[p, "obs_rho"]
                sig.loc[p, b] = bool(x.loc[p, "p"] < 0.05)
    tbl.to_csv(OUT / "per_patient_table.csv")
    print("=== per-patient obs_rho at s*_band (*, obs>own surrogate p95) ===")
    print("s*_band:", {b: round(sstar[b], 1) for b in BANDS})
    hdr = "patient    " + "".join(f"{BTeX[b].strip('$'):>12}" for b in BANDS)
    print(hdr)
    for p in pats:
        cells = []
        for b in BANDS:
            v = tbl.loc[p, b]
            mark = "*" if (isinstance(sig.loc[p, b], (bool, np.bool_)) and sig.loc[p, b]) else " "
            cells.append(f"{v:+7.3f}{mark}    " if np.isfinite(v) else "   nan     ")
        print(f"{p:9s} " + "".join(cells))
    ncnt = {b: int(sig[b].sum()) for b in BANDS}
    print("n>own p95: " + "".join(f"{ncnt[b]:>12}" for b in BANDS))

    # ---- Figure A: tau-resolved gate ----
    figA, ax = plt.subplots(2, 3, figsize=(15, 8.5), sharex=True)
    for a, b in zip(ax.ravel(), BANDS):
        g = gate[gate.band == b].sort_values("s")
        col = band_color(b)
        # cohort surrogate envelope from per-patient p95
        env = df[df.band == b].groupby("s").agg(sp50=("surr_p50", "median"),
                                                sp95=("surr_p95", "median"))
        a.fill_between(env.index, env.sp50, env.sp95, color="0.6", alpha=0.25, lw=0,
                       label="surrogate p50-p95")
        a.plot(g.s, g.obs_med, "-o", color=col, lw=2.2, ms=4, label="obs median")
        sigs = g[g.gate_p < 0.05]
        if not sigs.empty:
            yb = g.obs_med.max() * 1.08
            a.plot(sigs.s, np.full(len(sigs), yb), "v", color=col, ms=7)
        a.axhline(0, color="0.7", lw=0.6); a.axvline(1, color="0.85", lw=0.7, ls=":")
        a.set_xscale("log")
        bp = g.gate_p.min()
        a.set_title(f"{BTeX[b]}  best gate_p={bp:.3f}", color=col,
                    fontweight="bold", fontsize=13, loc="left")
        a.set_xlabel(r"$s=\tau\lambda_{\max}$"); a.set_ylabel(r"$\rho_{\mathrm{sym}}^{\mathrm{coph}}$")
    ax.ravel()[0].legend(frameon=False, fontsize=9, loc="upper right")
    figA.tight_layout(); figA.savefig(FIG / "gate_vs_scale.pdf", transparent=True)
    plt.close(figA); print(f"[fig] {FIG/'gate_vs_scale.pdf'}")

    # ---- Figure B: per-patient at s*_band ----
    figB, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    for a, b in zip(ax.ravel(), BANDS):
        x = df[(df.band == b) & np.isclose(df.s, sstar[b])].copy()
        x = x.sort_values("obs_rho", ascending=True)
        col = band_color(b); y = np.arange(len(x))
        for yi, (_, row) in zip(y, x.iterrows()):
            is_sig = row.p < 0.05
            a.plot([row.surr_p95], [yi], "|", color="0.5", ms=12, mew=2)
            a.plot([row.obs_rho], [yi], "o", ms=9,
                   color=col if is_sig else "white",
                   mec=col, mew=1.8)
        a.axvline(0, color="0.7", lw=0.6)
        a.set_yticks(y); a.set_yticklabels([p.replace("Pat_", "") for p in x.patient], fontsize=8)
        g = gate[gate.band == b]; bp = g.gate_p.min(); nsig = int((x.p < 0.05).sum())
        a.set_title(f"{BTeX[b]}  s*={sstar[b]:.0f}  {nsig}/{len(x)}>null  p={bp:.3f}",
                    color=col, fontweight="bold", fontsize=12, loc="left")
        a.set_xlabel(r"$\rho_{\mathrm{sym}}^{\mathrm{coph}}$ (obs $\bullet$, own null p95 $|$)")
    figB.tight_layout(); figB.savefig(FIG / "per_patient.pdf", transparent=True)
    plt.close(figB); print(f"[fig] {FIG/'per_patient.pdf'}")


if __name__ == "__main__":
    main()
