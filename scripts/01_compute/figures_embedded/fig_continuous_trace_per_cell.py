#!/usr/bin/env python3
"""Per-cell visual of the continuous trace.

Four panels for one (patient, band):
  A  Δ_task = D_test − D_pre   (N×N, diverging colormap)
  B  Δ_rest = D_post − D_pre   (N×N, same colormap, same vmin/vmax)
  C  σ(i,j) = sgn(Δ_task)·sgn(Δ_rest) ∈ {−1,0,+1}  (categorical map)
  D  scatter Δ_task vs Δ_rest  (hex-bin, with cohort median dashed
                                 line and ρ annotation)

Reads:  data/reports/imcoh_continuous_trace/per_cell_summary.csv
        data/reports/imcoh_continuous_trace/per_pair/{patient}_{band}.npz
Writes: data/reports/imcoh_continuous_trace/figures/per_cell/{patient}_{band}.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_continuous_trace"
OUT_DIR = IN_DIR / "figures" / "per_cell"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def render_one(patient: str, band: str) -> None:
    summary = pd.read_csv(IN_DIR / "per_cell_summary.csv")
    row = summary[(summary["patient"] == patient) & (summary["band"] == band)]
    if row.empty:
        raise SystemExit(f"No per-cell row for {patient} {band}")
    rho = float(row["rho"].iloc[0])
    N = int(row["N_p"].iloc[0])
    m = int(row["m_p"].iloc[0])
    frac_pos = float(row["frac_pos_sigma"].iloc[0])
    frac_neg = float(row["frac_neg_sigma"].iloc[0])

    pair_path = IN_DIR / "per_pair" / f"{patient}_{band}.npz"
    pp = np.load(pair_path)
    dD_task = pp["dD_task"]
    dD_rest = pp["dD_rest"]
    sigma = pp["sigma"]
    iu_i = pp["iu_i"]
    iu_j = pp["iu_j"]

    # Reconstruct N×N square forms
    M_task = np.zeros((N, N))
    M_rest = np.zeros((N, N))
    M_sigma = np.zeros((N, N), dtype=np.int8)
    M_task[iu_i, iu_j] = dD_task; M_task[iu_j, iu_i] = dD_task
    M_rest[iu_i, iu_j] = dD_rest; M_rest[iu_j, iu_i] = dD_rest
    M_sigma[iu_i, iu_j] = sigma; M_sigma[iu_j, iu_i] = sigma

    # Symmetric vmin/vmax for the diverging panels
    vmax = float(max(np.abs(dD_task).max(), np.abs(dD_rest).max()))

    fig, axes = plt.subplots(1, 4, figsize=(18.0, 4.6), dpi=150)

    # --- Panel A: Δ_task ---
    ax = axes[0]
    im = ax.imshow(M_task, cmap="RdBu_r", vmin=-vmax, vmax=+vmax,
                   aspect="equal", rasterized=True)
    ax.set_title(r"A — $\Delta_{\mathrm{task}} = D_{\mathrm{test}} - D_{\mathrm{pre}}$",
                 fontsize=11)
    ax.set_xlabel("contact j", fontsize=9)
    ax.set_ylabel("contact i", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.85)

    # --- Panel B: Δ_rest ---
    ax = axes[1]
    im = ax.imshow(M_rest, cmap="RdBu_r", vmin=-vmax, vmax=+vmax,
                   aspect="equal", rasterized=True)
    ax.set_title(r"B — $\Delta_{\mathrm{rest}} = D_{\mathrm{post}} - D_{\mathrm{pre}}$",
                 fontsize=11)
    ax.set_xlabel("contact j", fontsize=9)
    ax.set_yticklabels([])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.85)

    # --- Panel C: σ(i,j) ---
    ax = axes[2]
    sigma_cmap = ListedColormap(["#7e3a8a",   # -1: opposite (purple)
                                 "#cccccc",   #  0: zero (grey)
                                 "#f4c430"])  # +1: same (yellow)
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], sigma_cmap.N)
    im = ax.imshow(M_sigma, cmap=sigma_cmap, norm=norm,
                   aspect="equal", rasterized=True)
    ax.set_title(r"C — $\sigma = \mathrm{sgn}(\Delta_{\mathrm{task}}) \cdot "
                 r"\mathrm{sgn}(\Delta_{\mathrm{rest}})$", fontsize=11)
    ax.set_xlabel("contact j", fontsize=9)
    ax.set_yticklabels([])
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.85,
                      ticks=[-1, 0, +1])
    cb.set_ticklabels(["−1 (opp)", "0", "+1 (same)"])
    ax.text(0.02, 0.98,
            f"σ>0: {frac_pos:.2f}\nσ<0: {frac_neg:.2f}",
            transform=ax.transAxes, fontsize=9, va="top", ha="left",
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"))

    # --- Panel D: hexbin scatter ---
    ax = axes[3]
    ax.axhline(0, color="#888", lw=0.6, zorder=1)
    ax.axvline(0, color="#888", lw=0.6, zorder=1)
    lim = float(max(np.abs(dD_task).max(), np.abs(dD_rest).max()))
    ax.plot([-lim, lim], [-lim, lim], color="#444", lw=0.8, ls="--",
            zorder=1)
    hb = ax.hexbin(dD_task, dD_rest, gridsize=40, cmap="viridis",
                   bins="log", mincnt=1, rasterized=True)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel(r"$\Delta_{\mathrm{task}}(i,j)$", fontsize=10)
    ax.set_ylabel(r"$\Delta_{\mathrm{rest}}(i,j)$", fontsize=10)
    ax.set_title(rf"D — pair scatter, $\rho = {rho:+.3f}$ "
                 rf"on $m = {m}$ pairs", fontsize=11)
    ax.set_aspect("equal")
    fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.04, shrink=0.85,
                 label="log density")

    fig.tight_layout()
    out = OUT_DIR / f"{patient}_{band}.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default="Pat_06")
    ap.add_argument("--band", default="alpha")
    ap.add_argument("--all", action="store_true",
                    help="Render every (patient, band) cell.")
    args = ap.parse_args()

    if args.all:
        summary = pd.read_csv(IN_DIR / "per_cell_summary.csv")
        for _, r in summary.iterrows():
            render_one(str(r["patient"]), str(r["band"]))
    else:
        render_one(args.patient, args.band)


if __name__ == "__main__":
    main()
