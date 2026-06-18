#!/usr/bin/env python3
"""Per-pair cophenetic co-movement as a density — the tie-free way.

The earlier joint-density renders failed because ``Δ^coph`` is ~90% tied: the
non-mover atom piles on the rank-centre (= the diagonal) and fakes a trace in
every band. The fix is to **actually remove the non-movers** (not cancel them):

  1. per patient, drop the modal (non-mover) atom in each axis — keep only the
     pairs that genuinely moved (~82–85% of pairs);
  2. RE-RANK the survivors with random tie-breaking → clean uniform marginals;
  3. pool the within-patient ranks over the cohort and KDE them;
  4. plot the empirical-copula excess (density − 1), hyperzoomed.

The diagonal ridge that appears IS the per-pair rank co-movement the §5.3
Spearman ρ measures — shown directly, with the tie artifact gone. α, β and γ_l
show a clear diagonal; θ/δ/γ_h are flat.

SCOPE — read before using this figure. This is a POOLED density: it shows the
*magnitude* of co-movement, and it is dominated by the patients that carry it.
It therefore shows that γ_l's per-pair co-movement is as strong as α/β's — which
is TRUE at the per-pair level (Pat_05 ρ=0.85, Pat_02 0.76, Pat_06 0.59). It does
NOT and cannot deliver the cohort verdict: a density integrates out the patient
index, so it cannot express cohort COHERENCE (α/β 8–9/10 coherent; γ_l a 6/4
split). The verdict that separates α/β (trace) from γ_l (sub-threshold) lives in
the per-patient gate — the forest figure
(`preprint_23_bands_perpair_rho_null_forest`). Use the two together: this
density justifies that the correlation is a real diagonal; the forest carries
the cohort verdict.

Output: data/preprint/figures/all_bands/fig_bands_perpair_movers_density_coph.pdf
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
sys.path.insert(0, str(Path(__file__).resolve().parent))
p18 = importlib.import_module("preprint_18_bands_joint_density_rawfc")
use_lrg_style()

COHORT = p18.COHORT
BAND_ORDER = p18.BAND_ORDER
TRACE_CMAP = p18.TRACE_CMAP

GRID = 120
KDE_BW = 0.13
SUBSAMPLE = 24000          # per band, for KDE speed (pooled n ≈ 56k)
SEED = 20260605


def movers_uv(dt: np.ndarray, dr: np.ndarray, rng) -> tuple[np.ndarray, np.ndarray]:
    """Drop the non-mover (modal) atom in each axis, re-rank survivors → (u, v)."""
    dt = np.round(dt, 12)
    dr = np.round(dr, 12)
    vt, ct = np.unique(dt, return_counts=True)
    vr, cr = np.unique(dr, return_counts=True)
    mask = (dt != vt[ct.argmax()]) & (dr != vr[cr.argmax()])
    return (p18.rank_random_tiebreak(dt[mask], rng),
            p18.rank_random_tiebreak(dr[mask], rng))


def pooled_excess(band: str, rng) -> tuple[np.ndarray, float, float]:
    pairs = p18.coph_obs_pairs(band)
    U, V, fmov = [], [], []
    for dt, dr in pairs.values():
        u, v = movers_uv(dt, dr, rng)
        U.append(u); V.append(v); fmov.append(u.size / dt.size)
    U = np.concatenate(U); V = np.concatenate(V)
    if U.size > SUBSAMPLE:
        sel = rng.choice(U.size, SUBSAMPLE, replace=False)
        U, V = U[sel], V[sel]
    xy = np.linspace(0.0, 1.0, GRID)
    xx, yy = np.meshgrid(xy, xy)
    dens = gaussian_kde(np.vstack([U, V]), bw_method=KDE_BW)(
        np.vstack([xx.ravel(), yy.ravel()])).reshape(GRID, GRID)
    excess = dens - 1.0                                   # [v, u]
    c = (np.arange(GRID) + 0.5) / GRID
    UU, VV = np.meshgrid(c, c)
    on = np.abs(UU - VV) < 0.1
    diag = float(np.maximum(excess[on], 0).mean())
    return excess, float(np.mean(fmov)), diag


def main() -> Path:
    rng = np.random.default_rng(SEED)
    fields, diags = {}, {}
    print("Per-pair cophenetic co-movement (movers-removed, pooled empirical copula):")
    for band in BAND_ORDER:
        fields[band], fmov, diag = pooled_excess(band, rng)
        diags[band] = diag
        print(f"  {band:11s} movers={fmov:.2f}  mean diagonal excess (|u-v|<0.1) = {diag:+.3f}")
    vabs = float(np.percentile(
        np.concatenate([np.abs(f).ravel() for f in fields.values()]), 98.0))
    print(f"  colour limit (p98 |density-1|) = {vabs:.3f}")

    fig = plt.figure(figsize=(15.0, 10.0))
    outer = fig.add_gridspec(2, 3, wspace=0.16, hspace=0.20,
                             left=0.05, right=0.91, top=0.95, bottom=0.08)
    last = None
    for idx, band in enumerate(BAND_ORDER):
        rr, cc = idx // 3, idx % 3
        ax = fig.add_subplot(outer[rr, cc])
        last = ax.imshow(fields[band], origin="lower", extent=[0, 1, 0, 1],
                         cmap=TRACE_CMAP, vmin=-vabs, vmax=vabs,
                         aspect="auto", interpolation="bilinear")
        ax.plot([0, 1], [0, 1], color="white", lw=1.5, alpha=0.7, zorder=3)
        ax.plot([0, 1], [1, 0], color="white", lw=0.8, ls="--", alpha=0.35,
                zorder=3)
        ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "½", "1"]); ax.set_yticklabels(["0", "½", "1"])
        ax.set_box_aspect(1.0)
        for s in ax.spines.values():
            s.set_edgecolor("0.55"); s.set_linewidth(0.9)
        tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
        ax.text(0.06, 0.94, tex, transform=ax.transAxes, ha="left", va="top",
                fontsize=24, fontweight="bold", color="0.18")
        if cc == 0:
            ax.set_ylabel(r"$v$ = within-patient rank of "
                          r"$\Delta^{\mathrm{coph}}_{\mathrm{rest}}$ (movers)",
                          color="0.20")
        if rr == 1:
            ax.set_xlabel(r"$u$ = within-patient rank of "
                          r"$\Delta^{\mathrm{coph}}_{\mathrm{task}}$ (movers)",
                          color="0.20")

    cbar_ax = fig.add_axes([0.925, 0.18, 0.014, 0.60])
    cb = fig.colorbar(last, cax=cbar_ax, orientation="vertical", extend="both")
    cb.set_label(r"empirical-copula excess  (density $-$ 1)   "
                 r"(green diagonal = per-pair co-movement)",
                 rotation=270, labelpad=22)

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_perpair_movers_density_coph.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
