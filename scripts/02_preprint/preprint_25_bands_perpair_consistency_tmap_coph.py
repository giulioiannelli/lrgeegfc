#!/usr/bin/env python3
"""Per-pair cophenetic trace as a CONSISTENCY statistical map (not a magnitude
density). Resolves the two pathologies of a plain co-movement density:

  • ties      → removed (drop the non-mover atom, re-rank the movers);
  • outliers  → killed by a per-cell cohort t-statistic (÷ between-patient SE),
                so one ρ=0.9 patient cannot light a cell — only CONSISTENT
                excess across the cohort does;
  • too-fine  → coarse cells (default 6×6, ~200 movers/cell/patient) so β's
                distributed, moderate correlation clears per-cell noise instead
                of being buried by it.

Construction (per band, per cell of a GRID×GRID grid on the movers' rank square):
  per patient p:  E_p = (movers density, obs) − (movers density, matched-strength
                  surrogate median)             — excess over the patient's OWN null
  cohort map:     t(cell) = mean_p E_p / (std_p E_p / sqrt(n))   (n = 10)
Cells with t > t_{0.05, df=9} = 1.833 (one-sided) are marked: a CONSISTENT
cohort excess over the matched-strength null — i.e. the gate, resolved in the
joint-rank plane. Green diagonal = per-pair trace.

What it shows (cohort, n=10): α, β, γ_l carry a significant diagonal; γ_h, θ, δ
are flat. β sits firmly in the trace group, well above γ_h — because the
t-statistic rewards β's consistent-moderate effect and discounts the
few-extreme-patient bands (γ_h/γ_l) that a magnitude density over-weights.
The exact β-vs-γ_l order is at the per-pair noise level (they are genuinely
comparable per-pair); the cohort VERDICT that separates them (β trace, γ_l
sub-threshold) is the paired-Wilcoxon gate / forest figure
(`preprint_23_bands_perpair_rho_null_forest`). This map shows WHERE and HOW
CONSISTENTLY the co-movement sits; the forest carries the per-band verdict.

Output: data/preprint/figures/all_bands/fig_bands_perpair_consistency_tmap_coph.pdf
Cache:  data/cache/perpair_consistency_tmap/<band>_G{G}_NS{NS}.npz
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import zoom
from scipy.stats import t as student_t

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
sys.path.insert(0, str(Path(__file__).resolve().parent))
p18 = importlib.import_module("preprint_18_bands_joint_density_rawfc")
use_lrg_style()

COHORT = p18.COHORT
BAND_ORDER = p18.BAND_ORDER
TRACE_CMAP = p18.TRACE_CMAP

GRID = 6                 # coarse: ~200 movers/cell/patient → distributed signal clears noise
N_SURR = 100             # surrogate realizations for the per-patient null median
SEED = 20260605
TMAP_CACHE = CACHE_ROOT / "perpair_consistency_tmap"
T_CRIT = float(student_t.ppf(0.95, df=len(COHORT) - 1))   # one-sided p<0.05, df=9 → 1.833
VMAX = 4.0


def _movers_uv(dt, dr, rng):
    dt = np.round(dt, 12); dr = np.round(dr, 12)
    vt, ct = np.unique(dt, return_counts=True)
    vr, cr = np.unique(dr, return_counts=True)
    k = (dt != vt[ct.argmax()]) & (dr != vr[cr.argmax()])
    return (p18.rank_random_tiebreak(dt[k], rng),
            p18.rank_random_tiebreak(dr[k], rng))


def _density(u, v, edges):
    H, _, _ = np.histogram2d(u, v, bins=[edges, edges], density=True)
    return H


def band_tmap(band: str) -> np.ndarray:
    """Per-cell cohort t of (obs − own matched-strength-null) movers density."""
    cache = TMAP_CACHE / f"{band}_G{GRID}_NS{N_SURR}_seed{SEED}.npz"
    if cache.exists():
        with np.load(cache) as d:
            return d["tmap"]
    rng = np.random.default_rng(SEED)
    edges = np.linspace(0.0, 1.0, GRID + 1)
    pairs = p18.coph_obs_pairs(band)
    E = []
    for pat in COHORT:
        dt, dr = pairs[pat]
        obs = _density(*_movers_uv(dt, dr, rng), edges)
        dts, drs, _ = p18._surrogate_rho_and_pairs(pat, band, "coph")
        surr = [_density(*_movers_uv(s_dt, s_dr, rng), edges)
                for s_dt, s_dr in list(zip(dts, drs))[:N_SURR]]
        null = np.median(np.stack(surr), axis=0)
        E.append(obs - null)
    E = np.stack(E)                                   # (n, G, G)
    tmap = E.mean(0) / (E.std(0, ddof=1) / np.sqrt(E.shape[0]) + 1e-12)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, tmap=tmap)
    return tmap


def main() -> Path:
    edges = np.linspace(0.0, 1.0, GRID + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    cc, rr = np.meshgrid(centres, centres, indexing="ij")    # [iu, iv]
    on_diag = np.abs(cc - rr) < (1.2 / GRID)

    tmaps = {}
    print(f"Per-pair consistency t-map (GRID={GRID}, n={len(COHORT)}, "
          f"t_crit={T_CRIT:.2f}); diagonal summary:")
    for band in BAND_ORDER:
        t = band_tmap(band)
        tmaps[band] = t
        nsig = int((t[on_diag] > T_CRIT).sum())
        print(f"  {band:11s} mean diag t = {t[on_diag].mean():+.2f}  "
              f"sig diagonal cells (t>{T_CRIT:.2f}) = {nsig}/{int(on_diag.sum())}")

    fig = plt.figure(figsize=(15.0, 10.0))
    outer = fig.add_gridspec(2, 3, wspace=0.16, hspace=0.20,
                             left=0.05, right=0.91, top=0.95, bottom=0.08)
    last = None
    for idx, band in enumerate(BAND_ORDER):
        ax = fig.add_subplot(outer[idx // 3, idx % 3])
        t = tmaps[band]
        # smooth, vector display of the coarse-cell statistics: Gouraud shading
        # between the GRID×GRID cell-centre t-values (the stats stay coarse for
        # stability; only the rendering is interpolated).
        last = ax.pcolormesh(centres, centres, t.T, shading="gouraud",
                             cmap=TRACE_CMAP, vmin=-VMAX, vmax=VMAX, zorder=1)
        # smooth significance boundary: contour of t = t_crit (one-sided p<0.05),
        # interpolated from the coarse grid for a clean outline of the trace zone.
        tf = zoom(t.T, 10, order=1)
        gf = np.linspace(centres[0], centres[-1], tf.shape[0])
        xf, yf = np.meshgrid(gf, gf)
        ax.contour(xf, yf, tf, levels=[T_CRIT], colors=["white"], linewidths=2.2,
                   zorder=5)
        ax.contour(xf, yf, tf, levels=[T_CRIT], colors=["0.12"], linewidths=0.7,
                   zorder=5)
        ax.plot([0, 1], [0, 1], color="0.30", lw=1.0, ls="--", alpha=0.85, zorder=4)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "½", "1"]); ax.set_yticklabels(["0", "½", "1"])
        ax.set_box_aspect(1.0)
        for s in ax.spines.values():
            s.set_edgecolor("0.55"); s.set_linewidth(0.9)
        ax.text(0.06, 0.94, BRAIN_BAND_TEX_DICT.get(band, rf"${band}$"),
                transform=ax.transAxes, ha="left", va="top",
                fontsize=24, fontweight="bold", color="0.18")
        if idx % 3 == 0:
            ax.set_ylabel(r"$v$ = rank of $\Delta^{\mathrm{coph}}_{\mathrm{rest}}$"
                          r" (movers)", color="0.20")
        if idx // 3 == 1:
            ax.set_xlabel(r"$u$ = rank of $\Delta^{\mathrm{coph}}_{\mathrm{task}}$"
                          r" (movers)", color="0.20")

    cbar_ax = fig.add_axes([0.925, 0.18, 0.014, 0.60])
    cb = fig.colorbar(last, cax=cbar_ax, orientation="vertical", extend="both")
    cb.set_label(r"cohort consistency $t$  (per-pair excess over matched-strength"
                 r" null)" "\n" r"outlined zone = consistent cohort trace, $p<0.05$",
                 rotation=270, labelpad=30)

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_perpair_consistency_tmap_coph.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
