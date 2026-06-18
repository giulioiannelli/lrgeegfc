#!/usr/bin/env python3
r"""Dissociation map — per-pair cophenetic trace vs leading-subspace Grassmann.

A single 2-D plane, one point per band, that locates where the two LRG trace
probes agree and where they dissociate:

  * x = cophenetic ρ_split^coph effect, signed −log10(p_Wilcoxon) · sign(median).
    The novel per-pair, multiscale-aware backbone probe: does each pair keep
    its relative ultrametric (communication-distance) position task→rest?
  * y = Grassmann T_G* effect, −log10(cluster-mass permutation p). The
    "usual spectral analysis" cross-check: do the dominant collective
    eigenmodes (span{φ_2..φ_{k+1}}) rotate?

Significance thresholds (p = 0.05) are drawn on both axes, partitioning the
plane into four cells: NEITHER, COPH-ONLY (local trace, no global-mode trace),
GRASSMANN-ONLY (global trace, no local), BOTH (local + global).

γ_l is plotted TWICE: full-graph (Grassmann-only) and epi-excluded (the
cophenetic trace EMERGES under epileptic-node exclusion, p 0.116 → 0.024),
joined by an arrow — its dissociation is CONDITIONAL on the epileptic confound.

All numbers are surfaced from the locked matched-strength outputs (R = 200);
nothing is recomputed. Sources:
    cophenet   data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
    cophenet (epi-X)  data/audit/epi_stratified/cophenetic_cohort.csv
    Grassmann  data/audit/grassmann_cluster_extent/cohort_summary.csv
Cross-checked against .agents/preprint/locked/VERDICT_LEDGER.md.

No in-axes numeric text beyond axis labels + single-letter panel tag (project
rule); the full numeric dissociation table is printed to stdout.

Output: data/preprint/figures/all_bands/fig_bands_coph_grassmann_dissociation_map.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
ALPHA = 0.05
PMIN = 1e-3  # empirical-null floor 1/(R+1) for R=200 ⇒ -log10 = 2.305 ceiling

# Quadrant background tints (low alpha, white-safe saturated hues).
CELL_COLOR = {
    "BOTH": "#1f7a1f",
    "COPH-ONLY": "#2c5f8a",
    "GRASSMANN-ONLY": "#e67e22",
    "NEITHER": "#95a5a6",
}


def _signed_neglog10(p: float, sign: float) -> float:
    """Signed −log10(p), floored at the empirical null 1/(R+1)."""
    p = max(float(p), PMIN)
    return float(np.sign(sign)) * (-np.log10(p))


def main() -> None:
    coph_src = (ROOT / "data" / "audit"
                / "matched_strength_surrogate_split_baseline" / "cohort_summary.csv")
    grass_src = (ROOT / "data" / "audit"
                 / "grassmann_cluster_extent" / "cohort_summary.csv")
    epi_src = ROOT / "data" / "audit" / "epi_stratified" / "cophenetic_cohort.csv"
    for s in (coph_src, grass_src, epi_src):
        if not s.exists():
            raise SystemExit(f"[preprint_34] missing {s}")

    coph = pd.read_csv(coph_src).set_index("band")
    grass = pd.read_csv(grass_src).set_index("band")
    epi = pd.read_csv(epi_src)
    epi_full = epi[epi.config == "full"].set_index("band")
    epi_x = epi[epi.config == "exclude_epi"].set_index("band")

    thr = -np.log10(ALPHA)  # 1.301

    # ---- assemble per-band coordinates -------------------------------------
    rows = []
    for b in BANDS:
        cp = float(coph.loc[b, "paired_wilcoxon_p"])
        cmed = float(coph.loc[b, "obs_median_rho"])
        cratio = cmed / float(coph.loc[b, "surr_median_rho_median"])
        cn = str(coph.loc[b, "n_above_surrogate"])
        x = _signed_neglog10(cp, cmed)

        gp = float(grass.loc[b, "cluster_p_cluster_mass"])
        gloo = float(grass.loc[b, "cluster_p_mass_loo_max"])
        glr = int(grass.loc[b, "obs_longest_run"])
        gverd = str(grass.loc[b, "verdict_cluster_extent"])
        y = -np.log10(max(gp, PMIN))  # Grassmann is one-sided ≥0 by construction

        coph_sig = (cp < ALPHA) and (cmed > 0)
        # Grassmann strong-tier needs LOO < 0.05 too (Decision 12); for the
        # plane we use the cohort gate p_mass<0.05 as the axis threshold and
        # annotate LOO-fragility (δ) in the printout / marker edge.
        grass_sig = gp < ALPHA
        cell = (
            "BOTH" if coph_sig and grass_sig
            else "COPH-ONLY" if coph_sig and not grass_sig
            else "GRASSMANN-ONLY" if grass_sig and not coph_sig
            else "NEITHER"
        )
        loo_fragile = grass_sig and (gloo >= ALPHA)
        rows.append(dict(band=b, x=x, y=y, cp=cp, cmed=cmed, cratio=cratio,
                         cn=cn, gp=gp, gloo=gloo, glr=glr, gverd=gverd,
                         cell=cell, loo_fragile=loo_fragile))
    R = pd.DataFrame(rows).set_index("band")

    # γ_l epi-excluded cophenetic point (the conditional emergence)
    glx = epi_x.loc["low_gamma"]
    glf = epi_full.loc["low_gamma"]
    glx_p = float(glx["paired_wilcoxon_p"])
    glx_med = float(glx["obs_median"])
    glx_x = _signed_neglog10(glx_p, glx_med)
    gl_y = R.loc["low_gamma", "y"]  # Grassmann unchanged (full-graph gate)

    # ---- stdout dissociation table -----------------------------------------
    print("\n=== ρ^coph × Grassmann dissociation (n=10, matched-strength R=200) ===")
    hdr = (f"{'band':<11}{'coph_p':>9}{'ratio':>9}{'n/10':>7}{'coph?':>7}"
           f"{'grass_p':>9}{'grLOO':>8}{'LR':>4}{'grass?':>8}  cell")
    print(hdr)
    for b in BANDS:
        r = R.loc[b]
        cflag = "TRACE" if (r.cp < ALPHA and r.cmed > 0) else "—"
        gflag = ("weak*" if r.loo_fragile else "TRACE") if r.gp < ALPHA else "—"
        print(f"{b:<11}{r.cp:>9.4f}{r.cratio:>9.2f}{r.cn:>7}{cflag:>7}"
              f"{r.gp:>9.4f}{r.gloo:>8.3f}{int(r.glr):>4}{gflag:>8}  {r.cell}")
    print(f"\nlow_gamma  exclude_epi cophenetic: p={glx_p:.4f} "
          f"(full {glf['paired_wilcoxon_p']:.4f}), median {glx_med:+.3f}, "
          f"8/10 → cophenetic trace EMERGES; Grassmann CONTRACTS under epi-X "
          f"(ledger C5: mass 66→33, LR 13→10). Conditional BOTH (epi-X) "
          f"vs Grassmann-only (full).")
    print("  * δ Grassmann: cohort gate p_mass=0.005 but LOO max=0.055 "
          "(Pat_08) → weak (Decision 12).")

    # ---- figure ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.4, 5.0))

    xmax = max(2.6, float(np.nanmax(np.abs(R["x"].values))) + 0.5, abs(glx_x) + 0.5)
    ymax = max(2.6, float(np.nanmax(R["y"].values)) + 0.4)
    xlo = -1.4

    # quadrant background tints (only the four meaningful corners)
    ax.add_patch(plt.Rectangle((thr, thr), xmax - thr, ymax - thr,
                               facecolor=CELL_COLOR["BOTH"], alpha=0.08, zorder=0))
    ax.add_patch(plt.Rectangle((thr, 0), xmax - thr, thr,
                               facecolor=CELL_COLOR["COPH-ONLY"], alpha=0.08, zorder=0))
    ax.add_patch(plt.Rectangle((xlo, thr), thr - xlo, ymax - thr,
                               facecolor=CELL_COLOR["GRASSMANN-ONLY"], alpha=0.08, zorder=0))
    ax.add_patch(plt.Rectangle((xlo, 0), thr - xlo, thr,
                               facecolor=CELL_COLOR["NEITHER"], alpha=0.08, zorder=0))

    # significance guides
    ax.axhline(thr, color="0.45", lw=0.9, ls="--", zorder=1)
    ax.axvline(thr, color="0.45", lw=0.9, ls="--", zorder=1)
    ax.axvline(0.0, color="0.7", lw=0.7, ls=":", zorder=1)

    # quadrant labels (cell names — structural labels, not stats)
    qkw = dict(ha="center", va="center", fontsize=8, color="0.35", zorder=1)
    ax.text((thr + xmax) / 2, (thr + ymax) / 2, "both\n(local + global)", **qkw)
    ax.text((thr + xmax) / 2, thr / 2, "local-only\n(coph)", **qkw)
    ax.text((xlo + thr) / 2, (thr + ymax) / 2, "global-only\n(Grassmann)", **qkw)
    ax.text((xlo + thr) / 2, thr / 2, "neither", **qkw)

    # band points coloured by cell
    for b in BANDS:
        r = R.loc[b]
        col = CELL_COLOR[r.cell]
        edge = "crimson" if r.loo_fragile else "white"
        ax.scatter(r.x, r.y, s=190, c=col, edgecolors=edge,
                   linewidths=1.6, zorder=4)
        # band glyph just outside the marker
        ax.annotate(BRAIN_BAND_TEX_DICT.get(b, b), (r.x, r.y),
                    textcoords="offset points", xytext=(11, 6),
                    fontsize=11, zorder=5)

    # γ_l conditional emergence: arrow from full → epi-excluded cophenetic
    ax.annotate(
        "", xy=(glx_x, gl_y), xytext=(R.loc["low_gamma", "x"], gl_y),
        arrowprops=dict(arrowstyle="-|>", color=CELL_COLOR["GRASSMANN-ONLY"],
                        lw=1.4, ls=(0, (4, 2))), zorder=3)
    ax.scatter(glx_x, gl_y, s=150, marker="D",
               c=CELL_COLOR["BOTH"], edgecolors="white", linewidths=1.4, zorder=4)
    ax.annotate(BRAIN_BAND_TEX_DICT["low_gamma"] + r"$^{\mathrm{exc.epi}}$",
                (glx_x, gl_y), textcoords="offset points", xytext=(6, -16),
                fontsize=9, color=CELL_COLOR["BOTH"], zorder=5)

    ax.set_xlim(xlo, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel(r"per-pair cophenetic effect  "
                  r"$\mathrm{sign}(\tilde\rho^{\,\mathrm{coph}})\cdot"
                  r"(-\log_{10}p_{\mathrm{Wilcoxon}})$")
    ax.set_ylabel(r"Grassmann subspace effect  "
                  r"$-\log_{10}p_{\mathrm{mass}}$")
    ax.text(0.02, 0.97, "a", transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top")
    fig.tight_layout()

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_bands_coph_grassmann_dissociation_map.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"\n[preprint_34] -> {out}")


if __name__ == "__main__":
    main()
