"""Per-leaf recovery-slope ``s_TR(c)`` on the rsPost dendrogram, β band.

For each β cohort-backbone patient, draw the rest_post ``D_coph`` dendrogram
(UPGMA linkage on the propagator distance) with each leaf coloured by its
per-contact recovery slope

    s_TR(c) = Σ_{j≠c} Δ_task(c,j)·Δ_rest(c,j) / Σ_{j≠c} Δ_task(c,j)²

— the row-c restriction of the global ``s_TR``, read as "the fraction of the
task-induced distance shift at contact c's row that is recovered in rsPost".
Diverging palette centred at 0 (anti — grey — trace). This is the literal
per-leaf decomposition of the cohort slope statistic; averaging the leaf
colours lands near (not exactly at) the global ``s_TR`` because of per-row
magnitude weighting.

Caveat shown by the figure: because ``s_TR`` is magnitude-weighted, leaves on
high-distance rows dominate, and the cohort-level slope fails the matched-
strength gate at β (audit_73). The per-leaf map is a descriptive decomposition,
not a verified localisation.

PDF only, vector. No suptitle. Dendrogram log-y limits per project rule.
Output: data/audit/pair_trace_measures_comparison/figures/per_leaf_slope_<band>/<pat>.pdf
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.notebook import move_to_rootf
from lrg_eegfc.visuals.styles import use_lrg_style

move_to_rootf(pathname="lrgeegfc")
use_lrg_style()

_AUDIT = Path("scripts/01_compute/audit/audit_73_pair_trace_measures_comparison.py")
_spec = importlib.util.spec_from_file_location("audit_73", _AUDIT)
a73 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a73)

FIG_ROOT = Path("data/audit/pair_trace_measures_comparison/figures")
BACKBONE = ["Pat_02", "Pat_05", "Pat_06", "Pat_08"]
CMAP = "coolwarm"  # saturated diverging, readable on white (no near-white wash)


def dcoph_square(eig):
    _, D = a73.primitives_from_eig(*eig)
    Z = linkage(squareform(D, checks=False), method="average")
    return a73.cophenet_matrix(Z, condensed=False), D


def per_leaf_slope(dT_sq, dR_sq):
    """s_TR(c) per row c of the full square Δ matrices (diagonal excluded)."""
    np.fill_diagonal(dT_sq, 0.0)
    np.fill_diagonal(dR_sq, 0.0)
    num = np.sum(dT_sq * dR_sq, axis=1)
    den = np.sum(dT_sq * dT_sq, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(den > 0, num / den, np.nan)
    return s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--band", default="beta")
    ap.add_argument("--patients", default=",".join(BACKBONE))
    args = ap.parse_args()
    band = args.band
    patients = [p.strip() for p in args.patients.split(",") if p.strip()]

    out_dir = FIG_ROOT / f"per_leaf_slope_{band}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # shared symmetric colour scale across patients for comparability
    vspan = 1.5  # s_TR(c) clipped to [-vspan, +vspan] for colour
    norm = TwoSlopeNorm(vmin=-vspan, vcenter=0.0, vmax=vspan)
    cmap = plt.get_cmap(CMAP)

    for pat in patients:
        eigs = a73.observed_eigs(pat, band)
        N = eigs["rest_pre_A"][1].shape[0]
        Dc = {}
        Dpost_raw = None
        for ph in a73.PHASES:
            coph, Draw = dcoph_square(eigs[ph])
            Dc[ph] = coph
            if ph == "rest_post":
                Dpost_raw = Draw

        dT = Dc["task_test"] - Dc["rest_pre_A"]
        dR = Dc["rest_post"] - Dc["rest_pre_B"]
        s_leaf = per_leaf_slope(dT.copy(), dR.copy())

        # rsPost dendrogram structure
        Zpost = linkage(squareform(Dpost_raw, checks=False), method="average")

        fig, ax = plt.subplots(figsize=(7.5, 3.4))
        R = dendrogram(Zpost, ax=ax, no_labels=True, color_threshold=0,
                       above_threshold_color="0.75")
        merge_h = np.sort(Zpost[:, 2])
        tmin = merge_h[0] * 0.8
        tmax = merge_h[-1] * 1.05
        ax.set_yscale("log")
        ax.set_ylim(tmin, tmax)

        # leaf -> parent merge height (for the coloured vline length)
        leaf_parent_h = {}
        for row in Zpost:
            a_, b_, h = int(row[0]), int(row[1]), float(row[2])
            if a_ < N:
                leaf_parent_h.setdefault(a_, h)
            if b_ < N:
                leaf_parent_h.setdefault(b_, h)

        for i, lf in enumerate(R["leaves"]):
            sc = s_leaf[lf]
            if not np.isfinite(sc):
                col = "0.6"
            else:
                col = cmap(norm(np.clip(sc, -vspan, vspan)))
            x = 5.0 + 10.0 * i
            h_parent = leaf_parent_h.get(lf, tmin * 1.3)
            ax.vlines(x, tmin, h_parent, color=col, lw=1.6, zorder=4)

        ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=7)
        for sn in ("top", "right", "bottom"):
            ax.spines[sn].set_visible(False)
        ax.set_ylabel(r"$D_{\mathrm{coph}}$ merge height (rsPost)", fontsize=8)

        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02,
                          extend="both")
        cb.set_label(r"per-leaf $s_{\mathrm{TR}}(c)$", fontsize=8)
        cb.ax.axhline(0.0, color="0.2", lw=0.8)

        fig.tight_layout()
        out = out_dir / f"{pat}.pdf"
        fig.savefig(out, transparent=True)
        plt.close(fig)
        med = np.nanmedian(s_leaf)
        print(f"Wrote {out}  (median s_TR(c)={med:+.3f}, N={N})")


if __name__ == "__main__":
    main()
