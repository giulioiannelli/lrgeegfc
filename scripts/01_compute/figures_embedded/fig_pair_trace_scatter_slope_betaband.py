"""Canonical per-pair scatter for the asymmetric slope `s_TR` at the β band.

For each of the four β cohort-backbone patients (Pat_02, Pat_05, Pat_06,
Pat_08), three panels on the cophenet primitive ``D_coph``:
  A — Δ_task per-pair shift heatmap (D_coph(task) − D_coph(rsPre_A))
  B — Δ_rest per-pair shift heatmap (D_coph(rsPost) − D_coph(rsPre_B))
  C — scatter Δ_task (x) vs Δ_rest (y), one dot per upper-triangle pair, with
      the through-origin regression line of slope ``s_TR`` and the y=x "full
      trace" reference. Annotated with s_TR, R², ρ_Spearman, ρ_Pearson.

The scatter is the point of the figure: on ``D_coph`` the cloud is dominated by
a handful of large-magnitude merge-height pairs, so the magnitude-weighted
slope ``s_TR`` is leveraged by those few pairs (and reproduced by matched-
strength surrogates), whereas the rank correlation weights all pairs equally.
This is the visual mechanism behind the audit_73 verdict that ``s_TR`` loses
the β trace that ``ρ_Spearman`` keeps.

Band kwarg: ``--band beta`` (default) extends to other bands.

PDF only, vector. No suptitle. ``--band`` other than beta is supported.
Output: data/audit/pair_trace_measures_comparison/figures/scatter_slope_<band>/<pat>.pdf
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.notebook import move_to_rootf
from lrg_eegfc.utils.metrics.hypothesis import regression_slope_through_origin
from lrg_eegfc.visuals.styles import use_lrg_style
from scipy.stats import pearsonr, spearmanr

move_to_rootf(pathname="lrgeegfc")
use_lrg_style()

# Reuse the audit_73 compute helpers (single source of truth for the primitives)
_AUDIT = Path("scripts/01_compute/audit/audit_73_pair_trace_measures_comparison.py")
_spec = importlib.util.spec_from_file_location("audit_73", _AUDIT)
a73 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a73)

FIG_ROOT = Path("data/audit/pair_trace_measures_comparison/figures")
BACKBONE = ["Pat_02", "Pat_05", "Pat_06", "Pat_08"]
DIVERGING = "coolwarm"  # saturated endpoints, used on symmetric-limit heatmaps


def dcoph_square(eig):
    """Full N×N cophenet matrix for one phase eigendecomposition."""
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform
    _, D = a73.primitives_from_eig(*eig)
    Z = linkage(squareform(D, checks=False), method="average")
    return a73.cophenet_matrix(Z, condensed=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--band", default="beta")
    ap.add_argument("--patients", default=",".join(BACKBONE))
    args = ap.parse_args()
    band = args.band
    patients = [p.strip() for p in args.patients.split(",") if p.strip()]

    out_dir = FIG_ROOT / f"scatter_slope_{band}"
    out_dir.mkdir(parents=True, exist_ok=True)

    for pat in patients:
        eigs = a73.observed_eigs(pat, band)
        N = eigs["rest_pre_A"][1].shape[0]
        iu = np.triu_indices(N, k=1)
        Dc = {ph: dcoph_square(eigs[ph]) for ph in a73.PHASES}
        dT = Dc["task_test"][iu] - Dc["rest_pre_A"][iu]
        dR = Dc["rest_post"][iu] - Dc["rest_pre_B"][iu]

        slope, rsq = regression_slope_through_origin(dT, dR)
        rho_s = float(spearmanr(dT, dR).statistic)
        rho_p = float(pearsonr(dT, dR).statistic)

        dT_sq = Dc["task_test"] - Dc["rest_pre_A"]
        dR_sq = Dc["rest_post"] - Dc["rest_pre_B"]
        vmax = max(np.abs(dT_sq).max(), np.abs(dR_sq).max())

        fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.1))

        for ax, M, lab in (
            (axes[0], dT_sq, r"$\Delta_{\mathrm{task}}=D_{\mathrm{coph}}^{\,task}-D_{\mathrm{coph}}^{\,rsPre_A}$"),
            (axes[1], dR_sq, r"$\Delta_{\mathrm{rest}}=D_{\mathrm{coph}}^{\,rsPost}-D_{\mathrm{coph}}^{\,rsPre_B}$"),
        ):
            Mm = M.copy()
            np.fill_diagonal(Mm, np.nan)
            im = ax.imshow(Mm, cmap=DIVERGING, vmin=-vmax, vmax=vmax)
            ax.set_xlabel(lab, fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        ax = axes[2]
        ax.axhline(0, color="0.7", lw=0.6)
        ax.axvline(0, color="0.7", lw=0.6)
        lim = max(np.abs(dT).max(), np.abs(dR).max()) * 1.05
        # y = x full-trace reference
        ax.plot([-lim, lim], [-lim, lim], color="0.6", ls=":", lw=1.0,
                label=r"$y=x$ (full trace)")
        # through-origin s_TR fit
        ax.plot([-lim, lim], [-slope * lim, slope * lim], color="#d62728",
                lw=1.8, label=rf"$s_{{\mathrm{{TR}}}}={slope:+.2f}$")
        ax.scatter(dT, dR, s=8, alpha=0.35, color="#0b5394",
                   edgecolors="none", zorder=2)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$\Delta_{\mathrm{task}}$ (per pair)", fontsize=9)
        ax.set_ylabel(r"$\Delta_{\mathrm{rest}}$ (per pair)", fontsize=9)
        ax.legend(loc="upper left", fontsize=7, frameon=False)
        ax.text(0.97, 0.06,
                f"$R^2$={rsq:.2f}\n"
                rf"$\rho_{{\mathrm{{Spear}}}}$={rho_s:+.2f}"
                "\n"
                rf"$\rho_{{\mathrm{{Pear}}}}$={rho_p:+.2f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5)

        fig.tight_layout()
        out = out_dir / f"{pat}.pdf"
        fig.savefig(out, transparent=True)
        plt.close(fig)
        print(f"Wrote {out}  (s_TR={slope:+.3f} rho_S={rho_s:+.3f} R2={rsq:.3f})")


if __name__ == "__main__":
    main()
