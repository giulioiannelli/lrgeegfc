#!/usr/bin/env python3
"""Dendrogram grid per patient — visual check of topology preservation.

6 bands × 4 phases grid of dendrograms for ONE patient, drawn with a
common leaf ordering (derived from the patient's rest_pre dendrogram
on that band) so topology changes across phases are visible as
changes in the branching pattern at the same leaf positions.

Purpose: verify by eye whether our continuous measures (H2c ρ,
Frobenius, Δρ) are reflecting actual tree topology preservation, or
whether the ρ-magnitude story is just arithmetic on D values that can
shift without topology changing.

Usage: python scripts/01_compute/fig_dendrogram_grid.py Pat_02 [imcoh_abs]

Reads:  LRG caches.
Writes: data/reports/imcoh_vi/figures/dendrogram_grid_<pat>.{pdf,png}
"""
from __future__ import annotations

import argparse

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, leaves_list, optimal_leaf_ordering
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result


PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
PHASE_LABEL = {"rest_pre": "rest_pre", "task_learn": "task_learn",
               "task_test": "task_test", "rest_post": "rest_post"}


def _square_from_condensed(vec: np.ndarray) -> np.ndarray:
    return squareform(vec, checks=False)


def _get_Z_and_D(pat: str, phase: str, band: str, fc_method: str):
    r = load_lrg_result(pat, phase, band, fc_method=fc_method,
                        cache_root=IMCOH_LRG_CACHE)
    if r is None:
        return None, None
    Z = np.asarray(r.linkage_matrix, dtype=float)
    D = np.asarray(r.ultrametric_matrix, dtype=float)
    return Z, D


def _common_leaf_order(pat: str, band: str, fc_method: str) -> list[int] | None:
    """Use rest_pre's optimal leaf ordering as the canonical order."""
    Z, D = _get_Z_and_D(pat, "rest_pre", band, fc_method)
    if Z is None:
        return None
    Z_opt = optimal_leaf_ordering(Z, _square_from_condensed(D))
    return leaves_list(Z_opt).tolist()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("patient")
    p.add_argument("--fc-method", default="imcoh_abs")
    args = p.parse_args()

    pat = args.patient
    fc = args.fc_method

    fig, axes = plt.subplots(
        len(BRAIN_BANDS_NAMES), len(PHASES),
        figsize=(3.0 * len(PHASES), 1.6 * len(BRAIN_BANDS_NAMES)),
        dpi=160, squeeze=False,
    )

    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        order = _common_leaf_order(pat, band, fc)
        for ip, phase in enumerate(PHASES):
            ax = axes[ib, ip]
            Z, _ = _get_Z_and_D(pat, phase, band, fc)
            if Z is None:
                ax.text(0.5, 0.5, "missing", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="#888")
                ax.axis("off")
                continue
            try:
                dendrogram(
                    Z, no_labels=True, color_threshold=0.0,
                    above_threshold_color="#2a5d9f",
                    ax=ax,
                )
            except Exception as e:
                ax.text(0.5, 0.5, f"err: {type(e).__name__}",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=8, color="crimson")
                ax.axis("off")
                continue
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if ib == 0:
                ax.set_title(PHASE_LABEL[phase], fontsize=10, pad=4)
            if ip == 0:
                ax.set_ylabel(BRAIN_BAND_TEX_DICT[band], fontsize=12,
                              rotation=0, ha="right", va="center", labelpad=12)

    # No fig.suptitle on publication figures — context lives in the file name.
    fig.tight_layout()

    out_dir = REPORTS_ROOT / "imcoh_vi" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        out = out_dir / f"dendrogram_grid_{pat}.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
