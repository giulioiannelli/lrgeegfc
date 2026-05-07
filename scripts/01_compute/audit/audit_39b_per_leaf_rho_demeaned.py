#!/usr/bin/env python3
"""Audit 39b — Within-patient demeaned per-leaf rho + dendrogram link coloring.

Fixes the methodological issue in audit_39: per-leaf rho_ell of a high-overall-
rho patient produces ~uniformly positive per-leaf rho values because the
patient's global rank correlation dominates the per-leaf marginals. This is
NOT localization -- it is the patient being a global outlier whose signal
propagates to every leaf.

The proper localization measure demeans per-leaf rho_ell relative to the
patient-band's own mean:

    rho_demeaned(ell) = rho_ell(ell) - mean_{ell' in p,b}( rho_ell(ell') )
    rho_z(ell)        = (rho_ell(ell) - mean) / std        (within p,b)

A leaf with rho_demeaned > 0 has its neighborhood rotating in the trace
direction MORE than the patient's average -- that is genuine localization.
A leaf with rho_demeaned < 0 is reorganizing LESS than the patient's average
(or anti-trace within a globally-trace patient).

Visualization: per-(patient, band) dendrogram with BOTH links AND leaves
colored by rho_demeaned. A contiguous warm-colored branch indicates a
coherent localized trace-direction subtree above the patient's baseline.

Outputs
-------
``data/audit/per_leaf_rho_demeaned/leaf_rho_demeaned.csv``
``data/audit/per_leaf_rho_demeaned/cohort_summary.csv``
``data/outputs/figures/section_5_lrg_trace/per_leaf_rho_demeaned/{Pat_XX}_{band}.pdf``
``data/outputs/figures/section_5_lrg_trace/per_leaf_rho_demeaned/cohort_distribution.pdf``
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import dendrogram, leaves_list

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import _load_Z, BAND_ORDER, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE  # noqa: E402

LEAF_RHO_CSV = ROOT / "data" / "audit" / "per_leaf_sigma_aggregate" / "leaf_rho.csv"
OUT_DIR = ROOT / "data" / "audit" / "per_leaf_rho_demeaned"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "per_leaf_rho_demeaned"


def link_descendant_means(Z: np.ndarray, leaf_values: np.ndarray) -> dict[int, float]:
    """Mean of leaf_values across descendants of every internal node."""
    n_leaves = Z.shape[0] + 1
    desc_sum = np.zeros(Z.shape[0])
    desc_n = np.zeros(Z.shape[0], dtype=int)
    for r in range(Z.shape[0]):
        s = 0.0
        c = 0
        for col in (0, 1):
            child = int(Z[r, col])
            if child < n_leaves:
                s += leaf_values[child]
                c += 1
            else:
                idx = child - n_leaves
                s += desc_sum[idx]
                c += desc_n[idx]
        desc_sum[r] = s
        desc_n[r] = c
    return {n_leaves + r: desc_sum[r] / desc_n[r] for r in range(Z.shape[0])}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(LEAF_RHO_CSV)
    df["rho_mean_pb"] = df.groupby(["patient", "band"])["rho_leaf"].transform("mean")
    df["rho_std_pb"] = df.groupby(["patient", "band"])["rho_leaf"].transform("std")
    df["rho_demeaned"] = df["rho_leaf"] - df["rho_mean_pb"]
    df["rho_z"] = df["rho_demeaned"] / df["rho_std_pb"].replace(0, 1)
    df["is_local_peak_demeaned_03"] = (df["rho_demeaned"] >= 0.3).astype(int)
    df["is_local_peak_z_15"] = (df["rho_z"] >= 1.5).astype(int)
    df.to_csv(OUT_DIR / "leaf_rho_demeaned.csv", index=False)

    # Per-(patient, band) summary
    summary = (
        df.groupby(["patient", "band"])
        .agg(
            n_leaves=("leaf_id", "count"),
            rho_mean=("rho_leaf", "mean"),
            rho_std=("rho_leaf", "std"),
            n_local_peaks_demeaned=("is_local_peak_demeaned_03", "sum"),
            n_local_peaks_z=("is_local_peak_z_15", "sum"),
        )
        .reset_index()
    )
    summary["frac_local_peaks_demeaned"] = summary["n_local_peaks_demeaned"] / summary["n_leaves"]
    summary["frac_local_peaks_z"] = summary["n_local_peaks_z"] / summary["n_leaves"]
    summary.to_csv(OUT_DIR / "cohort_summary.csv", index=False)

    band_summary = (
        summary.groupby("band")
        .agg(
            mean_rho_pb=("rho_mean", "mean"),
            mean_n_local_peaks_demeaned=("n_local_peaks_demeaned", "mean"),
            mean_n_local_peaks_z=("n_local_peaks_z", "mean"),
            mean_frac_peaks_demeaned=("frac_local_peaks_demeaned", "mean"),
            n_pat_with_any_peak=("n_local_peaks_demeaned", lambda s: int((s > 0).sum())),
        )
        .reset_index()
    )
    band_summary.to_csv(OUT_DIR / "cohort_band_summary.csv", index=False)
    print(band_summary.to_string(index=False))

    # Per-(patient, band) figures: dendrogram with link + leaf coloring by rho_demeaned
    cmap = plt.get_cmap("RdBu_r")
    for pat in PATIENTS_4PHASE:
        for band in BAND_ORDER:
            sub = df[(df["patient"] == pat) & (df["band"] == band)]
            if sub.empty:
                continue
            Z = _load_Z(pat, "rest_post", band)
            if Z is None:
                continue
            rho = sub.sort_values("leaf_id")["rho_leaf"].values.astype(float)
            rho_dm = sub.sort_values("leaf_id")["rho_demeaned"].values.astype(float)
            n_leaves = Z.shape[0] + 1
            link_means = link_descendant_means(Z, rho_dm)
            vmax = max(0.20, float(np.nanmax(np.abs(rho_dm))))
            norm = Normalize(vmin=-vmax, vmax=vmax)

            def link_color(nid: int) -> str:
                return to_hex(cmap(norm(link_means.get(nid, 0.0))))

            order = leaves_list(Z)
            rho_ord = rho[order]
            rho_dm_ord = rho_dm[order]
            fig, axes = plt.subplots(3, 1, figsize=(11, 5.5),
                                     gridspec_kw={"height_ratios": [3.5, 0.35, 0.35], "hspace": 0.06})
            ax_d, ax_bar_dm, ax_bar_raw = axes
            dendrogram(Z, ax=ax_d, no_labels=True, color_threshold=0.0,
                       link_color_func=link_color, above_threshold_color="#888888")
            mean_rho = float(np.nanmean(rho))
            ax_d.set_title(
                rf"{pat} {BRAIN_BAND_TEX_DICT[band]}  --  "
                rf"$\overline{{\rho_\ell}}_{{p,b}} = {mean_rho:.2f}$  "
                rf"(links + leaves coloured by $\rho_\ell - \overline{{\rho_\ell}}_{{p,b}}$)",
                fontsize=10,
            )
            ax_bar_dm.imshow(rho_dm_ord.reshape(1, -1), aspect="auto",
                             cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            ax_bar_dm.set_yticks([])
            ax_bar_dm.set_xticks([])
            ax_bar_dm.set_ylabel(r"$\rho_\ell - \bar{\rho}$", rotation=0, ha="right", va="center", fontsize=8)
            vmax_raw = max(0.5, float(np.nanmax(np.abs(rho))))
            ax_bar_raw.imshow(rho_ord.reshape(1, -1), aspect="auto",
                              cmap="RdBu_r", vmin=-vmax_raw, vmax=vmax_raw)
            ax_bar_raw.set_yticks([])
            ax_bar_raw.set_xticks([])
            ax_bar_raw.set_ylabel(r"$\rho_\ell$ (raw)", rotation=0, ha="right", va="center", fontsize=8)
            # Two colorbars, one per bar
            sm_dm = ScalarMappable(norm=norm, cmap="RdBu_r")
            cb_dm = fig.colorbar(sm_dm, ax=ax_bar_dm, fraction=0.04, pad=0.01,
                                 location="right")
            cb_dm.set_label(r"demeaned", fontsize=8)
            sm_raw = ScalarMappable(norm=Normalize(vmin=-vmax_raw, vmax=vmax_raw), cmap="RdBu_r")
            cb_raw = fig.colorbar(sm_raw, ax=ax_bar_raw, fraction=0.04, pad=0.01,
                                  location="right")
            cb_raw.set_label(r"raw", fontsize=8)
            fig.tight_layout()
            fig.savefig(FIG_DIR / f"{pat}_{band}.pdf")
            plt.close(fig)

    # Cohort violin: rho_demeaned distributions per band
    fig, ax = plt.subplots(figsize=(9, 4.2))
    data = [df[df["band"] == b]["rho_demeaned"].values for b in BAND_ORDER]
    parts = ax.violinplot(data, positions=range(1, len(BAND_ORDER) + 1),
                          showmeans=True, showmedians=False)
    for pc in parts["bodies"]:
        pc.set_facecolor("#9ec6e5")
        pc.set_edgecolor("#1a4f73")
        pc.set_alpha(0.8)
    ax.axhline(0, color="0.4", lw=0.7, ls="--")
    ax.axhline(0.3, color="#d62728", lw=0.6, ls="--", label=r"$\rho_\ell - \bar{\rho} \geq 0.3$")
    ax.set_xticks(range(1, len(BAND_ORDER) + 1))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER])
    ax.set_ylabel(r"$\rho_\ell - \bar{\rho}_{(p, b)}$ -- positive = local peak above patient baseline")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cohort_distribution.pdf")
    plt.close(fig)
    print(f"[audit_39b] outputs at {OUT_DIR}; figures at {FIG_DIR}")


if __name__ == "__main__":
    main()
