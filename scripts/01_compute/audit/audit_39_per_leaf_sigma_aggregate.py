#!/usr/bin/env python3
"""Audit 39 — Per-leaf rho_ell localization (replaces TA-CBR).

Uses the existing controlled CTM per-pair Delta_task and Delta_rest data at
``data/reports/imcoh_continuous_trace/per_pair_split/{Pat_XX}_{band}.npz``.
For each leaf ell, define the per-leaf trace strength as the Spearman
rank correlation over the (N-1) pairs containing ell:

    rho_ell(p, b) = Spearman( Delta_task[ell, j], Delta_rest[ell, j] )
                    for j != ell

Positive rho_ell = leaf ell's per-pair distance shifts in task correlate
in rank with shifts in rest_post -- ell's neighbourhood reorganized in the
trace direction. Negative = anti-trace. Near zero = chance.

This is the per-leaf decomposition of the cohort-level rho_split that drives
the controlled CTM verdict in measure 05. Unlike the per-leaf sigma-mean
(which is just a rescaling of frac_pos_sigma and inflates when sign-agreement
exceeds 0.5 even at zero rank correlation), per-leaf rho agrees with the
cohort rank-correlation statistic.

Trace-leaf flag: rho_ell >= 0.30 (descriptive). Stricter rho_ell >= 0.50.

Outputs
-------
``data/audit/per_leaf_sigma_aggregate/leaf_rho.csv``
``data/audit/per_leaf_sigma_aggregate/cohort_summary.csv``
``data/outputs/figures/section_5_lrg_trace/per_leaf_sigma/{Pat_XX}_{band}.pdf``
``data/outputs/figures/section_5_lrg_trace/per_leaf_sigma/cohort_distribution.pdf``
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, leaves_list
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import _load_Z, BAND_ORDER, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE  # noqa: E402

CTM_PAIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"
OUT_DIR = ROOT / "data" / "audit" / "per_leaf_sigma_aggregate"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "per_leaf_sigma"


def per_leaf_rho(npz_path: Path) -> tuple[np.ndarray, int]:
    """Per-leaf Spearman rho between Delta_task and Delta_rest over pairs containing ell."""
    d = np.load(npz_path)
    dD_task = d["dD_task"]
    dD_rest = d["dD_rest"]
    iu_i = d["iu_i"].astype(int)
    iu_j = d["iu_j"].astype(int)
    n = int(max(iu_i.max(), iu_j.max())) + 1
    rho = np.zeros(n, dtype=float)
    for ell in range(n):
        mask = (iu_i == ell) | (iu_j == ell)
        if mask.sum() < 3:
            rho[ell] = np.nan
            continue
        a, b = dD_task[mask], dD_rest[mask]
        r, _ = spearmanr(a, b)
        rho[ell] = float(r) if r is not None and not np.isnan(r) else 0.0
    return rho, n


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    leaf_rows = []
    cohort_rows = []
    for pat in PATIENTS_4PHASE:
        for band in BAND_ORDER:
            npz = CTM_PAIR / f"{pat}_{band}.npz"
            if not npz.exists():
                continue
            rho, n = per_leaf_rho(npz)
            for ell in range(n):
                leaf_rows.append({
                    "patient": pat, "band": band, "leaf_id": ell,
                    "rho_leaf": float(rho[ell]) if not np.isnan(rho[ell]) else 0.0,
                    "is_trace_leaf_05": int(rho[ell] >= 0.5),
                    "is_trace_leaf_03": int(rho[ell] >= 0.3),
                })
            cohort_rows.append({
                "patient": pat, "band": band, "n_leaves": n,
                "rho_leaf_mean": float(np.nanmean(rho)),
                "rho_leaf_median": float(np.nanmedian(rho)),
                "rho_leaf_max": float(np.nanmax(rho)),
                "n_trace_05": int((rho >= 0.5).sum()),
                "n_trace_03": int((rho >= 0.3).sum()),
                "frac_trace_05": float((rho >= 0.5).mean()),
                "frac_trace_03": float((rho >= 0.3).mean()),
            })
            # Per-(patient, band) figure
            Z = _load_Z(pat, "rest_post", band)
            if Z is not None:
                fig, axes = plt.subplots(2, 1, figsize=(10, 4.5),
                                         gridspec_kw={"height_ratios": [4, 0.5], "hspace": 0.05})
                ax_d, ax_bar = axes
                vmax = max(0.3, float(np.nanmax(np.abs(rho))))
                dendrogram(Z, ax=ax_d, no_labels=True, color_threshold=0.0,
                           link_color_func=lambda nid: "#888888",
                           above_threshold_color="#888888")
                order = leaves_list(Z)
                rho_ord = rho[order]
                ax_bar.imshow(rho_ord.reshape(1, -1), aspect="auto",
                              cmap="RdBu_r", vmin=-vmax, vmax=vmax)
                ax_bar.set_yticks([])
                ax_bar.set_xticks([])
                ax_d.set_title(rf"{pat} {BRAIN_BAND_TEX_DICT[band]} -- per-leaf $\rho_\ell$ (T_RPost order)",
                               fontsize=11)
                ax_bar.set_xlabel(r"leaf order (rest_post dendrogram)", fontsize=9)
                cbar = fig.colorbar(ax_bar.images[0], ax=ax_bar, fraction=0.06, pad=0.02,
                                    orientation="horizontal", location="bottom")
                cbar.set_label(r"$\rho_\ell$ -- positive = trace direction", fontsize=9)
                fig.tight_layout()
                fig.savefig(FIG_DIR / f"{pat}_{band}.pdf")
                plt.close(fig)

    leaf_df = pd.DataFrame(leaf_rows)
    leaf_df.to_csv(OUT_DIR / "leaf_rho.csv", index=False)
    cohort_df = pd.DataFrame(cohort_rows)
    cohort_df.to_csv(OUT_DIR / "cohort_summary.csv", index=False)

    # Per-band cohort summary
    band_summary = (
        cohort_df.groupby("band")
        .agg(
            mean_rho_mean=("rho_leaf_mean", "mean"),
            median_rho_mean=("rho_leaf_mean", "median"),
            mean_frac_trace_05=("frac_trace_05", "mean"),
            median_frac_trace_05=("frac_trace_05", "median"),
            mean_frac_trace_03=("frac_trace_03", "mean"),
            n_pat_ge_5pct_trace_05=("frac_trace_05", lambda s: int((s >= 0.05).sum())),
            n_pat_ge_10pct_trace_05=("frac_trace_05", lambda s: int((s >= 0.10).sum())),
            n_pat_ge_20pct_trace_05=("frac_trace_05", lambda s: int((s >= 0.20).sum())),
        )
        .reset_index()
    )
    band_summary.to_csv(OUT_DIR / "cohort_band_summary.csv", index=False)
    print(band_summary.to_string(index=False))

    # Cohort violin
    fig, ax = plt.subplots(figsize=(9, 4.2))
    data = [leaf_df[leaf_df["band"] == b]["rho_leaf"].values for b in BAND_ORDER]
    parts = ax.violinplot(data, positions=range(1, len(BAND_ORDER) + 1),
                          showmeans=True, showmedians=False)
    for pc in parts["bodies"]:
        pc.set_facecolor("#9ec6e5")
        pc.set_edgecolor("#1a4f73")
        pc.set_alpha(0.8)
    ax.axhline(0, color="0.4", lw=0.7, ls="--")
    ax.axhline(0.5, color="#d62728", lw=0.6, ls="--", label=r"$\rho_\ell \geq 0.5$")
    ax.axhline(0.3, color="#7f7f7f", lw=0.6, ls=":", label=r"$\rho_\ell \geq 0.3$")
    ax.set_xticks(range(1, len(BAND_ORDER) + 1))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER])
    ax.set_ylabel(r"$\rho_\ell = \mathrm{Spearman}(\Delta_\mathrm{task}[\ell,\cdot], \Delta_\mathrm{rest}[\ell,\cdot])$")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cohort_distribution.pdf")
    plt.close(fig)
    print(f"[audit_39] outputs at {OUT_DIR}; figures at {FIG_DIR}")


if __name__ == "__main__":
    main()
