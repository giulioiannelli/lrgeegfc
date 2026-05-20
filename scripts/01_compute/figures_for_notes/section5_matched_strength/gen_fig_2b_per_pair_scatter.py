#!/usr/bin/env python3
"""§5.7 Figure 2b — Per-pair LRG-CTM scatter with matched-strength annotation.

Same layout as the manuscript per-pair scatter (`manuscript_ctm_per_pair_scatter.pdf`,
audit_55) at three cohort-median-anchored (patient, band) cells, plus
an annotation in each panel of the matched-strength cohort ratio and
paired-Wilcoxon p value for the band so the per-pair geometry stays
linked to the §5.7 controlled cohort claim.

Inputs:
    data/cache/lrg_ctm_pairs/Pat_NN_<band>.npz   (Δ_task, Δ_rest, iu_i, iu_j)
    data/audit/ctm_per_pair_scatter/tables/selection.csv
    data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv

Output:
    data/reports/section_5_matched_strength_refinement/figures/
        fig_2b_ctm_per_pair_scatter_matched_strength.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.lines as mlines
import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.probe import extract_probe_labels
from lrg_eegfc.visuals.styles import use_lrg_style

from _shared_ms import FIG_DIR, LRG_CTM_DIR


SEL_CSV = ROOT / "data" / "audit" / "ctm_per_pair_scatter" / "tables" / "selection.csv"
CTM_PAIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"


def load_pair_data(pat: str, band: str) -> dict:
    npz_path = CTM_PAIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    dD_task = np.asarray(d["dD_task"])
    dD_rest = np.asarray(d["dD_rest"])
    iu_i = d["iu_i"].astype(int)
    iu_j = d["iu_j"].astype(int)
    pat_dir = ROOT / "data" / "raw" / "stereoeeg_patients" / pat
    ch = pd.read_csv(pat_dir / "channel_labels.csv")
    labels = ch["label"].tolist()
    probes = np.array(extract_probe_labels(labels))
    n_contacts = int(max(iu_i.max(), iu_j.max())) + 1
    if probes.size != n_contacts:
        probes = probes[:n_contacts]
    same_probe = probes[iu_i] == probes[iu_j]
    return dict(dD_task=dD_task, dD_rest=dD_rest, same_probe=same_probe)


def main() -> Path:
    use_lrg_style()
    sel = pd.read_csv(SEL_CSV)
    cohort_ms = pd.read_csv(LRG_CTM_DIR / "cohort_summary.csv")

    panels = []
    for _, r in sel.iterrows():
        pair = load_pair_data(r.patient, r.band)
        rho_full, _ = spearmanr(pair["dD_task"], pair["dD_rest"])
        cp = ~pair["same_probe"]
        rho_x, _ = spearmanr(pair["dD_task"][cp], pair["dD_rest"][cp])
        ms_row = cohort_ms[cohort_ms.band == r.band].iloc[0]
        ms_ratio = (abs(ms_row.obs_median_rho)
                    / max(abs(ms_row.surr_median_rho_median), 1e-12))
        panels.append(dict(
            pat=r.patient, band=r.band,
            rho_full=float(rho_full), rho_x=float(rho_x),
            **pair,
            N=int(len(pair["dD_task"])),
            N_same=int(pair["same_probe"].sum()),
            N_cross=int((~pair["same_probe"]).sum()),
            ms_ratio=float(ms_ratio),
            ms_p=float(ms_row.paired_wilcoxon_p),
            ms_n_above=str(ms_row.n_above_surrogate),
        ))

    all_vals = np.concatenate(
        [p["dD_task"] for p in panels] + [p["dD_rest"] for p in panels]
    )
    q98 = float(np.quantile(np.abs(all_vals), 0.98))
    lim = q98 * 1.05

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.8))
    for ax, p in zip(axes, panels):
        sp = p["same_probe"]
        cp = ~sp
        ax.scatter(p["dD_task"][sp], p["dD_rest"][sp],
                   s=8, color="#bbbbbb", alpha=0.45,
                   edgecolors="none", zorder=1)
        ax.scatter(p["dD_task"][cp], p["dD_rest"][cp],
                   s=10, color="#1f3d6e", alpha=0.55,
                   edgecolors="none", zorder=2)
        ax.axhline(0, color="#cccccc", lw=0.6, zorder=0)
        ax.axvline(0, color="#cccccc", lw=0.6, zorder=0)
        ax.plot([-lim, lim], [-lim, lim], color="#888888", lw=0.6,
                ls="--", zorder=0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$\Delta_{\mathrm{task}}(i,j) = "
                      r"D_{\mathrm{task}} - D_{\mathrm{rsPre,A}}$")
        ax.set_ylabel(r"$\Delta_{\mathrm{rest}}(i,j) = "
                      r"D_{\mathrm{rsPost}} - D_{\mathrm{rsPre,B}}$")
        ax.spines[["top", "right"]].set_visible(False)

        band_tex = BRAIN_BAND_TEX_DICT[p["band"]]
        # Significance annotation colour: red if cohort-paired ≤ 0.05
        ms_col = "#c0392b" if p["ms_p"] <= 0.05 else "#222"
        ax.text(
            0.03, 0.97,
            f"{p['pat']}, {band_tex}\n"
            rf"$\rho_{{\rm split}} = {p['rho_full']:+.3f}$"
            "\n"
            rf"$\rho_{{\rm xprobe}} = {p['rho_x']:+.3f}$"
            "\n"
            f"$N$ = {p['N']} "
            f"({p['N_same']} same / {p['N_cross']} cross)",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=9.5,
            bbox=dict(facecolor="white", edgecolor="#bbbbbb",
                      boxstyle="round,pad=0.35"),
        )
        # Bottom-right: cohort matched-strength panel
        ax.text(
            0.97, 0.03,
            rf"cohort {band_tex}: ratio $={p['ms_ratio']:.1f}\times$"
            "\n"
            rf"$p_{{\rm MS}} = {p['ms_p']:.3f}$, "
            rf"$n_{{>}} = {p['ms_n_above']}$",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.5, color=ms_col,
            fontweight="bold" if p["ms_p"] <= 0.05 else "normal",
            bbox=dict(facecolor="white", edgecolor="#bbbbbb",
                      boxstyle="round,pad=0.32"),
        )

    handles = [
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=8,
                      markerfacecolor="#bbbbbb", markeredgecolor="none",
                      alpha=0.7, label="same-probe pairs"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=8,
                      markerfacecolor="#1f3d6e", markeredgecolor="none",
                      alpha=0.8, label="cross-probe pairs"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.04))
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    out = FIG_DIR / "fig_2b_ctm_per_pair_scatter_matched_strength.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
