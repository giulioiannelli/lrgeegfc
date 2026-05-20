#!/usr/bin/env python3
"""Audit 58 — Band-agnostic LRG trace, first-pass via band-stacked CTM proxy.

Question
--------
§5 of `notes_imcoh.tex` shows the cognitive-task trace at six canonical
EEG bands (δ θ α β γ_l γ_h). Three bands carry the per-pair correlation
trace at q=0.027 (α / β / γ_l) under the Run-A controlled CTM
formulation. Is band specificity essential, or does a single
band-agnostic LRG read recover the same trace direction?

Strategy
--------
A canonical broadband |ImCoh| pipeline would re-compute |ImCoh| on the
union spectrum, feed to LRG, and re-derive split-half ultrametric
shifts. That requires:
  (a) a broadband freq-resolved cache (does NOT exist; per-band only at
      `data/cache/imcoh/Pat_NN/{band}_{phase}_imcoh_freqresolved_nperseg-4096.npy`),
  (b) a half-baseline broadband LRG cache (does NOT exist),
  (c) ~30 minutes of LRG compute.

This script implements the FAST PROXY from the scope report:

  Stack the existing per-band Δ_task(i,j) and Δ_rest(i,j) vectors —
  produced by `continuous_trace_matrix.py --mode split-baseline` and
  cached at `data/reports/imcoh_continuous_trace/per_pair_split/` —
  across all 6 bands per patient, then compute a single broadband
  ρ_split per patient.

Two stacking conventions are reported:

  ρ_split^pool   Spearman on the concatenated 6 × m_p vector (pooled
                 rank-agreement across pair-and-band).
  ρ_split^band   Mean of per-band ρ_split (cohort-level "band-average").

A drift floor is reconstructed analogously by stacking the band-level
ρ_null_drift values from `controls_summary.csv`. The cohort-paired
Wilcoxon `ρ_split^pool > ρ_drift^pool` is the headline number.

Caveats (full list in scope report)
-----------------------------------
* The pool proxy is a *band-pool*, not a *spectral broadband*: it does
  not literally average spectral bins. It is a band-rank pooling.
* The canonical pipeline (re-compute |ImCoh| on the union spectrum and
  feed to LRG) is in §11 of the scope; runtime ≤ 30 minutes; not run
  here.
* 1/f weighting under canonical broadband would amplify δ; band-pool
  proxy here gives equal weight to each band.
* Pat_03 at 1024 Hz has half the high-γ coverage; flagged.

Outputs
-------
``data/audit/band_agnostic_lrg/cohort_summary.csv``
``data/audit/band_agnostic_lrg/per_patient_summary.csv``
``data/reports/notes_verification_2026-05-08/figures/band_agnostic_lrg_teaser.pdf``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z

CTM_PAIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"
CONTROLS = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "controls_summary.csv"
TD_CSV = ROOT / "data" / "audit" / "ctm_triangle" / "Td_per_patient_per_band.csv"

OUT = ROOT / "data" / "audit" / "band_agnostic_lrg"
FIG = ROOT / "data" / "reports" / "notes_verification_2026-05-08" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
            "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = list(BRAIN_BANDS_NAMES)
TRACE_BANDS_TEX = r"$\alpha$, $\beta$, $\gamma_{\mathrm{l}}$"

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
})


def stack_per_patient(pat: str) -> dict | None:
    """Concatenate per-band Δ_task / Δ_rest into a single broadband vector."""
    parts_task = []
    parts_rest = []
    band_lengths = {}
    for band in BANDS:
        npz_path = CTM_PAIR / f"{pat}_{band}.npz"
        if not npz_path.exists():
            return None
        d = np.load(npz_path)
        parts_task.append(np.asarray(d["dD_task"], dtype=float))
        parts_rest.append(np.asarray(d["dD_rest"], dtype=float))
        band_lengths[band] = parts_task[-1].size
    dD_task_pool = np.concatenate(parts_task)
    dD_rest_pool = np.concatenate(parts_rest)
    return dict(
        dD_task=dD_task_pool, dD_rest=dD_rest_pool,
        band_lengths=band_lengths,
    )


def main() -> None:
    print("[audit_58] band-stacked CTM proxy for broadband ρ_split")
    print(f"[audit_58] cohort = {len(PATIENTS)} patients × {len(BANDS)} bands")

    # Load per-band ρ table (Td_per_patient_per_band.csv has rho_split column)
    td = pd.read_csv(TD_CSV)
    ctrls = pd.read_csv(CONTROLS)

    rows = []
    for pat in PATIENTS:
        st = stack_per_patient(pat)
        if st is None:
            print(f"[audit_58] SKIP {pat}: missing per-band per-pair NPZ")
            continue
        rho_pool, p_pool = spearmanr(st["dD_task"], st["dD_rest"])
        # Per-patient mean of per-band ρ_split (cohort proxy 2)
        per_band_rho = td[td["patient"] == pat].set_index("band")["rho_split"]
        rho_band_mean = float(per_band_rho.reindex(BANDS).mean())
        # Band-pool drift: same construction but using (rest_pre_B - rest_pre_A)
        # vs (rest_post_B - rest_post_A) is not stored as per-pair NPZ; use the
        # CSV's already-aggregated rho_null_drift as the band-mean drift.
        per_band_drift = ctrls[ctrls["patient"] == pat].set_index("band")["rho_null_drift"]
        rho_drift_band_mean = float(per_band_drift.reindex(BANDS).mean())

        # Best-band ρ for this patient among the manuscript's α / β / γ_l set
        trace_bands = ["alpha", "beta", "low_gamma"]
        rho_best_trace = float(per_band_rho.reindex(trace_bands).max())
        best_trace_band = per_band_rho.reindex(trace_bands).idxmax()

        m_total = int(st["dD_task"].size)
        rows.append({
            "patient": pat,
            "rho_split_pool": float(rho_pool),
            "rho_split_band_mean": rho_band_mean,
            "rho_drift_band_mean": rho_drift_band_mean,
            "rho_best_trace_band": rho_best_trace,
            "best_trace_band": str(best_trace_band),
            "n_total_pairs": m_total,
        })

    per_pat = pd.DataFrame(rows)
    per_pat.to_csv(OUT / "per_patient_summary.csv", index=False)

    # Cohort summary
    n = len(per_pat)
    n_pool_pos = int((per_pat["rho_split_pool"] > 0).sum())
    n_band_mean_pos = int((per_pat["rho_split_band_mean"] > 0).sum())
    med_pool = float(per_pat["rho_split_pool"].median())
    med_band = float(per_pat["rho_split_band_mean"].median())
    med_drift = float(per_pat["rho_drift_band_mean"].median())

    z_pool, p_pool_gt0 = wilcoxon_z(per_pat["rho_split_pool"].values)
    diff_pool_drift = (per_pat["rho_split_pool"].values
                       - per_pat["rho_drift_band_mean"].values)
    z_pool_drift, p_pool_drift = wilcoxon_z(diff_pool_drift)
    diff_band_drift = (per_pat["rho_split_band_mean"].values
                       - per_pat["rho_drift_band_mean"].values)
    z_band_drift, p_band_drift = wilcoxon_z(diff_band_drift)

    n_pool_above_drift = int((diff_pool_drift > 0).sum())
    n_band_above_drift = int((diff_band_drift > 0).sum())

    cohort = pd.DataFrame([{
        "n_patients": n,
        "n_pool_pos": n_pool_pos,
        "n_band_mean_pos": n_band_mean_pos,
        "n_pool_above_drift": n_pool_above_drift,
        "n_band_above_drift": n_band_above_drift,
        "median_rho_pool": med_pool,
        "median_rho_band_mean": med_band,
        "median_rho_drift": med_drift,
        "wilcoxon_pool_gt_0_z": z_pool,
        "wilcoxon_pool_gt_0_p": p_pool_gt0,
        "wilcoxon_pool_gt_drift_z": z_pool_drift,
        "wilcoxon_pool_gt_drift_p": p_pool_drift,
        "wilcoxon_band_gt_drift_z": z_band_drift,
        "wilcoxon_band_gt_drift_p": p_band_drift,
    }])
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)

    print()
    print(f"[audit_58] cohort n = {n}")
    print(f"  rho_split_pool      median = {med_pool:+.3f}   "
          f"n_pos = {n_pool_pos}/{n}   "
          f"n_above_drift = {n_pool_above_drift}/{n}")
    print(f"  rho_split_band_mean median = {med_band:+.3f}   "
          f"n_pos = {n_band_mean_pos}/{n}   "
          f"n_above_drift = {n_band_above_drift}/{n}")
    print(f"  rho_drift_band_mean median = {med_drift:+.3f}")
    print(f"  Wilcoxon pool > 0     : z = {z_pool:+.3f}, p = {p_pool_gt0:.4f}")
    print(f"  Wilcoxon pool > drift : z = {z_pool_drift:+.3f}, p = {p_pool_drift:.4f}")
    print(f"  Wilcoxon band > drift : z = {z_band_drift:+.3f}, p = {p_band_drift:.4f}")

    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # Panel A — ρ_split^pool vs ρ_drift connected per patient
    axA = axes[0]
    pos = np.array([0, 1.0])
    for _, r in per_pat.iterrows():
        axA.plot(pos, [r["rho_drift_band_mean"], r["rho_split_pool"]],
                 color="#888888", lw=0.7, alpha=0.7, zorder=1)
    axA.scatter(np.full(n, pos[0]), per_pat["rho_drift_band_mean"].values,
                s=42, color="#bbbbbb", edgecolors="#444444", lw=0.6,
                zorder=2, label=r"$\bar{\rho}_{\mathrm{drift}}$ (band mean)")
    axA.scatter(np.full(n, pos[1]), per_pat["rho_split_pool"].values,
                s=46, color="#1f3d6e", edgecolors="white", lw=0.6,
                zorder=3, label=r"$\rho_{\mathrm{split}}^{\mathrm{pool}}$ (broadband)")
    axA.axhline(0, color="0.4", lw=0.7, ls="--", zorder=0)
    axA.set_xticks(pos)
    axA.set_xticklabels([r"drift", r"broadband"])
    axA.set_xlim(-0.4, 1.4)
    axA.set_ylabel(r"Spearman $\rho$ (positive = trace direction)")
    axA.set_title(f"Per-patient pairing — Wilcoxon "
                  rf"$p = {p_pool_drift:.4f}$,  $n_{{>}} = {n_pool_above_drift}/{n}$",
                  fontsize=10)
    axA.legend(loc="lower right", frameon=False, fontsize=9)
    axA.spines[["top", "right"]].set_visible(False)

    # Panel B — broadband vs best-trace-band
    axB = axes[1]
    bands_color = {"alpha": "#3aa57d", "beta": "#d6603a", "low_gamma": "#a070c0"}
    for _, r in per_pat.iterrows():
        c = bands_color.get(r["best_trace_band"], "#444444")
        axB.scatter(r["rho_best_trace_band"], r["rho_split_pool"],
                    s=58, color=c, edgecolors="white", lw=0.7, zorder=3)
        axB.annotate(r["patient"].replace("Pat_", ""),
                     (r["rho_best_trace_band"], r["rho_split_pool"]),
                     xytext=(4, 2), textcoords="offset points",
                     fontsize=7.5, color="#222222")
    lim = max(0.05,
              float(np.abs(per_pat[["rho_best_trace_band",
                                    "rho_split_pool"]]).max().max()) * 1.15)
    axB.plot([-lim, lim], [-lim, lim], color="#888888", lw=0.6, ls="--", zorder=0)
    axB.axhline(0, color="0.7", lw=0.5, zorder=0)
    axB.axvline(0, color="0.7", lw=0.5, zorder=0)
    axB.set_xlim(-lim, lim)
    axB.set_ylim(-lim, lim)
    axB.set_aspect("equal")
    axB.set_xlabel(r"per-patient best $\rho_{\mathrm{split}}$ "
                   rf"among ({TRACE_BANDS_TEX})")
    axB.set_ylabel(r"$\rho_{\mathrm{split}}^{\mathrm{pool}}$ (broadband)")
    axB.set_title("Best trace-band vs broadband", fontsize=10)
    axB.spines[["top", "right"]].set_visible(False)

    # Color legend for B
    handles = []
    import matplotlib.lines as mlines
    for b in ["alpha", "beta", "low_gamma"]:
        handles.append(mlines.Line2D(
            [], [], marker="o", linestyle="None",
            markersize=8, markerfacecolor=bands_color[b],
            markeredgecolor="white",
            label=BRAIN_BAND_TEX_DICT[b],
        ))
    axB.legend(handles=handles, loc="lower right", frameon=False,
               fontsize=9, title="best trace band")

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig_path = FIG / "band_agnostic_lrg_teaser.pdf"
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[audit_58] wrote {fig_path}")
    print(f"[audit_58] outputs at {OUT}")


if __name__ == "__main__":
    main()
