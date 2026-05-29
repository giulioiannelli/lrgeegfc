#!/usr/bin/env python3
"""Audit 33 — CTM (continuous-trace matrix) sigma/rho re-mining.

Re-mines the existing controlled CTM outputs at
``data/reports/imcoh_continuous_trace/`` into the triangle form used by
measures 02-04. CTM produces, per (patient, band):

    rho_split = Spearman(Delta_task, Delta_rest)
                where Delta_task = D^task_test - D^rest_pre,
                      Delta_rest = D^rest_post - D^rest_pre,
                computed with INDEPENDENT half-baselines D^rest_pre_A and
                D^rest_pre_B (Run A control: split-baseline).

CTM trace direction is rho_split > 0 (task-induced and rest-induced
distance shifts agree in sign). **Project-wide sign convention: T > 0 =
trace, T < 0 = anti-trace, T = 0 = no trace.** We report

    T_CTM(p, b) = rho_split(p, b)

so positive = trace. Cohort summary mirrors the existing
``controls_band_stats.md`` numbers in the (positive = trace) convention.

Outputs
-------
``data/audit/ctm_triangle/Td_per_patient_per_band.csv``
``data/audit/ctm_triangle/cohort_summary.csv``
``data/outputs/figures/section_5_lrg_trace/ctm_triangle/boxplots.pdf``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT

CTM_SPLIT = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_cell_summary_split.csv"
CTM_CTRL = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "controls_summary.csv"
OUT_DIR = ROOT / "data" / "audit" / "ctm_triangle"
FIG_DIR = ROOT / "data" / "outputs" / "figures" / "section_5_lrg_trace" / "ctm_triangle"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    split = pd.read_csv(CTM_SPLIT)
    ctrl = pd.read_csv(CTM_CTRL)

    df = split[["patient", "band", "rho", "frac_pos_sigma"]].copy()
    df = df.rename(columns={"rho": "rho_split", "frac_pos_sigma": "frac_pos_sigma_shared"})
    df = df.merge(
        ctrl[["patient", "band", "rho_shared", "rho_null_drift",
              "rho_split_same_probe", "rho_split_cross_probe"]],
        on=["patient", "band"],
        how="left",
    )
    df["T_CTM"] = df["rho_split"]                    # positive = trace
    df["passes_drift"] = (df["rho_split"] > df["rho_null_drift"]).astype(int)
    df.to_csv(OUT_DIR / "Td_per_patient_per_band.csv", index=False)

    # Per-band cohort summary: include controlled cohort scalars
    summary_rows = []
    for b in BRAIN_BANDS_NAMES:
        sub = df[df["band"] == b]
        v = sub["rho_split"].dropna().values
        v_drift = sub["rho_null_drift"].dropna().values
        v_xp = sub["rho_split_cross_probe"].dropna().values
        n_trace_split = int((v > 0).sum())                    # rho > 0 = trace
        n_trace_xp = int((v_xp > 0).sum())
        n_above_drift = int((v > v_drift).sum())
        try:
            _, p_one = wilcoxon(v, alternative="greater")
        except Exception:
            p_one = float("nan")
        try:
            _, p_drift = wilcoxon(v - v_drift, alternative="greater")
        except Exception:
            p_drift = float("nan")
        summary_rows.append({
            "band": b,
            "n_trace_split": f"{n_trace_split}/{len(v)}",
            "n_trace_split_int": n_trace_split,
            "rho_split_median": float(np.median(v)),
            "n_above_drift": f"{n_above_drift}/{len(v)}",
            "n_above_drift_int": n_above_drift,
            "n_trace_xprobe": f"{n_trace_xp}/{len(v_xp)}",
            "n_trace_xprobe_int": n_trace_xp,
            "rho_xprobe_median": float(np.median(v_xp)),
            "wilcoxon_split_gt_0_p": float(p_one),
            "wilcoxon_split_gt_drift_p": float(p_drift),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "cohort_summary.csv", index=False)
    print(summary.to_string(index=False))

    fig, ax = plt.subplots(figsize=(8, 4))
    data_split = [df[df["band"] == b]["rho_split"].dropna().values for b in BRAIN_BANDS_NAMES]
    data_drift = [df[df["band"] == b]["rho_null_drift"].dropna().values for b in BRAIN_BANDS_NAMES]
    pos = np.arange(1, len(BRAIN_BANDS_NAMES) + 1)
    bp_a = ax.boxplot(data_split, positions=pos - 0.18, widths=0.32,
                      patch_artist=True, boxprops=dict(facecolor="#4c92c3"))
    bp_b = ax.boxplot(data_drift, positions=pos + 0.18, widths=0.32,
                      patch_artist=True, boxprops=dict(facecolor="#cccccc"))
    for i, b in enumerate(BRAIN_BANDS_NAMES):
        v_a = df[df["band"] == b]["rho_split"].dropna().values
        v_b = df[df["band"] == b]["rho_null_drift"].dropna().values
        ax.scatter(np.full(len(v_a), pos[i] - 0.18), v_a, s=10, color="#1a4f73", alpha=0.7, zorder=3)
        ax.scatter(np.full(len(v_b), pos[i] + 0.18), v_b, s=10, color="#666666", alpha=0.7, zorder=3)
    ax.set_xticks(pos)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES])
    ax.axhline(0, color="0.4", lw=0.7, ls="--")
    ax.set_ylabel(r"$\rho$ (positive = trace direction)")
    ax.legend([bp_a["boxes"][0], bp_b["boxes"][0]],
              [r"$\rho_{\mathrm{split}}$ (Run A)", r"$\rho_{\mathrm{null\_drift}}$ (Run C)"],
              loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "boxplots.pdf")
    plt.close(fig)
    print(f"[audit_33] outputs at {OUT_DIR}; figure at {FIG_DIR}")


if __name__ == "__main__":
    main()
