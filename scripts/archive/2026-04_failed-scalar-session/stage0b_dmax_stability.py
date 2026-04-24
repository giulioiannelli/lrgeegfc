#!/usr/bin/env python3
"""Stage 0b — tree-height (dmax) stability across phases.

For each (patient, band) load the LRG linkage matrices for the 4 phases
under `imcoh_abs` and compare their top-merge heights (`dmax = Z[-1, 2]`).

Purpose: decide whether fractional-height cuts (`h / dmax`) are meaningful
across phases. If dmax varies by ≪10% between phases within a (pat, band),
fractional cuts correspond to comparable absolute structural levels. If it
varies wildly, any h_rel-based analysis downstream must normalize carefully.

Outputs
-------
- ``data/reports/imcoh_vi/stage0b_dmax_table.csv``  — long CSV
- ``data/reports/imcoh_vi/figures/stage0b_dmax.pdf`` — two panels:
  (1) dmax per phase, one box per phase, one color per band;
  (2) ratio `dmax(phase) / dmax(rest_pre)` as a diagnostic for how much
      the tree stretches/compresses across phases.
- ``data/reports/imcoh_vi/stage0b_dmax.md`` — summary with per-band medians.
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE, PHASE_LABELS,
)
from lrg_eegfc.config.paths import REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

FC_METHOD = "imcoh_abs"
PATIENTS = [p for p in PATIENTS_4PHASE if p != "Pat_14"]  # Pat_14 has corrupt task_test
OUT_DIR = REPORTS_ROOT / "imcoh_vi"
FIG_DIR = OUT_DIR / "figures"


def collect() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS:
        for band in BRAIN_BANDS_NAMES:
            for phase in PHASE_LABELS:
                lrg = load_lrg_result(pat, phase, band, FC_METHOD)
                if lrg is None:
                    print(f"  MISSING: {pat} {phase} {band}")
                    continue
                Z = lrg.linkage_matrix
                rows.append({
                    "patient": pat,
                    "band": band,
                    "phase": phase,
                    "dmax": float(Z[-1, 2]),
                    "dmin": float(Z[Z[:, 2] > 0, 2].min()) if np.any(Z[:, 2] > 0) else 0.0,
                    "n_nodes": lrg.n_nodes,
                })
    return pd.DataFrame(rows)


def plot(df: pd.DataFrame, out: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # Panel 1 — dmax per (phase, band)
    band_colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(BRAIN_BANDS_NAMES)))
    width = 0.12
    x = np.arange(len(PHASE_LABELS))
    for i, band in enumerate(BRAIN_BANDS_NAMES):
        meds = []
        for phase in PHASE_LABELS:
            sub = df[(df["band"] == band) & (df["phase"] == phase)]["dmax"].values
            meds.append(np.median(sub) if len(sub) else np.nan)
        ax1.bar(x + (i - 2.5) * width, meds, width=width,
                label=BRAIN_BAND_TEX_DICT[band], color=band_colors[i])
    ax1.set_xticks(x)
    ax1.set_xticklabels(PHASE_LABELS, rotation=15, ha="right")
    ax1.set_ylabel(r"median $d_{\max}$ across 9 patients")
    ax1.set_title(r"Top-merge height by phase & band (|ImCoh|)")
    ax1.legend(fontsize=8, ncol=2)
    ax1.grid(alpha=0.25, axis="y")

    # Panel 2 — ratio dmax(phase)/dmax(rest_pre), per (patient, band, phase)
    pivot = df.pivot_table(index=["patient", "band"], columns="phase", values="dmax")
    if "rest_pre" in pivot.columns:
        for ph in ["task_learn", "task_test", "rest_post"]:
            if ph not in pivot.columns:
                continue
            ratios = (pivot[ph] / pivot["rest_pre"]).dropna()
            ax2.scatter(np.full(len(ratios), PHASE_LABELS.index(ph)),
                        ratios.values, s=14, alpha=0.5,
                        label=f"{ph} / rest_pre")
    ax2.axhline(1.0, color="k", lw=0.8, ls="--")
    ax2.set_xticks(range(len(PHASE_LABELS)))
    ax2.set_xticklabels(PHASE_LABELS, rotation=15, ha="right")
    ax2.set_ylabel(r"$d_{\max}(\mathrm{phase}) \;/\; d_{\max}(\mathrm{rest\_pre})$")
    ax2.set_title("Phase-to-rest_pre ratio per (patient, band)")
    ax2.grid(alpha=0.25, axis="y")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def summarize(df: pd.DataFrame) -> str:
    out = ["# Stage 0b — dmax stability diagnostic", "", ""]
    out.append(f"FC: `{FC_METHOD}`; patients: n={len(PATIENTS)} ({', '.join(PATIENTS)}); "
               f"4 phases, 6 bands.")
    out.append(f"Rows loaded: {len(df)}; expected {len(PATIENTS)*4*6}={len(PATIENTS)*24}.")
    out.append("")
    out.append("## Median dmax per (band, phase)")
    out.append("")
    pivot = df.pivot_table(index="band", columns="phase", values="dmax",
                           aggfunc="median").reindex(BRAIN_BANDS_NAMES)[list(PHASE_LABELS)]
    out.append(pivot.round(4).to_markdown())
    out.append("")
    out.append("## Ratio dmax(phase) / dmax(rest_pre), median across patients")
    out.append("")
    rel = (
        df.pivot_table(index=["patient", "band"], columns="phase", values="dmax")
        .assign(tl_ratio=lambda d: d["task_learn"] / d["rest_pre"],
                tt_ratio=lambda d: d["task_test"] / d["rest_pre"],
                rp_ratio=lambda d: d["rest_post"] / d["rest_pre"])
        .reset_index()
    )
    ratio_tbl = rel.groupby("band")[["tl_ratio", "tt_ratio", "rp_ratio"]].median() \
                   .reindex(BRAIN_BANDS_NAMES)
    out.append(ratio_tbl.round(3).to_markdown())
    out.append("")
    out.append("## Interpretation")
    out.append("")
    out.append("- Ratios near 1.0 → fractional-height cuts (`h/dmax`) are comparable "
               "across phases; h_rel-based analyses are safe.")
    out.append("- Systematic ratio ≠ 1 in some bands → fractional cuts within that band "
               "reflect structurally different levels across phases.")
    out.append("- Downstream implication: if ratios are stable, Stage 1's scale-aware "
               "metrics (TopK, Weighted ARI/VI at fixed k) remain comparable; if not, "
               "Stage 4's band × scale map must normalize scales more carefully.")
    return "\n".join(out)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/3] Loading dmax from cached Z matrices…")
    df = collect()
    csv_path = OUT_DIR / "stage0b_dmax_table.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Wrote {csv_path} ({len(df)} rows)")

    print("[2/3] Plot…")
    fig_path = FIG_DIR / "stage0b_dmax.pdf"
    plot(df, fig_path)
    print(f"  Wrote {fig_path}")

    print("[3/3] Summary…")
    md = summarize(df)
    md_path = OUT_DIR / "stage0b_dmax.md"
    md_path.write_text(md)
    print(f"  Wrote {md_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
