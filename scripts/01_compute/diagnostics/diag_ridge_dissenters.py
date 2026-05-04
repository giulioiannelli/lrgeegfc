#!/usr/bin/env python3
"""Per-patient dissenter audit on the n=10 H2a Δ_VI ridges.

Three checks driven by the question: are the 1–2 ≥8/10-cohort dissenters
the same patient (Pat_03 1024 Hz outlier?) and if so, is the pattern
explainable by implant positioning rather than physiology?

  Check 1 — Per-patient ridge sign grid.
    For each headline ridge (band, k_range from the n=10 ≥8/10 unanimity
    audit), compute per-patient mean Δ_VI averaged over the ridge cells.
    Heatmap rows = (band, k_range) ridges, cols = patients, colour = mean
    Δ_VI. Negative cells = dissenters.

  Check 2 — Pat_03 in/out sensitivity.
    Recompute the n=10 ≥8/10 unanimity ridges with Pat_03 dropped (n=9).
    Report (a) which Δ_VI ridges survive at ≥7/9, (b) which strengthen,
    (c) which collapse. If the headline δ k=20–32 ridge survives without
    Pat_03 the cohort claim is robust to the 1024 Hz instrument; if it
    depends on Pat_03 the claim is fragile.

  Check 3 — Implant probe analysis for persistent dissenters.
    Identify any patient that dissents (mean Δ_VI < 0) on ≥ 3 of the 6
    primary headline ridges. For each, dump probe-prefix counts from
    channel_labels.csv vs the cohort distribution. If a single patient
    has unusually skewed probe coverage (one probe ≥ 30% of channels,
    or unusually low coverage of a region the cohort dominates) flag it.

Outputs:
  data/audit/ridge_diagnostics/per_patient_ridge_sign.pdf
  data/audit/ridge_diagnostics/per_patient_ridge_sign.csv
  data/audit/ridge_diagnostics/pat03_sensitivity.csv
  data/audit/ridge_diagnostics/dissenter_summary.md
"""
from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT, SEEG_DATAPATH


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

# Ridges from the n=10 ≥8/10 audit on Δ_VI (L >= 2 only).
RIDGES = [
    ("delta",       (10, 11), "δ k=10–11"),
    ("delta",       (20, 32), "δ k=20–32 ★"),
    ("alpha",       (20, 26), "α k=20–26"),
    ("beta",        ( 6,  7), "β k=6–7"),
    ("high_gamma",  (21, 23), "γ_h k=21–23"),
    ("high_gamma",  (30, 35), "γ_h k=30–35"),
]


def _diverging_cmap():
    return LinearSegmentedColormap.from_list(
        "ridge_div",
        [(0.05, 0.20, 0.55), (1.0, 1.0, 1.0), (0.70, 0.05, 0.05)],
        N=256,
    )


def _per_patient_ridge_means(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for band, (k_lo, k_hi), _ in RIDGES:
        sub = df[(df["band"] == band) & df["k"].between(k_lo, k_hi)]
        for pat in COHORT_N10:
            x = sub[sub["patient"] == pat]["d_VI"].to_numpy()
            x = x[np.isfinite(x)]
            rows.append({
                "band": band, "k_lo": k_lo, "k_hi": k_hi,
                "patient": pat,
                "mean_dVI": float(x.mean()) if x.size else np.nan,
                "frac_pos": float((x > 0).mean()) if x.size else np.nan,
                "n_cells": int(x.size),
            })
    return pd.DataFrame(rows)


def _check1_grid(df_means: pd.DataFrame, out_dir: Path) -> None:
    grid = (df_means
            .pivot_table(index=["band", "k_lo", "k_hi"], columns="patient",
                         values="mean_dVI", aggfunc="first")
            .reindex(columns=COHORT_N10))
    # Re-order rows to match RIDGES list
    row_keys = [(b, lo, hi) for (b, (lo, hi), _) in RIDGES]
    grid = grid.reindex(row_keys)
    row_labels = [lbl for (_, _, lbl) in RIDGES]

    fig, ax = plt.subplots(figsize=(11.0, 4.8), dpi=160)
    vmax = float(np.nanmax(np.abs(grid.values)))
    im = ax.imshow(grid.values, aspect="auto", cmap=_diverging_cmap(),
                   vmin=-vmax, vmax=+vmax, interpolation="nearest")
    ax.set_xticks(range(len(COHORT_N10)))
    ax.set_xticklabels(COHORT_N10, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=10)
    # Annotate each cell with sign + 3-decimal value
    for ir in range(grid.shape[0]):
        for ic in range(grid.shape[1]):
            v = grid.values[ir, ic]
            if not np.isfinite(v):
                continue
            colour = "white" if abs(v) > 0.55 * vmax else "black"
            ax.text(ic, ir, f"{v:+.2f}", ha="center", va="center",
                    fontsize=8, color=colour)
            if v < 0:
                ax.add_patch(plt.Rectangle((ic - 0.5, ir - 0.5), 1, 1,
                                            facecolor="none",
                                            edgecolor="black",
                                            linewidth=1.6, zorder=4))
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(r"per-patient mean $\Delta_{VI}$ over ridge cells",
                   fontsize=9)
    ax.set_title(
        "Per-patient ridge sign — black-bordered cells = dissenters "
        "(mean $\\Delta_{VI} < 0$)",
        loc="left", fontsize=10,
    )
    fig.tight_layout()
    out_pdf = out_dir / "per_patient_ridge_sign.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_pdf}")
    out_csv = out_dir / "per_patient_ridge_sign.csv"
    df_means.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")


def _check2_pat03_sensitivity(df: pd.DataFrame, out_dir: Path) -> dict:
    """Recompute ≥8/10 (n=10) and ≥7/9 (n=9, Pat_03 dropped) ridges."""
    rows = []
    for band, (k_lo, k_hi), label in RIDGES:
        sub = df[(df["band"] == band) & df["k"].between(k_lo, k_hi)]
        per_k_full = sub.pivot(index="patient", columns="k", values="d_VI")
        per_k_no03 = per_k_full.drop("Pat_03", errors="ignore")

        n_full_ridge_pos = ((per_k_full > 0).sum(axis=0) >= 8).sum()
        n_no03_ridge_pos = ((per_k_no03 > 0).sum(axis=0) >= 7).sum()

        full_mean = per_k_full.mean(axis=0).mean()
        no03_mean = per_k_no03.mean(axis=0).mean()
        full_pos_frac = ((per_k_full > 0).mean(axis=0)).mean()
        no03_pos_frac = ((per_k_no03 > 0).mean(axis=0)).mean()

        rows.append({
            "ridge": label,
            "k_lo": k_lo, "k_hi": k_hi, "L": k_hi - k_lo + 1,
            "cells_geq_8of10_n10": int(n_full_ridge_pos),
            "cells_geq_7of9_n9_no_p03": int(n_no03_ridge_pos),
            "mean_dVI_n10": float(full_mean),
            "mean_dVI_n9_no_p03": float(no03_mean),
            "pos_frac_n10": float(full_pos_frac),
            "pos_frac_n9_no_p03": float(no03_pos_frac),
        })
    df_sens = pd.DataFrame(rows)
    df_sens.to_csv(out_dir / "pat03_sensitivity.csv", index=False)
    print(f"wrote {out_dir / 'pat03_sensitivity.csv'}")
    return df_sens.to_dict(orient="records")


def _channel_probe_counts(pat: str) -> dict[str, int]:
    p = Path(SEEG_DATAPATH) / pat / "channel_labels.csv"
    if not p.exists():
        return {}
    labels = pd.read_csv(p).iloc[:, 0].astype(str)
    counts: Counter[str] = Counter()
    for lbl in labels:
        # Probe prefix = leading non-digit chars before first digit/space
        m = re.match(r"^([A-Za-z']+)", lbl.strip())
        if not m:
            continue
        counts[m.group(1).upper()] += 1
    return dict(counts)


def _check3_implant_dissenters(df_means: pd.DataFrame, out_dir: Path) -> str:
    """Identify patients dissenting on ≥3 ridges and dump probe-prefix counts."""
    grid = df_means.pivot_table(
        index="patient", columns=["band", "k_lo", "k_hi"],
        values="mean_dVI", aggfunc="first"
    ).reindex(COHORT_N10)
    n_dissent = (grid < 0).sum(axis=1)
    persistent = n_dissent[n_dissent >= 3].index.tolist()

    md = ["# Dissenter probe analysis (Δ_VI ridges, n=10)\n"]
    md.append(f"Total ridges checked: {len(RIDGES)}.\n")
    md.append("\n## Dissent count per patient\n\n")
    md.append("| patient | ridges with mean Δ_VI < 0 |\n|---|---:|\n")
    for pat in COHORT_N10:
        md.append(f"| {pat} | {int(n_dissent[pat])} |\n")
    md.append("\n")

    if not persistent:
        md.append("## Persistent dissenters (≥ 3 negative ridges)\n\n")
        md.append("**None.** No patient dissents on more than 2 of the 6 "
                   "headline ridges. The negative cells in Check 1 are "
                   "scattered across patients × ridges — no single patient "
                   "drives the ≥ 8/10 boundary.\n")
    else:
        md.append(f"## Persistent dissenters: {', '.join(persistent)}\n\n")
        cohort_probes: Counter[str] = Counter()
        for pat in COHORT_N10:
            for probe, n in _channel_probe_counts(pat).items():
                cohort_probes[probe] += n
        cohort_total = sum(cohort_probes.values())
        for pat in persistent:
            md.append(f"### {pat}\n\n")
            md.append("| probe | count | %_pat | %_cohort | skew |\n")
            md.append("|---|---:|---:|---:|---:|\n")
            cnt = _channel_probe_counts(pat)
            n_pat = sum(cnt.values()) or 1
            for probe in sorted(cnt, key=lambda x: -cnt[x]):
                pct = cnt[probe] / n_pat
                cohort_pct = cohort_probes.get(probe, 0) / cohort_total if cohort_total else 0
                skew = pct - cohort_pct
                md.append(f"| {probe} | {cnt[probe]} | {pct:.2%} | "
                           f"{cohort_pct:.2%} | {skew:+.2%} |\n")
            md.append("\n")
        md.append("**Skew interpretation.** Positive skew = this probe is "
                   "over-represented in the patient relative to cohort mean. "
                   "Skew > 10pp on a single probe is large; > 20pp warrants "
                   "investigation.\n")

    md_path = out_dir / "dissenter_summary.md"
    md_path.write_text("".join(md), encoding="utf-8")
    print(f"wrote {md_path}")
    return "".join(md)


def main() -> None:
    df = pd.read_csv(REPORTS_ROOT / "imcoh_vi" / "h2_partition_multiscale_raw.csv")
    out_dir = ROOT / "data" / "audit" / "ridge_diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)

    df_means = _per_patient_ridge_means(df)
    print("\n=== Check 1 — per-patient ridge sign grid ===")
    _check1_grid(df_means, out_dir)

    print("\n=== Check 2 — Pat_03 in/out sensitivity ===")
    sens = _check2_pat03_sensitivity(df, out_dir)
    for r in sens:
        print(f"  {r['ridge']:18s}  L={r['L']}  ≥8/10 cells [n=10] = "
              f"{r['cells_geq_8of10_n10']}/{r['L']};  "
              f"≥7/9 cells [no Pat_03] = {r['cells_geq_7of9_n9_no_p03']}/{r['L']};  "
              f"mean ΔVI: n=10 {r['mean_dVI_n10']:+.4f} → no-P03 {r['mean_dVI_n9_no_p03']:+.4f}")

    print("\n=== Check 3 — implant probe analysis ===")
    md = _check3_implant_dissenters(df_means, out_dir)
    # Print head of the summary so it's visible without opening the file
    head = "\n".join(md.splitlines()[:30])
    print(head)


if __name__ == "__main__":
    main()
