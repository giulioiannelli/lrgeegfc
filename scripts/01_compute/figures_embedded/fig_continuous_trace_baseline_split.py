#!/usr/bin/env python3
"""Overlay figure: Run A common-baseline control.

Compares per-(patient, band) ρ from two modes of the continuous-trace
matrix:

  shared          Δ_task = D_test − D_pre,    Δ_rest = D_post − D_pre
                  (single shared baseline; original mode)

  split-baseline  Δ_task = D_test − D_pre_A,  Δ_rest = D_post − D_pre_B
                  (independent half-baselines; D_pre_A ⊥ D_pre_B)

The shared mode contains a spurious-correlation contribution from
shared baseline noise. Run A removes it. The drop ρ_shared − ρ_split
quantifies the artefact magnitude per (patient, band); the residual
ρ_split is the noise-floor-aware estimate of the genuine task-induced
trace.

Reads:
    data/reports/imcoh_continuous_trace/per_cell_summary.csv
    data/reports/imcoh_continuous_trace/per_cell_summary_split.csv
Writes:
    data/reports/imcoh_continuous_trace/figures/baseline_split_overlay.pdf
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT


IN_DIR = REPORTS_ROOT / "imcoh_continuous_trace"
OUT_DIR = IN_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


BANDS = list(BRAIN_BANDS_NAMES)
TEX = BRAIN_BAND_TEX_DICT


def render() -> None:
    shared = pd.read_csv(IN_DIR / "per_cell_summary.csv")
    split = pd.read_csv(IN_DIR / "per_cell_summary_split.csv")
    merged = shared.merge(
        split[["patient", "band", "rho"]].rename(columns={"rho": "rho_split"}),
        on=["patient", "band"], how="inner",
    ).rename(columns={"rho": "rho_shared"})

    fig = plt.figure(figsize=(15.0, 5.0), dpi=150)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.6, 1.0, 1.4])
    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[0, 2])

    # --- Panel A: per-patient paired ρ across bands ---
    rng = np.random.default_rng(0)
    for j, band in enumerate(BANDS):
        sub = merged[merged["band"] == band]
        x_shared = j - 0.18 + rng.uniform(-0.06, 0.06, size=len(sub))
        x_split  = j + 0.18 + rng.uniform(-0.06, 0.06, size=len(sub))
        ax_A.scatter(x_shared, sub["rho_shared"], s=34, color="#1f77b4",
                     alpha=0.85, edgecolor="black", linewidth=0.3, zorder=3,
                     label="shared D_pre" if j == 0 else None)
        ax_A.scatter(x_split, sub["rho_split"], s=34, color="#d62728",
                     alpha=0.85, edgecolor="black", linewidth=0.3, zorder=3,
                     label=r"split D_pre (Run A)" if j == 0 else None)
        for xs, xp, rs, rp in zip(x_shared, x_split,
                                   sub["rho_shared"], sub["rho_split"]):
            ax_A.plot([xs, xp], [rs, rp], color="#999", lw=0.4, alpha=0.5,
                      zorder=2)
        med_s = sub["rho_shared"].median()
        med_p = sub["rho_split"].median()
        ax_A.hlines(med_s, j - 0.30, j - 0.06, color="#1f77b4", lw=2.4,
                    zorder=4)
        ax_A.hlines(med_p, j + 0.06, j + 0.30, color="#d62728", lw=2.4,
                    zorder=4)
    ax_A.axhline(0, color="black", lw=0.7)
    ax_A.set_xticks(range(len(BANDS)))
    ax_A.set_xticklabels([TEX[b] for b in BANDS], fontsize=11)
    ax_A.set_xlim(-0.5, len(BANDS) - 0.5)
    ax_A.set_ylim(-0.5, 1.0)
    ax_A.set_ylabel(r"$\rho = \mathrm{Spearman}(\Delta_{\mathrm{task}}, "
                    r"\Delta_{\mathrm{rest}})$", fontsize=11)
    ax_A.set_xlabel("band", fontsize=11)
    ax_A.legend(loc="lower right", frameon=False, fontsize=9)
    ax_A.set_title("A — Per-patient ρ: shared vs split-baseline", fontsize=11)

    # --- Panel B: scatter rho_shared vs rho_split (all 60 cells) ---
    cmap = plt.get_cmap("viridis")
    band_to_color = {b: cmap(j / max(1, len(BANDS) - 1))
                     for j, b in enumerate(BANDS)}
    for band in BANDS:
        sub = merged[merged["band"] == band]
        ax_B.scatter(sub["rho_shared"], sub["rho_split"],
                     s=46, color=band_to_color[band], alpha=0.85,
                     edgecolor="black", linewidth=0.3, label=TEX[band])
    lo = min(merged["rho_shared"].min(), merged["rho_split"].min()) - 0.05
    hi = max(merged["rho_shared"].max(), merged["rho_split"].max()) + 0.05
    ax_B.plot([lo, hi], [lo, hi], color="#444", lw=0.8, ls="--", zorder=1)
    ax_B.axhline(0, color="#888", lw=0.5)
    ax_B.axvline(0, color="#888", lw=0.5)
    ax_B.set_xlim(lo, hi); ax_B.set_ylim(lo, hi)
    ax_B.set_aspect("equal")
    ax_B.set_xlabel(r"$\rho_{\mathrm{shared}}$", fontsize=11)
    ax_B.set_ylabel(r"$\rho_{\mathrm{split}}$", fontsize=11)
    ax_B.set_title("B — Per-cell shift (60 cells)", fontsize=11)
    ax_B.legend(loc="lower right", frameon=False, fontsize=8)

    # --- Panel C: per-band cohort medians + IQR comparison ---
    xs = np.arange(len(BANDS))
    med_s = [merged[merged["band"] == b]["rho_shared"].median() for b in BANDS]
    med_p = [merged[merged["band"] == b]["rho_split"].median()  for b in BANDS]
    q1_s = [merged[merged["band"] == b]["rho_shared"].quantile(0.25) for b in BANDS]
    q3_s = [merged[merged["band"] == b]["rho_shared"].quantile(0.75) for b in BANDS]
    q1_p = [merged[merged["band"] == b]["rho_split"].quantile(0.25)  for b in BANDS]
    q3_p = [merged[merged["band"] == b]["rho_split"].quantile(0.75)  for b in BANDS]
    npos_s = [int((merged[merged["band"] == b]["rho_shared"] > 0).sum())
              for b in BANDS]
    npos_p = [int((merged[merged["band"] == b]["rho_split"]  > 0).sum())
              for b in BANDS]

    bar_w = 0.36
    ax_C.bar(xs - bar_w/2, med_s, bar_w, color="#1f77b4", alpha=0.85,
             edgecolor="black", linewidth=0.3, label="shared")
    ax_C.bar(xs + bar_w/2, med_p, bar_w, color="#d62728", alpha=0.85,
             edgecolor="black", linewidth=0.3, label="split (Run A)")
    ax_C.errorbar(xs - bar_w/2, med_s,
                  yerr=[np.array(med_s) - np.array(q1_s),
                        np.array(q3_s) - np.array(med_s)],
                  fmt="none", ecolor="#1f3f6f", lw=1.2, capsize=3, zorder=4)
    ax_C.errorbar(xs + bar_w/2, med_p,
                  yerr=[np.array(med_p) - np.array(q1_p),
                        np.array(q3_p) - np.array(med_p)],
                  fmt="none", ecolor="#5e1818", lw=1.2, capsize=3, zorder=4)
    for j, (b, ns, np_) in enumerate(zip(BANDS, npos_s, npos_p)):
        ax_C.text(j - bar_w/2, max(med_s[j], 0) + 0.04, f"{ns}/10",
                  ha="center", fontsize=8, color="#1f3f6f", fontweight="bold")
        ax_C.text(j + bar_w/2, max(med_p[j], 0) + 0.04, f"{np_}/10",
                  ha="center", fontsize=8, color="#5e1818", fontweight="bold")
    ax_C.axhline(0, color="black", lw=0.7)
    ax_C.set_xticks(xs)
    ax_C.set_xticklabels([TEX[b] for b in BANDS], fontsize=11)
    ax_C.set_ylim(-0.15, 0.85)
    ax_C.set_ylabel(r"cohort median $\rho$ (IQR bars)", fontsize=11)
    ax_C.set_xlabel("band", fontsize=11)
    ax_C.legend(loc="upper right", frameon=False, fontsize=9)
    ax_C.set_title("C — Cohort summary (numbers above bars: ρ>0 / 10)",
                   fontsize=11)

    fig.tight_layout()
    out = OUT_DIR / "baseline_split_overlay.pdf"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")
    plt.close(fig)

    # numeric summary written next to the figure
    summary_lines = ["# Run A — common-baseline control summary",
                     "",
                     "| band | n | shared med | split med | drop | "
                     "shared n_+ | split n_+ |",
                     "|------|--:|-----------:|----------:|-----:|"
                     "----------:|----------:|"]
    for b, ms, mp, ns, np_ in zip(BANDS, med_s, med_p, npos_s, npos_p):
        summary_lines.append(
            f"| {TEX[b]} | 10 | {ms:+.3f} | {mp:+.3f} | "
            f"{ms - mp:+.3f} | {ns}/10 | {np_}/10 |"
        )
    summary_lines += [
        "",
        "* `shared` = Δ_task and Δ_rest both reference D_pre (the same "
        "full-phase baseline). Vulnerable to spurious correlation through "
        "shared D_pre noise.",
        "* `split` = Run A. Δ_task − D_pre_A, Δ_rest − D_pre_B with "
        "D_pre_A ⊥ D_pre_B (independent half-of-rest_pre baselines from "
        "data/cache/imcoh_lrg_halves/). The shared-baseline-noise "
        "contribution to ρ is removed by construction.",
        "* `drop` = ρ_shared − ρ_split, the spurious-correlation magnitude.",
        "* The half-data D_pre_A/B uses halved Welch nperseg → noisier "
        "baseline → modestly inflates Var(Δ_task), Var(Δ_rest) and so "
        "deflates ρ_split slightly compared to the underlying population "
        "correlation. The right comparison is split vs the H2e split-half "
        "null (which lives in the same noise regime).",
        "",
    ]
    summary_md = OUT_DIR.parent / "run_A_baseline_split_summary.md"
    summary_md.write_text("\n".join(summary_lines), encoding="utf-8")
    print(f"saved {summary_md}")


if __name__ == "__main__":
    render()
