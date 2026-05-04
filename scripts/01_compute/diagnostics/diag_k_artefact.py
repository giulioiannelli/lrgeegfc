#!/usr/bin/env python3
"""k-artefact companion diagnostic to the n=10 IMCOH_ABS headline figure.

Per (patient, band, phase) — across k ∈ [2, 49] — compute three diagnostic
curves that surface the spurious extremes the user flagged:

  - n_eff(k)               = Simpson effective number of clusters (1/Σpᵢ²)
  - singleton_fraction(k)  = fraction of clusters of size 1
  - max_cluster_fraction(k)= largest cluster size / N

At k ≈ N the partition is singleton-dominated; at k ≈ O(1) one giant cluster
dominates. Both regimes deserve transparent visual treatment instead of
silent inclusion in the headline. We do NOT mask the headline figure (per
user choice 2026-04-25); this companion lets the reader see directly where
on the k axis those regimes live.

Inputs:  Pat_NN LRG imcoh_abs caches (4 phases × 6 bands × 10 patients).
Outputs:
  data/outputs/figures/section6/k_artefact_diagnostic_n10_imcoh_abs.pdf
  data/outputs/figures/section6/k_artefact_diagnostic_n10_imcoh_abs.md (sidecar)
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PHASE_LABELS,
)
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.utils.metrics.tree import simpson_neff, cluster_size_stats
from lrg_eegfc.workflow.lrg import load_lrg_result


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
K_RANGE = np.arange(2, 50)
PHASE_COLORS = {
    "rest_pre":  "#2c7bb6",
    "task_learn": "#fdae61",
    "task_test":  "#d7191c",
    "rest_post": "#1a9641",
}
PAT_LINE_ALPHA = 0.45


def _load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return np.asarray(r.linkage_matrix) if r is not None else None


def _curves_for(Z: np.ndarray) -> dict:
    """For a single linkage matrix, compute (n_eff, singleton, max_frac) over K_RANGE."""
    if Z is None:
        return {}
    n_eff = np.full(K_RANGE.size, np.nan)
    sfrac = np.full(K_RANGE.size, np.nan)
    mfrac = np.full(K_RANGE.size, np.nan)
    for jk, k in enumerate(K_RANGE):
        labels = fcluster(Z, k, criterion="maxclust")
        st = cluster_size_stats(labels)
        n_eff[jk] = st["n_eff"]
        sfrac[jk] = st["singleton_fraction"]
        mfrac[jk] = st["max_cluster_fraction"]
    return {"n_eff": n_eff, "singleton": sfrac, "max_frac": mfrac}


def main() -> None:
    out_dir = FIGURES_ROOT / "section6"
    out_dir.mkdir(parents=True, exist_ok=True)

    # data[pat][band][phase] = {n_eff, singleton, max_frac}
    data: dict = {}
    for pat in COHORT_N10:
        data[pat] = {}
        for band in BRAIN_BANDS_NAMES:
            data[pat][band] = {}
            for phase in PHASE_LABELS:
                Z = _load_Z(pat, phase, band)
                if Z is None:
                    print(f"  SKIP missing: {pat} {band} {phase}")
                    continue
                data[pat][band][phase] = _curves_for(Z)

    # Figure: 6 bands × 3 metric columns
    fig, axes = plt.subplots(
        len(BRAIN_BANDS_NAMES), 3,
        figsize=(13.0, 14.0), dpi=160,
        sharex=True,
        gridspec_kw={"hspace": 0.32, "wspace": 0.20},
    )

    metric_cols = [
        ("n_eff",     r"$n_{\mathrm{eff}}(k)$ — Simpson"),
        ("singleton", r"singleton fraction"),
        ("max_frac",  r"max cluster / $N$"),
    ]

    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        for jc, (mkey, mlabel) in enumerate(metric_cols):
            ax = axes[ib, jc]
            for phase in PHASE_LABELS:
                colour = PHASE_COLORS[phase]
                # patient curves (faint), then median curve (solid)
                stack = []
                for pat in COHORT_N10:
                    curves = data[pat][band].get(phase, None)
                    if curves is None or mkey not in curves:
                        continue
                    y = curves[mkey]
                    ax.plot(K_RANGE, y, color=colour,
                            alpha=PAT_LINE_ALPHA, linewidth=0.6)
                    stack.append(y)
                if stack:
                    arr = np.vstack(stack)
                    med = np.nanmedian(arr, axis=0)
                    ax.plot(K_RANGE, med, color=colour, linewidth=2.0,
                            label=phase if (ib == 0 and jc == 0) else None)
            if jc == 0:
                ax.set_ylabel(BRAIN_BAND_TEX_DICT[band], fontsize=10,
                              rotation=0, ha="right", va="center", labelpad=10)
            if ib == 0:
                ax.set_title(mlabel, fontsize=11)
            if ib == len(BRAIN_BANDS_NAMES) - 1:
                ax.set_xlabel("k", fontsize=10)
            ax.grid(alpha=0.25, linestyle=":")
            # Reference lines
            if mkey == "singleton":
                ax.axhline(0.5, color="0.5", linewidth=0.6, linestyle="--")
            elif mkey == "max_frac":
                ax.axhline(0.7, color="0.5", linewidth=0.6, linestyle="--")

    # Legend on the very first axis
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels,
                   loc="upper center", ncol=4, frameon=False, fontsize=10,
                   bbox_to_anchor=(0.5, 1.005))

    out_pdf = out_dir / "k_artefact_diagnostic_n10_imcoh_abs.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_pdf}")

    md = out_dir / "k_artefact_diagnostic_n10_imcoh_abs.md"
    md.write_text(
        "---\n"
        "name: k-artefact-diagnostic-n10\n"
        "type: figure-sidecar\n"
        "era: COHORT_N10\n"
        "status: current\n"
        "created: 2026-04-25\n"
        "updated: 2026-04-26\n"
        "pointers:\n"
        "  - .agents/reports/2026-04-25_task-trace-audit-and-recovery.md\n"
        "  - data/outputs/figures/section6/task_trace_band_k_n10_imcoh_abs.pdf\n"
        "---\n\n"
        "# k-artefact companion diagnostic (n=10, IMCOH_ABS)\n\n"
        "Per (patient, band, phase) curves of three k-axis diagnostics:\n"
        "`n_eff(k)` (Simpson), `singleton_fraction(k)`, `max_cluster_fraction(k)`.\n"
        "Faint lines = individual patients; solid line = phase median.\n"
        "Phase colors: `rest_pre`=blue, `task_learn`=orange, `task_test`=red,\n"
        "`rest_post`=green.\n\n"
        "## Reading rules\n\n"
        "- **Where `singleton_fraction(k) ≥ 0.5`** the partition is\n"
        "  singleton-dominated. VI(k) and partition-divergence statistics\n"
        "  computed at those k are noise-dominated.\n"
        "- **Where `max_cluster_fraction(k) ≥ 0.7`** one giant cluster\n"
        "  swallows the partition. Statistics there are insensitive.\n"
        "- **`n_eff(k)`** falls below ~5 in both regimes — useful as a\n"
        "  single-glance check.\n\n"
        "## Why this is a companion, not a mask\n\n"
        "User explicitly chose 2026-04-25 to NOT mask the headline figure's k\n"
        "extremes. This diagnostic surfaces the artefact transparently for the\n"
        "reader instead. Cross-check: any cell with `≥ 8/10` unanimity in the\n"
        "headline that lies inside a singleton-dominated or giant-dominated\n"
        "k-region should be treated as fragile.\n",
        encoding="utf-8",
    )
    print(f"wrote {md}")


if __name__ == "__main__":
    main()
