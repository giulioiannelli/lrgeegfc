#!/usr/bin/env python3
"""Diagnostic — does Spearman ρ on D matrices track tree-topology overlap?

For every patient × band × phase pair, compute BOTH:

  * ρ = Spearman(upper-tri(D^p1), upper-tri(D^p2))          (what we use)
  * bip_overlap = |bipartitions(Z^p1) ∩ bipartitions(Z^p2)| / (N-1)
    (what REAL tree topology similarity looks like: count of internal
     bipartitions that are identical sets-of-leaves in both trees)

Then show: are they correlated? Do bands or phase-pairs with high ρ
also have high bip_overlap, or does ρ tell a different story than
topology?

Also compute cross-phase averages per band, and split "within-type"
(rpre↔rpost, tlearn↔ttest) vs "cross-type" pairs.

Reads:  LRG imcoh_abs caches.
Writes: data/reports/imcoh_vi/diag_topology_vs_rho.{csv,md}
        data/reports/imcoh_vi/figures/diag_topology_vs_rho.{pdf,png}
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result


PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
TYPE = {"rest_pre": "rest", "rest_post": "rest",
        "task_learn": "task", "task_test": "task"}

BAND_COLORS = {
    "delta":      "#3b6e9c",
    "theta":      "#c44e4e",
    "alpha":      "#5ea85e",
    "beta":       "#b28ad1",
    "low_gamma":  "#e5a24b",
    "high_gamma": "#7d5a50",
}


def _upper_tri(D: np.ndarray) -> np.ndarray:
    D = np.asarray(D)
    return D if D.ndim == 1 else D[np.triu_indices(D.shape[0], k=1)]


def _bipartitions(Z: np.ndarray, n: int) -> set[frozenset[int]]:
    members: dict[int, set[int]] = {i: {i} for i in range(n)}
    out: set[frozenset[int]] = set()
    for i, (a, b, *_) in enumerate(Z):
        a, b = int(a), int(b)
        node_id = n + i
        members[node_id] = members[a] | members[b]
        out.add(frozenset(members[node_id]))
    return out


def main() -> None:
    rows = []
    for pat in PATIENTS_LIST:
        # Load all 4 phases' Z and D for each band
        for band in BRAIN_BANDS_NAMES:
            per_phase: dict[str, tuple[np.ndarray, np.ndarray, int]] = {}
            for phase in PHASES:
                r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
                if r is None:
                    continue
                Z = np.asarray(r.linkage_matrix)
                D = _upper_tri(r.ultrametric_matrix)
                N = Z.shape[0] + 1
                per_phase[phase] = (Z, D, N)

            for p1, p2 in combinations(per_phase, 2):
                Z1, D1, N1 = per_phase[p1]
                Z2, D2, N2 = per_phase[p2]
                if N1 != N2 or D1.shape != D2.shape:
                    continue
                rho, _ = stats.spearmanr(D1, D2)
                bps1 = _bipartitions(Z1, N1)
                bps2 = _bipartitions(Z2, N2)
                overlap = len(bps1 & bps2) / (N1 - 1)
                same_type = TYPE[p1] == TYPE[p2]
                rows.append({
                    "patient": pat, "band": band,
                    "p1": p1, "p2": p2,
                    "rho": float(rho),
                    "bip_overlap": float(overlap),
                    "pair_type": "within" if same_type else "cross",
                })

    df = pd.DataFrame(rows)
    out_dir = REPORTS_ROOT / "imcoh_vi"
    df.to_csv(out_dir / "diag_topology_vs_rho.csv", index=False)

    # ── Figure: two-panel diagnostic ─────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.5), dpi=150,
                              gridspec_kw={"width_ratios": [1.2, 1.0]})

    # Left: scatter ρ vs bipartition overlap, coloured by band, marker by type
    ax = axes[0]
    for band in BRAIN_BANDS_NAMES:
        sub = df[df["band"] == band]
        ax.scatter(
            sub[sub["pair_type"] == "within"]["rho"],
            sub[sub["pair_type"] == "within"]["bip_overlap"],
            marker="o", s=28, color=BAND_COLORS[band],
            edgecolor="black", linewidth=0.3, alpha=0.75,
            label=f"{BRAIN_BAND_TEX_DICT[band]} (within)" if band == "delta" else None,
        )
        ax.scatter(
            sub[sub["pair_type"] == "cross"]["rho"],
            sub[sub["pair_type"] == "cross"]["bip_overlap"],
            marker="x", s=28, color=BAND_COLORS[band],
            alpha=0.75,
            label=f"{BRAIN_BAND_TEX_DICT[band]} (cross)" if band == "delta" else None,
        )
    ax.plot([0, 1], [0, 1], color="#888", lw=0.6, linestyle="--")
    # Correlation of ρ vs bip_overlap across all rows
    r_all, p_all = stats.spearmanr(df["rho"], df["bip_overlap"])
    ax.set_xlabel("Spearman ρ on ultrametric distance matrices", fontsize=11)
    ax.set_ylabel("bipartition overlap |bps(Z1) ∩ bps(Z2)| / (N−1)", fontsize=11)
    ax.set_title(f"(a) ρ vs true tree bipartition overlap  "
                  f"(Spearman across all pairs = {r_all:+.3f}, p={p_all:.1e})",
                  fontsize=11, loc="left")
    ax.set_xlim(-0.2, 1.0)
    ax.set_ylim(0.0, 1.0)
    # Simple legend: dots = same-type phase pair, x = cross-type
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", color="grey", linestyle="none",
                markerfacecolor="grey", markersize=7, label="within-type pair"),
        Line2D([0], [0], marker="x", color="grey", linestyle="none",
                markersize=7, label="cross-type pair"),
    ]
    for band in BRAIN_BANDS_NAMES:
        legend_handles.append(
            Line2D([0], [0], marker="s", color=BAND_COLORS[band],
                    linestyle="none", markersize=7,
                    label=BRAIN_BAND_TEX_DICT[band])
        )
    ax.legend(handles=legend_handles, frameon=False, fontsize=8.5,
               loc="lower right", ncol=2)
    ax.grid(alpha=0.3)

    # Right: mean bip_overlap per band, within vs cross
    ax = axes[1]
    band_means = (df.groupby(["band", "pair_type"])["bip_overlap"].mean()
                  .unstack("pair_type"))
    band_means = band_means.reindex(index=BRAIN_BANDS_NAMES)
    x = np.arange(len(BRAIN_BANDS_NAMES))
    w = 0.36
    ax.bar(x - w/2, band_means["within"].values, w,
           color=[BAND_COLORS[b] for b in BRAIN_BANDS_NAMES],
           edgecolor="black", linewidth=0.5, label="within-type")
    ax.bar(x + w/2, band_means["cross"].values, w,
           color=[BAND_COLORS[b] for b in BRAIN_BANDS_NAMES],
           edgecolor="black", linewidth=0.5, alpha=0.55, label="cross-type")
    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                       fontsize=12)
    ax.set_ylabel("mean bipartition overlap", fontsize=11)
    ax.set_title("(b) Tree-topology similarity per band, within vs cross phase-type",
                  fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, 0.7)

    # Print annotations
    for xi, b in enumerate(BRAIN_BANDS_NAMES):
        wv = band_means.at[b, "within"]
        cv = band_means.at[b, "cross"]
        ax.text(xi, max(wv, cv) + 0.02, f"Δ={wv-cv:+.3f}",
                ha="center", va="bottom", fontsize=8.5)

    fig.tight_layout()
    out_png = out_dir / "figures" / "diag_topology_vs_rho.png"
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    print(f"saved {out_png}")
    plt.close(fig)

    # Report
    lines: list[str] = []
    ap = lines.append
    ap("# Diagnostic — Spearman ρ vs tree-bipartition overlap")
    ap("")
    ap(f"All pairs: Spearman ρ vs bipartition overlap = **{r_all:+.3f}** "
       f"(p = {p_all:.2e}).")
    ap("")
    ap("Interpretation: if ρ tracks topology, this correlation is high. "
       "If ρ is arithmetic-dominated, it could be near zero.")
    ap("")
    ap("## Mean bipartition overlap per band × pair-type")
    ap("")
    ap("| band | within-type | cross-type | Δ within−cross |")
    ap("|------|------------:|-----------:|---------------:|")
    for b in BRAIN_BANDS_NAMES:
        wv = band_means.at[b, "within"]
        cv = band_means.at[b, "cross"]
        ap(f"| {BRAIN_BAND_TEX_DICT[b]} | {wv:.3f} | {cv:.3f} | {wv-cv:+.3f} |")
    ap("")
    md = out_dir / "diag_topology_vs_rho.md"
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")
    print(f"\nGlobal Spearman(ρ, bip_overlap) = {r_all:+.3f}  (p = {p_all:.2e})")


if __name__ == "__main__":
    main()
