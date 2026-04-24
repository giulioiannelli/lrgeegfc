#!/usr/bin/env python3
"""H2-PART cluster-permutation — formalize per-band per-k significance.

For each (band, measure ∈ {VI, ARI, H, NMI}), run sign-flip cluster
permutation across k on the per-patient contrast values (from
`h2_partition_multiscale_raw.csv`). Identifies scale-coherent k ranges
per band where the trace is statistically significant beyond chance.

Reuses `cluster_perm` from `h2d_coactivation_persistence.py`.

Reads:  data/reports/imcoh_vi/h2_partition_multiscale_raw.csv
Writes: data/reports/imcoh_vi/h2_partition_cluster_perm.md
        data/reports/imcoh_vi/h2_partition_cluster_perm.csv
        data/reports/imcoh_vi/figures/h2_partition_cluster_perm.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT

from h2d_coactivation_persistence import cluster_perm


CONTRASTS = [
    ("d_VI",  r"$\Delta_{VI}$"),
    ("d_ARI", r"$\Delta_{ARI}$"),
    ("d_H",   r"$\Delta_H$"),
    ("d_NMI", r"$\Delta_{NMI}$"),
]

BAND_COLORS = {
    "delta":      "#3b6e9c", "theta":      "#c44e4e",
    "alpha":      "#5ea85e", "beta":       "#b28ad1",
    "low_gamma":  "#e5a24b", "high_gamma": "#7d5a50",
}


def _build_mat(sub: pd.DataFrame, col: str,
                patients: list[str], ks: list[int]) -> np.ndarray:
    mat = np.full((len(patients), len(ks)), np.nan)
    for ip, pat in enumerate(patients):
        s = sub[sub["patient"] == pat].set_index("k")
        for jk, k in enumerate(ks):
            if k in s.index:
                mat[ip, jk] = s.at[k, col]
    return mat


def main() -> None:
    raw = pd.read_csv(REPORTS_ROOT / "imcoh_vi" / "h2_partition_multiscale_raw.csv")
    patients = sorted(raw["patient"].unique())
    ks = sorted(raw["k"].unique())

    records = []
    all_clusters: dict[tuple[str, str], list[dict]] = {}
    for col, label in CONTRASTS:
        for band in BRAIN_BANDS_NAMES:
            sub = raw[raw["band"] == band]
            mat = _build_mat(sub, col, patients, ks)
            keep_row = ~np.isnan(mat).all(axis=1)
            mat = mat[keep_row]
            ok_col = ~np.isnan(mat).any(axis=0)
            if not ok_col.any():
                continue
            clusters = cluster_perm(mat[:, ok_col], np.array(ks)[ok_col])
            all_clusters[(col, band)] = clusters
            for c in clusters:
                records.append({
                    "measure": col, "band": band,
                    "k_start": c["k_start"], "k_end": c["k_end"],
                    "length": c["len"], "mass": c["mass"], "p": c["p"],
                })

    df = pd.DataFrame(records)
    out_dir = REPORTS_ROOT / "imcoh_vi"
    df.to_csv(out_dir / "h2_partition_cluster_perm.csv", index=False)

    # ── Markdown report ─────────────────────────────────────
    lines: list[str] = []
    ap = lines.append
    ap("# H2-PART cluster-permutation — significant k-ranges per band × measure")
    ap("")
    ap("Sign-flip cluster permutation (5,000 perms, primary |z|>1.96) on "
       "per-patient contrasts across k ∈ [2, 49], per band per measure.")
    ap("")
    for col, label in CONTRASTS:
        ap(f"## {label}")
        ap("")
        ap("| band | k-range | length | mass | p (cluster) |")
        ap("|------|:--------|------:|-----:|------------:|")
        any_found = False
        for band in BRAIN_BANDS_NAMES:
            cs = all_clusters.get((col, band), [])
            for c in cs:
                any_found = True
                sig = "★" if c["p"] < 0.05 else ""
                ap(f"| {BRAIN_BAND_TEX_DICT[band]} | k={c['k_start']}–{c['k_end']} | "
                   f"{c['len']} | {c['mass']:.1f} | {c['p']:.4f}{sig} |")
        if not any_found:
            ap("| — | no supra-threshold clusters | — | — | — |")
        ap("")

    md = out_dir / "h2_partition_cluster_perm.md"
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")
    print(f"wrote {out_dir / 'h2_partition_cluster_perm.csv'}")

    # ── Figure: band × k significant-cluster strip per measure ──
    fig, axes = plt.subplots(4, 1, figsize=(12.5, 7.5), dpi=160, sharex=True)
    for ax, (col, label) in zip(axes, CONTRASTS):
        ax.set_xlim(ks[0] - 0.5, ks[-1] + 0.5)
        ax.set_ylim(-0.5, len(BRAIN_BANDS_NAMES) - 0.5)
        ax.invert_yaxis()
        ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                            fontsize=11)
        ax.set_ylabel(label, fontsize=11, rotation=0, ha="right",
                       va="center", labelpad=20)
        ax.grid(axis="x", alpha=0.3)
        for ib, band in enumerate(BRAIN_BANDS_NAMES):
            for c in all_clusters.get((col, band), []):
                alpha = 0.85 if c["p"] < 0.05 else 0.30
                edgecolor = "black" if c["p"] < 0.05 else BAND_COLORS[band]
                lw = 1.2 if c["p"] < 0.05 else 0.4
                rect = Rectangle(
                    (c["k_start"] - 0.5, ib - 0.38),
                    c["k_end"] - c["k_start"] + 1, 0.76,
                    facecolor=BAND_COLORS[band], alpha=alpha,
                    edgecolor=edgecolor, linewidth=lw,
                )
                ax.add_patch(rect)
                if c["p"] < 0.05:
                    ax.text(
                        (c["k_start"] + c["k_end"]) / 2, ib,
                        f"p={c['p']:.3f}★", ha="center", va="center",
                        fontsize=8, fontweight="bold",
                    )
    axes[-1].set_xlabel("k (dendrogram cut scale)", fontsize=11)
    fig.suptitle(
        "Cluster-permutation on partition-level trace contrasts per (band, k).\n"
        "Filled: supra-threshold clusters (|z|>1.96). Outlined+text: p<0.05. "
        "Transparent: observed but non-significant under permutation.",
        fontsize=11, y=1.00,
    )
    fig.tight_layout()
    out_png = out_dir / "figures" / "h2_partition_cluster_perm.png"
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    print(f"saved {out_png}")
    plt.close(fig)


if __name__ == "__main__":
    main()
