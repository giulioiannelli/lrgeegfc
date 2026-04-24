#!/usr/bin/env python3
"""H2-PART consensus — single-score combination of VI/ARI/H/NMI contrasts.

Motivation: cluster-permutation per single measure at n=9 is conservative
(2^9 = 512 sign-flip null). Yet the unanimity landscape shows that four
INDEPENDENT partition-distance measures (VI, ARI, conditional H, NMI)
converge on the same band × k trace pattern. That convergence is
evidence we should leverage with a pooled statistic.

Combination: for each measure, z-score the per-patient contrast within
(band, k) across all patients × k values globally. Then for each
(patient, band, k) take the mean of the four z-scores as the
"consensus contrast". Positive consensus ⇒ all four measures agree
that P^rpost is closer to P^ttest than to P^rpre at that scale.

Cluster-permutation is then applied on the consensus across k per
band, using the same sign-flip null as the single-measure tests. This
is strictly stronger when measures correlate (and they do — all four
measure partition similarity).

Reads:  data/reports/imcoh_vi/h2_partition_multiscale_raw.csv
Writes: data/reports/imcoh_vi/h2_partition_consensus_raw.csv
        data/reports/imcoh_vi/h2_partition_consensus.{md,csv}
        data/reports/imcoh_vi/figures/h2_partition_consensus.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT

from h2d_coactivation_persistence import cluster_perm


CONTRAST_COLS = ["d_VI", "d_ARI", "d_H", "d_NMI"]


def _zscore(x: np.ndarray) -> np.ndarray:
    finite = np.isfinite(x)
    if finite.sum() < 2:
        return np.full_like(x, np.nan, dtype=float)
    mu = np.nanmean(x)
    sd = np.nanstd(x, ddof=1)
    return (x - mu) / sd if sd > 0 else np.zeros_like(x)


def build_consensus(raw: pd.DataFrame) -> pd.DataFrame:
    """Z-score each measure globally, then mean across measures per row."""
    df = raw.copy()
    for c in CONTRAST_COLS:
        df[f"z_{c}"] = _zscore(df[c].to_numpy())
    df["consensus"] = df[[f"z_{c}" for c in CONTRAST_COLS]].mean(axis=1)
    return df


def _band_mat(df: pd.DataFrame, band: str,
              patients: list[str], ks: list[int]) -> np.ndarray:
    bdf = df[df["band"] == band]
    mat = np.full((len(patients), len(ks)), np.nan)
    for ip, pat in enumerate(patients):
        s = bdf[bdf["patient"] == pat].set_index("k")
        for jk, k in enumerate(ks):
            if k in s.index:
                mat[ip, jk] = s.at[k, "consensus"]
    return mat


def _frac_pos(mat: np.ndarray) -> np.ndarray:
    finite = np.isfinite(mat)
    pos = finite & (mat > 0)
    n = finite.sum(axis=0)
    p = pos.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, p / n, np.nan)


BAND_COLORS = {
    "delta":      "#3b6e9c", "theta":      "#c44e4e",
    "alpha":      "#5ea85e", "beta":       "#b28ad1",
    "low_gamma":  "#e5a24b", "high_gamma": "#7d5a50",
}


def main() -> None:
    raw = pd.read_csv(REPORTS_ROOT / "imcoh_vi" / "h2_partition_multiscale_raw.csv")
    cons = build_consensus(raw)
    cons.to_csv(REPORTS_ROOT / "imcoh_vi" / "h2_partition_consensus_raw.csv",
                index=False)

    patients = sorted(cons["patient"].unique())
    ks = sorted(cons["k"].unique())

    # ── Cluster-permutation on consensus per band ────────────
    cluster_rows = []
    band_matrices: dict[str, np.ndarray] = {}
    band_fracs: dict[str, np.ndarray] = {}
    for band in BRAIN_BANDS_NAMES:
        mat = _band_mat(cons, band, patients, ks)
        band_matrices[band] = mat
        band_fracs[band] = _frac_pos(mat)
        keep_row = ~np.isnan(mat).all(axis=1)
        m = mat[keep_row]
        ok = ~np.isnan(m).any(axis=0)
        if not ok.any():
            continue
        clusters = cluster_perm(m[:, ok], np.array(ks)[ok])
        for c in clusters:
            cluster_rows.append({"band": band, "k_start": c["k_start"],
                                  "k_end": c["k_end"], "length": c["len"],
                                  "mass": c["mass"], "p": c["p"]})

    cluster_df = pd.DataFrame(cluster_rows)
    cluster_df.to_csv(
        REPORTS_ROOT / "imcoh_vi" / "h2_partition_consensus.csv", index=False,
    )

    # ── Markdown report ─────────────────────────────────────
    lines: list[str] = []
    ap = lines.append
    ap("# H2-PART consensus — pooled VI/ARI/H/NMI cluster-permutation")
    ap("")
    ap("Each partition-distance contrast is z-scored across the pooled "
       "(patient × band × k) population. The **consensus contrast** per "
       "(patient, band, k) is the mean of the four z-scored measures — "
       "robust to noise in any single measure, amplifies agreement.")
    ap("")
    ap("Sign-flip cluster-permutation across k per band, primary |z|>1.96, "
       "5,000 perms.")
    ap("")
    ap("## Significant consensus clusters (p < 0.05)")
    ap("")
    ap("| band | k-range | length | mass | p (cluster) |")
    ap("|------|:--------|------:|-----:|------------:|")
    sig_found = False
    for _, r in cluster_df.iterrows():
        if r["p"] < 0.05:
            sig_found = True
            ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | "
               f"k={int(r['k_start'])}–{int(r['k_end'])} | "
               f"{int(r['length'])} | {r['mass']:.1f} | {r['p']:.4f}★ |")
    if not sig_found:
        ap("| — | no significant clusters | — | — | — |")
    ap("")
    ap("## All supra-threshold clusters (|z|>1.96, including non-sig.)")
    ap("")
    ap("| band | k-range | length | mass | p (cluster) |")
    ap("|------|:--------|------:|-----:|------------:|")
    for _, r in cluster_df.iterrows():
        star = "★" if r["p"] < 0.05 else ""
        ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | "
           f"k={int(r['k_start'])}–{int(r['k_end'])} | "
           f"{int(r['length'])} | {r['mass']:.1f} | {r['p']:.4f}{star} |")
    if cluster_df.empty:
        ap("| — | no supra-threshold clusters | — | — | — |")
    ap("")

    md = REPORTS_ROOT / "imcoh_vi" / "h2_partition_consensus.md"
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {md}")

    # ── Figure: consensus unanimity landscape + sig-cluster strip ──
    fig, axes = plt.subplots(
        2, 1, figsize=(13.0, 6.0), dpi=160,
        gridspec_kw={"height_ratios": [1.2, 0.6], "hspace": 0.50},
    )

    # Top: per-band consensus unanimity map
    ax = axes[0]
    cmap = LinearSegmentedColormap.from_list(
        "trace",
        [(1.0, 1.0, 1.0), (1.0, 0.95, 0.7),
         (1.0, 0.6, 0.2), (0.7, 0.05, 0.05)],
        N=256,
    )
    frac_mat = np.array([band_fracs[b] for b in BRAIN_BANDS_NAMES])
    im = ax.imshow(frac_mat, aspect="auto", cmap=cmap,
                    vmin=0.0, vmax=1.0,
                    extent=[ks[0] - 0.5, ks[-1] + 0.5,
                            len(BRAIN_BANDS_NAMES) - 0.5, -0.5],
                    interpolation="nearest")
    ys, xs = np.where(frac_mat >= 7/9)
    ax.scatter(np.array(ks)[xs], ys, s=1.2, c="black", alpha=0.5)
    ys, xs = np.where(frac_mat >= 0.999)
    ax.scatter(np.array(ks)[xs], ys, s=5, c="black")
    ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=11)
    ax.set_xlabel("k (dendrogram cut scale)", fontsize=11)
    ax.set_title(
        "(a) Consensus partition-level trace per (band, k) — "
        "mean of z-scored VI, ARI, H, NMI contrasts. "
        "Colour = fraction of patients with positive consensus; "
        "dots: ≥ 7/9 (faint), 9/9 (black).",
        fontsize=10, loc="left",
    )
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, fraction=0.03)
    cbar.set_label("fraction of patients with positive consensus", fontsize=9)

    # Bottom: cluster-permutation strip
    ax = axes[1]
    ax.set_xlim(ks[0] - 0.5, ks[-1] + 0.5)
    ax.set_ylim(-0.5, len(BRAIN_BANDS_NAMES) - 0.5)
    ax.invert_yaxis()
    ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=11)
    ax.set_xlabel("k (dendrogram cut scale)", fontsize=11)
    ax.grid(axis="x", alpha=0.3)
    for _, r in cluster_df.iterrows():
        band = r["band"]
        ib = BRAIN_BANDS_NAMES.index(band)
        alpha = 0.85 if r["p"] < 0.05 else 0.35
        edgecolor = "black" if r["p"] < 0.05 else BAND_COLORS[band]
        lw = 1.3 if r["p"] < 0.05 else 0.5
        rect = Rectangle(
            (r["k_start"] - 0.5, ib - 0.38),
            r["k_end"] - r["k_start"] + 1, 0.76,
            facecolor=BAND_COLORS[band], alpha=alpha,
            edgecolor=edgecolor, linewidth=lw,
        )
        ax.add_patch(rect)
        if r["p"] < 0.05:
            ax.text(
                (r["k_start"] + r["k_end"]) / 2, ib,
                f"p={r['p']:.3f}★", ha="center", va="center",
                fontsize=8.5, fontweight="bold",
            )
    ax.set_title(
        "(b) Sign-flip cluster-permutation on the consensus score "
        "(5,000 perms, primary |z|>1.96).  "
        "Filled: supra-threshold clusters; outlined+text: p<0.05.",
        fontsize=10, loc="left",
    )

    fig.suptitle(
        "Consensus partition-level task trace per band (n=9).  "
        "Pooling VI + ARI + H + NMI leverages cross-measure agreement.",
        fontsize=12, y=1.00,
    )
    for ext in ("pdf", "png"):
        out = REPORTS_ROOT / "imcoh_vi" / "figures" / f"h2_partition_consensus.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)

    # Summary to console
    if cluster_df.empty:
        print("\nNo supra-threshold consensus clusters at |z|>1.96.")
    else:
        sig = cluster_df[cluster_df["p"] < 0.05]
        print(f"\nSignificant consensus clusters: {len(sig)}")
        for _, r in sig.iterrows():
            print(f"  {r['band']}  k={int(r['k_start'])}–{int(r['k_end'])}  "
                  f"L={int(r['length'])}  p={r['p']:.4f}")


if __name__ == "__main__":
    main()
