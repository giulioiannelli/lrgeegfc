#!/usr/bin/env python3
"""H2-CBR — cluster-birth-retention on rest_post internal tree nodes.

Primary new partition-level directed test. For every internal node η of
the rest_post LRG dendrogram, measure how well its leaf set matches a
subtree in task_test vs rest_pre:

    match_tt(η)   = max_{η'∈tree(Z^tt)}   |L(η)∩L(η')| / |L(η)∪L(η')|
    match_rpre(η) = max_{η'∈tree(Z^rpre)} |L(η)∩L(η')| / |L(η)∪L(η')|
    δ(η) = match_tt(η) − match_rpre(η)

A module with δ(η) > 0 is "present in rest_post, matches task, and did not
exist in rest_pre" — the task-born-retained object the paper claims. Size-
weighted across every internal node of rest_post:

    S(pat, band) = Σ_η |L(η)| · δ(η)  /  Σ_η |L(η)|

Primary statistic: S > 0 per band (one-sample Wilcoxon, FDR-BH m=6).
Scale-localization: S(pat, band, h) binned on a log grid of h_rel, sign-
flip cluster-permutation across h (reuses `cluster_perm`).

Reads:  LRG imcoh_abs caches for 3 phases.
Writes: data/reports/imcoh_vi/h2_cluster_birth_retention_raw.csv (per-node)
        data/reports/imcoh_vi/h2_cluster_birth_retention_patient.csv
        data/reports/imcoh_vi/h2_cluster_birth_retention.{md,csv}
        data/reports/imcoh_vi/figures/h2_cluster_birth_retention.{pdf,png}
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST,
)
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.utils.metrics import (
    tree_internal_nodes, jaccard_leafsets, h_log_grid,
)
from lrg_eegfc.workflow.lrg import load_lrg_result

from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr
from h2d_coactivation_persistence import cluster_perm


N_H = 60
H_MIN_FLOOR = 1e-3


# ──────────────────────────── data loaders ────────────────────────────

def _load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return np.asarray(r.linkage_matrix) if r is not None else None


# ──────────────────────────── per-node scoring ────────────────────────────

def _best_jaccard(leaves: frozenset[int], ref_nodes: list[dict]) -> float:
    """Max Jaccard of ``leaves`` against every internal subtree in ``ref_nodes``."""
    if not ref_nodes:
        return 0.0
    best = 0.0
    for rn in ref_nodes:
        j = jaccard_leafsets(leaves, rn["leaves"])
        if j > best:
            best = j
            if best == 1.0:
                break
    return best


def score_rpost_nodes(Z_rpre: np.ndarray, Z_tt: np.ndarray,
                      Z_rpost: np.ndarray) -> pd.DataFrame:
    """Per-node rows: h, h_rel, size, match_tt, match_rpre, delta."""
    rpost_nodes = tree_internal_nodes(Z_rpost)
    tt_nodes = tree_internal_nodes(Z_tt)
    rpre_nodes = tree_internal_nodes(Z_rpre)
    rows = []
    for η in rpost_nodes:
        m_tt = _best_jaccard(η["leaves"], tt_nodes)
        m_rpre = _best_jaccard(η["leaves"], rpre_nodes)
        rows.append({
            "h_rpost": η["h"],
            "h_rel_rpost": η["h_rel"],
            "size": η["size"],
            "match_tt": m_tt,
            "match_rpre": m_rpre,
            "delta": m_tt - m_rpre,
        })
    return pd.DataFrame(rows)


def _trees_compatible(Z_rpre, Z_tt, Z_rpost) -> bool:
    if Z_rpre is None or Z_tt is None or Z_rpost is None:
        return False
    n = Z_rpre.shape[0] + 1
    return (Z_tt.shape[0] + 1 == n) and (Z_rpost.shape[0] + 1 == n)


def collect_raw() -> pd.DataFrame:
    """One row per (patient, band, rest_post internal node)."""
    parts = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            Z_rpre = _load_Z(pat, "rest_pre", band)
            Z_tt = _load_Z(pat, "task_test", band)
            Z_rpost = _load_Z(pat, "rest_post", band)
            if not _trees_compatible(Z_rpre, Z_tt, Z_rpost):
                continue
            df = score_rpost_nodes(Z_rpre, Z_tt, Z_rpost)
            df.insert(0, "patient", pat)
            df.insert(1, "band", band)
            parts.append(df)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


# ──────────────────────────── aggregation ────────────────────────────

def patient_scalar(raw: pd.DataFrame) -> pd.DataFrame:
    """S(pat, band) = Σ size·δ / Σ size (size-weighted mean δ over nodes)."""
    g = raw.groupby(["patient", "band"], sort=False)
    num = g.apply(lambda d: float((d["size"] * d["delta"]).sum()))
    den = g["size"].sum()
    out = (num / den).reset_index(name="S")
    out["n_nodes"] = g.size().values
    return out


def h_resolved(raw: pd.DataFrame, h_edges: np.ndarray
                ) -> tuple[dict[str, np.ndarray], list[str], np.ndarray, np.ndarray]:
    """Matrix dict {band: (n_patients, n_hbins)} of size-weighted δ per h-bin.

    ``h_edges`` is the length-(n+1) array of bin boundaries. Bin i holds
    nodes with ``h_rel ∈ [h_edges[i], h_edges[i+1])``; the last bin
    additionally absorbs ``h_rel == 1.0``.

    Returns ``(mats, patients, centers, edges)`` where ``centers`` is the
    geometric midpoint per bin (for log-axis plotting) and ``edges`` is
    ``h_edges`` passed through for convenience.
    """
    patients = sorted(raw["patient"].unique())
    n_p = len(patients)
    n_b = len(h_edges) - 1
    centers = np.sqrt(h_edges[:-1] * h_edges[1:])  # geometric midpoint
    out = {}
    for band in BRAIN_BANDS_NAMES:
        M = np.full((n_p, n_b), np.nan)
        bdf = raw[raw["band"] == band]
        for ip, pat in enumerate(patients):
            pdf = bdf[bdf["patient"] == pat]
            if pdf.empty:
                continue
            h = pdf["h_rel_rpost"].to_numpy()
            w = pdf["size"].to_numpy(dtype=float)
            d = pdf["delta"].to_numpy()
            idx = np.clip(np.searchsorted(h_edges, h, side="right") - 1,
                          0, n_b - 1)
            num = np.bincount(idx, weights=w * d, minlength=n_b)
            den = np.bincount(idx, weights=w, minlength=n_b)
            with np.errstate(invalid="ignore", divide="ignore"):
                M[ip, :] = np.where(den > 0, num / den, np.nan)
        out[band] = M
    return out, patients, centers, h_edges


# ──────────────────────────── stats ────────────────────────────

def band_tests(scalar: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(0)
    for band in BRAIN_BANDS_NAMES:
        x = scalar[scalar["band"] == band]["S"].to_numpy(float)
        z, p = wilcoxon_z(x)
        rb = rank_biserial(x)
        mean, lo, hi = boot_ci_mean(x, rng=rng)
        rows.append({"band": band, "n": int(np.isfinite(x).sum()),
                     "mean_S": mean, "ci_lo": lo, "ci_hi": hi,
                     "z": z, "p": p, "rb": rb})
    df = pd.DataFrame(rows)
    df["q_fdr"] = bh_fdr(df["p"].fillna(1.0).tolist())
    return df


def h_cluster_perm(h_mats: dict[str, np.ndarray], centers: np.ndarray
                    ) -> pd.DataFrame:
    rows = []
    for band in BRAIN_BANDS_NAMES:
        mat = h_mats[band]
        keep_row = ~np.isnan(mat).all(axis=1)
        m = mat[keep_row]
        ok = ~np.isnan(m).any(axis=0)
        if not ok.any():
            continue
        clusters = cluster_perm(m[:, ok], centers[ok])
        for c in clusters:
            rows.append({"band": band,
                         "h_start": c["scale_start"],
                         "h_end": c["scale_end"],
                         "len": c["len"],
                         "mass": c["mass"], "p": c["p"]})
    return pd.DataFrame(rows)


# ──────────────────────────── report ────────────────────────────

BAND_COLORS = {
    "delta":     "#3b6e9c", "theta":      "#c44e4e",
    "alpha":     "#5ea85e", "beta":       "#b28ad1",
    "low_gamma": "#e5a24b", "high_gamma": "#7d5a50",
}


def write_markdown(scalar_tab: pd.DataFrame, cluster_tab: pd.DataFrame,
                   out_path) -> None:
    lines: list[str] = []
    ap = lines.append
    ap("# H2-CBR — cluster-birth-retention (primary partition-level directed test)")
    ap("")
    ap("For every internal node η of the rest_post dendrogram, δ(η) is the "
       "Jaccard overlap with the best-matching task_test subtree minus the "
       "Jaccard overlap with the best-matching rest_pre subtree. The "
       "patient-band scalar S is the size-weighted mean of δ over all "
       "rest_post internal nodes.")
    ap("")
    ap("Primary: one-sample Wilcoxon S > 0 per band, FDR-BH m=6. "
       "Secondary: sign-flip cluster-permutation on S(h) across a log-h grid "
       "(60 bins, |z|>1.96, 5,000 perms).")
    ap("")
    ap("## Per-band scalar test")
    ap("")
    ap("| band | n | mean S | 95% CI | z | p (raw) | q (BH) | rb |")
    ap("|------|--:|-------:|:-------|--:|--------:|-------:|---:|")
    for _, r in scalar_tab.iterrows():
        star = "★" if r["q_fdr"] < 0.05 else ""
        ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | {r['n']} | "
           f"{r['mean_S']:+.4f} | [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] | "
           f"{r['z']:+.2f} | {r['p']:.4f} | {r['q_fdr']:.4f}{star} | "
           f"{r['rb']:+.3f} |")
    ap("")
    ap("## Scale-resolved cluster-permutation (log h_rel)")
    ap("")
    ap("| band | h-range | length | mass | p |")
    ap("|------|:--------|------:|-----:|--:|")
    sig = cluster_tab[cluster_tab["p"] < 0.05] if not cluster_tab.empty else cluster_tab
    if sig.empty:
        ap("| — | no p<0.05 clusters | — | — | — |")
    else:
        for _, r in sig.iterrows():
            ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | "
               f"[{r['h_start']:.3f}, {r['h_end']:.3f}] | "
               f"{int(r['len'])} | {r['mass']:.1f} | {r['p']:.4f}★ |")
    ap("")
    ap("### All supra-threshold clusters (including non-sig.)")
    ap("")
    ap("| band | h-range | length | mass | p |")
    ap("|------|:--------|------:|-----:|--:|")
    if cluster_tab.empty:
        ap("| — | no supra-threshold clusters | — | — | — |")
    else:
        for _, r in cluster_tab.iterrows():
            star = "★" if r["p"] < 0.05 else ""
            ap(f"| {BRAIN_BAND_TEX_DICT[r['band']]} | "
               f"[{r['h_start']:.3f}, {r['h_end']:.3f}] | "
               f"{int(r['len'])} | {r['mass']:.1f} | {r['p']:.4f}{star} |")
    ap("")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def draw_figure(h_mats: dict[str, np.ndarray], centers: np.ndarray,
                edges: np.ndarray, cluster_tab: pd.DataFrame,
                out_base) -> None:
    # Unanimity map: fraction of patients with positive S per (band, h).
    frac = np.full((len(BRAIN_BANDS_NAMES), len(centers)), np.nan)
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        M = h_mats[band]
        finite = np.isfinite(M)
        pos = finite & (M > 0)
        n = finite.sum(axis=0)
        p = pos.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            frac[ib, :] = np.where(n > 0, p / n, np.nan)

    # Auto-crop x-axis to columns that any band populates.
    populated = np.isfinite(frac).any(axis=0)
    if populated.any():
        x_lo = edges[np.argmax(populated)]
        last_col = len(populated) - 1 - np.argmax(populated[::-1])
        x_hi = edges[last_col + 1]
    else:
        x_lo, x_hi = edges[0], edges[-1]

    fig, axes = plt.subplots(
        2, 1, figsize=(12.5, 7.2), dpi=160,
        gridspec_kw={"height_ratios": [1.2, 0.9], "hspace": 0.45},
    )

    # Top: log-h unanimity map
    ax = axes[0]
    cmap = LinearSegmentedColormap.from_list(
        "trace",
        [(1.0, 1.0, 1.0), (1.0, 0.95, 0.7),
         (1.0, 0.6, 0.2), (0.7, 0.05, 0.05)], N=256,
    )
    im = ax.pcolormesh(edges, np.arange(len(BRAIN_BANDS_NAMES) + 1) - 0.5,
                        frac, cmap=cmap, vmin=0.0, vmax=1.0, shading="flat")
    ys, xs = np.where(frac >= 7 / 9)
    ax.scatter(centers[xs], ys, s=10, c="black", alpha=0.55,
                edgecolors="white", linewidths=0.3)
    ys, xs = np.where(frac >= 0.999)
    ax.scatter(centers[xs], ys, s=35, c="black",
                edgecolors="white", linewidths=0.5)
    ax.set_xscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=12)
    ax.invert_yaxis()
    ax.set_xlabel(r"$h_{\mathrm{rel}} = h / D_{\max}(\mathrm{rpost})$",
                   fontsize=11)
    ax.set_title(
        "(a) Cluster-birth-retention landscape — size-weighted δ(η) per "
        "rest_post internal node, binned on log-h_rel. "
        "Colour = fraction of patients with positive S(h); "
        "dots: ≥7/9 (faint), 9/9 (black).",
        fontsize=10, loc="left",
    )
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, fraction=0.03)
    cbar.set_label("fraction of patients with S(h) > 0", fontsize=9)

    # Bottom: cluster-perm strip
    ax = axes[1]
    ax.set_xscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-0.5, len(BRAIN_BANDS_NAMES) - 0.5)
    ax.invert_yaxis()
    ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
    ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                        fontsize=12)
    ax.set_xlabel(r"$h_{\mathrm{rel}} = h / D_{\max}(\mathrm{rpost})$",
                   fontsize=11)
    ax.grid(axis="x", which="both", alpha=0.3)
    # Geometric half-step (so a single-bin cluster still shows up).
    step_ratio = (edges[1] / edges[0]) ** 0.5
    for _, r in cluster_tab.iterrows():
        band = r["band"]
        ib = BRAIN_BANDS_NAMES.index(band)
        sig = r["p"] < 0.05
        x0 = r["h_start"] / step_ratio
        x1 = r["h_end"] * step_ratio
        rect = Rectangle(
            (x0, ib - 0.38), x1 - x0, 0.76,
            facecolor=BAND_COLORS[band],
            alpha=0.85 if sig else 0.32,
            edgecolor="black" if sig else BAND_COLORS[band],
            linewidth=1.2 if sig else 0.4,
        )
        ax.add_patch(rect)
        if sig:
            ax.text(np.sqrt(x0 * x1), ib,
                     f"p={r['p']:.3f}★", ha="center", va="center",
                     fontsize=8.5, fontweight="bold")
    ax.set_title(
        "(b) Sign-flip cluster-permutation on S(pat, band, h) "
        "(5,000 perms, |z|>1.96). "
        "Filled: supra-threshold; outlined+text: p<0.05.",
        fontsize=10, loc="left",
    )

    for ext in ("pdf", "png"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight")
        print(f"saved {out_base}.{ext}")
    plt.close(fig)


# ──────────────────────────── main ────────────────────────────

def main() -> None:
    raw = collect_raw()
    if raw.empty:
        raise SystemExit("no LRG caches found — aborting")
    out_dir = REPORTS_ROOT / "imcoh_vi"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    raw.to_csv(out_dir / "h2_cluster_birth_retention_raw.csv", index=False)

    scalar = patient_scalar(raw)
    scalar.to_csv(out_dir / "h2_cluster_birth_retention_patient.csv",
                  index=False)

    h_edges = h_log_grid([1.0], n=N_H + 1, h_min_floor=H_MIN_FLOOR)
    h_mats, patients, centers, h_edges = h_resolved(raw, h_edges)

    scalar_tab = band_tests(scalar)
    scalar_tab.to_csv(out_dir / "h2_cluster_birth_retention.csv", index=False)

    cluster_tab = h_cluster_perm(h_mats, centers)
    cluster_tab.to_csv(out_dir / "h2_cluster_birth_retention_clusters.csv",
                       index=False)

    write_markdown(scalar_tab, cluster_tab,
                    out_dir / "h2_cluster_birth_retention.md")
    print(f"wrote {out_dir / 'h2_cluster_birth_retention.md'}")

    draw_figure(h_mats, centers, h_edges, cluster_tab,
                 out_dir / "figures" / "h2_cluster_birth_retention")

    print("\nPer-band scalar summary:")
    for _, r in scalar_tab.iterrows():
        star = "★" if r["q_fdr"] < 0.05 else ""
        print(f"  {r['band']:12s}  mean_S={r['mean_S']:+.4f}  "
              f"z={r['z']:+.2f}  p={r['p']:.4f}  q={r['q_fdr']:.4f}{star}")


if __name__ == "__main__":
    main()
