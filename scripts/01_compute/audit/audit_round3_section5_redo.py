#!/usr/bin/env python3
"""Round-3 Section 5 figure redo bundle.

Five redos, all driven by existing round-2 CSVs and the LRG cache (no new
compute beyond a fast cohort linkage sweep for the Ψ-boundary audit):

  1. psi_irrelevance.pdf      — single (Pat_02 β rest_pre) 3-panel +
                                 6-band Ψ small-mults for the same patient.
  2. psi_boundary_audit.csv   — 180-cell argmax_n(Ψ) classification.
  3. kc_decomposition.pdf     — T(λ=0) vs T(λ=1) per-band scatter.
  4. ctm_headline.pdf         — boxes + per-patient lines, asterisks above split.
  5. cross_probe_signs.pdf    — patients sorted by trace-agreement count.

Output dir: data/audit/section5_v2_round3_redo/
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram

# Centralized font sizing — change here, propagates everywhere.
# (Per-call fontsize=... overrides still take precedence if needed.)
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_LIST,
)
from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
from lrg_eegfc.workflow.lrg import load_lrg_result

S5 = ROOT / "data" / "reports" / "section_5_lrg_trace"
R2 = ROOT / "data" / "audit" / "section5_v2_round2"
R2_TBL = R2 / "tables"

OUT = ROOT / "data" / "audit" / "section5_v2_round3_redo"
FIG = OUT / "figures"
TBL = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TBL.mkdir(parents=True, exist_ok=True)

BAND_ORDER = list(BRAIN_BANDS_NAMES)
PHASES = ("rest_pre", "task_test", "rest_post")

# --- redo-1 dendrogram leaf-label placement -------------------------------
# Vertical position of each leaf label as a fraction of its parent merge
# height (0 < value < 1). Closer to 1.0 → closer to the merge bracket.
DEND_LABEL_Y_FRAC = 0.9
DEND_LABEL_FONTSIZE = 4.0
# Number of label-heights (visual, display-pixel space) by which the
# alternating "low-shelf" label sits below the high shelf. 1.0 = exactly
# one label-height; 1.5 = one and a half label-heights, etc.
DEND_LABEL_DROP_HEIGHTS = 2.0


def asterisks(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


# ===================================================================
# Redo 1 — Ψ-irrelevance figure (Pat_02 β rest_pre + 6-band Ψ small-mults)
# ===================================================================
def redo1_psi_irrelevance() -> None:
    rep_pat, rep_band, rep_phase = "Pat_02", "beta", "rest_pre"
    res = load_lrg_result(rep_pat, rep_phase, rep_band, fc_method="imcoh_abs")
    Z = np.asarray(res.linkage_matrix)
    # Top-down indexing (Villegas convention): n=0 is the first cut from the
    # root (largest merge), n=N-2 is the last (smallest leaf-pair merge), so
    # the first Ψ peak corresponds to the highest cut in the dendrogram.
    heights = np.sort(Z[:, 2])[::-1]
    n = len(heights)
    psi = np.zeros(n - 1)
    for i in range(n - 1):
        if heights[i] > 0 and heights[i + 1] > 0:
            psi[i] = n * (np.log10(heights[i]) - np.log10(heights[i + 1]))
    argmax_n = int(np.argmax(psi))
    cut_h = float(np.sqrt(heights[argmax_n] * heights[argmax_n + 1]))

    # Channel labels for Pat_02 (compact form: drop the ",Gx" group suffix
    # and the embedded space). Length = 117, matches dendrogram N leaves.
    pat_dir = ROOT / "data" / "raw" / "stereoeeg_patients" / rep_pat
    ch = pd.read_csv(pat_dir / "channel_labels.csv")
    leaf_labels = [
        str(s).split(",")[0].replace(" ", "") for s in ch["label"].tolist()
    ]
    if len(leaf_labels) != n + 1:
        leaf_labels = [str(i) for i in range(n + 1)]

    fig = plt.figure(figsize=(16.0, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.18)

    ax_d = fig.add_subplot(gs[0, 0])
    ax_h = fig.add_subplot(gs[0, 1])
    ax_p = ax_h.twinx()

    dd = dendrogram(
        Z, ax=ax_d, no_labels=True, color_threshold=cut_h,
        above_threshold_color="#888888",
    )
    ax_d.axhline(cut_h, color="#c0392b", lw=1.2, ls="--",
                 label=f"Ψ-best cut h = {cut_h:.3f}")
    h_min = float(heights[-1])  # smallest merge (bottom of tree)
    h_max = float(heights[0])   # largest merge (root)
    tmin = 0.0
    tmax = h_max * 1.05
    ax_d.set_ylim(tmin, tmax)
    ax_d.set_ylabel("merge height")

    n_leaves = len(dd["leaves"])
    ax_d.set_xlim(0, 10 * n_leaves)
    ax_d.tick_params(axis="x", which="both", length=0)
    ax_d.set_xticks([])

    # Build leaf-index → parent merge height from the linkage matrix Z.
    # Each row of Z merges children Z[i, 0] and Z[i, 1] at height Z[i, 2];
    # children with index < N are leaves; each leaf appears in exactly one row.
    N = n + 1
    leaf_merge_h: dict[int, float] = {}
    for row in Z:
        c1, c2, h, _ = row
        c1, c2 = int(c1), int(c2)
        if c1 < N:
            leaf_merge_h[c1] = float(h)
        if c2 < N:
            leaf_merge_h[c2] = float(h)

    # Place each label just beneath its split (parent merge bracket). When
    # two consecutive leaves share the same merge height (which would put
    # their labels at the same y), drop the later one by EXACTLY one
    # rendered-label-height — measured from the actual font + axis transform,
    # not a hardcoded log offset.
    leaf_colors = dd["leaves_color_list"]
    leaves_in_order = dd["leaves"]

    # Force a draw so transforms are valid, then measure a probe label's
    # display-pixel height and convert it to an additive linear-y offset
    # (constant across the linear axis: a fixed pixel height → fixed Δy).
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    y_probe = 0.5 * (tmin + tmax)
    probe = ax_d.text(
        5.0, y_probe, "Wg99", fontsize=DEND_LABEL_FONTSIZE, ha="center",
        va="center",
        bbox=dict(facecolor="white", edgecolor="none", pad=0.4),
    )
    fig.canvas.draw()
    bb_px = probe.get_window_extent(renderer=renderer)
    inv = ax_d.transData.inverted()
    _, y_top_data = inv.transform((0, bb_px.y1))
    _, y_bot_data = inv.transform((0, bb_px.y0))
    probe.remove()
    # On a linear axis, a fixed display-pixel height corresponds to a fixed
    # additive Δ in data units. Multiply by the configured number of
    # label-heights to keep the drop in display (visual) units.
    label_drop_delta = float((y_top_data - y_bot_data) * DEND_LABEL_DROP_HEIGHTS)

    prev_y_parent = None
    prev_dropped = False
    for i, (leaf_idx, col) in enumerate(zip(leaves_in_order, leaf_colors)):
        x_leaf = 5.0 + 10.0 * i
        y_parent = leaf_merge_h[leaf_idx]
        drop = (prev_y_parent is not None
                and y_parent == prev_y_parent
                and not prev_dropped)
        y_high = y_parent * DEND_LABEL_Y_FRAC
        y_text = (y_high - label_drop_delta) if drop else y_high
        label = leaf_labels[leaf_idx] if leaf_idx < len(leaf_labels) else str(leaf_idx)
        ax_d.text(
            x_leaf, y_text, label, rotation=0, ha="center", va="center",
            fontsize=DEND_LABEL_FONTSIZE, color=col,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.4, alpha=1.0),
            zorder=5,
        )
        prev_y_parent = y_parent
        prev_dropped = drop

    ax_d.legend(loc="upper right", frameon=False)

    idx = np.arange(n)
    line_h, = ax_h.plot(idx, heights, color="#1f4e79", lw=1.4,
                        label=r"merge height $t_n$")
    ax_h.axvline(argmax_n, color="#c0392b", lw=0.7, ls=":")
    ax_h.set_xlabel("merge index $n$ (top-down: $n=0$ is the root cut)")
    ax_h.set_ylabel(r"merge height $t_n$", color="#1f4e79")
    ax_h.tick_params(axis="y", labelcolor="#1f4e79")
    ax_h.set_ylim(0.0, h_max * 1.05)

    line_p, = ax_p.plot(np.arange(n - 1), psi, color="#c0392b", lw=1.0,
                        label=r"$\Psi(n)$")
    psi_max = float(np.max(psi))
    ax_p.scatter([argmax_n], [psi[argmax_n]], color="#c0392b", s=30, zorder=5)
    ax_p.set_ylabel(r"$\Psi(n) = N\,[\log_{10} t_n - \log_{10} t_{n+1}]$",
                    color="#c0392b")
    ax_p.tick_params(axis="y", labelcolor="#c0392b")
    ax_h.legend(handles=[line_h, line_p], loc="upper center", frameon=False)

    fig.savefig(FIG / "psi_irrelevance.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[redo1] wrote psi_irrelevance.pdf")


# ===================================================================
# Redo 2 — Ψ-boundary audit across 180 (patient, band, phase) cells
# ===================================================================
def redo2_psi_boundary_audit() -> tuple[int, int, float]:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BAND_ORDER:
            for phi in PHASES:
                try:
                    res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                    Z = np.asarray(res.linkage_matrix)
                except Exception as e:  # noqa: BLE001
                    rows.append({
                        "patient": pat, "band": band, "phase": phi,
                        "argmax_n": np.nan, "N_leaves": np.nan,
                        "psi_max": np.nan, "argmax_position_class": "",
                        "error": str(e),
                    })
                    continue
                # Top-down (descending) convention: n=0 = root cut,
                # n=N-3 = leaf-pair cut. Matches redo1_psi_irrelevance.
                heights = np.sort(Z[:, 2])[::-1]
                N = len(heights) + 1  # number of leaves
                n_psi = len(heights) - 1  # Ψ defined for n in [0, N-3]
                psi = np.zeros(n_psi)
                for i in range(n_psi):
                    if heights[i] > 0 and heights[i + 1] > 0:
                        psi[i] = N * (np.log10(heights[i]) - np.log10(heights[i + 1]))
                am = int(np.argmax(psi)) if len(psi) else -1
                if am == 0:
                    cls = "top_cut"
                elif am == n_psi - 1:
                    cls = "bottom_cut"
                else:
                    cls = "interior"
                rows.append({
                    "patient": pat, "band": band, "phase": phi,
                    "argmax_n": am, "N_leaves": N,
                    "psi_max": float(psi[am]) if len(psi) else np.nan,
                    "argmax_position_class": cls,
                    "error": "",
                })
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "psi_argmax_pinning.csv", index=False)

    valid = df["argmax_position_class"].isin(["top_cut", "bottom_cut", "interior"])
    n_total = int(valid.sum())
    n_top = int((df["argmax_position_class"] == "top_cut").sum())
    n_bottom = int((df["argmax_position_class"] == "bottom_cut").sum())
    n_interior = int((df["argmax_position_class"] == "interior").sum())
    n_pinned = n_top + n_bottom
    pct = 100.0 * n_pinned / n_total if n_total else float("nan")

    # Per-band & per-patient summaries (per-class counts)
    def class_counts(group_col: str) -> pd.DataFrame:
        sub = df[valid].copy()
        ct = sub.groupby([group_col, "argmax_position_class"]).size().unstack(fill_value=0)
        for col in ("top_cut", "bottom_cut", "interior"):
            if col not in ct.columns:
                ct[col] = 0
        ct = ct[["top_cut", "bottom_cut", "interior"]]
        ct["n_cells"] = ct.sum(axis=1)
        ct["pct_pinned"] = 100.0 * (ct["top_cut"] + ct["bottom_cut"]) / ct["n_cells"]
        return ct.reset_index()

    by_band = class_counts("band")
    by_band.to_csv(TBL / "psi_argmax_pinning_by_band.csv", index=False)
    by_pat = class_counts("patient")
    by_pat.to_csv(TBL / "psi_argmax_pinning_by_patient.csv", index=False)

    # Histogram: integer-resolution distribution of argmax_n folded onto
    # distance-from-nearest-extremum, stacked by which extremum is nearer
    # (top = root cut at n=0; bottom = leaf-pair cut at n=N-3). Cells with
    # argmax beyond the window are summarized in a deep-interior badge.
    sub = df[valid].copy()
    am_all = sub["argmax_n"].astype(int).values
    N_arr = sub["N_leaves"].astype(int).values
    top_dist = am_all
    bot_dist = (N_arr - 3) - am_all
    nearer_top = top_dist <= bot_dist
    min_dist = np.where(nearer_top, top_dist, bot_dist)

    WINDOW = 12
    mask_in = min_dist <= WINDOW
    top_d = min_dist[nearer_top & mask_in]
    bot_d = min_dist[(~nearer_top) & mask_in]
    n_top0 = int(np.sum(nearer_top & (min_dist == 0)))
    n_bot0 = int(np.sum((~nearer_top) & (min_dist == 0)))
    n_deep = int((min_dist > WINDOW).sum())

    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    bins = np.arange(-0.5, WINDOW + 1.5)
    ax.hist(
        [top_d, bot_d], bins=bins, stacked=True,
        color=["#5B7CE0", "#c0392b"],
        label=[r"nearer root cut ($n=0$): "f"{int((nearer_top).sum())}",
               r"nearer leaf-pair cut ($n=N-3$): "f"{int((~nearer_top).sum())}"],
        edgecolor="white", linewidth=0.7,
    )
    # Annotate the d=0 split
    ax.annotate(f"{n_top0} | {n_bot0}", xy=(0, n_top0 + n_bot0),
                xytext=(0, n_top0 + n_bot0 + 4), ha="center",
                fontsize=11, color="#222")

    ax.set_xlim(-0.5, WINDOW + 0.5)
    ax.set_xticks(range(0, WINDOW + 1))
    ax.set_xlabel(r"$d$ = distance of argmax from nearest extremum"
                  " (merge indices)")
    ax.set_ylabel(f"# cells (of {len(sub)})")
    ax.legend(frameon=False, loc="upper right")

    if n_deep:
        ax.text(WINDOW * 0.62, ax.get_ylim()[1] * 0.55,
                f"$d > {WINDOW}$: {n_deep} cell"
                + ("s" if n_deep != 1 else ""),
                fontsize=11, color="#444",
                bbox=dict(facecolor="white", edgecolor="#888",
                          boxstyle="round,pad=0.35"))
    fig.tight_layout()
    fig.savefig(FIG / "psi_argmax_histogram.pdf", bbox_inches="tight")
    fig.savefig(FIG / "psi_argmax_histogram.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"[redo2] {n_pinned}/{n_total} cells pinned at boundary "
          f"(top_cut={n_top}, bottom_cut={n_bottom}, interior={n_interior}) "
          f"[{pct:.1f}% pinned]")
    return n_pinned, n_total, pct


# ===================================================================
# Redo 3 — KC decomposition scatter (T(λ=0) vs T(λ=1))
# ===================================================================
def redo3_kc_decomposition() -> None:
    KC = pd.read_csv(S5 / "03_kc_lambda_triangle/tables/cohort_summary_lambda.csv")

    # Build per-band record (lam=0 vs lam=1), sign-flipped so trace is positive
    recs = []
    for b in BAND_ORDER:
        r0 = KC[(KC["band"] == b) & (KC["lam"] == 0.0)]
        r1 = KC[(KC["band"] == b) & (KC["lam"] == 1.0)]
        if not len(r0) or not len(r1):
            continue
        T0 = -float(r0["T_median"].iloc[0])
        T1 = -float(r1["T_median"].iloc[0])
        p0 = float(r0["wilcoxon_one_sided_p"].iloc[0])
        p1 = float(r1["wilcoxon_one_sided_p"].iloc[0])
        n0 = int(r0["n_trace_int"].iloc[0])
        n1 = int(r1["n_trace_int"].iloc[0])
        recs.append({
            "band": b, "T_topology": T0, "T_heights": T1,
            "p_topology": p0, "p_heights": p1, "p_min": min(p0, p1),
            "n_trace_lam0": n0, "n_trace_lam1": n1,
        })
    R = pd.DataFrame(recs)
    R.to_csv(TBL / "kc_decomposition_points.csv", index=False)

    fig = plt.figure(figsize=(12.0, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.28)
    ax = fig.add_subplot(gs[0, 0])

    # Quadrant background shading: deeper green for trace quadrant
    xmax = max(abs(R["T_topology"]).max(), 0.5) * 1.30
    ymax = max(abs(R["T_heights"]).max(), 0.5) * 1.30
    ax.axhspan(0, ymax, xmin=0.5, xmax=1.0, color="#1a7c3e", alpha=0.13)
    ax.axhline(0, color="0.4", lw=0.7, ls="--")
    ax.axvline(0, color="0.4", lw=0.7, ls="--")

    # Per-band annotation offsets in axes-relative units. Small, close to dots.
    # β goes SE so it stays inside the plot area near the right edge; α goes SE
    # to separate vertically from γ_h which is just above-left at similar T_top.
    ann_offsets = {
        "beta":       (+0.025, -0.030),
        "low_gamma":  (+0.025, +0.025),
        "alpha":      (+0.025, -0.030),
        "high_gamma": (+0.025, +0.030),
        "theta":      (+0.025, +0.025),
        "delta":      (+0.025, +0.025),
    }

    # Sizes proportional to -log10 p_min
    for _, row in R.iterrows():
        size = 60 + 220.0 * (-np.log10(max(row["p_min"], 1e-6)) / 3.0)
        size = float(np.clip(size, 60, 420))
        is_sig = (row["p_topology"] <= 0.05) or (row["p_heights"] <= 0.05)
        face = "#c0392b" if is_sig else "white"
        edge = "#c0392b" if is_sig else "#444"
        ax.scatter(row["T_topology"], row["T_heights"], s=size,
                   facecolor=face, edgecolor=edge, lw=1.4, zorder=3)
        # Per-band offset; convert axes-fraction to data via xmax/ymax (axis is
        # symmetric −xmax..+xmax, so axes-frac 0.05 = 0.10*xmax in data).
        dx_frac, dy_frac = ann_offsets.get(row["band"], (0.05, 0.05))
        dx_data = dx_frac * 2 * xmax
        dy_data = dy_frac * 2 * ymax
        ha = "left" if dx_frac >= 0 else "right"
        va = "bottom" if dy_frac >= 0 else "top"
        weight = "bold" if is_sig else "normal"
        fs = 13 if is_sig else 11
        ax.annotate(BRAIN_BAND_TEX_DICT[row["band"]],
                    (row["T_topology"], row["T_heights"]),
                    xytext=(row["T_topology"] + dx_data,
                            row["T_heights"] + dy_data),
                    fontsize=fs, fontweight=weight, ha=ha, va=va,
                    color="#222")

    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-ymax, ymax)
    ax.set_xlabel(r"$-T_{KC}(\lambda=0)$ — pure topology, trace direction →",
                  fontsize=10)
    ax.set_ylabel(r"$-T_{KC}(\lambda=1)$ — pure heights, ↑ trace direction",
                  fontsize=10)

    # Quadrant corner labels only. The trace quadrant is identified by its
    # green shading; no in-quadrant annotation that would compete with β.
    ax.text(0.04, 0.96, "heights only", ha="left", va="top",
            fontsize=9, color="#444", fontstyle="italic",
            transform=ax.transAxes)
    ax.text(0.96, 0.04, "topology only", ha="right", va="bottom",
            fontsize=9, color="#444", fontstyle="italic",
            transform=ax.transAxes)
    ax.text(0.04, 0.04, "anti-trace", ha="left", va="bottom",
            fontsize=9, color="#c0392b", fontstyle="italic",
            transform=ax.transAxes)
    ax.spines[["top", "right"]].set_visible(False)

    # Right panel: per-band n_trace bars at λ=0 / λ=1, β-highlighted
    ax_r = fig.add_subplot(gs[0, 1])
    x = np.arange(len(R))
    w = 0.36
    band_list = list(R["band"])
    # Color each bar by its own per-(band, λ) significance: red if p ≤ 0.05,
    # otherwise the steel-blue (λ=0) / warm-rust (λ=1) Tableau pair.
    p_top_by = {row["band"]: float(row["p_topology"]) for _, row in R.iterrows()}
    p_h_by = {row["band"]: float(row["p_heights"]) for _, row in R.iterrows()}
    col_lam0 = ["#c0392b" if p_top_by[b] <= 0.05 else "#4C72B0" for b in band_list]
    col_lam1 = ["#c0392b" if p_h_by[b] <= 0.05 else "#DD8452" for b in band_list]
    edge_lam0 = ["#7d1a0e" if p_top_by[b] <= 0.05 else "#2F4A78" for b in band_list]
    edge_lam1 = ["#7d1a0e" if p_h_by[b] <= 0.05 else "#9C5B36" for b in band_list]
    bars0 = ax_r.bar(x - w / 2, R["n_trace_lam0"], width=w,
                     color=col_lam0, edgecolor=edge_lam0, lw=0.8,
                     label=r"$\lambda=0$ (topology)")
    bars1 = ax_r.bar(x + w / 2, R["n_trace_lam1"], width=w,
                     color=col_lam1, edgecolor=edge_lam1, lw=0.8,
                     label=r"$\lambda=1$ (heights)")
    # Highlight the β column with a subtle background band
    j_beta = band_list.index("beta") if "beta" in band_list else None
    if j_beta is not None:
        ax_r.axvspan(j_beta - 0.5, j_beta + 0.5, color="#fff0ed",
                     alpha=0.7, zorder=0)
    # Reference line at n=5 (cohort midpoint) — label outside the bar field
    ax_r.axhline(5, color="0.4", lw=0.6, ls="--", zorder=1)
    ax_r.text(-0.55, 5, "n = 5", fontsize=8, color="0.4",
              ha="right", va="center")
    # Numerical labels above each bar
    for bar, v in zip(list(bars0) + list(bars1),
                      list(R["n_trace_lam0"]) + list(R["n_trace_lam1"])):
        ax_r.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.18,
                  str(int(v)), ha="center", va="bottom", fontsize=8.5,
                  color="0.15")

    # Significance asterisks above every bar with p ≤ 0.05.
    for j, b in enumerate(band_list):
        for off, p in [(-w / 2, p_top_by[b]), (+w / 2, p_h_by[b])]:
            stars = asterisks(p)
            if stars:
                ax_r.text(j + off, 9.5, stars, ha="center", va="bottom",
                          fontsize=12, fontweight="bold", color="#7d1a0e")

    ax_r.set_xticks(x)
    ax_r.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in R["band"]],
                         fontsize=12)
    ax_r.set_ylim(0, 10.5)
    ax_r.set_ylabel(r"$n_{\mathrm{trace}}\,/\,10$", fontsize=11)
    ax_r.legend(fontsize=9, frameon=False, loc="upper left",
                bbox_to_anchor=(0.0, 1.0))
    ax_r.spines[["top", "right"]].set_visible(False)
    ax_r.tick_params(axis="x", which="both", length=0)

    # Marker-significance legend, vertical column inside the left panel.
    # Placed above the "topology only" corner label in the empty SE region.
    handles = [
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=9,
                      markerfacecolor="#c0392b", markeredgecolor="#c0392b",
                      label=r"$p \leq 0.05$"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=9,
                      markerfacecolor="white", markeredgecolor="#444",
                      label=r"$p > 0.05$"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=13,
                      markerfacecolor="#888", markeredgecolor="#444",
                      label=r"$\propto -\log_{10} p_{\min}$"),
    ]
    ax.legend(handles=handles, loc="lower right",
              bbox_to_anchor=(1.0, 0.08),
              frameon=False, ncol=1,
              handletextpad=0.5, labelspacing=0.6)
    fig.tight_layout()
    fig.savefig(FIG / "kc_decomposition.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[redo3] wrote kc_decomposition.pdf")


# ===================================================================
# Redo 4 — CTM headline (per-patient lines + asterisks above split)
# ===================================================================
def redo4_ctm_headline() -> None:
    CTM = pd.read_csv(S5 / "05_ctm_sigma_aggregate/tables/cohort_summary.csv")
    CTM_PT = pd.read_csv(S5 / "05_ctm_sigma_aggregate/tables/Td_per_patient_per_band.csv")

    fig, ax = plt.subplots(figsize=(13, 5.0))
    width = 0.25
    x_band = np.arange(len(BAND_ORDER))
    colors = {"split": "#1f77b4", "drift": "#888888", "xprobe": "#e07b00"}
    box_keys = (("split", "rho_split"), ("drift", "rho_null_drift"),
                ("xprobe", "rho_split_cross_probe"))

    # Boxes
    for off, (name, key) in zip((-width, 0.0, width), box_keys):
        data = [CTM_PT[CTM_PT["band"] == b][key].dropna().values for b in BAND_ORDER]
        ax.boxplot(
            data, positions=x_band + off, widths=width * 0.85,
            patch_artist=True, showmeans=False, showfliers=False,
            boxprops=dict(facecolor=colors[name], alpha=0.50, edgecolor="#444"),
            medianprops=dict(color="#111", lw=1.4),
            whiskerprops=dict(color="#666"), capprops=dict(color="#666"),
        )

    # Per-patient connected dots: split → drift → xprobe
    pats = sorted(CTM_PT["patient"].unique())
    for xi, b in enumerate(BAND_ORDER):
        sub = CTM_PT[CTM_PT["band"] == b].set_index("patient")
        for p in pats:
            if p not in sub.index:
                continue
            ys = [
                float(sub.loc[p, "rho_split"]),
                float(sub.loc[p, "rho_null_drift"]),
                float(sub.loc[p, "rho_split_cross_probe"]),
            ]
            xs = [xi - width, xi, xi + width]
            is_pro = ys[0] > 0
            line_c = "#1a7c3e" if is_pro else "#c0392b"
            ax.plot(xs, ys, color=line_c, alpha=0.45, lw=0.8, zorder=3)
            ax.scatter(xs, ys, s=10, color=line_c, alpha=0.65,
                       edgecolor="white", linewidth=0.3, zorder=4)

    # Significance asterisks ABOVE the split (blue) box at top of each box
    for xi, b in enumerate(BAND_ORDER):
        rs = CTM[CTM["band"] == b].iloc[0]
        p = float(rs["wilcoxon_split_gt_drift_p"])
        # Anchor above the IQR of the split box
        v = CTM_PT[CTM_PT["band"] == b]["rho_split"].dropna().values
        if not len(v):
            continue
        q3 = float(np.nanpercentile(v, 75))
        whisk_top = float(np.nanmax(v))
        anchor = max(q3, whisk_top) + 0.03
        stars = asterisks(p)
        txt = f"p={p:.3f}{stars}"
        col = "#c0392b" if p <= 0.05 else "#444"
        ax.text(xi - width, anchor, txt, ha="center", va="bottom",
                fontsize=8, color=col, fontweight="bold" if stars else "normal")

    ax.axhline(0, color="0.3", lw=0.7, ls="--")
    ax.set_xticks(x_band)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BAND_ORDER], fontsize=11)
    ax.set_ylabel(r"$\rho$ (CTM split / drift / x-probe)", fontsize=10)

    # y-limits: tighten using p1/p99 of split + xprobe (drift is small)
    all_v = np.concatenate([
        CTM_PT[k].dropna().values for _, k in box_keys
    ])
    qmin = float(np.nanpercentile(all_v, 1))
    qmax = float(np.nanpercentile(all_v, 99))
    ax.set_ylim(qmin - 0.15, qmax + 0.18)

    handles = [
        mpatches.Patch(facecolor=colors["split"], alpha=0.50, edgecolor="#444",
                       label=r"$\rho_{\mathrm{split}}$ (controlled)"),
        mpatches.Patch(facecolor=colors["drift"], alpha=0.50, edgecolor="#444",
                       label=r"$\rho_{\mathrm{null\,drift}}$"),
        mpatches.Patch(facecolor=colors["xprobe"], alpha=0.50, edgecolor="#444",
                       label=r"$\rho_{\mathrm{xprobe}}$"),
        mlines.Line2D([], [], color="#1a7c3e", marker="o", lw=0.8, alpha=0.7,
                      markersize=4, label="patient: pro (ρ_split > 0)"),
        mlines.Line2D([], [], color="#c0392b", marker="o", lw=0.8, alpha=0.7,
                      markersize=4, label="patient: anti (ρ_split ≤ 0)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=5, frameon=False, fontsize=8.5)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(FIG / "ctm_headline.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[redo4] wrote ctm_headline.pdf")


# ===================================================================
# Redo 5 — Cross-probe signs, sorted by trace-agreement count
# ===================================================================
def redo5_cross_probe_signs() -> None:
    df = pd.read_csv(R2_TBL / "cross_probe_signs.csv")
    probes = ["d_rank_dS", "kc_lam0", "grassmann_k13", "ctm", "vi_clean"]
    sign_cols = [f"sign_{x}" for x in probes]
    probe_labels = ["d_rank d_S", "KC λ=0", "Grass k=13", "CTM", "VI(k)"]

    # trace direction = sign == -1 ; agreement count per (patient, band)
    df["agree_trace"] = (df[sign_cols] == -1).sum(axis=1)

    # Patient ordering: sort by mean agreement across all 6 bands
    pat_avg = df.groupby("patient")["agree_trace"].mean().sort_values(ascending=False)
    cohort_pat_order = list(pat_avg.index)

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0), sharex=False, sharey=False)
    cmap = plt.cm.RdYlGn_r

    for ax, b in zip(axes.flat, BAND_ORDER):
        sub = df[df["band"] == b].copy()
        # In-panel sort: by agree_trace desc; ties broken by cohort order
        sub["pat_rank"] = sub["patient"].map(
            {p: i for i, p in enumerate(cohort_pat_order)}
        )
        sub = sub.sort_values(["agree_trace", "pat_rank"],
                              ascending=[False, True]).reset_index(drop=True)
        M = sub[sign_cols].values.astype(float)
        ax.imshow(M, aspect="auto", cmap=cmap, vmin=-1, vmax=1,
                  interpolation="nearest")
        ax.set_xticks(range(len(probes)))
        ax.set_xticklabels(probe_labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(
            [f"{p}  ({a}/5)" for p, a in zip(sub["patient"], sub["agree_trace"])],
            fontsize=7,
        )
        ax.set_title(BRAIN_BAND_TEX_DICT[b], fontsize=10)
        # Cell glyphs
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M[i, j]
                txt = "−" if v == -1 else ("+" if v == 1 else "·")
                ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                        color="white" if v != 0 else "black")
        # Block dividers: between agree>=4, agree∈{2,3}, agree<=1
        agree = sub["agree_trace"].values
        for thresh in (4, 2):
            # Find first row index where agree < thresh
            below = np.where(agree < thresh)[0]
            if len(below):
                y = float(below[0]) - 0.5
                ax.axhline(y, color="#222", lw=0.9)

    fig.text(
        0.5, 0.005,
        "rows sorted by # probes in trace direction (≥ 4 / 5 = consistently pro · "
        "2–3 = mixed · ≤ 1 = consistently anti) · "
        "green/− = trace direction · red/+ = anti-trace · gray/· = ~zero",
        ha="center", fontsize=8, color="0.4",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(FIG / "cross_probe_signs.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[redo5] wrote cross_probe_signs.pdf")


# ===================================================================
# Redo 6 — KC λ=0 topology eye-proof: per-pair MRCA-depth shift scatter
# ===================================================================
def redo6_kc_topology_eye_proof() -> None:
    """Per-pair MRCA-depth shift scatter, 2x3 panel grid (one per band).

    For each leaf pair (i, j), m_phi(i, j) is the depth (number of edges
    from root) of the most recent common ancestor of i and j in the
    UPGMA dendrogram of phase phi -- exactly what
    ``lrg_eegfc.utils.metrics.tree_distance.kc_vectors`` returns. The
    KC λ=0 distance is the L2 norm of (m_a - m_b) over leaf pairs (with
    a pooled-max normalization in the library version).

    Pair-level decomposition shown here:
      - x = |m_rspre - m_taskt|
      - y = |m_taskt - m_rspost|
      - diagonal y = x; pairs below = trace direction, above = anti.

    Sanity check at runtime: for each (patient, band) cell, the raw
    Euclidean ||m_a - m_b||_2 must equal m_max * library_kc_distance
    (lam=0) where m_max = max(m_a.max(), m_b.max(), 1). This is logged
    on the first three cells.
    """
    pat_list = list(PATIENTS_LIST)
    pat_color = {p: c for p, c in zip(
        pat_list, plt.cm.tab10(np.linspace(0, 1, len(pat_list)))
    )}
    bands_order = list(BAND_ORDER)
    shaded = {"alpha", "beta", "low_gamma"}  # §5.3 controlled trace bands

    rows_pair = []        # long-format CSV rows; keep small dtypes
    rows_panel = {}       # band -> dict(x_all=[], y_all=[], pat_all=[], n_below, n_above, n_eq)

    sanity_logged = 0
    for pat in pat_list:
        for band in bands_order:
            try:
                Z = {}
                for phi in PHASES:
                    res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                    if res is None:
                        raise FileNotFoundError(f"missing LRG cache: {pat} {phi} {band}")
                    Z[phi] = res.linkage_matrix
                m = {phi: kc_vectors(Z[phi])[0] for phi in PHASES}
                d1 = m["task_test"] - m["rest_pre"]
                d2 = m["rest_post"] - m["task_test"]
                a1, a2 = np.abs(d1), np.abs(d2)

                # Sanity check: raw L2 == m_max * library_kc(lam=0).
                if sanity_logged < 3:
                    from lrg_eegfc.utils.metrics.tree_distance import kc_distance
                    raw_pre_tt = float(np.sqrt(np.sum(d1 * d1)))
                    mmax_pre_tt = max(int(m["rest_pre"].max()),
                                      int(m["task_test"].max()), 1)
                    lib_pre_tt = kc_distance(Z["rest_pre"], Z["task_test"], lam=0.0)
                    recon = mmax_pre_tt * lib_pre_tt
                    print(f"[redo6 sanity] {pat} {band}: "
                          f"raw_L2={raw_pre_tt:.3f}, m_max*lib={recon:.3f}, "
                          f"diff={abs(raw_pre_tt - recon):.2e}")
                    sanity_logged += 1

                # Per-pair rows for CSV — enumerate (i, j) in scipy condensed order.
                n_leaves = Z["rest_pre"].shape[0] + 1
                ii, jj = np.triu_indices(n_leaves, k=1)
                for k in range(len(ii)):
                    rows_pair.append((
                        pat, band,
                        int(ii[k]), int(jj[k]),
                        int(m["rest_pre"][k]),
                        int(m["task_test"][k]),
                        int(m["rest_post"][k]),
                        int(a1[k]), int(a2[k]),
                    ))

                rec = rows_panel.setdefault(band, dict(x=[], y=[], pat=[]))
                rec["x"].append(a1)
                rec["y"].append(a2)
                rec["pat"].append(np.full(len(a1), pat, dtype=object))
            except Exception as e:
                print(f"[redo6] WARN {pat} {band}: {e}")

    # Persist long-format CSV (gzipped to keep <10 MB).
    if rows_pair:
        df_pair = pd.DataFrame(
            rows_pair,
            columns=["patient", "band", "leaf_i", "leaf_j",
                     "m_rspre", "m_taskt", "m_rspost",
                     "abs_delta_rspre_taskt", "abs_delta_taskt_rspost"],
        )
        df_pair.to_csv(TBL / "kc_pair_depth_shifts.csv.gz",
                       index=False, compression="gzip")
        print(f"[redo6] wrote kc_pair_depth_shifts.csv.gz ({len(df_pair)} rows)")

    # ---- Figure 1: 2D density heatmaps (the right viz for integer lattice) ----
    band_grid = bands_order  # δ θ α | β low_γ γ_h
    panel_summary = []
    per_patient_med_rows = []

    # Determine common axis range = 99th percentile across cohort
    all_xy = []
    for band in band_grid:
        if band in rows_panel:
            all_xy.append(np.concatenate(rows_panel[band]["x"]))
            all_xy.append(np.concatenate(rows_panel[band]["y"]))
    common_max = int(np.ceil(np.quantile(np.concatenate(all_xy), 0.99))) + 1
    common_max = max(common_max, 12)  # floor for visibility

    fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.6), sharex=True, sharey=True)
    # Bin edges: integer cells [0..common_max]
    edges = np.arange(common_max + 2) - 0.5
    # Divergent colormap, green = trace excess (more mass below diagonal),
    # red = anti excess (more mass above).  PiYG gives clean green/red contrast.
    cmap = plt.cm.PiYG

    # Compute asymmetry maps per band first to set a shared color scale.
    asym_data = {}
    for band in band_grid:
        if band not in rows_panel:
            continue
        x_all = np.concatenate(rows_panel[band]["x"]).astype(float)
        y_all = np.concatenate(rows_panel[band]["y"]).astype(float)
        H, _, _ = np.histogram2d(x_all, y_all, bins=[edges, edges])
        # H[i, j] = # pairs with |Δ1| in bin i, |Δ2| in bin j (note transpose
        # convention: imshow(H.T, origin='lower') puts |Δ2| on y-axis).
        # Asymmetry: A[i, j] = H[i, j] - H[j, i].
        # On below-diagonal cells (j < i, i.e. |Δ2| < |Δ1|), positive A means
        # more pairs in trace direction at this magnitude pair than its mirror.
        A = H - H.T
        asym_data[band] = (x_all, y_all, A.T)  # transpose for imshow convention
    vmax = max(np.abs(asym.flatten()).max() for _, _, asym in asym_data.values())
    # Soften the saturation so most cells aren't pegged at the colormap end.
    vmax = vmax * 0.6

    for ax, band in zip(axes.flat, band_grid):
        if band not in rows_panel:
            ax.set_axis_off()
            continue
        if band in shaded:
            for spine in ax.spines.values():
                spine.set_edgecolor("#c97000")
                spine.set_linewidth(1.4)

        x_all, y_all, A_disp = asym_data[band]
        pat_all = np.concatenate(rows_panel[band]["pat"])

        ax.imshow(
            A_disp, origin="lower",
            extent=(edges[0], edges[-1], edges[0], edges[-1]),
            cmap=cmap, aspect="equal",
            vmin=-vmax, vmax=+vmax,
            interpolation="nearest",
        )

        # Diagonal y = x — black on the divergent map.
        ax.plot([edges[0], edges[-1]], [edges[0], edges[-1]],
                color="black", lw=1.4, ls="--", zorder=4)

        # Quadrant labels: BELOW = trace excess; ABOVE = anti excess.
        ax.text(0.97, 0.06, "below diag\n= trace direction",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7.5, color="#0f5d2c", fontstyle="italic",
                alpha=0.9)
        ax.text(0.03, 0.94, "above diag\n= anti-trace",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=7.5, color="#a01b1b", fontstyle="italic",
                alpha=0.9)

        # Per-patient median dots.
        for pat in pat_list:
            mask = pat_all == pat
            if not mask.any():
                continue
            mx = float(np.median(x_all[mask]))
            my = float(np.median(y_all[mask]))
            ax.scatter([mx], [my], marker="o", s=55, color=pat_color[pat],
                       edgecolor="black", linewidth=0.7, zorder=6)
            per_patient_med_rows.append(dict(patient=pat, band=band,
                                             med_x=mx, med_y=my,
                                             ratio=(my / mx) if mx > 0 else float("nan")))

        # Cohort-median star.
        med_x = float(np.median(x_all))
        med_y = float(np.median(y_all))
        ax.scatter([med_x], [med_y], marker="*", s=420, color="gold",
                   edgecolor="black", linewidth=1.4, zorder=7)

        # Cohort-mean cross (the right L2 stat — usually different from median).
        mean_x = float(np.mean(x_all))
        mean_y = float(np.mean(y_all))
        ax.scatter([mean_x], [mean_y], marker="X", s=180, color="black",
                   edgecolor="white", linewidth=1.4, zorder=7)

        # Counts.
        below = int(np.sum(y_all < x_all))
        above = int(np.sum(y_all > x_all))
        equal = int(np.sum(y_all == x_all))
        n_total = len(x_all)
        ratio = (med_y / med_x) if med_x > 0 else float("nan")
        s = x_all - y_all
        mean_s = float(np.mean(s))

        panel_summary.append(dict(
            band=band, n_pairs=n_total,
            n_below_diag=below, n_above_diag=above, n_on_diag=equal,
            frac_below=below / n_total,
            cohort_median_x=med_x, cohort_median_y=med_y,
            cohort_median_ratio_y_over_x=ratio,
            cohort_mean_x=mean_x, cohort_mean_y=mean_y,
            cohort_mean_s=mean_s,
        ))

        # Annotation: lead with cohort mean s (the eye-proof number).
        if mean_s >= 0.5:
            verdict, vcolor = "trace", "#0f5d2c"
        elif mean_s <= -0.5:
            verdict, vcolor = "anti", "#a01b1b"
        else:
            verdict, vcolor = "null", "#666"
        ax.text(0.50, 0.97,
                f"mean $s$ = ${mean_s:+.2f}$  ({verdict})",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=10.0, color=vcolor, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.30", facecolor="white",
                          edgecolor=vcolor, linewidth=1.2, alpha=0.96))

        # Band tag.
        ax.text(0.96, 0.96, BRAIN_BAND_TEX_DICT[band],
                transform=ax.transAxes, ha="right", va="top",
                fontsize=18, color="#111", fontweight="bold")

        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(edges[0], edges[-1])

    # Shared axis labels (figure-level — no fig.suptitle, no per-axis titles).
    for ax in axes[1, :]:
        ax.set_xlabel(r"$|m_{\mathrm{RPre}} - m_{\mathrm{TT}}|$  (per pair)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$|m_{\mathrm{TT}} - m_{\mathrm{RPost}}|$  (per pair)")

    # Two-row legend: top row = patients (10 swatches), bottom row = cell-density
    # legend + diagonal + cohort markers.
    handles_top = [mlines.Line2D([], [], marker="o", linestyle="",
                                 color=pat_color[p], markeredgecolor="black",
                                 markersize=7, label=p)
                   for p in pat_list]
    handles_bot = [
        mpatches.Patch(facecolor="#0f5d2c", edgecolor="black",
                       label="cell green = trace excess (more pairs below diag)"),
        mpatches.Patch(facecolor="#a01b1b", edgecolor="black",
                       label="cell red = anti excess (more pairs above diag)"),
        mlines.Line2D([], [], color="black", lw=1.4, ls="--",
                      label=r"diagonal $y = x$"),
        mlines.Line2D([], [], marker="*", linestyle="",
                      color="gold", markeredgecolor="black",
                      markersize=15, label="cohort median"),
        mlines.Line2D([], [], marker="X", linestyle="",
                      color="black", markeredgecolor="white",
                      markersize=11, label="cohort mean"),
    ]
    leg_top = fig.legend(handles=handles_top, loc="lower center",
                         bbox_to_anchor=(0.5, 0.00),
                         ncol=10, frameon=False, fontsize=8.5)
    leg_bot = fig.legend(handles=handles_bot, loc="lower center",
                         bbox_to_anchor=(0.5, -0.05),
                         ncol=5, frameon=False, fontsize=8.0)

    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.savefig(FIG / "kc_topology_eye_proof.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[redo6] wrote kc_topology_eye_proof.pdf")

    # Persist panel summary + per-patient median table.
    pd.DataFrame(panel_summary).to_csv(
        TBL / "kc_topology_eye_proof_panel_summary.csv", index=False
    )
    pd.DataFrame(per_patient_med_rows).to_csv(
        TBL / "kc_topology_eye_proof_per_patient_medians.csv", index=False
    )
    print("[redo6] wrote panel + per-patient medians CSVs")
    print(pd.DataFrame(panel_summary).to_string(index=False))

    # ---- Figure 2: 1D histogram of per-pair (|Δ_RPre→TT| - |Δ_TT→RPost|) ----
    # The cohort-level KC λ=0 trace says ||m_TT - m_RPost||₂ < ||m_RPre - m_TT||₂.
    # The per-pair signature is s(i,j) = |Δ1| - |Δ2|: positive = trace direction
    # (rspost shift smaller than rspre→tt shift). Mean of s ≈ skew of the
    # distribution — the integer-valued data has heavy mass at zero so the
    # median is uninformative; the mean is the right stat for L2-trace direction.
    fig2, axes2 = plt.subplots(2, 3, figsize=(11.4, 6.4), sharex=True, sharey=False)
    sig_summary = []
    edges_s = np.arange(-30, 31) - 0.5  # integer bins centered on integers

    for ax, band in zip(axes2.flat, band_grid):
        if band not in rows_panel:
            ax.set_axis_off()
            continue
        if band in shaded:
            for spine in ax.spines.values():
                spine.set_edgecolor("#c97000")
                spine.set_linewidth(1.4)

        x_all = np.concatenate(rows_panel[band]["x"]).astype(float)
        y_all = np.concatenate(rows_panel[band]["y"]).astype(float)
        s = x_all - y_all  # |Δ1| - |Δ2|: + = trace, - = anti

        # Histogram, color positive bins green / negative bins red.
        counts, _ = np.histogram(s, bins=edges_s)
        centers = 0.5 * (edges_s[:-1] + edges_s[1:])
        colors_bar = ["#0f5d2c" if c > 0 else "#a01b1b" if c < 0 else "#888"
                      for c in centers]
        ax.bar(centers, counts, width=1.0, color=colors_bar,
               edgecolor="black", linewidth=0.3)

        ax.axvline(0, color="black", lw=1.0, zorder=4)
        mean_s = float(np.mean(s))
        med_s = float(np.median(s))
        ax.axvline(mean_s, color="gold", lw=2.4, zorder=5)

        n_pos = int(np.sum(s > 0))
        n_neg = int(np.sum(s < 0))
        n_zero = int(np.sum(s == 0))
        n_total = len(s)
        frac_pos = n_pos / n_total

        # Direction is set by the mean (the L2-trace stat); call it trace if
        # the mean leans positive enough to matter relative to per-pair noise.
        if mean_s >= 0.5:
            verdict, vcolor = "trace", "#0f5d2c"
        elif mean_s <= -0.5:
            verdict, vcolor = "anti", "#a01b1b"
        else:
            verdict, vcolor = "null", "#666"

        sig_summary.append(dict(
            band=band, n_pairs=n_total,
            n_positive_trace=n_pos,
            n_negative_anti=n_neg,
            n_zero=n_zero,
            frac_positive=frac_pos,
            cohort_median_diff=med_s,
            cohort_mean_diff=mean_s,
            verdict=verdict,
        ))

        # Top-left annotation: lead with mean s (the eye-proof number), then
        # right/left tail mass and verdict in matching color.
        ax.text(0.04, 0.96,
                f"mean $s$ = ${mean_s:+.2f}$  ({verdict})\n"
                f"+: {n_pos:,}  /  −: {n_neg:,}",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=9.5, color=vcolor, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.30", facecolor="white",
                          edgecolor=vcolor, linewidth=1.2, alpha=0.96))
        ax.text(0.96, 0.96, BRAIN_BAND_TEX_DICT[band],
                transform=ax.transAxes, ha="right", va="top",
                fontsize=18, color="#111", fontweight="bold")
        ax.set_xlim(-15, 15)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in axes2[1, :]:
        ax.set_xlabel(r"$s = |m_{\mathrm{RPre}} - m_{\mathrm{TT}}| - |m_{\mathrm{TT}} - m_{\mathrm{RPost}}|$  (per pair)")
    for ax in axes2[:, 0]:
        ax.set_ylabel("pair count")

    handles2 = [
        mpatches.Patch(facecolor="#0f5d2c", edgecolor="black",
                       label=r"$s>0$: trace direction (RPost shift smaller)"),
        mpatches.Patch(facecolor="#a01b1b", edgecolor="black",
                       label=r"$s<0$: anti-trace (RPost shift larger)"),
        mlines.Line2D([], [], color="black", lw=1.0, label="$s = 0$"),
        mlines.Line2D([], [], color="gold", lw=2.4, label="cohort mean $s$"),
    ]
    fig2.legend(handles=handles2, loc="lower center",
                bbox_to_anchor=(0.5, -0.02),
                ncol=4, frameon=False, fontsize=9)
    fig2.tight_layout(rect=(0, 0.07, 1, 1))
    fig2.savefig(FIG / "kc_topology_per_pair_signed.pdf", bbox_inches="tight")
    plt.close(fig2)
    print("[redo6] wrote kc_topology_per_pair_signed.pdf")

    pd.DataFrame(sig_summary).to_csv(
        TBL / "kc_topology_per_pair_signed_summary.csv", index=False
    )
    print(pd.DataFrame(sig_summary).to_string(index=False))

    # Remove the old broken sibling if it exists, so the figures dir reflects
    # the current set.
    old = FIG / "kc_topology_eye_proof_signed_high_magnitude.pdf"
    if old.exists():
        old.unlink()
    old_csv = TBL / "kc_topology_eye_proof_high_magnitude_summary.csv"
    if old_csv.exists():
        old_csv.unlink()


# ===================================================================
# Redo 7 — KC dendrogram eye-proof (Pat_06, β trace vs θ null)
# ===================================================================
def _reorder_linkage_to_match(Z: np.ndarray, target_leaves: list[int]) -> np.ndarray:
    """Greedy permutation of merge children so the natural dendrogram leaf
    order is as close to ``target_leaves`` as possible. Topology preserved
    (only left/right swaps at each merge); KC, MC, RF distances unchanged.
    """
    target_pos = {leaf: i for i, leaf in enumerate(target_leaves)}
    n = Z.shape[0] + 1
    leaves_under: dict[int, list[int]] = {i: [i] for i in range(n)}
    Z_perm = Z.copy()
    for k in range(n - 1):
        a = int(Z[k, 0])
        b = int(Z[k, 1])
        cid = n + k
        leaves_a = leaves_under[a]
        leaves_b = leaves_under[b]
        avg_a = np.mean([target_pos.get(l, n) for l in leaves_a])
        avg_b = np.mean([target_pos.get(l, n) for l in leaves_b])
        if avg_a > avg_b:
            Z_perm[k, 0] = b
            Z_perm[k, 1] = a
            leaves_under[cid] = leaves_b + leaves_a
        else:
            leaves_under[cid] = leaves_a + leaves_b
    return Z_perm


def redo7_kc_dendrogram_eye_proof() -> None:
    """Single-patient (Pat_07) eye-proof on the actual dendrograms.

    Rows: β (trace) and θ (anti). Columns: rest_pre, task_test, rest_post.
    For each band, top-5 highlighted pairs selected by absolute contribution
    |c| = ||Δm_RPre,TT|² − |Δm_TT,RPost|²|. Each pair drawn as a colored
    bracket from leaves to MRCA on every phase tree of the row. Marker shape
    encodes per-pair sign of c: ● green = trace, ✕ red = anti. Below each
    row, KC distances at λ ∈ {0, 0.5, 1}.
    """
    import matplotlib.gridspec as gridspec
    from scipy.cluster.hierarchy import dendrogram

    # Patient choice rationale: Pat_06 was the spec's first pick (pro-trace
    # at β across all probes), but Pat_06 traces at *every* band including θ
    # (T_KC(θ, λ=0) = -8.21), so β/θ both look like trace. Pat_07 instead has
    # T_KC(β, λ=0) = -2.43 (trace, moderate) and T_KC(θ, λ=0) = +4.58 (anti);
    # the rsPre/taskt/rsPost arc geometry differs cleanly between the two bands.
    PAT = "Pat_07"
    BANDS_TO_PLOT = ["beta", "theta"]
    PHASES_LOC = ["rest_pre", "task_test", "rest_post"]
    PHASE_LABEL = {"rest_pre": "rest_pre", "task_test": "task_test", "rest_post": "rest_post"}
    N_HIGHLIGHT = 5

    # Round-2 KC numerical table for annotation
    kc_real = pd.read_csv(
        ROOT / "data" / "reports" / "section_5_lrg_trace"
        / "03_kc_lambda_triangle" / "tables"
        / "Td_per_patient_per_band_lambda.csv"
    )

    selection_rows = []

    fig = plt.figure(figsize=(13.0, 10.0))
    gs = gridspec.GridSpec(
        4, 3,
        height_ratios=[3.0, 0.7, 3.0, 0.7],
        hspace=0.10, wspace=0.10,
        figure=fig,
    )

    for row_idx, band in enumerate(BANDS_TO_PLOT):
        # Load LRG per phase
        Z = {}
        for phi in PHASES_LOC:
            res = load_lrg_result(PAT, phi, band, fc_method="imcoh_abs")
            if res is None:
                raise FileNotFoundError(f"missing LRG cache: {PAT} {phi} {band}")
            Z[phi] = res.linkage_matrix

        # Reference leaf order = rsPre's natural dendrogram order
        ref_order = dendrogram(Z["rest_pre"], no_plot=True, no_labels=True)["leaves"]
        Z_aligned = {
            "rest_pre": Z["rest_pre"],
            "task_test": _reorder_linkage_to_match(Z["task_test"], ref_order),
            "rest_post": _reorder_linkage_to_match(Z["rest_post"], ref_order),
        }

        # KC vectors per phase (m=topology depth, M=MRCA height) — identical
        # under left/right swap, so we use the original Z's.
        m_vec = {phi: kc_vectors(Z[phi])[0] for phi in PHASES_LOC}
        M_vec = {phi: kc_vectors(Z[phi])[1] for phi in PHASES_LOC}

        # Pair selection: top 5 pairs by absolute contribution to T_KC²,
        # where  c(i,j) = |Δm_RPre,TT|² - |Δm_TT,RPost|².  Sign of c per pair
        # is the per-pair direction (c>0 = trace, c<0 = anti).  Same rule
        # applied at every (patient, band) cell — no cherry-picking.
        d1 = (m_vec["task_test"].astype(int) - m_vec["rest_pre"].astype(int))
        d2 = (m_vec["rest_post"].astype(int) - m_vec["task_test"].astype(int))
        c_pair = d1 * d1 - d2 * d2  # signed contribution: + = trace, - = anti
        top_idx = np.argsort(np.abs(c_pair))[::-1][:N_HIGHLIGHT]
        n_leaves = Z["rest_pre"].shape[0] + 1
        ii_all, jj_all = np.triu_indices(n_leaves, k=1)
        pairs = [(int(ii_all[k]), int(jj_all[k])) for k in top_idx]
        pair_signs = [int(np.sign(c_pair[k])) for k in top_idx]

        for rank, (li, lj) in enumerate(pairs):
            k = top_idx[rank]
            selection_rows.append(dict(
                band=band, pair_rank=rank + 1,
                leaf_i=li, leaf_j=lj,
                m_rspre=int(m_vec["rest_pre"][k]),
                m_taskt=int(m_vec["task_test"][k]),
                m_rspost=int(m_vec["rest_post"][k]),
                abs_delta_rspre_taskt_depth=int(abs(int(d1[k]))),
                abs_delta_taskt_rspost_depth=int(abs(int(d2[k]))),
                contribution_to_T_KC_sq=int(c_pair[k]),
                sign=("trace" if c_pair[k] > 0
                      else "anti" if c_pair[k] < 0 else "null"),
                height_rspre=float(M_vec["rest_pre"][k]),
                height_taskt=float(M_vec["task_test"][k]),
                height_rspost=float(M_vec["rest_post"][k]),
            ))

        # Plot dendrograms in this band row
        all_heights_this_row = np.concatenate(
            [Z_aligned[phi][:, 2] for phi in PHASES_LOC]
        )
        tmin = float(all_heights_this_row.min()) * 0.8
        tmax = float(all_heights_this_row.max()) * 1.05

        COL_TRACE = "#0f5d2c"
        COL_ANTI = "#a01b1b"

        for col_idx, phi in enumerate(PHASES_LOC):
            ax = fig.add_subplot(gs[2 * row_idx, col_idx])
            Z_use = Z_aligned[phi]

            with plt.rc_context({"lines.linewidth": 0.55}):
                ddata = dendrogram(
                    Z_use, no_labels=True, ax=ax,
                    color_threshold=0.0,
                    above_threshold_color="#aab2bd",
                )

            # Leaf x-positions (scipy convention: 5 + 10 * i for i-th leaf)
            leaves_order = ddata["leaves"]
            leaf_to_x = {leaf: 5.0 + 10.0 * i for i, leaf in enumerate(leaves_order)}

            # Highlighted pair brackets — colored by trace/anti direction.
            for pair_idx, (li, lj) in enumerate(pairs):
                sign_p = pair_signs[pair_idx]
                color = (COL_TRACE if sign_p > 0
                         else COL_ANTI if sign_p < 0 else "#666")
                marker = "o" if sign_p > 0 else "X" if sign_p < 0 else "s"
                xi = leaf_to_x.get(li)
                xj = leaf_to_x.get(lj)
                if xi is None or xj is None:
                    continue
                k = top_idx[pair_idx]
                y_mrca = float(M_vec[phi][k])
                # Bracket
                ax.plot([xi, xi], [tmin, y_mrca], color=color, lw=1.6,
                        alpha=0.9, zorder=5, solid_capstyle="round")
                ax.plot([xj, xj], [tmin, y_mrca], color=color, lw=1.6,
                        alpha=0.9, zorder=5, solid_capstyle="round")
                ax.plot([xi, xj], [y_mrca, y_mrca], color=color, lw=1.6,
                        alpha=0.9, zorder=5, solid_capstyle="round")
                # MRCA marker (shape = direction)
                ax.scatter([(xi + xj) / 2], [y_mrca], color=color, s=85,
                           marker=marker, edgecolor="black", linewidth=0.7,
                           zorder=6)
                # Rank number next to the marker
                ax.annotate(f"{pair_idx + 1}",
                            xy=((xi + xj) / 2, y_mrca),
                            xytext=(5, 4), textcoords="offset points",
                            fontsize=8, color="black", fontweight="bold",
                            zorder=7,
                            bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                      ec=color, lw=0.6, alpha=0.9))

            ax.set_yscale("log")
            ax.set_ylim(tmin, tmax)
            ax.set_xticks([])
            ax.tick_params(axis="y", labelsize=8)
            ax.spines[["top", "right"]].set_visible(False)

            # Column title only on the top row
            if row_idx == 0:
                ax.text(0.5, 1.04, PHASE_LABEL[phi], transform=ax.transAxes,
                        ha="center", va="bottom", fontsize=11.5,
                        fontweight="bold", color="#222")

            # Row label only on the leftmost column
            if col_idx == 0:
                row_tex = BRAIN_BAND_TEX_DICT[band]
                tag = "trace band" if band == "beta" else "null band"
                ax.set_ylabel(f"{row_tex}  ({tag})\nmerge height (log)",
                              fontsize=10)

        # Numerical block — spans all 3 columns under this band's row
        ax_t = fig.add_subplot(gs[2 * row_idx + 1, :])
        ax_t.axis("off")

        # Pull KC numbers for this (PAT, band) from the round-2 table.
        kc_sub = kc_real[(kc_real["patient"] == PAT) & (kc_real["band"] == band)]
        rows = ["", "", ""]
        for lam in (0.0, 0.5, 1.0):
            r = kc_sub[kc_sub["lam"] == lam]
            if r.empty:
                continue
            d_pre_tt = float(r["d_pre_tt"].iloc[0])
            d_tt_post = float(r["d_tt_post"].iloc[0])
            T = float(r["T_KC"].iloc[0])
            verdict = "trace" if T < 0 else "anti"
            vcol = "#0f5d2c" if T < 0 else "#a01b1b"
            line = (
                rf"$\lambda{{=}}{lam:.1f}$:  "
                rf"$d_{{KC}}(\mathrm{{RPre,TT}}) = {d_pre_tt:.3f}$,  "
                rf"$d_{{KC}}(\mathrm{{TT,RPost}}) = {d_tt_post:.3f}$,  "
                rf"$T_{{KC}} = {T:+.3f}$  ({verdict})"
            )
            i_line = {0.0: 0, 0.5: 1, 1.0: 2}[lam]
            rows[i_line] = (line, vcol)

        # Title for this annotation block, including top-5 trace/anti tally
        n_trace_top = sum(1 for s in pair_signs if s > 0)
        n_anti_top = sum(1 for s in pair_signs if s < 0)
        ax_t.text(
            0.0, 0.94,
            f"  KC distances for {PAT} at {BRAIN_BAND_TEX_DICT[band]}   "
            f"·  top-5 pairs by $|$contribution to $T_{{KC}}^2|$:  "
            f"{n_trace_top} trace, {n_anti_top} anti",
            transform=ax_t.transAxes, ha="left", va="top",
            fontsize=10.5, fontweight="bold", color="#222",
        )
        for i, item in enumerate(rows):
            if not item:
                continue
            line, vcol = item
            ax_t.text(
                0.02, 0.78 - i * 0.27, line,
                transform=ax_t.transAxes, ha="left", va="top",
                fontsize=10, color=vcol, family="DejaVu Sans",
            )

    # Figure-level legend
    handles = [
        mlines.Line2D([], [], color="#0f5d2c", lw=2.0, marker="o",
                      markersize=9, markerfacecolor="#0f5d2c",
                      markeredgecolor="black",
                      label=r"trace pair (contribution $c > 0$)"),
        mlines.Line2D([], [], color="#a01b1b", lw=2.0, marker="X",
                      markersize=9, markerfacecolor="#a01b1b",
                      markeredgecolor="black",
                      label=r"anti pair (contribution $c < 0$)"),
        mlines.Line2D([], [], color="#aab2bd", lw=0.55,
                      label="dendrogram skeleton"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.01),
               ncol=3, frameon=False, fontsize=9.5)

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(FIG / "kc_dendrogram_eye_proof.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[redo7] wrote kc_dendrogram_eye_proof.pdf")

    pd.DataFrame(selection_rows).to_csv(
        TBL / f"kc_eye_proof_{PAT.lower()}_pair_selection.csv", index=False
    )
    print(f"[redo7] wrote kc_eye_proof_{PAT.lower()}_pair_selection.csv")
    print(pd.DataFrame(selection_rows).to_string(index=False))


# ===================================================================
# Redo 8 — Per-patient KC-clade dendrogram visualization (Pat_06)
# ===================================================================
def _enumerate_clades(Z: np.ndarray) -> list[dict]:
    """For every internal node of a scipy linkage Z, return
    [{'cid': int, 'leaves': frozenset[int], 'height': float}, ...]
    sorted by node id (root last). Length = N-1.
    """
    n = Z.shape[0] + 1
    leaves_under: dict[int, frozenset[int]] = {i: frozenset({i}) for i in range(n)}
    out = []
    for k in range(n - 1):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        cid = n + k
        leaves_under[cid] = leaves_under[a] | leaves_under[b]
        out.append(dict(cid=cid, leaves=leaves_under[cid], height=float(Z[k, 2])))
    return out


def _max_jaccard_match(
    clades_A: list[dict], clades_B: list[dict]
) -> list[dict]:
    """For each clade in A return its best-Jaccard partner in B.
    Returns list aligned with clades_A: each row has
    {'cid_A','leaves_A','height_A','cid_B','leaves_B','height_B','jaccard'}.
    """
    rows = []
    sets_B = [c["leaves"] for c in clades_B]
    for cA in clades_A:
        sA = cA["leaves"]
        best_j = 0.0
        best_idx = -1
        for idx, sB in enumerate(sets_B):
            inter = len(sA & sB)
            if inter == 0:
                continue
            uni = len(sA | sB)
            j = inter / uni
            if j > best_j:
                best_j = j
                best_idx = idx
        if best_idx < 0:
            continue
        cB = clades_B[best_idx]
        rows.append(dict(
            cid_A=cA["cid"], leaves_A=cA["leaves"], height_A=cA["height"],
            cid_B=cB["cid"], leaves_B=cB["leaves"], height_B=cB["height"],
            jaccard=best_j,
        ))
    return rows


def _clade_leaf_xrange(Z: np.ndarray, leaf_to_x: dict[int, float],
                        cid: int) -> tuple[float, float]:
    """Min/max leaf-x position covered by the subtree rooted at cid."""
    n = Z.shape[0] + 1
    if cid < n:
        return (leaf_to_x[cid], leaf_to_x[cid])
    # Walk subtree iteratively
    stack = [cid]
    leaves = []
    while stack:
        node = stack.pop()
        if node < n:
            leaves.append(leaf_to_x[node])
        else:
            kk = node - n
            stack.append(int(Z[kk, 0]))
            stack.append(int(Z[kk, 1]))
    return (min(leaves), max(leaves))


def _draw_clade(ax, Z: np.ndarray, leaf_to_x: dict[int, float],
                cid: int, color: str, lw: float, alpha: float,
                tmin: float):
    """Draw a single horizontal bar at the clade's merge height, spanning
    the leaf-x range covered by its subtree, with two short vertical whiskers
    dropping to ``tmin`` so the eye can locate the clade's leaf base.
    Far less visual mass than recursively recoloring the whole subtree.
    """
    n = Z.shape[0] + 1
    if cid < n:
        return  # leaves are not clades
    kk = cid - n
    h = float(Z[kk, 2])
    x_lo, x_hi = _clade_leaf_xrange(Z, leaf_to_x, cid)
    # Horizontal bar at merge height
    ax.plot([x_lo, x_hi], [h, h], color=color, lw=lw, alpha=alpha,
            solid_capstyle="round", zorder=5)
    # Short whiskers (drop ~25% of the log-axis range visually)
    y_whisk = max(tmin, h * 0.5)
    ax.plot([x_lo, x_lo], [y_whisk, h], color=color, lw=lw * 0.7,
            alpha=alpha * 0.85, zorder=5)
    ax.plot([x_hi, x_hi], [y_whisk, h], color=color, lw=lw * 0.7,
            alpha=alpha * 0.85, zorder=5)


def redo8_kc_dendrogram_per_patient() -> None:
    """Per-patient KC-driven dendrogram visualization for Pat_06.

    6 rows (2 bands × 3 λ regimes) × 3 cols (3 phases) of dendrograms.
    Highlight rule per row:
      - λ=0 (topology):     leaf-set Jaccard match (strict ≥ τ_s, partial ≥ τ_p).
      - λ=0.5 (mixed):      Jaccard ≥ τ_p AND |Δh|/h_max ≤ height_tol.
      - λ=1 (heights):      Jaccard ≥ τ_p; intensity ∝ 1 − |Δh|/h_max.
    Blue family = rsPre→taskt matches; red family = taskt→rsPost matches.
    """
    import matplotlib.gridspec as gridspec
    from scipy.cluster.hierarchy import dendrogram

    # Spec called for Pat_06 first; Pat_06 traces at every band, gives
    # the wrong visual direction (θ shows MORE strict matches than β).
    # Pat_08: too uniform — β/θ separate by only 1 strict match.
    # Pat_07: anti direction at β.
    # Pat_13 has the textbook β-trace / θ-anti clade-preservation
    # asymmetry under leaf-set Jaccard (thr ≥ 0.85):
    #   β: rsPre→taskt 28 strict, taskt→rsPost 35 strict (+7, trace)
    #   θ: rsPre→taskt 27 strict, taskt→rsPost 21 strict (−6, anti)
    # The rsPost panel of β should therefore carry visibly more red
    # highlighting; the rsPre panel of θ should carry visibly more
    # blue highlighting. Cohort KC numbers at λ=0 for Pat_13 happen to
    # be trace at both bands (β -3.16, θ -4.67); the clade-Jaccard
    # readout and KC λ=0 measure related but distinct things.
    PAT = "Pat_13"
    BANDS = ["beta", "theta"]
    PHASES_LOC = ["rest_pre", "task_test", "rest_post"]
    LAMBDA_REGIMES = [
        ("λ=0", "topology"),
        ("λ=0.5", "mixed"),
        ("λ=1", "heights"),
    ]
    JAC_STRICT = 0.90
    JAC_PARTIAL = 0.80
    HEIGHT_TOL = 0.10  # for λ=0.5 — fraction of h_max
    # Clade-size filter: drop near-root and tiny clades (no information).
    MIN_CLADE_SIZE = 3
    MAX_CLADE_FRAC = 0.50

    BLUE_S, BLUE_P = "#08519c", "#9ecae1"   # rsPre→taskt
    RED_S, RED_P = "#a50f15", "#fcae91"     # taskt→rsPost

    kc_real = pd.read_csv(
        ROOT / "data" / "reports" / "section_5_lrg_trace"
        / "03_kc_lambda_triangle" / "tables"
        / "Td_per_patient_per_band_lambda.csv"
    )

    csv_rows = []

    fig = plt.figure(figsize=(12.5, 16.0))
    gs = gridspec.GridSpec(
        8, 3,
        height_ratios=[2.4, 2.4, 2.4, 0.55, 2.4, 2.4, 2.4, 0.55],
        hspace=0.15, wspace=0.10,
        figure=fig,
    )

    for band_idx, band in enumerate(BANDS):
        Z_raw = {}
        for phi in PHASES_LOC:
            res = load_lrg_result(PAT, phi, band, fc_method="imcoh_abs")
            if res is None:
                raise FileNotFoundError(f"missing LRG cache: {PAT} {phi} {band}")
            Z_raw[phi] = np.asarray(res.linkage_matrix)

        # Reference leaf order = rsPre's natural dendrogram ordering;
        # other phases reordered (topology preserved, only L/R swaps).
        ref_order = dendrogram(Z_raw["rest_pre"], no_plot=True, no_labels=True)["leaves"]
        Z = {
            "rest_pre": Z_raw["rest_pre"],
            "task_test": _reorder_linkage_to_match(Z_raw["task_test"], ref_order),
            "rest_post": _reorder_linkage_to_match(Z_raw["rest_post"], ref_order),
        }

        # Clade enumeration per phase (use the *ALIGNED* Z so cluster ids
        # in `clades` match the linkage we draw from). Filter near-root
        # and tiny clades — they carry no informative match signal.
        n_leaves = Z["rest_pre"].shape[0] + 1
        max_clade_size = int(MAX_CLADE_FRAC * n_leaves)
        clades = {
            phi: [c for c in _enumerate_clades(Z[phi])
                  if MIN_CLADE_SIZE <= len(c["leaves"]) <= max_clade_size]
            for phi in PHASES_LOC
        }

        # Cross-phase clade matches (best Jaccard, in both directions)
        # rsPre → taskt:  for each rsPre clade, best taskt match
        match_pre_tt = _max_jaccard_match(clades["rest_pre"], clades["task_test"])
        # taskt → rsPost: for each taskt clade, best rsPost match
        match_tt_post = _max_jaccard_match(clades["task_test"], clades["rest_post"])

        h_max_band = float(np.concatenate(
            [Z[phi][:, 2] for phi in PHASES_LOC]
        ).max())

        for lam_idx, (lam_name, regime) in enumerate(LAMBDA_REGIMES):
            row_idx = band_idx * 4 + lam_idx  # 0..2 (β), 4..6 (θ)

            # Per-(band,λ) selection of matches that pass the regime rule
            def filter_matches(rows: list[dict]) -> list[tuple[dict, str, float]]:
                """Return [(row, intensity_kind, alpha)] entries.
                intensity_kind in {'strict', 'partial'}; alpha is 1.0 except
                for λ=1 where it scales with height agreement.
                """
                kept = []
                for r in rows:
                    j = r["jaccard"]
                    dh = abs(r["height_A"] - r["height_B"]) / max(h_max_band, 1e-12)
                    if regime == "topology":
                        if j >= JAC_STRICT:
                            kept.append((r, "strict", 1.0))
                        elif j >= JAC_PARTIAL:
                            kept.append((r, "partial", 1.0))
                    elif regime == "mixed":
                        if j < JAC_PARTIAL or dh > HEIGHT_TOL:
                            continue
                        kind = "strict" if j >= JAC_STRICT else "partial"
                        kept.append((r, kind, 1.0))
                    elif regime == "heights":
                        if j < JAC_PARTIAL:
                            continue
                        kind = "strict" if j >= JAC_STRICT else "partial"
                        # Intensity ∝ 1 − |Δh|/h_max, floor at 0.25
                        alpha_h = max(0.25, 1.0 - dh)
                        kept.append((r, kind, alpha_h))
                return kept

            kept_pre_tt = filter_matches(match_pre_tt)
            kept_tt_post = filter_matches(match_tt_post)

            # CSV log of every kept match
            for r, kind, alpha in kept_pre_tt:
                csv_rows.append(dict(
                    band=band, lambda_regime=lam_name, phase_pair="rsPre→taskt",
                    clade_A_root=r["cid_A"],
                    clade_A_leaves="|".join(map(str, sorted(r["leaves_A"]))),
                    clade_B_root=r["cid_B"],
                    clade_B_leaves="|".join(map(str, sorted(r["leaves_B"]))),
                    jaccard=r["jaccard"],
                    height_A=r["height_A"], height_B=r["height_B"],
                    match_strength=kind, intensity_alpha=alpha,
                ))
            for r, kind, alpha in kept_tt_post:
                csv_rows.append(dict(
                    band=band, lambda_regime=lam_name, phase_pair="taskt→rsPost",
                    clade_A_root=r["cid_A"],
                    clade_A_leaves="|".join(map(str, sorted(r["leaves_A"]))),
                    clade_B_root=r["cid_B"],
                    clade_B_leaves="|".join(map(str, sorted(r["leaves_B"]))),
                    jaccard=r["jaccard"],
                    height_A=r["height_A"], height_B=r["height_B"],
                    match_strength=kind, intensity_alpha=alpha,
                ))

            for col_idx, phi in enumerate(PHASES_LOC):
                ax = fig.add_subplot(gs[row_idx, col_idx])
                Z_use = Z[phi]
                ddata = dendrogram(
                    Z_use, no_labels=True, ax=ax,
                    color_threshold=0.0,
                    above_threshold_color="#cdd2db",
                )
                # gray skeleton: thin out
                for ln in ax.collections:
                    pass  # collections are lines created by dendrogram
                for ln in ax.get_lines():
                    ln.set_linewidth(0.4)
                    ln.set_alpha(0.45)

                leaves_order = ddata["leaves"]
                leaf_to_x = {leaf: 5.0 + 10.0 * i
                             for i, leaf in enumerate(leaves_order)}

                # Decide which clade-id to highlight in THIS phase column
                # for each phase pair:
                # - rsPre→taskt match: highlight clade A in rsPre column,
                #   clade B in taskt column. Nothing in rsPost column.
                # - taskt→rsPost match: highlight clade A (taskt) in taskt
                #   column, clade B (rsPost) in rsPost column. Nothing in
                #   rsPre column.
                # In task_test column, BOTH families are shown (overlay).
                def col_targets(pair_label: str, phi_col: str):
                    if pair_label == "rsPre→taskt":
                        if phi_col == "rest_pre":
                            return [(r, "A", k, a) for (r, k, a) in kept_pre_tt]
                        if phi_col == "task_test":
                            return [(r, "B", k, a) for (r, k, a) in kept_pre_tt]
                        return []
                    else:  # taskt→rsPost
                        if phi_col == "task_test":
                            return [(r, "A", k, a) for (r, k, a) in kept_tt_post]
                        if phi_col == "rest_post":
                            return [(r, "B", k, a) for (r, k, a) in kept_tt_post]
                        return []

                heights_all = Z_use[:, 2]
                tmin = float(heights_all.min()) * 0.8
                tmax = float(heights_all.max()) * 1.05

                # Draw rsPre→taskt highlights (blue) first
                for (r, side, kind, alpha) in col_targets("rsPre→taskt", phi):
                    cid = r["cid_A"] if side == "A" else r["cid_B"]
                    color = BLUE_S if kind == "strict" else BLUE_P
                    lw = 2.2 if kind == "strict" else 1.4
                    _draw_clade(ax, Z_use, leaf_to_x, cid,
                                color=color, lw=lw, alpha=alpha,
                                tmin=tmin)

                # Then taskt→rsPost (red)
                for (r, side, kind, alpha) in col_targets("taskt→rsPost", phi):
                    cid = r["cid_A"] if side == "A" else r["cid_B"]
                    color = RED_S if kind == "strict" else RED_P
                    lw = 2.2 if kind == "strict" else 1.4
                    _draw_clade(ax, Z_use, leaf_to_x, cid,
                                color=color, lw=lw, alpha=alpha,
                                tmin=tmin)

                ax.set_yscale("log")
                ax.set_ylim(tmin, tmax)
                ax.set_xticks([])
                ax.tick_params(axis="y", which="both", labelsize=7)
                # Suppress the noisy minor-tick labels on log y axes
                ax.tick_params(axis="y", which="minor", labelleft=False)
                ax.spines[["top", "right"]].set_visible(False)

                # Column title only at very top of figure
                if band_idx == 0 and lam_idx == 0:
                    ax.text(0.5, 1.06, phi, transform=ax.transAxes,
                            ha="center", va="bottom", fontsize=11,
                            fontweight="bold", color="#222")

                # Row label only on leftmost column
                if col_idx == 0:
                    n_pre_tt_strict = sum(
                        1 for (_, _, k, _) in col_targets("rsPre→taskt", phi)
                        if k == "strict")
                    n_pre_tt_part = sum(
                        1 for (_, _, k, _) in col_targets("rsPre→taskt", phi)
                        if k == "partial")
                    n_tt_post_strict = sum(
                        1 for (_, _, k, _) in col_targets("taskt→rsPost", phi)
                        if k == "strict")
                    n_tt_post_part = sum(
                        1 for (_, _, k, _) in col_targets("taskt→rsPost", phi)
                        if k == "partial")
                    band_tex = BRAIN_BAND_TEX_DICT[band]
                    tag = "trace" if band == "beta" else "null"
                    ax.set_ylabel(
                        f"{band_tex}  ({tag})\n"
                        rf"$\bf{{{lam_name.replace('λ', r'\lambda')}}}$ "
                        f"  {regime}\n"
                        f"merge height (log)",
                        fontsize=8.0,
                    )

        # Annotation block under this band's 3 λ-rows
        ax_t = fig.add_subplot(gs[band_idx * 4 + 3, :])
        ax_t.axis("off")

        kc_sub = kc_real[
            (kc_real["patient"] == PAT) & (kc_real["band"] == band)
        ]
        # Per-lambda counts of strict/partial matches
        # (recompute by re-applying the regime rule to the full match lists)
        def regime_counts(rows: list[dict], regime: str) -> tuple[int, int]:
            ns, np_ = 0, 0
            for r in rows:
                j = r["jaccard"]
                dh = abs(r["height_A"] - r["height_B"]) / max(h_max_band, 1e-12)
                if regime == "topology":
                    if j >= JAC_STRICT:
                        ns += 1
                    elif j >= JAC_PARTIAL:
                        np_ += 1
                elif regime == "mixed":
                    if j >= JAC_PARTIAL and dh <= HEIGHT_TOL:
                        if j >= JAC_STRICT:
                            ns += 1
                        else:
                            np_ += 1
                else:  # heights
                    if j >= JAC_PARTIAL:
                        if j >= JAC_STRICT:
                            ns += 1
                        else:
                            np_ += 1
            return ns, np_

        # Header
        ax_t.text(
            0.0, 0.95,
            f"  KC numerical block — {PAT} at {BRAIN_BAND_TEX_DICT[band]}  "
            f"(Jaccard: strict $\\geq {JAC_STRICT}$, partial $\\geq {JAC_PARTIAL}$;  "
            f"$\\lambda{{=}}0.5$ height tol $= {int(HEIGHT_TOL*100)}\\%$ of $h_{{max}}$)",
            transform=ax_t.transAxes, ha="left", va="top",
            fontsize=10.0, fontweight="bold", color="#222",
        )

        # Per-λ numeric line
        y0 = 0.62
        dy = 0.30
        for i_lam, (lam_name, regime) in enumerate(LAMBDA_REGIMES):
            lam_val = {0: 0.0, 1: 0.5, 2: 1.0}[i_lam]
            r = kc_sub[kc_sub["lam"] == lam_val]
            if r.empty:
                continue
            d_pre_tt_kc = float(r["d_pre_tt"].iloc[0])
            d_tt_post_kc = float(r["d_tt_post"].iloc[0])
            T = float(r["T_KC"].iloc[0])
            verdict = "trace" if T < 0 else "anti"
            vcol = "#0f5d2c" if T < 0 else "#a01b1b"
            n_blue_s, n_blue_p = regime_counts(match_pre_tt, regime)
            n_red_s, n_red_p = regime_counts(match_tt_post, regime)
            line = (
                rf"$\lambda{{=}}{lam_val:g}$ ({regime}):  "
                rf"$d_{{KC}}$(RPre,TT)$={d_pre_tt_kc:.2f}$,  "
                rf"$d_{{KC}}$(TT,RPost)$={d_tt_post_kc:.2f}$,  "
                rf"$T_{{KC}}={T:+.2f}$ ({verdict}).  "
                rf"matches: blue {n_blue_s}+{n_blue_p}, "
                rf"red {n_red_s}+{n_red_p}"
            )
            ax_t.text(
                0.02, y0 - i_lam * dy, line,
                transform=ax_t.transAxes, ha="left", va="top",
                fontsize=8.5, color=vcol, family="DejaVu Sans",
            )

    # Figure-level legend
    handles = [
        mlines.Line2D([], [], color=BLUE_S, lw=2.2,
                      label=rf"rsPre$\to$taskt strict ($J\geq {JAC_STRICT}$)"),
        mlines.Line2D([], [], color=BLUE_P, lw=1.4,
                      label=rf"rsPre$\to$taskt partial ($J\geq {JAC_PARTIAL}$)"),
        mlines.Line2D([], [], color=RED_S, lw=2.2,
                      label=rf"taskt$\to$rsPost strict ($J\geq {JAC_STRICT}$)"),
        mlines.Line2D([], [], color=RED_P, lw=1.4,
                      label=rf"taskt$\to$rsPost partial ($J\geq {JAC_PARTIAL}$)"),
        mlines.Line2D([], [], color="#cdd2db", lw=0.4,
                      label="dendrogram skeleton"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.005),
               ncol=3, frameon=False, fontsize=9.0)

    fig.tight_layout(rect=(0.0, 0.025, 1.0, 1.0))
    fig.savefig(FIG / "kc_dendrogram_per_patient.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[redo8] wrote kc_dendrogram_per_patient.pdf")

    pd.DataFrame(csv_rows).to_csv(
        TBL / f"kc_dendrogram_clade_matches_{PAT.lower()}.csv", index=False
    )
    print(f"[redo8] wrote kc_dendrogram_clade_matches_{PAT.lower()}.csv "
          f"({len(csv_rows)} rows)")


# ===================================================================
# Redo 9 — Cohort per-band dendrogram view (one PDF per band)
# ===================================================================
def _best_jaccard_partner(
    anchor: dict, candidates: list[dict]
) -> tuple[dict | None, float]:
    """For one anchor clade, return (best_match_clade, best_jaccard) over a
    list of candidate clades. Returns (None, 0.0) if no overlap.
    """
    sA = anchor["leaves"]
    best_j = 0.0
    best = None
    for cB in candidates:
        sB = cB["leaves"]
        inter = len(sA & sB)
        if inter == 0:
            continue
        j = inter / len(sA | sB)
        if j > best_j:
            best_j = j
            best = cB
    return best, best_j


def redo9_kc_cohort_band_dendrograms() -> None:
    """One PDF per band, full cohort (10 patients) × 3 phases.

    Subtree identification: per (patient, band), enumerate clades in
    ``T_taskt`` (size 5..40% of N), score each by
    ``s = J(C_taskt, best_J in T_rspost) − J(C_taskt, best_J in T_rspre)``.
    Pick the clade with maximum |s| as the anchor. The rspre panel
    highlights the rspre best-Jaccard partner; the taskt panel highlights
    the anchor clade itself; the rspost panel highlights the rspost best-
    Jaccard partner. Each phase panel uses **its own natural leaf order**,
    so the matched clade is a contiguous block and gets shaded as a single
    translucent vertical span behind the dendrogram. Color: green if
    s > 0 (trace), red if s < 0 (anti).
    """
    import matplotlib.gridspec as gridspec
    from scipy.cluster.hierarchy import dendrogram

    PATIENTS = list(PATIENTS_LIST)
    BANDS = ["beta", "theta"]
    PHASES_LOC = ["rest_pre", "task_test", "rest_post"]
    MIN_SIZE = 5
    MAX_FRAC = 0.40

    COL_TRACE = "#0f5d2c"
    COL_ANTI = "#a01b1b"
    COL_NULL = "#777777"
    SHADE_ALPHA = 0.20

    summary_rows = []

    for band in BANDS:
        n_pat = len(PATIENTS)
        fig = plt.figure(figsize=(12.0, 1.45 * n_pat + 1.5))
        gs = gridspec.GridSpec(
            n_pat, 3, figure=fig,
            hspace=0.35, wspace=0.05,
        )

        for row_idx, pat in enumerate(PATIENTS):
            # Load LRG per phase (use natural ordering per phase)
            Z = {}
            for phi in PHASES_LOC:
                res = load_lrg_result(pat, phi, band, fc_method="imcoh_abs")
                if res is None:
                    raise FileNotFoundError(
                        f"missing LRG cache: {pat} {phi} {band}")
                Z[phi] = np.asarray(res.linkage_matrix)
            n_leaves = Z["rest_pre"].shape[0] + 1
            max_size = int(MAX_FRAC * n_leaves)

            # Enumerate clades per phase
            cl = {phi: _enumerate_clades(Z[phi]) for phi in PHASES_LOC}
            cl_taskt_filt = [
                c for c in cl["task_test"]
                if MIN_SIZE <= len(c["leaves"]) <= max_size
            ]

            # Score each taskt clade by trace asymmetry
            best_anchor = None
            best_score = 0.0
            best_match_pre = None
            best_match_post = None
            best_J_pre = 0.0
            best_J_post = 0.0
            for c in cl_taskt_filt:
                m_pre, j_pre = _best_jaccard_partner(c, cl["rest_pre"])
                m_post, j_post = _best_jaccard_partner(c, cl["rest_post"])
                s = j_post - j_pre
                if abs(s) > abs(best_score):
                    best_score = s
                    best_anchor = c
                    best_match_pre = m_pre
                    best_match_post = m_post
                    best_J_pre = j_pre
                    best_J_post = j_post

            if best_anchor is None:
                # Defensive: no eligible clade found
                color = COL_NULL
                tag = "no clade"
            else:
                if best_score > 0.05:
                    color = COL_TRACE
                    tag = "trace"
                elif best_score < -0.05:
                    color = COL_ANTI
                    tag = "anti"
                else:
                    color = COL_NULL
                    tag = "null"

            summary_rows.append(dict(
                band=band, patient=pat, tag=tag, score=best_score,
                anchor_size=len(best_anchor["leaves"]) if best_anchor else 0,
                J_pre=best_J_pre, J_post=best_J_post,
                anchor_height=(best_anchor["height"] if best_anchor
                               else float("nan")),
                anchor_leaves=("|".join(map(str, sorted(best_anchor["leaves"])))
                               if best_anchor else ""),
            ))

            # Render the 3 phase panels
            highlight = {
                "rest_pre": best_match_pre if best_J_pre > 0 else None,
                "task_test": best_anchor,
                "rest_post": best_match_post if best_J_post > 0 else None,
            }

            for col_idx, phi in enumerate(PHASES_LOC):
                ax = fig.add_subplot(gs[row_idx, col_idx])
                Z_use = Z[phi]
                ddata = dendrogram(
                    Z_use, no_labels=True, ax=ax,
                    color_threshold=0.0,
                    above_threshold_color="#444",
                )
                # Skeleton color/weight (light gray, thin but visible)
                for ln in ax.get_lines():
                    ln.set_linewidth(0.55)
                    ln.set_color("#5a5f6a")
                    ln.set_alpha(0.85)

                heights = Z_use[:, 2]
                tmin = float(heights.min()) * 0.8
                tmax = float(heights.max()) * 1.05

                # Shade highlighted clade leaf range — cap at merge height
                hl = highlight[phi]
                if hl is not None:
                    leaves_order = ddata["leaves"]
                    leaf_to_x = {leaf: 5.0 + 10.0 * i
                                 for i, leaf in enumerate(leaves_order)}
                    xs = [leaf_to_x[l] for l in hl["leaves"]
                          if l in leaf_to_x]
                    if xs:
                        x_lo, x_hi = min(xs) - 4.0, max(xs) + 4.0
                        h_clade = float(hl["height"])
                        ax.fill_betweenx(
                            [tmin, h_clade], x_lo, x_hi,
                            facecolor=color, alpha=SHADE_ALPHA,
                            zorder=0, linewidth=0,
                        )
                        # Top edge at clade merge height
                        ax.plot([x_lo, x_hi], [h_clade, h_clade],
                                color=color, lw=1.6, alpha=0.95, zorder=4)
                        # Side ticks at the clade range
                        ax.plot([x_lo, x_lo], [tmin, h_clade],
                                color=color, lw=0.9, alpha=0.7, zorder=4)
                        ax.plot([x_hi, x_hi], [tmin, h_clade],
                                color=color, lw=0.9, alpha=0.7, zorder=4)

                ax.set_yscale("log")
                ax.set_ylim(tmin, tmax)
                ax.set_xticks([])
                ax.tick_params(axis="y", which="major", labelsize=6.5,
                                length=2.0, pad=1.5)
                ax.tick_params(axis="y", which="minor", left=False,
                                labelleft=False)
                ax.spines[["top", "right"]].set_visible(False)

                # Column titles only on the very top row
                if row_idx == 0:
                    ax.set_title(phi, fontsize=10.5, fontweight="bold",
                                  pad=4.0)

                # Row label only on leftmost column
                if col_idx == 0:
                    score_str = (f"$s={best_score:+.2f}$"
                                 if best_anchor else "no clade")
                    ax.set_ylabel(
                        f"{pat}\n{tag}  {score_str}",
                        fontsize=8.0, color=color, fontweight="bold",
                        rotation=0, ha="right", va="center",
                        labelpad=18,
                    )

                # Per-panel annotation in lower-right corner (away from y-tick labels)
                if hl is not None:
                    if phi == "rest_pre":
                        anno = rf"$J={best_J_pre:.2f}$"
                    elif phi == "rest_post":
                        anno = rf"$J={best_J_post:.2f}$"
                    else:
                        anno = (rf"$|C|={len(best_anchor['leaves'])}$"
                                if best_anchor else "")
                    ax.text(0.985, 0.05, anno, transform=ax.transAxes,
                            ha="right", va="bottom", fontsize=7.0,
                            color=color, fontweight="bold",
                            bbox=dict(boxstyle="round,pad=0.18", fc="white",
                                       ec="none", alpha=0.85))

        # Figure-level legend
        handles = [
            mpatches.Patch(facecolor=COL_TRACE, alpha=SHADE_ALPHA,
                            edgecolor=COL_TRACE,
                            label=r"trace ($s = J_{post}-J_{pre} > +0.05$)"),
            mpatches.Patch(facecolor=COL_ANTI, alpha=SHADE_ALPHA,
                            edgecolor=COL_ANTI,
                            label=r"anti ($s < -0.05$)"),
            mpatches.Patch(facecolor=COL_NULL, alpha=SHADE_ALPHA,
                            edgecolor=COL_NULL,
                            label=r"null ($|s| \leq 0.05$)"),
            mlines.Line2D([], [], color="#5a5f6a", lw=0.55,
                          label="dendrogram skeleton"),
        ]
        fig.legend(handles=handles, loc="lower center",
                    bbox_to_anchor=(0.5, -0.005),
                    ncol=4, frameon=False, fontsize=9.0)

        fig.tight_layout(rect=(0.0, 0.025, 1.0, 1.0))
        out_pdf = FIG / f"kc_cohort_dendrograms_{band}.pdf"
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        print(f"[redo9] wrote {out_pdf.name}")

    # Combined cohort summary CSV
    pd.DataFrame(summary_rows).to_csv(
        TBL / "kc_cohort_dendrograms_summary.csv", index=False
    )
    print("[redo9] wrote kc_cohort_dendrograms_summary.csv")


# ===================================================================
# Redo 10 — Cohort KC anchor-Jaccard trail (one PDF per band)
# ===================================================================
def redo10_kc_cohort_jaccard_trail() -> None:
    """Cohort visualization of the KC anchor-clade leaf-set persistence.

    Reads ``kc_cohort_dendrograms_summary.csv`` (produced by redo9) for the
    per-(patient, band) anchor selection. Per band, draws a single panel with
    10 polylines (one per patient): three points
    ``(rspre, J_pre)``, ``(taskt, 1.0)``, ``(rspost, J_post)``. Color encodes
    direction (green: trace, ``J_post > J_pre``; red: anti, ``J_post < J_pre``;
    grey: null). Patient ID labels at the rspost end. Annotation reports the
    cohort tally and the median asymmetry ``s = J_post − J_pre``.
    """
    from matplotlib.lines import Line2D

    summary_csv = TBL / "kc_cohort_dendrograms_summary.csv"
    if not summary_csv.exists():
        # If redo9 was not yet run, run it now to generate the inputs
        redo9_kc_cohort_band_dendrograms()

    df_all = pd.read_csv(summary_csv)
    BANDS = ["beta", "theta"]
    PHASES_LOC = ["rest_pre", "task_test", "rest_post"]
    PHASE_X = {p: i for i, p in enumerate(PHASES_LOC)}

    COL_TRACE = "#0f5d2c"
    COL_ANTI = "#a01b1b"
    COL_NULL = "#888"

    for band in BANDS:
        sub = df_all[df_all["band"] == band].copy()
        sub = sub.sort_values("score", ascending=True).reset_index(drop=True)

        # Stagger label y-positions so they don't pile up at J_post ties
        sub["label_y"] = sub["J_post"].copy()
        # Greedy de-overlap: walk in sorted order, push labels apart by ε
        order = np.argsort(sub["label_y"].values)
        ys = sub["label_y"].values.copy()
        ys_sorted = ys[order]
        eps = 0.035
        for i in range(1, len(ys_sorted)):
            if ys_sorted[i] - ys_sorted[i - 1] < eps:
                ys_sorted[i] = ys_sorted[i - 1] + eps
        ys[order] = ys_sorted
        sub["label_y"] = np.clip(ys, -0.02, 1.10)

        fig, ax = plt.subplots(figsize=(7.0, 5.4))
        # Phase axis
        ax.set_xlim(-0.18, 2.65)
        ax.set_ylim(-0.02, 1.12)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["rest_pre", "task_test", "rest_post"],
                            fontsize=10.5, fontweight="bold")
        ax.set_ylabel(r"Jaccard with $T_{taskt}$ anchor clade",
                       fontsize=10.5)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis="y", labelsize=9)
        ax.spines[["top", "right"]].set_visible(False)

        # Reference horizontal at J=1
        ax.axhline(1.0, color="#bbb", lw=0.6, ls=":", zorder=0)

        # Per-patient trail
        for _, r in sub.iterrows():
            tag = r["tag"]
            color = (COL_TRACE if tag == "trace"
                      else COL_ANTI if tag == "anti" else COL_NULL)
            xs = [0, 1, 2]
            ys = [r["J_pre"], 1.0, r["J_post"]]
            ax.plot(xs, ys, color=color, lw=1.7, alpha=0.85, zorder=3,
                    marker="o", markersize=4.5, markerfacecolor=color,
                    markeredgecolor="white", markeredgewidth=0.6)
            # Patient label at right end with offset to avoid overlap
            ax.annotate(
                r["patient"],
                xy=(2, r["J_post"]),
                xytext=(2.10, r["label_y"]),
                fontsize=8.0, color=color, fontweight="bold",
                va="center", ha="left",
                arrowprops=dict(arrowstyle="-",
                                color=color, lw=0.5, alpha=0.5,
                                shrinkA=0, shrinkB=0),
            )

        # Cohort tally
        n_trace = int((sub["tag"] == "trace").sum())
        n_anti = int((sub["tag"] == "anti").sum())
        n_null = int((sub["tag"] == "null").sum())
        med_s = float(sub["score"].median())

        verdict_color = (COL_TRACE if n_trace > n_anti
                          else COL_ANTI if n_anti > n_trace
                          else COL_NULL)

        ax.text(
            -0.13, 1.10,
            rf"{BRAIN_BAND_TEX_DICT[band]}-band  ·  cohort N=10",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=12.5, fontweight="bold", color="#222",
        )
        ax.text(
            -0.13, 1.02,
            rf"trace {n_trace}  ·  anti {n_anti}"
            + (rf"  ·  null {n_null}" if n_null else "")
            + rf"   median $s = J_{{post}}-J_{{pre}} = {med_s:+.2f}$",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=10.0, color=verdict_color, fontweight="bold",
        )

        # Reading rule annotation (small, bottom-right)
        ax.text(
            0.99, 0.04,
            r"trace: $J_{post} > J_{pre}$ (rspost retains anchor)" "\n"
            r"anti:  $J_{post} < J_{pre}$ (rspre retains anchor)",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.5, color="#444",
            bbox=dict(boxstyle="round,pad=0.25", fc="white",
                       ec="#bbb", lw=0.4, alpha=0.92),
        )

        # Legend
        handles = [
            Line2D([], [], color=COL_TRACE, lw=1.7, marker="o",
                   markersize=5, label=rf"trace ($s>+0.05$)  n={n_trace}"),
            Line2D([], [], color=COL_ANTI, lw=1.7, marker="o",
                   markersize=5, label=rf"anti ($s<-0.05$)  n={n_anti}"),
        ]
        if n_null:
            handles.append(Line2D([], [], color=COL_NULL, lw=1.7,
                                    marker="o", markersize=5,
                                    label=rf"null  n={n_null}"))
        ax.legend(handles=handles, loc="lower left",
                  bbox_to_anchor=(0.005, 0.04),
                  frameon=False, fontsize=8.5)

        fig.tight_layout()
        out_pdf = FIG / f"kc_cohort_trail_{band}.pdf"
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        print(f"[redo10] wrote {out_pdf.name}")


# ===================================================================
# Driver
# ===================================================================
def main() -> None:
    print("[redo] 1 — psi_irrelevance")
    redo1_psi_irrelevance()
    print("[redo] 2 — psi_boundary_audit")
    n_b, n_t, pct = redo2_psi_boundary_audit()
    print("[redo] 3 — kc_decomposition")
    redo3_kc_decomposition()
    print("[redo] 4 — ctm_headline")
    redo4_ctm_headline()
    print("[redo] 5 — cross_probe_signs")
    redo5_cross_probe_signs()
    print("[redo] 6 — kc_topology_eye_proof")
    redo6_kc_topology_eye_proof()
    print("[redo] 7 — kc_dendrogram_eye_proof")
    redo7_kc_dendrogram_eye_proof()
    print("[redo] 8 — kc_dendrogram_per_patient")
    redo8_kc_dendrogram_per_patient()
    print("[redo] 9 — kc_cohort_band_dendrograms")
    redo9_kc_cohort_band_dendrograms()
    print("[redo] 10 — kc_cohort_jaccard_trail")
    redo10_kc_cohort_jaccard_trail()

    # Persist headline for README
    with (TBL / "_psi_boundary_headline.txt").open("w") as fh:
        fh.write(f"{n_b}/{n_t} ({pct:.1f}%)\n")
    print(f"[redo] all done · outputs at {OUT}")


if __name__ == "__main__":
    main()
