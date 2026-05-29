#!/usr/bin/env python3
"""Audit round-3 Section-5 redos — psi-irrelevance + Psi-boundary audit.

Split out of audit_round3_section5_redo.py on 2026-05-29 (Phase 4-B
split 2/7). Imports shared constants + helpers from
_audit_round3_shared and runs: redo1, redo2.
"""
from __future__ import annotations

import matplotlib  # noqa: F401
import matplotlib.lines as mlines  # noqa: F401
import matplotlib.patches as mpatches  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram  # noqa: F401

from _audit_round3_shared import (
    BAND_ORDER, BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, DEND_LABEL_DROP_HEIGHTS,
    DEND_LABEL_FONTSIZE, DEND_LABEL_Y_FRAC, FIG, OUT, PATIENTS_LIST,
    PHASES, R2, R2_TBL, ROOT, S5, TBL, asterisks, kc_vectors,
    load_lrg_result,
)


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


def main() -> None:
    print("[redo] 1 — psi_irrelevance")
    redo1_psi_irrelevance()
    print("[redo] 2 — psi_boundary_audit")
    n_b, n_t, pct = redo2_psi_boundary_audit()
    with (TBL / "_psi_boundary_headline.txt").open("w") as fh:
        fh.write(f"{n_b}/{n_t} ({pct:.1f}%)\n")
    print(f"[psi] done · outputs at {OUT}")


if __name__ == "__main__":
    main()
