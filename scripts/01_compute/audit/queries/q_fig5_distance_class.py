#!/usr/bin/env python3
"""q_bundle figure fig5 -- distance-class illustration.

Split out of q_bundle_figures.py on 2026-05-29 (Phase 4-B split 3/7).
"""
from __future__ import annotations

import sys
from pathlib import Path

# Folder-local imports
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # noqa: F401
from matplotlib.patches import Ellipse, Patch, Rectangle, Arc  # noqa: F401
from matplotlib.lines import Line2D  # noqa: F401
from matplotlib.colors import Normalize  # noqa: F401
from scipy.stats import pearsonr, spearmanr  # noqa: F401

from _q_bundle_shared import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, NODE_SIZE_DEFAULT, N_COHORT,
    PAIR_COLORS, PAIR_LABELS, PAIR_ORDER, PATIENTS_4PHASE,
    PHASE_SHORT, ROOT,
    compute_network_layout, draw_network_edges, _probe_color_map,
    extract_probe_labels, imshow_colorbar_caxdivider, load_bundle_data,
    load_channel_labels, load_fc_matrix, lookup_pair, resolve_substrate,
)


def main() -> None:
    SUBSTRATE, SUFFIX, SRC, OUT = resolve_substrate()
    FC_METHOD_FIG5 = SUBSTRATE
    print(f"=== q_fig5 for substrate: {SUBSTRATE} ===")
    distances_long, td, contrast = load_bundle_data(SRC)



    # ----------------------------------------------------------------------
    # Fig 5 — Distance-class illustration (real-pair examples)
    # ----------------------------------------------------------------------
    def _common_giant_indices(adj_phases):
        """Common giant-component node indices across phases (per-pair V*)."""
        sets = []
        for adj in adj_phases.values():
            G = nx.from_numpy_array(np.abs(adj))
            comps = list(nx.connected_components(G))
            sets.append(set(max(comps, key=len)) if comps else set())
        common = set.intersection(*sets) if sets else set()
        return np.array(sorted(common), dtype=int)


    from matplotlib.colors import LogNorm, LinearSegmentedColormap
    from scipy.optimize import minimize as _spo_minimize
    from matplotlib.ticker import (
        LogLocator,
        LogFormatterMathtext,
        NullLocator,
        NullFormatter,
    )

    # Probe-sort + outline helpers from the same archive module
    from _shared import (  # noqa: E402
        probe_sort_indices,
        probe_boundaries,
        draw_probe_outlines,
    )

    EXAMPLES = [
        {
            "patient": "Pat_14",
            "band": "low_gamma",
            "phase_A": "task_learn",
            "phase_B": "rest_post",
            "label": r"high $d_S$, low $d_P$  —  topological reorganization",
        },
        {
            "patient": "Pat_02",
            "band": "high_gamma",
            "phase_A": "task_test",
            "phase_B": "rest_post",
            "label": r"low $d_S$, high $d_P$  —  magnitude redistribution",
        },
    ]

    # Three super-columns so the heatmap pair + network pair can each be
    # packed tight (wspace ~ 0.05 inside the pair) while the scatter
    # in the middle has its own breathing room from the colorbar at the
    # right edge of the heatmap pair.
    fig = plt.figure(figsize=(14.6, 6.6))
    gs_outer = fig.add_gridspec(
        2, 3,
        width_ratios=[2.0, 1.05, 2.0],
        height_ratios=[1, 1],
        wspace=0.10,
        hspace=0.55,
        top=0.90, bottom=0.13, left=0.04, right=0.985,
    )

    # axes_by_row[r] = [heatA, heatB, scatter, netA, netB]
    axes_by_row = [[None] * 5 for _ in EXAMPLES]

    for r, ex in enumerate(EXAMPLES):
        gs_heat = gs_outer[r, 0].subgridspec(1, 2, wspace=0.06)
        ax_a = fig.add_subplot(gs_heat[0, 0])
        ax_b = fig.add_subplot(gs_heat[0, 1])
        ax_sc = fig.add_subplot(gs_outer[r, 1])
        gs_net = gs_outer[r, 2].subgridspec(1, 2, wspace=0.04)
        ax_n1 = fig.add_subplot(gs_net[0, 0])
        ax_n2 = fig.add_subplot(gs_net[0, 1])
        axes_by_row[r] = [ax_a, ax_b, ax_sc, ax_n1, ax_n2]

    for r, ex in enumerate(EXAMPLES):
        ax_a, ax_b, ax_sc, ax_n1, ax_n2 = axes_by_row[r]
        p, b = ex["patient"], ex["band"]
        phA, phB = ex["phase_A"], ex["phase_B"]

        # Load FC + restrict to per-pair V*
        adj = {ph: np.asarray(load_fc_matrix(p, ph, b, FC_METHOD_FIG5),
                              dtype=np.float64)
               for ph in (phA, phB)}
        common_idx = _common_giant_indices(adj)
        A_a = adj[phA][np.ix_(common_idx, common_idx)]
        A_b = adj[phB][np.ix_(common_idx, common_idx)]
        np.fill_diagonal(A_a, 0.0)
        np.fill_diagonal(A_b, 0.0)

        # Channel labels + probe info (restricted to V*)
        ch_full = load_channel_labels(p)
        ch = [ch_full[i] for i in common_idx]
        probes = extract_probe_labels(ch)
        # Sort everything by sEEG probe → block structure on the heatmaps
        sort_idx = probe_sort_indices(ch)
        A_a = A_a[np.ix_(sort_idx, sort_idx)]
        A_b = A_b[np.ix_(sort_idx, sort_idx)]
        ch = [ch[i] for i in sort_idx]
        probes = [probes[i] for i in sort_idx]
        probe_cmap = _probe_color_map(probes)
        node_colors = [probe_cmap[pl] for pl in probes]

        # Probe-block boundaries → tick positions at probe block centers
        bnd = probe_boundaries(probes)
        edges = [0] + list(bnd) + [len(probes)]
        block_centers = [(edges[k] + edges[k + 1]) / 2 - 0.5
                         for k in range(len(edges) - 1)]
        block_names = [probes[edges[k]] for k in range(len(edges) - 1)]

        # Distance triplet
        ut_a = A_a[np.triu_indices_from(A_a, k=1)]
        ut_b = A_b[np.triu_indices_from(A_b, k=1)]
        d_P = float(1.0 - pearsonr(ut_a, ut_b).statistic)
        d_S = float(1.0 - spearmanr(ut_a, ut_b).statistic)
        d_F = float(np.linalg.norm(A_a - A_b, ord="fro") /
                    np.sqrt(np.linalg.norm(A_a, ord="fro")
                            * np.linalg.norm(A_b, ord="fro")))

        # Log-scale shared range across both heatmaps and the scatter.
        # vmin/vmax come straight from the data (NOT snapped to integer
        # decades). Tick labels are placed only at integer powers of 10
        # via LogLocator(subs=(1.0,)) + LogFormatterMathtext(labelOnlyBase),
        # so a colorbar whose data tops at 0.15 only shows 10**-2 / 10**-1,
        # never 10**0. Minor ticks at 2..9 × 10**k are drawn as tick marks
        # for scale density without their own labels.
        pos_vals = np.concatenate([ut_a[ut_a > 0], ut_b[ut_b > 0]])
        lo = float(np.percentile(pos_vals, 1))
        hi = float(np.percentile(pos_vals, 99.5))
        norm = LogNorm(vmin=lo, vmax=hi)

        def _apply_log_ticks(axis, max_majors=4):
            """Major ticks at 10**k labelled; minor ticks at 2..9 × 10**k
            without labels (just tick marks)."""
            axis.set_major_locator(
                LogLocator(base=10, subs=(1.0,), numticks=max_majors)
            )
            axis.set_major_formatter(
                LogFormatterMathtext(base=10, labelOnlyBase=True)
            )
            axis.set_minor_locator(
                LogLocator(base=10, subs=np.arange(2, 10), numticks=80)
            )
            axis.set_minor_formatter(NullFormatter())

        def _apply_probe_ticks(ax, axis="both"):
            ax.set_xticks(block_centers)
            ax.set_xticklabels(block_names, fontsize=5.0, rotation=90,
                               ha="center")
            ax.set_yticks(block_centers)
            ax.set_yticklabels(block_names, fontsize=5.0)
            ax.tick_params(axis="both", which="both",
                           length=2, pad=1.5, color="#888")

        # ---- Heatmap A (log scale, probe-sorted, probe-tick labels) ----
        ax = ax_a
        im_a = ax.imshow(np.where(A_a > 0, A_a, np.nan), cmap="magma",
                         norm=norm, aspect="equal", interpolation="nearest")
        _apply_probe_ticks(ax)
        draw_probe_outlines(ax, ch, color="#ffffff", lw=0.4, alpha=0.4)
        ax.set_title(rf"$A^{{\mathrm{{{PHASE_SHORT[phA]}}}}}$", fontsize=10,
                     pad=4)

        # ---- Heatmap B + per-row colorbar (log scale, same sort/blocks) -
        ax = ax_b
        im_b = ax.imshow(np.where(A_b > 0, A_b, np.nan), cmap="magma",
                         norm=norm, aspect="equal", interpolation="nearest")
        _apply_probe_ticks(ax)
        draw_probe_outlines(ax, ch, color="#ffffff", lw=0.4, alpha=0.4)
        # Don't repeat the y-tick labels on the right heatmap
        ax.set_yticklabels([""] * len(block_centers))
        ax.set_title(rf"$A^{{\mathrm{{{PHASE_SHORT[phB]}}}}}$", fontsize=10,
                     pad=4)
        _div, _cax, cbar = imshow_colorbar_caxdivider(im_b, ax,
                                                      size="5%", pad=0.04)
        _apply_log_ticks(cbar.ax.yaxis, max_majors=4)
        cbar.ax.tick_params(which="major", labelsize=6.5, pad=1, length=3)
        cbar.ax.tick_params(which="minor", length=1.5)
        cbar.set_label(r"$|\mathrm{ImCoh}|$  (log)", fontsize=7.5,
                       labelpad=2)

        # ---- Scatter: log-log density + identity + power-law fit --------
        ax = ax_sc
        cmap_sc = plt.get_cmap("cividis")
        ax.set_facecolor(cmap_sc(0.0))     # empty cells blend with cmap
        # Limits match the actual data extent (no spurious upper white band)
        s_lo = float(min(ut_a.min(), ut_b.min()))
        s_hi = float(max(ut_a.max(), ut_b.max()))
        s_lo = max(s_lo, 1e-7)
        bins = np.geomspace(s_lo, s_hi, 60)
        H, xedges, yedges = np.histogram2d(
            np.clip(ut_a, s_lo, s_hi), np.clip(ut_b, s_lo, s_hi),
            bins=[bins, bins],
        )
        H_log = np.log10(H + 1)
        ax.pcolormesh(xedges, yedges, H_log.T, cmap=cmap_sc,
                      shading="auto", zorder=1)
        ax.plot([s_lo, s_hi], [s_lo, s_hi], color="white", lw=1.2, ls="--",
                alpha=0.9, zorder=3, label="identity")
        lx = np.log10(np.clip(ut_a, s_lo, s_hi))
        ly = np.log10(np.clip(ut_b, s_lo, s_hi))
        slope, intercept = np.polyfit(lx, ly, 1)
        xline = np.geomspace(s_lo, s_hi, 100)
        yline = 10 ** intercept * xline ** slope
        ax.plot(xline, yline, color="#ff5252", lw=1.4,
                label=rf"power fit (slope={slope:.2f})", zorder=4)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(s_lo, s_hi); ax.set_ylim(s_lo, s_hi)
        ax.set_aspect("equal")
        # Same log-tick policy as the colorbar: 10**k labelled, 2..9 × 10**k
        # as unlabelled minor tick marks for scale density.
        _apply_log_ticks(ax.xaxis, max_majors=4)
        _apply_log_ticks(ax.yaxis, max_majors=4)
        ax.set_xlabel(rf"$\mathrm{{triu}}(A^{{\mathrm{{{PHASE_SHORT[phA]}}}}})$",
                      fontsize=8.5, labelpad=2)
        ax.set_ylabel(rf"$\mathrm{{triu}}(A^{{\mathrm{{{PHASE_SHORT[phB]}}}}})$",
                      fontsize=8.5, labelpad=2)
        ax.tick_params(labelsize=7, pad=1.5)
        ax.legend(fontsize=6.5, loc="upper left", frameon=True,
                  facecolor="white", edgecolor="#bbb", framealpha=0.92,
                  borderpad=0.3, handlelength=1.6, handletextpad=0.4)
        ax.text(0.97, 0.05,
                (rf"$d_S$ = {d_S:.2f}" "\n"
                 rf"$d_P$ = {d_P:.2f}" "\n"
                 rf"$d_F$ = {d_F:.2f}"),
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="#999", alpha=0.92))

        # ---- LRG-MDS-layout networks (probe-coloured) -------------------
        # Pass the probe-sorted A so positions are in sorted indexing.
        pos = compute_network_layout(
            0.5 * (A_a + A_b),
            fc_method=FC_METHOD_FIG5,
            patient=p, phase=phA, band=b,
        )
        for net_ax, mat, ph in [(ax_n1, A_a, phA), (ax_n2, A_b, phB)]:
            # Canonical defaults (rank-based scaling, per-fc_method
            # gamma/wmin/wmax/amin/amax from config.const.EDGE_RANK_*).
            # NEVER pass scaling="value" with a heatmap-derived value_range
            # — the heavy-tailed FC distribution then gives a few huge
            # edges and mostly invisible ones.
            draw_network_edges(
                net_ax, pos, mat, fc_method=FC_METHOD_FIG5,
                probe_labels=probes,
                highlight_same_probe=True,
            )
            net_ax.scatter(pos[:, 0], pos[:, 1], c=node_colors,
                           s=NODE_SIZE_DEFAULT * 0.5,
                           edgecolors="white", linewidths=0.5, zorder=5)
            net_ax.set_xticks([]); net_ax.set_yticks([])
            net_ax.set_aspect("equal")
            net_ax.margins(0.06)
            net_ax.set_title(
                rf"network  $A^{{\mathrm{{{PHASE_SHORT[ph]}}}}}$",
                fontsize=10, pad=4,
            )
            for spine in net_ax.spines.values():
                spine.set_visible(False)

    # Per-row title placed in the hspace gap, computed from axis positions
    for r, ex in enumerate(EXAMPLES):
        p, b = ex["patient"], ex["band"]
        phA, phB = ex["phase_A"], ex["phase_B"]
        pos_left = axes_by_row[r][0].get_position()
        pos_right = axes_by_row[r][-1].get_position()
        y = pos_left.y1 + 0.045
        x_centre = 0.5 * (pos_left.x0 + pos_right.x1)
        fig.text(
            x_centre, y,
            f"{p} · {BRAIN_BAND_TEX_DICT[b]} · "
            f"{PHASE_SHORT[phA]} $\\leftrightarrow$ {PHASE_SHORT[phB]}    "
            f"—    {ex['label']}",
            ha="center", va="bottom", fontsize=10.5, fontweight="bold",
        )

    # Footer-bar legend explaining the network + scatter conventions
    handles_f5 = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="#bbb", markeredgecolor="white",
               markersize=8, label="node = electrode (colour = sEEG probe)"),
        Line2D([0], [0], color="#666", lw=2.0,
               label="cross-probe edge (gray, weight-scaled)"),
        Line2D([0], [0], color="#1f77b4", lw=2.4,
               label="same-probe edge (probe colour, full alpha)"),
        Line2D([0], [0], color="white", lw=1.4, ls="--",
               label="identity (scatter)"),
        Line2D([0], [0], color="#ff5252", lw=1.4,
               label="power-law fit (scatter)"),
    ]
    fig.legend(handles=handles_f5, loc="lower center",
               bbox_to_anchor=(0.5, 0.015), ncol=len(handles_f5),
               frameon=False, fontsize=7.8)

    fig.savefig(OUT / f"fig5_distance_class_examples{SUFFIX}.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'fig5_distance_class_examples{SUFFIX}.pdf'}")


if __name__ == "__main__":
    main()
