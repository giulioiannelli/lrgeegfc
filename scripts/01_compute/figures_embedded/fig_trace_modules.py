#!/usr/bin/env python3
"""Trace-modules figures — cohort summary + per-patient mini-dendrograms.

Three visual outputs grounding the headline (band × k) cells in the actual
(J_pre, J_post) candidate cloud and best-match leafsets:

  1. data/outputs/figures/section6/trace_modules_summary_n10_imcoh_abs.pdf
     - Top: (J_pre, J_post) candidate scatter per (band, k_bin) cell, with
       J_min thresholds rendered. Visually shows whether the cohort sits on
       the identity diagonal (no trace) or in the T quadrant
       (¬present_pre ∧ present_post).
     - Middle: J_min sweep — bar chart of T-module counts per (band, k_bin)
       at J_min ∈ {0.70, 0.85, 0.90, 0.95}.
     - Bottom: per-patient T-module count grid (band × patient) at primary
       J_min = 0.90 — surfaces single-patient dominance.

  2. data/outputs/figures/section6/trace_modules_candidates_n10_imcoh_abs.pdf
     - Per (band, k_bin) cell: 4-panel grid (one panel per J_min) showing the
       same candidate cloud with cell-by-cell threshold-quadrant colouring.

  3. data/audit/trace_modules/Pat_NN_band_<band>_k<k>.pdf — per-patient
     mini-dendrograms (3 stacked: rest_pre, task_test, rest_post) for every
     patient that contributes ≥ 1 T module at primary J_min in the cell. The
     T module's leaves are coloured consistently across phases. Cap by the
     CSV — currently ≤ 5 PDFs since strict-J trace is near-null.

Reads:
  data/audit/trace_modules/trace_subtrees_n10_imcoh_abs{,_jmin{70,85,95}}.csv
  data/audit/trace_modules/trace_subtrees_n10_imcoh_abs_candidates.csv
  Pat_NN LRG imcoh_abs caches (for mini-dendrograms only).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PRIMARY_J = 0.90
J_SWEEP = [0.70, 0.85, 0.90, 0.95]
PHASE_ORDER = ["rest_pre", "task_test", "rest_post"]
T_COLOURS = ["#d7191c", "#2c7bb6", "#fdae61", "#1a9641", "#6a51a3"]


def _audit_dir() -> Path:
    p = Path("data/audit/trace_modules")
    if not p.is_absolute():
        p = ROOT / p
    return p


def _load_jmin_csv(j_min: float) -> pd.DataFrame:
    base = _audit_dir() / "trace_subtrees_n10_imcoh_abs"
    if abs(j_min - PRIMARY_J) < 1e-6:
        path = base.with_suffix(".csv")
    else:
        path = base.parent / (base.name + f"_jmin{int(j_min*100):02d}.csv")
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _load_candidates() -> pd.DataFrame:
    p = _audit_dir() / "trace_subtrees_n10_imcoh_abs_candidates.csv"
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def _cell_label(band: str, k_lo: int, k_hi: int) -> str:
    return f"{BRAIN_BAND_TEX_DICT.get(band, band)}, k={k_lo}–{k_hi}"


def _summary_figure(cands: pd.DataFrame, dfs_by_jmin: dict[float, pd.DataFrame]) -> None:
    """Cohort summary PDF — three rows of panels."""
    cells = (cands[["band", "k_lo", "k_hi"]]
             .drop_duplicates()
             .sort_values(["band", "k_lo"])
             .to_records(index=False))
    n_cells = len(cells)
    if n_cells == 0:
        print("no cells in candidates — skipping summary figure")
        return

    fig = plt.figure(figsize=(6 + 4 * n_cells, 12.5), dpi=160)
    gs = fig.add_gridspec(3, max(n_cells, 1), height_ratios=[1.4, 1.0, 1.0],
                          hspace=0.42, wspace=0.30)

    # ── Row 1: candidate (J_pre, J_post) cloud per cell ───────────────
    band_palette = {b: c for b, c in zip(
        sorted(cands["band"].unique()),
        plt.cm.tab10.colors,
    )}
    for ic, row in enumerate(cells):
        ax = fig.add_subplot(gs[0, ic])
        sub = cands[(cands["band"] == row.band)
                    & (cands["k_lo"] == row.k_lo)
                    & (cands["k_hi"] == row.k_hi)]
        for pat, sp in sub.groupby("patient"):
            ax.scatter(sp["J_pre"], sp["J_post"],
                       s=22, edgecolor="white", linewidth=0.4,
                       color=band_palette.get(row.band, "0.4"),
                       alpha=0.65, label=pat if ic == 0 else None)
        for jm in J_SWEEP:
            ax.axhline(jm, color="0.6", linewidth=0.5, linestyle=":")
            ax.axvline(jm, color="0.6", linewidth=0.5, linestyle=":")
        # Highlight the T quadrant at primary J_min
        ax.add_patch(plt.Rectangle((0, PRIMARY_J), PRIMARY_J, 1 - PRIMARY_J,
                                    facecolor=(0.95, 0.6, 0.2), alpha=0.10,
                                    edgecolor="#d7191c", linewidth=1.2,
                                    linestyle="--", zorder=0))
        ax.text(PRIMARY_J / 2, (1 + PRIMARY_J) / 2, "T",
                ha="center", va="center", fontsize=22, color="#d7191c",
                alpha=0.30, weight="bold", zorder=0)
        ax.plot([0, 1], [0, 1], color="0.5", linewidth=0.6, alpha=0.6)
        ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02)
        ax.set_xlabel(r"$J^*$ vs rest_pre", fontsize=10)
        if ic == 0:
            ax.set_ylabel(r"$J^*$ vs rest_post", fontsize=10)
        ax.set_title(_cell_label(row.band, row.k_lo, row.k_hi),
                     fontsize=11, loc="left")
        ax.grid(alpha=0.2, linestyle=":")

    # ── Row 2: J_min sweep T-module counts ──────────────────────────
    ax_sweep = fig.add_subplot(gs[1, :])
    bars_x = np.arange(n_cells)
    width = 0.22
    for ij, jm in enumerate(J_SWEEP):
        df = dfs_by_jmin.get(jm, pd.DataFrame())
        counts = []
        for row in cells:
            if df.empty:
                counts.append(0)
                continue
            n = ((df["band"] == row.band)
                 & (df["k_lo"] == row.k_lo)
                 & (df["k_hi"] == row.k_hi)).sum()
            counts.append(int(n))
        offset = (ij - (len(J_SWEEP) - 1) / 2) * width
        ax_sweep.bar(bars_x + offset, counts, width=width,
                     label=f"$J_{{\\min}} = {jm}$",
                     edgecolor="black" if abs(jm - PRIMARY_J) < 1e-6 else "none",
                     linewidth=1.5 if abs(jm - PRIMARY_J) < 1e-6 else 0)
    ax_sweep.set_xticks(bars_x)
    ax_sweep.set_xticklabels([_cell_label(r.band, r.k_lo, r.k_hi) for r in cells],
                              fontsize=10)
    ax_sweep.set_ylabel("T-module count\n(pooled across n=10)", fontsize=10)
    ax_sweep.legend(ncol=4, fontsize=9, frameon=False, loc="upper right")
    ax_sweep.grid(axis="y", alpha=0.25, linestyle=":")
    ax_sweep.set_title(
        r"$J_{\min}$ sweep — strict subtree-identity T-regime is near-null at "
        r"$J_{\min} \geq 0.85$ (matches the MRL $\leftrightarrow$ CBR reconciliation)",
        fontsize=11, loc="left",
    )

    # ── Row 3: per-patient T count grid at primary J ────────────────
    df_primary = dfs_by_jmin.get(PRIMARY_J, pd.DataFrame())
    ax_grid = fig.add_subplot(gs[2, :])
    grid = np.zeros((n_cells, len(COHORT_N10)), dtype=int)
    for ic, row in enumerate(cells):
        for ip, pat in enumerate(COHORT_N10):
            if df_primary.empty:
                continue
            grid[ic, ip] = ((df_primary["band"] == row.band)
                            & (df_primary["k_lo"] == row.k_lo)
                            & (df_primary["k_hi"] == row.k_hi)
                            & (df_primary["patient"] == pat)).sum()
    im = ax_grid.imshow(grid, aspect="auto", cmap="OrRd",
                         vmin=0, vmax=max(1, int(grid.max())),
                         interpolation="nearest")
    ax_grid.set_xticks(range(len(COHORT_N10)))
    ax_grid.set_xticklabels(COHORT_N10, fontsize=9, rotation=45, ha="right")
    ax_grid.set_yticks(range(n_cells))
    ax_grid.set_yticklabels([_cell_label(r.band, r.k_lo, r.k_hi) for r in cells],
                              fontsize=10)
    for ic in range(n_cells):
        for ip in range(len(COHORT_N10)):
            v = grid[ic, ip]
            if v > 0:
                ax_grid.text(ip, ic, str(v), ha="center", va="center",
                              fontsize=10, color="white" if v >= 2 else "black")
    cbar = fig.colorbar(im, ax=ax_grid, shrink=0.7, pad=0.02)
    cbar.set_label(f"T-modules at $J_{{\\min}} = {PRIMARY_J}$", fontsize=9)
    ax_grid.set_title(
        f"Per-patient T-module count grid at primary $J_{{\\min}} = {PRIMARY_J}$",
        fontsize=11, loc="left",
    )

    out = FIGURES_ROOT / "section6" / "trace_modules_summary_n10_imcoh_abs.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def _mini_dendro_for_T(df_primary: pd.DataFrame) -> None:
    """Per-patient 3-stacked dendrogram PDFs for every (patient, band, k_bin)
    cell with at least one T-module at primary J_min.
    """
    if df_primary.empty:
        print("no T-modules at primary J_min — skipping mini-dendrograms")
        return
    audit_dir = _audit_dir()
    audit_dir.mkdir(parents=True, exist_ok=True)
    for (pat, band, k_lo, k_hi), grp in df_primary.groupby(
            ["patient", "band", "k_lo", "k_hi"]):
        Z = {}
        for phase in PHASE_ORDER:
            try:
                r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
                Z[phase] = np.asarray(r.linkage_matrix)
            except Exception:
                print(f"  skip {pat} {phase} {band}")
                Z[phase] = None
        if any(z is None for z in Z.values()):
            continue
        # Build a colour map per T-module
        modules = []  # list[set[int]]
        for _, row in grp.iterrows():
            leaves = set(int(x) for x in str(row["leaf_indices"]).split(";") if x)
            modules.append(leaves)
        leaf_to_colour: dict[int, str] = {}
        for im, mod in enumerate(modules):
            colour = T_COLOURS[im % len(T_COLOURS)]
            for ℓ in mod:
                leaf_to_colour[ℓ] = colour

        fig, axes = plt.subplots(3, 1, figsize=(11.0, 9.0), dpi=160,
                                  sharex=False)
        for ax, phase in zip(axes, PHASE_ORDER):
            link_z = Z[phase]
            heights = link_z[:, 2]
            tmin = float(heights[0]) * 0.8
            tmax = float(heights[-1]) * 1.05
            d = dendrogram(
                link_z,
                ax=ax,
                no_labels=True,
                color_threshold=0,  # all branches grey
                above_threshold_color="0.55",
                leaf_font_size=6,
            )
            # Recolor coloured leaves
            leaves_order = d["leaves"]
            xs = np.arange(len(leaves_order)) * 10 + 5  # scipy default leaf spacing
            for x, ℓ in zip(xs, leaves_order):
                colour = leaf_to_colour.get(int(ℓ), None)
                if colour is None:
                    continue
                ax.scatter([x], [tmin * 1.02], color=colour, s=18, marker="|",
                           clip_on=False, linewidths=2.5, zorder=5)
            ax.set_yscale("log")
            ax.set_ylim(tmin, tmax)
            ax.set_title(
                f"{pat} | {BRAIN_BAND_TEX_DICT.get(band, band)} | "
                f"k={k_lo}–{k_hi} | {phase} | "
                f"{len(modules)} T-module(s)",
                fontsize=10, loc="left",
            )
            ax.set_xlabel("")  # leaves are anonymous on this view
        out = audit_dir / f"{pat}_band_{band}_k{k_lo}-{k_hi}.pdf"
        fig.tight_layout()
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        print(f"saved {out}")


def main() -> None:
    cands = _load_candidates()
    if cands.empty:
        raise SystemExit("trace_subtrees_*_candidates.csv missing — run audit_15 first")
    dfs_by_jmin = {jm: _load_jmin_csv(jm) for jm in J_SWEEP}
    _summary_figure(cands, dfs_by_jmin)
    _mini_dendro_for_T(dfs_by_jmin[PRIMARY_J])


if __name__ == "__main__":
    main()
