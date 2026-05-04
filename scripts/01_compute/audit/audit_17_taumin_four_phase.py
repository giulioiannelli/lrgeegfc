#!/usr/bin/env python3
"""Audit Step 17 — τ=1/λ_max ultrametric dendrograms across the 4 phases.

For every (patient, band) in the n=10 cohort × 6 bands grid, render two
PDFs from the cached imcoh_abs LRG analysis (which uses
``compute_laplacian_properties(G, tau=None)`` → τ defaults to 1/λ_max):

  * <band>_dendrograms.pdf — 4-axis figure (one per phase) with the
                              Ψ-optimal cut drawn at the matching height.
  * <band>_psi.pdf          — Ψ(n; τ=1/λ_max) curves overlaid for the
                              4 phases plus a vertical mark at n*.

A CSV ledger ``ledger.csv`` records (patient, band, phase, n_nodes,
psi_n*, psi_max, optimal_threshold).

Outputs
-------
data/audit/four_phase_taumin/<Patient>/<band>_dendrograms.pdf
data/audit/four_phase_taumin/<Patient>/<band>_psi.pdf
data/audit/four_phase_taumin/ledger.csv
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, to_tree

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BAND_TEX_DICT,
    BRAIN_BANDS_NAMES,
    PATIENTS_4PHASE,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import DATA_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.lrg import (
    _load_channel_labels,
    _psi_threshold_from_linkage,
    compute_partition_stability_index,
    plot_lrg_dendrogram_shaft_colored,
)
from lrg_eegfc.workflow.lrg import load_lrg_result


FC_METHOD = "imcoh_abs"
OUT_DIR = DATA_ROOT / "audit" / "four_phase_taumin"
LEDGER_PATH = OUT_DIR / "ledger.csv"
#: How many top-Ψ cut heights to overlay on each dendrogram (the largest is
#: drawn as a thicker dashed line, the rest as dotted secondary cuts).
TOP_M_CUTS: int = 5
#: Reference phase whose K-cluster cut defines the leaf colour palette
#: that is projected onto every phase's dendrogram.
REF_PHASE: str = "rest_post"
#: K values to scan for the reference cross-phase figure.
K_REF_LIST: tuple = (10, 20, 40)

PHASE_TEX = {
    "rest_pre":   r"rest$_\mathrm{pre}$",
    "task_learn": r"task$_\mathrm{learn}$",
    "task_test":  r"task$_\mathrm{test}$",
    "rest_post":  r"rest$_\mathrm{post}$",
}
PHASE_COLOR = {
    "rest_pre":   "#4477aa",
    "task_learn": "#ccbb44",
    "task_test":  "#ee6677",
    "rest_post":  "#228833",
}


# ─── plotting ──────────────────────────────────────────────────────────────


def _shared_ylim(results: Dict[str, "object"]) -> Optional[tuple]:
    heights = []
    for r in results.values():
        if r is None:
            continue
        h = r.linkage_matrix[:, 2]
        h = h[h > 0]
        if h.size:
            heights.append(h)
    if not heights:
        return None
    merged = np.concatenate(heights)
    # No padding: bottom == global smallest merge, top == global root merge.
    # Anything looser stretches the leaves below their actual height and
    # leaves visible empty space above the root.
    return float(merged.min()), float(merged.max())


def _top_m_cut_thresholds(Z: np.ndarray, m: int) -> List[Tuple[int, float]]:
    """Top *m* Ψ-cuts ranked by log10 gap, as ``[(n_communities, height), ...]``
    in descending Ψ order (largest first)."""
    psi_v, psi_nc = compute_partition_stability_index(Z)
    if psi_v.size == 0:
        return []
    order = np.argsort(psi_v)[::-1][:m]
    n_nodes = Z.shape[0] + 1
    fallback = float(np.median(Z[:, 2]))
    out: List[Tuple[int, float]] = []
    for idx in order:
        n_target = int(psi_nc[int(idx)])
        thresh = _psi_threshold_from_linkage(Z, n_target, n_nodes, fallback)
        out.append((n_target, float(thresh)))
    return out


def _plot_dendrograms(
    pat: str, band: str, results: Dict[str, "object"], labels: List[str],
    out_dir: Path = OUT_DIR,
) -> Path:
    out_path = out_dir / pat / f"{band}_dendrograms.pdf"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ylim = _shared_ylim(results)

    fig, axes = plt.subplots(
        1, len(PHASE_LABELS), figsize=(4.5 * len(PHASE_LABELS), 4.5),
        sharey=False,  # see note below — sharey=True breaks the log-axis
                       # rendering: scipy's `dendrogram` internally calls
                       # ax.set_ylim([0, dvw]); when sharey is True the first
                       # axis is already log by the time the 2nd–4th panels
                       # are drawn, which raises the non-positive-ylim warning.
    )
    axes = np.atleast_1d(axes)

    for ax_idx, (ax, phase) in enumerate(zip(axes, PHASE_LABELS)):
        r = results.get(phase)
        if r is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes, color="#999")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(PHASE_TEX[phase], fontsize=11)
            continue
        # Bypass the library helper's "prefer larger-n peak" heuristic and
        # force the cut at the *raw* argmax of Ψ (= biggest log10 gap between
        # consecutive merge heights). On log y-axis this is the visibly
        # biggest empty band in the dendrogram.
        psi_v, psi_nc = compute_partition_stability_index(r.linkage_matrix)
        n_target = (
            int(psi_nc[int(np.argmax(psi_v))]) if psi_v.size else 2
        )
        info = plot_lrg_dendrogram_shaft_colored(
            ax,
            r,
            labels,
            show_xlabels=(r.n_nodes <= 80),
            leaf_font_size=3.5,
            ylim=ylim,
            min_n_communities=n_target,
            max_n_communities=n_target,
            show_psi_cut=False,  # we draw a thicker solid cut below
            above_threshold_color="0.25",
            title=f"{PHASE_TEX[phase]}  ·  n*={'?'}",  # patched below
        )
        # Top-M Ψ cuts (ranked by log10 gap). The largest is drawn as a
        # thicker dashed line at full opacity; the secondary cuts as dotted
        # transparent lines so they read as a "where would alternative cuts
        # land" overlay. All drawn behind the tree (zorder=0).
        top_cuts = _top_m_cut_thresholds(r.linkage_matrix, TOP_M_CUTS)
        for rank, (_, h) in enumerate(top_cuts):
            if rank == 0:
                ax.axhline(h, color="black", ls="--", lw=1.4, alpha=1.0,
                           zorder=0)
            else:
                ax.axhline(h, color="black", ls=":", lw=0.6, alpha=0.9,
                           zorder=0)
        ax.set_title(
            f"{PHASE_TEX[phase]}  ·  n*={info['psi_n']}  ·  "
            rf"cut={info['psi_threshold']:.3g}",
            fontsize=10,
        )
        if ax_idx > 0:
            plt.setp(ax.get_yticklabels(), visible=False)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _build_distinct_palette(k: int) -> Dict[int, str]:
    """Return ``{cluster_id (1..k): hex_color}`` with *k* visually distinct
    colours. Uses tab10 / tab20 for small k, then concatenates tab20 +
    tab20b + tab20c (60 curated colours), then falls back to evenly-spaced
    HSV samples for very large k.
    """
    if k <= 0:
        return {}
    if k <= 10:
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i) for i in range(k)]
    elif k <= 20:
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i) for i in range(k)]
    elif k <= 60:
        bases = [plt.get_cmap(name) for name in ("tab20", "tab20b", "tab20c")]
        pool = [c(i) for c in bases for i in range(20)]
        colors = pool[:k]
    else:
        cmap = plt.get_cmap("hsv")
        colors = [cmap(i / k) for i in range(k)]
    return {
        i + 1: plt.matplotlib.colors.to_hex(c)
        for i, c in enumerate(colors)
    }


def _first_merge_heights_per_leaf(Z: np.ndarray) -> np.ndarray:
    """For each leaf id (0..N-1), the height of the first merge that consumes
    that leaf — i.e. the top of its singleton stem in the dendrogram."""
    n_leaves = Z.shape[0] + 1
    heights = np.full(n_leaves, np.nan, dtype=float)
    for row in range(Z.shape[0]):
        c0 = int(Z[row, 0])
        c1 = int(Z[row, 1])
        h = float(Z[row, 2])
        if c0 < n_leaves and np.isnan(heights[c0]):
            heights[c0] = h
        if c1 < n_leaves and np.isnan(heights[c1]):
            heights[c1] = h
    return heights


def _link_color_uniform(
    Z: np.ndarray, leaf_cl: np.ndarray, palette: Dict[int, tuple], mixed: str,
) -> Dict[int, tuple]:
    """For each internal node id, return the cluster colour iff every leaf
    below it shares one cluster id, else ``mixed``."""
    out: Dict[int, tuple] = {}

    def walk(node):
        if node.is_leaf():
            return int(leaf_cl[node.id])
        left = walk(node.get_left())
        right = walk(node.get_right())
        cluster = left if (left == right and left is not None) else None
        out[int(node.id)] = (
            palette[cluster] if cluster is not None else mixed
        )
        return cluster

    walk(to_tree(Z, rd=False))
    return out


def _plot_taskref_panel(
    ax,
    Z: np.ndarray,
    leaf_cl: np.ndarray,
    palette: Dict[int, tuple],
    ylim,
    title: str,
    show_cut_at: Optional[float] = None,
    mixed_color: str = "0.25",
) -> None:
    """Render dendrogram with branches/leaves coloured by ``leaf_cl``.

    Branches whose entire subtree shares one cluster id get that cluster's
    palette colour; mixed branches stay ``mixed_color``. Short coloured
    stems are overdrawn at the bottom so each leaf's own cluster identity
    is visible even when its parent merge is ``mixed``.
    """
    color_map = _link_color_uniform(Z, leaf_cl, palette, mixed_color)
    dend = dendrogram(
        Z,
        ax=ax,
        no_labels=True,
        color_threshold=0,
        above_threshold_color=mixed_color,
        link_color_func=lambda kk: color_map.get(int(kk), mixed_color),
    )
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_yscale("log")

    leaf_order = dend["leaves"]
    y_lo, _ = ax.get_ylim()
    first_h = _first_merge_heights_per_leaf(Z)
    for i, lid in enumerate(leaf_order):
        cl = int(leaf_cl[lid])
        if cl not in palette:
            continue
        h_top = float(first_h[lid])
        if not np.isfinite(h_top) or h_top <= y_lo:
            continue
        ax.plot(
            [5 + 10 * i, 5 + 10 * i],
            [y_lo, h_top],
            color=palette[cl],
            lw=0.9,
            zorder=4,
            solid_capstyle="butt",
        )
    if show_cut_at is not None:
        ax.axhline(show_cut_at, color="black", ls="--", lw=1.4,
                   alpha=1.0, zorder=0)
    ax.set_title(title, fontsize=10)


REF_PHASE_TAGS = {
    "task_test": "tasktestref",
    "rest_post": "restpostref",
}


def _plot_dendrograms_phaseref(
    pat: str,
    band: str,
    results: Dict[str, "object"],
    k: int,
    ref_phase: str,
    out_dir: Path = OUT_DIR,
) -> Optional[Path]:
    """Cross-phase figure where leaves are coloured by the cluster id obtained
    from cutting *ref_phase*'s dendrogram at K = ``k`` communities."""
    ref_r = results.get(ref_phase)
    if ref_r is None:
        return None
    Z_ref = ref_r.linkage_matrix
    n_ref = ref_r.n_nodes

    # Skip when phases disagree on giant-component size — leaf-id alignment
    # would be undefined.
    for phase, r in results.items():
        if r is not None and r.n_nodes != n_ref:
            print(f"  [skip {ref_phase}-ref] {pat} {band}: {phase} has "
                  f"{r.n_nodes} nodes vs {ref_phase} {n_ref}")
            return None

    leaf_cl = fcluster(Z_ref, t=k, criterion="maxclust")
    palette = _build_distinct_palette(int(leaf_cl.max()))

    heights = Z_ref[:, 2]
    cut_idx = n_ref - k
    if 0 < cut_idx < len(heights):
        cut_h = float((heights[cut_idx - 1] + heights[cut_idx]) / 2)
    else:
        cut_h = float(heights[-1])

    tag = REF_PHASE_TAGS.get(ref_phase, f"{ref_phase}ref")
    out_path = out_dir / pat / f"{band}_dendrogram_{tag}_k={k}.pdf"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ylim = _shared_ylim(results)
    fig, axes = plt.subplots(
        1, len(PHASE_LABELS), figsize=(4.5 * len(PHASE_LABELS), 4.5),
        sharey=False,  # see _plot_dendrograms note: sharey=True triggers
                       # scipy's non-positive-ylim warning on log axes.
    )
    axes = np.atleast_1d(axes)

    for ax_idx, (ax, phase) in enumerate(zip(axes, PHASE_LABELS)):
        r = results.get(phase)
        if r is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes, color="#999")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(PHASE_TEX[phase], fontsize=11)
            continue
        title = f"{PHASE_TEX[phase]}  ·  K_ref={k}"
        if phase == ref_phase:
            title += f"  ·  cut={cut_h:.3g}"
        _plot_taskref_panel(
            ax,
            r.linkage_matrix,
            leaf_cl,
            palette,
            ylim,
            title=title,
            show_cut_at=(cut_h if phase == ref_phase else None),
        )
        if ax_idx > 0:
            plt.setp(ax.get_yticklabels(), visible=False)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_psi(
    pat: str, band: str, results: Dict[str, "object"],
    out_dir: Path = OUT_DIR,
) -> Path:
    out_path = out_dir / pat / f"{band}_psi.pdf"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    plotted = False
    for phase in PHASE_LABELS:
        r = results.get(phase)
        if r is None:
            continue
        psi_values, n_communities = compute_partition_stability_index(
            r.linkage_matrix
        )
        if n_communities.size == 0:
            continue
        ax.plot(n_communities, psi_values, "-",
                color=PHASE_COLOR[phase], lw=1.4, label=PHASE_TEX[phase])
        argmax_n = int(n_communities[np.argmax(psi_values)])
        ax.axvline(argmax_n, color=PHASE_COLOR[phase], ls=":", lw=0.9, alpha=0.5)
        plotted = True

    if not plotted:
        ax.text(0.5, 0.5, "all phases missing", ha="center", va="center",
                transform=ax.transAxes, color="#999")

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    ax.set_xlabel(r"$n$ (number of communities)")
    ax.set_ylabel(r"$\Psi(n;\tau{=}1/\lambda_{\max})$")
    ax.set_title(f"{pat} — {band_tex}", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ─── main ──────────────────────────────────────────────────────────────────


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: List[dict] = []
    n_dend = n_psi = n_skip = 0

    for pat in PATIENTS_4PHASE:
        labels = _load_channel_labels(pat, SEEG_DATAPATH) or []
        print(f"{pat} (channels={len(labels)})", flush=True)

        for band in BRAIN_BANDS_NAMES:
            results: Dict[str, object] = {}
            for phase in PHASE_LABELS:
                r = load_lrg_result(pat, phase, band, FC_METHOD)
                results[phase] = r
                if r is None:
                    n_skip += 1
                    continue
                psi_values, n_comms = compute_partition_stability_index(
                    r.linkage_matrix
                )
                psi_max = float(psi_values.max()) if psi_values.size else float("nan")
                argmax_n = int(n_comms[np.argmax(psi_values)]) if psi_values.size else -1
                rows.append(dict(
                    patient=pat,
                    band=band,
                    phase=phase,
                    n_nodes=int(r.n_nodes),
                    psi_argmax_n=argmax_n,
                    psi_max=psi_max,
                    optimal_threshold=float(r.optimal_threshold),
                ))

            if any(v is not None for v in results.values()):
                dp = _plot_dendrograms(pat, band, results, labels)
                pp = _plot_psi(pat, band, results)
                ref_paths: List[Path] = []
                for k in K_REF_LIST:
                    rp = _plot_dendrograms_phaseref(
                        pat, band, results, k, REF_PHASE
                    )
                    if rp is not None:
                        ref_paths.append(rp)
                n_dend += 1
                n_psi += 1
                msg = f"  {band}: wrote {dp} + {pp}"
                for rp in ref_paths:
                    msg += f" + {rp}"
                print(msg)
            else:
                print(f"  {band}: skipped (no LRG cache for any phase)")

    with LEDGER_PATH.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "patient", "band", "phase", "n_nodes",
            "psi_argmax_n", "psi_max", "optimal_threshold",
        ])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {LEDGER_PATH}  "
          f"({len(rows)} rows, {n_dend} dend pdfs, {n_psi} psi pdfs, "
          f"{n_skip} cells skipped)")


if __name__ == "__main__":
    main()
