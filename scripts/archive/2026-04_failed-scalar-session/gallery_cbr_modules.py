#!/usr/bin/env python3
"""CBR-module gallery: sanity-check examples of four cross-phase patterns.

For each (patient, band, rpost-internal-node) we compute leaf-set Jaccard
against the best-matching internal node of each other phase tree and
classify the node into one of four patterns:

  * TRACE     — task-born, retained. J_rpre ≤ 0.5 ; J_tl , J_tt ≥ 0.75.
  * RESET     — pre-existing, dissolved in task, reappeared at rest_post.
                J_rpre ≥ 0.75 ; min(J_tl, J_tt) ≤ 0.5.
  * PERSIST   — stable module across all four phases. All J ≥ 0.8.
  * REARRANGE — module unique to rest_post. All three J ≤ 0.4.

For each pattern we rank candidates globally by
``score = size·(size−1)·(1 − size/N)·π·(pattern-specific directional gap)``,
select the top few diverse examples (different patients, bands, scales),
and write one 4-panel dendrogram figure per module under
``data/outputs/figures/cbr_gallery/{pattern}/``.

This is a sanity check: does the phenomenon hold visually across patients
and bands, at multiple scales, before we commit to a production metric.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, to_tree

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST,
)
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, OUTPUTS_ROOT
from lrg_eegfc.utils.metrics import jaccard_leafsets
from lrg_eegfc.workflow.lrg import load_lrg_result


PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_TEX = {"rest_pre":  r"rest$_\mathrm{pre}$",
             "task_learn": r"task$_\mathrm{learn}$",
             "task_test":  r"task$_\mathrm{test}$",
             "rest_post":  r"rest$_\mathrm{post}$"}
SIZE_MIN = 5           # demand non-trivial mesoscopic modules
SIZE_MAX = 60          # modules must not span (almost) the whole tree
BACKGROUND_GRAY = "#bfbfbf"
STRIP_GRAY = "#ededed"


# ──────────────────────────── tree utilities ────────────────────────────

def nodes_with_prominence(Z: np.ndarray, h_floor_frac: float = 1e-6
                           ) -> list[dict]:
    n_leaves = Z.shape[0] + 1
    n_int = Z.shape[0]
    root = to_tree(Z, rd=False)
    nodes: list[dict | None] = [None] * n_int

    def walk(nd) -> frozenset[int]:
        if nd.is_leaf():
            return frozenset([nd.id])
        left = walk(nd.get_left())
        right = walk(nd.get_right())
        leaves = left | right
        row = nd.id - n_leaves
        nodes[row] = {"row": row, "id": nd.id, "h": float(nd.dist),
                       "size": len(leaves), "leaves": leaves}
        return leaves
    walk(root)

    dmax = float(Z[-1, 2])
    h_floor = h_floor_frac * dmax

    for j in range(n_int):
        for col in (0, 1):
            child = int(Z[j, col])
            if child >= n_leaves:
                nodes[child - n_leaves]["parent_row"] = j
    nodes[n_int - 1]["parent_row"] = -1

    for nd in nodes:
        p = nd["parent_row"]
        if p < 0:
            nd["pi"] = np.nan
            nd["h_parent"] = np.nan
        else:
            nd["h_parent"] = float(Z[p, 2])
            nd["pi"] = np.log(nd["h_parent"]) - np.log(max(nd["h"], h_floor))
        nd["h_rel"] = nd["h"] / dmax if dmax > 0 else 0.0
    return nodes


def best_match(leafset: frozenset[int], ref_nodes: list[dict]
                ) -> tuple[int, float]:
    best_j, best_row = 0.0, -1
    for rn in ref_nodes:
        j = jaccard_leafsets(leafset, rn["leaves"])
        if j > best_j:
            best_j, best_row = j, rn["row"]
    return best_row, best_j


def highlight_subset_rows(Z: np.ndarray, target_leaves: frozenset[int]
                           ) -> list[int]:
    n_leaves = Z.shape[0] + 1
    rows_leaves: list[frozenset[int]] = [frozenset()] * Z.shape[0]
    for r in range(Z.shape[0]):
        out = set()
        for col in (0, 1):
            child = int(Z[r, col])
            if child < n_leaves:
                out.add(child)
            else:
                out |= rows_leaves[child - n_leaves]
        rows_leaves[r] = frozenset(out)
    parent_row = [-1] * Z.shape[0]
    for r in range(Z.shape[0]):
        for col in (0, 1):
            child = int(Z[r, col])
            if child >= n_leaves:
                parent_row[child - n_leaves] = r
    out_rows = []
    for r in range(Z.shape[0]):
        if rows_leaves[r].issubset(target_leaves) and len(rows_leaves[r]) >= 2:
            p = parent_row[r]
            if p < 0 or not rows_leaves[p].issubset(target_leaves):
                out_rows.append(r)
    return out_rows


def phase_style(Z: np.ndarray, highlight_rows: dict[int, str]
                 ) -> dict[int, str]:
    n_leaves = Z.shape[0] + 1
    color_by_id: dict[int, str] = {}
    def paint(row, color):
        stack = [row]
        while stack:
            r = stack.pop()
            color_by_id[n_leaves + r] = color
            for col in (0, 1):
                child = int(Z[r, col])
                if child >= n_leaves:
                    stack.append(child - n_leaves)
    for row, color in highlight_rows.items():
        paint(row, color)
    return color_by_id


def plot_one_dendrogram(ax, Z, highlight_rows, title,
                         y_bottom: float,
                         target_leaves: frozenset[int] | None = None,
                         color: str | None = None) -> dict:
    color_by_id = phase_style(Z, highlight_rows)
    dend = dendrogram(
        Z, ax=ax, no_labels=True, color_threshold=0,
        above_threshold_color=BACKGROUND_GRAY,
        link_color_func=lambda k: color_by_id.get(k, BACKGROUND_GRAY),
    )
    ax.set_title(title, fontsize=11)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=9)

    # Overdraw leaf stems for target leaves whose parent is a mixed node
    # (so scipy coloured their stem gray via link_color_func on the parent).
    # Without this, isolated target leaves are visually invisible in the
    # dendrogram even though they belong to the module.
    if target_leaves and color:
        n_leaves = Z.shape[0] + 1
        parent_h = np.full(n_leaves, np.nan)
        for r in range(Z.shape[0]):
            h = float(Z[r, 2])
            for col in (0, 1):
                child = int(Z[r, col])
                if child < n_leaves:
                    parent_h[child] = h
        for i, lid in enumerate(dend["leaves"]):
            if lid in target_leaves:
                h_top = parent_h[lid]
                if np.isfinite(h_top):
                    ax.plot([5 + 10 * i, 5 + 10 * i],
                             [y_bottom, h_top], color=color,
                             linewidth=1.1, zorder=5,
                             solid_capstyle="butt")
    return dend


def plot_module_strip(ax, leaves_order: list[int],
                       target_leaves: frozenset[int], color: str) -> None:
    """Colored bar under a dendrogram marking which leaves are in the module.

    ``leaves_order`` is ``dend['leaves']`` — the display-order list of leaf
    ids from ``scipy.cluster.hierarchy.dendrogram``. scipy places each leaf
    at x = 5 + 10*i for i in display order, so the strip's x-axis spans
    [0, 10·n_leaves] to align with the dendrogram above.
    """
    n = len(leaves_order)
    is_target = np.array([lid in target_leaves for lid in leaves_order])
    ax.set_xlim(0, 10 * n)
    ax.set_ylim(0, 1)
    ax.set_facecolor(STRIP_GRAY)
    # One rectangle per module leaf, width 10 units, centred on 5+10·i.
    for i, is_mod in enumerate(is_target):
        if is_mod:
            ax.axvspan(10 * i, 10 * (i + 1), color=color, alpha=0.95, lw=0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


# ──────────────────────────── pattern logic ────────────────────────────

def classify(J_rpre: float, J_tl: float, J_tt: float) -> str | None:
    # Ordered most-specific → least-specific.
    if J_rpre >= 0.80 and J_tl >= 0.80 and J_tt >= 0.80:
        return "persist"
    if J_rpre <= 0.50 and J_tl >= 0.75 and J_tt >= 0.75:
        return "trace"
    if J_rpre >= 0.75 and min(J_tl, J_tt) <= 0.50:
        return "reset"
    if J_rpre <= 0.40 and J_tl <= 0.40 and J_tt <= 0.40:
        return "rearrange"
    return None


def pattern_score(pattern: str, pi: float, size: int, n_leaves: int,
                   J_rpre: float, J_tl: float, J_tt: float) -> float:
    """size·(size−1)·(1 − size/N) · π · directional gap."""
    bal = size * (size - 1) * (1.0 - size / n_leaves)
    if pattern == "trace":
        gap = (J_tl + J_tt) / 2.0 - J_rpre
    elif pattern == "reset":
        gap = J_rpre - (J_tl + J_tt) / 2.0
    elif pattern == "persist":
        gap = min(J_rpre, J_tl, J_tt)
    elif pattern == "rearrange":
        gap = 1.0 - max(J_rpre, J_tl, J_tt)
    else:
        return 0.0
    return max(pi, 0.0) * bal * max(gap, 0.0)


# ──────────────────────────── candidate scanning ────────────────────────────

def _load_Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return np.asarray(r.linkage_matrix) if r is not None else None


def scan_candidates() -> pd.DataFrame:
    """One row per (patient, band, rpost-node) classified into a pattern."""
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            Zs = {ph: _load_Z(pat, ph, band) for ph in PHASES}
            if any(v is None for v in Zs.values()):
                continue
            N = Zs["rest_post"].shape[0] + 1
            if any(Z.shape[0] + 1 != N for Z in Zs.values()):
                continue
            nodes_rp = nodes_with_prominence(Zs["rest_post"])
            nodes_tt = nodes_with_prominence(Zs["task_test"])
            nodes_tl = nodes_with_prominence(Zs["task_learn"])
            nodes_rpre = nodes_with_prominence(Zs["rest_pre"])
            for η in nodes_rp:
                if η["size"] < SIZE_MIN or η["size"] > min(SIZE_MAX, N - 2):
                    continue
                if not np.isfinite(η["pi"]):
                    continue
                _, J_rpre = best_match(η["leaves"], nodes_rpre)
                _, J_tl = best_match(η["leaves"], nodes_tl)
                _, J_tt = best_match(η["leaves"], nodes_tt)
                pat_cls = classify(J_rpre, J_tl, J_tt)
                if pat_cls is None:
                    continue
                score = pattern_score(pat_cls, η["pi"], η["size"], N,
                                       J_rpre, J_tl, J_tt)
                rows.append({
                    "pattern": pat_cls, "patient": pat, "band": band,
                    "row": η["row"], "size": η["size"],
                    "h_rel": η["h_rel"], "pi": η["pi"],
                    "J_rpre": J_rpre, "J_tl": J_tl, "J_tt": J_tt,
                    "score": score, "n_leaves": N,
                })
    return pd.DataFrame(rows)


def pick_per_patient(df: pd.DataFrame, pattern: str,
                       k_per_patient: int = 2) -> pd.DataFrame:
    """Up to ``k_per_patient`` best examples for each patient of this pattern.

    Prefers band diversity within a patient: we pick the top candidate in the
    strongest band, then the top candidate in a different band (if any), etc.
    Ensures every patient who has the pattern at all is represented.
    """
    sub = df[df["pattern"] == pattern].copy()
    if sub.empty:
        return sub
    # Collapse to at most one row per (patient, band).
    sub = sub.sort_values("score", ascending=False)
    sub = sub.drop_duplicates(subset=["patient", "band"], keep="first")

    picks = []
    for pat, pdf in sub.groupby("patient", sort=False):
        pdf_sorted = pdf.sort_values("score", ascending=False)
        seen_bands: set[str] = set()
        pat_picks = []
        for _, row in pdf_sorted.iterrows():
            if row["band"] in seen_bands:
                continue
            pat_picks.append(row)
            seen_bands.add(row["band"])
            if len(pat_picks) >= k_per_patient:
                break
        picks.extend(pat_picks)

    if not picks:
        return pd.DataFrame(columns=sub.columns)
    out = pd.DataFrame(picks)
    out = out.sort_values(["patient", "score"],
                           ascending=[True, False]).reset_index(drop=True)
    return out


# ──────────────────────────── figure per module ────────────────────────────

PATTERN_COLORS = {
    "trace":     "#2a9d8f",
    "reset":     "#e76f51",
    "persist":   "#1d3557",
    "rearrange": "#8338ec",
}
PATTERN_LABELS = {
    "trace":     "TRACE (task-born, retained)",
    "reset":     "RESET (pre-existing, dissolved in task, reappeared)",
    "persist":   "PERSIST (stable across all phases)",
    "rearrange": "REARRANGE (unique to rest_post)",
}


def plot_module_figure(row: pd.Series, out_path) -> None:
    pat = row["patient"]; band = row["band"]; pattern = row["pattern"]
    Zs = {ph: _load_Z(pat, ph, band) for ph in PHASES}
    nodes_rp = nodes_with_prominence(Zs["rest_post"])
    target = nodes_rp[int(row["row"])]
    leaves = target["leaves"]
    color = PATTERN_COLORS[pattern]

    # Always use subset-rows: colour only internal nodes whose leafset is a
    # strict subset of the module. Never pull in extra leaves. Isolated
    # target leaves whose parent is mixed get no internal-link colour, but
    # they ARE shown in the module strip below — so the count of module
    # leaves is consistent across panels by construction.
    highlights = {ph: {r: color for r in highlight_subset_rows(Zs[ph], leaves)}
                  for ph in PHASES}

    fig = plt.figure(figsize=(18.0, 5.2), dpi=160)
    gs = fig.add_gridspec(
        2, 4,
        height_ratios=[10.0, 0.45],
        hspace=0.02, wspace=0.10,
    )
    tree_axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    strip_axes = [fig.add_subplot(gs[1, i]) for i in range(4)]

    # Pooled y-range across all 4 phase trees so sharey + tight bottom are
    # compatible. The gap between y_bottom and the lowest merge is ≈2 % of
    # the height span — enough breathing room without wasting vertical
    # real estate on empty stems.
    all_h = np.concatenate([Z[:, 2] for Z in Zs.values()])
    h_positive = all_h[all_h > 0]
    hmin = float(h_positive.min()) if h_positive.size else 0.0
    hmax = float(all_h.max())
    span = max(hmax - hmin, hmax * 0.01)
    y_bottom = max(hmin - 0.02 * span, 0.0)
    y_top = hmax + 0.02 * span

    # Share y across the dendrogram row only.
    for ax in tree_axes[1:]:
        ax.sharey(tree_axes[0])

    for i, ph in enumerate(PHASES):
        dend = plot_one_dendrogram(
            tree_axes[i], Zs[ph], highlights[ph], PHASE_TEX[ph],
            y_bottom=y_bottom,
            target_leaves=leaves, color=color,
        )
        plot_module_strip(strip_axes[i], dend["leaves"], leaves, color)
        if i > 0:
            tree_axes[i].tick_params(axis="y", labelleft=False)

    tree_axes[0].set_ylim(y_bottom, y_top)
    tree_axes[0].set_ylabel(r"merge height $h$", fontsize=10)
    strip_axes[0].set_ylabel(f"module\n(|L|={len(leaves)})", fontsize=8,
                              rotation=0, ha="right", va="center", labelpad=32)

    fig.suptitle(
        f"{pat}, band = {BRAIN_BAND_TEX_DICT[band]}   —   "
        f"{PATTERN_LABELS[pattern]}   |   "
        f"size = {int(row['size'])}, "
        f"$h_\\mathrm{{rel}}$(rpost) = {row['h_rel']:.2f}, "
        f"$\\pi$ = {row['pi']:.2f}\n"
        f"$J_\\mathrm{{rpre}}$ = {row['J_rpre']:.2f},  "
        f"$J_\\mathrm{{tl}}$ = {row['J_tl']:.2f},  "
        f"$J_\\mathrm{{tt}}$ = {row['J_tt']:.2f}   "
        r"(bar below each tree marks the exact same " f"{len(leaves)} "
        r"leaves at their dendrogram display positions)",
        fontsize=10, y=1.02,
    )
    fig.savefig(f"{out_path}.pdf", bbox_inches="tight")
    plt.close(fig)


# ──────────────────────────── main ────────────────────────────

def main() -> None:
    print("Scanning candidates across all patients × bands...")
    df = scan_candidates()
    if df.empty:
        raise SystemExit("no candidates found")

    counts = df["pattern"].value_counts().to_dict()
    print(f"Total candidates per pattern: {counts}")

    gallery_root = OUTPUTS_ROOT / "figures" / "cbr_gallery"
    gallery_root.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for pattern in ("trace", "reset", "persist", "rearrange"):
        picks = pick_per_patient(df, pattern, k_per_patient=4)
        print(f"\n=== {pattern.upper()}  ({len(picks)} examples across "
              f"{picks['patient'].nunique() if not picks.empty else 0} "
              f"patients, of {counts.get(pattern, 0)} total candidates) ===")
        if picks.empty:
            print("  (no candidates)")
            continue
        print(picks[["patient", "band", "size", "h_rel", "pi",
                     "J_rpre", "J_tl", "J_tt", "score"]].to_string(index=False))
        # Also report which patients are MISSING this pattern.
        missing = sorted(set(PATIENTS_LIST) - set(picks["patient"]))
        if missing:
            print(f"  patients with no {pattern} module at size ∈ "
                  f"[{SIZE_MIN}, {SIZE_MAX}]: {missing}")
        # Rank across patients by score for file ordering.
        ranked = picks.sort_values("score", ascending=False).reset_index(drop=True)
        rank_map = {(r["patient"], r["band"], r["row"]): i + 1
                     for i, r in ranked.iterrows()}
        for _, row in picks.iterrows():
            rank = rank_map[(row["patient"], row["band"], row["row"])]
            pat_dir = gallery_root / pattern / row["patient"]
            pat_dir.mkdir(parents=True, exist_ok=True)
            name = (f"{row['band']}_size{int(row['size']):03d}_"
                    f"h{row['h_rel']:.2f}")
            out = pat_dir / name
            plot_module_figure(row, out)
            print(f"  saved {out}.pdf")
            manifest_rows.append({
                "pattern": pattern, "rank_in_pattern": rank,
                "patient": row["patient"], "band": row["band"],
                "size": int(row["size"]), "h_rel": row["h_rel"],
                "pi": row["pi"], "J_rpre": row["J_rpre"],
                "J_tl": row["J_tl"], "J_tt": row["J_tt"],
                "score": row["score"],
                "file": f"{pattern}/{row['patient']}/{name}.pdf",
            })

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = gallery_root / "manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    print(f"\nmanifest: {manifest_path}")
    df.to_csv(gallery_root / "all_candidates.csv", index=False)
    print(f"all candidates: {gallery_root / 'all_candidates.csv'}")


if __name__ == "__main__":
    main()
