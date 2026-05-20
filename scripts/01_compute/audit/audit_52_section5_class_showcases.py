#!/usr/bin/env python3
"""Audit 52 — Section 5.6 single-class network reprojection showcases.

For each of the five cross-phase classes (trace / reset / rearrange /
anchor / diffuse), render a single (patient, band) showcase figure
proving that the dendrogram-level class identity projects onto a
coherent (or scattered) block in the LRG ultrametric distance matrix
AND in the FC network drawing.

Per-class chosen cells (locked 2026-05-07):

  trace      Pat_07 beta        6 modules sizes [4, 4, 5, 5, 6, 24]
  reset      Pat_10 low_gamma   5 modules
  rearrange  Pat_03 high_gamma  12 rearr modules, 1 small trace, 0 reset, 0 anchor
  anchor     Pat_02 beta        9 anchor modules
  diffuse    Pat_15 beta        52 diffuse leaves of 118 total

Each figure layout (3 rows x 3 columns):

  row 0     dendrograms at rest_pre / task_test / rest_post; class
            leaves highlighted as colored leaf bars + intra-class
            internal links; per-module shading (largest -> darkest)
  row 1     LRG ultrametric distance matrix D-hat(tau_min) reordered
            with class leaves first (by module-size descending, then
            by leaf-id), complement second; shared vmin/vmax across
            the three phase panels; thin black separator between
            class-block and complement-block; sequential viridis_r
            colormap (lower distance = brighter)
  row 2     FC network at each phase using a shared LRG-SFDP layout
            (computed once on the task_test |ImCoh| matrix); nodes
            colored by class-module membership (per-module shade);
            intra-class edges highlighted in module color with edge
            alpha/width scaled by phase-specific rank-percentile; all
            other edges drawn as faint gray background.

  tau_min = 1 / lambda_max from the cached LRG eigendecomposition.

Outputs
-------
  data/audit/section5_v2_round3_redo/figures/taxonomy_class_<class>.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.visuals.network_templates import compute_layout
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.diagnostics import (
    _compute_propagator,
    _ultrametric_from_propagator,
)

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
import audit_47_kc_trace_network_view as a47
import audit_48_kc_reset_module_view as a48
import audit_49_kc_rearrangement_module_view as a49
import audit_50_kc_anchor_module_view as a50
from audit_47_kc_trace_network_view import find_trace_pairs
from audit_48_kc_reset_module_view import find_reset_pairs
from audit_49_kc_rearrangement_module_view import find_rearrangement_modules
from audit_50_kc_anchor_module_view import find_anchor_triples

# Reuse each audit's draw_network as-is; their per-class
# `module_leaves_for_network` rules are part of their canonical visual
# convention and we want this script to match the originals.
DRAW_NETWORK = {
    "trace":     a47.draw_network,
    "reset":     a48.draw_network,
    "rearrange": a49.draw_network,
    "anchor":    a50.draw_network,
}


PHASES = ("rest_pre", "task_test", "rest_post")
PHASE_LABELS = {
    "rest_pre": "rest pre",
    "task_test": "task",
    "rest_post": "rest post",
}

CLASS_COLOR = {
    "trace":     "#2c7fb8",
    "reset":     "#41ab5d",
    "rearrange": "#fc9272",
    "anchor":    "#807dba",
    "diffuse":   "#8c8c8c",
}
CLASS_CMAPS = {
    "trace":     "Blues",
    "reset":     "Greens",
    "rearrange": "Oranges",
    "anchor":    "Purples",
    "diffuse":   "Greys",
}
GRAY_LEAF = "#c8c8c8"
GRAY_LINK = "#bababa"
LAYOUT_K = 20  # n communities for the LRG-SFDP layout (matches audit_47)

# Locked (patient, band) per class.
SHOWCASES = [
    {"cls": "trace",     "patient": "Pat_07", "band": "beta"},
    {"cls": "reset",     "patient": "Pat_10", "band": "low_gamma"},
    {"cls": "rearrange", "patient": "Pat_03", "band": "high_gamma"},
    {"cls": "anchor",    "patient": "Pat_02", "band": "beta"},
    {"cls": "diffuse",   "patient": "Pat_15", "band": "beta"},
]

OUT_DIR = ROOT / "data/audit/section5_v2_round3_redo/figures"
MATRIX_CMAP = "viridis_r"   # low D = bright


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_phase_lrg(patient: str, phase: str, band: str):
    """Return (Z, eigenvalues, eigenvectors) for the cached LRG result."""
    res = load_lrg_result(patient, phase, band, fc_method="imcoh_abs")
    if res is None or res.linkage_matrix is None:
        raise FileNotFoundError(f"missing LRG result {patient} {phase} {band}")
    if res.eigenvalues is None or res.eigenvectors is None:
        raise FileNotFoundError(
            f"LRG result {patient} {phase} {band} has no eigendecomposition"
        )
    return (
        np.asarray(res.linkage_matrix),
        np.asarray(res.eigenvalues),
        np.asarray(res.eigenvectors),
    )


def D_hat_at_tau_min(eigvals: np.ndarray, eigvecs: np.ndarray) -> np.ndarray:
    tau_min = 1.0 / eigvals[-1]
    K = _compute_propagator(eigvals, eigvecs, tau_min)
    return _ultrametric_from_propagator(K)


# ---------------------------------------------------------------------------
# Module detection (mirrors audit_51b conventions)
# ---------------------------------------------------------------------------

def detect_modules(Z_pre, Z_tt, Z_post, n):
    return {
        "trace":     find_trace_pairs(Z_tt, Z_post, Z_pre, n)[0.0],
        "reset":     find_reset_pairs(Z_pre, Z_post, Z_tt, n)[0.0],
        "rearrange": [m for m in find_rearrangement_modules(Z_pre, Z_tt, Z_post, n)
                      if m["size"] >= 3],
        "anchor":    find_anchor_triples(Z_pre, Z_tt, Z_post, n)[0.0],
    }


def strict_module_leaves(cls: str, module: dict) -> set:
    if cls == "trace":
        return module["leaves_tt"] & module["leaves_post"]
    if cls == "reset":
        return module["leaves_pre"] & module["leaves_post"]
    if cls == "anchor":
        return module["leaves_pre"] & module["leaves_tt"] & module["leaves_post"]
    if cls == "diffuse":
        return module["leaves"]
    return module.get("leaves_post_residual", module["leaves_post"])


def shades_for_class(cls: str, n: int) -> list:
    if n == 0:
        return []
    cmap = plt.get_cmap(CLASS_CMAPS[cls])
    if n == 1:
        return [to_hex(cmap(0.75))]
    return [to_hex(cmap(0.95 - 0.45 * i / (n - 1))) for i in range(n)]


def assemble_class_modules(modules_by_class: dict, n: int):
    non_rearrange_leaves: set = set()
    for cls in ("trace", "reset", "anchor"):
        for m in modules_by_class[cls]:
            non_rearrange_leaves |= strict_module_leaves(cls, m)
    residual_modules = []
    for m in modules_by_class["rearrange"]:
        residual = m["leaves_post"] - non_rearrange_leaves
        if len(residual) >= 3:
            mm = {**m, "leaves_post_residual": residual, "size": len(residual)}
            residual_modules.append(mm)
    modules_by_class["rearrange"] = residual_modules

    classified: set = set()
    for cls in ("trace", "reset", "rearrange", "anchor"):
        for m in modules_by_class[cls]:
            classified |= strict_module_leaves(cls, m)
    diffuse_leaves = set(range(n)) - classified
    modules_by_class["diffuse"] = (
        [{"leaves": diffuse_leaves, "size": len(diffuse_leaves)}]
        if diffuse_leaves else []
    )

    for cls in modules_by_class:
        mods = modules_by_class[cls]
        mods_by_size = sorted(mods, key=lambda m: -len(strict_module_leaves(cls, m)))
        shades = shades_for_class(cls, len(mods_by_size))
        for m, sh in zip(mods_by_size, shades):
            m["_shade"] = sh
        modules_by_class[cls] = mods_by_size

    return modules_by_class


# ---------------------------------------------------------------------------
# Dendrogram coloring (mirrors audit_51b)
# ---------------------------------------------------------------------------

def assign_leaf_for_cell(cls, modules, n):
    out = [None] * n
    # smallest first so the largest module wins (overwrites) per leaf
    for m in sorted(modules, key=lambda x: x.get("size", 0)):
        leaves = strict_module_leaves(cls, m)
        for lf in leaves:
            if out[lf] is None:
                out[lf] = m
    # second pass: any leaf claimed by a larger module overrides smaller assignments
    for m in sorted(modules, key=lambda x: -x.get("size", 0)):
        leaves = strict_module_leaves(cls, m)
        for lf in leaves:
            out[lf] = m
    return out


def make_link_color_func(Z, modules, cls, n):
    leaf_m = assign_leaf_for_cell(cls, modules, n)
    sorted_mods = sorted(modules, key=lambda x: x.get("size", 0))

    subtree_leaves = {i: {i} for i in range(n)}
    link_color = {}
    for i, row in enumerate(Z):
        a, b = int(row[0]), int(row[1])
        leaves = subtree_leaves[a] | subtree_leaves[b]
        new_id = n + i
        subtree_leaves[new_id] = leaves
        link_color[new_id] = GRAY_LINK
        for m in sorted_mods:
            mod_leaves = strict_module_leaves(cls, m)
            if mod_leaves and leaves <= mod_leaves:
                link_color[new_id] = m["_shade"]
                break

    leaf_color = {}
    for lf in range(n):
        leaf_color[lf] = leaf_m[lf]["_shade"] if leaf_m[lf] is not None else GRAY_LEAF

    def fn(node_id):
        if node_id < n:
            return leaf_color[node_id]
        return link_color.get(node_id, GRAY_LINK)

    return fn, leaf_m


def draw_dendrogram(ax, Z, modules, cls, n, *, with_yaxis=True):
    fn, leaf_m = make_link_color_func(Z, modules, cls, n)
    R = dendrogram(
        Z, ax=ax, no_labels=True,
        color_threshold=0,
        above_threshold_color=GRAY_LINK,
        link_color_func=fn,
    )
    merge_heights = sorted(set(Z[:, 2]))
    tmin = 0.0
    tmax = merge_heights[-1] * 1.05
    ax.set_ylim(tmin, tmax)
    ax.set_xticks([])
    if with_yaxis:
        ax.tick_params(axis="y", labelsize=6)
    else:
        ax.tick_params(axis="y", labelleft=False)
    for sn in ("top", "right", "bottom"):
        ax.spines[sn].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)

    leaf_parent_h = {}
    for i, row in enumerate(Z):
        a, b, h = int(row[0]), int(row[1]), float(row[2])
        if a < n:
            leaf_parent_h.setdefault(a, h)
        if b < n:
            leaf_parent_h.setdefault(b, h)

    leaf_order = R["leaves"]
    for i, lf_idx in enumerate(leaf_order):
        m = leaf_m[lf_idx]
        if m is None:
            continue
        x = 5.0 + 10.0 * i
        h_parent = leaf_parent_h.get(lf_idx, tmax * 0.05)
        ax.vlines(x, tmin, h_parent, color=m["_shade"], lw=1.4, zorder=4)


# ---------------------------------------------------------------------------
# Matrix panel
# ---------------------------------------------------------------------------

def class_block_order(cls, modules, n):
    """Return (perm, n_class_block).

    Class leaves first (modules ordered largest -> smallest, leaves
    sorted by id within each module, no duplicates across modules),
    complement leaves second (sorted by id).
    """
    seen: set = set()
    perm: list = []
    for m in modules:  # already sorted largest first by assemble_class_modules
        leaves = sorted(strict_module_leaves(cls, m))
        for lf in leaves:
            if lf not in seen:
                perm.append(lf)
                seen.add(lf)
    n_block = len(perm)
    complement = sorted(set(range(n)) - seen)
    perm.extend(complement)
    return np.asarray(perm, dtype=int), n_block


def draw_matrix(ax, D, perm, n_block, vmin, vmax, *, cmap_name=MATRIX_CMAP):
    n = D.shape[0]
    Dp = D[np.ix_(perm, perm)]
    im = ax.imshow(Dp, cmap=cmap_name, vmin=vmin, vmax=vmax,
                   interpolation="nearest", origin="upper")
    ax.set_xticks([])
    ax.set_yticks([])
    for sn in ("top", "right", "bottom", "left"):
        ax.spines[sn].set_linewidth(0.4)
        ax.spines[sn].set_color("#888")
    if 0 < n_block < n:
        ax.axhline(n_block - 0.5, color="black", lw=0.7, zorder=5)
        ax.axvline(n_block - 0.5, color="black", lw=0.7, zorder=5)
    return im


# ---------------------------------------------------------------------------
# Network panel — delegates to audit_47/48/49/50.draw_network so the §5.6
# showcase inherits each class's existing visual convention. We just
# attach (a) `color_idx` per module so audit_X.draw_network can index
# the palette, and (b) `leaves_post = leaves_post_residual` for rearrange
# residual modules so the network reflects the residual (not raw) set.
# ---------------------------------------------------------------------------

def network_modules_and_palette(cls, modules, n):
    """Return (modules_for_network, palette) compatible with audit_X.draw_network."""
    if cls == "diffuse":
        # No audit_diffuse — synthesize a single 'module' whose
        # `leaves_post` is the diffuse set, then route through
        # audit_49.draw_network (which reads `leaves_post`).
        if not modules:
            return [], {}
        m = modules[0]
        diff_leaves = strict_module_leaves(cls, m)
        synth = {
            "leaves_post": diff_leaves,
            "size": len(diff_leaves),
            "color_idx": 0,
        }
        pal = {0: m["_shade"]}
        return [synth], pal

    pal = {}
    out_mods = []
    for i, m in enumerate(modules):
        m = dict(m)
        m["color_idx"] = i
        if cls == "rearrange" and "leaves_post_residual" in m:
            # audit_49.draw_network reads `leaves_post`; route the
            # residual set through that key for §5.6 consistency.
            m["leaves_post"] = m["leaves_post_residual"]
        pal[i] = m["_shade"]
        out_mods.append(m)
    return out_mods, pal


def draw_class_network(ax, A, pos, cls, modules, n):
    """Render the class's FC-network panel by delegating to the
    corresponding audit_X.draw_network. Diffuse routes through audit_49."""
    fn = DRAW_NETWORK.get(cls, a49.draw_network)
    net_modules, pal = network_modules_and_palette(cls, modules, n)
    fn(ax, A, pos, net_modules, pal, n)


# ---------------------------------------------------------------------------
# Main render
# ---------------------------------------------------------------------------

def render_showcase(cls: str, patient: str, band: str) -> dict:
    Z_pre, ev_pre, V_pre = load_phase_lrg(patient, "rest_pre", band)
    Z_tt,  ev_tt,  V_tt  = load_phase_lrg(patient, "task_test", band)
    Z_post, ev_post, V_post = load_phase_lrg(patient, "rest_post", band)

    A_pre  = np.asarray(load_fc_matrix(patient, "rest_pre",  band, fc_method="imcoh_abs"))
    A_tt   = np.asarray(load_fc_matrix(patient, "task_test", band, fc_method="imcoh_abs"))
    A_post = np.asarray(load_fc_matrix(patient, "rest_post", band, fc_method="imcoh_abs"))

    n = min(V_pre.shape[0], V_tt.shape[0], V_post.shape[0],
            A_pre.shape[0], A_tt.shape[0], A_post.shape[0])

    Z_by_phase = {
        "rest_pre":  Z_pre,
        "task_test": Z_tt,
        "rest_post": Z_post,
    }
    A_by_phase = {
        "rest_pre":  A_pre[:n, :n],
        "task_test": A_tt[:n, :n],
        "rest_post": A_post[:n, :n],
    }
    D_by_phase = {
        "rest_pre":  D_hat_at_tau_min(ev_pre, V_pre)[:n, :n],
        "task_test": D_hat_at_tau_min(ev_tt,  V_tt)[:n, :n],
        "rest_post": D_hat_at_tau_min(ev_post, V_post)[:n, :n],
    }

    modules_by_class = detect_modules(Z_pre, Z_tt, Z_post, n)
    modules_by_class = assemble_class_modules(modules_by_class, n)
    modules = modules_by_class[cls]
    n_modules = len(modules)
    class_leaves = set()
    for m in modules:
        class_leaves |= strict_module_leaves(cls, m)
    n_class_leaves = len(class_leaves)

    perm, n_block = class_block_order(cls, modules, n)

    # Joint percentile clip across the three phases for visual contrast.
    # Off-diagonal only; D-hat blows up where K(tau) ~ 0 (near-disconnected
    # pairs in fully-connected weighted graphs), so the raw min/max wastes
    # the colormap on extreme tail values.
    pooled = []
    for ph in PHASES:
        D = D_by_phase[ph]
        offdiag = D[~np.eye(n, dtype=bool)]
        pooled.append(offdiag)
    pooled = np.concatenate(pooled)
    vmin = float(np.percentile(pooled, 5.0))
    vmax = float(np.percentile(pooled, 95.0))

    # Shared LRG-SFDP layout, computed once on task_test (matches audit_47).
    labels_layout = fcluster(Z_tt, t=LAYOUT_K, criterion="maxclust")
    pos = compute_layout(
        "lrg_sfdp",
        A_by_phase["task_test"],
        community_labels=labels_layout.tolist(),
    )
    pos = np.asarray(pos)[:n]

    fig = plt.figure(figsize=(8.6, 8.4))
    gs = fig.add_gridspec(
        3, 4,
        width_ratios=[1, 1, 1, 0.04],
        height_ratios=[1.0, 1.05, 1.4],
        left=0.06, right=0.94, top=0.93, bottom=0.04,
        hspace=0.16, wspace=0.06,
    )

    # Row 0: dendrograms
    den_axes = []
    for c, ph in enumerate(PHASES):
        sharey = den_axes[0] if den_axes else None
        ax = fig.add_subplot(gs[0, c], sharey=sharey)
        draw_dendrogram(ax, Z_by_phase[ph], modules, cls, n,
                        with_yaxis=(c == 0))
        ax.set_title(PHASE_LABELS[ph], fontsize=10, pad=4)
        if c == 0:
            ax.set_ylabel("merge height", fontsize=8)
        den_axes.append(ax)

    # Row 1: D-hat matrices
    mat_axes = []
    last_im = None
    for c, ph in enumerate(PHASES):
        ax = fig.add_subplot(gs[1, c])
        last_im = draw_matrix(
            ax, D_by_phase[ph], perm, n_block, vmin, vmax,
        )
        if c == 0:
            ax.set_ylabel(r"$\hat D(\tau_{\min})$", fontsize=8)
        mat_axes.append(ax)

    # Colorbar for the matrix row
    cax = fig.add_subplot(gs[1, 3])
    cb = fig.colorbar(last_im, cax=cax)
    cb.ax.tick_params(labelsize=6)
    cb.set_label(r"$\hat D(\tau_{\min})$", fontsize=7)

    # Row 2: FC networks (delegates to audit_47/48/49/50.draw_network)
    for c, ph in enumerate(PHASES):
        ax = fig.add_subplot(gs[2, c])
        draw_class_network(ax, A_by_phase[ph], pos, cls, modules, n)
        if c == 0:
            ax.set_ylabel("FC network", fontsize=8)

    # Class label + cell info
    cls_color = CLASS_COLOR[cls]
    fig.text(
        0.07, 0.965, cls,
        ha="left", va="center",
        fontsize=14, color=cls_color, fontweight="bold",
    )
    if cls == "diffuse":
        info = f"{patient}  {band}   {n_class_leaves}/{n} diffuse leaves"
    else:
        info = (f"{patient}  {band}   {n_modules} module"
                f"{'s' if n_modules != 1 else ''}  ·  "
                f"{n_class_leaves}/{n} class leaves")
    fig.text(0.93, 0.965, info, ha="right", va="center",
             fontsize=8.5, color="#333")

    fig.text(
        0.5, 0.02,
        "ordering: class leaves (block, top-left) | complement (sorted by leaf-id)",
        ha="center", va="center", fontsize=7, color="#666",
    )

    out_path = OUT_DIR / f"taxonomy_class_{cls}.pdf"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)

    return {
        "cls": cls,
        "patient": patient,
        "band": band,
        "n": n,
        "n_modules": n_modules,
        "n_class_leaves": n_class_leaves,
        "n_block": n_block,
        "vmin": float(vmin),
        "vmax": float(vmax),
        "out_path": out_path,
    }


def main() -> int:
    rows = []
    for case in SHOWCASES:
        info = render_showcase(case["cls"], case["patient"], case["band"])
        print(f"  {info['cls']:9s}  {info['patient']}  {info['band']:10s}  "
              f"n={info['n']:3d}  n_modules={info['n_modules']:2d}  "
              f"n_class_leaves={info['n_class_leaves']:3d}/{info['n']}  "
              f"vmin={info['vmin']:.3g}  vmax={info['vmax']:.3g}  "
              f"-> {info['out_path'].name}")
        rows.append(info)

    import pandas as pd
    out_csv = ROOT / "data/audit/section5_v2_round3_redo/tables/taxonomy_class_showcases.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame([
        {k: r[k] for k in ("cls", "patient", "band", "n", "n_modules",
                           "n_class_leaves", "n_block", "vmin", "vmax")}
        for r in rows
    ])
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
