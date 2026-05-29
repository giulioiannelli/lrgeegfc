#!/usr/bin/env python3
"""Per-patient 3-phase LRG circular-dendrogram persistence triptych.

Standalone version of panel (b) from
``preprint_07_beta_rho_split_figure_test2``.  Renders the three circular
dendrograms (rest_pre, task_test, rest_post) for one exemplar patient,
with a chosen reference phase that sets:

* the rainbow leaf colour ordering (continuous 1D HLS rainbow over the
  reference phase's ``leaves_list`` slot index);
* the optimal-leaf-ordering target for the other two phases (their
  child swaps minimise reference-slot distance, no topology change);
* the per-leaf cophenet correlation gate used to grey out non-preserved
  leaves in the other two phases.  Leaves with
  ``ρ_c = Spearman(D_coph_ref[c, ·], D_coph_X[c, ·]) > 0.35`` keep their
  reference-slot hue; the rest are greyed.

Per-pair cohort similarity reported to stdout is
``ρ^coph(ref, X) = Spearman(D_coph_ref[triu], D_coph_X[triu])``
— the project canonical continuous cophenet correlation
(see [[feedback-no-partition-metrics-use-rho-coph]]); partition-cut
metrics (ARI/NMI/Fowlkes-Mallows) are NOT used.

Usage
-----
    python preprint_15_dendrogram_persistence_triptych.py \\
        --patient Pat_05 --band beta --reference task
    python preprint_15_dendrogram_persistence_triptych.py \\
        --patient Pat_05 --band beta --reference post

Output (PDF only, no PNG sibling)
---------------------------------
data/preprint/figures/<band>/dendrogram_persistence/
    fig_<band>_dendrogram_persistence_<Pat>_<taskref|postref>.pdf
"""
from __future__ import annotations

import argparse
import colorsys
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from scipy.cluster.hierarchy import (
    cophenet, leaves_list, optimal_leaf_ordering,
)
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_patient_metadata
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

PHASES = ("rest_pre", "task_test", "rest_post")
PHASE_TEX = {
    "rest_pre":  r"rest$_{\mathrm{pre}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}

REFERENCE_TO_PHASE = {"task": "task_test", "post": "rest_post"}
REFERENCE_TAG = {"task": "taskref", "post": "postref"}

GATE_MODES = ("unilateral", "exclusive", "continuous", "discount")
GATE_TAG = {
    "unilateral": "uni",
    "exclusive":  "excl",
    "continuous": "cont",
    "discount":   "disc",
}

METRICS = ("spearman", "pearson")
METRIC_TAG = {"spearman": "sp", "pearson": "pe"}
METRIC_FN = {"spearman": spearmanr, "pearson": pearsonr}

PERSIST_THRESHOLD_DEFAULT = 0.35
GRAY_RGBA = (0.80, 0.80, 0.80, 1.0)


def _threshold_tag(theta: float) -> str:
    """Slug-safe filename tag for a [-1, 1] correlation threshold."""
    sign = "n" if theta < 0 else "p"
    return f"t{sign}{abs(theta):.2f}".replace(".", "")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _parse_contact_label(raw: str) -> str:
    """Compact ``<letter><number>`` label — strips comma annotation +
    apostrophes + mojibake.  Matches ``preprint_07`` / panel (b) helper.
    """
    head = str(raw).split(",")[0]
    return "".join(c for c in head if c.isascii() and c.isalnum())


def load_contact_labels(patient: str) -> list[str]:
    meta = load_patient_metadata(patient, SEEG_DATAPATH)
    if meta is None or "label" not in meta.columns:
        return []
    return [_parse_contact_label(s) for s in meta["label"].tolist()]


def _build_rainbow_cmap() -> LinearSegmentedColormap:
    """Custom fixed-luminance HLS rainbow (L=0.42, S=0.90) avoiding the
    near-white yellow band that bleaches subtree-mean branch colours
    against the paper background.  See [[feedback-no-near-white-cmaps]].
    """
    anchors = [colorsys.hls_to_rgb(h, 0.42, 0.90)
               for h in np.linspace(0.0, 0.83, 24)]
    return LinearSegmentedColormap.from_list("rainbow_fixed_L", anchors, N=256)


def _draw_circular_dendrogram(
    ax: plt.Axes, Z: np.ndarray, *,
    leaf_colors=None, leaf_labels=None,
    h_max_override: float | None = None,
    r_outer: float = 1.0, r_inner: float = 0.08,
    line_width: float = 0.95,
    node_size: float = 18.0,
    label_fontsize: float = 4.2,
    label_offset: float = 1.04,
) -> np.ndarray:
    """Circular dendrogram with size-weighted subtree-RGB branch colouring.

    Equiangular leaves in ``leaves_list(Z)`` order around the unit
    circle; internal nodes at ``r = r_outer · (1 − h / h_max)``; each
    merge = two radial spokes from children to merge radius + one CCW
    arc joining them at that radius.  Branch colour at every segment
    is the size-weighted RGB mean of the descendant leaves' colours:
    pure subtrees keep a vibrant hue, mixed subtrees grey out.
    """
    Z = np.asarray(Z)
    N = int(Z.shape[0]) + 1
    h_max = float(h_max_override) if h_max_override is not None else float(Z[-1, 2])
    n_total = 2 * N - 1

    leaf_order = [int(x) for x in leaves_list(Z)]
    n_leaves = len(leaf_order)
    pos_in_order = {lf: idx for idx, lf in enumerate(leaf_order)}

    theta_of_leaf = np.zeros(N)
    for slot, lf in enumerate(leaf_order):
        theta_of_leaf[lf] = 2 * math.pi * slot / n_leaves - math.pi / 2

    node_r = np.zeros(n_total)
    node_theta = np.zeros(n_total)
    idx_min = np.zeros(n_total, dtype=int)
    idx_max = np.zeros(n_total, dtype=int)

    for i in range(N):
        node_r[i] = r_outer
        node_theta[i] = theta_of_leaf[i]
        idx_min[i] = pos_in_order[i]
        idx_max[i] = pos_in_order[i]

    def _angle(slot_f: float) -> float:
        return 2 * math.pi * slot_f / n_leaves - math.pi / 2

    for r in range(Z.shape[0]):
        node_id = N + r
        h = float(Z[r, 2])
        c1, c2 = int(Z[r, 0]), int(Z[r, 1])
        node_r[node_id] = r_outer - (h / max(h_max, 1e-12)) * (r_outer - r_inner)
        idx_min[node_id] = min(idx_min[c1], idx_min[c2])
        idx_max[node_id] = max(idx_max[c1], idx_max[c2])
        node_theta[node_id] = _angle(0.5 * (idx_min[node_id] + idx_max[node_id]))

    if leaf_colors is None:
        leaf_rgb_arr = np.full((N, 3), 0.35)
    else:
        leaf_rgb_arr = np.array([to_rgb(c) for c in leaf_colors])
    node_rgb = np.zeros((n_total, 3))
    node_size_w = np.zeros(n_total)
    for i in range(N):
        node_rgb[i] = leaf_rgb_arr[i]
        node_size_w[i] = 1.0
    for r in range(Z.shape[0]):
        node_id = N + r
        c1, c2 = int(Z[r, 0]), int(Z[r, 1])
        s1, s2 = float(node_size_w[c1]), float(node_size_w[c2])
        node_rgb[node_id] = (s1 * node_rgb[c1] + s2 * node_rgb[c2]) / (s1 + s2)
        node_size_w[node_id] = s1 + s2

    for r in range(Z.shape[0]):
        node_id = N + r
        c1, c2 = int(Z[r, 0]), int(Z[r, 1])
        if idx_min[c1] <= idx_min[c2]:
            c_left, c_right = c1, c2
        else:
            c_left, c_right = c2, c1
        merge_r = node_r[node_id]
        t_left = node_theta[c_left]
        t_right = node_theta[c_right]
        for c in (c_left, c_right):
            cr = node_r[c]; th = node_theta[c]
            rgb_seg = tuple(node_rgb[c])
            ax.plot([cr * math.cos(th), merge_r * math.cos(th)],
                    [cr * math.sin(th), merge_r * math.sin(th)],
                    color=rgb_seg, lw=line_width, zorder=3,
                    solid_capstyle="round")
        dth = (t_right - t_left) % (2 * math.pi)
        arc_th = t_left + np.linspace(0.0, dth, 48)
        rgb_arc = tuple(node_rgb[node_id])
        ax.plot(merge_r * np.cos(arc_th), merge_r * np.sin(arc_th),
                color=rgb_arc, lw=line_width, zorder=3,
                solid_capstyle="round")

    leaf_x = r_outer * np.cos(theta_of_leaf)
    leaf_y = r_outer * np.sin(theta_of_leaf)
    if leaf_colors is None:
        leaf_colors = ["0.35"] * N
    ax.scatter(leaf_x, leaf_y, c=leaf_colors,
               s=node_size, edgecolors="white", linewidths=0.55,
               zorder=5)

    if leaf_labels is not None:
        for i in range(min(N, len(leaf_labels))):
            if not leaf_labels[i]:
                continue
            th = float(theta_of_leaf[i])
            lr = r_outer * label_offset
            x = lr * math.cos(th); y = lr * math.sin(th)
            ang = math.degrees(th)
            if -90 < ang <= 90:
                rot, ha = ang, "left"
            else:
                rot, ha = ang + 180, "right"
            ax.text(x, y, leaf_labels[i],
                    rotation=rot, rotation_mode="anchor",
                    ha=ha, va="center",
                    fontsize=label_fontsize,
                    color=leaf_colors[i],
                    zorder=6, clip_on=False)

    ax.set_xlim(-r_outer * 1.13, r_outer * 1.13)
    ax.set_ylim(-r_outer * 1.13, r_outer * 1.13)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return theta_of_leaf


# ---------------------------------------------------------------------------
# Triptych
# ---------------------------------------------------------------------------
def render_triptych(patient: str, band: str, reference: str,
                    gate_mode: str = "unilateral",
                    threshold: float = PERSIST_THRESHOLD_DEFAULT,
                    metric: str = "spearman") -> Path:
    """Render the 3-phase circular-dendrogram triptych for one patient.

    Parameters
    ----------
    patient : ``"Pat_05"``, …
    band : ``"beta"``, …
    reference : ``"task"`` (color ring = task_test) or ``"post"``
        (color ring = rest_post).  Sets the OLO target + per-leaf
        cophenet correlation gate against which the other two phases
        are coloured / greyed.
    gate_mode : ``"unilateral"`` (default — colour leaf c in panel X iff
        ρ_c(ref, X) > θ, ignoring the other comparison) or ``"exclusive"``
        (colour leaf c in panel X iff ρ_c(ref, X) > θ AND
        ρ_c(ref, Y) ≤ θ for the other comparison Y, isolating leaves
        unique to X — removes the "always preserved" background that
        otherwise washes out the contrast between the two panels).
    threshold : per-leaf cophenet-correlation cut θ (default 0.35).
        Counts depend on this — leaves with ρ_c > θ keep their colour.
    metric : ``"spearman"`` (default — rank-only on the 117 cophenet
        distances per leaf, immune to absolute-height scaling; same
        family as the §5.3 ``ρ_split^coph`` headline) or ``"pearson"``
        (uses the raw cophenet heights, picks up rank agreement AND
        merge-height scaling).
    """
    if reference not in REFERENCE_TO_PHASE:
        raise ValueError(
            f"reference must be 'task' or 'post', got {reference!r}"
        )
    if gate_mode not in GATE_MODES:
        raise ValueError(
            f"gate_mode must be in {GATE_MODES}, got {gate_mode!r}"
        )
    if metric not in METRICS:
        raise ValueError(
            f"metric must be in {METRICS}, got {metric!r}"
        )
    ref_phase = REFERENCE_TO_PHASE[reference]
    ref_tag = REFERENCE_TAG[reference]
    gate_tag = GATE_TAG[gate_mode]
    thr_tag = _threshold_tag(threshold)
    metric_tag = METRIC_TAG[metric]
    corr_fn = METRIC_FN[metric]
    other_phases = [p for p in PHASES if p != ref_phase]

    Zs: dict[str, np.ndarray] = {}
    for ph in PHASES:
        lrg = load_lrg_result(patient, ph, band, "imcoh_abs")
        if lrg is None or lrg.linkage_matrix is None:
            raise FileNotFoundError(
                f"No LRG cache for {patient}/{band}/{ph}/imcoh_abs"
            )
        Zs[ph] = np.asarray(lrg.linkage_matrix, dtype=float)

    # (1) OLO reference against its own cophenet → smooth rainbow ring.
    Zs[ref_phase] = optimal_leaf_ordering(Zs[ref_phase],
                                          cophenet(Zs[ref_phase]))
    ref_leaf_order = [int(x) for x in leaves_list(Zs[ref_phase])]
    N = len(ref_leaf_order)
    ref_slot_of = np.zeros(N, dtype=int)
    for slot, c in enumerate(ref_leaf_order):
        if c < N:
            ref_slot_of[c] = slot

    # (2) OLO each non-reference phase against reference's slot
    # distance ``|slot_ref(i) - slot_ref(j)|`` → fewer crossings.
    ref_dist_2d = np.abs(
        ref_slot_of[:, None] - ref_slot_of[None, :]
    ).astype(float)
    np.fill_diagonal(ref_dist_2d, 0.0)
    y_ref_slot = squareform(ref_dist_2d, checks=False)
    for ph in other_phases:
        Zs[ph] = optimal_leaf_ordering(Zs[ph], y_ref_slot)

    # (3) Reference-slot continuous rainbow (fixed-luminance HLS).
    cmap_seq = _build_rainbow_cmap()
    leaf_colors_ref = [cmap_seq(ref_slot_of[c] / max(N - 1, 1))
                       for c in range(N)]

    # (4) Per-leaf cophenet correlation against the reference, computed
    # once per non-reference phase.  Used by the exclusive gate below.
    # Metric is either Spearman (rank-only, default, matches §5.3) or
    # Pearson (uses raw cophenet heights).
    Dc = {ph: squareform(cophenet(Zs[ph])) for ph in PHASES}
    rho_per_leaf: dict[str, np.ndarray] = {}
    for ph in other_phases:
        rhos = np.zeros(N)
        for c in range(N):
            mask = np.ones(N, dtype=bool)
            mask[c] = False
            rho_c, _ = corr_fn(Dc[ref_phase][c, mask],
                               Dc[ph][c, mask])
            rhos[c] = float(rho_c) if np.isfinite(rho_c) else 0.0
        rho_per_leaf[ph] = rhos

    # (5) Per-leaf colour assignment in each non-reference panel.
    # Four modes:
    #
    # * 'unilateral' (binary threshold): colour leaf c in panel X iff
    #   ρ_c(ref, X) > θ.  Legacy panel-(b) behaviour.
    # * 'exclusive' (binary threshold + AND-NOT): colour in X iff
    #   ρ_c(ref, X) > θ AND ρ_c(ref, Y) ≤ θ for the other comparison Y.
    #   Isolates leaves unique to X.
    # * 'continuous' (no threshold): leaf colour in panel X is the
    #   reference rainbow colour blended toward grey by
    #   1 - clip(ρ_c(ref, X), 0, 1).  Preserved leaves stay full
    #   rainbow; weakly-preserved fade; ρ_c ≤ 0 → fully grey.  No
    #   cherry-picked cut.  ``threshold`` is ignored in this mode.
    # * 'discount' (continuous + baseline-subtracted): leaf colour in
    #   panel X is the reference rainbow colour blended toward grey by
    #   1 - clip(ρ_c(ref, X) - ρ_c(ref, Y), 0, 1).  Anything similar to
    #   the reference in BOTH non-reference panels (trivial structural
    #   baseline) cancels to grey; only leaves where the reference is
    #   genuinely closer to X than to Y (i.e. memory specifically
    #   carried by X) retain colour.  Most-honest visual for "X carries
    #   trace of ref that Y does not".  ``threshold`` is ignored.
    leaf_colors_per_phase: dict[str, list] = {ref_phase: leaf_colors_ref}
    n_colored: dict[str, int] = {ref_phase: N}
    n_both = 0
    n_neither = 0
    for c in range(N):
        x_close = [rho_per_leaf[ph][c] > threshold
                   for ph in other_phases]
        if all(x_close):
            n_both += 1
        elif not any(x_close):
            n_neither += 1

    gray_rgb = np.array(GRAY_RGBA[:3])

    for ph_x in other_phases:
        ph_y = [p for p in other_phases if p != ph_x][0]
        cs_ph: list = []
        n_kept = 0
        for c in range(N):
            if gate_mode == "continuous":
                # Blend each leaf's rainbow colour toward grey by its
                # own ρ_c.  No threshold; rho_c ≤ 0 hits full grey,
                # rho_c = 1 stays full rainbow.
                rho_c = rho_per_leaf[ph_x][c]
                f = float(max(0.0, min(1.0, rho_c)))
                rgb_rainbow = np.array(
                    matplotlib.colors.to_rgb(leaf_colors_ref[c])
                )
                blended = f * rgb_rainbow + (1 - f) * gray_rgb
                cs_ph.append((blended[0], blended[1], blended[2], 1.0))
                if f > 0:
                    n_kept += 1
            elif gate_mode == "discount":
                # Subtract baseline ρ_c(ref, Y) from ρ_c(ref, X) so any
                # leaf similarly close to ref in BOTH non-ref panels
                # cancels to grey; only the differential — "ref is
                # closer to X than to Y at leaf c" — keeps colour.
                delta = (rho_per_leaf[ph_x][c]
                         - rho_per_leaf[ph_y][c])
                f = float(max(0.0, min(1.0, delta)))
                rgb_rainbow = np.array(
                    matplotlib.colors.to_rgb(leaf_colors_ref[c])
                )
                blended = f * rgb_rainbow + (1 - f) * gray_rgb
                cs_ph.append((blended[0], blended[1], blended[2], 1.0))
                if f > 0:
                    n_kept += 1
            else:
                x_close = rho_per_leaf[ph_x][c] > threshold
                y_close = rho_per_leaf[ph_y][c] > threshold
                if gate_mode == "exclusive":
                    keep = x_close and not y_close
                else:  # 'unilateral'
                    keep = x_close
                if keep:
                    cs_ph.append(leaf_colors_ref[c])
                    n_kept += 1
                else:
                    cs_ph.append(GRAY_RGBA)
        leaf_colors_per_phase[ph_x] = cs_ph
        n_colored[ph_x] = n_kept

    contact_labels = load_contact_labels(patient)[:N]

    # Unified h_max across phases → radii directly comparable.
    h_max_unified = float(max(Zs[ph][-1, 2] for ph in PHASES))

    # ---------------------------------------------------------------------
    # Figure — 3 circular dendrograms side-by-side, no extra panels.
    # ---------------------------------------------------------------------
    fig = plt.figure(figsize=(15.0, 5.4))
    gs = fig.add_gridspec(
        1, 3, wspace=0.04,
        left=0.02, right=0.98, top=0.93, bottom=0.04,
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_map = dict(zip(PHASES, axes))

    for ph in PHASES:
        _draw_circular_dendrogram(
            ax_map[ph], Zs[ph],
            leaf_colors=leaf_colors_per_phase[ph],
            leaf_labels=contact_labels,
            h_max_override=h_max_unified,
            r_outer=1.0, r_inner=0.08,
            line_width=0.95,
            node_size=18.0,
            label_fontsize=4.2,
            label_offset=1.04,
        )
        # Star the reference panel so the reader sees which phase
        # is setting the colour map.
        title = PHASE_TEX[ph]
        if ph == ref_phase:
            title = title + r"$^{\,\star}$"
        ax_map[ph].set_title(title, fontsize=12.0, pad=4)

    # Numerics → stdout (no in-axes text — keeps the rings big).  The
    # global per-pair ρ^coph stays Spearman regardless of the per-leaf
    # metric, because §5.3 ρ_split^coph is the project canonical metric;
    # the per-leaf gate is the choice the CLI exposes.
    iu = np.triu_indices(N, k=1)
    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    print(f"  {patient} {band_tex} reference={ref_phase} (★) "
          f"gate={gate_mode} metric={metric}:")
    for ph in other_phases:
        rho_pair, _ = spearmanr(Dc[ref_phase][iu], Dc[ph][iu])
        if gate_mode == "continuous":
            rho_mean = float(np.mean(rho_per_leaf[ph]))
            rho_pos = int((rho_per_leaf[ph] > 0).sum())
            print(f"    ρ^coph({ref_phase}, {ph}) = {rho_pair:+.4f}, "
                  f"mean ρ_c[{metric}] = {rho_mean:+.3f}, "
                  f"ρ_c > 0 in {rho_pos}/{N} leaves "
                  f"(continuous fade — no threshold)")
        elif gate_mode == "discount":
            ph_other = [p for p in other_phases if p != ph][0]
            delta = rho_per_leaf[ph] - rho_per_leaf[ph_other]
            delta_mean = float(np.mean(delta))
            delta_pos = int((delta > 0).sum())
            delta_max = float(np.max(delta))
            print(f"    ρ^coph({ref_phase}, {ph}) = {rho_pair:+.4f}, "
                  f"mean Δρ_c[{metric}] vs {ph_other} = {delta_mean:+.3f}, "
                  f"Δρ_c > 0 in {delta_pos}/{N} leaves "
                  f"(max Δρ_c = {delta_max:+.3f}, discount — no threshold)")
        else:
            print(f"    ρ^coph({ref_phase}, {ph}) = {rho_pair:+.4f}, "
                  f"colored {n_colored[ph]}/{N} leaves "
                  f"(ρ_c[{metric}] > {threshold})")
    if gate_mode not in ("continuous", "discount"):
        print(f"    Venn at ρ_c[{metric}] > {threshold}: "
              f"both = {n_both}/{N}, neither = {n_neither}/{N}")

    out_dir = (ROOT / "data" / "preprint" / "figures" / band
               / "dendrogram_persistence")
    out_dir.mkdir(parents=True, exist_ok=True)
    # continuous/discount modes don't use the threshold — omit it from
    # the filename to avoid implying it's part of the rendering.
    if gate_mode in ("continuous", "discount"):
        fname = (f"fig_{band}_dendrogram_persistence_{patient}_"
                 f"{ref_tag}_{gate_tag}_{metric_tag}.pdf")
    else:
        fname = (f"fig_{band}_dendrogram_persistence_{patient}_"
                 f"{ref_tag}_{gate_tag}_{metric_tag}_{thr_tag}.pdf")
    out = out_dir / fname
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


def main() -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient", default="Pat_05", choices=COHORT)
    parser.add_argument("--band", default="beta",
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    parser.add_argument(
        "--reference", default="task",
        choices=["task", "post"],
        help="Reference phase that sets the rainbow ring + per-leaf "
             "cophenet correlation gate.  'task' = task_test (default, "
             "matches preprint_07 panel b); 'post' = rest_post.",
    )
    parser.add_argument(
        "--gate-mode", default="unilateral",
        choices=list(GATE_MODES),
        help="Per-leaf colour gate.  'unilateral' (default, legacy): "
             "colour leaf c in panel X iff ρ_c(ref, X) > θ.  "
             "'exclusive': colour leaf c in panel X iff "
             "ρ_c(ref, X) > θ AND ρ_c(ref, Y) ≤ θ for the other "
             "comparison Y — isolates leaves unique to X by removing "
             "the 'preserved everywhere' background.  "
             "'continuous': blend each leaf toward grey by "
             "1 - clip(ρ_c(ref, X), 0, 1), no threshold.  "
             "'discount': blend each leaf toward grey by "
             "1 - clip(ρ_c(ref, X) - ρ_c(ref, Y), 0, 1) — subtracts "
             "the baseline similarity present in BOTH non-ref panels, "
             "so only the differential memory carried by X retains "
             "colour.  Most-honest visual for 'X carries trace of ref "
             "that Y does not'.  Threshold is ignored.",
    )
    parser.add_argument(
        "--threshold", type=float, default=PERSIST_THRESHOLD_DEFAULT,
        help=(f"Per-leaf cophenet-correlation cut θ "
              f"(default {PERSIST_THRESHOLD_DEFAULT}).  "
              "Leaves with ρ_c > θ keep their colour; the rest grey "
              "out.  Both kept-count and Venn breakdown shift with θ; "
              "filename embeds the value so you can sweep it."),
    )
    parser.add_argument(
        "--metric", default="spearman",
        choices=list(METRICS),
        help="Per-leaf correlation metric.  'spearman' (default, "
             "rank-only on the 117 cophenet distances, immune to "
             "absolute-height scaling, matches §5.3 ρ_split^coph "
             "family).  'pearson' (raw cophenet heights, picks up "
             "rank agreement AND merge-height scaling).",
    )
    args = parser.parse_args()
    return render_triptych(args.patient, args.band, args.reference,
                           gate_mode=args.gate_mode,
                           threshold=args.threshold,
                           metric=args.metric)


if __name__ == "__main__":
    main()
