#!/usr/bin/env python3
"""Section 5.2 — brain-based unified β figure for the preprint.

Six tiles arranged as 4 columns × 4 rows.  Rows 1–3 are *phase-resolved*
visualisations of the same exemplar patient (Pat_06, cohort-median ρ_split
at β); row 4 is a cohort-level anatomical summary panel.

Layout
------
Row A (4 cols)    LRG ultrametric distance matrices D̂(τ) per phase.
                  Same leaf ordering: hierarchical clustering of the
                  task_test matrix.  Temporal-lobe leaves marked on the
                  axes.

Row B (4 cols)    Dendrograms per phase, same leaf ordering, leaves
                  coloured by Desikan-Killany region group
                  (temporal-lobe / other-cortex / Wm-Unk).

Row C (4 cols)    Leading-eigenvector amplitude on a 2D axial scatter
                  of the patient's contacts (electrode positions).
                  Same colour scale across phases.

Row D (full-width) Cohort-level glass-brain projections (axial
                  + sagittal MIP) showing β trace-flagged contacts
                  pooled across patients.  Two layers: per-pair (CTM)
                  trace-leaves and hierarchical (KC λ=0 topology)
                  trace-leaves.

Reads
-----
- data/cache/imcoh_lrg/Pat_06/beta_{rest_pre,task_learn,task_test,rest_post}_lrg_imcoh-abs.npz
- data/raw/stereoeeg_patients/Pat_06/{channel_labels.csv, implant_pat_06.csv}
- data/audit/lrg_localization_anatomy/per_trace_leaf.csv  (per-pair leaves)
- data/kc_topology_anatomy/per_topology_trace_leaf.csv  (hierarchical leaves)

Writes
------
data/reports/section_5_lrg_trace/headline/figures/section5_2_beta_unified.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, leaves_list
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()


# --------------------------------------------------------------------------
# Anatomy helpers (lift from audit_57_kc_topology_anatomy)
# --------------------------------------------------------------------------
TEMPORAL_LOBE = {
    "ctx-lh-superiortemporal", "ctx-rh-superiortemporal",
    "ctx-lh-middletemporal", "ctx-rh-middletemporal",
    "ctx-lh-inferiortemporal", "ctx-rh-inferiortemporal",
    "ctx-lh-fusiform", "ctx-rh-fusiform",
    "ctx-lh-bankssts", "ctx-rh-bankssts",
    "ctx-lh-entorhinal", "ctx-rh-entorhinal",
    "ctx-lh-parahippocampal", "ctx-rh-parahippocampal",
    "ctx-lh-temporalpole", "ctx-rh-temporalpole",
    "ctx-lh-transversetemporal", "ctx-rh-transversetemporal",
    "Hip",
}


def _normalize_label(raw: str) -> str:
    """'A 1,G2' -> 'A1'; drop reference suffix and collapse spaces.
    Mirrors `audit_44_anatomical_mapping.normalize_label`."""
    return str(raw).split(",")[0].replace(" ", "").upper()


def _parse_dk_dominant(dk) -> str:
    if pd.isna(dk):
        return "Unk"
    return str(dk).split(",")[0].strip()


def load_anatomy(patient: str) -> pd.DataFrame:
    """Return per-leaf anatomy with [leaf_id, channel, x, y, z, region].
    Coordinates are converted from custom 1000x units to mm by /1000 — same
    convention used by `audit_51_anatomy_2d_disclosed`."""
    from lrg_eegfc.utils.io.patient import PATIENT_CHANNEL_DROP
    pat_dir = ROOT / "data" / "raw" / "stereoeeg_patients" / patient
    labels = pd.read_csv(pat_dir / "channel_labels.csv")
    impl = pd.read_csv(pat_dir / f"implant_pat_{patient[-2:]}.csv")
    dk_col = next(c for c in impl.columns if c.startswith("Desikan"))
    impl["norm_label"] = impl["label"].astype(str).apply(_normalize_label)
    impl_map = impl.set_index("norm_label")

    label_drops = PATIENT_CHANNEL_DROP.get(patient, {}).get("__labels__", [])
    if label_drops:
        labels = labels.drop(index=list(label_drops)).reset_index(drop=True)

    rows = []
    for leaf_id, raw in enumerate(labels["label"].astype(str).tolist()):
        ch = _normalize_label(raw)
        rec = {"leaf_id": int(leaf_id), "channel": ch, "channel_raw": raw,
               "x": np.nan, "y": np.nan, "z": np.nan, "region": "Unk"}
        if ch in impl_map.index:
            r = impl_map.loc[ch]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            try:
                # Implant CSVs use 1000x scaled mm units;
                # divide by 1000 → MNI-mm-like.
                rec["x"] = float(str(r["x"]).replace(",", ".")) / 1000.0
                rec["y"] = float(str(r["y"]).replace(",", ".")) / 1000.0
                rec["z"] = float(str(r["z"]).replace(",", ".")) / 1000.0
            except (ValueError, TypeError):
                pass
            rec["region"] = _parse_dk_dominant(r[dk_col])
        rows.append(rec)
    return pd.DataFrame(rows)


def region_group(region: str) -> str:
    if region in TEMPORAL_LOBE:
        return "temporal"
    if region.startswith("ctx-"):
        return "other_cortex"
    return "wm_unk"


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
EXEMPLAR = "Pat_06"
BAND = "beta"
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_LABEL = {
    "rest_pre": r"rsPre",
    "task_learn": r"task$_{\rm learn}$",
    "task_test": r"task$_{\rm test}$",
    "rest_post": r"rsPost",
}

PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

LRG_CACHE = ROOT / "data" / "cache" / "imcoh_lrg"
PER_PAIR_LEAVES = ROOT / "data" / "audit" / "lrg_localization_anatomy" / "per_trace_leaf.csv"
KC_TOPOLOGY_LEAVES = ROOT / "data" / "kc_topology_anatomy" / "per_topology_trace_leaf.csv"

OUT = ROOT / "data" / "reports" / "section_5_lrg_trace" / "headline" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT / "section5_2_beta_unified.pdf"


# --------------------------------------------------------------------------
# Plotting setup
# --------------------------------------------------------------------------
COL_TEMPORAL = "#c0392b"   # temporal lobe
COL_OTHER_CORTEX = "#7f8c8d"
COL_WM_UNK = "#dddddd"

COL_PERPAIR = "#1f3d6e"
COL_HIERARCH = "#7d3c98"
COL_BG = "#cccccc"

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def load_lrg(patient: str, phase: str, band: str) -> dict:
    """Load LRG cache for one (patient, phase, band).

    The cache stores ρ̂(τ) as a flat upper-triangular vector in
    ``ultrametric_matrix``.  We expose both ρ̂ and the dissimilarity
    1−ρ̂ for plotting (1−ρ̂ ∈ [0,1], close = 0 = dark)."""
    f = LRG_CACHE / patient / f"{band}_{phase}_lrg_imcoh-abs.npz"
    d = np.load(f)
    n = int(d["n_nodes"])
    rho = d["ultrametric_matrix"]
    R = squareform(rho, checks=False)
    np.fill_diagonal(R, 1.0)
    Z = d["linkage_matrix"]
    eig_vals = d["eigenvalues"]
    eig_vecs = d["eigenvectors"]
    return dict(n=n, R=R, Z=Z, eig_vals=eig_vals, eig_vecs=eig_vecs)


def draw_custom_dendrogram(ax, Z, leaf_order, leaf_color, leaf_y, line_color="#888888"):
    """Render a dendrogram with leaves in the specified order.

    `leaf_order` is an array of leaf indices giving the desired
    left-to-right placement. `leaf_color` is an array of length n with
    the colour to draw each leaf marker.  `leaf_y` is the y-position
    at which leaf markers are drawn (must be > 0 if axis is log-scale).
    Internal merges are drawn at the average x of their children;
    crossings mean the imposed order differs from the natural linkage
    order — that is the desired diagnostic.
    """
    n = Z.shape[0] + 1
    leaf_x = np.zeros(n)
    for x, leaf in enumerate(leaf_order):
        leaf_x[int(leaf)] = float(x)

    node_x = np.zeros(2 * n - 1)
    node_y = np.zeros(2 * n - 1)
    node_x[:n] = leaf_x
    node_y[:n] = leaf_y
    for k in range(n - 1):
        a = int(Z[k, 0]); b = int(Z[k, 1])
        node_y[n + k] = float(Z[k, 2])
        node_x[n + k] = (node_x[a] + node_x[b]) / 2.0

    for k in range(n - 1):
        a = int(Z[k, 0]); b = int(Z[k, 1])
        h = float(Z[k, 2])
        ya = node_y[a]; yb = node_y[b]
        xa = node_x[a]; xb = node_x[b]
        ax.plot([xa, xa], [ya, h], color=line_color, lw=0.45, zorder=2)
        ax.plot([xb, xb], [yb, h], color=line_color, lw=0.45, zorder=2)
        ax.plot([xa, xb], [h, h], color=line_color, lw=0.45, zorder=2)

    for leaf in range(n):
        ax.plot([leaf_x[leaf]], [leaf_y], marker="|",
                color=leaf_color[leaf],
                markersize=4, markeredgewidth=1.2, zorder=4)
    ax.set_xticks([])
    ax.set_xlim(-1, n)


# --------------------------------------------------------------------------
# Row A: D̂(τ) per phase
# --------------------------------------------------------------------------
def panel_row_A(axes, phase_data, leaf_order, anatomy: pd.DataFrame) -> None:
    """ρ̂(τ) heatmap reordered by the task_test linkage so block structure
    aligns visually under task_test. ρ̂ ∈ [0, 1]; bright = strong
    communication (low diffusion distance)."""
    # Shared color scale: percentile clip at 99% off-diagonal across phases.
    vmax_candidates = []
    for d in phase_data.values():
        R = d["R"].copy()
        np.fill_diagonal(R, np.nan)
        vmax_candidates.append(float(np.nanpercentile(R, 99)))
    vmax = max(vmax_candidates)
    vmin = 0.0
    temporal_mask = anatomy.region.apply(lambda r: r in TEMPORAL_LOBE).values
    for ax, phase in zip(axes, PHASES):
        R = phase_data[phase]["R"][np.ix_(leaf_order, leaf_order)]
        im = ax.imshow(R, cmap="magma_r", vmin=vmin, vmax=vmax,
                       interpolation="nearest", aspect="equal")
        ax.set_xticks([])
        ax.set_yticks([])
        # Tick the temporal-lobe leaves on left axis as red strokes.
        for li in np.where(temporal_mask[leaf_order])[0]:
            ax.plot([-3.0, -0.5], [li, li], color=COL_TEMPORAL,
                    lw=0.7, zorder=3, clip_on=False)
        ax.set_title(PHASE_LABEL[phase], fontsize=10)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    div = make_axes_locatable(axes[-1])
    cax = div.append_axes("right", size="5%", pad=0.10)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(r"$\hat\rho(\tau)$", fontsize=8)
    cbar.ax.tick_params(labelsize=7)


# --------------------------------------------------------------------------
# Row B: dendrograms per phase
# --------------------------------------------------------------------------
def panel_row_B(axes, phase_data, leaf_order, anatomy: pd.DataFrame) -> None:
    """Custom dendrograms with leaves pinned to the task_test order.
    Crossings = reorganisation relative to task_test."""
    region_arr = anatomy.region.apply(region_group).values
    color_for_group = {"temporal": COL_TEMPORAL,
                       "other_cortex": COL_OTHER_CORTEX,
                       "wm_unk": COL_WM_UNK}
    leaf_color = np.array([color_for_group[g] for g in region_arr], dtype=object)

    y_max = max(d["Z"][:, 2].max() for d in phase_data.values())
    y_min = max(1e-3, min(d["Z"][:, 2].min() for d in phase_data.values()) * 0.8)

    leaf_y = y_min  # place leaf markers at axis floor
    for ax, phase in zip(axes, PHASES):
        Z = phase_data[phase]["Z"]
        draw_custom_dendrogram(ax, Z, leaf_order, leaf_color, leaf_y=leaf_y)
        ax.set_yscale("log")
        ax.set_ylim(y_min, y_max * 1.05)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(PHASE_LABEL[phase], fontsize=10)
        n = Z.shape[0] + 1
        ax.set_xlim(-1, n)
    axes[0].set_ylabel(r"merge height (log)")


# --------------------------------------------------------------------------
# Row C: leading-eigenvector amplitude on contact positions (axial)
# --------------------------------------------------------------------------
def panel_row_C(axes, phase_data, anatomy: pd.DataFrame) -> None:
    """Fiedler eigenvector projected onto axial contact scatter.
    Sign-align all phases to task_test then standardize per phase so the
    spatial *pattern* (not the global amplitude) is what's shown. Marker
    size encodes |u_1|; marker color encodes sign(u_1) on a diverging
    scale.
    """
    xs = anatomy.x.values
    ys = anatomy.y.values

    evec_test = phase_data["task_test"]["eig_vecs"][:, 1]
    aligned = {}
    for phase in PHASES:
        v = phase_data[phase]["eig_vecs"][:, 1]
        if np.dot(v, evec_test[:len(v)]) < 0:
            v = -v
        # Standardize per phase: unit L2 norm so the colour scale is
        # comparable across phases despite raw amplitude differences.
        nrm = np.linalg.norm(v)
        if nrm > 0:
            v = v / nrm
        aligned[phase] = v

    vmax = max(np.nanmax(np.abs(v)) for v in aligned.values())
    s_min = 8.0
    s_max = 80.0

    for ax, phase in zip(axes, PHASES):
        evec = aligned[phase]
        n_eff = min(len(xs), len(evec))
        valid = (~np.isnan(xs[:n_eff])) & (~np.isnan(ys[:n_eff]))
        x_v = xs[:n_eff][valid]
        y_v = ys[:n_eff][valid]
        u_v = evec[:n_eff][valid]
        sizes = s_min + (s_max - s_min) * (np.abs(u_v) / max(vmax, 1e-12))
        sc = ax.scatter(
            x_v, y_v, c=u_v, cmap="RdBu_r",
            vmin=-vmax, vmax=vmax,
            s=sizes, edgecolors="black", linewidths=0.4,
            zorder=3,
        )
        ax.axhline(0, color="0.85", lw=0.4, zorder=0)
        ax.axvline(0, color="0.85", lw=0.4, zorder=0)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("L  ←  x  →  R", fontsize=8)
        ax.set_aspect("equal")
        ax.set_title(PHASE_LABEL[phase], fontsize=10)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    div = make_axes_locatable(axes[-1])
    cax = div.append_axes("right", size="5%", pad=0.10)
    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_label(r"$u_1$ (unit-norm, sign-aligned)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    axes[0].set_ylabel(r"P  ←  y  →  A", fontsize=8)


# --------------------------------------------------------------------------
# Row D: communication-distance arcs on brain (top-K ρ̂ pairs) per phase
# --------------------------------------------------------------------------
def panel_row_arcs(axes, phase_data, anatomy: pd.DataFrame, top_k: int = 60) -> None:
    """For each phase, draw the top-K contact pairs by ρ̂(τ) as quadratic
    Bezier arcs over the patient's axial contact scatter.  ρ̂ is the
    inverse communication distance — strong arcs = "tight" communication
    pathways."""
    xs = anatomy.x.values
    ys = anatomy.y.values
    n_eff = min(len(xs), phase_data[PHASES[0]]["R"].shape[0])
    x_v = xs[:n_eff]
    y_v = ys[:n_eff]
    valid_mask = (~np.isnan(x_v)) & (~np.isnan(y_v))

    # Pre-extract top-K pairs per phase using upper-triangular indices.
    iu, ju = np.triu_indices(n_eff, k=1)

    # Use shared color scale across phases on ρ̂ values of the displayed
    # pairs.
    rho_vals_all = []
    for phase in PHASES:
        rho = phase_data[phase]["R"][:n_eff, :n_eff][iu, ju]
        rho_vals_all.append(rho)
    rho_global = np.concatenate(rho_vals_all)
    vmax = float(np.nanpercentile(rho_global, 99.5))
    vmin = float(np.nanpercentile(rho_global, 99.5 - 100 * (top_k / len(iu))))
    if vmax <= vmin:
        vmax = vmin + 1e-3

    cmap = plt.get_cmap("magma_r")

    for ax, phase in zip(axes, PHASES):
        rho = phase_data[phase]["R"][:n_eff, :n_eff][iu, ju]
        # Pick top-K pair indices (highest ρ̂) but exclude NaN-position pairs.
        valid_pairs = valid_mask[iu] & valid_mask[ju]
        rho_filt = rho.copy()
        rho_filt[~valid_pairs] = -np.inf
        order = np.argsort(rho_filt)[::-1][:top_k]
        # Background contacts (light grey)
        ax.scatter(x_v[valid_mask], y_v[valid_mask], color="#dddddd",
                   s=10, alpha=0.7, edgecolors="none", zorder=1)
        # Draw arcs sorted by strength so weaker arcs are drawn first.
        sub_rho = rho[order]
        order_sorted = order[np.argsort(sub_rho)]
        for k in order_sorted:
            i, j = iu[k], ju[k]
            xi, yi = x_v[i], y_v[i]
            xj, yj = x_v[j], y_v[j]
            r = float(rho[k])
            # Color-encode ρ̂ on the magma_r scale.
            col = cmap((r - vmin) / max(vmax - vmin, 1e-9))
            # Quadratic bezier with control point above the midpoint.
            mx = 0.5 * (xi + xj)
            my = 0.5 * (yi + yj)
            dist = float(np.hypot(xj - xi, yj - yi))
            # Bend amount scales with chord length so long arcs aren't
            # crushed and short ones aren't oversized.
            bend = 0.4
            # Orient bend along the perpendicular of the chord.
            if dist > 1e-6:
                px = -(yj - yi) / dist
                py = (xj - xi) / dist
            else:
                px, py = 0, 1
            cx = mx + bend * dist * px
            cy = my + bend * dist * py
            ts = np.linspace(0, 1, 24)
            bx = (1 - ts) ** 2 * xi + 2 * (1 - ts) * ts * cx + ts ** 2 * xj
            by = (1 - ts) ** 2 * yi + 2 * (1 - ts) * ts * cy + ts ** 2 * yj
            ax.plot(bx, by, color=col, lw=0.6,
                    alpha=min(0.95, 0.35 + 0.6 * (r - vmin) / max(vmax - vmin, 1e-9)),
                    zorder=2 + (r > vmin))
        # Highlight contact positions of the displayed pairs.
        endpoints = sorted(set(int(iu[k]) for k in order_sorted)
                           | set(int(ju[k]) for k in order_sorted))
        ax.scatter([x_v[p] for p in endpoints],
                   [y_v[p] for p in endpoints],
                   color="black", s=14, edgecolors="white",
                   linewidths=0.4, zorder=3)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("L  ←  x  →  R", fontsize=8)
        ax.set_aspect("equal")
        ax.set_title(PHASE_LABEL[phase], fontsize=10)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    div = make_axes_locatable(axes[-1])
    cax = div.append_axes("right", size="5%", pad=0.10)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label(rf"$\hat\rho(\tau)$  (top {top_k} pairs)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    axes[0].set_ylabel(r"P  ←  y  →  A", fontsize=8)


# --------------------------------------------------------------------------
# Row E: cohort glass-brain — axial + sagittal scatter overlays
# --------------------------------------------------------------------------
def panel_row_D(axes, anatomy_pool: pd.DataFrame, perpair_leaves: pd.DataFrame,
                kc_leaves: pd.DataFrame) -> None:
    # Two axes: axial (x, y) and sagittal (y, z) — pooled cohort anatomy.
    pool = anatomy_pool.copy()
    pool["is_perpair"] = list(zip(pool.patient, pool.leaf_id))
    pool["is_kc"] = list(zip(pool.patient, pool.leaf_id))
    perpair_set = set(zip(perpair_leaves.patient, perpair_leaves.leaf_id))
    kc_set = set(zip(kc_leaves.patient, kc_leaves.leaf_id))
    pool["is_perpair"] = pool["is_perpair"].isin(perpair_set)
    pool["is_kc"] = pool["is_kc"].isin(kc_set)
    # Drop rows with no anatomy.
    pool = pool.dropna(subset=["x", "y", "z"]).reset_index(drop=True)

    bg = pool[~pool.is_perpair & ~pool.is_kc]
    pp = pool[pool.is_perpair & ~pool.is_kc]
    kc = pool[pool.is_kc & ~pool.is_perpair]
    both = pool[pool.is_perpair & pool.is_kc]

    titles = ["axial (z = MIP)", "sagittal (x = MIP)"]
    coords_xy = [
        ("x", "y"),  # axial
        ("y", "z"),  # sagittal
    ]
    axis_labels = [("L  ←  x  →  R", "P  ←  y  →  A"),
                   ("P  ←  y  →  A", "I  ←  z  →  S")]
    for ax, (cx, cy), title, alab in zip(axes, coords_xy, titles, axis_labels):
        ax.scatter(bg[cx], bg[cy], color=COL_BG, s=8, alpha=0.45,
                   edgecolors="none", zorder=1, label=None)
        ax.scatter(pp[cx], pp[cy], color=COL_PERPAIR, s=22, alpha=0.85,
                   edgecolors="white", linewidths=0.5, zorder=2,
                   label=None)
        ax.scatter(kc[cx], kc[cy], color=COL_HIERARCH, s=24, alpha=0.85,
                   edgecolors="white", linewidths=0.5, zorder=3,
                   label=None)
        ax.scatter(both[cx], both[cy], color="#e74c3c", s=32, alpha=0.95,
                   edgecolors="white", linewidths=0.6, zorder=4,
                   label=None)
        ax.set_xlabel(alab[0], fontsize=8)
        ax.set_ylabel(alab[1], fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.set_aspect("equal")

    handles = [
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=7,
                      markerfacecolor=COL_BG, markeredgecolor="none",
                      alpha=0.6, label="cortical contact"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=7,
                      markerfacecolor=COL_PERPAIR, markeredgecolor="white",
                      label="per-pair (CTM) trace-leaf"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=7,
                      markerfacecolor=COL_HIERARCH, markeredgecolor="white",
                      label=r"hierarchical (KC $\lambda{=}0$) trace-leaf"),
        mlines.Line2D([], [], marker="o", linestyle="None", markersize=7,
                      markerfacecolor="#e74c3c", markeredgecolor="white",
                      label="both layers"),
    ]
    axes[1].legend(handles=handles, loc="lower right",
                   bbox_to_anchor=(1.0, -0.18), ncol=4, frameon=False,
                   fontsize=8)


# --------------------------------------------------------------------------
# Cohort anatomy pool
# --------------------------------------------------------------------------
def load_anatomy_pool() -> pd.DataFrame:
    parts = []
    for p in PATIENTS:
        df = load_anatomy(p)
        df.insert(0, "patient", p)
        parts.append(df)
    return pd.concat(parts, ignore_index=True)


# --------------------------------------------------------------------------
# Composite assembly
# --------------------------------------------------------------------------
def main() -> None:
    # Load exemplar phase data.
    phase_data = {ph: load_lrg(EXEMPLAR, ph, BAND) for ph in PHASES}
    n = phase_data[PHASES[0]]["n"]
    anatomy_pat = load_anatomy(EXEMPLAR).iloc[:n].reset_index(drop=True)

    # Common leaf order: ordered by hierarchical clustering of task_test.
    Z_test = phase_data["task_test"]["Z"]
    leaf_order = leaves_list(Z_test)

    # Cohort-pool anatomy.
    anatomy_pool = load_anatomy_pool()

    # Trace-leaf tables.
    perpair = pd.read_csv(PER_PAIR_LEAVES)
    perpair_beta = perpair[perpair.band == BAND]
    kc = pd.read_csv(KC_TOPOLOGY_LEAVES)
    kc_beta = kc[kc.band == BAND]

    # 5 rows × 4 cols.  Row 5 has 2 wide axes (axial + sagittal cohort).
    fig = plt.figure(figsize=(13.5, 17.0))
    gs = fig.add_gridspec(
        nrows=5, ncols=4,
        height_ratios=[1.0, 0.85, 0.95, 0.95, 1.0],
        wspace=0.25, hspace=0.40,
        left=0.06, right=0.97, top=0.97, bottom=0.04,
    )

    axA = [fig.add_subplot(gs[0, j]) for j in range(4)]
    axB = [fig.add_subplot(gs[1, j]) for j in range(4)]
    axC = [fig.add_subplot(gs[2, j]) for j in range(4)]
    axD = [fig.add_subplot(gs[3, j]) for j in range(4)]
    axE1 = fig.add_subplot(gs[4, 0:2])
    axE2 = fig.add_subplot(gs[4, 2:4])

    panel_row_A(axA, phase_data, leaf_order, anatomy_pat)
    panel_row_B(axB, phase_data, leaf_order, anatomy_pat)
    panel_row_C(axC, phase_data, anatomy_pat)
    panel_row_arcs(axD, phase_data, anatomy_pat)
    panel_row_D([axE1, axE2], anatomy_pool, perpair_beta, kc_beta)

    axA[0].text(-0.32, 0.5, r"A — LRG  $\hat\rho(\tau)$",
                transform=axA[0].transAxes, ha="center", va="center",
                rotation=90, fontsize=11, weight="bold")
    axB[0].text(-0.32, 0.5, "B — dendrogram",
                transform=axB[0].transAxes, ha="center", va="center",
                rotation=90, fontsize=11, weight="bold")
    axC[0].text(-0.32, 0.5, "C — leading eigenmode",
                transform=axC[0].transAxes, ha="center", va="center",
                rotation=90, fontsize=11, weight="bold")
    axD[0].text(-0.32, 0.5, "D — communication arcs",
                transform=axD[0].transAxes, ha="center", va="center",
                rotation=90, fontsize=11, weight="bold")
    axE1.text(-0.18, 0.5, "E — cohort anatomy",
              transform=axE1.transAxes, ha="center", va="center",
              rotation=90, fontsize=11, weight="bold")

    # Region-group legend for row B.
    handles_b = [
        mlines.Line2D([], [], marker="|", linestyle="None", markersize=10,
                      color=COL_TEMPORAL, markeredgewidth=2,
                      label="temporal-lobe leaf"),
        mlines.Line2D([], [], marker="|", linestyle="None", markersize=10,
                      color=COL_OTHER_CORTEX, markeredgewidth=2,
                      label="other-cortex leaf"),
        mlines.Line2D([], [], marker="|", linestyle="None", markersize=10,
                      color=COL_WM_UNK, markeredgewidth=2,
                      label="white matter / unknown"),
    ]
    fig.legend(handles=handles_b, loc="upper center",
               bbox_to_anchor=(0.5, 0.815), ncol=3, frameon=False,
               fontsize=8)

    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig_section5_2_beta_unified] wrote {OUT_PDF}")
    print(f"  exemplar = {EXEMPLAR}, β, {len(PHASES)} phases")
    print(f"  cohort: {len(anatomy_pool)} contacts pooled")
    print(f"  per-pair β trace-leaves: {len(perpair_beta)}")
    print(f"  hierarchical β trace-leaves: {len(kc_beta)}")


if __name__ == "__main__":
    main()
