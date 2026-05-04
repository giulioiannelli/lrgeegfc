"""Polished bundle figures for the writing-agent handoff.

Rebuilds the three headline figures into
`.agents/writing-bundles/raw-fc/` with publication-style visuals
(cohort mean + 1σ ellipse for the trace scatter, paired strip+IQR
box for the phase geometry, numbered markers + trace-zone shading
for the d_S × d_P convergence). Uses TRACE terminology throughout
(see .agents/guides/01_project/terminology.md).

Outputs (overwrite):
- fig1_trace_scatter_dS{SUFFIX}.pdf
- fig2_phase_geometry{SUFFIX}.pdf
- fig3_dS_dP_convergence{SUFFIX}.pdf
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Patch, Rectangle, Arc
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from scipy.stats import pearsonr, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    PATIENTS_4PHASE,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    NODE_SIZE_DEFAULT,
)
from lrg_eegfc.config.paths import DATA_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.probe import extract_probe_labels
from lrgsglib.plotlib import imshow_colorbar_caxdivider

# Network helpers — promoted location TBD; archived copy still authoritative
sys.path.insert(0, str(ROOT / "scripts/archive/2026-02_imcoh-dev-notes"))
from _shared import (  # noqa: E402
    load_channel_labels,
    compute_network_layout,
    draw_network_edges,
    _probe_color_map,
)

PHASE_SHORT = {
    "rest_pre":   "RPre",
    "task_learn": "TL",
    "task_test":  "TT",
    "rest_post":  "RPost",
}

N_COHORT = len(PATIENTS_4PHASE)

# Substrate selection — CLI arg picks `imcoh_abs` (default) or `imcoh_sq`.
# Output filenames get a suffix so both substrates can coexist in the bundle.
SUBSTRATE = sys.argv[1] if len(sys.argv) > 1 else "imcoh_abs"
assert SUBSTRATE in ("imcoh_abs", "imcoh_sq"), \
    f"unknown substrate {SUBSTRATE}"
if SUBSTRATE == "imcoh_abs":
    SRC = DATA_ROOT / "audit" / "raw_fc_phase_distance"
else:
    SRC = DATA_ROOT / "audit" / "raw_fc_phase_distance_imcoh_sq"
SUFFIX = f"_{SUBSTRATE}"
FC_METHOD_FIG5 = SUBSTRATE   # used by load_fc_matrix in fig 5
OUT = ROOT / ".agents" / "writing-bundles" / "raw-fc"
OUT.mkdir(parents=True, exist_ok=True)
print(f"=== bundle figures for substrate: {SUBSTRATE} ===")
print(f"    SRC = {SRC}")
print(f"    OUT suffix = {SUFFIX}")

# ----------------------------------------------------------------------
# Phase-pair semantic palette (consistent across figures)
# ----------------------------------------------------------------------
PAIR_COLORS = {
    ("task_learn", "task_test"): "#1f77b4",  # blue   — within-task cohesion
    ("task_test", "rest_post"):  "#d62728",  # red    — THE trace pair
    ("rest_pre",  "task_test"):  "#ff7f0e",  # orange — pre→task baseline shift
    ("rest_pre",  "rest_post"):  "#7f7f7f",  # gray   — within-rest drift floor
}
PAIR_LABELS = {
    ("task_learn", "task_test"): "TL$\\leftrightarrow$TT",
    ("task_test", "rest_post"):  "TT$\\leftrightarrow$RPost",
    ("rest_pre",  "task_test"):  "RPre$\\leftrightarrow$TT",
    ("rest_pre",  "rest_post"):  "RPre$\\leftrightarrow$RPost",
}
# Sorted small → large by typical cohort median: task-tight, trace, baseline shift, drift floor
PAIR_ORDER = [
    ("task_learn", "task_test"),
    ("task_test",  "rest_post"),
    ("rest_pre",   "task_test"),
    ("rest_pre",   "rest_post"),
]


def lookup_pair(df, band, A, B, distance="S"):
    """Return per-patient `d_obs` Series indexed by patient for a phase pair."""
    sub = df[
        (df.distance == distance)
        & (df.band == band)
        & (((df.phase_A == A) & (df.phase_B == B))
           | ((df.phase_A == B) & (df.phase_B == A)))
    ]
    return sub.groupby("patient").d_obs.first()


# ----------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------
distances_long = []
for p in PATIENTS_4PHASE:
    f = SRC / p / "distance_4phase.csv"
    if not f.exists():
        print(f"[{p}] missing distance_4phase.csv, skipping")
        continue
    distances_long.append(pd.read_csv(f))
distances_long = pd.concat(distances_long, ignore_index=True)

td = pd.read_csv(SRC / "Td_per_patient_per_band.csv")
contrast = pd.read_csv(SRC / "Td_dS_vs_dF_band_contrast.csv").set_index("band")


# ----------------------------------------------------------------------
# Fig 1 — Trace scatter on d_S
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 6, figsize=(15.6, 2.95), squeeze=False)
sub_S = distances_long[distances_long.distance == "S"]

for c, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes[0, c]
    x_s = lookup_pair(sub_S, band, "rest_pre",  "task_test", "S")
    y_s = lookup_pair(sub_S, band, "task_test", "rest_post", "S")
    common = sorted(set(x_s.index) & set(y_s.index))
    if not common:
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
        continue
    x = x_s.loc[common].values
    y = y_s.loc[common].values

    span = max(x.max(), y.max()) - min(x.min(), y.min())
    lo = min(x.min(), y.min()) - 0.05 * span
    hi = max(x.max(), y.max()) + 0.05 * span

    # Trace zone (below identity) — light red shade
    poly_x = [lo, hi, hi]
    poly_y = [lo, lo, hi]
    ax.fill(poly_x, poly_y, color="#d62728", alpha=0.06, zorder=0)

    # Identity line
    ax.plot([lo, hi], [lo, hi], color="k", lw=0.7, ls="--",
            alpha=0.55, zorder=1)

    # Cohort 1σ covariance ellipse
    mu_x, mu_y = float(x.mean()), float(y.mean())
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2.0 * np.sqrt(np.maximum(vals, 0.0))
    if width > 0 and height > 0:
        ell = Ellipse(
            (mu_x, mu_y), width, height, angle=angle,
            facecolor="#888", alpha=0.18, edgecolor="#444",
            linewidth=0.7, zorder=2,
        )
        ax.add_patch(ell)

    # Per-patient circles (Pat_03 distinct)
    pat_codes = np.asarray(common)
    is_p03 = pat_codes == "Pat_03"
    ax.scatter(x[~is_p03], y[~is_p03], c="white",
               edgecolor="#1f77b4", linewidth=1.1, s=34, zorder=3)
    if is_p03.any():
        ax.scatter(x[is_p03], y[is_p03], marker="^", c="white",
                   edgecolor="#ff7f0e", linewidth=1.2, s=46, zorder=4)

    # Cohort mean star (atop everything)
    ax.scatter([mu_x], [mu_y], marker="*", s=180, c="#2ca02c",
               edgecolor="k", linewidth=0.9, zorder=5)

    # n_trace annotation
    n_trace = int((y < x).sum())
    n_total = len(common)
    ax.text(0.05, 0.95, f"trace: {n_trace}/{n_total}",
            transform=ax.transAxes, fontsize=8.5, va="top",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#999", alpha=0.92))

    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
    ax.set_xlabel(r"$d_S$(RPre, TT)", fontsize=8.5)
    if c == 0:
        ax.set_ylabel(r"$d_S$(TT, RPost)", fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(color="#eee", lw=0.4, zorder=0)
    ax.set_axisbelow(True)

handles_f1 = [
    Line2D([0], [0], marker="*", color="none",
           markerfacecolor="#2ca02c", markeredgecolor="k",
           markersize=12, label="cohort mean"),
    Patch(facecolor="#888", alpha=0.18, edgecolor="#444",
          label=r"cohort 1$\sigma$ ellipse"),
    Patch(facecolor="#d62728", alpha=0.10,
          label="trace zone (below identity)"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor="#1f77b4",
           markersize=8, label=f"patient (n={N_COHORT - 1})"),
    Line2D([0], [0], marker="^", color="none",
           markerfacecolor="white", markeredgecolor="#ff7f0e",
           markersize=9, label="Pat_03 (1024 Hz)"),
]
fig.legend(handles=handles_f1, loc="lower center",
           bbox_to_anchor=(0.5, -0.03), ncol=len(handles_f1),
           frameon=False, fontsize=8.2)
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(OUT / f"fig1_trace_scatter_dS{SUFFIX}.pdf", bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT / f'fig1_trace_scatter_dS{SUFFIX}.pdf'}")


# ----------------------------------------------------------------------
# Fig 2 — Phase geometry (paired strip + IQR box, sorted phase pairs)
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 6, figsize=(15.6, 3.6), squeeze=False)

for c, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes[0, c]
    pair_data = {pr: lookup_pair(sub_S, band, pr[0], pr[1], "S")
                 for pr in PAIR_ORDER}
    common = sorted(set.intersection(*[set(d.index) for d in pair_data.values()]))
    # Per-band sort: small → large by cohort median
    band_order = sorted(
        PAIR_ORDER,
        key=lambda pr: float(np.median(pair_data[pr].loc[common].values)),
    )
    xs = np.arange(len(band_order))

    # Per-patient connecting lines (in band-sorted order)
    for p in common:
        ys = [pair_data[pr][p] for pr in band_order]
        ax.plot(xs, ys, color="#bbb", lw=0.45, alpha=0.55, zorder=1)

    # Per-pair: IQR box + median bar + strip dots
    for i, pr in enumerate(band_order):
        vals = pair_data[pr].loc[common].values
        col = PAIR_COLORS[pr]
        q1, med, q3 = np.percentile(vals, [25, 50, 75])

        # IQR box (translucent fill)
        ax.add_patch(Rectangle(
            (i - 0.20, q1), 0.40, q3 - q1,
            facecolor=col, alpha=0.20, edgecolor=col,
            linewidth=0.7, zorder=2,
        ))
        # Median line
        ax.plot([i - 0.24, i + 0.24], [med, med],
                color=col, lw=2.2, zorder=4)

        # Strip dots (deterministic jitter)
        rs = np.random.RandomState(42 + i + c * 13)
        jitter = rs.uniform(-0.11, 0.11, size=len(vals))
        is_p03 = np.asarray(common) == "Pat_03"
        ax.scatter(np.full_like(vals, i)[~is_p03] + jitter[~is_p03],
                   vals[~is_p03], c="white",
                   edgecolor=col, linewidth=1.0, s=26,
                   alpha=0.95, zorder=3)
        if is_p03.any():
            ax.scatter(np.full_like(vals, i)[is_p03] + jitter[is_p03],
                       vals[is_p03], marker="^", c="white",
                       edgecolor=col, linewidth=1.0, s=34,
                       alpha=0.95, zorder=3)

    ax.set_xticks(xs)
    ax.set_xticklabels([PAIR_LABELS[pr] for pr in band_order],
                       fontsize=7.2, rotation=22, ha="right")
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
    ax.tick_params(axis="y", labelsize=7)
    if c == 0:
        ax.set_ylabel(rf"$d_S$  (n={N_COHORT} patients)", fontsize=9)
    ax.grid(axis="y", color="#eee", lw=0.5, zorder=0)
    ax.set_axisbelow(True)

handles_f2 = [Patch(facecolor=PAIR_COLORS[pr], alpha=0.55,
                    edgecolor=PAIR_COLORS[pr],
                    label=PAIR_LABELS[pr])
              for pr in PAIR_ORDER]
handles_f2 += [
    Line2D([0], [0], color="#bbb", lw=0.7, label="patient (paired)"),
    Line2D([0], [0], color="k", lw=2.0, label="cohort median"),
]
fig.legend(handles=handles_f2, loc="lower center",
           bbox_to_anchor=(0.5, -0.02), ncol=len(handles_f2),
           frameon=False, fontsize=8.2)
fig.tight_layout(rect=[0, 0.07, 1, 1])
fig.savefig(OUT / f"fig2_phase_geometry{SUFFIX}.pdf", bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT / f'fig2_phase_geometry{SUFFIX}.pdf'}")


# ----------------------------------------------------------------------
# Fig 3 — d_S × d_P convergence (numbered markers, trace-zone shading)
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 6, figsize=(15.6, 3.0), squeeze=False)

for c, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes[0, c]
    sub = td[td.band == band].dropna(subset=["S", "P"])
    if sub.empty:
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11); continue
    xs_v = sub.S.values
    ys_v = sub.P.values
    pats = sub.patient.values

    pad_pct = 0.10
    xmin, xmax = float(xs_v.min()), float(xs_v.max())
    ymin, ymax = float(ys_v.min()), float(ys_v.max())
    lo = min(xmin, ymin) - pad_pct * max(abs(xmin), abs(xmax),
                                         abs(ymin), abs(ymax), 0.05)
    hi = max(xmax, ymax) + pad_pct * max(abs(xmin), abs(xmax),
                                         abs(ymin), abs(ymax), 0.05)
    # Symmetrize a bit to keep zero centered-ish
    extent = max(abs(lo), abs(hi))
    lo, hi = -extent, extent

    # Trace zone (lower-left quadrant)
    ax.add_patch(Rectangle((lo, lo), -lo, -lo,
                           facecolor="#2ca02c", alpha=0.07, zorder=0))
    # Crosshair at origin
    ax.axhline(0, color="#999", lw=0.6, ls="--", zorder=1)
    ax.axvline(0, color="#999", lw=0.6, ls="--", zorder=1)
    # Identity line
    ax.plot([lo, hi], [lo, hi], color="#444", lw=0.7, ls=":",
            alpha=0.7, zorder=1)

    # Cohort mean (star)
    mu_x, mu_y = float(xs_v.mean()), float(ys_v.mean())
    ax.scatter([mu_x], [mu_y], marker="*", s=240, c="#1f77b4",
               edgecolor="k", linewidth=0.9, zorder=5, alpha=0.95)

    # Numbered patient circles
    for x, y, p in zip(xs_v, ys_v, pats):
        nid = p.replace("Pat_", "")
        in_trace = (x < 0) and (y < 0)
        is_p03 = (p == "Pat_03")
        if is_p03:
            edge = "#ff7f0e"
        elif in_trace:
            edge = "#1b5e20"
        else:
            edge = "#666"
        ax.scatter([x], [y], marker="o", s=320, c="white",
                   edgecolor=edge, linewidth=1.1, zorder=3)
        ax.text(x, y, nid, ha="center", va="center",
                fontsize=7.2, fontweight="bold", color=edge,
                zorder=4)

    n_trace = int(((xs_v < 0) & (ys_v < 0)).sum())
    n_total = len(xs_v)
    ax.text(0.05, 0.95, f"trace: {n_trace}/{n_total}",
            transform=ax.transAxes, fontsize=8.5, va="top",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#999", alpha=0.92))

    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
    ax.set_xlabel(r"$T_d^{(d_S)}$  rank-only", fontsize=8.5)
    if c == 0:
        ax.set_ylabel(r"$T_d^{(d_P)}$  mag-weighted", fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(color="#eee", lw=0.4, zorder=0)
    ax.set_axisbelow(True)

handles_f3 = [
    Line2D([0], [0], marker="*", color="none",
           markerfacecolor="#1f77b4", markeredgecolor="k",
           markersize=14, label="cohort mean"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor="#1b5e20",
           markersize=12, label="patient (in trace zone)"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor="#666",
           markersize=12, label="patient (off zone)"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor="#ff7f0e",
           markersize=12, label="Pat_03 (1024 Hz)"),
    Patch(facecolor="#2ca02c", alpha=0.13,
          label=r"trace zone ($T_d^{(d_S)}<0$ AND $T_d^{(d_P)}<0$)"),
]
fig.legend(handles=handles_f3, loc="lower center",
           bbox_to_anchor=(0.5, -0.02), ncol=len(handles_f3),
           frameon=False, fontsize=8.0)
fig.tight_layout(rect=[0, 0.07, 1, 1])
fig.savefig(OUT / f"fig3_dS_dP_convergence{SUFFIX}.pdf", bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT / f'fig3_dS_dP_convergence{SUFFIX}.pdf'}")


# ----------------------------------------------------------------------
# Fig 4 — Structural vs drift (T_d^(d_S) vs T_d^(d_F))
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 6, figsize=(15.6, 3.0), squeeze=False)

for c, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes[0, c]
    sub = td[td.band == band].dropna(subset=["S", "F"])
    if sub.empty:
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11); continue
    xs_v = sub.S.values
    ys_v = sub.F.values
    pats = sub.patient.values

    # Symmetric square axes
    extent = max(np.abs(xs_v).max(), np.abs(ys_v).max())
    extent = extent * 1.15  # padding
    lo, hi = -extent, extent

    # Trace zone (lower-left quadrant, both T_d < 0)
    ax.add_patch(Rectangle((lo, lo), -lo, -lo,
                           facecolor="#2ca02c", alpha=0.07, zorder=0))
    # Crosshair at origin
    ax.axhline(0, color="#999", lw=0.6, ls="--", zorder=1)
    ax.axvline(0, color="#999", lw=0.6, ls="--", zorder=1)
    # Identity line — dotted
    ax.plot([lo, hi], [lo, hi], color="#444", lw=0.8, ls=":",
            alpha=0.75, zorder=1)

    # Cohort 1σ covariance ellipse — shows the elongation along
    # identity == strong rank/amplitude correlation
    mu_x, mu_y = float(xs_v.mean()), float(ys_v.mean())
    cov = np.cov(xs_v, ys_v)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2.0 * np.sqrt(np.maximum(vals, 0.0))
    if width > 0 and height > 0:
        ell = Ellipse(
            (mu_x, mu_y), width, height, angle=angle,
            facecolor="#888", alpha=0.18, edgecolor="#444",
            linewidth=0.7, zorder=2,
        )
        ax.add_patch(ell)

    # Cohort mean star (atop ellipse)
    ax.scatter([mu_x], [mu_y], marker="*", s=240, c="#1f77b4",
               edgecolor="k", linewidth=0.9, zorder=5, alpha=0.95)

    # Numbered patient circles (last 2 digits)
    for x, y, p in zip(xs_v, ys_v, pats):
        nid = p.replace("Pat_", "")
        in_trace = (x < 0) and (y < 0)
        is_p03 = (p == "Pat_03")
        if is_p03:
            edge = "#ff7f0e"
        elif in_trace:
            edge = "#1b5e20"
        else:
            edge = "#666"
        ax.scatter([x], [y], marker="o", s=320, c="white",
                   edgecolor=edge, linewidth=1.1, zorder=3)
        ax.text(x, y, nid, ha="center", va="center",
                fontsize=7.2, fontweight="bold", color=edge,
                zorder=4)

    # Inset: per-band Spearman ρ + sign agreement (the load-bearing
    # numbers for the "this is not drift" call)
    rho = float(contrast.loc[band, "spearman_rho_S_F"])
    sign_str = str(contrast.loc[band, "sign_agree_S_F"])
    n_trace = int(((xs_v < 0) & (ys_v < 0)).sum())
    n_total = len(xs_v)
    ax.text(0.05, 0.95,
            (f"trace: {n_trace}/{n_total}\n"
             rf"$\rho_S$ = {rho:+.2f}" + "\n"
             f"sign: {sign_str}"),
            transform=ax.transAxes, fontsize=7.5, va="top",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#999", alpha=0.92))

    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
    ax.set_xlabel(r"$T_d^{(d_S)}$  rank-only", fontsize=8.5)
    if c == 0:
        ax.set_ylabel(r"$T_d^{(d_F)}$  amplitude-only", fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(color="#eee", lw=0.4, zorder=0)
    ax.set_axisbelow(True)

handles_f4 = [
    Line2D([0], [0], marker="*", color="none",
           markerfacecolor="#1f77b4", markeredgecolor="k",
           markersize=14, label="cohort mean"),
    Patch(facecolor="#888", alpha=0.18, edgecolor="#444",
          label=r"cohort 1$\sigma$ ellipse (elongation $\sim$ correlation)"),
    Line2D([0], [0], color="#444", lw=0.9, ls=":",
           label="identity (perfect rank/amp agreement)"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor="#1b5e20",
           markersize=12, label="patient (in trace zone)"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor="#666",
           markersize=12, label="patient (off zone)"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor="#ff7f0e",
           markersize=12, label="Pat_03 (1024 Hz)"),
]
fig.legend(handles=handles_f4, loc="lower center",
           bbox_to_anchor=(0.5, -0.02), ncol=3,
           frameon=False, fontsize=8.0)
fig.tight_layout(rect=[0, 0.10, 1, 1])
fig.savefig(OUT / f"fig4_structural_vs_drift{SUFFIX}.pdf", bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT / f'fig4_structural_vs_drift{SUFFIX}.pdf'}")


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


# ----------------------------------------------------------------------
# Fig 6 — 4-phase geometry as chord/arc diagrams (no redundant matrix)
# ----------------------------------------------------------------------
# A 4-phase distance "matrix" has 4×4 = 16 cells but only 6 unique pair
# distances (4 choose 2). The diagonal is trivial (=0) and the lower
# triangle is symmetric to the upper. Replacing the matrix with a
# chord diagram per (band, distance) cell:
#   - 4 phase markers on a baseline preserve the temporal order
#     RPre → TL → TT → RPost
#   - 6 arcs above connect each pair
#   - arc COLOUR encodes the cohort-median distance value
#   - arc HEIGHT encodes the temporal gap (adjacent / skip-1 / skip-2)
geom = pd.read_csv(SRC / "cohort_geometry_4phase_summary.csv")

PHASES_ORDER = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_X = {ph: i for i, ph in enumerate(PHASES_ORDER)}
PHASE_LABEL_SHORT = {"rest_pre": "RPre", "task_learn": "TL",
                     "task_test": "TT", "rest_post": "RPost"}
DISTANCE_KEYS = ["S", "P", "F"]
DISTANCE_LABEL = {
    "S": r"$d_S$  Spearman",
    "P": r"$d_P$  Pearson",
    "F": r"$d_F$  Frobenius",
}

fig = plt.figure(figsize=(15.4, 5.8))
gs = fig.add_gridspec(
    3, 6,
    width_ratios=[1, 1, 1, 1, 1, 1],
    height_ratios=[1, 1, 1],
    wspace=0.32, hspace=0.32,
    top=0.91, bottom=0.10, left=0.05, right=0.985,
)
REF_LINES = [0.0, 0.25, 0.5, 0.75, 1.0]

# Custom 2-colour colormap: vivid green (close) → neutral gray midpoint
# → vivid red (far). NO yellow waypoint, so the perceived semantic is
# binary "good/bad" with a clean neutral middle.
cmap_chord = LinearSegmentedColormap.from_list(
    "close_gray_far",
    [(0.00, "#1a9850"),   # vivid green = close
     (0.50, "#9e9e9e"),   # neutral gray midpoint (no yellow)
     (1.00, "#d73027")],  # vivid red   = far
    N=256,
)


def _solve_1d_spring_layout(panel_distances, n_phases=4):
    """Weighted 1D MDS: place `n_phases` phase points on a line so the
    pairwise distances |p_i - p_j| match `panel_distances` as closely
    as possible, with **spring stiffness w_ij = 1 / d_ij**² (closer
    pairs pull MUCH harder than far ones — amplifies the difference
    between near and distant pairs). Monotonic order p_0 < p_1 < ...
    < p_{n-1} is enforced via a positive-delta parameterisation.

    panel_distances: list of (i, j, d_ij) with i < j.
    Returns: array of n_phases positions, p_0 = 0.
    """
    def _stress(deltas):
        steps = np.abs(deltas)
        p = np.concatenate([[0.0], np.cumsum(steps)])
        s = 0.0
        for i, j, d in panel_distances:
            w = 1.0 / max(d, 1e-3) ** 2     # 1/d^2 stiffness
            s += w * (abs(p[i] - p[j]) - d) ** 2
        return s

    mean_d = np.mean([d for _, _, d in panel_distances])
    x0 = np.full(n_phases - 1, mean_d)
    result = _spo_minimize(
        _stress, x0, method="Nelder-Mead",
        options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 8000},
    )
    deltas = np.abs(result.x)
    return np.concatenate([[0.0], np.cumsum(deltas)])


from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable

for r_idx, dist in enumerate(DISTANCE_KEYS):
    sub = geom[geom.distance == dist]

    for c_idx, band in enumerate(BRAIN_BANDS_NAMES):
        ax = fig.add_subplot(gs[r_idx, c_idx])
        panel = sub[sub.band == band]

        panel_distances = []
        for _, row in panel.iterrows():
            i = PHASE_X[row.phase_A]
            j = PHASE_X[row.phase_B]
            if i > j:
                i, j = j, i
            panel_distances.append((i, j, float(row["median"])))

        # Per-PANEL colour normalisation (relative within panel —
        # distances are not directly comparable across cells).
        d_values = np.array([d for _, _, d in panel_distances])
        vmin = float(d_values.min())
        vmax = float(d_values.max())
        if vmax <= vmin:
            vmax = vmin + 1e-6
        norm_d = Normalize(vmin=vmin, vmax=vmax)

        # 1-D spring layout (stiffness 1/d^2) — rescaled to [0, 1].
        raw_pos = _solve_1d_spring_layout(panel_distances, n_phases=4)
        if raw_pos.max() > 0:
            positions = raw_pos / raw_pos.max()
        else:
            positions = np.linspace(0, 1, 4)

        # Reference dashed grid at 0, .25, .5, .75, 1
        for ref_x in REF_LINES:
            ax.axvline(ref_x, color="#cccccc", lw=0.45, ls=(0, (2, 2)),
                       alpha=0.7, zorder=0)

        # Arcs
        for i, j, d in panel_distances:
            x_i, x_j = positions[i], positions[j]
            x_lo, x_hi = sorted([x_i, x_j])
            arc_w = x_hi - x_lo
            arc_h = max(0.55 * arc_w, 0.05)
            colour = cmap_chord(float(norm_d(d)))
            arc = Arc(
                ((x_lo + x_hi) / 2, 0),
                width=arc_w, height=2 * arc_h,
                theta1=0, theta2=180,
                edgecolor=colour, linewidth=2.2,
                zorder=3,
            )
            ax.add_patch(arc)
            ax.text(
                (x_lo + x_hi) / 2, arc_h + 0.015,
                f"{d:.2f}",
                ha="center", va="bottom",
                fontsize=5.5, color=colour, fontweight="bold",
                zorder=4,
            )

        # Phase markers + labels at their spring-equilibrated positions
        for px_i, ph in enumerate(PHASES_ORDER):
            ax.scatter([positions[px_i]], [0], s=38, c="white",
                       edgecolor="k", linewidth=0.9, zorder=5)
        if r_idx == len(DISTANCE_KEYS) - 1:
            for px_i, ph in enumerate(PHASES_ORDER):
                ax.text(positions[px_i], -0.10,
                        PHASE_LABEL_SHORT[ph],
                        ha="center", va="top", fontsize=7)

        # Per-row label (rotated) on the leftmost panel only
        if c_idx == 0:
            ax.text(
                -0.32, 0.35, DISTANCE_LABEL[dist],
                rotation=90, ha="right", va="center",
                fontsize=10, fontweight="bold",
                transform=ax.transData,
            )

        if r_idx == 0:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10, pad=4)

        ax.set_xlim(-0.10, 1.10)
        ax.set_ylim(-0.18, 0.78)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Per-PANEL colorbar — small, just shows vmin/vmax for that
        # cell so each panel's colour gradient is interpretable on its
        # own scale.
        divider = _make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.04)
        sm = plt.cm.ScalarMappable(cmap=cmap_chord, norm=norm_d)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels([f"{vmin:.2f}", f"{vmax:.2f}"])
        cbar.minorticks_off()
        cbar.ax.tick_params(labelsize=5.5, pad=0.5, length=1.5)

# Footer legend explaining the encoding
handles_f6 = [
    Line2D([0], [0], color="#1a9850", lw=2.6,
           label="arc COLOUR = within-panel distance (green=close, red=far)"),
    Line2D([0], [0], color="#888", lw=2.6,
           label=r"phase POSITIONS = 1-D spring layout (stiffness = $1/d^{2}$)"),
    Line2D([0], [0], color="#bbb", lw=0.8, ls="--",
           label="reference grid at 0, .25, .5, .75, 1"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="white", markeredgecolor="k",
           markersize=8, label="phase node"),
]
fig.legend(handles=handles_f6, loc="lower center",
           bbox_to_anchor=(0.5, 0.005), ncol=len(handles_f6),
           frameon=False, fontsize=8.0)

fig.savefig(OUT / f"fig6_4phase_geometry{SUFFIX}.pdf", bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT / f'fig6_4phase_geometry{SUFFIX}.pdf'}")


# ----------------------------------------------------------------------
# Fig S1 — Z inflation diagnostic (only when null artefacts exist;
# audit_25 within-RPre split-half null is currently only run on
# `imcoh_abs`, so figS1 is skipped for `imcoh_sq`).
# ----------------------------------------------------------------------
_per_patient_csv = SRC / "per_patient_with_scale.csv"
if not _per_patient_csv.exists():
    print(f"skipping figS1 (no {_per_patient_csv.name} for {SUBSTRATE})")
else:
    per_patient = pd.read_csv(_per_patient_csv)

    PAIR_S1 = [
        ("rest_pre", "task_test"),
        ("rest_pre", "rest_post"),
        ("task_test", "rest_post"),
    ]
    PAIR_S1_COLOUR = {
        ("rest_pre", "task_test"):  "#d62728",   # red
        ("rest_pre", "rest_post"):  "#1f77b4",   # blue
        ("task_test", "rest_post"): "#2ca02c",   # green
    }
    PAIR_S1_LABEL = {
        ("rest_pre", "task_test"):  r"$d_{\mathrm{obs}}$  RPre$\rightarrow$TT",
        ("rest_pre", "rest_post"):  r"$d_{\mathrm{obs}}$  RPre$\rightarrow$RPost",
        ("task_test", "rest_post"): r"$d_{\mathrm{obs}}$  TT$\rightarrow$RPost",
    }

    fig, axes = plt.subplots(3, 6, figsize=(15.6, 7.0), squeeze=False)
    fig.subplots_adjust(top=0.93, bottom=0.10, left=0.06, right=0.985,
                        wspace=0.18, hspace=0.32)

    for r_idx, dist in enumerate(DISTANCE_KEYS):
        for c_idx, band in enumerate(BRAIN_BANDS_NAMES):
            ax = axes[r_idx, c_idx]
            sub = per_patient[(per_patient.distance == dist)
                              & (per_patient.band == band)]
            for p_idx, p in enumerate(PATIENTS_4PHASE):
                pdata = sub[sub.patient == p]
                if pdata.empty:
                    continue
                null_q1 = float(pdata.iloc[0].null_q1)
                null_q3 = float(pdata.iloc[0].null_q3)
                null_median = float(pdata.iloc[0].null_median)
                ax.add_patch(Rectangle(
                    (null_q1, p_idx - 0.32), null_q3 - null_q1, 0.64,
                    facecolor="#cccccc", edgecolor="#888",
                    linewidth=0.4, alpha=0.75, zorder=1,
                ))
                ax.plot([null_median, null_median],
                        [p_idx - 0.38, p_idx + 0.38],
                        color="k", lw=0.7, zorder=2)
                for _, row in pdata.iterrows():
                    pair = (row.phase_A, row.phase_B)
                    ax.scatter([row.d_obs], [p_idx], s=22,
                               c=PAIR_S1_COLOUR.get(pair, "#444"),
                               edgecolor="k", linewidth=0.35,
                               zorder=3)
            ax.set_yticks(np.arange(len(PATIENTS_4PHASE)))
            if c_idx == 0:
                ax.set_yticklabels([p.replace("Pat_", "")
                                    for p in PATIENTS_4PHASE],
                                   fontsize=6)
            else:
                ax.set_yticklabels([])
            ax.set_ylim(-0.6, len(PATIENTS_4PHASE) - 0.4)
            ax.invert_yaxis()
            ax.tick_params(labelsize=6.5)
            if r_idx == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
            ax.set_xlabel(rf"$d_{dist}$", fontsize=8, labelpad=2)
            if c_idx == 0:
                ax.text(-0.40, 0.5, DISTANCE_LABEL[dist],
                        rotation=90, ha="right", va="center",
                        fontsize=10, fontweight="bold",
                        transform=ax.transAxes)
            ax.grid(axis="x", color="#eee", lw=0.4, zorder=0)
            ax.set_axisbelow(True)

    handles_s1 = [
        Patch(facecolor="#cccccc", edgecolor="#888",
              label=r"within-RPre null  Q1$-$Q3 (50 splits)"),
        Line2D([0], [0], color="k", lw=0.9, label="null median"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=PAIR_S1_COLOUR[("rest_pre", "task_test")],
               markeredgecolor="k", markersize=7,
               label=PAIR_S1_LABEL[("rest_pre", "task_test")]),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=PAIR_S1_COLOUR[("rest_pre", "rest_post")],
               markeredgecolor="k", markersize=7,
               label=PAIR_S1_LABEL[("rest_pre", "rest_post")]),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=PAIR_S1_COLOUR[("task_test", "rest_post")],
               markeredgecolor="k", markersize=7,
               label=PAIR_S1_LABEL[("task_test", "rest_post")]),
    ]
    fig.legend(handles=handles_s1, loc="lower center",
               bbox_to_anchor=(0.5, 0.015), ncol=len(handles_s1),
               frameon=False, fontsize=8.0)

    fig.savefig(OUT / f"figS1_z_inflation{SUFFIX}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'figS1_z_inflation{SUFFIX}.pdf'}")


# ----------------------------------------------------------------------
# Fig S2 — T_d swarm per band, INDEPENDENT y-axes
# ----------------------------------------------------------------------
td_swarm = pd.read_csv(SRC / "Td_per_patient_per_band.csv")

fig, axes = plt.subplots(1, 6, figsize=(15.4, 3.6), squeeze=False)
fig.subplots_adjust(top=0.88, bottom=0.18, left=0.05, right=0.985,
                    wspace=0.30)

for c_idx, band in enumerate(BRAIN_BANDS_NAMES):
    ax = axes[0, c_idx]
    sub = td_swarm[td_swarm.band == band]
    pat_arr = sub.patient.values
    is_p03 = pat_arr == "Pat_03"
    for x_pos, dist in enumerate(["S", "P", "F"]):
        vals = sub[dist].values
        rs = np.random.RandomState(11 + c_idx * 7 + x_pos * 13)
        jitter = rs.uniform(-0.13, 0.13, size=len(vals))
        # In-pool patients (not Pat_03)
        in_pool_mask = ~is_p03
        pool_vals = vals[in_pool_mask]
        pool_x = np.full(in_pool_mask.sum(), x_pos) + jitter[in_pool_mask]
        ax.scatter(
            pool_x, pool_vals,
            c=["#2ca02c" if v < 0 else "#d62728" for v in pool_vals],
            edgecolor="k", linewidth=0.4, s=42, zorder=3, alpha=0.95,
        )
        # Pat_03 — orange triangle
        if is_p03.any():
            p03_vals = vals[is_p03]
            p03_x = np.full(is_p03.sum(), x_pos) + jitter[is_p03]
            ax.scatter(p03_x, p03_vals, marker="^", c="#ff7f0e",
                       edgecolor="k", linewidth=0.4, s=58, zorder=4)
        # Cohort median bar
        med = float(np.median(vals))
        ax.plot([x_pos - 0.27, x_pos + 0.27], [med, med],
                color="k", lw=2.2, zorder=5)
        # n_neg / n_tot annotation, anchored just inside the top
        n_neg = int((vals < 0).sum())
        n_tot = len(vals)
        ax.annotate(f"{n_neg}/{n_tot}",
                    xy=(x_pos, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(0, -2), textcoords="offset points",
                    ha="center", va="top",
                    fontsize=7.2, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.18",
                              facecolor="white", edgecolor="#bbb",
                              alpha=0.92))
    ax.axhline(0, color="#999", lw=0.6, ls="--", zorder=1)
    ax.set_xlim(-0.55, 2.55)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([r"$d_S$", r"$d_P$", r"$d_F$"], fontsize=9)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10, pad=6)
    if c_idx == 0:
        ax.set_ylabel(r"$T_d$  (negative $=$ trace)",
                      fontsize=9, labelpad=2)
    ax.tick_params(labelsize=7)
    ax.grid(axis="y", color="#eee", lw=0.5, zorder=0)
    ax.set_axisbelow(True)

handles_s2 = [
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="#2ca02c", markeredgecolor="k",
           markersize=8, label=r"patient with $T_d<0$ (trace)"),
    Line2D([0], [0], marker="o", color="none",
           markerfacecolor="#d62728", markeredgecolor="k",
           markersize=8, label=r"patient with $T_d>0$"),
    Line2D([0], [0], marker="^", color="none",
           markerfacecolor="#ff7f0e", markeredgecolor="k",
           markersize=9, label="Pat_03 (1024 Hz)"),
    Line2D([0], [0], color="k", lw=2.2, label="cohort median"),
    Line2D([0], [0], color="#999", lw=0.8, ls="--", label="zero line"),
]
fig.legend(handles=handles_s2, loc="lower center",
           bbox_to_anchor=(0.5, 0.005), ncol=len(handles_s2),
           frameon=False, fontsize=8.0)

fig.savefig(OUT / f"figS2_Td_swarm_all_distances{SUFFIX}.pdf",
            bbox_inches="tight")
plt.close(fig)
print(f"wrote {OUT / f'figS2_Td_swarm_all_distances{SUFFIX}.pdf'}")


# ----------------------------------------------------------------------
# Fig S3 — Significance counts (stoplight heatmap, instant read)
# Skipped for substrates without audit_25 cohort_summary.csv.
# ----------------------------------------------------------------------
_cohort_sum_csv = SRC / "cohort_summary.csv"
if not _cohort_sum_csv.exists():
    print(f"skipping figS3 (no {_cohort_sum_csv.name} for {SUBSTRATE})")
else:
    cohort_sum = pd.read_csv(_cohort_sum_csv)

    PAIR_S3 = [
        ("rest_pre", "task_test"),
        ("rest_pre", "rest_post"),
        ("task_test", "rest_post"),
    ]
    PAIR_S3_LABEL = {
        ("rest_pre", "task_test"):  r"RPre$\rightarrow$TT",
        ("rest_pre", "rest_post"):  r"RPre$\rightarrow$RPost",
        ("task_test", "rest_post"): r"TT$\rightarrow$RPost",
    }
    DIST_ORDER_S3 = ["S", "P", "F"]

    n_rows = len(PAIR_S3) * len(DIST_ORDER_S3)   # 9
    n_cols = len(BRAIN_BANDS_NAMES)              # 6

    matrix = np.full((n_rows, n_cols), np.nan)
    row_labels = []
    for pi, pair in enumerate(PAIR_S3):
        for di, d in enumerate(DIST_ORDER_S3):
            row_labels.append(rf"$d_{d}$  ·  {PAIR_S3_LABEL[pair]}")
            for bi, band in enumerate(BRAIN_BANDS_NAMES):
                sel = cohort_sum[
                    (cohort_sum.distance == d) & (cohort_sum.band == band)
                    & (cohort_sum.phase_A == pair[0])
                    & (cohort_sum.phase_B == pair[1])
                ]
                if not sel.empty:
                    matrix[pi * 3 + di, bi] = float(sel.iloc[0].n_plus)

    # Stoplight cmap with hard thresholds: <5 red, 5-7 amber, >=8 green
    cmap_st = LinearSegmentedColormap.from_list(
        "stoplight",
        [(0.0, "#c0392b"),
         (0.49, "#c0392b"),
         (0.50, "#f1c40f"),
         (0.79, "#f1c40f"),
         (0.80, "#27ae60"),
         (1.00, "#27ae60")],
        N=256,
    )
    norm_st = Normalize(vmin=0, vmax=10)

    fig, ax = plt.subplots(figsize=(11.5, 5.6))
    fig.subplots_adjust(top=0.92, bottom=0.18, left=0.22, right=0.98)
    im = ax.imshow(matrix, cmap=cmap_st, norm=norm_st, aspect="auto",
                   interpolation="nearest")

    # Annotate each cell with the count
    for i in range(n_rows):
        for j in range(n_cols):
            val = int(matrix[i, j])
            col_text = "white" if val < 5 or val >= 8 else "#222"
            ax.text(j, i, f"{val}/10",
                    ha="center", va="center",
                    fontsize=10, fontweight="bold", color=col_text)

    # Black separator lines between phase-pair groups
    for k in [3, 6]:
        ax.axhline(k - 0.5, color="k", lw=1.4)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                       fontsize=10)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=8.5)
    ax.tick_params(axis="x", labelsize=10, pad=4)
    ax.tick_params(axis="y", labelsize=8.5, pad=2)
    ax.set_title(r"Patients with $Z > 2$ vs the within-RPre split-half null"
                 r"  (out of 10)  —  cohort screening only, NOT the trace test",
                 fontsize=10.5, pad=8)

    # Stoplight legend at bottom
    handles_s3 = [
        Patch(facecolor="#27ae60",
              label=r"$n_+ \geq 8/10$  passes cohort screen"),
        Patch(facecolor="#f1c40f",
              label=r"$5 \leq n_+ \leq 7/10$  marginal"),
        Patch(facecolor="#c0392b",
              label=r"$n_+ \leq 4/10$  fails"),
    ]
    fig.legend(handles=handles_s3, loc="lower center",
               bbox_to_anchor=(0.5, 0.015), ncol=len(handles_s3),
               frameon=False, fontsize=9)

    fig.savefig(OUT / f"figS3_significance_counts{SUFFIX}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT / f'figS3_significance_counts{SUFFIX}.pdf'}")
