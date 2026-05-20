#!/usr/bin/env python3
"""β preprint Figure 3 — D(τ') phase mosaic across four phases.

For each patient in the n=10 β cohort, produce three PDF heatmap-mosaics
of the LRG communication-distance matrix D(τ') across the four phases
(rsPre, taskL, taskT, rsPost). One PDF per row/column ordering:

    - native  : channel-label order (no permutation)
    - fiedler : ascending Fiedler-vector value at rsPre (spectral order)
    - linkage : UPGMA leaf order on D(τ')^(rsPre) (block-structure-friendly)

The same permutation is applied to all four phase panels per file and to
the per-channel labels so that probe-family ticks align with the matrix
indices.

D convention: option (C) of
``.agents/preprint/established_results/00_open_methodology_question_lrg_D_convention.md``.
Raw propagator distance, max-symmetrization, no UPGMA cophenetic wrap.

    D(τ')_ij = 1 / ρ̂_ij(τ')          off-diagonal
             = 0                       diagonal
    ρ̂(τ)    = exp(-τL) / Tr exp(-τL)
    τ'      = 1 / λ_max (per-phase per-patient)

Layout: uses ``plot_fc_adjacency_row`` from ``lrg_eegfc.visuals.fc_templates``
to obtain (a) probe-family x/y tick labels via ``tick_labels="chnames"``
and (b) the canonical pre-padded last-panel gridspec so the shared
colorbar at the right does NOT steal width from the last imshow.

Outputs (PDF only, no PNG; no suptitle):

    data/preprint/figures/beta/D_tau_phase_mosaic/
        README.md
        summary.csv
        Pat_NN/
            Pat_NN_beta_D_tau_phase_mosaic_native.pdf
            Pat_NN_beta_D_tau_phase_mosaic_fiedler.pdf
            Pat_NN_beta_D_tau_phase_mosaic_linkage.pdf
            Pat_NN_beta_D_tau_phase_mosaic.md
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency_row

use_lrg_style()

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
BAND = "beta"
COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
PHASE_LABELS = {
    "rest_pre":   r"rsPre",
    "task_learn": r"taskL",
    "task_test":  r"taskT",
    "rest_post":  r"rsPost",
}
ORDERINGS = ("native", "fiedler", "linkage")

CACHE = ROOT / "data" / "cache" / "imcoh_lrg"
OUT_BASE = ROOT / "data" / "preprint" / "figures" / "beta" / "D_tau_phase_mosaic"
COLORBAR_LABEL = r"$D(\tau')_{ij} = 1/\hat{\rho}_{ij}(\tau')$"


# ----------------------------------------------------------------------------
# Cache + channel-label loading
# ----------------------------------------------------------------------------
def load_eigs(pat: str, phase: str) -> tuple[np.ndarray, np.ndarray]:
    f = CACHE / pat / f"{BAND}_{phase}_lrg_imcoh-abs.npz"
    npz = np.load(f)
    return npz["eigenvalues"], npz["eigenvectors"]


def load_channel_labels(pat: str) -> list[str]:
    """Load contact labels in FC-matrix order from channel_labels.csv.

    The file is already in the canonical channel ordering matching the
    LRG cache's N_nodes for every patient (Pat_10's task rows
    [53, 54, 55] are dropped at TIME-SERIES load time, but
    channel_labels.csv is already the canonical 113-channel list, so
    no additional drop is needed here).
    """
    path = SEEG_DATAPATH / pat / "channel_labels.csv"
    df = pd.read_csv(path)
    return df.iloc[:, 0].astype(str).tolist()


# ----------------------------------------------------------------------------
# D(τ') from cached LRG eigendecomposition
# ----------------------------------------------------------------------------
def compute_D_raw(eigvals: np.ndarray, eigvecs: np.ndarray) -> np.ndarray:
    """Raw LRG communication distance at τ = 1/λ_max.

    Matches ``audit_63.lrg_ultrametric_condensed`` MINUS the
    ``cophenet(linkage(...))`` wrap. Returns a symmetric (N, N) array
    with zero diagonal and strictly positive off-diagonal entries.
    """
    lam_max = float(eigvals.max())
    tau = 1.0 / lam_max
    diag_exp = np.exp(-tau * eigvals)
    rho = (eigvecs * diag_exp) @ eigvecs.T
    rho /= np.trace(rho)

    with np.errstate(divide="ignore", invalid="ignore"):
        Trho = 1.0 / rho

    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0.0)

    finite_pos = np.isfinite(Trho) & (Trho > 0.0)
    pos_mask = finite_pos | np.eye(Trho.shape[0], dtype=bool)
    if not pos_mask.all():
        cap = float(Trho[finite_pos].max()) if finite_pos.any() else 1e6
        Trho = np.where(pos_mask, Trho, cap)
        off = ~np.eye(Trho.shape[0], dtype=bool)
        Trho = np.where(off & (Trho <= 0.0), cap, Trho)
    return Trho


# ----------------------------------------------------------------------------
# Orderings
# ----------------------------------------------------------------------------
def native_order(N: int) -> np.ndarray:
    return np.arange(N)


def fiedler_order(eigvals: np.ndarray, eigvecs: np.ndarray) -> np.ndarray:
    order_eig = np.argsort(eigvals)
    fiedler_vec = eigvecs[:, order_eig[1]]
    return np.argsort(fiedler_vec)


def linkage_leaf_order(D: np.ndarray) -> np.ndarray:
    cond = squareform(D, checks=False)
    Z = linkage(cond, method="average")
    return np.asarray(dendrogram(Z, no_plot=True)["leaves"], dtype=int)


# ----------------------------------------------------------------------------
# Plotting (uses fc_templates.plot_fc_adjacency_row)
# ----------------------------------------------------------------------------
def plot_phase_mosaic(
    pat: str,
    Ds: dict,
    chan_labels: list[str],
    order: np.ndarray,
    ordering_name: str,
    out_pdf: Path,
) -> dict:
    """Plot 1×4 D(τ') row using fc_templates.plot_fc_adjacency_row.

    Probe-family x/y ticks via ``tick_labels="chnames"``. Shared
    LogNorm colorbar via the helper's pre-padded last-column trick
    (no steal from the last imshow).
    """
    # Reorder matrices + channel labels with the same permutation
    D_reord = {ph: Ds[ph][np.ix_(order, order)] for ph in PHASES}
    chan_reord = [chan_labels[i] for i in order]

    matrices = [D_reord[ph] for ph in PHASES]
    titles = [PHASE_LABELS[ph] for ph in PHASES]

    # Pre-clip extreme values (99th percentile across the row) so the
    # LogNorm vmax inside plot_fc_adjacency_row uses a robust scale
    # instead of being dragged by a few outlier pairs.
    all_pos = np.concatenate([
        D[~np.eye(D.shape[0], dtype=bool)][D[~np.eye(D.shape[0], dtype=bool)] > 0.0]
        for D in matrices
    ])
    vmax_clip = float(np.percentile(all_pos, 99.0))
    matrices_clipped = [np.minimum(D, vmax_clip) for D in matrices]

    fig, axes = plot_fc_adjacency_row(
        matrices_clipped,
        titles=titles,
        fc_method=None,
        band=BAND,
        tick_labels="chnames",
        channel_labels=chan_reord,
        log_scale=True,
        shared_scale=True,
        cmap="viridis",
        panel_size=2.8,
    )

    # Replace the auto colorbar label (empty when fc_method=None) with
    # the D(τ') label. The cax is the last axes appended to fig.
    panel_set = set(id(a) for a in axes)
    cax_candidates = [a for a in fig.axes if id(a) not in panel_set]
    if cax_candidates:
        # cax_candidates[-1] is the divider-attached colorbar axis.
        cax_candidates[-1].set_ylabel(COLORBAR_LABEL, fontsize=8)

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)

    # Per-phase value-range info for the companion .md
    per_phase_min_off = {
        ph: float(D_reord[ph][~np.eye(D_reord[ph].shape[0], dtype=bool)].min())
        for ph in PHASES
    }
    per_phase_max_off = {
        ph: float(D_reord[ph][~np.eye(D_reord[ph].shape[0], dtype=bool)].max())
        for ph in PHASES
    }

    return {
        "N": int(order.size),
        "ordering": ordering_name,
        "vmax_99th_pct_clip": vmax_clip,
        "per_phase_min_off": per_phase_min_off,
        "per_phase_max_off": per_phase_max_off,
    }


# ----------------------------------------------------------------------------
# Companion .md
# ----------------------------------------------------------------------------
def write_companion_md(
    pat: str, info_per_order: dict, lam_max_per_phase: dict, out_md: Path
) -> None:
    lines = [
        "---",
        f"name: preprint-fig-D-tau-mosaic-{pat.lower()}-{BAND}",
        "era: IMCOH_ABS_COHORT_N10",
        "status: current",
        "kind: figure-companion",
        f"band: {BAND}",
        f"patient: {pat}",
        "scope: D(tau') phase mosaic for the β preprint paragraph (Figure 3)",
        "---",
        "",
        f"# {pat} — β D(τ') phase mosaic",
        "",
        "## What is plotted",
        "",
        r"`D(τ')_ij = 1/ρ̂_ij(τ')` at `τ' = 1/λ_max`, with",
        r"`ρ̂(τ) = exp(-τL) / Tr exp(-τL)`.",
        "Max-symmetrization, diagonal = 0. Raw propagator distance, no UPGMA",
        "cophenetic wrap (option C of",
        "`.agents/preprint/established_results/00_open_methodology_question_lrg_D_convention.md`).",
        "",
        "Probe-family x/y tick labels are placed at each probe-group midpoint",
        "with minor ticks at probe boundaries (see `_apply_tick_labels` mode",
        "`chnames` in `lrg_eegfc.visuals.fc_templates`). When channels are",
        "reordered (`fiedler`, `linkage`), the same permutation is applied to",
        "the label sequence and probe groups are recomputed on the new order;",
        "if a probe's contacts end up non-consecutive, that probe label",
        "appears multiple times along the axis.",
        "",
        "## Per-phase τ' = 1/λ_max",
        "",
    ]
    for ph in PHASES:
        lines.append(
            f"- {ph}: λ_max = {lam_max_per_phase[ph]:.4f}, τ' = {1.0/lam_max_per_phase[ph]:.6f}"
        )

    lines.extend([
        "",
        "## Three orderings (same permutation applied to matrix axes + channel labels)",
        "",
        "- `native`  — channel-label order from `channel_labels.csv` (probes consecutive)",
        "- `fiedler` — ascending Fiedler-vector value at rsPre (L^(rsPre) second-smallest eigenvector)",
        "- `linkage` — UPGMA leaf order on D(τ')^(rsPre) condensed",
        "",
        "## Files in this folder",
        "",
        "| Ordering | File | N | vmax (99th %ile, shared across phases) |",
        "|---|---|---|---|",
    ])
    for order_name in ORDERINGS:
        info = info_per_order[order_name]
        fname = f"{pat}_{BAND}_D_tau_phase_mosaic_{order_name}.pdf"
        lines.append(
            f"| {order_name} | `{fname}` | {info['N']} "
            f"| {info['vmax_99th_pct_clip']:.3e} |"
        )

    lines.extend([
        "",
        "## Per-phase off-diagonal range (pre-clipping)",
        "",
        "| Phase | min(D off-diag) | max(D off-diag) |",
        "|---|---|---|",
    ])
    info_native = info_per_order["native"]
    for ph in PHASES:
        lines.append(
            f"| {ph} | {info_native['per_phase_min_off'][ph]:.3e} "
            f"| {info_native['per_phase_max_off'][ph]:.3e} |"
        )

    lines.extend([
        "",
        "## Provenance",
        "",
        "- Build script: `scripts/02_preprint/preprint_01_beta_D_tau_phase_mosaic.py`",
        f"- LRG cache: `data/cache/imcoh_lrg/{pat}/{BAND}_{{phase}}_lrg_imcoh-abs.npz`",
        "  (reads `eigenvalues`, `eigenvectors` only)",
        f"- Channel labels: `data/raw/stereoeeg_patients/{pat}/channel_labels.csv`",
        "- Layout helper: `lrg_eegfc.visuals.fc_templates.plot_fc_adjacency_row`",
        "  (probe-family ticks via `tick_labels=chnames`; shared colorbar via",
        "  pre-padded last-column gridspec, so the last imshow keeps its width)",
        "- D recipe: top-level docstring of the build script",
        "- D convention: `.agents/preprint/established_results/00_open_methodology_question_lrg_D_convention.md` (option C)",
        "",
    ])

    out_md.write_text("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> None:
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict] = []

    for pat in COHORT:
        print(f"[{pat}] loading 4 phases...", flush=True)
        try:
            eigs = {ph: load_eigs(pat, ph) for ph in PHASES}
            chan_labels = load_channel_labels(pat)
        except FileNotFoundError as e:
            print(f"  SKIP: missing input: {e}")
            continue

        Ds = {ph: compute_D_raw(*eigs[ph]) for ph in PHASES}
        N = Ds[PHASES[0]].shape[0]
        if any(D.shape != (N, N) for D in Ds.values()):
            print(f"  SKIP {pat}: phase shape mismatch")
            continue
        if len(chan_labels) != N:
            print(f"  SKIP {pat}: channel-label count {len(chan_labels)} != N {N}")
            continue

        lam_max_per_phase = {ph: float(eigs[ph][0].max()) for ph in PHASES}

        eigvals_pre, eigvecs_pre = eigs["rest_pre"]
        orders = {
            "native":  native_order(N),
            "fiedler": fiedler_order(eigvals_pre, eigvecs_pre),
            "linkage": linkage_leaf_order(Ds["rest_pre"]),
        }

        out_dir = OUT_BASE / pat
        info_per_order = {}
        for order_name in ORDERINGS:
            out_pdf = out_dir / f"{pat}_{BAND}_D_tau_phase_mosaic_{order_name}.pdf"
            info = plot_phase_mosaic(
                pat, Ds, chan_labels, orders[order_name], order_name, out_pdf,
            )
            info_per_order[order_name] = info
            summary_rows.append({
                "patient": pat, "ordering": order_name, "N": info["N"],
                "vmax_99th_pct_clip": info["vmax_99th_pct_clip"],
            })
            print(f"  wrote {out_pdf.relative_to(ROOT)}")

        out_md = out_dir / f"{pat}_{BAND}_D_tau_phase_mosaic.md"
        write_companion_md(pat, info_per_order, lam_max_per_phase, out_md)

    # Top-level README + summary CSV
    readme = OUT_BASE / "README.md"
    readme.write_text(
        "# β D(τ') phase mosaic — preprint figure set\n"
        "\n"
        f"One folder per patient ({len(COHORT)} total). Each folder contains:\n"
        "\n"
        "- Three PDFs (one per ordering: `native`, `fiedler`, `linkage`)\n"
        "- One companion `.md` with shapes, value ranges, λ_max per phase, provenance\n"
        "\n"
        "## D convention\n"
        "\n"
        "Raw LRG propagator distance `D(τ')_ij = 1/ρ̂_ij(τ')` at `τ' = 1/λ_max`,\n"
        "with `ρ̂(τ) = exp(-τL)/Tr exp(-τL)`, max-symmetrized, diagonal = 0,\n"
        "no UPGMA cophenetic wrap.\n"
        "\n"
        "Reference: `.agents/preprint/established_results/`\n"
        "`00_open_methodology_question_lrg_D_convention.md` option (C).\n"
        "\n"
        "## Layout\n"
        "\n"
        "Uses `lrg_eegfc.visuals.fc_templates.plot_fc_adjacency_row` for:\n"
        "\n"
        "- Probe-family x/y ticks at probe-group midpoints (`tick_labels=chnames`)\n"
        "- Shared LogNorm colorbar at the row's right edge, with the last\n"
        "  panel's gridspec column pre-padded so the colorbar does not steal\n"
        "  width from the last imshow.\n"
        "\n"
        "## Build script\n"
        "\n"
        "`scripts/02_preprint/preprint_01_beta_D_tau_phase_mosaic.py`\n"
    )

    csv_path = OUT_BASE / "summary.csv"
    with csv_path.open("w") as f:
        writer = csv.DictWriter(
            f, fieldnames=["patient", "ordering", "N", "vmax_99th_pct_clip"]
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\ndone. {len(summary_rows)} PDFs generated.")
    print(f"summary CSV: {csv_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
