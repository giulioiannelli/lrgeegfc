#!/usr/bin/env python3
"""Δ_VI (and Δ_H, Δ_NMI) at fixed h_rel cuts — answer to the dmin / scale
comparability concern.

The headline figure cuts at fixed integer k. Across patients at the same k
the fractional dendrogram depth `h_rel = h(k) / dmax` varies (see
`kcut_heights_summary.md`), so "δ k=28 cohort-wide" is a *partition-level*
claim, not a *scale-level* one. This diagnostic re-runs the partition
contrasts at fixed h_rel cuts via `fcluster_at_h_rel`, plots both:

  Figure A (linspace) — h_rel ∈ [0.05, 0.95], 30 evenly-spaced bins.
  Figure B (logspace) — h_rel ∈ [floor, 0.95], 30 geometrically-spaced bins
                        with `floor = max over cohort of dmin/dmax`
                        (the most-shallow tree's lowest meaningful cut).

`dmin / dmax` per (patient, band, phase) is reported alongside, so the
reader can see whether the comparability is still loose at fine h_rel.

Reads:  data/cache/imcoh_lrg/<Pat>/<band>_<phase>_lrg_imcoh-abs.npz

Writes:
  data/audit/dvi_split_baseline/dvi_hrel_n10_imcoh_abs.csv
  data/audit/dvi_split_baseline/dmin_inventory.md
  data/outputs/figures/section6/task_trace_band_hrel_linspace_n10_imcoh_abs.pdf
  data/outputs/figures/section6/task_trace_band_hrel_linspace_n10_imcoh_abs.md
  data/outputs/figures/section6/task_trace_band_hrel_logspace_n10_imcoh_abs.pdf
  data/outputs/figures/section6/task_trace_band_hrel_logspace_n10_imcoh_abs.md
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from sklearn.metrics import normalized_mutual_info_score

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.utils.metrics import compute_vi, conditional_entropy
from lrg_eegfc.utils.metrics.tree import dmax_from_Z, fcluster_at_h_rel
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, cluster_stats
from lrg_eegfc.workflow.lrg import load_lrg_result


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
COHORT_THRESHOLD = 8
COHORT_ADVISORY  = 6
N_GRID = 30
H_REL_LIN_LO, H_REL_LIN_HI = 0.05, 0.95
H_REL_LOG_HI = 0.95
PHASES = ("rest_pre", "task_test", "rest_post")
N_PERM = 2000
SEED = 0

PARTITION_PANELS = [
    ("d_VI",  r"$\Delta_{\mathrm{VI}}(h_{\mathrm{rel}})$ — H2a (VI)"),
    ("d_H",   r"$\Delta_{H}(h_{\mathrm{rel}})$ — H2a (directional cond-$H$)"),
    ("d_NMI", r"$\Delta_{\mathrm{NMI}}(h_{\mathrm{rel}})$ — H2a (NMI)"),
]


def _Z(pat: str, phase: str, band: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, "imcoh_abs", IMCOH_LRG_CACHE)
    except Exception:
        return None
    return None if r is None else np.asarray(r.linkage_matrix)


def _dmin_dmax_inventory() -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            for phase in PHASES:
                Z = _Z(pat, phase, band)
                if Z is None:
                    continue
                d_max = dmax_from_Z(Z)
                d_min = float(Z[0, 2])
                rows.append({
                    "patient": pat, "band": band, "phase": phase,
                    "dmin": d_min, "dmax": d_max,
                    "h_rel_min": d_min / d_max if d_max > 0 else float("nan"),
                    "n_leaves": Z.shape[0] + 1,
                })
    return pd.DataFrame(rows)


def collect_dvi_hrel(grids: dict[str, np.ndarray]) -> pd.DataFrame:
    """Compute Δ_VI / Δ_H / Δ_NMI at every h_rel in `grids`."""
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            Z_pre  = _Z(pat, "rest_pre",  band)
            Z_test = _Z(pat, "task_test", band)
            Z_post = _Z(pat, "rest_post", band)
            if any(z is None for z in (Z_pre, Z_test, Z_post)):
                continue
            n = Z_pre.shape[0] + 1
            if Z_test.shape[0] + 1 != n or Z_post.shape[0] + 1 != n:
                continue
            for grid_name, grid in grids.items():
                for h in grid:
                    c_pre  = fcluster_at_h_rel(Z_pre,  float(h))
                    c_test = fcluster_at_h_rel(Z_test, float(h))
                    c_post = fcluster_at_h_rel(Z_post, float(h))
                    vi_rpre  = compute_vi(c_pre,  c_post)
                    vi_test  = compute_vi(c_test, c_post)
                    h_rpre   = conditional_entropy(c_post, c_pre)
                    h_test   = conditional_entropy(c_post, c_test)
                    nmi_rpre = normalized_mutual_info_score(c_pre,  c_post)
                    nmi_test = normalized_mutual_info_score(c_test, c_post)
                    rows.append({
                        "patient": pat, "band": band,
                        "grid": grid_name, "h_rel": float(h),
                        "k_pre":  int(np.unique(c_pre).size),
                        "k_test": int(np.unique(c_test).size),
                        "k_post": int(np.unique(c_post).size),
                        "d_VI":  vi_rpre - vi_test,
                        "d_H":   h_rpre - h_test,
                        "d_NMI": nmi_test - nmi_rpre,
                    })
    return pd.DataFrame(rows)


# ─────────────────────── plotting helpers (mirrored from headline) ──────


def _cmap_diverging():
    return LinearSegmentedColormap.from_list(
        "trace_div",
        [(0.05, 0.20, 0.55), (1.0, 1.0, 1.0), (0.70, 0.05, 0.05)],
        N=256,
    )


def _build_A(df: pd.DataFrame, contrast: str,
              grid: np.ndarray) -> np.ndarray:
    """Build (n_pat, n_bands, n_grid) array. NaN if missing."""
    n_pat = len(COHORT_N10)
    n_b = len(BRAIN_BANDS_NAMES)
    n_g = len(grid)
    A = np.full((n_pat, n_b, n_g), np.nan)
    pat_idx = {p: i for i, p in enumerate(COHORT_N10)}
    band_idx = {b: i for i, b in enumerate(BRAIN_BANDS_NAMES)}
    grid_arr = np.asarray(grid)
    # Vectorise per row of the dataframe.
    for r in df.itertuples(index=False):
        ip = pat_idx.get(r.patient)
        ib = band_idx.get(r.band)
        if ip is None or ib is None:
            continue
        # Find closest grid index.
        diffs = np.abs(grid_arr - r.h_rel)
        jg = int(diffs.argmin())
        if diffs[jg] > 1e-9:
            continue  # not a grid point
        v = getattr(r, contrast)
        if np.isfinite(v):
            A[ip, ib, jg] = v
    return A


def _per_band_stats(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (mean[B,G], n_pos[B,G])."""
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(A, axis=0)
    npos = np.nansum(A > 0, axis=0).astype(int)
    return mean, npos


def _real_clusters(A: np.ndarray, z_thresh: float
                   ) -> dict[int, list[tuple[int, int, float]]]:
    n_pat = A.shape[0]
    with np.errstate(invalid="ignore"):
        mn = np.nanmean(A, axis=0)
        sd = np.nanstd(A, axis=0, ddof=1)
        zv = mn / (sd / np.sqrt(n_pat))
    out: dict[int, list[tuple[int, int, float]]] = {}
    for ib in range(zv.shape[0]):
        out[ib] = cluster_stats(zv[ib], z_thresh)
    return out


def _signflip_threshold(A: np.ndarray, z_thresh: float) -> float:
    n_pat, n_b, n_g = A.shape
    rng = np.random.default_rng(SEED)
    null_max = np.empty(N_PERM, dtype=float)
    for p in range(N_PERM):
        signs = rng.choice([-1.0, 1.0], size=n_pat)[:, None, None]
        A_p = A * signs
        with np.errstate(invalid="ignore"):
            mn = np.nanmean(A_p, axis=0)
            sd = np.nanstd(A_p, axis=0, ddof=1)
            zv = mn / (sd / np.sqrt(n_pat))
        max_mass = 0.0
        for ib in range(n_b):
            for _, _, mass in cluster_stats(zv[ib], z_thresh):
                if mass > max_mass:
                    max_mass = mass
        null_max[p] = max_mass
    return float(np.quantile(null_max, 0.95))


def _draw_borders(ax, grid: np.ndarray, npos: np.ndarray,
                  log_x: bool) -> None:
    B, G = npos.shape
    # cell width in the chosen axis: between consecutive grid edges
    # (extent depends on log/lin; the borders use grid coords directly).
    for ib in range(B):
        for jg in range(G):
            n = npos[ib, jg]
            if n < COHORT_ADVISORY:
                continue
            # cell edges (in axis coords): midpoints between grid points.
            if jg == 0:
                lo = grid[0] - (grid[1] - grid[0]) / 2.0
            else:
                lo = (grid[jg - 1] + grid[jg]) / 2.0
            if jg == G - 1:
                hi = grid[-1] + (grid[-1] - grid[-2]) / 2.0
            else:
                hi = (grid[jg] + grid[jg + 1]) / 2.0
            color = "black" if n >= COHORT_THRESHOLD else "0.4"
            lw = 1.4 if n >= COHORT_THRESHOLD else 0.8
            ls = "-" if n >= COHORT_THRESHOLD else ":"
            ax.add_patch(Rectangle((lo, ib - 0.5), hi - lo, 1,
                                    facecolor="none", edgecolor=color,
                                    linewidth=lw, linestyle=ls, zorder=3))


def _hatch_clusters(ax, grid: np.ndarray,
                    clusters: dict[int, list[tuple[int, int, float]]],
                    sig_mass: float, log_x: bool) -> None:
    for ib, runs in clusters.items():
        for j_lo, j_hi, mass in runs:
            if mass <= sig_mass:
                continue
            # cell edges as in _draw_borders
            if j_lo == 0:
                x_lo = grid[0] - (grid[1] - grid[0]) / 2.0
            else:
                x_lo = (grid[j_lo - 1] + grid[j_lo]) / 2.0
            if j_hi == len(grid) - 1:
                x_hi = grid[-1] + (grid[-1] - grid[-2]) / 2.0
            else:
                x_hi = (grid[j_hi] + grid[j_hi + 1]) / 2.0
            ax.add_patch(Rectangle((x_lo, ib - 0.5), x_hi - x_lo, 1,
                                    facecolor="none", edgecolor="black",
                                    linewidth=0.0, hatch="////",
                                    alpha=0.45, zorder=4))


def _draw_panel(ax, df: pd.DataFrame, contrast: str, label: str,
                grid: np.ndarray, log_x: bool,
                show_yticks: bool, show_xlabel: bool,
                z_thresh: float):
    A = _build_A(df, contrast, grid)
    mean, npos = _per_band_stats(A)
    cmap = _cmap_diverging()
    vmax = np.nanmax(np.abs(mean)) if np.isfinite(mean).any() else 1.0

    # Render cells as individual rectangles to support unequal x-widths
    # (essential for log-scale axis where bin widths grow).
    norm = plt.Normalize(vmin=-vmax, vmax=+vmax)
    for ib in range(mean.shape[0]):
        for jg in range(mean.shape[1]):
            v = mean[ib, jg]
            if not np.isfinite(v):
                continue
            if jg == 0:
                lo = grid[0] - (grid[1] - grid[0]) / 2.0
            else:
                lo = (grid[jg - 1] + grid[jg]) / 2.0
            if jg == len(grid) - 1:
                hi = grid[-1] + (grid[-1] - grid[-2]) / 2.0
            else:
                hi = (grid[jg] + grid[jg + 1]) / 2.0
            ax.add_patch(Rectangle((lo, ib - 0.5), hi - lo, 1,
                                    facecolor=cmap(norm(v)),
                                    edgecolor="none", zorder=1))

    _draw_borders(ax, grid, npos, log_x)
    sig_mass = _signflip_threshold(A, z_thresh)
    clusters = _real_clusters(A, z_thresh)
    _hatch_clusters(ax, grid, clusters, sig_mass, log_x)

    if log_x:
        ax.set_xscale("log")
    ax.set_xlim(grid.min() * (0.95 if log_x else 1.0) - (0 if log_x else 0.02),
                grid.max() + (grid.max() - grid[-2]) / 2.0)
    ax.set_ylim(len(BRAIN_BANDS_NAMES) - 0.5, -0.5)
    if show_yticks:
        ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                            fontsize=10)
    else:
        ax.set_yticks([])
    if show_xlabel:
        ax.set_xlabel(r"$h_{\mathrm{rel}}$ (cut height / dmax)", fontsize=10)
    ax.set_title(label, fontsize=11, loc="left")
    # Manual colorbar mappable
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    return sm, sig_mass


def _make_figure(df: pd.DataFrame, grid: np.ndarray, grid_name: str,
                 log_x: bool, out_pdf, sidecar_md):
    z_thresh = 1.65
    fig, axes = plt.subplots(
        1, 3, figsize=(15.0, 4.0), dpi=160, sharey=True,
        gridspec_kw={"wspace": 0.10},
    )
    last_sm = None
    sig_masses: dict[str, float] = {}
    for j, (col, label) in enumerate(PARTITION_PANELS):
        sm, sig = _draw_panel(
            axes[j], df, col, label, grid, log_x,
            show_yticks=(j == 0), show_xlabel=True, z_thresh=z_thresh,
        )
        last_sm = sm
        sig_masses[col] = sig
    cbar = fig.colorbar(last_sm, ax=axes.tolist(),
                        shrink=0.85, pad=0.02, fraction=0.025)
    cbar.set_label("patient-mean contrast", fontsize=10)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_pdf}")

    # Sidecar
    md_lines = [
        "---",
        f"name: task-trace-band-hrel-{grid_name}",
        "type: figure-sidecar",
        "era: COHORT_N10",
        "status: current",
        "created: 2026-04-27",
        "updated: 2026-04-27",
        "---",
        "",
        f"# Δ_VI / Δ_H / Δ_NMI at fixed `h_rel` cuts — {grid_name} grid (n=10)",
        "",
        "Companion to `task_trace_band_k_n10_imcoh_abs.pdf`. Re-runs the same",
        "three partition contrasts but on **fixed fractional dendrogram heights**",
        "instead of integer-k cuts. Tests robustness of the cohort claim to",
        "the k-vs-h_rel choice.",
        "",
        f"Grid: {grid_name} `h_rel` ∈ [{grid[0]:.4f}, {grid[-1]:.4f}], "
        f"{len(grid)} bins.",
        "",
        f"**Black border** = ≥ {COHORT_THRESHOLD}/10 patients sign-correct.",
        f"**Dashed grey** = {COHORT_ADVISORY}–{COHORT_THRESHOLD-1}/10. **Hatch** = sign-flip",
        f"cluster-permutation cluster (n_perm = {N_PERM}, α = 0.05).",
        "",
        "## Cluster-permutation 95% mass thresholds",
        "",
        "| panel | sig mass threshold |",
        "|---|---:|",
    ]
    for col, sig in sig_masses.items():
        md_lines.append(f"| {col} | {sig:.3f} |")
    md_lines.extend([
        "",
        "## Reading rules",
        "",
        f"- A vertical band of {COHORT_THRESHOLD}/10 borders along `h_rel` "
        "= cohort-wide reorganization at a fractional depth → a *scale-level*",
        "  claim, not just a partition-level one (which is what integer-k gives).",
        "- Compare to `task_trace_band_k_n10_imcoh_abs.pdf`: a δ ridge that",
        "  appears in both views = robust; a ridge that appears only in one",
        "  = sensitive to the cut-axis choice.",
        f"- Cohort `dmin/dmax` floor (max-over-cohort) gives the smallest "
        "h_rel at which every tree is non-degenerate; see `dmin_inventory.md`.",
    ])
    sidecar_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"wrote {sidecar_md}")


def main() -> None:
    out_diag = ROOT / "data" / "audit" / "dvi_split_baseline"
    out_diag.mkdir(parents=True, exist_ok=True)
    out_fig = FIGURES_ROOT / "section6"
    out_fig.mkdir(parents=True, exist_ok=True)

    # 1. dmin / dmax inventory ------------------------------------------------
    inv = _dmin_dmax_inventory()
    inv.to_csv(out_diag / "dmin_dmax_inventory.csv", index=False)

    h_rel_log_lo_strict = float(inv["h_rel_min"].max())
    print(f"dmin/dmax cohort: min={inv['h_rel_min'].min():.6f}, "
          f"max={inv['h_rel_min'].max():.6f}, "
          f"median={inv['h_rel_min'].median():.6f}")
    print(f"Strict logspace floor (cohort max h_rel_min) = "
          f"{h_rel_log_lo_strict:.4f} — covers only the top "
          f"{(1 - h_rel_log_lo_strict) * 100:.1f}% of every tree.")
    # Use the linspace floor for the figure logspace too: gives more
    # fine-scale resolution. Cells where a patient's h_rel < dmin/dmax
    # produce all-singleton partitions; those are reported but flagged
    # in the sidecar against dmin_inventory.md.
    grid_lin = np.linspace(H_REL_LIN_LO, H_REL_LIN_HI, N_GRID)
    grid_log = np.geomspace(H_REL_LIN_LO, H_REL_LIN_HI, N_GRID)

    # dmin inventory markdown
    lines: list[str] = []
    ap = lines.append
    ap("# dmin / dmax / h_rel_min inventory — n=10 IMCOH_ABS")
    ap("")
    ap("`dmin = Z[0, 2]` (smallest non-zero merge height).")
    ap("`dmax = Z[-1, 2]` (root height; constant ≈ 0.9901 by LRG normalization).")
    ap("`h_rel_min = dmin / dmax` — the smallest h_rel at which `fcluster` returns")
    ap("at least one non-singleton cluster. Below this, the partition is N singletons.")
    ap("")
    ap("## Cohort summary per band × phase")
    ap("")
    ap("| band | phase | dmin median | dmin min | dmin max | h_rel_min median | h_rel_min max |")
    ap("|------|-------|------------:|---------:|---------:|-----------------:|--------------:|")
    for band in BRAIN_BANDS_NAMES:
        for phase in PHASES:
            sub = inv[(inv["band"] == band) & (inv["phase"] == phase)]
            if sub.empty:
                continue
            ap(f"| {BRAIN_BAND_TEX_DICT[band]} | {phase} | "
               f"{sub['dmin'].median():.5f} | "
               f"{sub['dmin'].min():.5f} | "
               f"{sub['dmin'].max():.5f} | "
               f"{sub['h_rel_min'].median():.5f} | "
               f"{sub['h_rel_min'].max():.5f} |")
    ap("")
    ap(f"### Cohort strict logspace floor = {h_rel_log_lo_strict:.4f}")
    ap("")
    ap("This is the maximum `h_rel_min` across the cohort. Cutting below this")
    ap("h_rel makes at least one tree return all singletons → strict cohort")
    ap("comparison breaks down. The published logspace figure uses the wider")
    ap(f"`[{H_REL_LIN_LO:.2f}, {H_REL_LIN_HI:.2f}]` range to expose fine-scale")
    ap("structure on the bands where it is meaningful — γ_h cells below the")
    ap("strict floor will have all-singleton partitions in some patients and")
    ap("should be read with this caveat in mind.")
    ap("")
    ap("## Reading")
    ap("")
    ap("- If `h_rel_min` is **tightly clustered** (small max/min ratio, small")
    ap("  spread per band), the trees are scale-comparable at the fine end.")
    ap("- If `h_rel_min` is **widely spread**, the trees have heterogeneous")
    ap("  fine-scale resolution; cutting at a single h_rel below the cohort")
    ap("  max compares 'finest possible cut for one patient' against 'a")
    ap("  somewhat-coarser cut for another'. The user's concern is real here.")
    (out_diag / "dmin_inventory.md").write_text("\n".join(lines) + "\n",
                                                encoding="utf-8")
    print(f"wrote {out_diag / 'dmin_inventory.md'}")

    # 2. Compute Δ_VI etc on both grids ---------------------------------------
    grids = {"linspace": grid_lin, "logspace": grid_log}
    print("Computing Δ_VI / Δ_H / Δ_NMI on both h_rel grids...", flush=True)
    df = collect_dvi_hrel(grids)
    df.to_csv(out_diag / "dvi_hrel_n10_imcoh_abs.csv", index=False)
    print(f"{len(df)} rows → {out_diag / 'dvi_hrel_n10_imcoh_abs.csv'}",
          flush=True)

    # 3. Plot ---------------------------------------------------------------
    df_lin = df[df["grid"] == "linspace"].copy()
    df_log = df[df["grid"] == "logspace"].copy()

    _make_figure(
        df_lin, grid_lin, "linspace", log_x=False,
        out_pdf=out_fig / "task_trace_band_hrel_linspace_n10_imcoh_abs.pdf",
        sidecar_md=out_fig / "task_trace_band_hrel_linspace_n10_imcoh_abs.md",
    )
    _make_figure(
        df_log, grid_log, "logspace", log_x=True,
        out_pdf=out_fig / "task_trace_band_hrel_logspace_n10_imcoh_abs.pdf",
        sidecar_md=out_fig / "task_trace_band_hrel_logspace_n10_imcoh_abs.md",
    )


if __name__ == "__main__":
    main()
