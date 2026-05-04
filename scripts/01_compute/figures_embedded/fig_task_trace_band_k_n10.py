#!/usr/bin/env python3
"""Headline 5-panel band×k task-trace map at n=10 (IMCOH_ABS).

Surfaces existing per-(band, k) cohort evidence — does not re-test, does
not invent measures. The 5 panels are 5 operationalizations of the same
residual claim "rest_post − rest_pre alignment with task − rest_pre":

  Top row    — partition-level Δ_VI, Δ_H (directional cond-H), Δ_NMI
  Bottom row — block-pair coactivation Δρ_H2d(k) / continuous ρ_H2c (no k)

Note: Δ_ARI is intentionally excluded per `feedback_no_ari_in_figures.md`
(2026-04-26). The CSV retains the `d_ARI` column for future reference;
this figure does not display it.

For each (band, k) cell:
  - colour = patient-mean contrast
  - black solid border  → ≥ 8/10 patients sign-correct (cohort-wide)
  - dashed grey border  → 6–7/10 (advisory)
  - diagonal hatch      → cell sits inside a sign-flip cluster-permutation
                           cluster significant at α = 0.05 (5 000 permutations)

H2c collapses k into a per-band scalar; rendered as a 2-column band-strip
{task, learn} on the bottom-right of the grid.

Inputs (read-only, post-Pat_14 backfill):
  data/reports/imcoh_vi/h2_partition_multiscale_raw.csv   (n=10, 2880 rows)
  data/reports/imcoh_vi/h2d_persistence_raw.csv           (n=10)
  data/reports/imcoh_vi/h2c_ultrametric_drift_raw.csv     (n=10)

Output:
  data/outputs/figures/section6/task_trace_band_k_n10_imcoh_abs.pdf
  data/outputs/figures/section6/task_trace_band_k_n10_imcoh_abs.md  (sidecar)

Era: IMCOH_ABS × COHORT_N10 (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15).
Cohort threshold ≥ 8/10 (preserves the prior 7/9 ≈ 78% fraction).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import REPORTS_ROOT, FIGURES_ROOT
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, cluster_stats


COHORT_N10 = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
              "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
COHORT_THRESHOLD = 8  # 8/10 (matches CLAUDE.md current era)
COHORT_ADVISORY  = 6  # 6/10 dashed grey border
K_RANGE = list(range(2, 50))
N_PERM  = 5000
SEED    = 0

PARTITION_PANELS = [
    ("d_VI",  r"$\Delta_{\mathrm{VI}}(k)$ — H2a (VI)"),
    ("d_H",   r"$\Delta_{H}(k)$ — H2a (directional cond-$H$)"),
    ("d_NMI", r"$\Delta_{\mathrm{NMI}}(k)$ — H2a (NMI)"),
]


def _cmap_diverging():
    return LinearSegmentedColormap.from_list(
        "trace_div",
        [(0.05, 0.20, 0.55), (1.0, 1.0, 1.0), (0.70, 0.05, 0.05)],
        N=256,
    )


def _per_band_stats_k(df: pd.DataFrame, contrast: str
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (mean, frac_pos, n_pos) per (band, k).

    mean[b,k]      — patient-mean contrast
    frac_pos[b,k]  — fraction of patients with contrast > 0
    n_pos[b,k]     — count of patients with contrast > 0  (out of n=10)
    """
    B, K = len(BRAIN_BANDS_NAMES), len(K_RANGE)
    mean = np.full((B, K), np.nan)
    frac = np.full((B, K), np.nan)
    npos = np.full((B, K), 0, dtype=int)
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        for jk, k in enumerate(K_RANGE):
            x = df[(df["band"] == band) & (df["k"] == k)][contrast].to_numpy()
            x = x[np.isfinite(x)]
            if x.size == 0:
                continue
            mean[ib, jk] = x.mean()
            frac[ib, jk] = (x > 0).sum() / x.size
            npos[ib, jk] = int((x > 0).sum())
    return mean, frac, npos


def _z_per_k(df: pd.DataFrame, contrast: str) -> np.ndarray:
    """Per (band, k) Wilcoxon z-statistic across patients (one-sided '> 0')."""
    B, K = len(BRAIN_BANDS_NAMES), len(K_RANGE)
    Z = np.full((B, K), np.nan)
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        for jk, k in enumerate(K_RANGE):
            x = df[(df["band"] == band) & (df["k"] == k)][contrast].to_numpy()
            x = x[np.isfinite(x)]
            if x.size < 3:
                continue
            z, _ = wilcoxon_z(x)
            Z[ib, jk] = z
    return Z


def _signflip_cluster_threshold(df: pd.DataFrame, contrast: str,
                                z_thresh: float, n_perm: int = N_PERM,
                                seed: int = SEED) -> float:
    """Sign-flip null: max-cluster mass distribution across (band, k) jointly.

    Maris-Oostenveld 2007. Per permutation we flip the sign of each patient's
    full (band, k) contrast vector; recompute per-(band, k) Wilcoxon z;
    extract supra-threshold clusters along k for each band; record the largest
    cluster-mass observed across all (band, perm). Threshold = 95th percentile
    of the resulting null distribution.
    """
    rng = np.random.default_rng(seed)
    bands = BRAIN_BANDS_NAMES
    pivot = (df.pivot_table(index=["patient"], columns=["band", "k"],
                            values=contrast, aggfunc="first")
             .reindex(index=COHORT_N10))
    A = pivot.values  # (n_patients, n_bands * n_k) — wide matrix, may have NaN

    null_max = np.empty(n_perm, dtype=float)
    n_pat = A.shape[0]
    for p in range(n_perm):
        signs = rng.choice([-1.0, 1.0], size=n_pat)[:, None]
        A_p = A * signs
        # Recompute per (band, k) one-sample Wilcoxon-like z under sign-flip.
        # Use the simple t-style z = mean / (std/sqrt(n)) as the cluster
        # statistic — under sign-flip null this is exchangeable and
        # numerically stable for n=10. The headline z (true data) is also
        # recomputed in the same form to keep the null comparable.
        with np.errstate(invalid="ignore"):
            mean = np.nanmean(A_p, axis=0)
            std  = np.nanstd(A_p, axis=0, ddof=1)
            z_p  = mean / (std / np.sqrt(n_pat))
        z_grid = z_p.reshape(len(bands), -1)  # (B, K)
        # Largest cluster mass across bands.
        max_mass = 0.0
        for ib in range(len(bands)):
            cls = cluster_stats(z_grid[ib], z_thresh)
            for _, _, mass in cls:
                if mass > max_mass:
                    max_mass = mass
        null_max[p] = max_mass
    return float(np.quantile(null_max, 0.95))


def _real_clusters_t_z(df: pd.DataFrame, contrast: str, z_thresh: float
                       ) -> tuple[np.ndarray, dict[int, list[tuple[int, int, float]]]]:
    """Compute real-data t-style z and cluster runs for hatching.

    Returns (Z_real (B,K), {band_idx: [(k_start_idx, k_end_idx, mass), ...]}).
    """
    bands = BRAIN_BANDS_NAMES
    pivot = (df.pivot_table(index=["patient"], columns=["band", "k"],
                            values=contrast, aggfunc="first")
             .reindex(index=COHORT_N10))
    A = pivot.values  # (n_pat, B*K)
    n_pat = A.shape[0]
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(A, axis=0)
        std  = np.nanstd(A, axis=0, ddof=1)
        z    = mean / (std / np.sqrt(n_pat))
    Z = z.reshape(len(bands), -1)
    out: dict[int, list[tuple[int, int, float]]] = {}
    for ib in range(len(bands)):
        out[ib] = cluster_stats(Z[ib], z_thresh)
    return Z, out


def _hatch_cells(ax, clusters: dict[int, list[tuple[int, int, float]]],
                 sig_mass_thresh: float, k_offset: float = 0.0) -> None:
    """Overlay diagonal-hatch on cells inside significant clusters (mass > thresh).

    Cells are rendered as transparent rectangles with `////` hatch.
    """
    for ib, runs in clusters.items():
        for k_lo, k_hi, mass in runs:
            if mass <= sig_mass_thresh:
                continue
            # Each cell at (band ib, k = K_RANGE[jk]) is at imshow coord
            # (k = K_RANGE[jk]+k_offset, y=ib). With aspect="auto",
            # extent is set in data coords below. We use the same extent.
            x_lo = K_RANGE[k_lo] - 0.5
            x_hi = K_RANGE[k_hi] + 0.5
            y_lo = ib - 0.5
            y_hi = ib + 0.5
            rect = Rectangle((x_lo, y_lo), x_hi - x_lo, y_hi - y_lo,
                             facecolor="none", edgecolor="black",
                             linewidth=0.0, hatch="////", alpha=0.45,
                             zorder=4)
            ax.add_patch(rect)


def _draw_borders(ax, npos: np.ndarray) -> None:
    """Per-cell border encoding cohort unanimity.

    Black solid border  → ≥ COHORT_THRESHOLD (8/10)
    Dashed grey border  → COHORT_ADVISORY (6) ≤ npos < COHORT_THRESHOLD
    """
    B, K = npos.shape
    for ib in range(B):
        for jk in range(K):
            n = npos[ib, jk]
            if n >= COHORT_THRESHOLD:
                ax.add_patch(Rectangle((K_RANGE[jk] - 0.5, ib - 0.5), 1, 1,
                                        facecolor="none", edgecolor="black",
                                        linewidth=1.4, zorder=3))
            elif n >= COHORT_ADVISORY:
                ax.add_patch(Rectangle((K_RANGE[jk] - 0.5, ib - 0.5), 1, 1,
                                        facecolor="none", edgecolor="0.4",
                                        linewidth=0.8, linestyle=":",
                                        zorder=3))


def _draw_partition_panel(ax, df: pd.DataFrame, contrast: str, label: str,
                          show_yticks: bool, show_xlabel: bool,
                          z_thresh: float):
    mean, frac, npos = _per_band_stats_k(df, contrast)
    cmap = _cmap_diverging()
    vmax = np.nanmax(np.abs(mean)) if np.isfinite(mean).any() else 1.0
    im = ax.imshow(
        mean, aspect="auto", cmap=cmap, vmin=-vmax, vmax=+vmax,
        extent=[K_RANGE[0] - 0.5, K_RANGE[-1] + 0.5,
                len(BRAIN_BANDS_NAMES) - 0.5, -0.5],
        interpolation="nearest", rasterized=True,
    )
    _draw_borders(ax, npos)
    Z_real, clusters = _real_clusters_t_z(df, contrast, z_thresh)
    sig_mass = _signflip_cluster_threshold(df, contrast, z_thresh)
    _hatch_cells(ax, clusters, sig_mass)
    if show_yticks:
        ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
                           fontsize=10)
    else:
        ax.set_yticks([])
    if show_xlabel:
        ax.set_xlabel("k (dendrogram cut)", fontsize=10)
    ax.set_title(label, fontsize=11, loc="left")
    return im, sig_mass


def _draw_h2c_strip(ax, h2c_raw: pd.DataFrame):
    """H2c summary strip: rows = bands, 2 columns = {task, learn}.

    Cell color = patient-mean ρ (Spearman of cophenetic-distance directional drift).
    Black border if ≥ 8/10 patients have ρ > 0.
    """
    B = len(BRAIN_BANDS_NAMES)
    cols = ["task", "learn"]
    mean = np.full((B, 2), np.nan)
    npos = np.full((B, 2), 0, dtype=int)
    for ib, band in enumerate(BRAIN_BANDS_NAMES):
        sub = h2c_raw[h2c_raw["band"] == band]
        for jc, col in enumerate(cols):
            x = sub[f"rho_{col}"].to_numpy()
            x = x[np.isfinite(x)]
            if x.size == 0:
                continue
            mean[ib, jc] = x.mean()
            npos[ib, jc] = int((x > 0).sum())
    cmap = _cmap_diverging()
    vmax = np.nanmax(np.abs(mean))
    ax.imshow(
        mean, aspect="auto", cmap=cmap, vmin=-vmax, vmax=+vmax,
        extent=[-0.5, 1.5, B - 0.5, -0.5],
        interpolation="nearest", rasterized=True,
    )
    for ib in range(B):
        for jc in range(2):
            n = npos[ib, jc]
            if n >= COHORT_THRESHOLD:
                ax.add_patch(Rectangle((jc - 0.5, ib - 0.5), 1, 1,
                                        facecolor="none", edgecolor="black",
                                        linewidth=1.4, zorder=3))
            elif n >= COHORT_ADVISORY:
                ax.add_patch(Rectangle((jc - 0.5, ib - 0.5), 1, 1,
                                        facecolor="none", edgecolor="0.4",
                                        linewidth=0.8, linestyle=":",
                                        zorder=3))
            ax.text(jc, ib, f"{mean[ib, jc]:.2f}",
                    ha="center", va="center", fontsize=8,
                    color="white" if abs(mean[ib, jc]) > 0.45 * vmax else "black")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["task_test", "task_learn"], fontsize=9)
    ax.set_yticks([])  # already labelled on the H2d panel to the left
    ax.set_title(r"$\rho_{\mathrm{H2c}}$ (no $k$ axis)", fontsize=11, loc="left")


def main() -> None:
    out_dir = FIGURES_ROOT / "section6"
    out_dir.mkdir(parents=True, exist_ok=True)

    p_part = REPORTS_ROOT / "imcoh_vi" / "h2_partition_multiscale_raw.csv"
    p_h2d  = REPORTS_ROOT / "imcoh_vi" / "h2d_persistence_raw.csv"
    p_h2c  = REPORTS_ROOT / "imcoh_vi" / "h2c_ultrametric_drift_raw.csv"

    df_part = pd.read_csv(p_part)
    df_h2d  = pd.read_csv(p_h2d).rename(columns={"delta_rho": "delta_rho"})
    df_h2c  = pd.read_csv(p_h2c)

    # Sanity asserts (era + cohort guard).
    pat_part = set(df_part["patient"].unique())
    pat_h2d  = set(df_h2d["patient"].unique())
    pat_h2c  = set(df_h2c["patient"].unique())
    expected = set(COHORT_N10)
    assert pat_part == expected, f"h2_partition cohort mismatch: {pat_part - expected} extra / {expected - pat_part} missing"
    assert pat_h2d  == expected, f"h2d cohort mismatch: {pat_h2d - expected} extra / {expected - pat_h2d} missing"
    assert pat_h2c  == expected, f"h2c cohort mismatch: {pat_h2c - expected} extra / {expected - pat_h2c} missing"

    # Cluster threshold for hatch — use t-style z = 1.65 (one-sided p<0.05).
    z_thresh = 1.65

    fig, axes = plt.subplots(
        2, 3, figsize=(15.0, 7.5), dpi=160,
        gridspec_kw={"hspace": 0.42, "wspace": 0.18,
                     "width_ratios": [1.0, 1.0, 0.55]},
    )
    sig_masses: dict[str, float] = {}
    last_im = None

    # Top row — three H2a partition contrasts: Δ_VI, Δ_H, Δ_NMI
    for ax, (col, label) in zip(axes[0], PARTITION_PANELS):
        im, sig = _draw_partition_panel(
            ax, df_part, col, label,
            show_yticks=(ax is axes[0, 0]),
            show_xlabel=False,
            z_thresh=z_thresh,
        )
        sig_masses[col] = sig
        last_im = im

    # Bottom-left: H2d block-pair coactivation Δρ(k)
    im_h2d, sig_h2d = _draw_partition_panel(
        axes[1, 0], df_h2d, "delta_rho",
        r"$\Delta\rho_{\mathrm{H2d}}(k)$ — block-pair persistence",
        show_yticks=True, show_xlabel=True, z_thresh=z_thresh,
    )
    sig_masses["delta_rho"] = sig_h2d

    # Bottom-middle: leave a placeholder — H2c band strip lives at (1,2)
    # to keep proportions clean.
    axes[1, 1].axis("off")

    # Bottom-right: H2c band-strip
    _draw_h2c_strip(axes[1, 2], df_h2c)

    # Single shared colorbar for the partition-style panels.
    cbar = fig.colorbar(last_im, ax=axes.ravel().tolist(),
                        shrink=0.65, pad=0.02, fraction=0.025)
    cbar.set_label("patient-mean contrast", fontsize=10)

    # Save PDF only (per CLAUDE.md never-PNG rule).
    out_pdf = out_dir / "task_trace_band_k_n10_imcoh_abs.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_pdf}")

    # Sidecar markdown — every cached row used + cluster-mass thresholds.
    n_pat_part = df_part["patient"].nunique()
    n_pat_h2d  = df_h2d["patient"].nunique()
    n_pat_h2c  = df_h2c["patient"].nunique()
    md = out_dir / "task_trace_band_k_n10_imcoh_abs.md"
    md_lines = [
        "---",
        "name: task-trace-band-k-n10",
        "type: figure-sidecar",
        "era: COHORT_N10",
        "status: current",
        "created: 2026-04-25",
        "updated: 2026-04-26",
        "pointers:",
        "  - .agents/reports/2026-04-25_task-trace-audit-and-recovery.md",
        "  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md",
        "---",
        "",
        "# Headline 5-panel band×k task-trace map (n=10, IMCOH_ABS)",
        "",
        "Surfaces existing per-(band, k) cohort evidence at n=10. The 5 panels",
        "are 5 operationalizations of the same residual H2a/H2d claim",
        "(rest_post − rest_pre alignment with task − rest_pre):",
        "",
        "- **Top:** three H2a partition-level contrasts — Δ_VI, Δ_H, Δ_NMI.",
        "- **Bottom-left:** Δρ_H2d block-pair coactivation persistence.",
        "- **Bottom-right:** H2c continuous ultrametric drift (band × {task, learn}, no k axis).",
        "",
        "**Δ_ARI is intentionally excluded** per `feedback_no_ari_in_figures.md`",
        "(2026-04-26 user decision); CSV column `d_ARI` retained for reference.",
        "",
        f"**Black border** = ≥ {COHORT_THRESHOLD}/{len(COHORT_N10)} patients sign-correct (cohort-wide).",
        f"**Dashed grey** = {COHORT_ADVISORY}–{COHORT_THRESHOLD-1}/{len(COHORT_N10)} (advisory). **Hatch** = sign-flip",
        f"cluster-permutation significant cluster (n_perm = {N_PERM}, α = 0.05).",
        "",
        "## Inputs",
        "",
        f"- `data/reports/imcoh_vi/h2_partition_multiscale_raw.csv` — {len(df_part)} rows, n_patients={n_pat_part}",
        f"- `data/reports/imcoh_vi/h2d_persistence_raw.csv` — {len(df_h2d)} rows, n_patients={n_pat_h2d}",
        f"- `data/reports/imcoh_vi/h2c_ultrametric_drift_raw.csv` — {len(df_h2c)} rows, n_patients={n_pat_h2c}",
        "",
        "## Cohort",
        "",
        f"`{', '.join(COHORT_N10)}` (n={len(COHORT_N10)}). Pat_03 1024 Hz outlier; Pat_14 vendor-replaced 2026-04-25.",
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
        "- **Black border** in a (band, k) cell → cohort-wide unanimity at that resolution.",
        "- **Hatched k-ranges** → cluster-permutation significant cluster.",
        "- **Compare across the 6 panels**: a cell that survives in 3+ panels is",
        "  triangulated evidence; a cell only in 1 panel is single-operationalization.",
        "- **k-axis extremes**: see companion `k_artefact_diagnostic_n10_imcoh_abs.pdf`",
        "  for `n_eff(k)`, singleton-fraction, max-cluster-fraction. The headline",
        "  figure does NOT mask extremes (per user choice 2026-04-25).",
        "",
        "## What this figure does NOT show",
        "",
        "- Per-patient breakdown — see `data/audit/trace_modules/Pat_NN_*.pdf`",
        "  (Stage 4 of audit-and-recovery).",
        "- Direction of subtree identity — see canonical reformalization regime",
        "  classifier in `2026-04-25_task-trace-canonical.md` (T regime).",
        "- Pat_03 in/out sensitivity — flagged as a known caveat; spot-check",
        "  per-band before citing any cell where Pat_03 is load-bearing.",
        "",
    ])
    md.write_text("\n".join(md_lines), encoding="utf-8")
    print(f"wrote {md}")


if __name__ == "__main__":
    main()
