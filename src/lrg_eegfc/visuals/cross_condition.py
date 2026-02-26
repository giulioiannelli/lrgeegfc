"""Cross-condition MSLCD diagnostic figures (Section 4.8).

Generates publication-quality figures comparing LRG multiscale community
detection diagnostics across patients, phases, and frequency bands.

Part B: 6 cross-condition figures (B1-B6)
Part C: Thresholding analysis figures (C1-C3)
Part D: LaTeX diagnostic table
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PHASE_LABELS

# ---------------------------------------------------------------------------
# Colour palettes (consistent across all figures)
# ---------------------------------------------------------------------------
PHASE_COLORS = {
    "rsPre": "#1f77b4", "taskLearn": "#ff7f0e",
    "taskTest": "#2ca02c", "rsPost": "#d62728",
}
BAND_COLORS = {
    "delta": "#1b9e77", "theta": "#d95f02", "alpha": "#7570b3",
    "beta": "#e7298a", "low_gamma": "#66a61e", "high_gamma": "#e6ab02",
}
PATIENT_MARKERS = {
    "Pat_02": "o", "Pat_03": "s", "Pat_05": "D",
    "Pat_07": "^", "Pat_08": "v",
}

_RC = {
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 14,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10,
    "text.usetex": False, "mathtext.fontset": "cm",
}

BANDS_ORDERED = list(BRAIN_BANDS_NAMES)
PHASES_ORDERED = list(PHASE_LABELS)

# Threshold percentiles for Part C
THRESHOLD_PERCENTILES = [0, 50, 70, 80, 85, 90, 95, 97, 99]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _save(fig, path: Path, dpi: int = 300):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    pdf_path = path.with_suffix(".pdf")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return path, pdf_path


def _text_color(value, vmin, vmax, cmap_obj):
    """Return 'white' or 'black' based on cell background luminance."""
    if np.isnan(value):
        return "black"
    norm = (value - vmin) / (vmax - vmin) if vmax > vmin else 0.5
    norm = np.clip(norm, 0.0, 1.0)
    rgba = cmap_obj(norm)
    lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "white" if lum < 0.5 else "black"


def _annotated_heatmap(ax, pivot, cmap_name, vmin=None, vmax=None,
                       fmt=".2f", fontsize=8):
    """Draw a heatmap with square cells and luminance-adaptive annotations."""
    data = pivot.values
    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)
    cmap_obj = cm.get_cmap(cmap_name)
    im = ax.imshow(data, aspect="equal", cmap=cmap_obj, vmin=vmin, vmax=vmax)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isnan(v):
                continue
            color = _text_color(v, vmin, vmax, cmap_obj)
            txt = f"{int(v)}" if fmt == "d" else f"{v:{fmt}}"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=fontsize, color=color)
    return im


def _setup_heatmap_axes(ax, bands=True, phases=True):
    """Standard tick setup for band-vs-phase heatmaps."""
    if phases:
        ax.set_xticks(range(len(PHASES_ORDERED)))
        ax.set_xticklabels(PHASES_ORDERED, rotation=45, ha="right")
    else:
        ax.set_xticks(range(len(PHASES_ORDERED)))
        ax.set_xticklabels([])
    if bands:
        ax.set_yticks(range(len(BANDS_ORDERED)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS_ORDERED])
    else:
        ax.set_yticks(range(len(BANDS_ORDERED)))
        ax.set_yticklabels([])


# ======================================================================
# B1: Coarsening trajectories -- beta band, all patients
# ======================================================================

def plot_coarsening_beta(df: pd.DataFrame,
                         trajectories_path: Path,
                         output_dir: Path,
                         dpi: int = 300) -> list[Path]:
    """n_max(tau) step-plots for beta band across all patients.

    NO title.  Common legend placed below the subplot grid.
    Consistent x-axis across subplots.
    """
    data = np.load(trajectories_path, allow_pickle=True)
    patients = sorted(df["patient"].unique())

    with plt.rc_context(_RC):
        ncols = min(len(patients), 3)
        nrows = int(np.ceil(len(patients) / ncols))
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(5.5 * ncols, 4 * nrows + 0.6),
            constrained_layout=True,
            sharex=True, sharey=True,
        )
        if len(patients) == 1:
            axes = np.array([axes])
        axes = np.atleast_2d(axes)

        lines_for_legend = {}

        for idx, patient in enumerate(patients):
            ax = axes.ravel()[idx]
            for phase in PHASES_ORDERED:
                key = f"{patient}__{phase}__beta"
                tau_key = key + "__tau"
                nmax_key = key + "__nmax"
                if tau_key not in data or nmax_key not in data:
                    continue
                tau = data[tau_key]
                nmax = data[nmax_key]
                line, = ax.step(np.log10(tau), nmax, where="mid", lw=1.8,
                                color=PHASE_COLORS[phase])
                if phase not in lines_for_legend:
                    lines_for_legend[phase] = line

            ax.text(0.04, 0.94, patient, transform=ax.transAxes,
                    fontsize=12, fontweight="bold", va="top")
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.3)
            if idx % ncols == 0:
                ax.set_ylabel(r"$n_{\max}(\tau)$")
            if idx >= (nrows - 1) * ncols:
                ax.set_xlabel(r"$\log_{10}\,\tau$")

        for idx in range(len(patients), nrows * ncols):
            axes.ravel()[idx].set_visible(False)

        if lines_for_legend:
            fig.legend(
                lines_for_legend.values(), lines_for_legend.keys(),
                loc="lower center", ncol=len(lines_for_legend),
                bbox_to_anchor=(0.5, -0.01), fontsize=11,
                frameon=True, fancybox=True,
            )

        p, _ = _save(fig, output_dir / "fig_B1_coarsening_beta.png", dpi)
    return [p]


# ======================================================================
# B2: Coarsening trajectories -- Pat_02, all bands
# ======================================================================

def plot_coarsening_patient(df: pd.DataFrame,
                            trajectories_path: Path,
                            output_dir: Path,
                            patient: str = "Pat_02",
                            dpi: int = 300) -> list[Path]:
    """n_max(tau) step-plots for one patient across all 6 bands.

    NO title.  Common legend below grid.  Consistent axes.
    """
    data = np.load(trajectories_path, allow_pickle=True)

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(
            2, 3, figsize=(16, 9.6),
            constrained_layout=True,
            sharex=True, sharey=True,
        )

        lines_for_legend = {}

        for k, band in enumerate(BANDS_ORDERED):
            ax = axes.ravel()[k]
            for phase in PHASES_ORDERED:
                key = f"{patient}__{phase}__{band}"
                tau_key = key + "__tau"
                nmax_key = key + "__nmax"
                if tau_key not in data or nmax_key not in data:
                    continue
                tau = data[tau_key]
                nmax = data[nmax_key]
                line, = ax.step(np.log10(tau), nmax, where="mid", lw=1.8,
                                color=PHASE_COLORS[phase])
                if phase not in lines_for_legend:
                    lines_for_legend[phase] = line

            band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
            ax.text(0.04, 0.94, band_tex, transform=ax.transAxes,
                    fontsize=13, fontweight="bold", va="top")
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.3)
            if k >= 3:
                ax.set_xlabel(r"$\log_{10}\,\tau$")
            if k % 3 == 0:
                ax.set_ylabel(r"$n_{\max}(\tau)$")

        if lines_for_legend:
            fig.legend(
                lines_for_legend.values(), lines_for_legend.keys(),
                loc="lower center", ncol=len(lines_for_legend),
                bbox_to_anchor=(0.5, -0.01), fontsize=11,
                frameon=True, fancybox=True,
            )

        p, _ = _save(fig, output_dir / f"fig_B2_coarsening_{patient}.png", dpi)
    return [p]


# ======================================================================
# B3: Metastability heatmaps (mu_mean, mu_max only)
# ======================================================================

def plot_metastability_heatmaps(df: pd.DataFrame,
                                output_dir: Path,
                                dpi: int = 300) -> list[Path]:
    r"""2-row heatmap grid: :math:`\bar\mu` and :math:`\mu_{\max}`.

    Square cells, sequential colormap, luminance-adaptive text.
    """
    patients = sorted(df["patient"].unique())
    n_pat = len(patients)
    stats = ["mu_mean", "mu_max"]
    stat_labels = [r"$\bar\mu$", r"$\mu_{\max}$"]
    cmap_name = "magma"

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(
            len(stats), n_pat,
            figsize=(3.5 * n_pat + 1.2, 3.8 * len(stats)),
            constrained_layout=True,
        )
        if n_pat == 1:
            axes = axes[:, np.newaxis]

        for row, (stat, label) in enumerate(zip(stats, stat_labels)):
            all_vals = df[stat].dropna().values
            vmin = float(all_vals.min()) if len(all_vals) else 0
            vmax = float(all_vals.max()) if len(all_vals) else 1

            for col, patient in enumerate(patients):
                ax = axes[row, col]
                sub = df[df["patient"] == patient]
                pivot = sub.pivot_table(
                    index="band", columns="phase",
                    values=stat, aggfunc="first",
                )
                pivot = pivot.reindex(index=BANDS_ORDERED, columns=PHASES_ORDERED)
                im = _annotated_heatmap(ax, pivot, cmap_name,
                                        vmin=vmin, vmax=vmax)
                _setup_heatmap_axes(
                    ax,
                    bands=(col == 0),
                    phases=(row == len(stats) - 1),
                )
                if row == 0:
                    ax.text(0.5, 1.08, patient, transform=ax.transAxes,
                            ha="center", fontsize=12, fontweight="bold")

            cb = fig.colorbar(im, ax=axes[row, :].tolist(),
                              fraction=0.02, pad=0.03)
            cb.set_label(label, fontsize=12)

        p, _ = _save(fig, output_dir / "fig_B3_metastability_heatmaps.png", dpi)
    return [p]


# ======================================================================
# B4: Metastability boxplots
# ======================================================================

def plot_metastability_boxplots(df: pd.DataFrame,
                                output_dir: Path,
                                dpi: int = 300) -> list[Path]:
    r"""Boxplots of :math:`\bar\mu` by phase, one subplot per band.

    Individual patient values overlaid as shaped markers.
    Legend outside the subplot grid explains both phase colours and
    patient marker shapes.
    """
    patients = sorted(df["patient"].unique())

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(
            1, len(BANDS_ORDERED),
            figsize=(3 * len(BANDS_ORDERED) + 2, 5),
            constrained_layout=True, sharey=True,
        )

        for k, band in enumerate(BANDS_ORDERED):
            ax = axes[k]
            band_data = df[df["band"] == band]
            positions = []
            data_list = []
            colors_list = []

            for j, phase in enumerate(PHASES_ORDERED):
                vals = band_data[
                    band_data["phase"] == phase
                ]["mu_mean"].dropna().values
                if len(vals) > 0:
                    data_list.append(vals)
                    positions.append(j)
                    colors_list.append(PHASE_COLORS[phase])

            if data_list:
                bp = ax.boxplot(
                    data_list, positions=positions, widths=0.55,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=2),
                )
                for patch, color in zip(bp["boxes"], colors_list):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.5)

            # Overlay individual patient values
            for j, phase in enumerate(PHASES_ORDERED):
                for patient in patients:
                    row = band_data[
                        (band_data["phase"] == phase)
                        & (band_data["patient"] == patient)
                    ]
                    if len(row) == 0 or pd.isna(row["mu_mean"].values[0]):
                        continue
                    ax.scatter(
                        j, row["mu_mean"].values[0],
                        marker=PATIENT_MARKERS.get(patient, "o"),
                        color=PHASE_COLORS[phase], s=50, zorder=5,
                        edgecolors="black", linewidths=0.5,
                    )

            ax.set_xticks(range(len(PHASES_ORDERED)))
            ax.set_xticklabels(PHASES_ORDERED, rotation=45,
                               ha="right", fontsize=9)
            band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
            ax.text(0.5, 1.02, band_tex, transform=ax.transAxes,
                    ha="center", fontsize=11, fontweight="bold")
            ax.grid(axis="y", alpha=0.3)
            if k == 0:
                ax.set_ylabel(r"$\bar\mu$")

        # Combined legend: phase colours + patient markers
        phase_handles = [
            Line2D([0], [0], color=PHASE_COLORS[p], lw=8, alpha=0.5, label=p)
            for p in PHASES_ORDERED
        ]
        patient_handles = [
            Line2D(
                [0], [0], marker=PATIENT_MARKERS[p], color="w",
                markerfacecolor="gray", markersize=8,
                markeredgecolor="black", markeredgewidth=0.5, label=p,
            )
            for p in patients if p in PATIENT_MARKERS
        ]
        all_handles = phase_handles + patient_handles
        all_labels = (
            list(PHASES_ORDERED)
            + [p for p in patients if p in PATIENT_MARKERS]
        )
        fig.legend(
            all_handles, all_labels,
            loc="center right", bbox_to_anchor=(1.13, 0.5),
            fontsize=9, frameon=True, title="Phase / Patient",
        )

        p, _ = _save(fig, output_dir / "fig_B4_metastability_boxplots.png", dpi)
    return [p]


# ======================================================================
# B5: Partition richness heatmap
# ======================================================================

def plot_partition_richness(df: pd.DataFrame,
                            output_dir: Path,
                            dpi: int = 300) -> list[Path]:
    """N_sens heatmap per patient.  Square cells, luminance-adaptive text."""
    patients = sorted(df["patient"].unique())
    n_pat = len(patients)
    cmap_name = "YlGnBu"

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(
            1, n_pat,
            figsize=(3.5 * n_pat + 1.2, 4.5),
            constrained_layout=True,
        )
        if n_pat == 1:
            axes = [axes]

        all_vals = df["N_sens"].dropna().values
        vmin = float(all_vals.min()) if len(all_vals) else 0
        vmax = float(all_vals.max()) if len(all_vals) else 1

        for idx, patient in enumerate(patients):
            ax = axes[idx]
            sub = df[df["patient"] == patient]
            pivot = sub.pivot_table(
                index="band", columns="phase",
                values="N_sens", aggfunc="first",
            )
            pivot = pivot.reindex(index=BANDS_ORDERED, columns=PHASES_ORDERED)
            im = _annotated_heatmap(ax, pivot, cmap_name,
                                    vmin=vmin, vmax=vmax, fmt="d")
            _setup_heatmap_axes(ax, bands=(idx == 0), phases=True)
            ax.text(0.5, 1.05, patient, transform=ax.transAxes,
                    ha="center", fontsize=12, fontweight="bold")

        cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.03)
        cb.set_label(r"$N_{\mathrm{sens}}$", fontsize=12)

        p, _ = _save(fig, output_dir / "fig_B5_partition_richness.png", dpi)
    return [p]


# ======================================================================
# B6: Community size balance heatmap
# ======================================================================

def plot_community_balance(df: pd.DataFrame,
                           output_dir: Path,
                           dpi: int = 300) -> list[Path]:
    """H_size heatmap per patient.  Square cells, luminance-adaptive text."""
    patients = sorted(df["patient"].unique())
    n_pat = len(patients)
    cmap_name = "RdYlGn"

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(
            1, n_pat,
            figsize=(3.5 * n_pat + 1.2, 4.5),
            constrained_layout=True,
        )
        if n_pat == 1:
            axes = [axes]

        for idx, patient in enumerate(patients):
            ax = axes[idx]
            sub = df[df["patient"] == patient]
            pivot = sub.pivot_table(
                index="band", columns="phase",
                values="H_size", aggfunc="first",
            )
            pivot = pivot.reindex(index=BANDS_ORDERED, columns=PHASES_ORDERED)
            im = _annotated_heatmap(ax, pivot, cmap_name,
                                    vmin=0.0, vmax=1.0)
            _setup_heatmap_axes(ax, bands=(idx == 0), phases=True)
            ax.text(0.5, 1.05, patient, transform=ax.transAxes,
                    ha="center", fontsize=12, fontweight="bold")

        cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.03)
        cb.set_label(r"$H_{\mathrm{size}}$", fontsize=12)

        p, _ = _save(fig, output_dir / "fig_B6_community_balance.png", dpi)
    return [p]


# ======================================================================
# Part C: Threshold analysis figures (C1, C2, C3)
# ======================================================================

def _load_threshold_data(threshold_path: Path) -> dict:
    """Load precomputed threshold analysis NPZ."""
    raw = np.load(threshold_path, allow_pickle=True)
    return dict(raw)


def plot_threshold_susceptibility(threshold_path: Path,
                                  output_dir: Path,
                                  dpi: int = 300) -> list[Path]:
    """C1: Entropic susceptibility C(tau) at different edge-weight thresholds.

    2x3 grid, one subplot per band.  Curves coloured by threshold percentile.
    """
    data = _load_threshold_data(threshold_path)
    percentiles = data["percentiles"]
    saved = []

    # Colour map for thresholds: light -> dark
    thresh_cmap = cm.get_cmap("plasma", len(percentiles))

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(
            2, 3, figsize=(16, 10),
            constrained_layout=True,
        )

        for k, band in enumerate(BANDS_ORDERED):
            ax = axes.ravel()[k]
            tau_key = f"{band}__tau"
            if tau_key not in data:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes)
                continue

            tau = data[tau_key]

            for ti, pct in enumerate(percentiles):
                c_key = f"{band}__C__{pct}"
                if c_key not in data:
                    continue
                C = data[c_key]
                ax.plot(tau, C, lw=1.5,
                        color=thresh_cmap(ti / max(len(percentiles) - 1, 1)),
                        label=f"{pct}%")

            ax.set_xscale("log")
            band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
            ax.text(0.04, 0.94, band_tex, transform=ax.transAxes,
                    fontsize=13, fontweight="bold", va="top")
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.3, which="both")
            if k >= 3:
                ax.set_xlabel(r"$\tau$")
            if k % 3 == 0:
                ax.set_ylabel(r"$C(\tau)$")

        # Shared colorbar for threshold percentile
        sm = cm.ScalarMappable(
            cmap=thresh_cmap,
            norm=plt.Normalize(vmin=percentiles[0], vmax=percentiles[-1]),
        )
        sm.set_array([])
        cb = fig.colorbar(sm, ax=axes.ravel().tolist(),
                          fraction=0.02, pad=0.02)
        cb.set_label("Edge-weight percentile threshold (%)", fontsize=11)

        p, _ = _save(
            fig, output_dir / "fig_C1_threshold_susceptibility.png", dpi)
        saved.append(p)

    return saved


def plot_threshold_spectrum(threshold_path: Path,
                            output_dir: Path,
                            dpi: int = 300) -> list[Path]:
    """C2: Laplacian eigenvalue spectrum at different thresholds.

    Three panels:
    (a) Sorted eigenvalues for a representative band at key thresholds
    (b) lambda_2 vs threshold for all bands
    (c) Spectral gap ratio (lambda_2 / lambda_max) vs threshold
    """
    data = _load_threshold_data(threshold_path)
    percentiles = data["percentiles"]
    saved = []

    # Key thresholds to highlight in panel (a)
    highlight_pcts = [0, 70, 90, 95, 99]
    ref_band = "beta"

    thresh_cmap = cm.get_cmap("plasma", len(highlight_pcts))

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(
            1, 3, figsize=(18, 5),
            constrained_layout=True,
        )

        # (a) Sorted eigenvalues for reference band
        ax = axes[0]
        for ti, pct in enumerate(highlight_pcts):
            eig_key = f"{ref_band}__eigenvalues__{pct}"
            if eig_key not in data:
                continue
            eig = data[eig_key]
            ax.plot(
                range(len(eig)), eig, lw=1.2,
                color=thresh_cmap(ti / max(len(highlight_pcts) - 1, 1)),
                label=f"{pct}%",
            )
        ax.set_yscale("log")
        ax.set_xlabel("Eigenvalue index")
        ax.set_ylabel(r"$\lambda_i$")
        ax.text(0.04, 0.94, r"(a) $\beta$ eigenvalue spectrum",
                transform=ax.transAxes, fontsize=11, va="top")
        ax.legend(fontsize=8, title="Threshold")
        ax.grid(alpha=0.3)

        # x-axis: edges remaining (%) on log scale
        edges_remaining = np.maximum(100.0 - percentiles, 0.5)

        # (b) lambda_2 (LCC) vs threshold for all bands
        ax = axes[1]
        for band in BANDS_ORDERED:
            lam2_vals = []
            for pct in percentiles:
                eig_key = f"{band}__eigenvalues__{pct}"
                if eig_key not in data:
                    lam2_vals.append(np.nan)
                    continue
                eig = data[eig_key]
                lam2_vals.append(float(eig[1]))
            ax.plot(edges_remaining, lam2_vals, "o-", lw=1.5, ms=4,
                    color=BAND_COLORS[band],
                    label=BRAIN_BAND_TEX_DICT.get(band, band))
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_yscale("log")
        ax.set_xlabel("Edges remaining (%)")
        ax.set_ylabel(r"$\lambda_2$ (LCC)")
        ax.text(0.04, 0.94, r"(b) algebraic connectivity",
                transform=ax.transAxes, fontsize=11, va="top")
        ax.legend(fontsize=7, ncol=2)
        ax.grid(alpha=0.3, which="both")

        # (c) Spectral gap ratio vs threshold
        ax = axes[2]
        for band in BANDS_ORDERED:
            gap_vals = []
            for pct in percentiles:
                eig_key = f"{band}__eigenvalues__{pct}"
                if eig_key not in data:
                    gap_vals.append(np.nan)
                    continue
                eig = data[eig_key]
                lam2 = float(eig[1])
                lam_max = float(eig[-1])
                gap_vals.append(lam2 / lam_max if lam_max > 0 else np.nan)
            ax.plot(edges_remaining, gap_vals, "o-", lw=1.5, ms=4,
                    color=BAND_COLORS[band],
                    label=BRAIN_BAND_TEX_DICT.get(band, band))
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xlabel("Edges remaining (%)")
        ax.set_ylabel(r"$\lambda_2 / \lambda_{\max}$")
        ax.text(0.04, 0.94, r"(c) spectral gap ratio",
                transform=ax.transAxes, fontsize=11, va="top")
        ax.legend(fontsize=7, ncol=2)
        ax.grid(alpha=0.3, which="both")

        p, _ = _save(
            fig, output_dir / "fig_C2_threshold_spectrum.png", dpi)
        saved.append(p)

    return saved


def plot_threshold_connectivity(threshold_path: Path,
                                output_dir: Path,
                                dpi: int = 300) -> list[Path]:
    """C2: Network connectivity metrics vs edge-weight threshold.

    2x2 grid: connected components, LCC fraction, LCC density, transitivity.
    One curve per band.
    """
    data = _load_threshold_data(threshold_path)
    percentiles = data["percentiles"]
    saved = []

    metric_labels = [
        (0, r"Connected components"),
        (1, r"LCC fraction $|C_1|/N$"),
        (2, r"LCC density"),
        (3, r"Transitivity"),
    ]

    # x-axis: edges remaining (%) on log scale
    edges_remaining = np.maximum(100.0 - percentiles, 0.5)

    with plt.rc_context(_RC):
        fig, axes = plt.subplots(
            2, 2, figsize=(12, 9),
            constrained_layout=True,
        )

        for ax, (col_idx, ylabel) in zip(axes.ravel(), metric_labels):
            for band in BANDS_ORDERED:
                conn_key = f"{band}__connectivity"
                if conn_key not in data:
                    continue
                conn = data[conn_key]  # shape (n_thresholds, 4)
                ax.plot(
                    edges_remaining, conn[:, col_idx], "o-", lw=1.5, ms=4,
                    color=BAND_COLORS[band],
                    label=BRAIN_BAND_TEX_DICT.get(band, band),
                )
            ax.set_xscale("log")
            ax.invert_xaxis()
            ax.set_xlabel("Edges remaining (%)")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.3, which="both")

        # Single legend outside
        handles, labels = axes.ravel()[0].get_legend_handles_labels()
        fig.legend(
            handles, labels,
            loc="lower center", ncol=len(BANDS_ORDERED),
            bbox_to_anchor=(0.5, -0.02), fontsize=10,
            frameon=True, fancybox=True,
        )

        p, _ = _save(
            fig, output_dir / "fig_C2_threshold_connectivity.png", dpi)
        saved.append(p)

    return saved


# ======================================================================
# Part D: LaTeX table (trimmed columns)
# ======================================================================

def generate_latex_table(df: pd.DataFrame,
                         output_path: Path,
                         patient: str = "Pat_02") -> Path:
    """Generate booktabs LaTeX table with trimmed diagnostics.

    Columns: N_sens, n_max(tau'), n_max(1/lambda_2), mu_mean, mu_max, H_size.
    """
    sub = df[df["patient"] == patient].copy()
    if sub.empty:
        raise ValueError(f"No data for {patient}")

    cols = ["N_sens", "n_start", "n_end", "mu_mean", "mu_max", "H_size"]
    col_labels = [
        r"$N_{\text{sens}}$",
        r"$n_{\max}(\tau')$",
        r"$n_{\max}(1/\lambda_2)$",
        r"$\bar\mu$",
        r"$\mu_{\max}$",
        r"$H_{\text{size}}$",
    ]

    lines = [
        "% Auto-generated by lrg-eegfc compute diagnostics",
        r"\begin{table}[htbp]",
        r"\centering",
        f"\\caption{{MSLCD diagnostics for {patient.replace('_', chr(92) + '_')}.}}",
        f"\\label{{tab:mslcd-diagnostics-{patient.lower().replace('_', '')}}}",
        "\\begin{tabular}{ll" + "r" * len(cols) + "}",
        r"\toprule",
        "Band & Phase & " + " & ".join(col_labels) + r" \\",
        r"\midrule",
    ]

    for band in BANDS_ORDERED:
        band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
        first = True
        for phase in PHASES_ORDERED:
            row = sub[(sub["band"] == band) & (sub["phase"] == phase)]
            if len(row) == 0:
                vals = ["--"] * len(cols)
            else:
                r = row.iloc[0]
                vals = []
                for c in cols:
                    v = r[c]
                    if pd.isna(v):
                        vals.append("--")
                    elif c in ("N_sens", "n_start", "n_end"):
                        vals.append(f"{int(v)}")
                    else:
                        vals.append(f"{v:.3f}")
            band_cell = (
                f"\\multirow{{4}}{{*}}{{{band_tex}}}" if first else ""
            )
            lines.append(
                f"{band_cell} & {phase} & " + " & ".join(vals) + r" \\"
            )
            first = False
        lines.append(r"\midrule")

    # Replace last midrule with bottomrule
    lines[-1] = r"\bottomrule"
    lines += [
        r"\end{tabular}",
        r"\end{table}",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines))
    return output_path


# ======================================================================
# Master dispatcher
# ======================================================================

FIGURE_TYPES = {
    "coarsening-beta": "B1 -- Coarsening trajectories (beta, all patients)",
    "coarsening-patient": "B2 -- Coarsening trajectories (Pat_02, all bands)",
    "metastability-heatmaps": "B3 -- Metastability heatmaps",
    "metastability-boxplots": "B4 -- Metastability boxplots",
    "partition-richness": "B5 -- Partition richness heatmap",
    "community-balance": "B6 -- Community size balance heatmap",
    "threshold-susceptibility": "C1 -- C(tau) at different thresholds",
    "threshold-connectivity": "C2 -- Connectivity metrics vs threshold",
    "latex-table": "Part D -- LaTeX diagnostic table",
    "all": "Generate all figures and table",
}


def generate_cross_condition_figure(
    figure_type: str,
    csv_path: Path = Path("data/tables/mslcd_diagnostics_master.csv"),
    trajectories_path: Path = Path(
        "data/tables/mslcd_coarsening_trajectories.npz"
    ),
    threshold_path: Path = Path("data/tables/threshold_analysis_Pat02.npz"),
    output_dir: Path = Path(
        "data/figures/report_mslcd_section/cross_condition"
    ),
    threshold_output_dir: Path = Path(
        "data/figures/report_mslcd_section/threshold_analysis"
    ),
    latex_output: Path = Path("data/tables/mslcd_diagnostics_Pat02.tex"),
    dpi: int = 300,
) -> list[Path]:
    """Generate one or all cross-condition figures.

    Returns list of saved file paths.
    """
    saved: list[Path] = []

    # Determine which figure types to generate
    non_all = [k for k in FIGURE_TYPES if k != "all"]
    types_to_run = non_all if figure_type == "all" else [figure_type]

    # Load CSV only once (needed by B1-B6 and latex-table)
    df = None
    needs_csv = {
        "coarsening-beta", "coarsening-patient",
        "metastability-heatmaps", "metastability-boxplots",
        "partition-richness", "community-balance", "latex-table",
    }
    if needs_csv & set(types_to_run):
        df = pd.read_csv(csv_path)

    for ft in types_to_run:
        if ft == "coarsening-beta":
            saved.extend(
                plot_coarsening_beta(df, trajectories_path, output_dir, dpi))
        elif ft == "coarsening-patient":
            saved.extend(
                plot_coarsening_patient(
                    df, trajectories_path, output_dir, dpi=dpi))
        elif ft == "metastability-heatmaps":
            saved.extend(
                plot_metastability_heatmaps(df, output_dir, dpi))
        elif ft == "metastability-boxplots":
            saved.extend(
                plot_metastability_boxplots(df, output_dir, dpi))
        elif ft == "partition-richness":
            saved.extend(
                plot_partition_richness(df, output_dir, dpi))
        elif ft == "community-balance":
            saved.extend(
                plot_community_balance(df, output_dir, dpi))
        elif ft == "threshold-susceptibility":
            saved.extend(
                plot_threshold_susceptibility(
                    threshold_path, threshold_output_dir, dpi))
        elif ft == "threshold-connectivity":
            saved.extend(
                plot_threshold_connectivity(
                    threshold_path, threshold_output_dir, dpi))
        elif ft == "latex-table":
            saved.append(generate_latex_table(df, latex_output))

    return saved
