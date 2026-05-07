#!/usr/bin/env python3
"""Section 2 figure C2 (spectral distribution per patient):
reproduction of the original MSC fig_msc_distribution_all_patients_rsPre.pdf
style for both MSC and ImCoh, side-by-side.

For each patient, shows the distribution across all off-diagonal pairs of
the coherence metric as a function of frequency:

  - black line   : median C_{ij}(f)
  - dark band    : inter-quartile range (25th-75th percentile)
  - light band   : 5th-95th percentile range
  - color bands  : brain-band shaded backgrounds (δ, θ, α, β, γₗ, γₕ)
  - per-band bar : median in that band, in the band's color

Produced for multiple phases (rest_pre, task_learn, rest_post).

Run:
  python scripts/10_notes_imcoh/fig_C2_spectral_distribution.py [-v]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt

from _shared import (
    ALL_PATIENTS, BANDS, BRAIN_BAND_TEX_DICT,
    apply_pub_style, save_fig, SECTION2_ROOT,
    set_memory_limit, checkpoint_cleanup,
    load_channel_labels,
)
from lrg_eegfc.config.const import BRAIN_BANDS, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH, CACHE_ROOT
from lrg_eegfc.utils.io.patient import load_mat_pat_data
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch

#: Where the per-patient percentile curves are cached so re-runs are instant.
PCTL_CACHE_DIR = CACHE_ROOT / "spectral_pctl"


#: Band background colours (soft versions of BAND_COLORS for shading).
BAND_BG_COLOR = {
    "delta":      "#4C72B0",
    "theta":      "#55A868",
    "alpha":      "#C44E52",
    "beta":       "#8172B3",
    "low_gamma":  "#CCB974",
    "high_gamma": "#64B5CD",
}

METHOD_LABELS = {"msc": "MSC", "imcoh": r"$|\mathrm{ImCoh}|$"}
#: Y-axis labels per method (use |ImCoh| since imaginary coherence can be negative)
Y_AXIS_LABEL = {
    "msc": r"MSC  $C_{ij}(f)$",
    "imcoh": r"|ImCoh|  $|\mathrm{Im}\,C_{ij}(f)|$",
}


def _extract_timeseries(patient: str, phase: str):
    """Return (X [N, L], fs) for a patient/phase, or (None, None) if missing."""
    try:
        mat_data = load_mat_pat_data(patient, phase, SEEG_DATAPATH)
    except (FileNotFoundError, OSError):
        return None, None
    X = None
    for k in mat_data:
        if not k.startswith('_'):
            v = mat_data[k]
            if hasattr(v, 'shape') and len(v.shape) == 2 and min(v.shape) > 10:
                X = np.array(v, dtype=np.float64)
                break
    if X is None:
        return None, None
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = 2048.0 if patient != "Pat_03" else 1024.0
    return X, fs


def _compute_or_load_percentiles(patient: str, phase: str, method: str):
    """Return (freqs, pctls_dict, N, n_pairs) using a disk cache.

    pctls_dict keys: 'p5', 'p25', 'p50', 'p75', 'p95' and 'band_med_<band>'.
    Cache key: one .npz per (patient, phase, method).
    When a file exists it is loaded directly — no Welch recomputation.
    """
    PCTL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = PCTL_CACHE_DIR / f"{patient}_{phase}_{method}.npz"
    if cache_path.exists():
        data = np.load(cache_path)
        return (data["freqs"],
                {k: data[k] for k in data.files if k != "freqs"
                 and k not in ("N", "n_pairs")},
                int(data["N"]), int(data["n_pairs"]))

    # Compute from raw timeseries
    X, fs = _extract_timeseries(patient, phase)
    if X is None:
        return None, None, 0, 0
    nperseg = nperseg_for_fs(fs)
    freqs, Coh = compute_msc_welch(X, fs, nperseg=nperseg, metric=method)

    N = Coh.shape[0]
    iu = np.triu_indices(N, k=1)
    vals = np.abs(Coh[iu[0], iu[1], :])
    p5, p25, p50, p75, p95 = np.percentile(vals, [5, 25, 50, 75, 95], axis=0)

    # Per-band median (over all pair-frequency samples inside the band)
    band_medians = {}
    for band_name, (flo, fhi) in BRAIN_BANDS.items():
        mask = (freqs >= flo) & (freqs <= fhi)
        band_medians[f"band_med_{band_name}"] = (
            float(np.median(vals[:, mask])) if mask.any() else np.nan
        )

    n_pairs = N * (N - 1) // 2
    np.savez_compressed(
        cache_path,
        freqs=freqs, p5=p5, p25=p25, p50=p50, p75=p75, p95=p95,
        N=N, n_pairs=n_pairs,
        **band_medians,
    )

    # Free the tensor immediately
    del X, Coh, vals
    import gc; gc.collect()

    pctls = {"p5": p5, "p25": p25, "p50": p50, "p75": p75, "p95": p95,
             **band_medians}
    return freqs, pctls, N, n_pairs


def _plot_patient_panel(ax, freqs, pctls, N, n_pairs, patient, method,
                        show_ylabel: bool, show_xlabel: bool,
                        show_legend: bool):
    """Draw one patient panel from precomputed percentiles."""
    p5, p25, p50, p75, p95 = (pctls["p5"], pctls["p25"], pctls["p50"],
                               pctls["p75"], pctls["p95"])

    # Background band shading
    for band_name, (flo, fhi) in BRAIN_BANDS.items():
        ax.axvspan(flo, fhi, color=BAND_BG_COLOR[band_name], alpha=0.10, zorder=0)
        # Band tex label at top
        ax.text(np.sqrt(flo * fhi), 0.97, BRAIN_BAND_TEX_DICT[band_name],
                transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=10,
                color=BAND_BG_COLOR[band_name], fontweight="bold")

    # 5th-95th percentile range
    ax.fill_between(freqs, p5, p95, color="0.80", alpha=0.55,
                    linewidth=0, zorder=1, label="5th-95th pctl")
    # IQR
    ax.fill_between(freqs, p25, p75, color="0.55", alpha=0.65,
                    linewidth=0, zorder=2, label="IQR")
    # Median line
    ax.plot(freqs, p50, color="black", lw=1.2, zorder=3, label="Median")

    # Per-band median horizontal bar in the band's color
    for band_name, (flo, fhi) in BRAIN_BANDS.items():
        band_med = pctls.get(f"band_med_{band_name}", np.nan)
        if np.isfinite(band_med):
            ax.hlines(band_med, flo, fhi,
                      colors=BAND_BG_COLOR[band_name], lw=2.2, zorder=4)

    ax.set_xscale("log")
    ax.set_xticks([1, 4, 8, 13, 30, 80, 300])
    ax.set_xticklabels([1, 4, 8, 13, 30, 80, 300])
    ax.set_xlim(freqs[1], 300)
    # Log-y for both methods: coherence distributions are heavily
    # right-skewed and a linear y-axis flattens the bulk into a sliver.
    ax.set_yscale("log")
    p5_pos = p5[p5 > 0]
    if p5_pos.size > 0:
        ax.set_ylim(bottom=max(p5_pos.min() * 0.5, 1e-4))
    if show_ylabel:
        ax.set_ylabel(Y_AXIS_LABEL[method], fontsize=11)
    if show_xlabel:
        ax.set_xlabel("Frequency (Hz)", fontsize=11)
    ax.set_title(f"{patient} (N={N}, {n_pairs} pairs)",
                 fontsize=11, fontweight="bold")
    if show_legend:
        ax.legend(loc="upper left", fontsize=7, framealpha=0.8)


def fig_all_patients(
    method: str,
    phase: str,
    patients: list[str],
    output_dir: Path,
    verbose: bool = False,
):
    """Replicate the fig_msc_distribution style for the given method & phase."""
    n = len(patients)
    cols = 3
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.2 * rows))
    axes = axes.flatten()

    for idx, patient in enumerate(patients):
        ax = axes[idx]
        freqs, pctls, N, n_pairs = _compute_or_load_percentiles(
            patient, phase, method,
        )
        if freqs is None:
            ax.text(0.5, 0.5, f"{patient}\nN/A", ha="center", va="center",
                    transform=ax.transAxes)
            ax.axis("off")
            continue

        show_legend = (idx == 0)
        show_ylabel = (idx % cols == 0)
        show_xlabel = (idx // cols == rows - 1)

        _plot_patient_panel(ax, freqs, pctls, N, n_pairs, patient, method,
                            show_ylabel, show_xlabel, show_legend)

        if verbose:
            print(f"    {patient}: done")

    for i in range(n, rows * cols):
        axes[i].axis("off")

    label = METHOD_LABELS[method]
    fig.tight_layout()
    save_fig(fig,
             output_dir / f"fig_C2_spectral_distribution_{method}_{phase}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=SECTION2_ROOT)
    parser.add_argument("--patients", nargs="+", default=ALL_PATIENTS)
    parser.add_argument("--phases", nargs="+",
                        default=["rest_pre", "task_learn", "rest_post"])
    parser.add_argument("--methods", nargs="+", default=["msc", "imcoh"])
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()
    set_memory_limit(12.0)
    out = args.output_dir / "fig_C"

    for method in args.methods:
        for phase in args.phases:
            print(f"--- {method.upper()} distribution, {phase} ---")
            fig_all_patients(method, phase, args.patients, out, args.verbose)
            checkpoint_cleanup(f"{method} {phase}")

    print("\nDone.")


if __name__ == "__main__":
    main()
