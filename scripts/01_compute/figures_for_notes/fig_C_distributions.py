#!/usr/bin/env python3
"""Section 2 Figure C: weight distributions (MSC vs ImCoh).

C1 — Per-band same/cross-probe distributions for Pat_02, rest_pre.
C2 — Raw spectral example: coherence(f) for one same-probe and one cross-probe pair.
C3 — KDE overlay across all patients for beta band.

Run:
  python scripts/10_notes_imcoh/fig_C_distributions.py [-v]
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
from scipy.stats import gaussian_kde

from _shared import (
    ALL_PATIENTS, BANDS, BAND_TEX, BRAIN_BAND_TEX_DICT,
    CLR_SAME, CLR_CROSS, CLR_MSC, CLR_IMCOH, PATIENT_COLORS,
    REPR_PATIENTS, REPR_BANDS, REPR_PHASES,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_fig, SECTION2_ROOT,
    build_probe_mask, probe_weight_distributions,
    set_memory_limit, checkpoint_cleanup,
    SECTION2_METHODS as METHODS,
    SECTION2_METHOD_LABELS as METHOD_LABELS
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import BRAIN_BANDS, nperseg_for_fs



# ---------------------------------------------------------------------------
# Figure C1 — Per-band same/cross-probe distributions
# ---------------------------------------------------------------------------

def fig_c1_weight_distributions(
    patient: str = "Pat_02",
    phase: str = "rest_pre",
    output_dir: Path = SECTION2_ROOT / "fig_C",
    verbose: bool = False,
):
    """2 subfigures (MSC, ImCoh), each with 6 band panels showing
    same-probe (orange) vs cross-probe (blue) histograms."""

    for method in METHODS:
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        axes_flat = axes.ravel()

        for k, band in enumerate(BANDS):
            ax = axes_flat[k]
            mat = load_fc_matrix(patient, phase, band, method)
            if mat is None:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes)
                continue

            ch = load_channel_labels(patient)
            sp_vals, cp_vals = probe_weight_distributions(mat, ch)

            if len(sp_vals) == 0 or len(cp_vals) == 0:
                continue

            # Log-scale x-axis bins
            all_vals = np.concatenate([sp_vals, cp_vals])
            all_pos = all_vals[all_vals > 0]
            if len(all_pos) == 0:
                continue

            xmin = max(all_pos.min() * 0.5, 1e-5)
            xmax = all_pos.max() * 1.5
            bins = np.logspace(np.log10(xmin), np.log10(xmax), 40)

            ax.hist(cp_vals[cp_vals > 0], bins=bins, density=True, alpha=0.5,
                    color=CLR_CROSS, label="cross-probe")
            ax.hist(sp_vals[sp_vals > 0], bins=bins, density=True, alpha=0.5,
                    color=CLR_SAME, label="same-probe")

            # Mean lines
            ax.axvline(cp_vals.mean(), color=CLR_CROSS, ls="--", lw=1.5)
            ax.axvline(sp_vals.mean(), color=CLR_SAME, ls="--", lw=1.5)

            ax.set_xscale("log")
            ax.set_title(BAND_TEX[k], fontsize=12, fontweight="bold")
            if k == 0:
                ax.legend(fontsize=8)
            if k >= 3:
                ax.set_xlabel("|weight|", fontsize=10)

            ratio = sp_vals.mean() / (cp_vals.mean() + 1e-30)
            ax.text(0.98, 0.95, f"ratio={ratio:.1f}x",
                    transform=ax.transAxes, ha="right", va="top", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

        label = METHOD_LABELS.get(method, method.upper())
        fig.tight_layout()
        tag = "a" if method == "msc" else "b"
        fname_stem = "msc" if method == "msc" else "imcoh"
        save_fig(fig, output_dir / f"fig_C1{tag}_weight_dist_{fname_stem}_{patient}_{phase}")


# ---------------------------------------------------------------------------
# Figure C2 — Raw spectral example
# ---------------------------------------------------------------------------

def fig_c2_spectral_example(
    patient: str = "Pat_02",
    phase: str = "rest_pre",
    output_dir: Path = SECTION2_ROOT / "fig_C",
    verbose: bool = False,
):
    """Show coherence(f) for one same-probe and one cross-probe pair."""
    from lrg_eegfc.utils.io.patient import load_mat_pat_data
    from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch

    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    N = len(ch)

    # Load timeseries
    mat_data = load_mat_pat_data(patient, phase, SEEG_DATAPATH)
    X = None
    for k in mat_data:
        if not k.startswith('_'):
            v = mat_data[k]
            if hasattr(v, 'shape') and len(v.shape) == 2 and min(v.shape) > 10:
                X = np.array(v, dtype=np.float64)
                break
    if X is None:
        print("  SKIP C2: could not load timeseries")
        return
    if X.shape[0] > X.shape[1]:
        X = X.T

    fs = 2048.0 if patient != "Pat_03" else 1024.0
    nperseg = nperseg_for_fs(fs)

    # Compute MSC and ImCoh spectral tensors
    freqs_msc, Coh_msc = compute_msc_welch(X, fs, nperseg=nperseg, metric="msc")
    freqs_im, Coh_im = compute_msc_welch(X, fs, nperseg=nperseg, metric="imcoh")

    # Find best same-probe pair (highest MSC in beta range)
    beta_mask = (freqs_msc >= 13) & (freqs_msc <= 30)
    probe_mask = build_probe_mask(ch)
    np.fill_diagonal(probe_mask, False)

    # Average beta MSC per pair
    beta_msc = Coh_msc[:, :, beta_mask].mean(axis=2)
    np.fill_diagonal(beta_msc, 0)

    # Best same-probe pair
    beta_sp = beta_msc.copy()
    beta_sp[~probe_mask] = 0
    i_sp, j_sp = np.unravel_index(beta_sp.argmax(), beta_sp.shape)

    # Best cross-probe pair with comparable MSC
    beta_cp = beta_msc.copy()
    beta_cp[probe_mask] = 0
    target_msc = beta_sp[i_sp, j_sp]
    # Find cross-probe pair closest to same MSC value
    diff = np.abs(beta_cp - target_msc)
    diff[probe_mask] = np.inf
    diff[diff == 0] = np.inf
    np.fill_diagonal(diff, np.inf)
    i_cp, j_cp = np.unravel_index(diff.argmin(), diff.shape)

    if verbose:
        print(f"    Same-probe pair: {ch[i_sp]}-{ch[j_sp]} (probe {probes[i_sp]})")
        print(f"    Cross-probe pair: {ch[i_cp]}-{ch[j_cp]} "
              f"(probes {probes[i_cp]}, {probes[j_cp]})")

    # Plot: 2x2 grid
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)

    pairs = [(i_sp, j_sp, "Same-probe"), (i_cp, j_cp, "Cross-probe")]
    spectra = [
        ("MSC", freqs_msc, Coh_msc),
        ("ImCoh", freqs_im, Coh_im),
    ]

    # Extract only the spectra we need (keep copies), then free the tensors
    sp_data = {
        ("Same-probe", "MSC"):   (freqs_msc, Coh_msc[i_sp, j_sp, :].copy()),
        ("Same-probe", "ImCoh"): (freqs_im,  Coh_im[i_sp, j_sp, :].copy()),
        ("Cross-probe", "MSC"):   (freqs_msc, Coh_msc[i_cp, j_cp, :].copy()),
        ("Cross-probe", "ImCoh"): (freqs_im,  Coh_im[i_cp, j_cp, :].copy()),
    }
    # Free the full (N, N, F) tensors immediately — they are huge
    del Coh_msc, Coh_im, beta_msc, beta_sp, beta_cp, diff, X, mat_data
    import gc; gc.collect()

    for row, (i, j, pair_label) in enumerate(pairs):
        for col, (method_label, _fq, _coh) in enumerate(spectra):
            ax = axes[row, col]
            freqs, spectrum = sp_data[(pair_label, method_label)]

            ax.plot(freqs, spectrum, color=CLR_MSC if col == 0 else CLR_IMCOH,
                    lw=1.2)

            # Shade frequency bands
            for band_name, (flo, fhi) in BRAIN_BANDS.items():
                ax.axvspan(flo, fhi, alpha=0.05, color="gray")

            ax.set_ylabel(method_label if col == 0 else "", fontsize=11)
            ax.set_xlim(0.5, 100)
            ax.set_xscale("log")

            title_parts = [f"{pair_label}: {ch[i]}–{ch[j]}"]
            if row == 0:
                title_parts.insert(0, method_label)
            ax.set_title(" — ".join(title_parts), fontsize=11)

            if row == 1:
                ax.set_xlabel("Frequency (Hz)", fontsize=10)

            ax.grid(True, alpha=0.2)

    fig.tight_layout()
    save_fig(fig, output_dir / f"fig_C2_spectral_example_{patient}_{phase}")


# ---------------------------------------------------------------------------
# Figure C3 — Weight distribution overlay across patients
# ---------------------------------------------------------------------------

def fig_c3_weight_overlay(
    patients: list[str] | None = None,
    band: str = "beta",
    phase: str = "rest_pre",
    output_dir: Path = SECTION2_ROOT / "fig_C",
    verbose: bool = False,
):
    """2 panels (MSC, ImCoh): KDE curves for all patients, same vs cross-probe."""
    if patients is None:
        patients = ALL_PATIENTS

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Titles use |ImCoh| because fc_method="imcoh_abs" loads the
    # band-averaged absolute value (Ewald 2012 / Bastos-Schoffelen 2016
    # convention). Use r"$|\mathrm{ImCoh}|$" everywhere we display that
    # quantity; r"$\mathrm{ImCoh}$" (signed, [-1,1]) is reserved for
    # fc_method="imcoh" plots.
    for ax, method, title in zip(
        axes, METHODS, ["MSC", r"$|\mathrm{ImCoh}|$"]
    ):
        for pat in patients:
            mat = load_fc_matrix(pat, phase, band, method)
            if mat is None:
                continue

            ch = load_channel_labels(pat)
            sp_vals, cp_vals = probe_weight_distributions(mat, ch)

            color = PATIENT_COLORS[pat]

            # KDE for cross-probe (solid) and same-probe (dashed)
            for vals, ls, label_suffix in [
                (cp_vals, "-", "cross"),
                (sp_vals, "--", "same"),
            ]:
                vals_pos = vals[vals > 0]
                if len(vals_pos) < 20:
                    continue
                log_vals = np.log10(vals_pos)
                kde = gaussian_kde(log_vals, bw_method=0.2)
                x_log = np.linspace(log_vals.min() - 0.5, log_vals.max() + 0.5, 200)
                x_lin = 10**x_log
                label = f"{pat}" if label_suffix == "cross" else None
                ax.plot(x_lin, kde(x_log), ls=ls, color=color, lw=1.5,
                        alpha=0.8, label=label)

        ax.set_xscale("log")
        ax.set_xlabel("|weight|", fontsize=11)
        ax.set_ylabel("KDE density (log scale)", fontsize=11)
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(True, alpha=0.2)

    # Add dashed/solid legend
    from matplotlib.lines import Line2D
    fig.legend(
        handles=[
            Line2D([0], [0], color="gray", ls="-", lw=1.5, label="cross-probe"),
            Line2D([0], [0], color="gray", ls="--", lw=1.5, label="same-probe"),
        ],
        loc="lower center", fontsize=10, ncol=2, bbox_to_anchor=(0.5, -0.02),
    )

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.tight_layout()
    save_fig(fig, output_dir / f"fig_C3_weight_overlay_all_patients_{band}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Section 2 weight distribution figures.")
    parser.add_argument("--output-dir", type=Path, default=SECTION2_ROOT)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()
    set_memory_limit(12.0)  # cap at 12 GB to prevent runaway
    out = args.output_dir / "fig_C"

    # C1
    for patient in REPR_PATIENTS:
        for phase in REPR_PHASES:
            print(f"--- Figure C1: {patient} {phase} ---")
            fig_c1_weight_distributions(patient, phase, out, args.verbose)
            checkpoint_cleanup(f"C1 {patient} {phase}")

    # C2 — spectral tensor per (patient, phase) is the heaviest step
    for patient in REPR_PATIENTS:
        for phase in REPR_PHASES:
            print(f"\n--- Figure C2: Spectral example {patient} {phase} ---")
            fig_c2_spectral_example(patient, phase, out, args.verbose)
            checkpoint_cleanup(f"C2 {patient} {phase}")

    # C3
    for band in REPR_BANDS:
        for phase in REPR_PHASES:
            print(f"\n--- Figure C3: Weight overlay ({band}, {phase}) ---")
            fig_c3_weight_overlay(band=band, phase=phase,
                                  output_dir=out, verbose=args.verbose)
            checkpoint_cleanup(f"C3 {band} {phase}")

    print("\nDone.")


if __name__ == "__main__":
    main()
