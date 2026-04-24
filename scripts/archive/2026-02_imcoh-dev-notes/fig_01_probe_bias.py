#!/usr/bin/env python3
"""Section 2.2 figure: same-probe bias comparison (MSC vs ImCoh).

Produces three panels motivating the estimator switch from MSC to ImCoh:

  A) Enrichment ratio heatmap (patients x bands, side-by-side MSC/ImCoh)
  B) Weight distribution comparison for a representative band
  C) Community-probe enrichment vs LRG scale

Uses the shared probe utilities from :mod:`lrg_eegfc.utils.probe`.

Run:
  python scripts/10_notes_imcoh/fig_01_probe_bias.py [-v]
  python scripts/10_notes_imcoh/fig_01_probe_bias.py --band beta --output-dir data/outputs/figures/section2/fig_01
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import (
    setup_script_env,
    save_figure,
    write_report,
    load_all_fc,
    load_all_lrg,
)

ROOT = setup_script_env()

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.probe import (
    extract_probe_labels,
    compute_probe_weight_ratio,
    probe_weight_distributions,
    compute_enrichment_vs_scale,
)

# ---------------------------------------------------------------------------
# Defaults (all from config, never hardcoded)
# ---------------------------------------------------------------------------
ALL_PATIENTS = list(PATIENTS_4PHASE) + ["Pat_06"]
BANDS = list(BRAIN_BANDS_NAMES)
BAND_TEX = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]
DEFAULT_OUTPUT_DIR = FIGURES_ROOT / "section2" / "fig_01"
DEFAULT_PHASE = "rest_pre"
DEFAULT_REPR_BAND = "beta"
N_COMMUNITIES = [3, 5, 10, 15, 20, 30]

# Plotting colours
CLR_SAME = "#E64A19"   # orange-red for same-probe
CLR_CROSS = "#1565C0"  # blue for cross-probe
CLR_MSC = "#2166AC"
CLR_IMCOH = "#B2182B"


# ---------------------------------------------------------------------------
# Channel label loading
# ---------------------------------------------------------------------------

def load_channel_labels(patient: str) -> list[str]:
    """Load cleaned channel labels for *patient*.

    Handles varying CSV formats: with/without header, quoted labels,
    reference suffixes (e.g. ``,G2``).
    """
    for ext in ("csv", "txt"):
        fpath = SEEG_DATAPATH / patient / f"channel_labels.{ext}"
        if not fpath.exists():
            continue
        labels = []
        with open(fpath) as f:
            for i, line in enumerate(f):
                line = line.strip().strip('"')
                if not line:
                    continue
                # Skip header row
                if i == 0 and line.lower() == "label":
                    continue
                # Take first field before comma (removes ,G2 reference)
                label = line.split(",")[0].strip().strip('"').replace(" ", "")
                if label:
                    labels.append(label)
        if labels:
            return labels
    raise FileNotFoundError(f"No channel_labels for {patient}")


# ---------------------------------------------------------------------------
# Figure A — Enrichment ratio heatmap
# ---------------------------------------------------------------------------

def fig_a_enrichment_heatmap(
    patients: list[str],
    bands: list[str],
    phase: str,
    output_dir: Path,
    verbose: bool = False,
) -> pd.DataFrame:
    """2-panel heatmap: MSC vs ImCoh same-probe/cross-probe mean ratio."""
    methods = ["msc", "imcoh_abs"]
    # Collect ratios: (patient, band, method) -> ratio
    rows = []

    for method in methods:
        fc_data = load_all_fc(
            patients=patients, bands=bands, phases=[phase],
            fc_method=method, verbose=verbose,
        )
        for (pat, ph, band), W in fc_data.items():
            ch = load_channel_labels(pat)
            ratio = compute_probe_weight_ratio(W, ch)
            rows.append({
                "patient": pat, "band": band, "method": method, "ratio": ratio,
            })

    df = pd.DataFrame(rows)
    if df.empty:
        print("  WARNING: no data for enrichment heatmap")
        return df

    fig, axes = plt.subplots(1, 2, figsize=(14, 0.6 * len(patients) + 2),
                             sharey=True)

    for ax, method, title in zip(
        axes, methods, ["MSC enrichment ratio", r"$|\mathrm{ImCoh}|$ enrichment ratio"]
    ):
        sub = df[df["method"] == method]
        piv = sub.pivot(index="patient", columns="band", values="ratio")
        piv = piv.reindex(index=patients, columns=bands)

        vmax = max(piv.max().max(), 2.0)
        norm = TwoSlopeNorm(vmin=0, vcenter=1.0, vmax=vmax)
        im = ax.imshow(
            piv.values, aspect="auto", cmap="RdBu_r", norm=norm,
            interpolation="nearest",
        )

        # Annotate cells
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = piv.values[i, j]
                if np.isfinite(v):
                    color = "white" if v > vmax * 0.7 else "black"
                    ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                            fontsize=8, color=color)

        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels(BAND_TEX, fontsize=10)
        ax.set_yticks(range(len(patients)))
        ax.set_yticklabels(patients, fontsize=10)
        ax.set_title(title, fontsize=12)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="ratio")

    fig.suptitle(
        f"Same-probe / cross-probe mean weight ratio  ({phase})",
        fontsize=14, y=1.02,
    )
    fig.tight_layout()
    save_figure(fig, output_dir / "fig_01A_probe_enrichment_heatmap",
                what="Heatmap of same-probe/cross-probe mean weight ratio for MSC and ImCoh.",
                proves="MSC has 2-8x inflation on same-probe pairs; ImCoh is near 1x.")
    return df


# ---------------------------------------------------------------------------
# Figure B — Weight distribution comparison
# ---------------------------------------------------------------------------

def fig_b_weight_distributions(
    patients: list[str],
    band: str,
    phase: str,
    output_dir: Path,
    verbose: bool = False,
):
    """2-row x N-col KDE: same-probe vs cross-probe weight distributions."""
    methods = ["msc", "imcoh_abs"]
    n_pats = len(patients)

    fig, axes = plt.subplots(2, n_pats, figsize=(3.5 * n_pats, 6), squeeze=False)

    for row, method in enumerate(methods):
        fc_data = load_all_fc(
            patients=patients, bands=[band], phases=[phase],
            fc_method=method, verbose=verbose,
        )
        for col, pat in enumerate(patients):
            ax = axes[row, col]
            key = (pat, phase, band)
            if key not in fc_data:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=12, color="gray")
                continue

            W = fc_data[key]
            ch = load_channel_labels(pat)
            sp_vals, cp_vals = probe_weight_distributions(W, ch)

            # Plot histograms
            bins = np.linspace(0, max(sp_vals.max(), cp_vals.max()) * 1.05, 50)
            ax.hist(cp_vals, bins=bins, density=True, alpha=0.6,
                    color=CLR_CROSS, label="cross-probe")
            ax.hist(sp_vals, bins=bins, density=True, alpha=0.6,
                    color=CLR_SAME, label="same-probe")

            # Mean lines
            ax.axvline(cp_vals.mean(), color=CLR_CROSS, ls="--", lw=1.5)
            ax.axvline(sp_vals.mean(), color=CLR_SAME, ls="--", lw=1.5)

            ratio = sp_vals.mean() / (cp_vals.mean() + 1e-30)
            ax.set_title(f"{pat}\nratio={ratio:.1f}x", fontsize=10)
            if row == 0 and col == 0:
                ax.legend(fontsize=7, loc="upper right")
            if col == 0:
                ax.set_ylabel(f"{method.upper()}\ndensity", fontsize=10)
            ax.set_xlabel("|weight|" if row == 1 else "", fontsize=9)
            ax.tick_params(labelsize=8)

    band_tex = BRAIN_BAND_TEX_DICT[band]
    fig.suptitle(
        f"Weight distributions: same-probe vs cross-probe  "
        f"({band_tex}, {phase})",
        fontsize=14, y=1.02,
    )
    fig.tight_layout()
    save_figure(fig, output_dir / f"fig_01B_weight_distributions_{band}",
                what=f"Same-probe (orange) vs cross-probe (blue) weight distributions for {band} band.",
                proves="MSC same-probe distribution is shifted right; ImCoh distributions overlap.")


# ---------------------------------------------------------------------------
# Figure C — Enrichment ratio vs community scale
# ---------------------------------------------------------------------------

def fig_c_enrichment_vs_scale(
    patients: list[str],
    band: str,
    phase: str,
    output_dir: Path,
    n_communities: list[int] | None = None,
    verbose: bool = False,
):
    """Line plot: community-probe enrichment at different LRG scales."""
    if n_communities is None:
        n_communities = N_COMMUNITIES
    methods = ["msc", "imcoh_abs"]

    fig, ax = plt.subplots(figsize=(9, 6))

    # Patient colors
    cmap = plt.get_cmap("tab10")

    for pi, pat in enumerate(patients):
        ch = load_channel_labels(pat)
        probes = extract_probe_labels(ch)
        color = cmap(pi % 10)

        for method, ls, marker in [("msc", "--", "o"), ("imcoh_abs", "-", "s")]:
            lrg_data = load_all_lrg(
                patients=[pat], bands=[band], phases=[phase],
                fc_method=method, verbose=False,
            )
            key = (pat, phase, band)
            if key not in lrg_data:
                if verbose:
                    print(f"  SKIP {pat} {method}: no LRG")
                continue

            result = lrg_data[key]
            enrich = compute_enrichment_vs_scale(
                result.linkage_matrix, probes,
                n_communities=n_communities,
                n_nodes=result.n_nodes,
            )
            if not enrich:
                continue

            xs = sorted(enrich.keys())
            ys = [enrich[x] for x in xs]

            label = f"{pat} {method.upper()}" if method == "msc" else fr"{pat} $|\mathrm{{ImCoh}}|$"
            ax.plot(xs, ys, ls=ls, marker=marker, color=color,
                    markersize=5, lw=1.5, alpha=0.85, label=label)

    ax.axhline(1.0, color="gray", ls=":", lw=1, label="no bias")
    ax.axhline(1.5, color="red", ls=":", lw=0.8, alpha=0.5, label="1.5x threshold")
    ax.set_xscale("log")
    ax.set_xlabel("Number of communities (n)", fontsize=12)
    ax.set_ylabel("Same-community same-probe enrichment", fontsize=12)
    ax.set_xticks(n_communities)
    ax.set_xticklabels([str(n) for n in n_communities])

    band_tex = BRAIN_BAND_TEX_DICT[band]
    ax.set_title(
        f"Community-probe enrichment vs scale  ({band_tex}, {phase})\n"
        fr"Dashed = MSC, Solid = $|\mathrm{{ImCoh}}|$",
        fontsize=13,
    )
    ax.legend(fontsize=7, loc="upper right", ncol=2)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    save_figure(fig, output_dir / f"fig_01C_enrichment_vs_scale_{band}",
                what="Community-probe enrichment ratio across LRG scales for MSC vs ImCoh.",
                proves="MSC enrichment grows at coarse scales; ImCoh stays closer to 1.")


# ---------------------------------------------------------------------------
# Stats report
# ---------------------------------------------------------------------------

def write_stats_report(
    enrichment_df: pd.DataFrame,
    patients: list[str],
    band: str,
    phase: str,
    output_dir: Path,
):
    """Write a summary stats report."""
    lines = [
        f"Phase: {phase}",
        f"Patients: {', '.join(patients)}",
        "",
        "## Enrichment ratio table (same-probe mean / cross-probe mean)",
        "",
    ]

    if not enrichment_df.empty:
        # Pivot for MSC
        _DISPLAY = {"msc": "MSC", "imcoh_abs": "|ImCoh|"}
        for method in ["msc", "imcoh_abs"]:
            sub = enrichment_df[enrichment_df["method"] == method]
            piv = sub.pivot(index="patient", columns="band", values="ratio")
            piv = piv.reindex(index=patients, columns=BANDS)
            lines.append(f"### {_DISPLAY[method]}")
            lines.append("")
            lines.append(piv.to_markdown(floatfmt=".2f"))
            lines.append("")

        # Per-patient: band with highest MSC enrichment
        lines.append("## Highest MSC enrichment per patient")
        lines.append("")
        msc_df = enrichment_df[enrichment_df["method"] == "msc"]
        for pat in patients:
            psub = msc_df[msc_df["patient"] == pat]
            if not psub.empty:
                worst = psub.loc[psub["ratio"].idxmax()]
                lines.append(f"- {pat}: {worst['band']} ({worst['ratio']:.2f}x)")
        lines.append("")

        # ImCoh mean
        imcoh_df = enrichment_df[enrichment_df["method"] == "imcoh_abs"]
        if not imcoh_df.empty:
            lines.append("## ImCoh enrichment summary")
            lines.append(f"- Mean: {imcoh_df['ratio'].mean():.2f}x")
            lines.append(f"- Max: {imcoh_df['ratio'].max():.2f}x "
                         f"({imcoh_df.loc[imcoh_df['ratio'].idxmax(), 'patient']} "
                         f"{imcoh_df.loc[imcoh_df['ratio'].idxmax(), 'band']})")
            lines.append("")

    write_report(lines, output_dir / "fig_01_stats.md",
                 title="Figure 01 — Probe bias stats")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate Section 2.2 probe bias comparison figures."
    )
    parser.add_argument("--patients", nargs="+", default=ALL_PATIENTS,
                        help="Patient IDs (default: all 6)")
    parser.add_argument("--phase", default=DEFAULT_PHASE,
                        help=f"Recording phase (default: {DEFAULT_PHASE})")
    parser.add_argument("--band", default=DEFAULT_REPR_BAND,
                        help=f"Representative band for panels B,C (default: {DEFAULT_REPR_BAND})")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Figure A: enrichment heatmap (all patients x bands)
    print("\n--- Figure A: Enrichment ratio heatmap ---")
    df = fig_a_enrichment_heatmap(
        args.patients, BANDS, args.phase, output_dir, args.verbose,
    )

    # Figure B: weight distributions (representative band)
    print(f"\n--- Figure B: Weight distributions ({args.band}) ---")
    fig_b_weight_distributions(
        args.patients, args.band, args.phase, output_dir, args.verbose,
    )

    # Figure C: enrichment vs LRG scale (representative band)
    print(f"\n--- Figure C: Enrichment vs scale ({args.band}) ---")
    fig_c_enrichment_vs_scale(
        args.patients, args.band, args.phase, output_dir, verbose=args.verbose,
    )

    # Stats report
    print("\n--- Stats report ---")
    write_stats_report(df, args.patients, args.band, args.phase, output_dir)

    print(f"\nDone. Figures saved to {output_dir}")


if __name__ == "__main__":
    main()
