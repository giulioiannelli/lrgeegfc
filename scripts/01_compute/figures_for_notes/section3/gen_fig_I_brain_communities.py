#!/usr/bin/env python3
"""Section 3 Figure I — Brain-space community plots MSC vs |ImCoh|.

Two 2x3 figures:
  fig_I1: fc_method=msc        rows = (n=5, n=10), cols = Pat_02, Pat_03, Pat_05
  fig_I2: fc_method=imcoh_abs  same layout
Phase = rest_pre, band = beta.

Each panel: nilearn glass-brain axial projection. Nodes at MNI coords,
coloured by community at the fixed n cut, sized by node strength. Top 10%
strongest edges shown. Same-shaft enrichment ratio annotated per panel.

Companion stats md:
  fig_I_community_enrichment_stats.md  — per-(fc_method, n, patient) ratio
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt

from _shared import apply_pub_style, save_fig
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.spatial import plot_brain_connectome_at_n


SECTION3_ROOT = FIGURES_ROOT / "section3"
FIG_I_DIR = SECTION3_ROOT / "fig_I"
WRITING_DIR = SECTION3_ROOT / "for_writing_agent"

PATIENTS = ["Pat_02", "Pat_03", "Pat_05"]
N_VALUES = [5, 20]
PHASE = "rest_pre"
BAND = "beta"

FC_METHODS = [("msc", r"$\mathrm{MSC}$", "MSC"),
              ("imcoh_abs", r"$|\mathrm{ImCoh}|$", "ImCoh")]


def render_method(fc_method: str, label_tex: str, label_short: str,
                  out_pdf: Path, verbose: bool) -> dict:
    n_rows, n_cols = len(N_VALUES), len(PATIENTS)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.0 * n_cols, 2.5 * n_rows),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.04, right=0.99, top=0.94, bottom=0.02,
                        wspace=0.0, hspace=0.05)
    enrichments: dict = {}

    for row, n in enumerate(N_VALUES):
        for col, pat in enumerate(PATIENTS):
            ax = axes[row, col]
            try:
                info = plot_brain_connectome_at_n(
                    ax, pat, PHASE, BAND, fc_method, n,
                    edge_top_pct=10.0,
                    display_mode="z",
                    palette="tab20",
                    title=None,
                    annotate_enrichment=False,
                    zoom_to_nodes=True,
                )
                enrichments[(pat, n)] = info["enrichment"]
                if verbose:
                    print(f"  [{label_short} {pat} n={n}] "
                          f"enrichment={info['enrichment']:.2f}, "
                          f"n_nodes={info['n_nodes']}")
            except Exception as exc:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"{pat}\nfailed:\n{exc}",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=8, color="red")
                enrichments[(pat, n)] = float("nan")

            if row == 0:
                ax.set_title(pat, fontsize=10, fontweight="bold")
            if col == 0:
                ax.text(-0.02, 0.5, f"n = {n}",
                        transform=ax.transAxes, rotation=90,
                        ha="right", va="center",
                        fontsize=11, fontweight="bold")

    save_fig(fig, out_pdf)
    shutil.copy(out_pdf, WRITING_DIR / out_pdf.name)
    return enrichments


def write_stats(all_enrichments: dict, out_md: Path) -> None:
    lines = []
    lines.append("# Figure I — Same-shaft enrichment by fc_method × n × patient\n")
    lines.append(
        f"Phase: `{PHASE}`, band: `{BAND}`. Enrichment = "
        "P(same-probe | same-community) / P(same-probe). "
        "1.0 = no probe bias; >>1 = communities follow electrode geometry.\n"
    )
    for fc_method, _, label_short in FC_METHODS:
        lines.append(f"## {label_short}\n")
        header = "| Patient | " + " | ".join(f"n={n}" for n in N_VALUES) + " |"
        sep = "|---" * (1 + len(N_VALUES)) + "|"
        lines.append(header)
        lines.append(sep)
        for pat in PATIENTS:
            cells = []
            for n in N_VALUES:
                v = all_enrichments[fc_method].get((pat, n), float("nan"))
                cells.append("nan" if not np.isfinite(v) else f"{v:.2f}")
            lines.append(f"| {pat} | " + " | ".join(cells) + " |")
        lines.append("")

    # Mean per (fc_method, n)
    lines.append("## Mean enrichment across patients\n")
    lines.append("| fc_method | " + " | ".join(f"n={n}" for n in N_VALUES) + " |")
    lines.append("|---" * (1 + len(N_VALUES)) + "|")
    for fc_method, _, label_short in FC_METHODS:
        cells = []
        for n in N_VALUES:
            vals = [all_enrichments[fc_method].get((pat, n), float("nan"))
                    for pat in PATIENTS]
            vals = [v for v in vals if np.isfinite(v)]
            cells.append("nan" if not vals else f"{np.mean(vals):.2f}")
        lines.append(f"| {label_short} | " + " | ".join(cells) + " |")
    lines.append("")

    lines.append(
        "**Expected:** MSC enrichment > 2× (probe-bias dominated); "
        "|ImCoh| enrichment ~ 1× (volume-conduction-immune). Pat_03 was "
        "recorded at 1024 Hz (others 2048 Hz).\n"
    )
    out_md.write_text("\n".join(lines))
    print(f"  Saved: {out_md}")


def main(verbose: bool = False) -> None:
    apply_pub_style()
    FIG_I_DIR.mkdir(parents=True, exist_ok=True)
    WRITING_DIR.mkdir(parents=True, exist_ok=True)

    n_tag = "_".join(f"n{n}" for n in N_VALUES)
    out_paths = {
        "msc": FIG_I_DIR / f"fig_I1_brain_communities_MSC_{n_tag}.pdf",
        "imcoh_abs": FIG_I_DIR / f"fig_I2_brain_communities_ImCoh_{n_tag}.pdf",
    }

    all_enrichments: dict = {}
    for fc_method, label_tex, label_short in FC_METHODS:
        all_enrichments[fc_method] = render_method(
            fc_method, label_tex, label_short,
            out_paths[fc_method], verbose,
        )

    stats_md = FIG_I_DIR / "fig_I_community_enrichment_stats.md"
    write_stats(all_enrichments, stats_md)
    shutil.copy(stats_md, WRITING_DIR / stats_md.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(verbose=args.verbose)
