#!/usr/bin/env python3
r"""fig (talk slide 13): raw FC reads the same in every phase — three matrices.

One patient, one band, the three cross-phase states side by side
(rest_pre | task_test | rest_post) as raw |ImCoh| adjacency matrices on a
shared log colour scale. This is the RAW-EDGE exemplar for "how the bands
manifest": you can see the network reorganize at task and lean back at rest,
but at the edge level this blunt drift looks the same in every band — the
band structure only emerges once the same weights are read through the
diffusion hierarchy (the multiscale exemplars).

Reads : cached imcoh_abs FC via workflow.fc.load_fc_matrix
Writes: data/outputs/figures/talk/fig_raw_fc_three_phase_<pat>_<band>.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency_row
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.fc import load_fc_matrix

use_lrg_style()

PHASES = ["rest_pre", "task_test", "rest_post"]
PHASE_TEX = {
    "rest_pre": r"rest$_{\mathrm{pre}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--patient", default="Pat_08")
    ap.add_argument("--band", default="beta")
    ap.add_argument("--fc-method", default="imcoh_abs")
    args = ap.parse_args()

    mats = []
    for ph in PHASES:
        M = load_fc_matrix(patient=args.patient, phase=ph,
                           band=args.band, fc_method=args.fc_method)
        if M is None:
            raise SystemExit(f"no cached FC for {args.patient}/{args.band}/{ph}")
        mats.append(M)

    chn = pd.read_csv(SEEG_DATAPATH / args.patient / "channel_labels.csv")
    chnames = chn["label"].astype(str).tolist()

    fig, _ = plot_fc_adjacency_row(
        mats,
        titles=[PHASE_TEX[p] for p in PHASES],
        fc_method=args.fc_method,
        band=args.band,
        tick_labels="chnames",
        channel_labels=chnames,
        log_scale=True,
        shared_scale=True,
    )

    out = (Path(FIGURES_ROOT) / "talk"
           / f"fig_raw_fc_three_phase_{args.patient}_{args.band}.pdf")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
