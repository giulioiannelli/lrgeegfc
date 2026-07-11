#!/usr/bin/env python3
r"""Slide 03 (sEEG dataset) — a few raw sEEG traces as a wide horizontal banner.

The "here is the raw signal" panel: a HANDFUL of intracranial channels shown as stacked
voltage traces over a long window, so the slide conveys what a contact records and at what
temporal resolution (2048 Hz — the wiggles are real). A right-pointing arrow at the end of
each trace shows the recording continues beyond the shown window.

Channels: a few random non-SOZ contacts (dark) + a couple random SEIZURE-ONSET-ZONE contacts
(RED — the clinical annotation, tying to the double-probe framing). Common gain keeps relative
amplitudes honest; a time + amplitude scale bar replaces numeric axes. Wide aspect for a slide.

Amplitude is in the recording's native export units (labelled a.u.; pass --uv if the vendor
export is confirmed microvolts). Purely illustrative — no analysis.

Reads : raw timeseries (utils.io.patient.load_timeseries), channel labels, SOZ labels.
Writes: data/outputs/figures/talk/seeg_traces_stacked.pdf   (+ --qa white-matte PNG)
Usage : python scripts/07_figures/gen_seeg_traces.py [--patient Pat_05] [--phase rest_pre]
                 [--start 60] [--dur 8] [--n-ctx 3] [--n-soz 2] [--seed 7] [--uv] [--qa out.png]
"""
from __future__ import annotations

import argparse

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt

from lrg_eegfc.config.paths import SEEG_DATAPATH, FIGURES_ROOT
from lrg_eegfc.config.const import FS_OVERRIDES, DEFAULT_SAMPLE_RATE
from lrg_eegfc.utils.io.patient import (
    load_timeseries, load_channel_labels, load_epileptic_nodes,
)
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()

DARK, RED, GRAY, INK = "#1f2a37", "#c0392b", "#6f757c", "#1a1a1a"


def _nice(v, opts=(10, 20, 50, 100, 200, 500, 1000, 2000, 5000)):
    return min(opts, key=lambda o: abs(o - v))


def select_channels(labels, soz, n_ctx, n_soz, rng):
    """A few RANDOM non-SOZ channels + a few random SOZ channels, in montage order."""
    non = [i for i, l in enumerate(labels) if l not in soz]
    soz_i = [i for i, l in enumerate(labels) if l in soz]
    ctx = rng.choice(non, size=min(n_ctx, len(non)), replace=False).tolist() if non else []
    shw = rng.choice(soz_i, size=min(n_soz, len(soz_i)), replace=False).tolist() if soz_i else []
    return sorted(set(ctx) | set(shw))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default="Pat_05")
    ap.add_argument("--phase", default="rest_pre")
    ap.add_argument("--start", type=float, default=60.0, help="window start (s)")
    ap.add_argument("--dur", type=float, default=8.0, help="window length (s)")
    ap.add_argument("--n-ctx", type=int, default=3, help="random non-SOZ channels (dark)")
    ap.add_argument("--n-soz", type=int, default=2, help="random SOZ channels (red)")
    ap.add_argument("--seed", type=int, default=7, help="RNG seed for the channel draw")
    ap.add_argument("--uv", action="store_true", help="label amplitude in microvolts")
    ap.add_argument("--qa", help="also write a white-matte QA PNG here")
    args = ap.parse_args()

    fs = FS_OVERRIDES.get(args.patient, DEFAULT_SAMPLE_RATE)
    ts = load_timeseries(args.patient, args.phase, SEEG_DATAPATH)
    if ts.shape[0] < ts.shape[1]:                      # -> (n_samples, n_channels)
        ts = ts.T
    n_samp = ts.shape[0]
    labels = list(load_channel_labels(args.patient, SEEG_DATAPATH))
    try:
        soz = set(load_epileptic_nodes(args.patient, SEEG_DATAPATH))
    except Exception:
        soz = set()

    n = int(args.dur * fs)
    s0 = int(min(args.start, max(0.0, n_samp / fs - args.dur)) * fs)   # clamp to recording
    seg = ts[s0:s0 + n].astype(float)
    seg = seg - seg.mean(0, keepdims=True)
    t = np.arange(n) / fs
    D = args.dur

    rng = np.random.default_rng(args.seed)
    sel = select_channels(labels, soz, args.n_ctx, args.n_soz, rng)
    M = len(sel)
    spacing = 6.0 * float(np.median(seg[:, sel].std(0)))
    amp_bar = _nice(spacing * 0.5)
    unit = "µV" if args.uv else "a.u."

    fig, ax = plt.subplots(figsize=(13.5, 0.5 * M + 1.1))     # wide horizontal banner
    for k, ci in enumerate(sel):
        y0 = (M - 1 - k) * spacing
        is_soz = labels[ci] in soz
        col = RED if is_soz else DARK
        ax.plot(t, seg[:, ci] + y0, lw=(0.8 if is_soz else 0.6), color=col,
                solid_capstyle="round")
        ax.text(-0.014 * D, y0, labels[ci], ha="right", va="center",
                fontsize=11, color=col, family="monospace")
        # continuation arrow at the end of each trace -> the recording goes on
        ax.annotate("", xy=(D * 1.055, y0), xytext=(D * 1.004, y0),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.4,
                                    shrinkA=0, shrinkB=0))
    ax.text(D * 1.028, (M - 0.5) * spacing, "continues", ha="left", va="center",
            fontsize=9.5, color=GRAY, style="italic")

    # scale bars (amplitude + time), bottom-right of the data area
    x0, yb = D * 0.86, -1.5 * spacing
    ax.plot([x0, x0 + 1.0], [yb, yb], color=INK, lw=2.6)
    ax.text(x0 + 0.5, yb - 0.42 * spacing, "1 s", ha="center", va="top", fontsize=11.5, color=INK)
    ax.plot([x0, x0], [yb, yb + amp_bar], color=INK, lw=2.6)
    ax.text(x0 - 0.008 * D, yb + amp_bar / 2, f"{amp_bar} {unit}", ha="right",
            va="center", fontsize=11.5, color=INK)

    # SOZ key (only if any shown) + resolution note, above the traces
    if soz & {labels[i] for i in sel}:
        ax.text(-0.014 * D, (M - 0.25) * spacing, "red = seizure-onset zone (SOZ)",
                ha="left", va="bottom", fontsize=11, color=RED)
    ax.text(D, (M - 0.25) * spacing, f"intracranial sEEG · {int(fs)} Hz",
            ha="right", va="bottom", fontsize=11, color=GRAY, style="italic")

    ax.set_xlim(-0.075 * D, D * 1.11)
    ax.set_ylim(yb - 0.8 * spacing, (M + 0.05) * spacing)
    ax.axis("off")

    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "seeg_traces_stacked.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    shown = [labels[i] for i in sel]
    print(f"{args.patient} {args.phase}  window {s0/fs:.1f}-{s0/fs+D:.1f}s  "
          f"channels={M} {shown} (SOZ={[labels[i] for i in sel if labels[i] in soz]})  "
          f"spacing={spacing:.0f} amp_bar={amp_bar}{unit}")
    print("wrote", out)
    if args.qa:
        fig.savefig(args.qa, facecolor="white", bbox_inches="tight", dpi=130)
        print("wrote QA", args.qa)
    plt.close(fig)


if __name__ == "__main__":
    main()
