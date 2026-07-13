#!/usr/bin/env python3
r"""Cool 4-panel gallery of consolidation DYNAMICS archetypes, on the recovered
mst@0.20 mesoscale scheme -- the "3D blobs + wall density-contours + glowing flow"
style of fig_arc_attractor_3d, but showing FOUR qualitatively different behaviours,
one clean example of each:

    anchor      -- the state barely moves across phases (task never displaces it)
    reset       -- task displaces the state, rest_post returns to baseline
    reorganize  -- task displaces it, rest_post lands in a NEW place (not baseline,
                   not the task attractor)
    trace       -- task displaces it, rest_post HOLDS the task displacement

Reads the cached embeddings written by 16_windowed_trajectory_mst020.py and reuses
the rendering helpers from fig_arc_attractor_3d (no forks). Picks are auto-selected
from dynamics_discriminants.csv (cleanest score per archetype) unless --cells given.

IMPORTANT (honesty): these are archetypes of the RAW windowed trajectory. The raw
"hold" is uncorrelated with the matched-strength trace (Spearman ~0), so a panel
labelled "trace" here means the state-space HELD, not that the cell cleared the gate.
Each panel is annotated with the cell's real T_test for transparency.

Writes: data/preprint/figures/new_results_sec2/dynamics_gallery_s{scale}.pdf
"""
from __future__ import annotations

import argparse
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
import fig_arc_attractor_3d as att  # noqa: E402  (reuse: blobs, wall contours, smoothing)

from lrg_eegfc.visuals.styles import use_lrg_style  # noqa: E402

use_lrg_style()

PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
PHLAB = {"rest_pre": r"rest$_\mathrm{pre}$", "task_learn": r"task$_\mathrm{learn}$",
         "task_test": r"task$_\mathrm{test}$", "rest_post": r"rest$_\mathrm{post}$"}
TCMAP = att.TCMAP

EMB = ROOT / "data" / "sparsified_arc" / "windowed_traj_mst020" / "emb"
DISC = ROOT / "data" / "sparsified_arc" / "windowed_traj_mst020" / "dynamics_discriminants.csv"
PER = ROOT / "data" / "sparsified_arc" / "enc_inf_arc_mst020" / "per_cell.csv"
OUT = ROOT / "data" / "preprint" / "figures" / "new_results_sec2"
OUT.mkdir(parents=True, exist_ok=True)


def load_embedding(pat, band, scale):
    d = np.load(EMB / f"{pat}_{band}_s{scale:g}.npz", allow_pickle=True)
    P, ph = [], []
    for p in PHASES:
        x, y, z = d[f"{p}_x"], d[f"{p}_y"], d[f"{p}_z"]
        P.append(np.column_stack([x, y, z]))
        ph += [p] * len(x)
    meta = json.loads(str(d["meta"]))
    return np.vstack(P), np.asarray(ph), meta


def phase_dists(P, ph):
    """Rest_post distance to HOME (rest_pre) and to the TASK attractor (task_test),
    between phase centroids, normalized by pooled within-phase scatter -- the direct
    reset/reorganize/trace discriminant, measured in the plotted 3-D space."""
    C = {p: P[ph == p].mean(0) for p in PHASES}
    sig = np.sqrt(np.mean([np.mean(np.sum((P[ph == p] - C[p]) ** 2, axis=1))
                           for p in PHASES])) + 1e-9
    return (float(np.linalg.norm(C["rest_post"] - C["rest_pre"]) / sig),
            float(np.linalg.norm(C["rest_post"] - C["task_test"]) / sig))


def autopick(scale):
    """Cleanest cell per archetype from the discriminants table."""
    d = pd.read_csv(DISC)
    d = d[d.scale == scale].copy()
    if PER.exists():
        d = d.merge(pd.read_csv(PER)[["patient", "band", "T_test__obs_s1", "T_test__obs_max"]],
                    on=["patient", "band"], how="left")
    dt, dp, al = d["D_task"], d["D_post"], d["align"]
    med = dt.median()
    d["anchor"] = -(dt + dp)
    d["reset"] = np.where(dt > med, dt - 3 * dp, -9)
    d["trace"] = np.where(dt > med, dp * al.clip(0), -9)
    d["reorganize"] = np.where((dt > med) & (dp > dp.median()), dp * (1 - al.clip(-1, 1)), -9)
    picks = {}
    used = set()
    for arch in ["anchor", "reset", "trace", "reorganize"]:
        cand = d.sort_values(arch, ascending=False)
        for _, r in cand.iterrows():
            key = (r.patient, r.band)
            if key not in used:
                picks[arch] = r
                used.add(key)
                break
    return picks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", type=float, default=5.6)
    ap.add_argument("--cells", nargs="+", default=None,
                    help="override picks, e.g. anchor=Pat_15:theta reset=Pat_10:beta ...")
    args = ap.parse_args()

    order = ["anchor", "reset", "reorganize", "trace"]
    label = {"anchor": "ANCHOR  (barely moves)", "reset": "RESET  (moves, returns home)",
             "reorganize": "REORGANIZE  (moves, new place)", "trace": "TRACE  (moves, holds)"}

    if args.cells:
        pc = pd.read_csv(PER) if PER.exists() else None
        chosen = {}
        for tok in args.cells:
            arch, cell = tok.split("=")
            pat, band = cell.split(":")
            _, _, meta = load_embedding(pat, band, args.scale)
            if pc is not None:
                row = pc[(pc.patient == pat) & (pc.band == band)]
                meta["T_test__obs_s1"] = float(row["T_test__obs_s1"].iloc[0]) if len(row) else np.nan
            chosen[arch] = (pat, band, meta)
    else:
        picks = autopick(args.scale)
        chosen = {a: (picks[a].patient, picks[a].band, picks[a].to_dict()) for a in order}

    fig = plt.figure(figsize=(10.4, 9.4))
    for i, arch in enumerate(order, 1):
        pat, band, meta = chosen[arch]
        P, ph, _ = load_embedding(pat, band, args.scale)
        d_home, d_task = phase_dists(P, ph)
        ax = fig.add_subplot(2, 2, i, projection="3d")
        tt = meta.get("T_test__obs_s1", float("nan"))
        sub = (rf"{pat} $\cdot$ {band}    $d_\mathrm{{home}}$={d_home:.2f}  "
               rf"$d_\mathrm{{task}}$={d_task:.2f}    real $T$={tt:+.2f}")
        att.render_portrait(ax, P, ph, np.arange(len(P)), label[arch], subtitle=sub)

    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=TCMAP(c),
                      markeredgecolor="0.4", markersize=9, label=PHLAB[p])
               for p, c in zip(PHASES, [0.0, 0.34, 0.66, 1.0])]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               ncol=4, frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    out = OUT / f"dynamics_gallery_s{args.scale:g}.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    print(f"wrote {out}")
    for arch in order:
        pat, band, meta = chosen[arch]
        print(f"  {label[arch]:11s} = {pat}/{band}  "
              f"D_task={meta['D_task']:.2f} D_post={meta['D_post']:.2f} align={meta['align']:+.2f}")


if __name__ == "__main__":
    main()
