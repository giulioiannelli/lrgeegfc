#!/usr/bin/env python3
r"""Windowed consolidation trajectory in the encoding/inference/residual 3-D space,
on the RECOVERED scale-dependent scheme (mst@0.20, mesoscale).

This is the scale-dependent redo of fig_arc_attractor_3d (which was DENSE, tau_min).
Slide a 30 s window across the whole recording (rest_pre -> task_learn -> task_test
-> rest_post), read the LRG cophenetic state per window ON THE mst@0.20 BACKBONE AT A
CHOSEN DIMENSIONLESS SCALE s = tau*lambda_max, and draw the time-ordered trajectory
through the 3-D state space whose axes are the two things the task did:

    encoding x = state . (task_learn - rest_pre)          (the leap the task makes)
    inference y = state . (task_test  - task_learn)_|x     (the drift within the task)
    residual  z = leading residual direction

Purpose: SHOW, per (patient, band), whether rest_post HOLDS the displaced state
(tracer) or RETURNS to baseline (resetter) IN THE NEW SCHEME -- eyeball band
selectivity directly, do not trust a dense-data retention number. Prints a
hold/return diagnostic per (patient, band); nulls (theta/low_gamma) are computed
for the diagnostic but need not be plotted.

NOT a test: overlapping windows are autocorrelated; the verdict of record is the
matched-strength gate (05_enc_inf_arc / 13_matched_strength_mst020). This is the
SHAPE of that verified effect on the recovered backbone.

Writes: data/preprint/figures/new_results_sec2/traj_{patient}_{band}_s{scale}.pdf
        data/sparsified_arc/windowed_traj_mst020/hold_return.csv
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.ndimage import gaussian_filter1d

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import segment_ffts, imcoh_abs_cube, band_abs_average
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = list(BRAIN_BANDS)
PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
PHLAB = {"rest_pre": r"rest$_\mathrm{pre}$", "task_learn": r"task$_\mathrm{learn}$",
         "task_test": r"task$_\mathrm{test}$", "rest_post": r"rest$_\mathrm{post}$"}

W_SEC = 30.0
OVERLAP = 0.25
FMAX = max(hi for _lo, hi in BRAIN_BANDS.values()) + 5.0
MIN_SEG = 4
FRAC = 0.20                                # mst@0.20 recovered backbone

FIG_OUT = ROOT / "data" / "preprint" / "figures" / "new_results_sec2"
DATA_OUT = ROOT / "data" / "sparsified_arc" / "windowed_traj_mst020"
EMB_OUT = DATA_OUT / "emb"
FIG_OUT.mkdir(parents=True, exist_ok=True)
EMB_OUT.mkdir(parents=True, exist_ok=True)


def save_embedding(pat, band, scale, coords, disc, hold):
    """Cache the 3-D embedding (per-phase x/y/z) + discriminants for later cool-plotting."""
    d = {}
    for ph in PHASES:
        x, y, z = coords[ph]
        d[f"{ph}_x"], d[f"{ph}_y"], d[f"{ph}_z"] = x, y, z
    meta = dict(patient=pat, band=band, scale=float(scale), **disc, **hold)
    np.savez_compressed(EMB_OUT / f"{pat}_{band}_s{scale:g}.npz",
                        meta=np.array(json.dumps(meta)), **d)

TCMAP = cm.plasma


def _load_xt(pat: str, phase: str) -> np.ndarray:
    X = np.asarray(load_timeseries(pat, phase, SEEG_DATAPATH), float)
    return X if X.shape[0] <= X.shape[1] else X.T


def _z(v: np.ndarray) -> np.ndarray:
    s = v.std()
    return (v - v.mean()) / s if s > 0 else v * 0.0


def windowed_coph(pat: str, fs: float, nperseg: int, bands: list[str],
                  scale: float, verbose: bool):
    """Per-window mst@0.20 cophenetic state at dimensionless scale s, per band.

    Returns states[band][phase] = (n_win, P) array of condensed cophenetic vectors."""
    wlen = int(round(W_SEC * fs))
    step = int(round(W_SEC * (1.0 - OVERLAP) * fs))
    states = {b: {ph: [] for ph in PHASES} for b in bands}
    for ph in PHASES:
        X = _load_xt(pat, ph)
        starts = list(range(0, X.shape[1] - wlen + 1, step))
        kept = 0
        for s0 in starts:
            freqs, ff = segment_ffts(X[:, s0:s0 + wlen], fs, nperseg, nperseg // 2, FMAX)
            if ff.shape[1] < MIN_SEG:
                continue
            cube = imcoh_abs_cube(ff)
            for b in bands:
                Wb = band_abs_average(cube, freqs, *BRAIN_BANDS[b])
                try:
                    Bk = mst_union_top_fraction(np.asarray(Wb, float), FRAC)
                    ev, V = laplacian_eig(Bk)
                    coph = cophenetic_at_scale(ev, V, scale)
                    states[b][ph].append(coph)
                except Exception:
                    states[b][ph].append(None)
            kept += 1
            del cube
        if verbose:
            print(f"    {ph:11s} {kept:3d} windows", flush=True)
    # drop failed windows, stack
    for b in bands:
        for ph in PHASES:
            rows = [r for r in states[b][ph] if r is not None]
            states[b][ph] = np.asarray(rows) if rows else np.empty((0, 0))
    return states


def embed_eir(states_b: dict, n_pc: int = 15):
    """Encoding/inference/residual 3-D embedding of one band's window cloud.

    Each window cophenetic vector is z-scored (unit pattern), then the cloud is
    PCA-denoised to its top ``n_pc`` components (kills per-window high-frequency
    noise, exactly as the audit_164 embedding did) before the fixed task-contrast
    axes are built. rest_pre = origin. Returns per-phase (x, y, z) window arrays and
    the phase-mean anchor points."""
    zrows = {ph: np.array([_z(states_b[ph][i]) for i in range(len(states_b[ph]))])
             if len(states_b[ph]) else np.empty((0, 0)) for ph in PHASES}
    if any(zrows[ph].size == 0 for ph in PHASES):
        return None
    counts = {ph: len(zrows[ph]) for ph in PHASES}
    allrows = np.vstack([zrows[ph] for ph in PHASES])
    mu = allrows.mean(0)
    Ac = allrows - mu
    K = int(min(n_pc, Ac.shape[0] - 1))
    _, _, Vt = np.linalg.svd(Ac, full_matrices=False)
    PC = Vt[:K]                                  # (K, P) principal directions
    red = Ac @ PC.T                              # (W, K) denoised coordinates
    rph, idx = {}, 0
    for ph in PHASES:
        rph[ph] = red[idx:idx + counts[ph]]
        idx += counts[ph]

    pre = rph["rest_pre"].mean(0)
    learn = rph["task_learn"].mean(0)
    test = rph["task_test"].mean(0)
    a1 = learn - pre
    a1 /= np.linalg.norm(a1) + 1e-9
    f = test - learn
    f -= (f @ a1) * a1
    a2 = f / (np.linalg.norm(f) + 1e-9)
    S0 = red - pre
    Res = S0 - np.outer(S0 @ a1, a1) - np.outer(S0 @ a2, a2)
    Resc = Res - Res.mean(0)
    _, _, Vt2 = np.linalg.svd(Resc, full_matrices=False)
    a3 = Vt2[0]

    def proj(M):
        d = M - pre
        return d @ a1, d @ a2, d @ a3

    coords, anchors = {}, {}
    for ph in PHASES:
        coords[ph] = proj(rph[ph])
        anchors[ph] = tuple(float(v.mean()) for v in coords[ph])

    # --- absolute (scatter-normalized) inter-phase separations, in the denoised
    #     PCA space -- NOT self-referential to the encoding/inference axes, so an
    #     "anchor" (task never moved the state) is detectable. ---------------------
    cen = {ph: rph[ph].mean(0) for ph in PHASES}
    within = np.mean([np.mean(np.sum((rph[ph] - cen[ph]) ** 2, axis=1))
                      for ph in PHASES])
    sigma = float(np.sqrt(within)) + 1e-9
    vt = cen["task_test"] - cen["rest_pre"]
    vp = cen["rest_post"] - cen["rest_pre"]
    disc = dict(
        D_learn=float(np.linalg.norm(cen["task_learn"] - cen["rest_pre"]) / sigma),
        D_task=float(np.linalg.norm(vt) / sigma),
        D_post=float(np.linalg.norm(vp) / sigma),
        align=float(np.dot(vt, vp) / (np.linalg.norm(vt) * np.linalg.norm(vp) + 1e-9)),
    )
    return coords, anchors, disc


def _hold(anchors: dict) -> dict:
    """Fraction of the task displacement retained at rest_post (>~1 hold, ~0 return)."""
    lx, ly, _ = anchors["task_learn"]
    tx, ty, _ = anchors["task_test"]
    px, py, _ = anchors["rest_post"]
    task = np.array([tx, ty])
    ext = np.linalg.norm(task) + 1e-9
    hold_task = float(np.dot([px, py], task / ext) / ext)
    hold_enc = float(px / lx) if abs(lx) > 1e-9 else np.nan
    hold_inf = float(py / ty) if abs(ty) > 1e-9 else np.nan
    return dict(hold_task=hold_task, hold_enc=hold_enc, hold_inf=hold_inf,
                enc_learn=lx, inf_test=ty, enc_post=px, inf_post=py)


CVAL = {"rest_pre": 0.0, "task_learn": 0.34, "task_test": 0.66, "rest_post": 1.0}
MK = {"rest_pre": "o", "task_learn": "^", "task_test": "s", "rest_post": "*"}


def _draw_traj(ax, coords, anchors, hold, sigma=1.2):
    """Render one time-coloured trajectory + phase anchors onto a 3-D axis."""
    order = PHASES
    xs = np.concatenate([coords[ph][0] for ph in order])
    ys = np.concatenate([coords[ph][1] for ph in order])
    zs = np.concatenate([coords[ph][2] for ph in order])
    t = np.arange(len(xs))
    norm = Normalize(t.min(), t.max())
    # faint raw window cloud
    ax.scatter(xs, ys, zs, c=t, cmap=TCMAP, norm=norm, s=6, alpha=0.28,
               depthshade=False, linewidths=0)
    # smoothed flow line
    xf = gaussian_filter1d(xs, sigma); yf = gaussian_filter1d(ys, sigma)
    zf = gaussian_filter1d(zs, sigma)
    pts = np.column_stack([xf, yf, zf]).reshape(-1, 1, 3)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = Line3DCollection(segs, cmap=TCMAP, norm=norm, linewidth=2.2, alpha=0.9)
    lc.set_array(t[:-1])
    ax.add_collection3d(lc)
    # phase-mean anchors
    for ph in order:
        ax_, ay_, az_ = anchors[ph]
        ax.scatter([ax_], [ay_], [az_], marker=MK[ph], s=150, color=TCMAP(CVAL[ph]),
                   edgecolor="k", linewidth=0.9, depthshade=False, zorder=6)
        ax.text(ax_, ay_, az_, "  " + PHLAB[ph], fontsize=8, zorder=7)
    ax.set_xlabel("encoding", fontsize=11, labelpad=-4)
    ax.set_ylabel("inference", fontsize=11, labelpad=-4)
    ax.set_zlabel("residual", fontsize=10, labelpad=-6)
    ax.tick_params(labelsize=6, pad=-2)
    ax.view_init(elev=18, azim=-60)


def _verdict(ht):
    return "HOLDS (tracer)" if ht > 0.5 else "RETURNS (resetter)" if ht < 0.25 else "partial"


def plot_traj(coords, anchors, pat, band, scale, hold):
    """One 3-D encoding/inference/residual trajectory, time-coloured."""
    fig = plt.figure(figsize=(5.2, 4.6))
    ax = fig.add_subplot(111, projection="3d")
    _draw_traj(ax, coords, anchors, hold)
    ht = hold["hold_task"]
    ax.set_title(rf"{pat} $\cdot$ {band}  (s={scale:g})" "\n"
                 rf"rest$_\mathrm{{post}}$ retains {ht:+.2f} of task drift $\to$ {_verdict(ht)}",
                 fontsize=9.5)
    fig.tight_layout()
    out = FIG_OUT / f"traj_{pat}_{band}_s{scale:g}.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    return out


def plot_compare(embA, embB, labelA, labelB, pat, scale):
    """Two bands side by side in the same grammar (e.g. beta tracer vs low_gamma null)."""
    fig = plt.figure(figsize=(9.6, 4.6))
    for i, (emb, lab) in enumerate([(embA, labelA), (embB, labelB)], 1):
        coords, anchors, hold = emb
        ax = fig.add_subplot(1, 2, i, projection="3d")
        _draw_traj(ax, coords, anchors, hold)
        ht = hold["hold_task"]
        ax.set_title(rf"{pat} $\cdot$ {lab}   rest$_\mathrm{{post}}$ retains {ht:+.2f}"
                     "\n" rf"$\to$ {_verdict(ht)}", fontsize=10)
    fig.tight_layout()
    out = FIG_OUT / f"compare_{pat}_{labelA}_vs_{labelB}_s{scale:g}.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=BANDS)
    ap.add_argument("--scale", type=float, default=5.6, help="dimensionless s = tau*lambda_max")
    ap.add_argument("--plot", action="store_true",
                    help="also emit the quick per-band trajectory PDFs (batch sweep skips them)")
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    N = len(args.patients)
    for pi, pat in enumerate(args.patients):
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        tp = time.time()
        print(f"[{pi+1}/{N}] {pat} (fs={fs:.0f} nperseg={nperseg} s={args.scale:g})", flush=True)
        states = windowed_coph(pat, fs, nperseg, args.bands, args.scale, verbose=True)
        embs = {}
        for b in args.bands:
            emb = embed_eir(states[b])
            if emb is None:
                print(f"    {b:11s} SKIP (empty phase)", flush=True)
                continue
            coords, anchors, disc = emb
            hold = _hold(anchors)
            embs[b] = (coords, anchors, hold)
            save_embedding(pat, b, args.scale, coords, disc, hold)
            rows.append(dict(patient=pat, band=b, scale=args.scale, **disc, **hold))
            if args.plot:
                plot_traj(coords, anchors, pat, b, args.scale, hold)
            print(f"    {b:11s} D_task={disc['D_task']:.2f} D_post={disc['D_post']:.2f} "
                  f"align={disc['align']:+.2f}  hold_task={hold['hold_task']:+.2f}", flush=True)
        if args.plot and "beta" in embs and "low_gamma" in embs:
            plot_compare(embs["beta"], embs["low_gamma"], "beta", "low_gamma", pat, args.scale)
        el = time.time() - tp
        print(f"    done {el:.1f}s  ETA {el*(N-pi-1)/60:.1f} min", flush=True)

    df = pd.DataFrame(rows)
    csv = DATA_OUT / "dynamics_discriminants.csv"
    if csv.exists():
        old = pd.read_csv(csv)
        df = pd.concat([old, df]).drop_duplicates(["patient", "band", "scale"], keep="last")
    df = df.sort_values(["scale", "band", "patient"])
    df.to_csv(csv, index=False)
    classify_and_report(df[df.scale == args.scale].copy())
    print(f"\nwrote {csv}\nembeddings -> {EMB_OUT}\ntotal {(time.time()-t0)/60:.1f} min", flush=True)


def classify_and_report(d):
    """Score every (patient, band) cell for the four dynamical archetypes and print
    the cleanest candidates. Annotate each with its REAL matched-strength trace
    (T_test) so the pick is honest about whether the raw dynamic matches the gate."""
    per = ROOT / "data" / "sparsified_arc" / "enc_inf_arc_mst020" / "per_cell.csv"
    if per.exists():
        pc = pd.read_csv(per)[["patient", "band", "T_test__obs_s1", "T_test__obs_max"]]
        d = d.merge(pc, on=["patient", "band"], how="left")
    else:
        d["T_test__obs_s1"] = d["T_test__obs_max"] = np.nan
    dt, dp, al = d["D_task"], d["D_post"], d["align"]
    d["anchor_score"] = -(dt + dp)                       # nothing moves
    d["reset_score"] = dt - dp                            # big task move, returns
    d["trace_score"] = np.where(dt > dt.median(), dp * al.clip(0), -9)      # holds, aligned
    d["reorg_score"] = np.where(dt > dt.median(), dp * (1 - al.clip(-1, 1)), -9)  # new place
    print("\n================  DYNAMICAL ARCHETYPES (raw trajectory; NOT the gate)  ================")
    print("  D_task=task moved state (sigmas)  D_post=rest_post dist from baseline  align=post·task")
    print("  T_test = the REAL matched-strength trace of that cell (for honesty)\n")
    for name, col in [("ANCHOR  (never moves)", "anchor_score"),
                      ("RESET   (moves, returns)", "reset_score"),
                      ("REORGANIZE (moves, new place)", "reorg_score"),
                      ("TRACE   (moves, holds aligned)", "trace_score")]:
        top = d.sort_values(col, ascending=False).head(3)
        print(f"--- {name} ---")
        for _, r in top.iterrows():
            print(f"    {r['patient']}/{r['band']:11s} D_task={r['D_task']:.2f} "
                  f"D_post={r['D_post']:.2f} align={r['align']:+.2f} | "
                  f"real T_test@s1={r['T_test__obs_s1']:+.2f} "
                  f"@max={r['T_test__obs_max']:+.2f}", flush=True)
        print()


if __name__ == "__main__":
    main()
