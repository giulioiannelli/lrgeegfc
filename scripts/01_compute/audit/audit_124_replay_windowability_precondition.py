#!/usr/bin/env python3
"""Audit 124 — N4 replay-states PRECONDITION: |ImCoh| windowability + N1 anchor.

This is STEP 0 of the replay-states verification brief
(`.agents/guides/task-persistence-investigation/
2026-06-22_replay-states-multiscale-reinstatement.md`). NO replay claim is made
here. Two things, both gating preconditions for everything downstream:

(A) CONTINUITY ANCHOR. The windowed |ImCoh| recipe used by the replay pipeline
    must reduce to the canonical N1 `imcoh_abs` in the full-phase limit, else the
    time-resolved measure is not the time-resolved version of the N1 trace. We
    recompute full-phase |ImCoh| from raw timeseries with the windowed recipe
    (per-bin abs, then band-mean; nperseg = nperseg_for_fs(fs)) and Spearman-
    correlate the pair vector against the cached `load_fc_matrix(..., imcoh_abs)`.
    Expect rho ~ 1.0.

(B) WINDOWABILITY SWEEP. |ImCoh| needs enough Welch segments to be a stable
    estimate; short windows + low bands may be unusable (preamble alternative C:
    estimator artifact manufactures spurious configuration excursions). For each
    (patient, phase, band) and candidate window length L, we measure the
    WITHIN-WINDOW split-half reliability of the |ImCoh| edge vector: split the
    window's Welch segments into even/odd interleaved subsets, compute
    |ImCoh|_band on each subset independently, and Spearman-correlate the two
    upper-triangular pair vectors. High r(L) => the connectivity structure is
    reproducible from independent halves of the data inside the window =>
    estimator-stable. (Interleaving removes slow drift, so this isolates
    estimator noise, not non-stationarity.) Split-half uses HALF the segments, so
    it is a conservative lower bound on full-window reliability.

We then read off, per band, the shortest L whose reliability plateaus -> L*[b],
and the usable band set B*. NO acceptance gate is pre-registered
(feedback_no_pre_registered_acceptance): we plot r(L,band) and pick the knee with
the PI.

Outputs
-------
    data/audit/replay_windowability/
        anchor.csv                 # full-phase recipe vs cached imcoh_abs
        reliability.csv            # r(L, band, patient, phase) split-half
        figures/reliability.pdf    # r(L) curves per band
        README.md
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BAND_TEX_DICT, FS_OVERRIDES, nperseg_for_fs,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.coherence import (
    compute_imcoh, segment_ffts, imcoh_abs_cube, band_abs_average,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.visuals.styles import use_lrg_style

use_lrg_style()

OUT = ROOT / "data" / "audit" / "replay_windowability"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

BANDS = list(BRAIN_BANDS.keys())          # all six, to honestly report usability
L_GRID = [8.0, 12.0, 16.0, 20.0, 30.0, 45.0, 60.0]
FMAX_KEEP = 305.0                          # cap freq bins computed (high_gamma->300)


# ---------------------------------------------------------------------------
# windowed |ImCoh| primitives (segment-level, so even/odd split is possible)
# ---------------------------------------------------------------------------
def full_phase_imcoh(X: np.ndarray, fs: float, lo: float, hi: float, nperseg: int) -> np.ndarray:
    """Full-phase windowed-recipe |ImCoh|_band (the anchor object), via the
    library's memory-batched compute_imcoh."""
    freqs, signed = compute_imcoh(X, fs, nperseg=nperseg, noverlap=nperseg // 2)
    mask = (freqs >= lo) & (freqs <= hi)
    W = np.abs(signed[:, :, mask]).mean(axis=-1)
    W = 0.5 * (W + W.T)
    np.fill_diagonal(W, 0.0)
    return W


def triu(W: np.ndarray) -> np.ndarray:
    iu = np.triu_indices(W.shape[0], k=1)
    return W[iu]


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=["Pat_05", "Pat_02"])
    ap.add_argument("--phases", nargs="+", default=["rest_post", "task_test"])
    ap.add_argument("--n-windows", type=int, default=12,
                    help="windows sampled per (patient,phase,L) for reliability")
    ap.add_argument("--anchor-band", default="beta")
    args = ap.parse_args()

    t0 = time.time()
    anchor_rows = []
    rel_rows = []

    for pat in args.patients:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)
        print(f"\n=== {pat} (fs={fs:.0f}, nperseg={nperseg}) ===")

        # ---- (A) continuity anchor on a couple of phases/bands ----
        for ph in args.phases:
            try:
                X = load_timeseries(pat, ph, SEEG_DATAPATH)
            except Exception as e:
                print(f"  [load fail] {pat}/{ph}: {e}")
                continue
            X = np.asarray(X, dtype=np.float64)
            if X.shape[0] > X.shape[1]:
                X = X.T
            N, T = X.shape

            for band in [args.anchor_band, "alpha"]:
                lo, hi = BRAIN_BANDS[band]
                cached = load_fc_matrix(pat, ph, band, fc_method="imcoh_abs")
                if cached is None:
                    print(f"  [anchor] no cache {pat}/{ph}/{band}")
                    continue
                cached = np.asarray(cached, dtype=float)
                W_recipe = full_phase_imcoh(X, fs, lo, hi, nperseg)
                if W_recipe.shape != cached.shape:
                    print(f"  [anchor] shape mismatch {pat}/{ph}/{band}: "
                          f"{W_recipe.shape} vs {cached.shape}")
                    continue
                rho, _ = spearmanr(triu(W_recipe), triu(cached))
                pear = np.corrcoef(triu(W_recipe), triu(cached))[0, 1]
                anchor_rows.append(dict(patient=pat, phase=ph, band=band,
                                        N=N, spearman=rho, pearson=pear))
                print(f"  [anchor] {ph}/{band}: spearman={rho:.4f} pearson={pear:.4f}")

            # ---- (B) windowability sweep on this phase ----
            for L in L_GRID:
                wlen = int(round(L * fs))
                if wlen > T:
                    continue
                starts = np.linspace(0, T - wlen, args.n_windows, dtype=int)
                rs = {b: [] for b in BANDS}
                nseg_used = None
                for s in starts:
                    win = X[:, s:s + wlen]
                    freqs, ff = segment_ffts(win, fs, nperseg, nperseg // 2, FMAX_KEEP)
                    nseg = ff.shape[1]
                    nseg_used = nseg
                    if nseg < 4:
                        continue
                    even = np.arange(0, nseg, 2)
                    odd = np.arange(1, nseg, 2)
                    cube_e = imcoh_abs_cube(ff, even)
                    cube_o = imcoh_abs_cube(ff, odd)
                    if cube_e is None or cube_o is None:
                        continue
                    for band in BANDS:
                        lo, hi = BRAIN_BANDS[band]
                        We = band_abs_average(cube_e, freqs, lo, hi)
                        Wo = band_abs_average(cube_o, freqs, lo, hi)
                        if We is None or Wo is None:
                            continue
                        r, _ = spearmanr(triu(We), triu(Wo))
                        if np.isfinite(r):
                            rs[band].append(r)
                for band in BANDS:
                    if rs[band]:
                        rel_rows.append(dict(
                            patient=pat, phase=ph, band=band, L_sec=L,
                            n_seg=nseg_used, n_win=len(rs[band]),
                            r_mean=float(np.mean(rs[band])),
                            r_std=float(np.std(rs[band])),
                            r_min=float(np.min(rs[band])),
                        ))
                print(f"  [L={L:4.0f}s nseg={nseg_used:3d}] " +
                      "  ".join(f"{b[:5]}={np.mean(rs[b]):.2f}" if rs[b] else f"{b[:5]}=  -"
                               for b in BANDS))

    import pandas as pd
    anchor = pd.DataFrame(anchor_rows)
    rel = pd.DataFrame(rel_rows)
    anchor.to_csv(OUT / "anchor.csv", index=False)
    rel.to_csv(OUT / "reliability.csv", index=False)

    print("\n=== ANCHOR (windowed recipe vs cached imcoh_abs) ===")
    print(anchor.to_string(index=False))
    print("\n=== RELIABILITY (split-half Spearman of |ImCoh| edge vector) ===")
    piv = rel.pivot_table(index=["band"], columns="L_sec", values="r_mean", aggfunc="mean")
    piv = piv.reindex(BANDS)
    print(piv.to_string(float_format=lambda v: f"{v:.2f}"))

    # ---- figure: r(L) per band, averaged over patients/phases ----
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    cmap = plt.get_cmap("turbo")
    for i, band in enumerate(BANDS):
        sub = rel[rel.band == band]
        if sub.empty:
            continue
        g = sub.groupby("L_sec")["r_mean"].mean()
        ax.plot(g.index, g.values, "-o", color=cmap(i / (len(BANDS) - 1)),
                label=BRAIN_BAND_TEX_DICT.get(band, band), ms=4)
    ax.axhline(0.8, color="0.5", ls="--", lw=0.8)
    ax.set_xlabel("window length L (s)")
    ax.set_ylabel(r"split-half reliability $r$ of $|\mathrm{ImCoh}|$")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "reliability.pdf", transparent=True)
    plt.close(fig)

    print(f"\n[audit_124] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
