#!/usr/bin/env python3
"""W0-B rungs N3 (block phase-label permutation) and N4 (within-rest placebo).

The only rungs that reach session nonstationarity, slow drift, artifact epochs
and the phase segmentation itself. The session ``rest_pre + task_test +
rest_post`` is concatenated in true temporal order, cut into contiguous blocks
much longer than ``nperseg``, and the blocks are reassigned to pseudo-phases
{A, B, task_test, rest_post} preserving the true per-phase BLOCK COUNTS -- so
every surrogate has the same spectral degrees of freedom as the observed and the
comparison is never a segment-count artifact. FC and ``rho_sym`` are then
recomputed from scratch for each surrogate.

Variants
--------
n3_free   uniformly random block -> pseudo-phase assignment. Destroys temporal
          order entirely: no drift, no artifact clustering, no session shape.
n3_order  ★ order-preserving. The session is treated as a ring and the phase
          labels are rotated along it by ``r`` blocks; every pseudo-phase stays a
          contiguous interval of the true length, local temporal order survives
          everywhere but the single seam, and any smooth session trend is
          RETAINED. Only the placement of the cut points moves. The realization
          set is the rotation group, so all ``n_blocks`` rotations are enumerated
          and the test is EXACT rather than Monte-Carlo.
n4        within-rest placebo: all four pseudo-phases drawn from ``rest_pre``
          alone, block counts scaled to preserve the true duration RATIOS. Its
          absolute segment counts are ~4x smaller than the observed, which makes
          each pseudo-FC noisier, so n4 is reported against a DURATION-MATCHED
          observed comparator (the real phases truncated to the same block
          counts) computed through the identical code path. Without that
          comparator n4 would be confounded by spectral degrees of freedom.

Block length. 30 s = 15 x nperseg at 2048 Hz, ~29 Welch segments per block: long
enough that a block's own cross-spectrum is estimable, short enough that a ~2350 s
session yields ~78 blocks, i.e. a large free-permutation entropy and a 78-point
exact rotation test. Segments straddling a block boundary are dropped (marked
-1), so no segment ever mixes two blocks and no surrogate can inherit a spurious
cross-block correspondence.

================== 5-point critical preamble (pre-registered) ==================

(1) CLAIM. ``rho_sym`` computed on the TRUE phase partition exceeds what the
    same session yields under alternative partitions -- i.e. the geometry is tied
    to the task/rest boundaries, not to the session.

(2) NULL. "Any partition of this session into four intervals of these durations
    produces this geometry." n3_free removes temporal structure; n3_order keeps
    it and moves only the boundaries; n4 removes the task entirely.

(3) STRONGEST ALTERNATIVE. A slow session-level nonstationarity -- electrode
    drift, impedance change, arousal decline, a run of artifact epochs -- makes
    temporally adjacent windows more similar than distant ones. Because
    ``task_test`` sits BETWEEN ``rest_pre`` and ``rest_post``, such a gradient
    alone produces a positive ``rho_sym`` with no task-induced reorganization
    whatsoever. Neither matched-strength nor N1/N2 can touch this: all three keep
    the true phase boundaries and only perturb what is inside them.

(4) DOES IT CONTROL FOR IT, BY MECHANISM -- AND WHAT IT CANNOT REJECT.
    n3_order is the mechanism-matched control: it preserves the session's actual
    temporal profile, its drift and its artifact distribution, and preserves the
    contiguity and durations of the four intervals; the ONLY thing it changes is
    where the boundaries fall. If a smooth gradient were sufficient, a rotated
    partition would score as highly as the true one and the observed value would
    sit inside the rotation null. n3_free is the looser bracket.
    This is NOT the retired drift null. That one tested whether a monotone trend
    could manufacture a distance ASYMMETRY and was degenerate because a
    directional task makes drift and trace the same object; the statistic here is
    ``rho_sym``, a symmetric split-half correlation, and the manipulation is the
    boundary placement, not a detrending.
    CANNOT REJECT: (i) a nonstationarity whose change point HAPPENS to coincide
    with the true task onset -- rotation moves the boundary away from it, so such
    an event is indistinguishable from a task effect by design, and no
    within-session null can separate them; (ii) block quantization means
    boundaries move in 30 s steps, so a rotation of +-1 block is nearly the true
    partition and inflates the null's upper tail (conservative);
    (iii) n3 changes the phase CONTENT, so it is not a lag test -- it cannot say
    whether the trace is carried by lag or by magnitude. That is N1's job.
    Read the ladder, not the rung.

(5) FALSIFICATION. If the observed ``rho_sym`` sits inside the n3_order rotation
    null for beta at all scales, the trace is a property of the session rather
    than of the task, and the central claim fails. Limitation regardless: n = 10,
    a coarse cohort p-grid, and an exact per-patient test with only ~78
    rotations (per-patient p floor ~1/79).

Outputs: data/paper_final/w0b_nulls/{n3_free,n3_order,n4}/...
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))

from lrg_eegfc.config.const import (
    BRAIN_BANDS, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales
from lrg_eegfc.utils.surrogate.timeseries_nulls import (
    block_partition_indices, blocks_to_segment_sets, coherency_from_csd,
    csd_from_segment_subset, imcoh_abs_from_coherency, segment_fft,
)
from audit_150_rho_sym_gate import COHORT, BANDS

PHASES = ("A", "B", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)
FRAC, BACKBONE = 0.20, "mst020"
R_FREE = 200
BLOCK_SECONDS = 30.0
BASE_SEED = 20260825
OUT_ROOT = ROOT / "data" / "paper_final" / "w0b_nulls"


def readout(W_by_phase):
    return rho_sym_over_scales(
        {ph: laplacian_eig(select_backbone(W_by_phase[ph], BACKBONE, frac=FRAC))
         for ph in PHASES}, SGRID)


def build_session(pat, band):
    """Concatenated-session per-segment FFTs + block labels + true phase sizes.

    Returns (F, scale, block_of_segment, sizes, n_blocks). ``block_of_segment``
    is -1 for any segment that straddles a block boundary (dropped).
    """
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs)
    step = nps // 2
    parts, lens = [], {}
    for ph in ("rest_pre", "task_test", "rest_post"):
        X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
        if X.shape[0] > X.shape[1]:
            X = X.T
        parts.append(X)
        lens[ph] = X.shape[1]
    Xs = np.concatenate(parts, axis=1)
    del parts
    freqs, F, scale = segment_fft(Xs, fs, nps, band=BRAIN_BANDS[band])
    del Xs
    n_seg = F.shape[1]

    blk = int(round(BLOCK_SECONDS * fs))
    seg_start = np.arange(n_seg) * step
    seg_end = seg_start + nps
    b0 = seg_start // blk
    b1 = (seg_end - 1) // blk
    block_of_segment = np.where(b0 == b1, b0, -1).astype(np.int64)   # drop straddlers
    n_blocks = int(block_of_segment.max()) + 1

    # true per-phase block counts (a block belongs to the phase its start is in)
    bnd_pre = lens["rest_pre"]
    bnd_task = bnd_pre + lens["task_test"]
    blk_start = np.arange(n_blocks) * blk
    n_pre = int((blk_start < bnd_pre).sum())
    n_task = int(((blk_start >= bnd_pre) & (blk_start < bnd_task)).sum())
    n_post = n_blocks - n_pre - n_task
    sizes = {"A": n_pre // 2, "B": n_pre - n_pre // 2,
             "task_test": n_task, "rest_post": n_post}
    return F, scale, block_of_segment, sizes, n_blocks


def rho_from_blocks(F, scale, block_of_segment, block_idx):
    segsets = blocks_to_segment_sets(block_of_segment, block_idx)
    W = {}
    for ph in PHASES:
        idx = segsets[ph]
        if idx.size < 8:
            return np.full(SGRID.size, np.nan)
        W[ph] = imcoh_abs_from_coherency(
            coherency_from_csd(csd_from_segment_subset(F, np.sort(idx), scale)))
    return readout(W)


def per_cell(job):
    idx, pat, band, variant = job
    t0 = time.time()
    try:
        F, scale, bos, sizes, nb = build_session(pat, band)
    except Exception as e:
        return None, f"{pat}/{band}/{variant}: build failed ({e})"
    rng = np.random.default_rng([BASE_SEED, idx])

    if variant == "n4":
        # all four pseudo-phases inside rest_pre, true duration RATIOS preserved
        n_pre = sizes["A"] + sizes["B"]
        tot = sum(sizes.values())
        frac = {k: sizes[k] / tot for k in sizes}
        s4 = {k: max(2, int(round(frac[k] * n_pre))) for k in sizes}
        while sum(s4.values()) > n_pre:
            s4[max(s4, key=s4.get)] -= 1
        while sum(s4.values()) < n_pre:
            s4[max(s4, key=s4.get)] += 1
        pre_blocks = np.arange(n_pre)
        # Duration-matched OBSERVED comparator: each TRUE phase truncated to the
        # n4 block count, so observed and placebo share spectral degrees of
        # freedom and the contrast is not a segment-count artifact.
        off, starts = 0, {}
        for k in PHASES:
            starts[k] = off
            off += sizes[k]
        obs_idx = {k: np.arange(starts[k], starts[k] + s4[k]) for k in PHASES}
        obs = rho_from_blocks(F, scale, bos, obs_idx)
        surr = np.full((R_FREE, SGRID.size), np.nan)
        for r in range(R_FREE):
            perm = rng.permutation(pre_blocks)
            k0, bi = 0, {}
            for k in PHASES:
                bi[k] = perm[k0:k0 + s4[k]]; k0 += s4[k]
            surr[r] = rho_from_blocks(F, scale, bos, bi)
        n_real = R_FREE
    else:
        obs = rho_from_blocks(
            F, scale, bos, block_partition_indices(nb, sizes, rng, mode="identity"))
        if variant == "n3_free":
            n_real = R_FREE
            surr = np.full((n_real, SGRID.size), np.nan)
            for r in range(n_real):
                surr[r] = rho_from_blocks(
                    F, scale, bos, block_partition_indices(nb, sizes, rng, mode="free"))
        elif variant == "n3_order":
            # EXACT: enumerate the whole rotation group, excluding the identity
            rots = [r for r in range(nb) if r != 0]
            n_real = len(rots)
            surr = np.full((n_real, SGRID.size), np.nan)
            for i, r in enumerate(rots):
                surr[i] = rho_from_blocks(
                    F, scale, bos,
                    block_partition_indices(nb, sizes, rng, mode="order_preserving",
                                            rotation=r))
        else:
            raise ValueError(variant)

    rows = []
    for j, s in enumerate(SGRID):
        col = surr[:, j][np.isfinite(surr[:, j])]
        o = obs[j]
        if variant == "n3_order":                       # exact test -> (k+1)/(n+1)
            p = float((np.sum(col >= o) + 1) / (col.size + 1)) if col.size and np.isfinite(o) else np.nan
        else:
            p = float(np.mean(col >= o)) if col.size and np.isfinite(o) else np.nan
        rows.append(dict(patient=pat, band=band, rung=variant, s=float(s),
                         n_blocks=int(nb), block_s=BLOCK_SECONDS,
                         obs_rho=float(o) if np.isfinite(o) else np.nan,
                         surr_p50=float(np.percentile(col, 50)) if col.size else np.nan,
                         surr_p95=float(np.percentile(col, 95)) if col.size else np.nan,
                         p=p, n_surr=int(col.size)))
    return rows, (f"{pat}/{band}/{variant} nb={nb} sizes={list(sizes.values())} "
                  f"obs[s=1]={obs[0]:+.3f} p50[s=1]={np.nanpercentile(surr[:,0],50):+.3f} "
                  f"({time.time()-t0:.0f}s)")


def cohort_gate(df):
    out = []
    for rung in sorted(df.rung.unique()):
        for band in sorted(df.band.unique()):
            for s in np.sort(df.s.unique()):
                x = df[(df.rung == rung) & (df.band == band) & np.isclose(df.s, s)]
                x = x.dropna(subset=["obs_rho", "surr_p50"])
                if len(x) < 5:
                    continue
                d = x.obs_rho.values - x.surr_p50.values
                try:
                    p = float(wilcoxon(d, alternative="greater")[1]) if np.any(d != 0) else np.nan
                except Exception:
                    p = np.nan
                out.append(dict(rung=rung, band=band, s=float(s), gate_p=p,
                                n_pat=len(x), n_above=int((x.p < 0.05).sum()),
                                obs_med=float(x.obs_rho.median()),
                                surr_med=float(x.surr_p50.median()),
                                margin_med=float(np.median(d))))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", default="n3_free,n3_order,n4")
    ap.add_argument("--bands", default="delta,theta,alpha,beta,low_gamma")
    ap.add_argument("--patients", default="")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("NL_WORKERS", 8)))
    a = ap.parse_args()

    variants = a.variants.split(",")
    bands = a.bands.split(",") if a.bands else BANDS
    pats = a.patients.split(",") if a.patients else COHORT
    jobs = [(i, p, b, v) for i, (v, b, p) in
            enumerate((v, b, p) for v in variants for b in bands for p in pats)]
    print(f"[w0b-n3] {len(jobs)} cells | variants={variants} bands={bands} "
          f"n_pat={len(pats)} block={BLOCK_SECONDS}s workers={a.workers}", flush=True)

    t0, rows = time.time(), []
    with Pool(a.workers) as pool:
        for i, (rl, msg) in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                rows.extend(rl)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {msg} | {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)

    df = pd.DataFrame(rows)
    if df.empty:
        print("[w0b-n3] no rows"); return
    for v in variants:
        d = df[df.rung == v]
        if d.empty:
            continue
        od = OUT_ROOT / v
        od.mkdir(parents=True, exist_ok=True)
        d.to_csv(od / "per_patient_scale.csv", index=False)
        g = cohort_gate(d)
        g.to_csv(od / "cohort_gate.csv", index=False)
        (od / "config.json").write_text(json.dumps(dict(
            rung=v, backbone=BACKBONE, frac=FRAC, block_seconds=BLOCK_SECONDS,
            R=R_FREE if v != "n3_order" else "exact (all rotations)",
            s_grid=[float(x) for x in SGRID], cohort=pats, bands=bands,
            base_seed=BASE_SEED), indent=2))
        print(f"\n=== {v}: cohort gate, per band (min gate_p over scales) ===", flush=True)
        if g.empty or "band" not in g.columns:
            print("  (needs >=5 patients)", flush=True); continue
        for band in bands:
            gb = g[g.band == band].dropna(subset=["gate_p"])
            if gb.empty:
                continue
            b = gb.loc[gb.gate_p.idxmin()]
            print(f"  {band:11s} best s={b.s:6.1f} gate_p={b.gate_p:.4f} "
                  f"clears {int((gb.gate_p<0.05).sum())}/{len(gb)} scales "
                  f"obs_med={b.obs_med:+.3f} surr_med={b.surr_med:+.3f}", flush=True)
    print(f"\n[w0b-n3] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT_ROOT}", flush=True)


if __name__ == "__main__":
    main()
