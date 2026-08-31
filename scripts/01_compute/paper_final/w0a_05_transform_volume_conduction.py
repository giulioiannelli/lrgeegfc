#!/usr/bin/env python3
"""W0-A / A1 extension — |ImCoh| vs ordinary coherence, with and without same-shaft pairs.

WHY THIS EXISTS. Sibling lane W0-B found that (a) the cross-phase trace fails
every lag-destroying, coherence-magnitude-preserving null, and (b) ordinary
coherence <|C|>_f yields a LARGER median rho_sym than <|Im C|>_f in both trace
bands (alpha +0.277 vs +0.212, beta +0.355 vs +0.281). Taken at face value that
says the project adopted a transform for volume-conduction immunity which is
weaker on its own headline statistic than the transform it rejected.

It must NOT be taken at face value, and this script is the control that decides
it. <|C|> is not VC-immune, and on sEEG it carries the 2-8x same-shaft bias that
is CLAUDE.md invariant 5. Anatomically fixed structure is trivially stable across
phases, so a VC-contaminated measure is EXPECTED to give a larger cross-phase
correlation for reasons that have nothing to do with the task. The larger number
may be the artifact rather than the improvement.

5-POINT CRITICAL PREAMBLE
1. Claim. The <|C|> advantage over <|Im C|> on the cross-phase trace is carried
   by same-shaft (volume-conduction / shared-reference) pairs, and does not
   survive their removal.
2. Null. H0: the advantage is genuine -- <|C|> still beats <|Im C|> after
   same-shaft pairs are removed from BOTH the graph and the statistic, i.e. the
   trace is carried by coupling magnitude rather than by lag.
3. Strongest plausible alternative the null should control for. That removing
   same-shaft pairs handicaps <|C|> for a reason unrelated to VC -- it simply
   deletes the edges where <|C|> has the most dynamic range, so ANY measure
   would lose. Controlled by applying the identical mask to all three transforms
   and comparing the DROP, not the level: the question is whether <|C|> loses
   more than <|Im C|> does, on the same pairs, the same nodes, the same
   surrogate draws.
4. Does the null control for it, by mechanism, and what it CANNOT reject.
   Same-shaft removal is applied at BOTH stages -- the edges are zeroed before
   sparsification (so the graph is not built from them) AND the pairs are
   dropped from the rho_sym correlation (so they cannot contribute through
   indirect cophenetic paths). All three transforms are derived from ONE Welch
   pass per phase, so the spectral estimate is bit-identical across arms and the
   only difference is the transform. What it CANNOT reject: same-shaft is a
   coarse proxy for volume conduction -- adjacent contacts on DIFFERENT shafts
   that are anatomically close are still included, so this is a lower bound on
   VC contamination, not its removal. It also cannot distinguish volume
   conduction from genuine short-range connectivity, which is the standard
   confound of every VC control.
5. Falsification / limits. The claim is falsified if <|C|> keeps its advantage
   after masking. If BOTH collapse, the trace depends on same-shaft structure,
   which would be the most consequential finding in the project and must be
   reported immediately rather than folded into a transform comparison. n = 10,
   matched-strength only; this is a transform comparison run under the incumbent
   null, not a re-test of the trace against the W0-B ladder.

DESIGN. Three transforms x two pair sets x three backbone fractions from the
locked contract window, five phases, four functionals, matched-strength null.
Every arm shares its Welch pass, its node set, its scale grid and (within an
arm) its surrogate draws, so the arms differ only in the transform and the mask.

Outputs (data/paper_final/w0a_substrate/a1b_volume_conduction/):
  per_patient_scale.csv   arm, transform, probe_mask, frac, patient, band, s, ...
  arm_summary.csv         cohort margin + cluster gate per arm x band x functional
  probe_bias.csv          same-shaft / cross-shaft weight ratio per transform
"""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, DEFAULT_SAMPLE_RATE, FS_OVERRIDES, PATIENTS_4PHASE, nperseg_for_fs,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction
from lrg_eegfc.utils.fc.heat_multiscale import (
    cross_phase_functionals_over_scales, laplacian_eig,
)
from lrg_eegfc.utils.io.patient import load_channel_labels, load_timeseries
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.surrogate.coherency_surrogate import complex_coherency_bands
from lrg_eegfc.workflow.substrate import CANONICAL, CANONICAL_SCALES

COHORT = list(PATIENTS_4PHASE)
BANDS = (os.environ.get("W0A_BANDS", "").split(",") if os.environ.get("W0A_BANDS")
         else ["theta", "alpha", "beta"])
PHASES = ("A", "B", "task_learn", "task_test", "rest_post")
TRANSFORMS = ("imcoh_abs", "imcoh_sq", "coh_abs")
MASKS = ("all", "xshaft")                       # all pairs / same-shaft removed
FRACS = tuple(CANONICAL.plateau_fracs[1:])      # (0.10, 0.14, 0.20) -- contract window
FUNCTIONALS = ("T_probe", "T_encode", "T_probespec", "T_probespec_pe")
R = int(os.environ.get("W0A_R", 100))
SWAP_FACTOR, W_MAX, BASE_SEED = 20, 1.0, 20260831
SGRID = CANONICAL_SCALES
OUT = ROOT / "data" / "paper_final" / "w0a_substrate" / "a1b_volume_conduction"


def band_adjacency(C: np.ndarray, transform: str) -> np.ndarray:
    """(F, N, N) complex coherency -> (N, N) adjacency, per-bin then band-averaged.

    ``imcoh_abs`` = mean|Im C| ; ``imcoh_sq`` = mean (Im C)^2 ; ``coh_abs`` =
    mean|C| (ordinary coherence magnitude, NOT volume-conduction immune).
    """
    if transform == "imcoh_abs":
        A = np.abs(C.imag).mean(0)
    elif transform == "imcoh_sq":
        A = (C.imag ** 2).mean(0)
    elif transform == "coh_abs":
        A = np.abs(C).mean(0)
    else:                                                       # pragma: no cover
        raise ValueError(transform)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return np.clip(A, 0.0, 1.0)


def per_patient(job):
    """Wrapper so one bad patient cannot kill a 20-minute cohort run."""
    try:
        return _per_patient(job)
    except Exception:
        import traceback
        print(f"  [{job[1]}] FAILED:\n{traceback.format_exc()}", flush=True)
        return None, None


def _per_patient(job):
    idx, pat = job
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nper_full, nper_half = nperseg_for_fs(fs), max(256, nperseg_for_fs(fs) // 2)
    bands = {b: BRAIN_BANDS[b] for b in BANDS}

    # --- one Welch pass per phase; all three transforms share it -------------- #
    C = {}
    try:
        X = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
        if X.shape[0] > X.shape[1]:
            X = X.T
        T = X.shape[1]
        C["A"] = complex_coherency_bands(np.ascontiguousarray(X[:, : T // 2]), fs,
                                         bands, nper_half)
        C["B"] = complex_coherency_bands(np.ascontiguousarray(X[:, T // 2:]), fs,
                                         bands, nper_half)
        del X
        for ph in ("task_learn", "task_test", "rest_post"):
            Y = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
            if Y.shape[0] > Y.shape[1]:
                Y = Y.T
            C[ph] = complex_coherency_bands(np.ascontiguousarray(Y), fs, bands, nper_full)
            del Y
    except Exception as exc:                                    # pragma: no cover
        print(f"  [{pat}] load failed: {exc}", flush=True)
        return None, None

    N = C["A"][BANDS[0]].shape[1]
    labels = load_channel_labels(pat)
    if labels is None or len(labels) != N:
        print(f"  [{pat}] channel_labels missing/misaligned (N={N}) — skip", flush=True)
        return None, None
    same = build_probe_mask(labels)                             # (N, N) same-shaft
    iu = np.triu_indices(N, 1)
    same_cond = same[iu]
    keep_cond = ~same_cond                                      # cross-shaft pairs

    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rows, bias_rows = [], []

    for band in BANDS:
        Wd = {tk: {ph: band_adjacency(C[ph][band], tk) for ph in PHASES}
              for tk in TRANSFORMS}
        # same-shaft bias, per transform (CLAUDE.md invariant 5)
        for tk in TRANSFORMS:
            w = Wd[tk]["task_test"][iu]
            bias_rows.append(dict(patient=pat, band=band, transform=tk,
                                  same_mean=float(w[same_cond].mean()),
                                  cross_mean=float(w[keep_cond].mean()),
                                  ratio=float(w[same_cond].mean() / w[keep_cond].mean())
                                  if w[keep_cond].mean() > 0 else np.nan,
                                  n_same=int(same_cond.sum()), n_cross=int(keep_cond.sum())))

        for tk in TRANSFORMS:
            for mask_name in MASKS:
                if mask_name == "all":
                    Wm = Wd[tk]
                    pmask = None
                else:                       # zero same-shaft EDGES and drop the PAIRS
                    Wm = {}
                    for ph in PHASES:
                        A = Wd[tk][ph].copy()
                        A[same] = 0.0
                        Wm[ph] = A
                    pmask = keep_cond
                obs = {f: cross_phase_functionals_over_scales(
                          {ph: laplacian_eig(mst_union_top_fraction(Wm[ph], f))
                           for ph in PHASES}, SGRID, pair_mask=pmask)
                       for f in FRACS}
                surr = {f: {k: np.full((R, SGRID.size), np.nan) for k in FUNCTIONALS}
                        for f in FRACS}
                # Deterministic arm seed. Python's hash() on str is salted per
                # process (PYTHONHASHSEED), so it must NOT appear in a seed --
                # it would make the surrogate draws irreproducible across runs
                # and inconsistent between workers.
                arm_id = (BANDS.index(band) * len(TRANSFORMS) * len(MASKS)
                          + TRANSFORMS.index(tk) * len(MASKS) + MASKS.index(mask_name))
                rs = np.random.default_rng(BASE_SEED + 7919 * idx + arm_id)
                for r in range(R):
                    Wsh = {ph: matched_strength_shuffle(Wm[ph], n_swaps, rs, W_MAX)
                           for ph in PHASES}
                    for f in FRACS:
                        so = cross_phase_functionals_over_scales(
                            {ph: laplacian_eig(mst_union_top_fraction(Wsh[ph], f))
                             for ph in PHASES}, SGRID, pair_mask=pmask)
                        for k in FUNCTIONALS:
                            surr[f][k][r] = so[k]
                for f in FRACS:
                    for k in FUNCTIONALS:
                        o_all, s_all = obs[f][k], surr[f][k]
                        for j, s in enumerate(SGRID):
                            col = s_all[:, j]
                            col = col[np.isfinite(col)]
                            o = o_all[j]
                            rows.append(dict(
                                patient=pat, band=band, transform=tk,
                                probe_mask=mask_name, frac=f, functional=k,
                                s=float(s), N=int(N),
                                obs_rho=float(o) if np.isfinite(o) else np.nan,
                                surr_p50=float(np.nanpercentile(col, 50)) if col.size else np.nan,
                                p=float(np.mean(col >= o)) if col.size and np.isfinite(o) else np.nan))
        del Wd
    return rows, bias_rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ncpu = int(os.environ.get("SA_WORKERS", 6))
    jobs = [(i, p) for i, p in enumerate(COHORT)]
    print(f"[a1b] {len(jobs)} patients x {len(BANDS)} bands x {len(TRANSFORMS)} transforms "
          f"x {len(MASKS)} masks x {len(FRACS)} fracs, R={R}, {ncpu} workers -> {OUT}",
          flush=True)
    print(f"[a1b] transforms={TRANSFORMS} masks={MASKS} fracs={FRACS} bands={BANDS}",
          flush=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    t0 = time.time()
    rows, bias = [], []
    with Pool(ncpu) as pool:
        for i, (rl, bl) in enumerate(pool.imap_unordered(per_patient, jobs), 1):
            if rl:
                rows.extend(rl); bias.extend(bl)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
            pd.DataFrame(rows).to_csv(OUT / "per_patient_scale.csv", index=False)
            pd.DataFrame(bias).to_csv(OUT / "probe_bias.csv", index=False)
    df, bdf = pd.DataFrame(rows), pd.DataFrame(bias)
    print(f"\n[a1b] {len(df)} rows in {time.time()-t0:.0f}s", flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="a1b_transform_volume_conduction", cohort=COHORT, bands=BANDS,
        transforms=list(TRANSFORMS), masks=list(MASKS), fracs=[float(f) for f in FRACS],
        R=R, phases=list(PHASES), s_grid=[float(x) for x in SGRID], base_seed=BASE_SEED,
        note="same-shaft removal applied to BOTH the graph (edges zeroed pre-sparsify) "
             "and the statistic (pairs dropped from rho_sym)"), indent=2))

    print("\n=== same-shaft bias per transform (CLAUDE.md invariant 5) ===", flush=True)
    for band in BANDS:
        for tk in TRANSFORMS:
            x = bdf[(bdf.band == band) & (bdf.transform == tk)]
            if x.empty:
                continue
            print(f"  {band:6s} {tk:10s} same/cross weight ratio = {x.ratio.median():.2f}x "
                  f"[{x.ratio.min():.2f}, {x.ratio.max():.2f}]", flush=True)

    print("\n=== the head-to-head: cohort median MARGIN (obs - surr_p50), "
          "median over scales and fractions ===", flush=True)
    out = []
    for k in FUNCTIONALS:
        for band in BANDS:
            for tk in TRANSFORMS:
                for mk in MASKS:
                    x = df[(df.functional == k) & (df.band == band)
                           & (df.transform == tk) & (df.probe_mask == mk)]
                    if x.empty:
                        continue
                    g = x.groupby(["patient", "s", "frac"]).first().reset_index()
                    m = g.obs_rho - g.surr_p50
                    out.append(dict(functional=k, band=band, transform=tk, probe_mask=mk,
                                    obs_med=float(g.obs_rho.median()),
                                    surr_med=float(g.surr_p50.median()),
                                    margin_med=float(m.median()),
                                    frac_pos=float((m > 0).mean())))
    sm = pd.DataFrame(out)
    sm.to_csv(OUT / "arm_summary.csv", index=False)
    for k in ("T_probe", "T_encode"):
        print(f"\n  --- [{k}] obs (raw rho_sym) | margin vs matched-strength ---", flush=True)
        for band in BANDS:
            print(f"    {band}", flush=True)
            for tk in TRANSFORMS:
                cells = []
                for mk in MASKS:
                    y = sm[(sm.functional == k) & (sm.band == band)
                           & (sm.transform == tk) & (sm.probe_mask == mk)]
                    if y.empty:
                        continue
                    cells.append(f"{mk:7s} obs={y.obs_med.iloc[0]:+.3f} "
                                 f"margin={y.margin_med.iloc[0]:+.3f}")
                print(f"      {tk:10s} " + " | ".join(cells), flush=True)
            # the decisive quantity: how much each transform LOSES to masking
            for tk in TRANSFORMS:
                a = sm[(sm.functional == k) & (sm.band == band) & (sm.transform == tk)
                       & (sm.probe_mask == "all")]
                b = sm[(sm.functional == k) & (sm.band == band) & (sm.transform == tk)
                       & (sm.probe_mask == "xshaft")]
                if a.empty or b.empty:
                    continue
                print(f"        {tk:10s} DROP on masking: obs "
                      f"{a.obs_med.iloc[0]:+.3f} -> {b.obs_med.iloc[0]:+.3f} "
                      f"({b.obs_med.iloc[0]-a.obs_med.iloc[0]:+.3f}); margin "
                      f"{a.margin_med.iloc[0]:+.3f} -> {b.margin_med.iloc[0]:+.3f} "
                      f"({b.margin_med.iloc[0]-a.margin_med.iloc[0]:+.3f})", flush=True)
    print(f"\n[a1b] done -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
