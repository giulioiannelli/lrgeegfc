#!/usr/bin/env python3
"""W0-B rungs N1 / N2 / N1b -- coherency- and segment-level nulls for the
cross-phase trace, per patient, per band, per scale.

Pipeline-identical to ``sparsified_arc/13_matched_strength_mst020.py`` from the
adjacency onward: dense ``imcoh_abs`` -> ``mst_union_top_fraction(frac=0.20)``
-> combinatorial Laplacian -> heat kernel at ``s = tau*lambda_max`` -> UPGMA ->
cophenetic -> ``rho_sym``. The ONLY difference is where the randomness enters:
script 13 shuffles the finished ``N x N`` matrix; this script perturbs the
complex band coherency ``C(f)`` or the Welch segment lattice, i.e. upstream of
the estimator.

The observed statistic is recomputed here from the timeseries through the same
door as the surrogates (never read from a CSV) and is checked against the
cached-FC value; the run aborts if they disagree beyond float32 precision.

================== 5-point critical preamble (pre-registered) ==================

(1) CLAIM. The cross-phase cophenetic trace ``rho_sym(s) > 0`` observed in
    ``alpha`` and ``beta`` at mst@0.20 reflects task-induced reorganization of
    *time-lagged* interaction that persists into ``rest_post``, and is not a
    by-product of the coherence-magnitude structure or of the estimator.

(2) NULLS.
    N1  lag-randomized coherency. Per channel one circular time shift, applied
        as its exact frequency-domain image (DFT shift theorem). ``|C_ij(f)|``
        and all per-channel spectra preserved EXACTLY; only cross-spectral
        phase -- the lag -- is randomized. Independent draws per phase.
    N2  per-(channel, frequency) uniform phase randomization (Theiler et al.
        1992, restricted to the cross-spectrum). Same invariants as N1, with
        the phase decorrelated across frequency bins as well.
    N1b segment-lattice circular shift: each channel's Welch segment sequence is
        circularly rolled. Each channel's Welch PSD is preserved exactly; the
        cross-channel temporal correspondence is destroyed. The spectral form of
        the shift predictor (Perkel, Gerstein & Moore 1967).

(3) STRONGEST ALTERNATIVE THE NULL SHOULD CONTROL FOR. That the trace is an
    artifact of *how much* two contacts cohere -- a volume-conduction /
    field-spread geometry that is stable across the session and, being shared by
    all four phases, manufactures a positive cross-phase Spearman without any
    task-induced change in lagged interaction. Matched-strength cannot address
    this: it acts after ``|ImCoh|`` has been formed, so it inherits whatever the
    estimator produced.

(4) DOES IT CONTROL FOR IT, BY MECHANISM -- AND WHAT IT CANNOT REJECT.
    N1/N2 hold ``|C_ij(f)|`` fixed at its observed, phase-specific value and
    destroy only the phase. Randomizing the phase drives ``<|Im C|>_f`` to
    ``(2/pi)<|C|>_f``: measured x4.6 inflation of the edge weights on Pat_05
    beta. This is harmless ONLY because the readout is invariant to a global
    positive rescale of W -- verified numerically at residual 0.0e+00 (script 02,
    V5): the backbone is rank-based and ``s = tau*lambda_max`` makes ``tau L``
    scale-free. So N1/N2 answer exactly: *does the coherence-magnitude structure
    alone, lag destroyed, reproduce the trace?*
    CANNOT REJECT: N1/N2 keep each phase's own ``|C|`` intact, so they cannot
    rule out an account in which the trace lives in cross-phase changes of
    coherence MAGNITUDE rather than of lag. They also cannot reach session
    nonstationarity, drift, artifact epochs, or the phase segmentation -- all
    four phases keep their true boundaries. That is N3's job.
    N1b additionally destroys coupling magnitude (measured: mean ``|C|`` 0.241 ->
    0.035), so it is an INDEPENDENCE null, not a lag null. It is the weakest rung
    and clearing it proves only that some genuine simultaneous coupling is
    required. It must not be reported as "the circular-shift null" without that
    qualification, because for a Welch-averaged estimator a whole-recording roll
    and a within-segment phase rotation are different operations.

(5) FALSIFICATION. If the cohort gate (one-sided Wilcoxon on
    ``obs - surr_p50`` across the 10 patients) fails at every scale for a band
    under N1, then that band's trace is not attributable to lagged interaction
    and the |ImCoh| framing of the result is wrong. If beta fails N1 the paper's
    central claim fails. Remaining limitation regardless of outcome: n = 10 gives
    a coarse Wilcoxon p-grid, and N1/N2 leave the magnitude channel untested.

Outputs: data/paper_final/w0b_nulls/{rung}/per_patient_scale.csv, cohort_gate.csv
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
    band_bin_frequencies, coherency_from_csd, csd_from_segment_subset,
    imcoh_abs_from_coherency, lag_randomized_coherency,
    phase_randomized_coherency, segment_fft, segment_shifted_csd,
)
from audit_150_rho_sym_gate import COHORT, BANDS, load_phase

PHASES = ("A", "B", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)          # identical to script 13
FRAC, BACKBONE = 0.20, "mst020"
R = 200
BASE_SEED = 20260825
OUT_ROOT = ROOT / "data" / "paper_final" / "w0b_nulls"


def eig_mst(W):
    return laplacian_eig(select_backbone(W, BACKBONE, frac=FRAC))


def readout(W_by_phase, masks=None):
    """rho_sym(s) through the production backbone + LRG path.

    ``masks`` pins the backbone EDGE SET to a supplied one instead of re-selecting
    it from the surrogate weights. This is what the ``n1_fixedbb`` rung needs.
    Rationale: N1 sends ``W`` toward ``(2/pi)|C|``, and the |C| edge ranking is not
    the |ImCoh| ranking -- on sEEG, |C| is dominated by same-probe pairs (a known
    2-8x bias, CLAUDE.md invariant 5). A plain N1 surrogate therefore lives on a
    DIFFERENT graph from the observed, and part of its high floor could be the
    stability of that anatomical scaffold across phases rather than anything about
    lag. Holding the observed backbone fixed removes that explanation, so the
    comparison isolates the lag content of the weights.
    """
    eig = {}
    for ph in PHASES:
        W = W_by_phase[ph]
        if masks is None:
            eig[ph] = eig_mst(W)
        else:
            eig[ph] = laplacian_eig(W * masks[ph])
    return rho_sym_over_scales(eig, SGRID)


def load_segment_ffts(pat, band):
    """Per-segment band FFTs for the four canonical phases.

    Segmentation matches the production FC exactly: task_test / rest_post at
    ``nperseg_for_fs(fs)``; the rest_pre halves A / B at ``nperseg_for_fs//2``,
    which is the setting the cached ``imcoh_halves_fc`` matrices were built with
    (audit_63 pre-flight). Using the production nperseg per phase is what makes
    the recomputed observed statistic reproduce the cached-FC one.
    """
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs)
    nps_half = max(256, nps // 2)
    bnd = BRAIN_BANDS[band]
    out = {}
    Xr = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    if Xr.shape[0] > Xr.shape[1]:
        Xr = Xr.T
    T = Xr.shape[1]
    for tag, sl in (("A", slice(0, T // 2)), ("B", slice(T // 2, T))):
        out[tag] = segment_fft(Xr[:, sl], fs, nps_half, band=bnd) + (fs, nps_half)
    del Xr
    for ph in ("task_test", "rest_post"):
        X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
        if X.shape[0] > X.shape[1]:
            X = X.T
        out[ph] = segment_fft(X, fs, nps, band=bnd) + (fs, nps)
        del X
    return out


def per_cell(job):
    idx, pat, band, rung = job
    t0 = time.time()
    try:
        segs = load_segment_ffts(pat, band)
    except Exception as e:
        return None, f"{pat}/{band}: load failed ({e})"

    # ---- observed, recomputed through the surrogate code path ----
    Cobs, Wobs = {}, {}
    for ph, (freqs, F, scale, fs, nps) in segs.items():
        S = csd_from_segment_subset(F, None, scale)
        Cobs[ph] = coherency_from_csd(S)
        Wobs[ph] = imcoh_abs_from_coherency(Cobs[ph])
    # n1_fixedbb pins the backbone to the OBSERVED one (see readout docstring)
    masks = None
    if rung == "n1_fixedbb":
        masks = {ph: (select_backbone(Wobs[ph], BACKBONE, frac=FRAC) > 0).astype(float)
                 for ph in PHASES}
    obs = readout(Wobs, masks)

    # ---- equivalence gate vs the cached-FC observed (script 13's statistic) ----
    try:
        Wc = {ph: load_phase(pat, ph, band) for ph in PHASES}
        obs_cached = readout(Wc, masks)     # same footing as obs for every rung
        d_obs = float(np.nanmax(np.abs(obs - obs_cached)))
    except Exception:
        obs_cached, d_obs = np.full(SGRID.size, np.nan), np.nan

    rng = np.random.default_rng([BASE_SEED, idx])
    surr = np.full((R, SGRID.size), np.nan)
    for r in range(R):
        Ws = {}
        for ph, (freqs, F, scale, fs, nps) in segs.items():
            if rung in ("n1", "n1_fixedbb"):
                Cs = lag_randomized_coherency(Cobs[ph], freqs, fs, rng, nperseg=nps)
                Ws[ph] = imcoh_abs_from_coherency(Cs)
            elif rung == "n2":
                Cs = phase_randomized_coherency(Cobs[ph], rng)
                Ws[ph] = imcoh_abs_from_coherency(Cs)
            elif rung == "n1b":
                S = segment_shifted_csd(F, rng, scale, min_shift=1)
                Ws[ph] = imcoh_abs_from_coherency(coherency_from_csd(S))
            else:
                raise ValueError(rung)
        surr[r] = readout(Ws, masks)

    rows = []
    for j, s in enumerate(SGRID):
        col = surr[:, j][np.isfinite(surr[:, j])]
        o = obs[j]
        rows.append(dict(
            patient=pat, band=band, rung=rung, s=float(s), N=int(Wobs["A"].shape[0]),
            obs_rho=float(o) if np.isfinite(o) else np.nan,
            obs_rho_cached=float(obs_cached[j]) if np.isfinite(obs_cached[j]) else np.nan,
            surr_p50=float(np.percentile(col, 50)) if col.size else np.nan,
            surr_p95=float(np.percentile(col, 95)) if col.size else np.nan,
            surr_mean=float(col.mean()) if col.size else np.nan,
            p=float(np.mean(col >= o)) if col.size and np.isfinite(o) else np.nan,
            n_surr=int(col.size)))
    return rows, (f"{pat}/{band}/{rung} obs_dev={d_obs:.2e} "
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
    ap.add_argument("--rungs", default="n1,n2,n1b",
                    help="n1 | n2 | n1b | n1_fixedbb (N1 on the observed backbone)")
    ap.add_argument("--bands", default="")
    ap.add_argument("--patients", default="")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("NL_WORKERS", 10)))
    a = ap.parse_args()

    rungs = a.rungs.split(",")
    bands = a.bands.split(",") if a.bands else BANDS
    pats = a.patients.split(",") if a.patients else COHORT
    jobs = [(i, p, b, rg) for i, (rg, b, p) in
            enumerate((rg, b, p) for rg in rungs for b in bands for p in pats)]
    print(f"[w0b] {len(jobs)} cells | rungs={rungs} bands={bands} n_pat={len(pats)} "
          f"R={R} scales={len(SGRID)} workers={a.workers}", flush=True)

    parts = OUT_ROOT / "_parts_rho"             # per-cell checkpoint; resumable
    parts.mkdir(parents=True, exist_ok=True)
    done = {p.stem for p in parts.glob("*.csv")}
    todo = [j for j in jobs if f"{j[1]}__{j[2]}__{j[3]}" not in done]
    print(f"[w0b] {len(done)} cells checkpointed; {len(todo)} to run", flush=True)

    t0 = time.time()
    with Pool(a.workers) as pool:
        for i, (rl, msg) in enumerate(pool.imap_unordered(per_cell, todo), 1):
            if rl:
                r0 = rl[0]
                pd.DataFrame(rl).to_csv(
                    parts / f"{r0['patient']}__{r0['band']}__{r0['rung']}.csv", index=False)
            el = time.time() - t0
            print(f"[{i}/{len(todo)}] {msg} | {el:.0f}s ETA {el/i*(len(todo)-i):.0f}s", flush=True)

    got = [pd.read_csv(p) for p in sorted(parts.glob("*.csv"))]
    df = pd.concat(got, ignore_index=True) if got else pd.DataFrame()
    if df.empty:
        print("[w0b] no rows"); return
    for rung in rungs:
        d = df[df.rung == rung]
        if d.empty:
            continue
        od = OUT_ROOT / rung
        od.mkdir(parents=True, exist_ok=True)
        # MERGE, never clobber: a band-subset re-run must extend the table, not
        # replace it. New rows win on (patient, band, rung, s).
        pp = od / "per_patient_scale.csv"
        if pp.exists():
            prev = pd.read_csv(pp)
            d = (pd.concat([prev, d], ignore_index=True)
                   .drop_duplicates(subset=["patient", "band", "rung", "s"], keep="last")
                   .sort_values(["band", "patient", "s"]))
        d.to_csv(pp, index=False)
        g = cohort_gate(d)
        g.to_csv(od / "cohort_gate.csv", index=False)
        (od / "config.json").write_text(json.dumps(dict(
            rung=rung, backbone=BACKBONE, frac=FRAC, R=R,
            s_grid=[float(x) for x in SGRID], cohort=pats, bands=bands,
            base_seed=BASE_SEED,
            note="timeseries/coherency-level null; readout identical to script 13",
        ), indent=2))
        print(f"\n=== {rung}: cohort gate, per band (min gate_p over scales) ===", flush=True)
        if g.empty or "band" not in g.columns:
            print("  (cohort gate needs >=5 patients; subset/timing run)", flush=True)
            continue
        for band in bands:
            gb = g[g.band == band].dropna(subset=["gate_p"])
            if gb.empty:
                continue
            b = gb.loc[gb.gate_p.idxmin()]
            n_clear = int((gb.gate_p < 0.05).sum())
            print(f"  {band:11s} best s={b.s:6.1f} gate_p={b.gate_p:.4f} "
                  f"clears {n_clear}/{len(gb)} scales  obs_med={b.obs_med:+.3f} "
                  f"surr_med={b.surr_med:+.3f}", flush=True)
    print(f"\n[w0b] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT_ROOT}", flush=True)


if __name__ == "__main__":
    main()
