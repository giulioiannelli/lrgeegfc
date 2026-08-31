#!/usr/bin/env python3
"""W0-B ★ NULL CALIBRATION of the five-phase cross-phase functionals.

THE POINT. A p-value is only meaningful if the statistic it gates returns ~0 when
fed an input that contains no signal. For a PARTIAL correlation that is not
automatic: the conditioning is entangled with the estimator, so a construction can
manufacture a positive conditional value out of nothing. This is documented, not
hypothetical -- ``.agents/preprint/supplementary/S2_drift_controls.md`` records a
windowed sham arc built entirely from pre-task ``rest_pre`` returning
``T_infspec_pe = +0.243`` for beta, ABOVE the real ``+0.091``, with sham-above-zero
``p = 0.007``. The conclusion drawn at the time ("inference is drift-confounded")
was withdrawn; the correct conclusion was that the NULL was invalid.

So: before any rung of the ladder is allowed to produce a p-value for
``T_infspec_pe`` (or ``T_learn``), this script measures what each functional returns
on a genuinely no-signal input, and by which mechanism any offset arises.

CONSTRUCTIONS (all five pseudo-phases drawn from ONE resting recording, where no
task and therefore no consolidation signal can exist):

  sham_ordered   five contiguous windows in true temporal order, sized in
                 proportion to the real phase durations. This reproduces the
                 audit_168 construction that broke, and RETAINS any within-session
                 drift.
  sham_shuffled  the same five window sizes, but assembled from randomly
                 reassigned 30 s blocks. Destroys temporal order, hence destroys
                 drift, while keeping every other property of the construction
                 (durations, spectral degrees of freedom, the estimator).

Contrasting the two SEPARATES the two candidate mechanisms, which the earlier work
could not do:
  * positive in ordered AND ~0 in shuffled  -> drift manufactures the offset;
  * positive in BOTH                        -> the estimator/construction itself
                                               manufactures it, and no temporal
                                               null can rescue the functional;
  * ~0 in both                              -> the functional is calibrated and
                                               its p-values are interpretable.

Run on ``rest_pre`` and, independently, on ``rest_post`` -- two disjoint no-signal
recordings, so a single-recording quirk cannot masquerade as a general property.

The REAL arc is recomputed here too, through the identical code path, so the
calibration offset can be compared against the effect it is supposed to gate.

================== 5-point critical preamble (pre-registered) ==================
(1) CLAIM under test: that ``T_infspec_pe`` (and ``T_learn``, ``T_infspec``,
    ``T_test``) measure a property of the task arc rather than of the way the arc
    is constructed.
(2) NULL: a five-phase arc carved out of a single resting recording. No encoding,
    no inference, no consolidation is possible inside it, so every functional must
    return ~0.
(3) STRONGEST ALTERNATIVE the check must control for: that a positive sham value
    is caused by session drift rather than by the estimator -- which would leave
    the functional itself sound and merely demand a drift-aware null.
(4) DOES IT CONTROL FOR IT, BY MECHANISM: yes, by the ordered/shuffled contrast
    above, which holds the estimator, the durations and the degrees of freedom
    fixed and toggles ONLY the presence of temporal order.
    CANNOT REJECT: a resting recording is not a perfect no-signal control -- it
    contains its own slow structure, and if the estimator responds to ANY temporal
    structure the shuffled arm removes it rather than proving its absence in the
    real arc. Nor does this settle whether the real value is large ENOUGH; it only
    settles whether zero is the right reference point.
(5) FALSIFICATION: if a functional is significantly positive on the shuffled sham
    (cohort Wilcoxon, two-sided, against 0), it is not zero-centred, its p-values
    against any rung are uninterpretable as stated, and it must be reported as
    uncalibrated rather than gated.

Outputs: data/paper_final/w0b_nulls/calibration/{per_cell.csv, cohort.csv}
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
from lrg_eegfc.utils.fc.heat_multiscale import (
    CROSS_PHASE_FUNCTIONALS, cross_phase_functionals_over_scales, laplacian_eig,
)
from lrg_eegfc.utils.surrogate.timeseries_nulls import (
    coherency_from_csd, csd_from_segment_subset, imcoh_abs_from_coherency, segment_fft,
)
from audit_150_rho_sym_gate import COHORT, BANDS

PHASES5 = ("A", "B", "task_learn", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)
FRAC, BACKBONE = 0.20, "mst020"
BLOCK_SECONDS = 30.0
N_SHUF = 40                      # shuffled-sham realizations per cell
BASE_SEED = 20260826
OUT = ROOT / "data" / "paper_final" / "w0b_nulls" / "calibration"
_DP = ROOT / "data" / "paper_final" / "w0b_nulls" / "durations.json"
if not _DP.exists():
    raise SystemExit(f"missing {_DP}; run 09_cache_durations.py first")
DURS = json.loads(_DP.read_text())


def eig(W):
    return laplacian_eig(select_backbone(W, BACKBONE, frac=FRAC))


def W_from_segments(F, idx, scale):
    return imcoh_abs_from_coherency(
        coherency_from_csd(csd_from_segment_subset(F, np.sort(idx), scale)))


def real_arcs(pat, band, keep_fracs):
    """REAL five-phase arcs at several truncation fractions, ONE load per phase.

    ``keep_fracs`` maps a tag to a fraction of each phase's leading segments
    (``None`` = full). Every requested arc is derived from the SAME per-segment
    coefficients, so the recording is read once -- reloading it per arc dominated
    the runtime (9.5 h ETA) for no scientific gain.

    The truncated arc is the DURATION-MATCHED comparator: the sham carves five
    pseudo-phases out of one ~600 s resting recording, so each sham phase holds
    roughly a quarter of the segments of its real counterpart. Comparing the sham
    against the FULL real arc would confound "no task" with "less data".
    """
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs); nps_h = max(256, nps // 2)
    bnd = BRAIN_BANDS[band]
    W = {t: {} for t in keep_fracs}

    def add(ph, F, sc):
        for t, kf in keep_fracs.items():
            n = F.shape[1] if kf is None else max(8, int(F.shape[1] * kf))
            W[t][ph] = W_from_segments(F, np.arange(n), sc)

    Xr = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    if Xr.shape[0] > Xr.shape[1]:
        Xr = Xr.T
    T = Xr.shape[1]
    for tag, sl in (("A", slice(0, T // 2)), ("B", slice(T // 2, T))):
        _, F, sc = segment_fft(Xr[:, sl], fs, nps_h, band=bnd)
        add(tag, F, sc)
        del F
    del Xr
    for ph in ("task_learn", "task_test", "rest_post"):
        X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
        if X.shape[0] > X.shape[1]:
            X = X.T
        _, F, sc = segment_fft(X, fs, nps, band=bnd)
        add(ph, F, sc)
        del X, F
    return {t: cross_phase_functionals_over_scales(
        {p: eig(W[t][p]) for p in PHASES5}, SGRID) for t in keep_fracs}


def sham_arcs(pat, band, source, rng):
    """Sham five-phase arcs carved out of ONE resting recording.

    Returns (ordered_result, [shuffled_result, ...]). Window sizes are in
    proportion to the REAL phase durations, so the sham reproduces the real arc's
    duration profile and its spectral degrees of freedom -- the only thing missing
    is the task.
    """
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs)
    bnd = BRAIN_BANDS[band]

    durs = DURS[pat]                       # cached shapes; no reload just for a shape
    frac = np.array([durs["rest_pre"] / 2, durs["rest_pre"] / 2, durs["task_learn"],
                     durs["task_test"], durs["rest_post"]], float)
    frac /= frac.sum()

    Xs = np.asarray(load_timeseries(pat, source, SEEG_DATAPATH), float)
    if Xs.shape[0] > Xs.shape[1]:
        Xs = Xs.T
    _, F, sc = segment_fft(Xs, fs, nps, band=bnd)
    del Xs
    n_seg = F.shape[1]
    step = nps // 2
    blk = int(round(BLOCK_SECONDS * fs))
    seg_start = np.arange(n_seg) * step
    b0, b1 = seg_start // blk, (seg_start + nps - 1) // blk
    bos = np.where(b0 == b1, b0, -1).astype(np.int64)
    nb = int(bos.max()) + 1
    sizes = np.maximum(2, np.floor(frac * nb).astype(int))
    while sizes.sum() > nb:
        sizes[int(np.argmax(sizes))] -= 1

    def score(block_sets):
        W = {}
        for ph, blocks in zip(PHASES5, block_sets):
            idx = np.concatenate([np.flatnonzero(bos == b) for b in blocks])
            if idx.size < 8:
                return None
            W[ph] = W_from_segments(F, idx, sc)
        return cross_phase_functionals_over_scales({p: eig(W[p]) for p in PHASES5}, SGRID)

    k, ordered_sets = 0, []
    for n in sizes:
        ordered_sets.append(np.arange(k, k + n)); k += n
    res_ord = score(ordered_sets)

    res_shuf = []
    for _ in range(N_SHUF):
        perm = rng.permutation(nb)
        k, sets = 0, []
        for n in sizes:
            sets.append(perm[k:k + n]); k += n
        r = score(sets)
        if r is not None:
            res_shuf.append(r)
    del F
    return res_ord, res_shuf, int(nb), sizes.tolist()


def per_cell(job):
    idx, pat, band, source = job
    t0 = time.time()
    rng = np.random.default_rng([BASE_SEED, idx])
    rows = []
    try:
        if source == "REAL":
            # rest_pre duration / whole-session duration = the fraction of each
            # phase the sham can afford; the matched real arc uses the same.
            durs = DURS[pat]
            kf = durs["rest_pre"] / float(durs["rest_pre"] + durs["task_learn"]
                                          + durs["task_test"] + durs["rest_post"])
            arcs = real_arcs(pat, band, {"real": None, "real_durmatched": kf})
            for tag, r in arcs.items():
                for k in CROSS_PHASE_FUNCTIONALS:
                    for j, s in enumerate(SGRID):
                        rows.append(dict(patient=pat, band=band, source="real",
                                         construction=tag, func=k, s=float(s),
                                         value=float(r[k][j])))
            msg = (f"{pat}/{band}/real T_test[s=1]={arcs['real']['T_test'][0]:+.3f} "
                   f"T_ispe[s=1]={arcs['real']['T_infspec_pe'][0]:+.3f} | "
                   f"durmatched(kf={kf:.2f}) "
                   f"T_test={arcs['real_durmatched']['T_test'][0]:+.3f} "
                   f"T_ispe={arcs['real_durmatched']['T_infspec_pe'][0]:+.3f}")
        else:
            ro, rs, nb, sz = sham_arcs(pat, band, source, rng)
            for k in CROSS_PHASE_FUNCTIONALS:
                for j, s in enumerate(SGRID):
                    if ro is not None:
                        rows.append(dict(patient=pat, band=band, source=source,
                                         construction="sham_ordered", func=k,
                                         s=float(s), value=float(ro[k][j])))
                    if rs:
                        v = np.array([x[k][j] for x in rs], float)
                        rows.append(dict(patient=pat, band=band, source=source,
                                         construction="sham_shuffled", func=k,
                                         s=float(s), value=float(np.nanmedian(v))))
            msg = (f"{pat}/{band}/{source} nb={nb} sizes={sz} "
                   f"ord T_ispe[s=1]={ro['T_infspec_pe'][0]:+.3f} "
                   f"shuf={np.nanmedian([x['T_infspec_pe'][0] for x in rs]):+.3f}")
    except Exception as e:
        return None, f"{pat}/{band}/{source}: FAILED ({type(e).__name__}: {e})", None
    return rows, f"{msg} ({time.time()-t0:.0f}s)", f"{pat}__{band}__{source}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", default="delta,theta,alpha,beta,low_gamma")
    ap.add_argument("--sources", default="REAL,rest_pre,rest_post")
    ap.add_argument("--patients", default="")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("NL_WORKERS", 6)))
    a = ap.parse_args()
    bands = a.bands.split(","); srcs = a.sources.split(",")
    pats = a.patients.split(",") if a.patients else COHORT
    jobs = [(i, p, b, sc) for i, (sc, b, p) in
            enumerate((sc, b, p) for sc in srcs for b in bands for p in pats)]
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[calib] {len(jobs)} cells | sources={srcs} bands={bands} "
          f"N_shuf={N_SHUF} workers={a.workers}", flush=True)

    parts = OUT / "_parts"                    # per-cell checkpoint; resumable
    parts.mkdir(parents=True, exist_ok=True)
    done = {p.stem for p in parts.glob("*.csv")}
    todo = [j for j in jobs if f"{j[1]}__{j[2]}__{j[3]}" not in done]
    print(f"[calib] {len(done)} cells checkpointed; {len(todo)} to run", flush=True)

    t0 = time.time()
    with Pool(a.workers) as pool:
        for i, (rl, msg, tag) in enumerate(pool.imap_unordered(per_cell, todo), 1):
            if rl and tag:
                pd.DataFrame(rl).to_csv(parts / f"{tag}.csv", index=False)
            el = time.time() - t0
            print(f"[{i}/{len(todo)}] {msg} | {el:.0f}s ETA {el/i*(len(todo)-i):.0f}s", flush=True)

    got = [pd.read_csv(p) for p in sorted(parts.glob("*.csv"))]
    df = pd.concat(got, ignore_index=True) if got else pd.DataFrame()
    if df.empty:
        print("[calib] no rows"); return
    pp = OUT / "per_cell.csv"
    df.to_csv(pp, index=False)

    # cohort: is each (construction, func, band, scale) different from ZERO?
    out = []
    for (src, con, fn, band), g in df.groupby(["source", "construction", "func", "band"]):
        for s in np.sort(g.s.unique()):
            v = g[np.isclose(g.s, s)].dropna(subset=["value"]).value.values
            if len(v) < 5:
                continue
            try:
                p2 = float(wilcoxon(v, alternative="two-sided")[1]) if np.any(v != 0) else np.nan
                pg = float(wilcoxon(v, alternative="greater")[1]) if np.any(v != 0) else np.nan
            except Exception:
                p2 = pg = np.nan
            out.append(dict(source=src, construction=con, func=fn, band=band, s=float(s),
                            n=len(v), median=float(np.median(v)),
                            iqr=float(np.percentile(v, 75) - np.percentile(v, 25)),
                            p_vs0_two=p2, p_vs0_greater=pg))
    co = pd.DataFrame(out)
    co.to_csv(OUT / "cohort.csv", index=False)

    print("\n" + "=" * 104)
    print("NULL CALIBRATION — value on a NO-SIGNAL input. A valid functional must sit at ~0.")
    print("=" * 104)
    if co.empty or "func" not in co.columns:
        print("  (cohort calibration needs >=5 patients; subset/timing run)")
        print(f"\nwrote {OUT/'per_cell.csv'}", flush=True)
        return
    for fn in CROSS_PHASE_FUNCTIONALS:
        print(f"\n--- {fn} ---")
        print(f"  {'band':11s} {'construction':15s} {'source':10s} "
              f"{'med(s=1)':>10s} {'med|over s|':>12s} {'scales p<.05 (>0)':>18s}")
        for band in bands:
            for con, src in (("real", "real"), ("real_durmatched", "real"),
                             ("sham_ordered", "rest_pre"),
                             ("sham_shuffled", "rest_pre"), ("sham_ordered", "rest_post"),
                             ("sham_shuffled", "rest_post")):
                g = co[(co.func == fn) & (co.band == band) & (co.construction == con)
                       & (co.source == src)].sort_values("s")
                if g.empty:
                    continue
                # NB: g["median"], never g.median -- the latter is the DataFrame
                # method and silently shadows the column of that name.
                med = g["median"].values
                print(f"  {band:11s} {con:15s} {src:10s} {med[0]:+10.3f} "
                      f"{np.median(med):+12.3f} "
                      f"{int((g.p_vs0_greater < 0.05).sum()):>13d}/{len(g)}")
    print(f"\nwrote {OUT/'per_cell.csv'} and {OUT/'cohort.csv'}", flush=True)
    print(f"[calib] {len(df)} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
