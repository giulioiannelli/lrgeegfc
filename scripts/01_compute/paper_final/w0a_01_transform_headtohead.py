#!/usr/bin/env python3
"""W0-A / A1 — |ImCoh| vs |ImCoh|^2 head to head on hypothesis-independent properties.

PRE-REGISTERED DECISION RULE R1 (frozen in w0a_00_preregistration.md BEFORE any
number was computed; reproduced verbatim so the script is self-contained):

  1. Rank-equivalence gate. If the two transforms are rank-equivalent on the
     upper-triangle edge vector (Spearman >= 0.9999 AND top-f edge-set Jaccard
     == 1 at f in {0.05, 0.10, 0.20}), the choice is downstream-vacuous for any
     rank-based backbone: keep the incumbent imcoh_abs and REPORT THE INVARIANCE
     as the justification. Rule stops.
  2. Split-half reliability. Otherwise prefer the transform with higher cohort
     SPEARMAN reliability between the two rest_pre halves' edge vectors
     (test-retest of the estimator itself). Decisive only if it wins in >= 4/6
     bands AND a Wilcoxon over the 60 (patient, band) cells is p < 0.05.
  3. Separation from the lag-destroyed floor (tie-break).
  4. Lower weight heterogeneity (final tie-break).

  Trace gates / marker AUCs under each transform are RECORDED but are NOT
  admissible at any step of R1.

5-POINT CRITICAL PREAMBLE
1. Claim. The band transform can be fixed without reference to the task-trace
   hypothesis, using (a) the estimator's own split-half test-retest and (b) the
   distance of the observed weights from a null that destroys lag but preserves
   coherence magnitude exactly.
2. Null. H0_A: the two transforms are rank-equivalent, so the choice cannot
   matter for any rank-based backbone and the 2026-04-15 debate is vacuous.
   H0_B (per edge): the observed |ImCoh| is what a lag-free process with the
   SAME coherence magnitude would produce -- i.e. the edge weight is volume
   conduction plus estimator bias, not lagged interaction.
3. Strongest plausible alternative. That any measured reliability difference is
   a units artifact: squaring is a monotone map on non-negative numbers, so a
   Pearson correlation changes under it while nothing about the estimator does.
4. Does the null control for it, by mechanism? Yes for the units artifact: the
   PRIMARY reliability statistic is Spearman, which is invariant under any
   strictly monotone map, so a Spearman difference can ONLY come from the fact
   that abs and sq are not related by a monotone map -- the transform is applied
   per frequency bin and THEN band-averaged, so mean(|x|) and mean(x^2) reorder
   edges (Jensen). Pearson is reported as secondary and is NOT used to decide.
   The circular-shift floor controls H0_B by mechanism: a circular time shift of
   channel c by t_c maps S_ij(f) -> S_ij(f) e^{-2*pi*i*f*(t_i - t_j)}, so
   |C_ij(f)| -- hence every bit of volume conduction and of the coupling
   MAGNITUDE -- is preserved exactly and only the phase (the lag) is randomized.
   What it CANNOT reject: it does not test the band split, the Welch settings,
   session nonstationarity, or whether the lag is neural rather than a
   referencing / filtering artifact; and the frequency-domain phase-ramp form
   used here is the idealized limit of a circular shift (Welch segmenting makes
   an exact time-domain circular shift only approximately a pure phase ramp).
5. Falsification / limits. R1 step 2 is falsified for a transform if it loses
   the reliability comparison; the whole approach is falsified if the two
   transforms turn out rank-equivalent (step 1 fires) -- in which case the
   correct report is "the choice does not matter", not a justification for one.
   Split-half reliability is a WITHIN-session test-retest: it cannot separate
   estimator noise from genuine within-session nonstationarity, and it is
   computed at half the Welch segment length, so it understates the reliability
   of the full-duration matrix.

Outputs (data/paper_final/w0a_substrate/a1_transform/):
  per_cell.csv        patient, band, transform: reliability, weight-shape,
                      lag-concentration vs circular-shift floor, backbone agreement
  summary.csv         cohort medians per (band, transform) + the R1 verdict inputs
  config.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew, spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    DEFAULT_SAMPLE_RATE,
    FS_OVERRIDES,
    PATIENTS_4PHASE,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.fc.backbone import backbone_density, mst_union_top_fraction
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.metrics.graph_descriptors import gini_coefficient
from lrg_eegfc.utils.surrogate.coherency_surrogate import complex_coherency_bands

COHORT = list(PATIENTS_4PHASE)
BANDS = list(BRAIN_BANDS)
TRANSFORMS = ("abs", "sq")
TOP_FRACS = (0.05, 0.10, 0.20)
N_SHIFT_DRAWS = 32
BASE_SEED = 20260825
OUT = ROOT / "data" / "paper_final" / "w0a_substrate" / "a1_transform"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"      # reference imcoh_abs halves (validation)

# Random-phase expectation of the lag concentration (analytic):
#   abs: E|sin(U)| = 2/pi ;  sq: E[sin^2(U)] = 1/2
RANDOM_PHASE_L = {"abs": 2.0 / np.pi, "sq": 0.5}


# --------------------------------------------------------------------------- #
# transforms + shape statistics
# --------------------------------------------------------------------------- #
def band_transform(C: np.ndarray, kind: str) -> np.ndarray:
    """(F, N, N) complex coherency -> (N, N) non-negative adjacency.

    ``abs`` = mean_f |Im C| (Ewald 2012 / Bastos-Schoffelen 2016);
    ``sq``  = mean_f (Im C)^2 (Ewald 2012). Transform per bin, THEN band-average
    -- the order the canonical loader uses.
    """
    im = C.imag
    A = np.abs(im).mean(0) if kind == "abs" else (im ** 2).mean(0)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


def coherence_magnitude(C: np.ndarray, kind: str) -> np.ndarray:
    """Band-averaged coherence MAGNITUDE in the same units as the transform.

    ``abs`` -> mean_f |C| ; ``sq`` -> mean_f |C|^2. This is the quantity a
    circular shift preserves exactly, hence the denominator of the lag
    concentration.
    """
    m = np.abs(C)
    A = m.mean(0) if kind == "abs" else (m ** 2).mean(0)
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    return A


def circular_shift_floor(C: np.ndarray, freqs: np.ndarray, duration_s: float,
                         rng: np.random.Generator, n_draw: int) -> dict:
    """Empirical circular-shift (lag-destroying) floor of both transforms.

    Shifting channel ``c`` by ``t_c`` seconds maps ``C_ij(f) -> C_ij(f) *
    exp(-2*pi*i*f*(t_i - t_j))``: ``|C_ij(f)|`` is preserved exactly, only the
    phase is randomized. Returns per-transform mean and sd (over draws) of the
    resulting adjacency, upper triangle only.
    """
    N = C.shape[1]
    r, c = np.triu_indices(N, 1)
    acc = {k: [] for k in TRANSFORMS}
    for _ in range(n_draw):
        t = rng.uniform(0.0, duration_s, N)
        ph = np.exp(-2j * np.pi * np.outer(freqs, t))            # (F, N)
        Ct = C * ph[:, :, None] * np.conj(ph)[:, None, :]
        im = Ct.imag
        acc["abs"].append(np.abs(im).mean(0)[r, c])
        acc["sq"].append((im ** 2).mean(0)[r, c])
    return {k: (np.mean(v, 0), np.std(v, 0)) for k, v in acc.items()}


def shape_stats(w: np.ndarray) -> dict:
    """Hypothesis-independent shape descriptors of an edge-weight vector."""
    w = np.asarray(w, float)
    med = float(np.median(w))
    return dict(
        w_mean=float(w.mean()), w_median=med, w_max=float(w.max()),
        cv=float(w.std() / w.mean()) if w.mean() > 0 else np.nan,
        gini=float(gini_coefficient(w)),
        skew=float(skew(w)), kurtosis=float(kurtosis(w)),
        dyn_range_p99_p50=float(np.percentile(w, 99) / med) if med > 0 else np.nan,
        dyn_range_max_p50=float(w.max() / med) if med > 0 else np.nan,
        # participation ratio: effective number of edges carrying the weight
        eff_n_edges=float(w.sum() ** 2 / (w ** 2).sum()) / w.size,
    )


# --------------------------------------------------------------------------- #
# per-patient worker
# --------------------------------------------------------------------------- #
def per_patient(job):
    idx, pat = job
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nper_full = nperseg_for_fs(fs)
    nper_half = max(256, nper_full // 2)
    try:
        X = load_timeseries(pat, "rest_pre", SEEG_DATAPATH)
    except Exception as exc:                                     # pragma: no cover
        print(f"  [{pat}] load failed: {exc}", flush=True)
        return None
    X = np.asarray(X, float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    N, T = X.shape

    rng = np.random.default_rng(BASE_SEED + idx)
    rows = []

    # --- halves (reliability) at the halved Welch segment ------------------- #
    Chalf = {}
    for tag, sl in (("A", slice(0, T // 2)), ("B", slice(T // 2, T))):
        Chalf[tag] = complex_coherency_bands(
            np.ascontiguousarray(X[:, sl]), fs, BRAIN_BANDS, nper_half
        )
    # --- full phase (canonical settings): floor + shape --------------------- #
    Cfull = complex_coherency_bands(np.ascontiguousarray(X), fs, BRAIN_BANDS, nper_full)
    del X

    r, c = np.triu_indices(N, 1)
    for band in BANDS:
        Cb = Cfull[band]
        # exact in-band frequency grid (Welch bins are k*fs/nperseg)
        allf = np.arange(nper_full // 2 + 1) * (fs / nper_full)
        m = (allf >= BRAIN_BANDS[band][0]) & (allf <= BRAIN_BANDS[band][1])
        freqs_b = allf[m]
        assert freqs_b.size == Cb.shape[0], (band, freqs_b.size, Cb.shape)

        floor = circular_shift_floor(Cb, freqs_b, T / fs, rng, N_SHIFT_DRAWS)

        obs_by_t, mag_by_t = {}, {}
        for tk in TRANSFORMS:
            W_obs = band_transform(Cb, tk)
            W_mag = coherence_magnitude(Cb, tk)
            obs_by_t[tk] = W_obs
            mag_by_t[tk] = W_mag
            o = W_obs[r, c]
            g = W_mag[r, c]
            f_mu, f_sd = floor[tk]

            # lag concentration: L = obs / magnitude, standardized against its
            # own random-phase value so abs and sq are on the same scale.
            L = np.divide(o, g, out=np.zeros_like(o), where=g > 0)
            L0 = RANDOM_PHASE_L[tk]
            lag_excess = (L - L0) / (1.0 - L0)

            with np.errstate(divide="ignore", invalid="ignore"):
                z = np.where(f_sd > 0, (o - f_mu) / f_sd, np.nan)

            # halves reliability (SAME transform, both halves)
            wa = band_transform(Chalf["A"][band], tk)[r, c]
            wb = band_transform(Chalf["B"][band], tk)[r, c]
            rel_s = float(spearmanr(wa, wb).statistic)
            rel_p = float(np.corrcoef(wa, wb)[0, 1])

            row = dict(patient=pat, band=band, transform=tk, N=int(N),
                       n_edges=int(o.size), duration_s=float(T / fs),
                       rel_spearman=rel_s, rel_pearson=rel_p,
                       floor_med=float(np.median(f_mu)),
                       obs_over_floor_med=float(np.median(o) / np.median(f_mu)),
                       lag_excess_med=float(np.median(lag_excess)),
                       lag_excess_p90=float(np.percentile(lag_excess, 90)),
                       frac_edges_above_floor=float(np.mean(o > f_mu)),
                       z_vs_floor_med=float(np.nanmedian(z)),
                       spearman_obs_vs_floor=float(spearmanr(o, f_mu).statistic),
                       spearman_obs_vs_magnitude=float(spearmanr(o, g).statistic))
            row.update(shape_stats(o))
            rows.append(row)

        # --- cross-transform agreement (the rank-equivalence gate) ---------- #
        oa = obs_by_t["abs"][r, c]
        os_ = obs_by_t["sq"][r, c]
        agree = dict(patient=pat, band=band, transform="_cross", N=int(N),
                     spearman_abs_vs_sq=float(spearmanr(oa, os_).statistic),
                     kendall_rank_changes=float(np.mean(
                         np.argsort(np.argsort(oa)) != np.argsort(np.argsort(os_)))))
        for f in TOP_FRACS:
            k = int(np.ceil(f * oa.size))
            ta = set(np.argsort(-oa)[:k].tolist())
            tb = set(np.argsort(-os_)[:k].tolist())
            agree[f"topjac_{f:.2f}"] = len(ta & tb) / len(ta | tb)
            # backbone-level agreement (mst-union, the incumbent family)
            Ba = mst_union_top_fraction(obs_by_t["abs"], f) > 0
            Bb = mst_union_top_fraction(obs_by_t["sq"], f) > 0
            ea, eb = Ba[r, c], Bb[r, c]
            inter = float((ea & eb).sum()); union = float((ea | eb).sum())
            agree[f"bbjac_{f:.2f}"] = inter / union if union else np.nan
            agree[f"bbdens_abs_{f:.2f}"] = backbone_density(
                mst_union_top_fraction(obs_by_t["abs"], f))
        rows.append(agree)

    # --- validation against the existing imcoh_abs halves cache ------------- #
    ref = HALVES_FC_CACHE / pat / f"beta_rest_pre_A_imcoh_abs.npy"
    if ref.exists():
        Wref = np.asarray(np.load(ref), float)
        Wnew = band_transform(Chalf["A"]["beta"], "abs")
        if Wref.shape == Wnew.shape:
            d = float(np.abs(Wref - Wnew).max())
            print(f"  [{pat}] halves-cache reproduction check (beta/A): max|diff| = {d:.3e}",
                  flush=True)
    return rows


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-only", action="store_true",
                    help="re-emit the R1 verdict from the cached per_cell.csv (no recompute)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.report_only:
        df = pd.read_csv(OUT / "per_cell.csv")
        print(f"[a1] report-only from {OUT/'per_cell.csv'} ({len(df)} rows)\n", flush=True)
    else:
        ncpu = int(os.environ.get("SA_WORKERS", 5))
        jobs = [(i, p) for i, p in enumerate(COHORT)]
        print(f"[a1] {len(jobs)} patients x {len(BANDS)} bands x 2 transforms, "
              f"{N_SHIFT_DRAWS} circular-shift draws, {ncpu} workers -> {OUT}", flush=True)
        t0 = time.time()
        rows = []
        with Pool(ncpu) as pool:
            for i, rl in enumerate(pool.imap_unordered(per_patient, jobs), 1):
                if rl:
                    rows.extend(rl)
                el = time.time() - t0
                print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el / i * (len(jobs) - i):.0f}s", flush=True)
        df = pd.DataFrame(rows)
        df.to_csv(OUT / "per_cell.csv", index=False)
        print(f"\n[a1] {len(df)} rows in {time.time() - t0:.0f}s -> {OUT}\n", flush=True)

    # ---------------- R1 step 1: rank-equivalence gate --------------------- #
    x = df[df["transform"] == "_cross"]
    print("=== R1 step 1 — rank-equivalence gate (abs vs sq) ===", flush=True)
    gate_rows = []
    for band in BANDS:
        b = x[x.band == band]
        if b.empty:
            continue
        sp = b.spearman_abs_vs_sq.median()
        jac = {f: b[f"topjac_{f:.2f}"].median() for f in TOP_FRACS}
        bb = {f: b[f"bbjac_{f:.2f}"].median() for f in TOP_FRACS}
        print(f"  {band:11s} spearman={sp:.5f}  top-jac "
              + "  ".join(f"{f:.2f}:{jac[f]:.3f}" for f in TOP_FRACS)
              + "  |  mst-union edge-set jac "
              + "  ".join(f"{f:.2f}:{bb[f]:.3f}" for f in TOP_FRACS), flush=True)
        gate_rows.append(dict(band=band, spearman_abs_vs_sq=sp,
                              **{f"topjac_{f:.2f}": jac[f] for f in TOP_FRACS},
                              **{f"bbjac_{f:.2f}": bb[f] for f in TOP_FRACS}))
    equivalent = (x.spearman_abs_vs_sq.min() >= 0.9999
                  and all((x[f"topjac_{f:.2f}"] >= 1.0 - 1e-12).all() for f in TOP_FRACS))
    print(f"  --> rank-equivalent: {equivalent}"
          f"  (gate fires and stops the rule only if True)", flush=True)

    # ---------------- R1 step 2: split-half reliability -------------------- #
    y = df[df["transform"].isin(TRANSFORMS)]
    piv = y.pivot_table(index=["patient", "band"], columns="transform",
                        values="rel_spearman")
    print("\n=== R1 step 2 — split-half reliability (Spearman, primary) ===", flush=True)
    wins = {"abs": 0, "sq": 0}
    for band in BANDS:
        b = y[y.band == band]
        ma = b[b["transform"] == "abs"].rel_spearman.median()
        ms = b[b["transform"] == "sq"].rel_spearman.median()
        pa = b[b["transform"] == "abs"].rel_pearson.median()
        ps = b[b["transform"] == "sq"].rel_pearson.median()
        w = "abs" if ma > ms else "sq"
        wins[w] += 1
        print(f"  {band:11s} spearman abs={ma:.4f} sq={ms:.4f}  (delta={ma-ms:+.4f}, win={w})"
              f"   [pearson abs={pa:.4f} sq={ps:.4f}]", flush=True)
    d = (piv["abs"] - piv["sq"]).dropna()
    try:
        p_two = float(wilcoxon(d.values)[1])
    except Exception:
        p_two = np.nan
    print(f"  cohort (n={len(d)} cells) median delta(abs-sq) = {np.median(d):+.5f}, "
          f"Wilcoxon two-sided p = {p_two:.4g};  band wins abs={wins['abs']} sq={wins['sq']}",
          flush=True)
    decisive2 = (max(wins.values()) >= 4) and (p_two < 0.05)
    print(f"  --> step 2 decisive: {decisive2}"
          f" (needs >=4/6 band wins AND p<0.05)", flush=True)

    # ---------------- R1 step 3: separation from the floor ----------------- #
    print("\n=== R1 step 3 — separation from the circular-shift (lag-destroyed) floor ===",
          flush=True)
    for band in BANDS:
        b = y[y.band == band]
        for tk in TRANSFORMS:
            t = b[b["transform"] == tk]
            print(f"  {band:11s} {tk:3s}  lag_excess_med={t.lag_excess_med.median():+.4f}  "
                  f"obs/floor={t.obs_over_floor_med.median():.4f}  "
                  f"frac_above={t.frac_edges_above_floor.median():.3f}  "
                  f"z_vs_floor={t.z_vs_floor_med.median():+.2f}  "
                  f"sp(obs,floor)={t.spearman_obs_vs_floor.median():.3f}", flush=True)

    # ---------------- R1 step 4: weight heterogeneity ---------------------- #
    print("\n=== R1 step 4 — weight heterogeneity (lower = further from the dense-degenerate regime) ===",
          flush=True)
    for band in BANDS:
        b = y[y.band == band]
        for tk in TRANSFORMS:
            t = b[b["transform"] == tk]
            print(f"  {band:11s} {tk:3s}  gini={t['gini'].median():.4f}  cv={t['cv'].median():.4f}  "
                  f"skew={t['skew'].median():+.3f}  p99/p50={t['dyn_range_p99_p50'].median():.2f}  "
                  f"eff_edge_frac={t['eff_n_edges'].median():.4f}", flush=True)

    # ------- POST-HOC DIAGNOSTIC (not part of the frozen rule) -------------- #
    # The strongest alternative reading of a reliability win: it is inherited
    # from coherence MAGNITUDE (anatomy / volume conduction), which is highly
    # reproducible across halves and has nothing to do with lag estimation. If
    # the winner's edge ranking is more magnitude-driven, its reliability
    # advantage is anatomy, not estimator quality. This diagnostic can only
    # WEAKEN the step-2 conclusion; it cannot create one.
    print("\n=== post-hoc diagnostic — is the reliability win inherited from coherence magnitude? ===",
          flush=True)
    for band in BANDS:
        b = y[y.band == band]
        va = b[b["transform"] == "abs"].spearman_obs_vs_magnitude.median()
        vs = b[b["transform"] == "sq"].spearman_obs_vs_magnitude.median()
        print(f"  {band:11s} spearman(obs, band-avg |C|): abs={va:+.4f}  sq={vs:+.4f}  "
              f"(delta={va-vs:+.4f}; positive => abs more magnitude-driven)", flush=True)

    summ = y.groupby(["band", "transform"]).median(numeric_only=True).reset_index()
    summ.to_csv(OUT / "summary.csv", index=False)
    pd.DataFrame(gate_rows).to_csv(OUT / "rank_equivalence_gate.csv", index=False)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="a1_transform_headtohead", cohort=COHORT, bands=BANDS,
        transforms=list(TRANSFORMS), n_shift_draws=N_SHIFT_DRAWS,
        base_seed=BASE_SEED, top_fracs=list(TOP_FRACS),
        rule="w0a_00_preregistration.md R1 (frozen before any number)",
        floor="frequency-domain circular shift: C_ij(f) -> C_ij(f) e^{-2 pi i f (t_i - t_j)}, "
              "t_c ~ U(0, T); preserves |C_ij(f)| exactly, randomizes lag only",
    ), indent=2))
    print(f"\n[a1] done -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
