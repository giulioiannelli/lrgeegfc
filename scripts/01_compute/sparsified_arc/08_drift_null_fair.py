#!/usr/bin/env python3
"""D2-drift-FAIR -- a fair, distribution-based, per-scale, SNR-matched drift null.

Supersedes 07_drift_null.py, whose scale-max gate was UNFAIR (it compared each
patient's best-scale obs to its best-scale drift, and the drift arm was a single
noisy 1/5-length-window realization whose max-over-16-scales is inflated). That
statistic wrongly reported "beta = drift". This rebuild fixes every unfairness.

5-point critical preamble
1. CLAIM: the whole-task cross-phase trace rho_sym(s) exceeds within-rest_pre
   temporal drift, per band, AT THE SCALE WHERE IT EMERGES (beta mesoscale s~17).
2. NULL: rho_sym computed on an arc whose four phases are all windows of a single
   rest_pre recording (no task, no consolidation) -> any rho_sym>0 is drift.
   Built as a DISTRIBUTION: R random ordered 4-window placements i<j<k<l mapped
   (i,j)->A,B (early baseline) and (k,l)->task,post (later) -> drift p50/p95 per s.
3. STRONGEST ALTERNATIVE the null must control: the trace is just monotonic FC
   nonstationarity (drift), which matched-strength does NOT remove.
4. DOES IT CONTROL IT -- and the fairness fixes vs 07:
   (a) PER-SCALE paired gate, never scale-max (the scale dependence IS the result);
   (b) DISTRIBUTION over R placements, not one realization (a proper null);
   (c) SNR-MATCHED: obs is recomputed at the SAME short window length L as drift
       (task_test/rest_post subsampled to L), so obs cannot win on data volume;
   (d) NON-ADJACENT: random ordered placements break the w3/w4 collinearity that
       inflated 07's floor.
   CANNOT reject: genuine slow task-independent state change on the SAME timescale
   as within-rest drift (indistinguishable by construction); does not re-test the
   matched-strength (strength) alternative -- a real trace must beat BOTH nulls.
5. FALSIFY: if obs_matched(s) does not exceed the drift p50 band at the trace's
   emergence scale for a band, that band's "trace" on the backbone is drift.

Optimization: per (patient,band) a window-FC bank is precomputed once
(backbone -> eig -> cophenetic(s) for every window); R=200 drift realizations and
the matched-obs draws then only resample cached cophenetic distances (Spearman),
so the whole null is one Welch pass per window. Reuses the D2 s-grid + full-length
obs (trace_arc/{band}/{pat}.npz) as the reference line.

Outputs:
    data/sparsified_arc/trace_arc/drift_fair/{band}/{patient}.npz
        (s, obs_full, obs_matched, drift_p50, drift_p95)
    data/sparsified_arc/trace_arc/drift_fair/gate_vs_s.csv
    data/sparsified_arc/trace_arc/drift_fair/verdict.csv
    data/sparsified_arc/trace_arc/drift_fair/config.json
"""
from __future__ import annotations
import argparse, itertools, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.fc.backbone import percolation_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale, rho_sym

ARC = ROOT / "data" / "sparsified_arc" / "trace_arc"
OUT = ARC / "drift_fair"

FRAC_L = 0.20        # window length L = FRAC_L * rest_pre  (the SNR unit)
NWIN_BANK = 24       # overlapping length-L windows tiling rest_pre (drift bank)
N_OBS = 8            # length-L windows per late phase (task_test, rest_post)
R_DRIFT = 200        # drift realizations (random ordered 4-tuples from the bank)
N_MATCH = 40         # matched-SNR obs draws (A,B in rest_pre; task,post real)
SEED_BASE = 20260712


def _band_abs(Coh, freqs, flo, fhi):        # verbatim audit_167 / 07
    m = (freqs >= flo) & (freqs <= fhi)
    A = np.abs(Coh[:, :, m]).mean(axis=-1)
    A = 0.5 * (A + A.T); np.fill_diagonal(A, 0.0)
    return np.clip(A, 0.0, 1.0)


def _coph_bank_allbands(X, starts, L, fs, nper, sgrids):
    """cophenetic(s) banks for EVERY band from ONE Welch pass per window.

    ``compute_msc_welch`` returns the full freq-resolved cross-spectrum; band-
    averaging is a cheap frequency mask, so we Welch each window ONCE and derive
    all 6 bands (6x fewer Welch than per-band). Returns ``{band: [D or None]}``
    with each D of shape ``(n_s_band, n_pairs)``. One transient C alive at a time."""
    out = {b: [] for b in BANDS}
    for a in starts:
        Xw = X[:, a:a + L]
        fr, C = compute_msc_welch(Xw, fs, nperseg=min(nper, Xw.shape[1]), metric="imcoh")
        for b in BANDS:
            flo, fhi = BRAIN_BANDS[b]
            W = _band_abs(C, fr, flo, fhi)
            try:
                ev, V = laplacian_eig(percolation_backbone(W)[0])
                out[b].append(np.array([cophenetic_at_scale(ev, V, s) for s in sgrids[b]]))
            except Exception:
                out[b].append(None)
        del C
    return out


def _rho_sym_s(DA, DB, Dt, Dp):
    """rho_sym at every scale from four (n_s, n_pairs) cophenetic banks."""
    n_s = DA.shape[0]
    out = np.full(n_s, np.nan)
    for si in range(n_s):
        try:
            out[si] = rho_sym(DA[si], DB[si], Dt[si], Dp[si])[0]
        except Exception:
            pass
    return out


def _band_cell(band, pat, s_grid, obs_full, Dpre, Dtask, Dpost, rng):
    """Drift distribution + matched-SNR obs for one band from precomputed banks."""
    n_s = s_grid.size
    pre_ok = [i for i, d in enumerate(Dpre) if d is not None]
    task_ok = [i for i, d in enumerate(Dtask) if d is not None]
    post_ok = [i for i, d in enumerate(Dpost) if d is not None]
    if len(pre_ok) < 4 or not task_ok or not post_ok:
        return None
    # drift distribution: random ordered 4-tuples i<j<k<l from the rest_pre bank
    combos = list(itertools.combinations(pre_ok, 4))
    sel = rng.choice(len(combos), size=min(R_DRIFT, len(combos)), replace=False)
    drift = np.full((sel.size, n_s), np.nan)
    for r, ci in enumerate(sel):
        i, j, k, l = combos[ci]
        drift[r] = _rho_sym_s(Dpre[i], Dpre[j], Dpre[k], Dpre[l])
    drift_p50 = np.nanpercentile(drift, 50, axis=0)
    drift_p95 = np.nanpercentile(drift, 95, axis=0)
    # matched-SNR obs: A,B two rest_pre windows; task,post REAL, all length L
    obs_m = np.full((N_MATCH, n_s), np.nan)
    for m in range(N_MATCH):
        a, b = rng.choice(pre_ok, 2, replace=False)
        t = int(rng.choice(task_ok)); p = int(rng.choice(post_ok))
        obs_m[m] = _rho_sym_s(Dpre[a], Dpre[b], Dtask[t], Dpost[p])
    obs_matched = np.nanmean(obs_m, axis=0)
    cell = OUT / band; cell.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell / f"{pat}.npz", s=s_grid, obs_full=obs_full,
                        obs_matched=obs_matched, drift_p50=drift_p50, drift_p95=drift_p95)
    return dict(patient=pat, band=band, L_frac=FRAC_L,
                s=";".join(f"{v:.4f}" for v in s_grid),
                obs_matched=";".join(f"{v:.4f}" for v in obs_matched),
                obs_full=";".join(f"{v:.4f}" for v in obs_full),
                drift_p50=";".join(f"{v:.4f}" for v in drift_p50),
                drift_p95=";".join(f"{v:.4f}" for v in drift_p95))


def per_patient(pat):
    """All 6 bands for one patient; ONE Welch bank per window shared across bands."""
    sgrids, obs_full = {}, {}
    for band in BANDS:
        cache = ARC / band / f"{pat}.npz"
        if not cache.exists():
            return []
        d2 = np.load(cache); sgrids[band] = d2["s"]; obs_full[band] = d2["obs"]
    fs = FS_OVERRIDES.get(pat, 2048.0); nper = nperseg_for_fs(fs)

    Xp = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
    Tp = Xp.shape[1]; L = int(Tp * FRAC_L)
    PRE = _coph_bank_allbands(Xp, np.linspace(0, Tp - L, NWIN_BANK).astype(int),
                              L, fs, nper, sgrids)
    del Xp
    Xt = np.asarray(load_timeseries(pat, "task_test", SEEG_DATAPATH), float)
    TASK = _coph_bank_allbands(Xt, np.linspace(0, Xt.shape[1] - L, N_OBS).astype(int),
                               L, fs, nper, sgrids)
    del Xt
    Xr = np.asarray(load_timeseries(pat, "rest_post", SEEG_DATAPATH), float)
    POST = _coph_bank_allbands(Xr, np.linspace(0, Xr.shape[1] - L, N_OBS).astype(int),
                               L, fs, nper, sgrids)
    del Xr

    rows = []
    for band in BANDS:
        rng = np.random.default_rng(SEED_BASE + COHORT.index(pat) * 13 + BANDS.index(band))
        r = _band_cell(band, pat, sgrids[band], obs_full[band],
                       PRE[band], TASK[band], POST[band], rng)
        if r:
            rows.append(r)
    return rows


def cohort(df):
    vs_s, verdict = [], []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        def M(col):
            return np.array([[float(v) for v in r.split(";")] for r in x[col]])
        S = M("s"); OBS = M("obs_matched"); P50 = M("drift_p50"); P95 = M("drift_p95")
        s_axis = np.nanmedian(S, 0); n_s = s_axis.size
        gate = np.full(n_s, np.nan)
        for j in range(n_s):
            d = OBS[:, j] - P50[:, j]; d = d[np.isfinite(d)]
            if d.size >= 5 and np.any(d != 0):
                try:
                    gate[j] = wilcoxon(d, alternative="greater")[1]
                except Exception:
                    pass
            vs_s.append(dict(band=band, s=float(s_axis[j]), gate_p=float(gate[j]),
                             obs_matched_med=float(np.nanmedian(OBS[:, j])),
                             drift_p50_med=float(np.nanmedian(P50[:, j])),
                             drift_p95_med=float(np.nanmedian(P95[:, j])),
                             n_pat_obs_gt_p50=int(np.nansum(OBS[:, j] > P50[:, j]))))
        fire = np.isfinite(gate) & (gate < 0.05)
        # emergence scale = argmax of cohort obs_matched
        s_star = int(np.nanargmax(np.nanmedian(OBS, 0)))
        verdict.append(dict(band=band, n=len(x),
                            n_scales_beat_drift=int(fire.sum()),
                            frac_beat_drift=round(float(fire.mean()), 3),
                            first_s_beat=float(s_axis[fire][0]) if fire.any() else np.nan,
                            s_emerge=float(s_axis[s_star]),
                            gate_at_emerge=float(gate[s_star]),
                            obs_at_emerge=float(np.nanmedian(OBS[:, s_star])),
                            drift_p50_at_emerge=float(np.nanmedian(P50[:, s_star])),
                            beats_drift=bool(fire.any())))
    return pd.DataFrame(vs_s), pd.DataFrame(verdict)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    pats = COHORT[:a.limit] if a.limit else COHORT
    ncpu = int(os.environ.get("SA_WORKERS", 4))
    print(f"[D2-drift-FAIR] {len(pats)} patients x {len(BANDS)} bands | dist R={R_DRIFT}, "
          f"bank={NWIN_BANK}, L={FRAC_L}*rest_pre, matched-SNR obs (1 Welch bank/patient) | "
          f"{ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, rlist in enumerate(pool.imap_unordered(per_patient, pats), 1):
            rows.extend(rlist)
            el = time.time() - t0
            print(f"[{i}/{len(pats)} pt] {el:.0f}s  ETA {el/i*(len(pats)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        print("[D2-drift-FAIR] no cells", flush=True); return
    vs_s, verdict = cohort(df)
    vs_s.to_csv(OUT / "gate_vs_s.csv", index=False)
    verdict.to_csv(OUT / "verdict.csv", index=False)
    print(f"\n[D2-drift-FAIR] {len(rows)} cells in {time.time()-t0:.0f}s\n", flush=True)
    print("=== FAIR drift null: does obs_matched(s) beat the drift p50 band per scale? ===", flush=True)
    print(verdict[["band", "n_scales_beat_drift", "frac_beat_drift", "first_s_beat",
                   "s_emerge", "gate_at_emerge", "obs_at_emerge", "drift_p50_at_emerge",
                   "beats_drift"]].to_string(index=False), flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="D2_drift_null_FAIR", supersedes="07_drift_null.py (scale-max unfair)",
        frac_L=FRAC_L, n_win_bank=NWIN_BANK, n_obs=N_OBS, r_drift=R_DRIFT, n_match=N_MATCH,
        gate="per-scale cohort Wilcoxon(obs_matched - drift_p50, greater); NO scale-max",
        fixes=["per-scale not scale-max", "distribution over placements", "SNR-matched obs",
               "non-adjacent random placements"], cohort=COHORT, bands=BANDS,
    ), indent=2))
    print(f"[D2-drift-FAIR] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
