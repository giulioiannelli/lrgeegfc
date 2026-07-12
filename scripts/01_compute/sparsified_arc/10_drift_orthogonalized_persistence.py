#!/usr/bin/env python3
"""Test 2 -- DRIFT-ORTHOGONALISED PERSISTENCE (drift-immune, keeps the monotonic signal).

Scope: .agents/guides/task-persistence-investigation/2026-07-12_drift-orthogonalised-persistence.md

The cross-phase design (rest_pre -> task -> rest_post) is time-ordered, so a trace is
MONOTONIC by construction -- and so is drift. Removing monotonic structure (the linear
detrend) removes the signal (null == alternative): void. Instead we remove the
SPONTANEOUS RESTING REPERTOIRE S (the directions the resting brain explores, which
CONTAIN drift), and ask whether the task pushed the network OUTSIDE S and rest_post
STAYED outside. The direction of change is NOT removed -> a monotonic trace survives.

  S       = span of top-r PCs of the rest_pre window trajectory   (the repertoire; drift in S)
  v       = x_task - xbar_pre                                      (task change)
  v_perp  = P_perp v ,  vhat = v_perp/||v_perp||                   (task change OUTSIDE S)
  T_perp  = < x_post - xbar_pre , vhat >                           (does post retain it?)

Fair CROSS-VALIDATED, SNR-MATCHED null (controls the under-sampling of S, caveat C-a):
  split rest_pre windows -> TRAIN (build S, xbar, vhat) and HELDOUT; obs uses POSTSUB
  (equal #windows). obs = <mean(POSTSUB)-xbar, vhat>; null = <mean(HELDOUT)-xbar, vhat>
  (a task-naive rest segment of the SAME size). Trace = obs > null, per scale r, cohort.

Two spaces (both reported): FC-edge (full imcoh_abs upper-tri; drift ~linear; drift-clean
check) and COPHENETIC@tau_min (LRG-on-backbone; the trace's native space). r-sweep so
drift-immunity is a curve, never a single picked value.

Outputs:
    data/sparsified_arc/trace_arc/orth_persist/{space}/{band}/{patient}.npz
        (r_grid, T_obs, T_null, novelty)
    data/sparsified_arc/trace_arc/orth_persist/gate_vs_r.csv
    data/sparsified_arc/trace_arc/orth_persist/verdict.csv
    data/sparsified_arc/trace_arc/orth_persist/config.json
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon, rankdata

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.fc.backbone import percolation_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale

OUT = ROOT / "data" / "sparsified_arc" / "trace_arc" / "orth_persist"
PHASES = ("rest_pre", "task_test", "rest_post")
SPACES = ("fc", "coph")
K = 12                 # windows per phase
B = 100                # CV splits (random TRAIN/HELDOUT/POSTSUB)
R_GRID = [0, 1, 2, 3, 4]   # repertoire dimension sweep (r < K//2)
SEED_BASE = 20260712
EPS = 1e-12


def _band_abs(Coh, freqs, flo, fhi):
    m = (freqs >= flo) & (freqs <= fhi)
    A = np.abs(Coh[:, :, m]).mean(axis=-1)
    A = 0.5 * (A + A.T); np.fill_diagonal(A, 0.0)
    return np.clip(A, 0.0, 1.0)


def _states_from_window(W):
    """(fc-edge vector, cophenetic@tau_min vector) for one band-window FC matrix."""
    iu = np.triu_indices_from(W, k=1)
    fc_vec = W[iu]
    try:
        ev, V = laplacian_eig(percolation_backbone(W)[0])
        coph_raw = cophenetic_at_scale(ev, V, 1.0)   # tau_min
        # rank-normalise: the 1/rho cophenetic magnitude explodes on near-disconnected
        # pairs and would dominate the linear projection; ranks are bounded and match
        # the Spearman/rank basis of rho_sym.
        coph_vec = rankdata(coph_raw) / coph_raw.size
    except Exception:
        coph_vec = None
    return fc_vec, coph_vec


def per_patient(pat):
    fs = FS_OVERRIDES.get(pat, 2048.0); nper = nperseg_for_fs(fs)
    # seqs[space][band][phase] = (K, m) array of state vectors
    seqs = {sp: {b: {ph: [] for ph in PHASES} for b in BANDS} for sp in SPACES}
    for ph in PHASES:
        X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
        T = X.shape[1]; L = T // K
        for w in range(K):
            fr, C = compute_msc_welch(X[:, w * L:(w + 1) * L], fs,
                                      nperseg=min(nper, L), metric="imcoh")
            for b in BANDS:
                flo, fhi = BRAIN_BANDS[b]
                fc_vec, coph_vec = _states_from_window(_band_abs(C, fr, flo, fhi))
                seqs["fc"][b][ph].append(fc_vec)
                seqs["coph"][b][ph].append(coph_vec)
            del C
        del X

    rng = np.random.default_rng(SEED_BASE + COHORT.index(pat))
    half = K // 2
    rows = []
    for sp in SPACES:
        for b in BANDS:
            pre = seqs[sp][b]["rest_pre"]; task = seqs[sp][b]["task_test"]; post = seqs[sp][b]["rest_post"]
            if any(v is None for v in pre + task + post):
                continue
            pre = np.asarray(pre); task = np.asarray(task); post = np.asarray(post)
            task_mean = task.mean(0)
            obs = np.full((B, len(R_GRID)), np.nan)
            null = np.full((B, len(R_GRID)), np.nan)
            nov = np.full((B, len(R_GRID)), np.nan)
            for s in range(B):
                perm = rng.permutation(K)
                tr, hd = perm[:half], perm[half:]
                ps = rng.choice(K, half, replace=False)      # POSTSUB (equal size)
                xbar = pre[tr].mean(0)
                Mc = pre[tr] - xbar
                # right singular vectors = repertoire directions in m-space
                _, _, Vt = np.linalg.svd(Mc, full_matrices=False)
                vfull = task_mean - xbar
                d_obs = post[ps].mean(0) - xbar
                d_null = pre[hd].mean(0) - xbar
                nvf = np.linalg.norm(vfull) + EPS
                for ri, r in enumerate(R_GRID):
                    if r == 0:
                        vperp = vfull
                    else:
                        Vr = Vt[:r]
                        vperp = vfull - Vr.T @ (Vr @ vfull)
                    nv = np.linalg.norm(vperp)
                    if nv < EPS:
                        continue
                    vhat = vperp / nv
                    obs[s, ri] = float(d_obs @ vhat)
                    null[s, ri] = float(d_null @ vhat)
                    nov[s, ri] = float(nv / nvf)
            T_obs = np.nanmean(obs, 0); T_null = np.nanmean(null, 0); NOV = np.nanmean(nov, 0)
            cell = OUT / sp / b; cell.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(cell / f"{pat}.npz", r_grid=np.array(R_GRID),
                                T_obs=T_obs, T_null=T_null, novelty=NOV)
            rows.append(dict(space=sp, band=b, patient=pat,
                             T_obs=";".join(f"{v:.5f}" for v in T_obs),
                             T_null=";".join(f"{v:.5f}" for v in T_null),
                             novelty=";".join(f"{v:.4f}" for v in NOV)))
    return rows


def cohort(df):
    gate, verdict = [], []
    for sp in SPACES:
        for b in BANDS:
            x = df[(df.space == sp) & (df.band == b)]
            if x.empty:
                continue
            OBS = np.array([[float(v) for v in r.split(";")] for r in x.T_obs])
            NUL = np.array([[float(v) for v in r.split(";")] for r in x.T_null])
            NOV = np.array([[float(v) for v in r.split(";")] for r in x.novelty])
            for ri, r in enumerate(R_GRID):
                d = OBS[:, ri] - NUL[:, ri]; d = d[np.isfinite(d)]
                p = np.nan
                if d.size >= 5 and np.any(d != 0):
                    try:
                        p = float(wilcoxon(d, alternative="greater")[1])
                    except Exception:
                        pass
                gate.append(dict(space=sp, band=b, r=r, gate_p=p,
                                 T_obs_med=float(np.nanmedian(OBS[:, ri])),
                                 T_null_med=float(np.nanmedian(NUL[:, ri])),
                                 novelty_med=float(np.nanmedian(NOV[:, ri])),
                                 n_pat_obs_gt_null=int(np.nansum(OBS[:, ri] > NUL[:, ri]))))
            g = pd.DataFrame([row for row in gate if row["space"] == sp and row["band"] == b])
            best = g.loc[g.gate_p.idxmin()] if g.gate_p.notna().any() else None
            verdict.append(dict(space=sp, band=b,
                                best_r=int(best.r) if best is not None else -1,
                                best_gate_p=float(best.gate_p) if best is not None else np.nan,
                                T_obs_at_best=float(best.T_obs_med) if best is not None else np.nan,
                                T_null_at_best=float(best.T_null_med) if best is not None else np.nan,
                                novelty_at_best=float(best.novelty_med) if best is not None else np.nan,
                                survives=bool(best is not None and best.gate_p < 0.05)))
    return pd.DataFrame(gate), pd.DataFrame(verdict)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    pats = COHORT[:a.limit] if a.limit else COHORT
    ncpu = int(os.environ.get("SA_WORKERS", 4))
    print(f"[Test2-orth] {len(pats)} pt x {len(BANDS)} bands x {len(SPACES)} spaces | "
          f"K={K} win, B={B} CV splits, r={R_GRID} | held-out-rest SNR-matched null | "
          f"{ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, rl in enumerate(pool.imap_unordered(per_patient, pats), 1):
            rows.extend(rl)
            el = time.time() - t0
            print(f"[{i}/{len(pats)} pt] {el:.0f}s ETA {el/i*(len(pats)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        print("[Test2-orth] no cells", flush=True); return
    gate, verdict = cohort(df)
    gate.to_csv(OUT / "gate_vs_r.csv", index=False)
    verdict.to_csv(OUT / "verdict.csv", index=False)
    print(f"\n[Test2-orth] {len(rows)} cells in {time.time()-t0:.0f}s\n", flush=True)
    print("=== DRIFT-ORTHOGONALISED PERSISTENCE: obs(post) vs held-out-rest null, best r ===", flush=True)
    print(verdict.to_string(index=False), flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="Test2_drift_orthogonalised_persistence",
        scope=".agents/guides/task-persistence-investigation/2026-07-12_drift-orthogonalised-persistence.md",
        K=K, B=B, r_grid=R_GRID, spaces=list(SPACES),
        null="cross-validated held-out-rest, SNR-matched (equal window counts)",
        note="removes the resting repertoire S, NOT the direction of change -> drift-immune, monotonic-safe",
        cohort=COHORT, bands=BANDS), indent=2))
    print(f"[Test2-orth] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
