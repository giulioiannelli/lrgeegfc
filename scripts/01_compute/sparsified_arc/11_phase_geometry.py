#!/usr/bin/env python3
"""Phase geometry -- the simplest, most interpretable cross-phase view. No LRG.

User observation (2026-07-12): the FC changes a lot rest_pre -> task, while
task_learn and task_test are very similar to each other. If learn~test DESPITE a
time gap, the task state is STABLE (a step at task onset, not a monotonic ramp) ->
the drift-ramp confound is empirically weak. Then the trace question is just:
does rest_post sit back at rest_pre (RESET) or hold the task state (TRACE)?

Measure = plain FC-state similarity. For each phase, imcoh_abs FC on two
non-overlapping halves -> upper-tri edge vectors. Cross-phase similarity uses
CROSS-HALF Spearman (p^A vs q^B, symmetrised) so no shared-half inflation;
within-phase split-half Spearman = the reliability ceiling. Build the 4x4 phase
geometry per (patient, band); cohort-average (Fisher-z).

Trace contrast (per patient, per band):
    toward_task = mean[ sim(post,learn), sim(post,test) ]
                - mean[ sim(pre ,learn), sim(pre ,test) ]
    > 0  => rest_post is more task-like than rest_pre was  (movement toward task)
Cohort Wilcoxon(toward_task, greater); band-specificity read off the per-band vector.

Outputs:
    data/sparsified_arc/trace_arc/phase_geom/{band}/{patient}.npz  (S 4x4, rel 4)
    data/sparsified_arc/trace_arc/phase_geom/geometry.csv          (cohort 4x4 + rel)
    data/sparsified_arc/trace_arc/phase_geom/toward_task.csv       (per band contrast+gate)
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch

OUT = ROOT / "data" / "sparsified_arc" / "trace_arc" / "phase_geom"
PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
PLAB = ("pre", "learn", "test", "post")


def _band_abs(Coh, freqs, flo, fhi):
    m = (freqs >= flo) & (freqs <= fhi)
    A = np.abs(Coh[:, :, m]).mean(axis=-1)
    A = 0.5 * (A + A.T); np.fill_diagonal(A, 0.0)
    return A


def _edge(W):
    return W[np.triu_indices_from(W, k=1)]


def _fz(r):   # Fisher z (clip for safety)
    return np.arctanh(np.clip(r, -0.999, 0.999))


def per_patient(pat):
    fs = FS_OVERRIDES.get(pat, 2048.0); nper = nperseg_for_fs(fs)
    # halves[band][phase] = (edgeA, edgeB)
    halves = {b: {} for b in BANDS}
    for ph in PHASES:
        try:
            X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
        except Exception:
            return []
        T = X.shape[1]; mid = T // 2
        vecs = {b: [] for b in BANDS}
        for a, z in ((0, mid), (mid, 2 * mid)):
            fr, C = compute_msc_welch(X[:, a:z], fs, nperseg=min(nper, z - a), metric="imcoh")
            for b in BANDS:
                flo, fhi = BRAIN_BANDS[b]
                vecs[b].append(_edge(_band_abs(C, fr, flo, fhi)))
            del C
        for b in BANDS:
            halves[b][ph] = vecs[b]
        del X

    rows = []
    for b in BANDS:
        H = halves[b]
        S = np.full((4, 4), np.nan); rel = np.full(4, np.nan)
        for i in range(4):
            rel[i] = spearmanr(H[PHASES[i]][0], H[PHASES[i]][1])[0]
            for j in range(4):
                if i == j:
                    S[i, j] = rel[i]
                else:
                    r1 = spearmanr(H[PHASES[i]][0], H[PHASES[j]][1])[0]
                    r2 = spearmanr(H[PHASES[i]][1], H[PHASES[j]][0])[0]
                    S[i, j] = 0.5 * (r1 + r2)
        cell = OUT / b; cell.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cell / f"{pat}.npz", S=S, rel=rel, labels=np.array(PLAB))
        post_task = 0.5 * (S[3, 1] + S[3, 2])   # sim(post, learn/test)
        pre_task = 0.5 * (S[0, 1] + S[0, 2])    # sim(pre , learn/test)
        rows.append(dict(patient=pat, band=b,
                         S=";".join(f"{v:.4f}" for v in S.ravel()),
                         rel=";".join(f"{v:.4f}" for v in rel),
                         sim_learn_test=float(S[1, 2]),
                         pre_task=float(pre_task), post_task=float(post_task),
                         toward_task=float(post_task - pre_task)))
    return rows


def cohort(df):
    geom, toward = [], []
    for b in BANDS:
        x = df[df.band == b]
        if x.empty:
            continue
        Ss = np.array([[float(v) for v in r.split(";")] for r in x.S]).reshape(-1, 4, 4)
        Sm = np.tanh(np.nanmean(_fz(Ss), axis=0))           # Fisher-z cohort mean 4x4
        rel = np.tanh(np.nanmean(_fz(np.array(
            [[float(v) for v in r.split(";")] for r in x.rel])), axis=0))
        for i in range(4):
            for j in range(4):
                geom.append(dict(band=b, row=PLAB[i], col=PLAB[j], sim=float(Sm[i, j])))
        tt = x.toward_task.values
        try:
            p = float(wilcoxon(tt, alternative="greater")[1])
        except Exception:
            p = np.nan
        toward.append(dict(band=b,
                           sim_learn_test=float(np.tanh(np.nanmean(_fz(x.sim_learn_test.values)))),
                           reliability=float(rel.mean()),
                           pre_task=float(np.tanh(np.nanmean(_fz(x.pre_task.values)))),
                           post_task=float(np.tanh(np.nanmean(_fz(x.post_task.values)))),
                           toward_task_med=float(np.median(tt)),
                           n_pos=int((tt > 0).sum()), gate_p=p,
                           trace=bool(np.isfinite(p) and p < 0.05)))
    return pd.DataFrame(geom), pd.DataFrame(toward)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    pats = COHORT[:a.limit] if a.limit else COHORT
    ncpu = int(os.environ.get("SA_WORKERS", 4))
    print(f"[phase-geom] {len(pats)} pt x {len(BANDS)} bands | 4-phase FC similarity, "
          f"split-half reliability | {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, rl in enumerate(pool.imap_unordered(per_patient, pats), 1):
            rows.extend(rl)
            print(f"[{i}/{len(pats)} pt] {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        print("[phase-geom] no cells", flush=True); return
    geom, toward = cohort(df)
    geom.to_csv(OUT / "geometry.csv", index=False)
    toward.to_csv(OUT / "toward_task.csv", index=False)
    print(f"\n[phase-geom] {len(rows)} cells in {time.time()-t0:.0f}s\n", flush=True)
    print("=== is the task state STABLE (learn~test vs reliability), and does POST hold it? ===", flush=True)
    print(toward[["band", "reliability", "sim_learn_test", "pre_task", "post_task",
                  "toward_task_med", "n_pos", "gate_p", "trace"]].to_string(index=False), flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="phase_geometry", measure="Spearman cross-half FC-edge similarity, 4 phases",
        phases=list(PHASES), note="no LRG/backbone; reliability = within-phase split-half.",
        cohort=COHORT, bands=BANDS), indent=2))
    print(f"[phase-geom] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
