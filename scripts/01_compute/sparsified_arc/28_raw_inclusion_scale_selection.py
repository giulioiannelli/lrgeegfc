#!/usr/bin/env python3
"""Raw FC as the fine-scale limit, and the scale-selection that turns it into the trace.

The point (user frame, 2026-07-23): raw FC is NOT a rival to the multiscale cophenetic
read -- it is provably its s->0 limit (e^{-tauL} ~ I - tauL, whose off-diagonal is W).
Raw FC registers reorganization NON-selectively (fires in every band); the diffusion,
read across scales, acts as a SELECTION that keeps only what is coherent across scales
(alpha, beta) and washes out the scale-local firing of the rest. This script shows the
inclusion + the morph + the selection, on the SAME mst@0.20 backbone + matched-strength
null as the trace gate, and it also brackets the GLOBAL methods (geodesic/resistance/
Grassmann; assembled from the ladder + probe table) as the COARSE limit of the same axis.

CHECK 1 -- inclusion + morph + selection:
  extend the trace SGRID BELOW s=1 (s in [0.05, 180]); per band, cophenetic rho_sym(s)
  under matched-strength, PLUS the raw-edge rho_sym (the s->0 reference). Read:
    (a) at fine s does the cophenetic band-pattern approach raw's all-band pattern?
    (b) does it morph continuously to alpha/beta as s grows (the selection made visible)?
    (c) at coarse s, which bands survive (the "global" limit to compare to Grassmann)?
CHECK 2 -- the selecting MECHANISM (descriptive, obs only, no null):
  per band, scale-coherence C of the per-pair task-reorganization fingerprint
  u(s) = D_task(s) - D_pre(s): mean cross-scale Spearman of u(s). Prediction: bands that
  survive to coarse scales (alpha, beta) are scale-coherent (high C); washed bands are
  scale-local (low C). Reported full-grid and s>=1. Cross-band corr(C, #scales-cleared).

5-point preamble
1. Claim: cophenetic rho_sym(s) -> raw rho_sym as s->0 (inclusion), morphs to alpha/beta
   as s grows (selection), and the selection tracks per-band scale-coherence (mechanism).
2. Null: the whole band-pattern is a node-strength artifact -> IDENTICAL matched-strength
   surrogate through the identical sparsify+LRG pipeline (audit_150 bit-identical shuffle).
3. Strongest alt it must beat: raw non-selectivity is real signal the hierarchy wrongly
   discards. Controlled by (a) reporting raw under the SAME null (raw fires >2 bands ->
   non-selective), and (b) the scale-coherence check -- if discarded bands are as coherent
   as alpha/beta, the selection is unprincipled and we DROP the multiscale claim.
4. Cannot: cophenetic is a deterministic (lossy) transform of raw, so coph carries NO
   information raw lacks (script 21 confirmed: coph|raw ~ 0). This is the point -- the
   contribution is SELECTION/precision, not new detection. Scale-coherence C over the full
   grid is partly self-referential (fine scales ~ raw); the s>=1 C and the task-side (no
   persistence) fingerprint are the less-circular reads.
5. Falsify: if fine-s cophenetic does NOT approach raw's pattern, or there is no smooth
   morph, or discarded bands are as scale-coherent as alpha/beta -> the inclusion/selection
   story fails and must not be written.

Outputs (data/sparsified_arc/raw_inclusion_mst020/):
  coph_scale.csv   per (patient, band, s): obs_rho, surr_p50/p95, p  (cophenetic sweep)
  raw_ref.csv      per (patient, band): obs_raw, surr_p50, surr_p95, p  (s->0 reference)
  coherence.csv    per (patient, band): C_full, C_meso (scale-coherence of task reorg)
  cohort_gate.csv  per (band, s): gate_p, n_above, obs_med, surr_med  (+ s=raw row)
  config.json
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import wilcoxon, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS, load_phase            # noqa: E402
from lrg_eegfc.utils.fc.backbone import select_backbone                 # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import (                        # noqa: E402
    laplacian_eig, rho_sym_over_scales, rho_sym, cophenetic_at_scale)
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle  # noqa: E402

BACKBONE = os.environ.get("SA_BACKBONE", "mst020")
FRAC = float(os.environ.get("SA_FRAC", "0.20"))
OUT = ROOT / "data" / "sparsified_arc" / "raw_inclusion_mst020"
PHASES = ("A", "B", "task_test", "rest_post")
# Extended grid: reach BELOW s=1 (raw limit) up to the coarse/global limit.
SGRID = np.logspace(np.log10(0.05), np.log10(180.0), 28)
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260723
_IU = {}


def _triu(M):
    N = M.shape[0]
    if N not in _IU:
        _IU[N] = np.triu_indices(N, 1)
    return M[_IU[N]]


def eig_bb(W):
    return laplacian_eig(select_backbone(W, BACKBONE, frac=FRAC))


def raw_rho(Ws):
    """Raw-edge cross-phase trace = the s->0 reference (dense |ImCoh| upper-tri)."""
    A = {ph: _triu(Ws[ph]) for ph in PHASES}
    r, _ = rho_sym(A["A"], A["B"], A["task_test"], A["rest_post"])
    return r


def scale_coherence(eig_obs, sgrid):
    """Cross-scale Spearman of the symmetric per-pair task-reorg fingerprint
    u(s) = 0.5[(D_task-D_A)+(D_task-D_B)] at scale s. Returns (C_full, C_meso)."""
    cols, keep = [], []
    for s in sgrid:
        try:
            Dt = cophenetic_at_scale(*eig_obs["task_test"], s)
            DA = cophenetic_at_scale(*eig_obs["A"], s)
            DB = cophenetic_at_scale(*eig_obs["B"], s)
            u = 0.5 * ((Dt - DA) + (Dt - DB))
        except Exception:
            u = None
        if u is not None and np.all(np.isfinite(u)) and np.std(u) > 0:
            cols.append(u); keep.append(s)
    if len(cols) < 2:
        return np.nan, np.nan
    U = np.array(cols); ks = np.array(keep)

    def _meanpair(idx):
        if idx.sum() < 2:
            return np.nan
        M = U[idx]
        rs = []
        for i in range(len(M)):
            for j in range(i + 1, len(M)):
                r, _ = spearmanr(M[i], M[j])
                if np.isfinite(r):
                    rs.append(r)
        return float(np.mean(rs)) if rs else np.nan

    return _meanpair(np.ones(len(ks), bool)), _meanpair(ks >= 1.0)


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    if any(v.shape[0] != N for v in Ws.values()):
        return None
    eig_obs = {ph: eig_bb(Ws[ph]) for ph in PHASES}
    obs = rho_sym_over_scales(eig_obs, SGRID)          # (nS,)
    obs_raw = raw_rho(Ws)
    C_full, C_meso = scale_coherence(eig_obs, SGRID)

    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    surr = np.full((R, SGRID.size), np.nan)
    surr_raw = np.full(R, np.nan)
    for r in range(R):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX) for ph in PHASES}
        try:
            surr[r] = rho_sym_over_scales({ph: eig_bb(Wsh[ph]) for ph in PHASES}, SGRID)
            surr_raw[r] = raw_rho(Wsh)
        except Exception:
            continue

    scale_rows = []
    for j, s in enumerate(SGRID):
        col = surr[:, j]; col = col[np.isfinite(col)]; o = obs[j]
        scale_rows.append(dict(
            patient=pat, band=band, s=float(s),
            obs_rho=float(o) if np.isfinite(o) else np.nan,
            surr_p50=float(np.nanpercentile(col, 50)) if col.size else np.nan,
            surr_p95=float(np.nanpercentile(col, 95)) if col.size else np.nan,
            p=float(np.mean(col >= o)) if col.size and np.isfinite(o) else np.nan))
    cr = surr_raw[np.isfinite(surr_raw)]
    raw_row = dict(patient=pat, band=band, obs_raw=float(obs_raw),
                   surr_p50=float(np.nanpercentile(cr, 50)) if cr.size else np.nan,
                   surr_p95=float(np.nanpercentile(cr, 95)) if cr.size else np.nan,
                   p=float(np.mean(cr >= obs_raw)) if cr.size and np.isfinite(obs_raw) else np.nan)
    coh_row = dict(patient=pat, band=band, C_full=C_full, C_meso=C_meso)
    return scale_rows, raw_row, coh_row


def _gate(df, valcol="obs_rho", surrcol="surr_p50"):
    d = df.dropna(subset=[valcol, surrcol])
    if len(d) < 5:
        return None
    diff = d[valcol].values - d[surrcol].values
    try:
        p = float(wilcoxon(diff, alternative="greater")[1]) if np.any(diff != 0) else np.nan
    except Exception:
        p = np.nan
    return dict(gate_p=p, n_above=int((d.p < 0.05).sum()), n_pat=len(d),
                obs_med=float(d[valcol].median()), surr_med=float(d[surrcol].median()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--bands", type=str, default="")
    ap.add_argument("--R", type=int, default=0)
    a = ap.parse_args()
    global R
    if a.R:
        R = a.R
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    pats = COHORT[:a.limit] if a.limit else COHORT
    bands = a.bands.split(",") if a.bands else list(BANDS)
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in bands for p in pats)]
    ncpu = int(os.environ.get("SA_WORKERS", 12))
    print(f"[raw-incl] {len(jobs)} cells, R={R}, {len(SGRID)} scales "
          f"(s={SGRID[0]:.3f}..{SGRID[-1]:.0f}), {ncpu} workers -> {OUT.name}", flush=True)
    t0 = time.time()
    scale_rows, raw_rows, coh_rows = [], [], []
    with Pool(ncpu) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if res:
                sr, rr, cr = res
                scale_rows.extend(sr); raw_rows.append(rr); coh_rows.append(cr)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
    ds = pd.DataFrame(scale_rows); dr = pd.DataFrame(raw_rows); dc = pd.DataFrame(coh_rows)
    if ds.empty:
        print("[raw-incl] no cells", flush=True); return
    ds.to_csv(OUT / "coph_scale.csv", index=False)
    dr.to_csv(OUT / "raw_ref.csv", index=False)
    dc.to_csv(OUT / "coherence.csv", index=False)

    gate = []
    for band in bands:
        rr = dr[dr.band == band]
        g = _gate(rr, "obs_raw", "surr_p50")
        if g:
            gate.append(dict(band=band, s=-1.0, kind="raw", **g))
        for s in np.sort(ds.s.unique()):
            x = ds[(ds.band == band) & np.isclose(ds.s, s)]
            g = _gate(x)
            if g:
                gate.append(dict(band=band, s=float(s), kind="coph", **g))
    gdf = pd.DataFrame(gate); gdf.to_csv(OUT / "cohort_gate.csv", index=False)
    if gdf.empty or "band" not in gdf.columns:
        print(f"\n[raw-incl] CSVs written; cohort summary needs >=5 patients "
              f"(subset/timing run) -> {OUT}", flush=True)
        return

    print("\n=== CHECK 1: band-pattern from raw (s->0) through the scale sweep ===", flush=True)
    print("   (cells = cohort gate_p; * = clears <0.05; raw row first, then s ascending)", flush=True)
    show_s = [-1.0] + [float(s) for s in np.sort(ds.s.unique())]
    hdr = "band        " + "".join(f"{('raw' if s < 0 else f'{s:.2f}'):>8s}" for s in show_s)
    print(hdr, flush=True)
    for band in bands:
        cells = []
        for s in show_s:
            row = gdf[(gdf.band == band) & np.isclose(gdf.s, s)]
            if row.empty or not np.isfinite(row.gate_p.iloc[0]):
                cells.append(f"{'--':>8s}")
            else:
                p = row.gate_p.iloc[0]
                cells.append(f"{('*' if p < 0.05 else ' ')}{p:6.3f} ")
        print(f"{band:11s} " + "".join(cells), flush=True)

    print("\n=== CHECK 2: scale-coherence of task reorganization (obs; mechanism) ===", flush=True)
    print("   C = mean cross-scale Spearman of per-pair (D_task - D_pre); higher = survives", flush=True)
    ncleared = {b: int((gdf[(gdf.band == b) & (gdf.kind == 'coph')].gate_p < 0.05).sum()) for b in bands}
    for band in bands:
        cc = dc[dc.band == band]
        print(f"  {band:11s} C_full={cc.C_full.median():+.3f}  C_meso(s>=1)={cc.C_meso.median():+.3f}"
              f"   scales_cleared={ncleared[band]}/{len(SGRID)}", flush=True)
    bvals = [(dc[dc.band == b].C_meso.median(), ncleared[b]) for b in bands
             if np.isfinite(dc[dc.band == b].C_meso.median())]
    if len(bvals) >= 3:
        cc, nn = zip(*bvals)
        rho, pp = spearmanr(cc, nn)
        print(f"  --> cross-band Spearman(C_meso, scales_cleared) = {rho:+.3f} (p={pp:.3f})", flush=True)

    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="raw_inclusion_scale_selection", backbone=BACKBONE, frac=FRAC,
        s_grid=[float(x) for x in SGRID], R=R, swap_factor=SWAP_FACTOR,
        cohort=pats, bands=bands, base_seed=BASE_SEED,
        null="matched-strength on dense FC -> select_backbone -> LRG cophenetic rho_sym(s); raw = dense upper-tri rho_sym"),
        indent=2))
    print(f"\n[raw-incl] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
