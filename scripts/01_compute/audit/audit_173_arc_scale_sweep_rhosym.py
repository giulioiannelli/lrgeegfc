#!/usr/bin/env python3
"""audit_173 — Part B: encoding/inference components over a DENSE scale sweep (rho_sym x s x MS).

Generalises audit_157 (rho_sym arc null at 3 tau: 1.0, 2.610, 6.813) to a dense
s-grid, exactly as audit_172 generalised the audit_150 trace gate. Same statistic
(the symmetric two-arm functionals T_learn = encoding echo, T_infspec_pe =
inference-specific | encoding), same 5-phase arc, same CACHED seed-20260511
matched-strength surrogate eigenpairs (cophenetic reformed at each s from cached
eigs, NO regeneration -> cheap). Question: do the encoding / inference components
have a characteristic scale, and does the beta inference component's mesoscale
preference (audit_157: s~2.6) hold up across the full collapse-free window?

s = tau*lam_max (= audit_157's tau_mult). Anchors: arm1 at s in {1.0,2.610,6.813}
reproduces audit_103d (rho_split); sym at s=1 reproduces audit_152 (rho_sym fine).
Read off the collapse-free window + SHAPE (peak = scale-tuned; monotone = collapse),
NOT the min-p over s. 5-point preamble = audit_157's (inherited); scope Part B:
.agents/guides/task-persistence-investigation/2026-07-11_ms-relative-scale-sweep.md

Does NOT edit audit_103*/152/157. Output: data/audit/arc_scale_sweep_rhosym/
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
import _wm_stratify as ws  # type: ignore
from audit_103_cophenetic_consolidation_arc import (  # type: ignore
    _arc_cell_rng, ensure_half_fcs, load_phase_fc,
)
from audit_103c_tau_sweep_arc import _cophenet_at_taumult, _eig_L  # type: ignore
from audit_103d_arc_mesoscale_null import TAU_MULTS  # type: ignore  (1.000,2.610,6.813)
from audit_152_consolidation_arc_rhosym import (  # type: ignore
    ARC_PHASES, BANDS, COHORT, FUNCTIONALS, _sym_functionals, _surr_functionals_sym,
)

OUT = ROOT / "data" / "audit" / "arc_scale_sweep_rhosym"
REF_MESO = ROOT / "data" / "audit" / "consolidation_arc" / "arc_mesoscale_null_R200.csv"
REF_ARC1 = ROOT / "data" / "audit" / "consolidation_arc_rhosym" / "arc_null_per_patient.csv"

# dense grid covering the collapse-free window; INCLUDES the 3 audit_157 anchor s
SGRID = np.unique(np.concatenate([[0.5, 0.7, 0.85],
                                  np.geomspace(1.0, 8.0, 12),
                                  list(TAU_MULTS)]))       # -> 1.0, 2.610, 6.813 present
S1 = float(SGRID[np.argmin(np.abs(SGRID - 1.0))])
HEAD = ("T_infspec_pe", "T_learn")


def _cell(job):
    pat, band = job
    R, sf = ws.N_SURROGATES, ws.SWAP_FACTOR
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}
        obs_eig = {ph: _eig_L(Ws[ph]) for ph in ARC_PHASES}
        eigs = {ph: load_or_compute_eigs_at_path(
            ws.surr_eig_path("full", pat, band, ph, R, sf), Ws[ph], R, sf,
            _arc_cell_rng(pat, band, ph), verbose=False) for ph in ARC_PHASES}
    except Exception as exc:  # noqa: BLE001
        print(f"[audit_173] FAIL {pat}/{band}: {type(exc).__name__}: {exc}", flush=True)
        return None
    rows = []
    for s in SGRID:
        D = {ph: _cophenet_at_taumult(*obs_eig[ph], s)[0] for ph in ARC_PHASES}
        sym_obs, arm1_obs = _sym_functionals(D)
        sc = {}
        for ph in ARC_PHASES:
            evals, evecs = eigs[ph]
            npair = Ws[ph].shape[0] * (Ws[ph].shape[0] - 1) // 2
            C = np.full((R, npair), np.nan)
            for r in range(R):
                if np.isfinite(evals[r]).all():
                    C[r] = _cophenet_at_taumult(evals[r], evecs[r], s)[0]
            sc[ph] = C
        surr = _surr_functionals_sym(sc, R)
        row = {"patient": pat, "band": band, "s": float(s), "tau_mult": float(s)}
        for k in FUNCTIONALS:
            sv = surr[k][np.isfinite(surr[k])]
            row[f"{k}_obs"] = sym_obs[k]
            row[f"{k}_arm1"] = arm1_obs[k]
            row[f"{k}_surr_p50"] = float(np.median(sv)) if sv.size else np.nan
            row[f"{k}_p"] = float(np.mean(sv >= sym_obs[k])) if sv.size else np.nan
        rows.append(row)
    return rows


def _anchors(df):
    print("=== CORRECTNESS ANCHORS ===")
    if REF_MESO.exists():
        ref = pd.read_csv(REF_MESO)
        m = df.merge(ref, on=["patient", "band", "tau_mult"], suffixes=("", "_ref"))
        d = [np.abs(m[f"{k}_arm1"] - m[f"{k}_obs_ref"]).max()
             for k in FUNCTIONALS if f"{k}_arm1" in m and f"{k}_obs_ref" in m]
        md = float(np.nanmax(d)) if d else np.nan
        print(f"  A1 arm1==audit_103d (rho_split @ 3 shared s): n={len(m)} "
              f"max|Δ|={md:.2e} → {'PASS' if md < 1e-9 else 'CHECK'}")
    if REF_ARC1.exists():
        ref = pd.read_csv(REF_ARC1).set_index(["patient", "band"])
        cur = df[np.isclose(df.s, 1.0)].set_index(["patient", "band"])
        d = [abs(r[f"{k}_obs"] - float(ref.loc[idx, f"{k}_obs"]))
             for k in FUNCTIONALS for idx, r in cur.iterrows()
             if idx in ref.index and f"{k}_obs" in ref.columns]
        md = float(np.nanmax(d)) if d else np.nan
        print(f"  A2 sym(s=1)==audit_152 (rho_sym fine): n={len(cur)} "
              f"max|Δ|={md:.2e} → {'PASS' if md < 1e-9 else 'CHECK'}")
    print(flush=True)


def _verdict(df):
    vrows = []
    print("=== COHORT MATCHED-STRENGTH VERDICT under rho_sym over s (positive = trace) ===")
    for k in HEAD:
        print(f"\n--- {k} : cohort Wilcoxon p(s) ---")
        for band in ("beta", "alpha", "low_gamma", "theta", "delta", "high_gamma"):
            line = f"  {band:<11}"
            for s in SGRID:
                d = df[(df.band == band) & np.isclose(df.s, s)].dropna(
                    subset=[f"{k}_obs", f"{k}_surr_p50"])
                if d.empty:
                    p = np.nan
                else:
                    n_above = int((d[f"{k}_p"] < 0.05).sum())
                    v = ws.cohort_verdict(d[f"{k}_obs"].values, d[f"{k}_surr_p50"].values, n_above)
                    p = v["wilcoxon_p"]
                    vrows.append({"band": band, "s": float(s), "functional": k,
                                  "med_obs": v["med_obs"], "med_surr": v["med_surr"],
                                  "n_above": n_above, "wilcoxon_p": p})
                if s in (S1, TAU_MULTS[1], TAU_MULTS[2]) or (band == "beta"):
                    line += f" {s:.2f}:{p:.3f}"
            print(line, flush=True)
    pd.DataFrame(vrows).to_csv(OUT / "arc_scale_cohort_verdict_s.csv", index=False)
    # headline: beta inference & alpha/beta encoding — s=1 vs collapse-free best
    print("\n=== HEADLINE: s=1 vs collapse-free-window best (s < ~3.3 Fiedler) ===")
    vd = pd.DataFrame(vrows)
    for k in HEAD:
        for band in ("beta", "alpha"):
            g = vd[(vd.functional == k) & (vd.band == band)]
            if g.empty:
                continue
            g1 = g[np.isclose(g.s, 1.0)].iloc[0]
            free = g[g.s < 3.3]
            gb = free.loc[free.wilcoxon_p.idxmin()]
            print(f"  {k:<13} {band:<6} s=1 p={g1.wilcoxon_p:.4f} | "
                  f"best(free) p={gb.wilcoxon_p:.4f} @s={gb.s:.2f}", flush=True)


def _run(cohort):
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, b) for b in BANDS for p in cohort]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_173] {len(jobs)} cells, R={ws.N_SURROGATES}, |s|={len(SGRID)} "
          f"(s in [{SGRID[0]:.2f},{SGRID[-1]:.1f}]), {ncpu} workers "
          f"(cached seed-20260511 eigs, reform-only)", flush=True)
    for pat in cohort:
        ensure_half_fcs(pat, BANDS)
    t0, rows = time.time(), []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(_cell, jobs), 1):
            if r:
                rows.extend(r)
                print(f"[{i}/{len(jobs)}] {r[0]['patient']}/{r[0]['band']} done", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient", "s"])
    df.to_csv(OUT / "arc_scale_null_s.csv", index=False)
    print(f"\n[audit_173] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT}\n", flush=True)
    _anchors(df)
    _verdict(df)


def _pilot():
    print(f"[pilot] Pat_08/beta, R={ws.N_SURROGATES}, |s|={len(SGRID)} ...", flush=True)
    ensure_half_fcs("Pat_08", ["beta"])
    t0 = time.time(); rows = _cell(("Pat_08", "beta")); dt = time.time() - t0
    if rows:
        r = pd.DataFrame(rows)[["s", "T_infspec_pe_obs", "T_infspec_pe_p",
                                "T_learn_obs", "T_learn_p"]]
        print(r.to_string(index=False))
    print(f"\n[pilot] {dt:.1f}s/cell -> 60 cells Pool(14) ~{dt*60/14/60:.1f} min", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "pilot":
        _pilot()
    else:
        _run(COHORT)
