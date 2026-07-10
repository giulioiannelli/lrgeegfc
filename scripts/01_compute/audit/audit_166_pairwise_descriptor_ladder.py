#!/usr/bin/env python3
"""audit_166 — pairwise-descriptor ladder: do the usual low-level network scalars
reproduce the multiscale cophenetic trace?

Runs a ladder of classical pairwise-graph descriptors through the IDENTICAL
cross-phase rho_sym estimator and the IDENTICAL matched-strength null as the LRG
cophenetic trace, ordered by structural order (raw edges -> node scalars ->
cophenetic). Every rung is scored on the four signatures that make the cophenetic
result a result: (C1) recovers the alpha whole-task gate, (C2) stays silent on the
clean-null band theta, (C3) survives matched-strength, (C5) reproduces the beta-only
inference-specific refinement T_infspec_pe. (C4 beta->OFC is a downstream
localization, not run here.)

Design (scope: .agents/guides/task-persistence-investigation/2026-07-09_pairwise-descriptor-ladder.md)
- Estimator: rho_sym arc functionals, copied verbatim from audit_152 (T_test =
  whole-task trace = rho_sym(A,B,TT,RP); T_learn, T_infspec, T_infspec_pe).
- Representation swap: only phi(W) changes across rungs (graph_descriptors + the
  LRG cophenetic). raw_fc/cophenetic are pair-level (N(N-1)/2); the classical
  scalars are node-level (N). rho_sym is representation-agnostic (Spearman).
- Null: MODE B. Reuse the EXACT flagship matched-strength ensemble on disk
  (seed 20260511, R=200, swap 20; cache 100% complete for all 5 phases). Per
  surrogate, reconstruct the adjacency from the cached Laplacian eigendecomposition
  (adjacency_from_laplacian_eigs) and recompute every descriptor on it; the
  cophenetic rung uses cophenetic_condensed_from_eigs directly. So every rung sees
  identical surrogate draws AND the cophenetic rung reproduces audit_152 bit-for-bit
  (anti-hallucination gate).

5-point critical preamble
1. Claim: no classical low-level pairwise descriptor reproduces the cophenetic
   discrimination (band-selectivity AND inference-beta-specificity) under matched
   strength; only the cophenetic rung passes.
2. Null: matched-strength (4-cycle +/-delta, strength preserving) reproduces any
   trace that is a function of node strength.
3. Strongest alternatives: (i) strength-slaving (our story = re-encoded strength;
   MS kills it); (ii) a non-strength scalar genuinely recovers the trace (would
   falsify "multiscale necessary" -> run the FULL menu); (iii) presence !=
   discrimination (raw FC is positive in all bands; score on selectivity).
4. Null's reach: (i) controlled BY CONSTRUCTION -- node strength is degenerate
   under MS (surrogate preserves strength exactly -> reproduces the strength-trace
   -> p~1); (ii) not auto-covered -> covered by running it; (iii) covered by the
   four-criterion score + the raw-FC broad-positive reference.
5. Falsify: any classical rung that survives MS AND is band-selective AND
   reproduces beta-only T_infspec_pe breaks the claim -> reported honestly.

Outputs: data/audit/pairwise_descriptor_ladder/
  per_patient_per_band.csv   per (descriptor, patient, band): obs + surrogate p, 4 functionals
  cohort_summary.csv         per (descriptor, band, functional): Wilcoxon gate
  strength_slaving.csv       per descriptor: cohort-median Spearman(phi, strength)
  README.md
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS
from lrg_eegfc.utils.surrogate.matched_strength import (
    adjacency_from_laplacian_eigs,
    cophenetic_condensed_from_eigs,
    load_or_compute_surrogate_eigs,
)
from lrg_eegfc.utils.metrics.graph_descriptors import (
    raw_edges,
    node_strength,
    weighted_clustering_onnela,
    eigenvector_centrality,
    pagerank,
    closeness_centrality,
    betweenness_centrality,
)

# Reuse the proven audit_63 loaders + LRG ultrametric (import, no fork) exactly as
# audit_152 does, so the observed cophenetic rung is bit-identical to the flagship.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
    lrg_ultrametric_condensed,
)

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = list(BRAIN_BANDS)
ARC_PHASES = ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")
FUNCTIONALS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")

R = 200
SWAP_FACTOR = 20
SEED = 20260511          # canonical cached ensemble (Mode B: exact flagship draws)
OUT = ROOT / "data" / "audit" / "pairwise_descriptor_ladder"

# Ladder rungs in structural order (bottom -> top). Classical node/edge descriptors;
# "cophenetic" is handled specially (from eigs). Betweenness is opt-in (--betweenness).
CLASSICAL = [
    ("raw_fc", raw_edges),
    ("strength", node_strength),
    ("clustering_onnela", weighted_clustering_onnela),
    ("eigcent", eigenvector_centrality),
    ("pagerank", pagerank),
    ("closeness", closeness_centrality),
]
RUNGS = [d for d, _ in CLASSICAL] + ["cophenetic"]


# --- rho_sym arc functionals (verbatim from audit_152) ---------------------------
def _rho(a: np.ndarray, b: np.ndarray) -> float:
    r, _ = spearmanr(a, b)
    return float(r)


def _partial_rho(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """First-order partial Spearman of (a, b) controlling c (rank world)."""
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    denom = np.sqrt(max(0.0, (1.0 - rac**2) * (1.0 - rbc**2)))
    return (rab - rac * rbc) / denom if denom > 0 else np.nan


def _arm_functionals(D_TL, D_TT, D_RP, D_bt, D_br) -> dict:
    e = D_TL - D_bt          # encoding
    g = D_TT - D_bt          # inference-online
    f = D_TT - D_TL          # inference-specific (arm-invariant)
    p = D_RP - D_br          # persistent
    return {"T_test": _rho(g, p), "T_learn": _rho(e, p),
            "T_infspec": _rho(f, p), "T_infspec_pe": _partial_rho(f, p, e)}


def _sym_functionals(D: dict) -> dict:
    """Mean over the two split-half arm assignments (rho_sym)."""
    arm1 = _arm_functionals(D["task_learn"], D["task_test"], D["rest_post"],
                            D["rest_pre_A"], D["rest_pre_B"])
    arm2 = _arm_functionals(D["task_learn"], D["task_test"], D["rest_post"],
                            D["rest_pre_B"], D["rest_pre_A"])
    return {k: 0.5 * (arm1[k] + arm2[k]) for k in FUNCTIONALS}


# --- representation bank ----------------------------------------------------------
def _phi(desc_name, desc_fn, W=None, eigs=None):
    """Descriptor vector for a phase. Classical descriptors act on W; 'cophenetic'
    acts on the Laplacian eigendecomposition (eigs=(evals, evecs)) or W."""
    if desc_name == "cophenetic":
        if eigs is not None:
            return cophenetic_condensed_from_eigs(eigs[0], eigs[1])
        return lrg_ultrametric_condensed(W)
    return desc_fn(W)


def per_cell(job):
    idx, pat, band, descriptors = job
    try:
        W = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}
    except Exception as exc:  # noqa: BLE001
        return {"error": f"{pat}/{band} load: {type(exc).__name__}: {exc}"}
    N = W["rest_pre_A"].shape[0]
    desc_map = dict(descriptors)  # name -> fn (classical); "cophenetic" special

    # ---- observed descriptor vectors + observed functionals ----
    obs_vec = {}
    for name in list(desc_map) + ["cophenetic"]:
        fn = desc_map.get(name)
        obs_vec[name] = {ph: _phi(name, fn, W=W[ph]) for ph in ARC_PHASES}
    strength_obs = {ph: node_strength(W[ph]) for ph in ARC_PHASES}

    obs_ff = {name: _sym_functionals(obs_vec[name]) for name in obs_vec}
    # strength-slaving: cohort will aggregate; here median-over-phases per descriptor
    slave = {}
    for name in obs_vec:
        if name in ("raw_fc", "cophenetic"):
            slave[name] = np.nan  # pair-level: not a node-vs-strength comparison
            continue
        rs = [abs(_rho(obs_vec[name][ph], strength_obs[ph])) for ph in ARC_PHASES]
        slave[name] = float(np.median(rs))

    # ---- surrogate draws (Mode B: cached eigs -> reconstruct W / cophenetic) ----
    dummy = np.random.default_rng(0)
    eigs = {}
    for ph in ARC_PHASES:
        ev, ec = load_or_compute_surrogate_eigs(
            pat, band, ph, W[ph], R, SWAP_FACTOR, SEED, dummy, fc_method="imcoh_abs")
        eigs[ph] = (ev, ec)

    surr_ff = {name: {k: [] for k in FUNCTIONALS} for name in obs_vec}
    for r in range(R):
        # validity: all 5 phases must have a finite eigendecomposition this row
        ok = all(np.isfinite(eigs[ph][0][r]).all() for ph in ARC_PHASES)
        if not ok:
            continue
        Wr = {ph: adjacency_from_laplacian_eigs(eigs[ph][0][r], eigs[ph][1][r])
              for ph in ARC_PHASES}
        for name in obs_vec:
            fn = desc_map.get(name)
            if name == "cophenetic":
                vec = {ph: cophenetic_condensed_from_eigs(eigs[ph][0][r], eigs[ph][1][r])
                       for ph in ARC_PHASES}
            else:
                vec = {ph: fn(Wr[ph]) for ph in ARC_PHASES}
            sff = _sym_functionals(vec)
            for k in FUNCTIONALS:
                surr_ff[name][k].append(sff[k])

    # ---- p-values ----
    rows = []
    for name in obs_vec:
        row = {"patient": pat, "band": band, "descriptor": name,
               "N_nodes": int(N), "slave_r": slave[name]}
        for k in FUNCTIONALS:
            s = np.asarray(surr_ff[name][k], float)
            s = s[np.isfinite(s)]
            o = obs_ff[name][k]
            row[f"{k}_obs"] = o
            row[f"{k}_surr_p50"] = float(np.median(s)) if s.size else np.nan
            row[f"{k}_p"] = float(np.mean(s >= o)) if s.size and np.isfinite(o) else np.nan
            row[f"{k}_nsurr"] = int(s.size)
        rows.append(row)
    return {"rows": rows}


def cohort_gate(df):
    """Per (descriptor, band, functional): one-sided Wilcoxon(obs - surr_p50)."""
    out = []
    for name in RUNGS:
        for band in BANDS:
            x = df[(df.descriptor == name) & (df.band == band)]
            if x.empty:
                continue
            for k in FUNCTIONALS:
                d = x.dropna(subset=[f"{k}_obs", f"{k}_surr_p50"])
                if d.empty:
                    continue
                diff = d[f"{k}_obs"].values - d[f"{k}_surr_p50"].values
                try:
                    _, p = wilcoxon(diff, alternative="greater")
                except Exception:
                    # all-zero differences (e.g. strength: MS preserves it exactly,
                    # so obs == surrogate) -> Wilcoxon undefined -> definitively
                    # non-significant (degenerate). Distinguish from a real NaN.
                    p = 1.0 if np.nanmax(np.abs(diff)) < 1e-6 else float("nan")
                out.append({
                    "descriptor": name, "band": band, "functional": k,
                    "med_obs": float(d[f"{k}_obs"].median()),
                    "med_surr": float(d[f"{k}_surr_p50"].median()),
                    "n_above": int((d[f"{k}_p"] < 0.05).sum()),
                    "n_pat": len(d), "gate_p": float(p),
                })
    return pd.DataFrame(out)


def cross_check(gate):
    """Cophenetic rung must reproduce the flagship verdict (anti-hallucination)."""
    print("\n=== CROSS-CHECK: cophenetic rung vs flagship (audit_150 / audit_152) ===")
    g = gate[gate.descriptor == "cophenetic"]
    tt = g[g.functional == "T_test"].set_index("band")["gate_p"]
    ip = g[g.functional == "T_infspec_pe"].set_index("band")["gate_p"]
    checks = [
        ("T_test alpha CLEAR (~0.024)", tt.get("alpha", np.nan), lambda p: p < 0.05),
        ("T_test beta  CLEAR (~0.032)", tt.get("beta", np.nan), lambda p: p < 0.05),
        ("T_test theta FAIL (~0.78)", tt.get("theta", np.nan), lambda p: p >= 0.20),
        ("T_infspec_pe beta ONLY (~0.0098)", ip.get("beta", np.nan), lambda p: p < 0.05),
    ]
    allok = True
    for label, p, ok in checks:
        status = "PASS" if (np.isfinite(p) and ok(p)) else "WARN"
        allok &= status == "PASS"
        print(f"  [{status}] {label:<34} gate_p={p:.4f}")
    # inference-specificity: beta the ONLY band that clears T_infspec_pe
    ip_clear = ip[ip < 0.05]
    print(f"  T_infspec_pe bands clearing: {list(ip_clear.index)} "
          f"({'beta-only PASS' if list(ip_clear.index) == ['beta'] else 'CHECK'})")
    print(f"  => cophenetic cross-check {'PASS' if allok else 'WARN (inspect before trusting)'}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--betweenness", choices=["off", "on"], default="off",
                    help="include Brandes betweenness rung (slow on dense graphs)")
    ap.add_argument("--limit", type=int, default=0, help="debug: first N cells only")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    descriptors = list(CLASSICAL)
    if args.betweenness == "on":
        descriptors.append(("betweenness", betweenness_centrality))
        RUNGS.insert(-1, "betweenness")

    for pat in COHORT:  # populate rest_pre halves cache (self-contained)
        ensure_half_fcs(pat, BANDS)

    jobs = [(i, p, b, descriptors)
            for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    if args.limit:
        jobs = jobs[:args.limit]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_166] {len(jobs)} cells, {len(RUNGS)} rungs, R={R}, "
          f"Mode B (seed {SEED}), {ncpu} workers")

    t0 = time.time()
    rows, errs = [], []
    with Pool(ncpu) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if "error" in res:
                errs.append(res["error"])
                print(f"[{i}/{len(jobs)}] ERR {res['error']}", flush=True)
                continue
            rows.extend(res["rows"])
            r0 = res["rows"][0]
            cop = next((r for r in res["rows"] if r["descriptor"] == "cophenetic"), r0)
            el = time.time() - t0
            eta = el / i * (len(jobs) - i)
            print(f"[{i}/{len(jobs)}] {r0['patient']}/{r0['band']} "
                  f"coph T_test={cop['T_test_obs']:+.3f}(p={cop['T_test_p']:.3f}) "
                  f"T_ispe={cop['T_infspec_pe_obs']:+.3f}(p={cop['T_infspec_pe_p']:.3f}) "
                  f"| {el:.0f}s ETA {eta:.0f}s", flush=True)

    df = pd.DataFrame(rows).sort_values(["descriptor", "band", "patient"])
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)
    gate = cohort_gate(df)
    gate.to_csv(OUT / "cohort_summary.csv", index=False)

    slav = (df[["descriptor", "slave_r"]].dropna()
            .groupby("descriptor")["slave_r"].median().reset_index())
    slav.to_csv(OUT / "strength_slaving.csv", index=False)

    rt = time.time() - t0
    print(f"\n[audit_166] {len(rows)} rows ({len(errs)} cell errors) in {rt:.0f}s")
    cross_check(gate)
    _print_ladder(gate, slav)
    write_readme(gate, slav, rt, args.betweenness == "on")
    print(f"[audit_166] outputs -> {OUT}")


def _print_ladder(gate, slav):
    slave = dict(zip(slav.descriptor, slav.slave_r))
    print("=== LADDER: whole-task gate p (T_test) + inference gate p (T_infspec_pe) ===")
    hdr = f"{'rung':<18}{'slaveR':>7}  " + "".join(f"{b[:5]:>7}" for b in BANDS) + "   | infspec_pe(beta)"
    print(hdr)
    for name in RUNGS:
        tt = gate[(gate.descriptor == name) & (gate.functional == "T_test")].set_index("band")["gate_p"]
        ip = gate[(gate.descriptor == name) & (gate.functional == "T_infspec_pe")].set_index("band")["gate_p"]
        sr = slave.get(name, np.nan)
        cells = "".join(f"{tt.get(b, np.nan):>7.3f}" for b in BANDS)
        print(f"{name:<18}{sr:>7.2f}  {cells}   | beta={ip.get('beta', np.nan):.4f}")


def write_readme(gate, slav, rt, with_bet):
    slave = dict(zip(slav.descriptor, slav.slave_r))
    L = ["---", "name: pairwise_descriptor_ladder",
         "scope: classical_pairwise_descriptors_vs_multiscale_cophenetic_trace",
         "date: 2026-07-09", "status: current", "---", "",
         "# Pairwise-descriptor ladder (do low-level scalars reproduce the trace?)", "",
         "**Head.** Classical low-level graph descriptors run through the IDENTICAL",
         "rho_sym estimator + matched-strength null as the LRG cophenetic trace.",
         "Mode B: exact flagship surrogate ensemble (seed 20260511), so the",
         "cophenetic rung reproduces audit_150/152. Scored per rung on whole-task",
         "band-selectivity (T_test) and inference-beta-specificity (T_infspec_pe).", "",
         "| rung | strengthR | T_test alpha | T_test beta | T_test theta | infspec_pe beta |",
         "|---|---|---|---|---|---|"]
    for name in RUNGS:
        tt = gate[(gate.descriptor == name) & (gate.functional == "T_test")].set_index("band")["gate_p"]
        ip = gate[(gate.descriptor == name) & (gate.functional == "T_infspec_pe")].set_index("band")["gate_p"]
        sr = slave.get(name, float("nan"))
        L.append(f"| {name} | {sr:.2f} | {tt.get('alpha', float('nan')):.4f} | "
                 f"{tt.get('beta', float('nan')):.4f} | {tt.get('theta', float('nan')):.4f} | "
                 f"{ip.get('beta', float('nan')):.4f} |")
    L += ["", "## Reading",
          "- **C1 recovers alpha / C2 silent on theta / C3 survives MS / C5 inference-beta-only.**",
          "- node **strength** is degenerate under MS by construction (p~1): a strength-",
          "  slaved descriptor cannot clear the null. strengthR = median |Spearman(phi, strength)|.",
          "- **cophenetic** is the anti-hallucination gate: it must clear alpha/beta and",
          "  the beta-only T_infspec_pe (audit_150 / audit_152).", "",
          "## Provenance",
          f"- estimator: rho_sym arc functionals (audit_152 verbatim); positive = trace.",
          f"- null: MODE B, matched-strength seed {SEED}, R={R}, swap {SWAP_FACTOR} (exact",
          "  flagship ensemble; surrogate adjacency = adjacency_from_laplacian_eigs).",
          f"- descriptors: src/lrg_eegfc/utils/metrics/graph_descriptors.py"
          f"{' (+ betweenness)' if with_bet else ''}.",
          f"- cohort gate: Wilcoxon(obs - surr_p50) one-sided greater, n=10.",
          f"- Build: scripts/01_compute/audit/audit_166_pairwise_descriptor_ladder.py",
          f"- runtime {rt:.0f}s."]
    (OUT / "README.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
