#!/usr/bin/env python3
"""Audit 86 — cluster-extent permutation gate on the WM-EXCLUDED Grassmann trace.

Resolves honesty-flag B on the C6 white-matter-exclusion control. C6-Grassmann
(audit_84) reported per-`k` cohort-Wilcoxon significant-cell COUNTS (e.g.
β 40→40) — the audit_66 basis. But the LOCKED Grassmann gate (C3/C5,
VERDICT_LEDGER Decision 6) is the audit_70 cluster-extent permutation **mass**
test `cluster_p_mass`, NOT the per-`k` count. So "β Grassmann survives WM-X" was
a sensitivity read, not a re-pass of the locked gate under WM exclusion.

This audit runs the audit_70 cluster-extent mass test on the `exclude_wm`
per-`k` T_G curves → `cluster_p_mass^wmX(b)`, judged on the SAME gate as C3
(< 0.05 weak, < 0.01 strong). The cluster-extent null is a phantom-surrogate
permutation, so it needs the per-surrogate T_G(k) stack, which audit_84 did not
persist (only summary percentiles). We regenerate it from the WM-excluded
matched-strength eigvec caches that audit_83/84 built
(`matched_strength_surrogate_wm_excluded_lrg/…_wmX_…npz`).

Critical preamble (per CLAUDE.md rule, before any code)
=======================================================
(1) **Claim.** The audit_66/C3 Grassmann subspace trace (β / γ_l load-bearing)
    survives white-matter exclusion **on the locked gate** — i.e.
    `cluster_p_mass^wmX < 0.05`, not merely "per-`k` counts unchanged".
(2) **Null.** Phantom-surrogate cluster-extent permutation on W_red =
    W[gray, gray]: R=200 matched-strength 4-cycle ±δ surrogates per phase; each
    surrogate r is treated as the phantom observation against the mean of the
    remaining R−1; cluster-mass = Σ −log10(p_k) over k with per-k cohort
    Wilcoxon p < ALPHA_K. Identical machinery to audit_70 (reused, not forked).
(3) **Strongest plausible alternative.** WM exclusion changes N (30–57 % of
    nodes), so the available slow-mode subspace dimension shrinks and the k-axis
    is only approximately comparable to the full-graph C3 gate. A change in
    cluster mass could be subspace-size, not tissue. Mitigated by: the gate is
    self-contained on the WM-excluded graph (obs and null share the same N per
    patient and the same NaN structure for k > N−1), so the p-value is honest
    *within* the WM-excluded geometry; only cross-comparison to the full-graph k
    range is approximate.
(4) **Null's mechanical reach.** Matched-strength fixes per-node strength on
    W_red, so a surviving cluster mass is not a strength artifact. It does NOT
    make the k-axis comparable across N (reported via per-band max k with ≥3
    contributing patients). NaN surrogate rows (rare strength-violation
    failures) are excluded via nanmean — a faithful extension of audit_70 (whose
    full-graph caches were NaN-free).
(5) **Falsification + limitations.** "β Grassmann survives WM-X on the locked
    gate" requires `cluster_p_mass^wmX < 0.05`. If it is ≥ 0.05 while the
    full-graph C3 was < 0.05, the locked gate does NOT survive WM exclusion and
    the C6-Grassmann "unchanged" read must be downgraded. Limitations: (i)
    k-axis not exactly comparable across configs; (ii) single WM definition;
    (iii) δ is expected to fail (it weakened on per-`k` already and is the
    cross-probe epi-biology channel, not the task trace).

Outputs
-------
``data/audit/wm_stratified/``
    grassmann_cluster_extent_wmX.csv     per-band obs mass, cluster_p_mass^wmX,
                                          LOO max, verdict vs C3 gate, full-graph
                                          C3 cross-reference
    grassmann_cluster_extent_wmX_per_k.csv   per-band per-k observed Wilcoxon p
    README_grassmann_cluster_extent.md
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# High-k cells where no patient has enough nodes are all-NaN columns; nanmean
# over them is intentionally NaN (dropped downstream by wilcoxon_per_k_greater).
warnings.filterwarnings("ignore", message="Mean of empty slice")

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.surrogate import load_or_compute_eigs_at_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
# Locked C3 gate machinery (reuse, no fork).
from audit_70_grassmann_cluster_extent import (  # type: ignore
    ALPHA_K,
    chordal_distances_for_k_grid,
    cluster_mass,
    longest_run_below,
    wilcoxon_per_k_greater,
)
# Grassmann + FC helpers (reuse, no fork).
from audit_67_grassmann_epi_exclusion import (  # type: ignore
    K_GRID,
    ensure_half_fcs,
    laplacian_eig,
    load_phase_fc_full,
)
import _wm_stratify as ws  # type: ignore


OUT = ROOT / "data" / "audit" / "wm_stratified"
OUT.mkdir(parents=True, exist_ok=True)
PHASES_3 = ws.PHASES_3  # ("rest_pre_A", "task_test", "rest_post")

# Full-graph C3 gate, for cross-reference (audit_70 output).
C3_CSV = ROOT / "data" / "audit" / "grassmann_cluster_extent" / "cohort_summary.csv"


def _verdict(cluster_p_mass: float) -> str:
    """Identical mapping to audit_70 (C3 gate, locked 2026-05-19)."""
    if not np.isfinite(cluster_p_mass):
        return "undefined"
    if cluster_p_mass < 0.01:
        return "strong"
    if cluster_p_mass < 0.05:
        return "weak"
    return "no_trace"


# ---------------------------------------------------------------------------
# Per-band: obs T_G (P, K) + surrogate T_G (P, R, K) on the WM-excluded graph
# ---------------------------------------------------------------------------
def _build_T_G(band: str, patients: list[str], config: str, k_grid: list[int],
               n_surr: int, swap_factor: int, verbose: bool
               ) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return (obs_T_G [P,K], surr_T_G [P,R,K], used_patients).

    obs eigvecs computed fresh from the WM-excluded submatrix; surrogate eigvecs
    loaded from the warm `wmX` matched-strength caches (computed if missing via
    the shared deterministic cell_rng). NaN where k > N_red−1 for a patient —
    identical NaN structure for obs and surrogates of that patient, so the
    phantom comparison stays fair within the WM-excluded geometry."""
    K = len(k_grid)
    obs_rows, surr_rows, used = [], [], []
    for pat in patients:
        try:
            Ws = {ph: load_phase_fc_full(pat, ph, band) for ph in PHASES_3}
        except Exception as e:
            if verbose:
                print(f"[audit_86] SKIP {pat}/{band}: {e}")
            continue
        N = Ws["rest_pre_A"].shape[0]
        if any(W.shape != (N, N) for W in Ws.values()):
            continue
        wm = ws.wm_mask_for(pat, n_expected=N)
        keep = ws.node_mask_for_config(wm, config)
        if keep is None:
            continue
        idx = np.ix_(keep, keep)
        Wk = {ph: np.ascontiguousarray(Ws[ph][idx]) for ph in PHASES_3}

        # Observed eigvecs (full N_red×N_red, chordal helper skips column 0).
        evec_obs = {ph: laplacian_eig(Wk[ph])[1] for ph in PHASES_3}
        d_pre_tt = chordal_distances_for_k_grid(
            evec_obs["rest_pre_A"], evec_obs["task_test"], k_grid)
        d_tt_post = chordal_distances_for_k_grid(
            evec_obs["task_test"], evec_obs["rest_post"], k_grid)
        obs_rows.append(d_pre_tt - d_tt_post)  # T_G > 0 = trace

        # Surrogate eigvecs from the warm wmX caches.
        evec_surr = {}
        for ph in PHASES_3:
            _, evecs = load_or_compute_eigs_at_path(
                ws.surr_eig_path(config, pat, band, ph, n_surr, swap_factor),
                Wk[ph], n_surr, swap_factor,
                ws.cell_rng(pat, band, ph, config), verbose=verbose)
            evec_surr[ph] = evecs  # (R, N_red, N_red)
        st = np.full((n_surr, K), np.nan)
        for r in range(n_surr):
            ea, eb, ec = (evec_surr["rest_pre_A"][r],
                          evec_surr["task_test"][r],
                          evec_surr["rest_post"][r])
            if not (np.isfinite(ea).all() and np.isfinite(eb).all()
                    and np.isfinite(ec).all()):
                continue
            dpt = chordal_distances_for_k_grid(ea, eb, k_grid)
            dtp = chordal_distances_for_k_grid(eb, ec, k_grid)
            st[r] = dpt - dtp
        surr_rows.append(st)
        used.append(pat)
        del evec_surr
        gc.collect()

    if not used:
        return np.empty((0, K)), np.empty((0, n_surr, K)), used
    return np.vstack(obs_rows), np.stack(surr_rows, axis=0), used


def _cluster_p_mass(obs_T_G: np.ndarray, surr_T_G: np.ndarray
                    ) -> tuple[float, float, float, int, np.ndarray]:
    """audit_70 cluster-extent mass on (P,K) obs vs (P,R,K) surrogates.

    Returns (cluster_p_mass, obs_mass, null_p95_mass, obs_longest_run, obs_p).
    Uses nanmean for the surrogate reference (wmX caches may carry occasional
    NaN surrogate rows; full-graph audit_70 caches were NaN-free)."""
    R = surr_T_G.shape[1]
    surr_mean = np.nanmean(surr_T_G, axis=1)               # (P, K)
    obs_p = wilcoxon_per_k_greater(obs_T_G, surr_mean)
    obs_mass = cluster_mass(obs_p, ALPHA_K)
    obs_LR = longest_run_below(obs_p, ALPHA_K)
    null_mass = np.zeros(R)
    for r in range(R):
        phantom = surr_T_G[:, r, :]
        mask = np.ones(R, dtype=bool)
        mask[r] = False
        ref = np.nanmean(surr_T_G[:, mask, :], axis=1)     # (P, K)
        null_p = wilcoxon_per_k_greater(phantom, ref)
        null_mass[r] = cluster_mass(null_p, ALPHA_K)
    cp_mass = (1 + int(np.sum(null_mass >= obs_mass))) / (R + 1)
    return (float(cp_mass), float(obs_mass),
            float(np.percentile(null_mass, 95)), int(obs_LR), obs_p)


def _loo_max(obs_T_G: np.ndarray, surr_T_G: np.ndarray, used: list[str]
             ) -> tuple[float, str]:
    """Worst-case leave-one-out cluster_p_mass (descriptive, never a gate)."""
    P = obs_T_G.shape[0]
    if P < 4:
        return float("nan"), "—"
    worst_p, worst_pat = -1.0, "—"
    for i in range(P):
        keep = np.ones(P, dtype=bool)
        keep[i] = False
        cp, *_ = _cluster_p_mass(obs_T_G[keep], surr_T_G[keep])
        if cp > worst_p:
            worst_p, worst_pat = cp, used[i]
    return float(worst_p), worst_pat


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(summary: pd.DataFrame, config: str, k_grid: list[int],
                 n_surr: int, swap_factor: int, runtime_s: float) -> None:
    L: list[str] = []
    a = L.append
    a("---")
    a("name: wm_grassmann_cluster_extent")
    a("scope: cluster_extent_mass_gate_on_wm_excluded_grassmann_trace")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: first_pass")
    a("build_script: scripts/01_compute/audit/audit_86_wm_grassmann_cluster_extent.py")
    a("gate: same as audit_70 / C3 (cluster_p_mass < 0.01 strong, < 0.05 weak)")
    a("---")
    a("")
    a("# WM-excluded Grassmann on the locked C3 cluster-extent gate (n=10)")
    a("")
    a("**Head.** Re-runs the locked audit_70 cluster-extent permutation **mass** "
      f"gate on the `{config}` per-`k` Grassmann T_G(k) curves, so the "
      "C6-Grassmann sensitivity is judged on the SAME gate as the locked C3 "
      "Grassmann verdict — not the lighter per-`k` Wilcoxon count audit_84 "
      "reported. `cluster_p_mass^wmX < 0.05` ⇒ the subspace trace survives WM "
      "exclusion on the locked gate.")
    a("")
    a("## Per-band verdict (WM-excluded, on the C3 gate)")
    a("")
    a("| band | obs mass | null p95 mass | obs run | cluster_p_mass^wmX "
      "| LOO max p (drop) | verdict^wmX | C3 full verdict | C3 full p |")
    a("|---|---|---|---|---|---|---|---|---|")
    for _, r in summary.iterrows():
        a(f"| {BRAIN_BAND_TEX_DICT.get(r['band'], r['band'])} "
          f"| {r['obs_cluster_mass']:.2f} | {r['null_p95_mass']:.2f} "
          f"| {int(r['obs_longest_run'])} | {r['cluster_p_mass_wmX']:.4f} "
          f"| {r['loo_max_p_mass']:.4f} ({r['loo_argmax_patient']}) "
          f"| **{r['verdict_wmX']}** | {r['c3_full_verdict']} "
          f"| {r['c3_full_cluster_p_mass']} |")
    a("")
    a("## Gate (identical to audit_70 / C3, locked 2026-05-19)")
    a("- `cluster_p_mass < 0.01` → strong")
    a("- `0.01 ≤ cluster_p_mass < 0.05` → weak")
    a("- `≥ 0.05` → no_trace")
    a("")
    a("## Reading")
    a("- `verdict^wmX` is the WM-excluded Grassmann verdict ON THE LOCKED GATE. "
      "Compare to `C3 full verdict` (the locked full-graph C3 Grassmann "
      "verdict). `persist on gate` = both < 0.05; the C6-Grassmann 'unchanged' "
      "read is upgraded from a per-`k` count to a re-pass of the locked gate.")
    a("- δ is expected to fail (it weakened on per-`k` already; it is the "
      "cross-probe epi-biology channel, not the task trace).")
    a("")
    a("## Caveats")
    a(f"- {config} changes N (30–57 % of nodes are WM) → the k-axis is only "
      "approximately comparable to the full-graph C3 gate; the p-value is honest "
      "*within* the WM-excluded geometry (obs and null share N and NaN "
      "structure per patient).")
    a("- NaN surrogate rows excluded via nanmean (faithful extension of "
      "audit_70, whose full-graph caches were NaN-free).")
    a("- Single WM definition (atlas-dominant `Wm`).")
    a("")
    a("## Provenance")
    a(f"- config = {config}; R = {n_surr}, SWAP_FACTOR = {swap_factor}; "
      f"k_grid = {k_grid[0]}..{k_grid[-1]} (n={len(k_grid)}, self-limits per "
      "patient via NaN).")
    a("- Surrogate eigvecs: warm `matched_strength_surrogate_wm_excluded_lrg/` "
      "caches (seed 20260608, built by audit_83/84).")
    a("- Phases: " + ", ".join(PHASES_3) + " (split-baseline rsPre leg = "
      "rest_pre_A).")
    a("- Cohort: " + ", ".join(ws.COHORT))
    a("- WM = dominant Desikan-Killiany tissue == 'Wm' (atlas argmax).")
    a(f"- Wall-clock: {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `grassmann_cluster_extent_wmX.csv` — per-band gate verdict + C3 xref")
    a("- `grassmann_cluster_extent_wmX_per_k.csv` — per-band per-k observed p")
    (OUT / "README_grassmann_cluster_extent.md").write_text(
        "\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(ws.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(ws.COHORT))
    ap.add_argument("--config", default="exclude_wm",
                    choices=["exclude_wm", "wm_only"])
    ap.add_argument("--n-surrogates", type=int, default=ws.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=ws.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    # C3 full-graph cross-reference.
    c3 = pd.read_csv(C3_CSV) if C3_CSV.exists() else pd.DataFrame()

    print("[audit_86] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    summary_rows, per_k_rows = [], []
    for band in args.bands:
        tb = time.time()
        obs_T_G, surr_T_G, used = _build_T_G(
            band, args.patients, args.config, K_GRID,
            args.n_surrogates, args.swap_factor, args.verbose)
        if not used:
            print(f"[audit_86] {band}: no usable patients")
            continue
        cp_mass, obs_mass, null_p95, obs_LR, obs_p = _cluster_p_mass(
            obs_T_G, surr_T_G)
        loo_p, loo_pat = _loo_max(obs_T_G, surr_T_G, used)

        c3_row = c3[c3.band == band] if not c3.empty else pd.DataFrame()
        c3_v = (str(c3_row.iloc[0]["verdict_cluster_extent"])
                if not c3_row.empty else "missing")
        c3_p = (f"{float(c3_row.iloc[0]['cluster_p_cluster_mass']):.4f}"
                if not c3_row.empty else "—")

        summary_rows.append({
            "band": band, "config": args.config, "n_patients": len(used),
            "obs_cluster_mass": obs_mass, "null_p95_mass": null_p95,
            "obs_longest_run": obs_LR,
            "cluster_p_mass_wmX": cp_mass,
            "verdict_wmX": _verdict(cp_mass),
            "loo_max_p_mass": loo_p, "loo_argmax_patient": loo_pat,
            "c3_full_verdict": c3_v, "c3_full_cluster_p_mass": c3_p,
            "max_k_ge3_patients": int(
                max((k for ki, k in enumerate(K_GRID)
                     if np.isfinite(obs_T_G[:, ki]).sum() >= 3), default=0)),
        })
        for ki, k in enumerate(K_GRID):
            per_k_rows.append({"band": band, "k": k,
                               "obs_p_one_sided_greater": float(obs_p[ki]),
                               "n_patients_finite": int(
                                   np.isfinite(obs_T_G[:, ki]).sum())})
        print(f"[audit_86] {band}: cp_mass^wmX={cp_mass:.4f} "
              f"({_verdict(cp_mass)}) obs_mass={obs_mass:.2f} "
              f"LOO_max={loo_p:.4f}({loo_pat}) C3={c3_v} "
              f"({time.time()-tb:.1f}s)")

    runtime = time.time() - t0
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "grassmann_cluster_extent_wmX.csv", index=False)
    pd.DataFrame(per_k_rows).to_csv(
        OUT / "grassmann_cluster_extent_wmX_per_k.csv", index=False)
    write_readme(summary, args.config, K_GRID, args.n_surrogates,
                 args.swap_factor, runtime)
    print(f"[audit_86] {len(summary)} bands, {runtime:.1f}s -> {OUT}")
    if not summary.empty:
        print(summary[["band", "cluster_p_mass_wmX", "verdict_wmX",
                       "c3_full_verdict", "loo_max_p_mass"]].to_string(
            index=False))


if __name__ == "__main__":
    main()
