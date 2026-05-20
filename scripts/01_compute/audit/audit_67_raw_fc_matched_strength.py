#!/usr/bin/env python3
"""Audit 67 — raw-FC matched-strength surrogate null for ρ_split^raw.

(File numbered as audit_67; the running background process was started
under the prior name audit_66 — its log prefix [audit_66] and the
output directory `data/audit/raw_fc_matched_strength/` are preserved
for that single-shot result. Renamed because audit_66 was independently
used by the parallel-session Grassmann matched-strength surrogate
script. Future invocations use the audit_67 name.)

Apples-to-apples companion to audit_63. Computes the SAME split-baseline
cross-phase coherence statistic on **raw imcoh_abs edges** (no LRG
ultrametric step) under the same R=200 4-cycle ±δ matched-strength
surrogate.

    ρ_split^raw = Spearman(Δ_task^raw, Δ_rest^raw)
                  Δ_task^raw = triu(A^tt) − triu(A^pre_A)
                  Δ_rest^raw = triu(A^post) − triu(A^pre_B)

Decision logic:
- If ρ_split^raw survives matched-strength at α/β with effect size
  comparable to or exceeding the LRG ρ_split (audit_63: α p=0.002 z=+54,
  β p=0.005 z=+52), then LRG adds nothing beyond a coordinate transform.
- If ρ_split^raw dies under matched-strength while LRG ρ_split survives,
  LRG is a real structural amplifier of cross-phase coherence.
- If both survive, LRG adds in a graded sense (compare effect sizes).

Output: data/audit/raw_fc_matched_strength/
            cohort_summary.csv
            per_patient_per_band.csv
            figures/cohort_distribution.pdf
            figures/per_patient_panel.pdf
            README.md

This script is intentionally a near-clone of audit_63 so the comparison
is clean. Differences vs audit_63 are isolated to:
- `lrg_ultrametric_condensed` replaced by `raw_fc_condensed` (squareform).
- Output directory + figure labels.
- Surrogate cells are MUCH faster (no eigendecomposition per surrogate).
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BAND_TEX_DICT, FS_OVERRIDES, nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.workflow.fc import load_fc_matrix

# Reuse the existing half-FC helper (no copy)
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # type: ignore


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["alpha", "beta", "low_gamma"]
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
N_SURROGATES = 200
SWAP_FACTOR = 20

OUT = ROOT / "data" / "audit" / "raw_fc_matched_strength"
FIG = OUT / "figures"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
HALVES_FC_CACHE.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Half-FC cache (lazy populate) — identical to audit_63
# ---------------------------------------------------------------------------
def _half_fc_path(pat: str, band: str, half: str) -> Path:
    return (HALVES_FC_CACHE / pat
            / f"{band}_rest_pre_{half}_imcoh_abs.npy")


def ensure_half_fcs(pat: str, bands: list[str]) -> None:
    missing = [(b, h) for b in bands for h in ("A", "B")
               if not _half_fc_path(pat, b, h).exists()]
    if not missing:
        return
    print(f"[audit_67] {pat}: caching {len(missing)} missing half FCs")
    X = load_timeseries(pat, "rest_pre", SEEG_DATAPATH)
    X = np.asarray(X, dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = FS_OVERRIDES.get(pat, 2048.0)
    nperseg_half = max(256, nperseg_for_fs(fs) // 2)
    halves_fc = compute_imcoh_abs_halves(X, fs, nperseg_half, BRAIN_BANDS)
    del X
    gc.collect()
    pat_dir = HALVES_FC_CACHE / pat
    pat_dir.mkdir(parents=True, exist_ok=True)
    for (band, half), A in halves_fc.items():
        if band not in bands:
            continue
        path = _half_fc_path(pat, band, half)
        np.save(path, np.asarray(A, dtype=np.float32))


def load_phase_fc(pat: str, phase: str, band: str) -> np.ndarray:
    if phase in ("rest_pre_A", "rest_pre_B"):
        half = phase[-1]
        W = np.load(_half_fc_path(pat, band, half))
    else:
        W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
    W = np.asarray(W, dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    W = 0.5 * (W + W.T)
    return W


# ---------------------------------------------------------------------------
# Raw-FC condensed triu (replaces lrg_ultrametric_condensed)
# ---------------------------------------------------------------------------
def raw_fc_condensed(W: np.ndarray) -> np.ndarray:
    """Return the condensed (off-diagonal upper triangle) raw-FC vector.

    No LRG transform. The condensed vector has length N(N-1)/2 and is in
    the canonical scipy `squareform` ordering, matching what
    `lrg_ultrametric_condensed` returns (so downstream ρ_split is on
    the same edge ordering).
    """
    W = W.copy()
    np.fill_diagonal(W, 0.0)
    W = np.maximum(W, W.T)  # ensure exact symmetry for squareform
    return squareform(W, checks=False)


# ---------------------------------------------------------------------------
# Strength-preserving 4-cycle ±δ rewiring — identical to audit_63
# ---------------------------------------------------------------------------
def strength_preserving_shuffle(W: np.ndarray, n_swaps: int,
                                rng: np.random.Generator,
                                w_max: float = 1.0) -> np.ndarray:
    W = W.copy()
    N = W.shape[0]
    samples = rng.integers(0, N, size=(n_swaps, 4))
    fracs = rng.uniform(0.0, 1.0, size=n_swaps)
    for i in range(n_swaps):
        a, b, c, d = samples[i]
        if a == b or a == c or a == d or b == c or b == d or c == d:
            continue
        w1 = W[a, b]; w2 = W[c, d]; w3 = W[a, d]; w4 = W[c, b]
        lo = -w1 if -w1 > -w2 else -w2
        if w3 - w_max > lo:
            lo = w3 - w_max
        if w4 - w_max > lo:
            lo = w4 - w_max
        hi = w_max - w1 if w_max - w1 < w_max - w2 else w_max - w2
        if w3 < hi:
            hi = w3
        if w4 < hi:
            hi = w4
        if lo >= hi:
            continue
        delta = lo + fracs[i] * (hi - lo)
        n1 = w1 + delta; n2 = w2 + delta
        n3 = w3 - delta; n4 = w4 - delta
        W[a, b] = n1; W[b, a] = n1
        W[c, d] = n2; W[d, c] = n2
        W[a, d] = n3; W[d, a] = n3
        W[c, b] = n4; W[b, c] = n4
    return W


def verify_strengths(W_obs: np.ndarray, W_surr: np.ndarray,
                      tol: float = 1e-6) -> bool:
    s_obs = W_obs.sum(axis=1)
    s_surr = W_surr.sum(axis=1)
    return bool(np.max(np.abs(s_obs - s_surr)) < tol)


# ---------------------------------------------------------------------------
# Per-cell pipeline
# ---------------------------------------------------------------------------
def rho_split_from_phases(A_pre_A: np.ndarray, A_pre_B: np.ndarray,
                           A_tt: np.ndarray, A_post: np.ndarray) -> float:
    dA_task = A_tt - A_pre_A
    dA_rest = A_post - A_pre_B
    rho, _ = spearmanr(dA_task, dA_rest)
    return float(rho)


def per_cell(pat: str, band: str, n_surr: int, swap_factor: int,
             rng: np.random.Generator, verbose: bool = False) -> dict | None:
    Ws: dict[str, np.ndarray] = {}
    for phase in PHASES:
        try:
            Ws[phase] = load_phase_fc(pat, phase, band)
        except Exception as e:
            if verbose:
                print(f"[audit_67] SKIP {pat}/{band}/{phase}: {e}")
            return None

    N = Ws[PHASES[0]].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_67] SKIP {pat}/{band}: phase shape mismatch")
        return None

    # Observed (raw-FC condensed, no LRG)
    A_obs = {phase: raw_fc_condensed(Ws[phase]) for phase in PHASES}
    Lc = A_obs[PHASES[0]].size
    rho_obs = rho_split_from_phases(
        A_obs["rest_pre_A"], A_obs["rest_pre_B"],
        A_obs["task_test"], A_obs["rest_post"])

    # Surrogates
    n_swaps = swap_factor * (N * (N - 1)) // 2
    rho_surr = np.empty(n_surr, dtype=float)
    for r in range(n_surr):
        A_surr = {}
        ok = True
        for phase in PHASES:
            W_s = strength_preserving_shuffle(Ws[phase], n_swaps, rng)
            if not verify_strengths(Ws[phase], W_s, tol=1e-4):
                ok = False
                break
            A_surr[phase] = raw_fc_condensed(W_s)
            if A_surr[phase].size != Lc:
                ok = False
                break
        if not ok:
            rho_surr[r] = np.nan
            continue
        rho_surr[r] = rho_split_from_phases(
            A_surr["rest_pre_A"], A_surr["rest_pre_B"],
            A_surr["task_test"], A_surr["rest_post"])

    surr = rho_surr[np.isfinite(rho_surr)]
    if surr.size == 0:
        return None
    surr_mean = float(np.mean(surr))
    surr_std = float(np.std(surr, ddof=1))
    z = (rho_obs - surr_mean) / surr_std if surr_std > 0 else float("nan")
    p_one = float(np.mean(surr >= rho_obs))

    return {
        "patient": pat,
        "band": band,
        "N_nodes": int(N),
        "n_pairs": int(Lc),
        "n_swaps_per_surrogate": int(n_swaps),
        "n_surrogates": int(surr.size),
        "obs_rho": rho_obs,
        "surr_mean_rho": surr_mean,
        "surr_std_rho": surr_std,
        "surr_p5": float(np.quantile(surr, 0.05)),
        "surr_p25": float(np.quantile(surr, 0.25)),
        "surr_p50": float(np.quantile(surr, 0.50)),
        "surr_p75": float(np.quantile(surr, 0.75)),
        "surr_p95": float(np.quantile(surr, 0.95)),
        "surr_p99": float(np.quantile(surr, 0.99)),
        "obs_z": z,
        "obs_p_one_sided": p_one,
        "_surr_array": surr,
    }


# ---------------------------------------------------------------------------
# Cohort summary + figures — identical layout to audit_63 for direct
# visual comparison against the LRG verdict.
# ---------------------------------------------------------------------------
def cohort_summary(rows: list[dict]) -> pd.DataFrame:
    out = []
    for band in TARGET_BANDS:
        sub = [r for r in rows if r["band"] == band]
        if not sub:
            continue
        rho_obs = np.array([r["obs_rho"] for r in sub])
        surr_med_per_pat = np.array([r["surr_p50"] for r in sub])
        n_above = int(sum(1 for r in sub if r["obs_p_one_sided"] < 0.05))
        try:
            wz, wp = wilcoxon(rho_obs - surr_med_per_pat,
                              alternative="greater")
            cohort_z = float(wz)
            cohort_p = float(wp)
        except Exception:
            cohort_z = float("nan")
            cohort_p = float("nan")
        med_obs = float(np.median(rho_obs))
        med_surr_med = float(np.median(surr_med_per_pat))
        p95_surr_med = float(np.quantile(surr_med_per_pat, 0.95))

        if (cohort_z > 2 and abs(med_surr_med) < 0.05
                and n_above >= 8):
            verdict = "separated"
        elif med_surr_med > 0.5 * med_obs and med_obs > 0:
            verdict = "also_positive"
        else:
            verdict = "intermediate"
        out.append({
            "band": band,
            "n_patients": len(sub),
            "obs_median_rho": med_obs,
            "surr_median_rho_median": med_surr_med,
            "surr_median_rho_p95": p95_surr_med,
            "paired_wilcoxon_z": cohort_z,
            "paired_wilcoxon_p": cohort_p,
            "n_above_surrogate": f"{n_above}/{len(sub)}",
            "verdict": verdict,
        })
    return pd.DataFrame(out)


def make_cohort_figure(rows: list[dict], cohort: pd.DataFrame,
                        out_path: Path) -> None:
    fig, axes = plt.subplots(1, len(TARGET_BANDS),
                             figsize=(4.5 * len(TARGET_BANDS), 4.6),
                             sharey=False)
    if len(TARGET_BANDS) == 1:
        axes = [axes]
    for ax, band in zip(axes, TARGET_BANDS):
        sub = [r for r in rows if r["band"] == band]
        if not sub:
            ax.set_axis_off()
            continue
        sub_sorted = sorted(sub, key=lambda r: r["obs_rho"])
        n = len(sub_sorted)
        ys = np.arange(n)
        for y, r in zip(ys, sub_sorted):
            ax.hlines(y, r["surr_p5"], r["surr_p95"],
                      color="#aaaaaa", lw=4, alpha=0.85)
            ax.plot(r["surr_p50"], y, marker="|", color="#444444",
                    markersize=10, mew=1.2)
            ax.plot(r["obs_rho"], y, marker="o", color="#1f3d6e",
                    markersize=7, mec="white", mew=0.8)
        ax.set_yticks(ys)
        ax.set_yticklabels([r["patient"].replace("Pat_", "P")
                            for r in sub_sorted], fontsize=8)
        ax.axvline(0, color="0.6", lw=0.7, ls="--", zorder=0)
        ax.set_xlabel(r"$\rho_{\mathrm{split}}^{\mathrm{raw}}$")
        cohort_row = cohort[cohort.band == band].iloc[0]
        verdict = str(cohort_row["verdict"])
        wp = float(cohort_row["paired_wilcoxon_p"])
        ax.set_title(
            f"{BRAIN_BAND_TEX_DICT[band]} — verdict: {verdict}\n"
            rf"Wilcoxon $p={wp:.4f}$,  $n_{{>}} = {cohort_row['n_above_surrogate']}$",
            fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_per_patient_panel(rows: list[dict], out_path: Path) -> None:
    fig, axes = plt.subplots(len(COHORT), len(TARGET_BANDS),
                             figsize=(3.4 * len(TARGET_BANDS),
                                      1.55 * len(COHORT)),
                             sharex=False)
    for i, pat in enumerate(COHORT):
        for j, band in enumerate(TARGET_BANDS):
            ax = axes[i, j]
            sub = [r for r in rows
                   if r["patient"] == pat and r["band"] == band]
            if not sub:
                ax.set_axis_off()
                continue
            r = sub[0]
            ax.hist(r["_surr_array"], bins=30, color="#cccccc",
                    edgecolor="#888888", alpha=0.85)
            ax.axvline(r["obs_rho"], color="#1f3d6e", lw=1.6,
                       label=f"obs={r['obs_rho']:+.3f}")
            ax.axvline(0, color="#888888", lw=0.6, ls="--")
            if i == len(COHORT) - 1:
                ax.set_xlabel(r"$\rho_{\mathrm{split}}^{\mathrm{raw}}$")
            if j == 0:
                ax.set_ylabel(pat.replace("Pat_", "P"), fontsize=8)
            if i == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
            ax.tick_params(labelsize=7)
            ax.legend(loc="upper left", frameon=False, fontsize=7)
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def write_readme(cohort: pd.DataFrame, n_surr: int, swap_factor: int,
                 runtime_s: float) -> None:
    lines = [
        "---",
        "name: raw_fc_matched_strength",
        "scope: section_5_baseline_lrg_vs_raw_fc_comparison",
        "date: 2026-05-11",
        "status: first_pass",
        "---",
        "",
        "# Raw-FC matched-strength surrogate null for ρ_split^raw",
        "",
        f"**Head.** R={n_surr} strength-preserving (4-cycle ±δ) "
        f"surrogate FCs per (patient, band, phase ∈ {{pre_A, pre_B, tt, "
        f"post}}) for trace bands {{α, β, γ_l}}. Same algorithm as "
        f"audit_63 but the cross-phase coherence statistic is computed "
        f"on **raw imcoh_abs edges** (`triu(A)`) instead of LRG "
        f"ultrametric distances. Direct apples-to-apples comparison: "
        f"does LRG add structural amplification to the cross-phase trace "
        f"signal, or is it a coordinate transform?",
        "",
        "## Cohort verdict",
        "",
        "| band | n | obs median ρ^raw | surr median (per-pat med) | n above own surr | Wilcoxon z | p | verdict |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for _, r in cohort.iterrows():
        lines.append(
            f"| {r['band']} | {r['n_patients']} "
            f"| {r['obs_median_rho']:+.3f} "
            f"| {r['surr_median_rho_median']:+.3f} "
            f"| {r['n_above_surrogate']} "
            f"| {r['paired_wilcoxon_z']:+.2f} "
            f"| {r['paired_wilcoxon_p']:.4f} "
            f"| **{r['verdict']}** |"
        )
    lines.extend([
        "",
        "## Decision rule for the LRG-vs-raw comparison",
        "",
        "Compare side-by-side against audit_63 (LRG ρ_split):",
        "",
        "| outcome | manuscript implication |",
        "|---|---|",
        "| raw dies, LRG survives | LRG is a real structural amplifier of cross-phase coherence; 'LRG deeper than raw FC' claim is defensible. |",
        "| both survive, similar effect | LRG is a coordinate transform; the headline becomes raw FC, LRG is a decorative sanity check. |",
        "| both survive, LRG stronger (higher z, lower p) | LRG amplifies in a graded sense; defensible but quieter claim. |",
        "| raw survives, LRG dies | LRG step is destroying signal; do not headline LRG. |",
        "",
        "## Verdict labels",
        "",
        "Same gate as audit_63 (cohort Wilcoxon z > 2 AND |surr median| < 0.05 "
        "AND n_above ≥ 8/10 → 'separated'). The strict per-patient gate is "
        "preserved for direct comparison.",
        "",
        "## Algorithm caveat (same as audit_63)",
        "",
        "Independent per-phase rewiring decorrelates cross-phase edge "
        "identity by construction. A 'separated' verdict confirms the "
        "observed cross-phase ρ_split^raw is not reproducible from "
        "independent per-phase strength-only randomization — does NOT "
        "isolate shape-driven from strength-driven cross-phase coherence. "
        "But this is the same null applied to the LRG statistic, so the "
        "comparison is fair.",
        "",
        "## Provenance",
        "",
        f"- N_surrogates = {n_surr} per cell",
        f"- n_swaps_per_surrogate = SWAP_FACTOR ({swap_factor}) × N(N-1)/2",
        "- Cohort: " + ", ".join(COHORT),
        "- Bands: " + ", ".join(TARGET_BANDS),
        "- Phases: " + ", ".join(PHASES),
        "- FC method: imcoh_abs",
        "- Statistic: ρ_split^raw = Spearman(Δ_task^raw, Δ_rest^raw) on triu(A)",
        f"- Wall-clock runtime: {runtime_s:.1f} s",
        "- Build script: `scripts/01_compute/audit/audit_67_raw_fc_matched_strength.py`",
        "",
        "## Files",
        "",
        "- `per_patient_per_band.csv` — observed + per-cell surrogate stats",
        "- `cohort_summary.csv` — per-band cohort verdict",
        "- `figures/cohort_distribution.pdf` — per-band cohort overlay",
        "- `figures/per_patient_panel.pdf` — per-(patient, band) histograms",
        "",
    ])
    (OUT / "README.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=TARGET_BANDS)
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--n-surrogates", type=int, default=N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=SWAP_FACTOR)
    ap.add_argument("--seed", type=int, default=20260511)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    bands = args.bands
    patients = args.patients
    n_surr = args.n_surrogates
    swap_factor = args.swap_factor

    print("[audit_67] pre-flight: ensure half FCs cached")
    for pat in patients:
        ensure_half_fcs(pat, bands)

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    rows: list[dict] = []
    for band in bands:
        for pat in patients:
            t_cell = time.time()
            r = per_cell(pat, band, n_surr, swap_factor, rng,
                          verbose=args.verbose)
            if r is None:
                continue
            rows.append(r)
            dt = time.time() - t_cell
            print(f"[audit_67] {pat}/{band}: obs={r['obs_rho']:+.3f} "
                  f"surr_mean={r['surr_mean_rho']:+.3f} z={r['obs_z']:+.2f} "
                  f"p={r['obs_p_one_sided']:.3f} ({dt:.1f}s)")

    runtime = time.time() - t0
    print(f"[audit_67] total: {runtime:.1f}s for {len(rows)} cells")

    per_pat = pd.DataFrame([{k: v for k, v in r.items()
                              if k != "_surr_array"} for r in rows])
    per_pat.to_csv(OUT / "per_patient_per_band.csv", index=False)

    cohort = cohort_summary(rows)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)
    print(cohort.to_string(index=False))

    make_cohort_figure(rows, cohort, FIG / "cohort_distribution.pdf")
    make_per_patient_panel(rows, FIG / "per_patient_panel.pdf")

    write_readme(cohort, n_surr, swap_factor, runtime)
    print(f"[audit_67] outputs at {OUT}")


if __name__ == "__main__":
    main()
