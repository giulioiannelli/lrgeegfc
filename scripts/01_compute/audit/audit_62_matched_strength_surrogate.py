#!/usr/bin/env python3
"""Audit 62 — matched-strength Laplacian surrogate null for §5.3 ρ.

C5 of the 2026-05-10 coding-agent batch (file 15 §4 of the verification
package). For each (patient, band) cell with band ∈ {α, β, γ_l},
generate R=100 strength-preserving surrogate FC matrices per phase,
run the LRG → ultrametric pipeline on each surrogate, and compute the
surrogate ρ_3phase = Spearman(Δ_task, Δ_rest)[upper triangular] under
the SHARED-baseline construction (see "Test choice" below).

Strength-preserving randomization (4-cycle ±δ perturbation)
----------------------------------------------------------
For a fully-connected weighted graph the canonical max-entropy
randomization preserving per-node strengths exactly is the 4-cycle
weight-shift (Squartini et al. 2011, MacMahon & Garlaschelli 2015):
pick four distinct nodes a, b, c, d, take w1=W[a,b], w2=W[c,d],
w3=W[a,d], w4=W[c,b], pick δ uniformly in
[max(-w1,-w2,w3-1,w4-1), min(1-w1,1-w2,w3,w4)] (here weights are in
[0, 1] for |ImCoh|), and set w1+=δ, w2+=δ, w3-=δ, w4-=δ. Per-node
strengths s_a, s_b, s_c, s_d are exactly preserved.

Test choice
-----------
The §5.3 split-baseline ρ uses D_pre_A and D_pre_B (LRG ultrametric
matrices from disjoint halves of the rest_pre time series). Building
matched-strength surrogates of the half-baseline FCs would require
re-cutting and re-FC'ing the source time series for every realization.
This script tests the SHARED-baseline ρ instead:

    Δ_task = D_test - D_pre,   Δ_rest = D_post - D_pre
    ρ_shared = Spearman(Δ_task, Δ_rest)

The shared-baseline ρ is *artifactually inflated* relative to the
split-baseline ρ (numerators share D_pre), but the inflation acts on
both the observed and the surrogate, so the **observed-vs-surrogate
contrast** still tests "is the trace direction shape-driven (rejected
by surrogate) or could it be reproduced by matched-strength noise
(not rejected)?". This is exactly the C5 question; the §5.3 framing
is unaffected by the choice of baseline as long as both observed and
surrogate use the same one.

Outputs
-------
    data/audit/matched_strength_surrogate/per_patient_per_band.csv
    data/audit/matched_strength_surrogate/cohort_summary.csv
    data/audit/matched_strength_surrogate/figures/cohort_distribution.pdf
    data/audit/matched_strength_surrogate/figures/per_patient_panel.pdf
    data/audit/matched_strength_surrogate/README.md

Surrogate cache reuse (post-2026-05-11)
---------------------------------------
This script regenerates surrogate adjacency matrices on every invocation
because it predates the disk cache introduced in audit_66. New audits
that need the same matched-strength ensemble MUST instead import from
``lrg_eegfc.utils.surrogate.matched_strength`` and call
``load_or_compute_surrogate_eigs(pat, band, phase, W, R, SF, seed, rng)``,
which returns ``(eigvals[R, N], eigvecs[R, N, N])`` cached at
``data/cache/matched_strength_surrogate_lrg/Pat_NN/``
``{band}_{phase}_R{R}_swap{SF}_seed{S}_imcoh_abs.npz``.
Any downstream statistic (KC, Grassmann, ρ propagator, eigenmode embedding)
on the SAME (R, swap_factor, seed) ensemble runs in minutes from cache
instead of hours from regeneration. The canonical R=200, SWAP_FACTOR=20
ensemble for trace bands {α, β, γ_l} is keyed on seed=20260511. This
script's R=100, SWAP_FACTOR=5 ensemble is NOT cached; cross-comparison
with the canonical ensemble requires re-running with the canonical
parameters or accepting independent realizations.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.fc import load_fc_matrix


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
TARGET_BANDS = ["alpha", "beta", "low_gamma"]
PHASES = ("rest_pre", "task_test", "rest_post")
FC_METHOD = "imcoh_abs"
N_SURROGATES = 100
SWAP_FACTOR = 5  # n_swaps_per_surrogate = SWAP_FACTOR * N(N-1)/2

OUT = ROOT / "data" / "audit" / "matched_strength_surrogate"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# LRG primitives — bypass entropy() loop (we don't need the entropy curve)
# ---------------------------------------------------------------------------
def lrg_ultrametric_condensed(W: np.ndarray) -> np.ndarray:
    """Run LRG on adjacency W → return condensed cophenet (ultrametric).

    Mirrors `compute_lrg_analysis` minus the 400-step entropy loop and
    minus the networkx wrapper (ImCoh adjacency is fully connected by
    construction; no need to extract a giant component).
    """
    deg = W.sum(axis=1)
    L = np.diag(deg) - W
    eigvals, eigvecs = np.linalg.eigh(L)
    lam_max = eigvals[-1]  # eigh returns ascending
    tau = 1.0 / lam_max
    diag_exp = np.exp(-tau * eigvals)
    rho = (eigvecs * diag_exp) @ eigvecs.T
    rho /= np.trace(rho)
    with np.errstate(divide="ignore"):
        Trho = 1.0 / rho
    Trho = np.maximum(Trho, Trho.T)
    np.fill_diagonal(Trho, 0.0)
    finite = np.isfinite(Trho)
    if not finite.all():
        cap = np.nanmax(Trho[finite]) if finite.any() else 1e6
        Trho = np.where(finite, Trho, cap)
    Trho_condensed = squareform(Trho, checks=False)
    Z = linkage(Trho_condensed, method="average")
    return cophenet(Z)


# ---------------------------------------------------------------------------
# Strength-preserving randomization — 4-cycle ±δ perturbation
# ---------------------------------------------------------------------------
def strength_preserving_shuffle(W: np.ndarray, n_swaps: int,
                                rng: np.random.Generator,
                                w_max: float = 1.0) -> np.ndarray:
    """Apply n_swaps 4-cycle ±δ perturbations to W (in-place on a copy).

    For each swap: pick 4 distinct nodes (a,b,c,d), perturb the four
    edges (a,b),(c,d),(a,d),(c,b) by (+δ,+δ,-δ,-δ) for δ uniform in
    its feasible range. Per-node strengths are preserved exactly.

    Vectorized: random nodes + δ-fractions are pre-batched in two
    numpy calls to avoid per-iteration rng overhead. Duplicates among
    the 4 nodes are detected inline and the swap is skipped (cheap to
    over-sample n_swaps slightly to compensate; for N>>4 the duplicate
    rate is < 5%).
    """
    W = W.copy()
    N = W.shape[0]
    samples = rng.integers(0, N, size=(n_swaps, 4))
    fracs = rng.uniform(0.0, 1.0, size=n_swaps)
    accepted = 0
    for i in range(n_swaps):
        a, b, c, d = samples[i]
        # Skip rows with duplicate nodes (rare for N >> 4)
        if a == b or a == c or a == d or b == c or b == d or c == d:
            continue
        w1 = W[a, b]
        w2 = W[c, d]
        w3 = W[a, d]
        w4 = W[c, b]
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
        new1 = w1 + delta
        new2 = w2 + delta
        new3 = w3 - delta
        new4 = w4 - delta
        W[a, b] = new1; W[b, a] = new1
        W[c, d] = new2; W[d, c] = new2
        W[a, d] = new3; W[d, a] = new3
        W[c, b] = new4; W[b, c] = new4
        accepted += 1
    return W


# ---------------------------------------------------------------------------
# Per-cell pipeline
# ---------------------------------------------------------------------------
def compute_rho_shared(D_pre_cond: np.ndarray, D_test_cond: np.ndarray,
                        D_post_cond: np.ndarray) -> float:
    dD_task = D_test_cond - D_pre_cond
    dD_rest = D_post_cond - D_pre_cond
    rho, _ = spearmanr(dD_task, dD_rest)
    return float(rho)


def per_cell(pat: str, band: str, n_surr: int,
             rng: np.random.Generator) -> dict | None:
    """Observed + surrogate ρ_shared for one (patient, band) cell."""
    Ws = {}
    for phase in PHASES:
        try:
            W = load_fc_matrix(pat, phase, band, fc_method=FC_METHOD)
        except Exception as e:
            print(f"[audit_62] SKIP {pat}/{band}/{phase}: {e}")
            return None
        if W is None:
            return None
        W = np.asarray(W, dtype=float)
        np.fill_diagonal(W, 0.0)
        # Clip into [0, 1] for safety (imcoh_abs is in [0, 1] but
        # numerical drift can push slightly outside)
        W = np.clip(W, 0.0, 1.0)
        Ws[phase] = W

    # Observed
    D_obs = {}
    for phase in PHASES:
        D_obs[phase] = lrg_ultrametric_condensed(Ws[phase])
    # Confirm condensed lengths match
    Lc = D_obs[PHASES[0]].size
    if any(D.size != Lc for D in D_obs.values()):
        print(f"[audit_62] SKIP {pat}/{band}: ultrametric size mismatch "
              f"(perhaps disconnected giant components per phase)")
        return None
    rho_obs = compute_rho_shared(D_obs["rest_pre"], D_obs["task_test"],
                                 D_obs["rest_post"])

    # Surrogates — per-phase independent randomization
    N = Ws["rest_pre"].shape[0]
    n_swaps = SWAP_FACTOR * (N * (N - 1)) // 2
    rho_surr = np.empty(n_surr, dtype=float)
    for r in range(n_surr):
        D_surr = {}
        ok = True
        for phase in PHASES:
            W_surr = strength_preserving_shuffle(Ws[phase], n_swaps, rng)
            d = lrg_ultrametric_condensed(W_surr)
            if d.size != Lc:
                ok = False
                break
            D_surr[phase] = d
        if not ok:
            rho_surr[r] = np.nan
            continue
        rho_surr[r] = compute_rho_shared(
            D_surr["rest_pre"], D_surr["task_test"], D_surr["rest_post"])

    surr_finite = rho_surr[np.isfinite(rho_surr)]
    if surr_finite.size == 0:
        return None
    surr_mean = float(np.mean(surr_finite))
    surr_sd = float(np.std(surr_finite, ddof=1))
    z_obs = (rho_obs - surr_mean) / surr_sd if surr_sd > 0 else float("nan")
    p_one = float(np.mean(surr_finite >= rho_obs))  # one-sided

    return {
        "patient": pat,
        "band": band,
        "N_nodes": int(N),
        "n_pairs": int(Lc),
        "n_swaps_per_surrogate": int(n_swaps),
        "n_surrogates": int(surr_finite.size),
        "rho_observed_shared": rho_obs,
        "surr_mean": surr_mean,
        "surr_sd": surr_sd,
        "surr_q05": float(np.quantile(surr_finite, 0.05)),
        "surr_q95": float(np.quantile(surr_finite, 0.95)),
        "z_observed_vs_surr": z_obs,
        "p_one_sided_obs_geq_surr": p_one,
        "_surr_array": surr_finite,  # used for plotting; stripped before CSV
    }


# ---------------------------------------------------------------------------
# Cohort + figures
# ---------------------------------------------------------------------------
def cohort_summary(rows: list[dict]) -> pd.DataFrame:
    out = []
    for band in TARGET_BANDS:
        sub = [r for r in rows if r["band"] == band]
        if not sub:
            continue
        rho_obs = np.array([r["rho_observed_shared"] for r in sub])
        rho_surr_mean = np.array([r["surr_mean"] for r in sub])
        n_pos = int((rho_obs > 0).sum())
        n_above_surr = int((rho_obs > rho_surr_mean).sum())
        # Cohort z-score: mean obs vs cohort-aggregated surrogate spread
        all_surr = np.concatenate([r["_surr_array"] for r in sub])
        cohort_z = ((np.mean(rho_obs) - np.mean(all_surr))
                    / (np.std(all_surr, ddof=1) if np.std(all_surr) > 0
                       else float("nan")))
        # Verdict
        med_obs = float(np.median(rho_obs))
        med_surr = float(np.median(all_surr))
        if abs(med_surr) < 0.05 and cohort_z > 2.0:
            verdict = "separated"
        elif med_surr > 0.05 and abs(med_obs - med_surr) < 0.05:
            verdict = "also_positive"
        else:
            verdict = "intermediate"
        out.append({
            "band": band,
            "n_patients": len(sub),
            "median_rho_obs": med_obs,
            "median_rho_surr": med_surr,
            "n_pos_obs": n_pos,
            "n_obs_above_surr_mean": n_above_surr,
            "cohort_z_obs_vs_surr": float(cohort_z),
            "verdict": verdict,
        })
    return pd.DataFrame(out)


def make_cohort_figure(rows: list[dict], out_path: Path) -> None:
    fig, axes = plt.subplots(1, len(TARGET_BANDS),
                             figsize=(4.2 * len(TARGET_BANDS), 4.4),
                             sharey=False)
    if len(TARGET_BANDS) == 1:
        axes = [axes]
    for ax, band in zip(axes, TARGET_BANDS):
        sub = [r for r in rows if r["band"] == band]
        if not sub:
            ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} — no data",
                         fontsize=10)
            continue
        all_surr = np.concatenate([r["_surr_array"] for r in sub])
        rho_obs = np.array([r["rho_observed_shared"] for r in sub])
        ax.hist(all_surr, bins=40, color="#cccccc", edgecolor="#888888",
                alpha=0.85, label="surrogate (cohort pooled)")
        for r in sub:
            ax.axvline(r["rho_observed_shared"], color="#1f3d6e",
                       lw=0.8, alpha=0.55)
        ax.axvline(np.median(rho_obs), color="#d6603a",
                   lw=2.0, label=f"median obs = {np.median(rho_obs):+.3f}")
        ax.axvline(np.mean(all_surr), color="#444444", ls="--", lw=1.0,
                   label=f"mean surr = {np.mean(all_surr):+.3f}")
        ax.set_xlabel(r"$\rho_{\mathrm{shared}}$")
        ax.set_ylabel("count")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=11)
        ax.legend(loc="upper left", frameon=False, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def make_per_patient_panel(rows: list[dict], out_path: Path) -> None:
    fig, axes = plt.subplots(len(COHORT), len(TARGET_BANDS),
                             figsize=(3.4 * len(TARGET_BANDS),
                                      1.5 * len(COHORT)),
                             sharex=False)
    for i, pat in enumerate(COHORT):
        for j, band in enumerate(TARGET_BANDS):
            ax = axes[i, j] if len(COHORT) > 1 else axes[j]
            sub = [r for r in rows
                   if r["patient"] == pat and r["band"] == band]
            if not sub:
                ax.set_axis_off()
                continue
            r = sub[0]
            ax.hist(r["_surr_array"], bins=30, color="#cccccc",
                    edgecolor="#888888", alpha=0.85)
            ax.axvline(r["rho_observed_shared"], color="#1f3d6e", lw=1.6,
                       label=f"obs={r['rho_observed_shared']:+.3f}")
            ax.set_xlabel(r"$\rho_{\mathrm{shared}}$" if i == len(COHORT) - 1
                          else "")
            if j == 0:
                ax.set_ylabel(pat.replace("Pat_", "P"), fontsize=8)
            if i == 0:
                ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=10)
            ax.tick_params(labelsize=7)
            ax.legend(loc="upper left", frameon=False, fontsize=7)
            ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def write_readme(cohort: pd.DataFrame, n_surr: int, n_swap_factor: int,
                  per_pat: pd.DataFrame, runtime_s: float) -> None:
    lines = [
        "---",
        "name: matched_strength_surrogate",
        "scope: section_5_3_lrg_rho_robustness",
        "date: 2026-05-10",
        "status: first_pass",
        "---",
        "",
        "# Matched-strength Laplacian surrogate null for §5.3 ρ",
        "",
        "**Head.** R=" + str(n_surr) + " strength-preserving (4-cycle "
        "±δ) surrogate Laplacians per (patient, band, phase) for the three "
        "trace bands {α, β, γ_l}, fed through the LRG → cophenet pipeline "
        "and re-correlated under the SHARED-baseline ρ (Δ_task vs Δ_rest "
        "with a common D_pre). Verdict per band lives in `cohort_summary.csv`.",
        "",
        "## Cohort verdict",
        "",
        "| band | n | median obs ρ | median surr ρ | n obs > surr mean | cohort z | verdict |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, r in cohort.iterrows():
        lines.append(
            f"| {r['band']} | {r['n_patients']} | "
            f"{r['median_rho_obs']:+.3f} | {r['median_rho_surr']:+.3f} "
            f"| {r['n_obs_above_surr_mean']}/{r['n_patients']} "
            f"| {r['cohort_z_obs_vs_surr']:+.2f} | **{r['verdict']}** |"
        )
    lines.extend([
        "",
        "## Verdict labels",
        "",
        "- **separated**: surrogate cohort pool centered near zero "
        "(|median_surr| < 0.05) AND cohort z > 2 → the trace direction "
        "is shape-driven; matched-strength noise does not reproduce it.",
        "- **also_positive**: surrogate cohort pool also positive AND "
        "observed ≈ surrogate median → the trace direction is reproducible "
        "by matched-strength noise; the test is sensitive to edge-strength "
        "heterogeneity rather than to specific edge-identity structure.",
        "- **intermediate**: anything else (typically: surrogate centered "
        "near zero but observed only modestly above; or partial overlap).",
        "",
        "## What this rules out and what it does not",
        "",
        "Rejecting the matched-strength null says: the per-pair LRG "
        "ultrametric pattern that drives ρ is NOT recoverable from the "
        "edge-strength sequence alone. Specific edge-identity structure "
        "carries the signal. This is a stricter null than the §5.3 "
        "within-baseline split-half null (which permutes time-window "
        "halves, preserving full network structure).",
        "",
        "Failing to reject (verdict 'also_positive') would mean the "
        "trace-direction observation is consistent with an artifact of "
        "the per-node strength sequence under the LRG → cophenet "
        "pipeline. The §5.3 framing would need a strength-explained "
        "caveat.",
        "",
        "## Test choice (see script docstring)",
        "",
        "Surrogate uses SHARED-baseline ρ rather than §5.3 SPLIT-baseline "
        "ρ for tractability (avoiding re-FC of half-time-series for every "
        "surrogate). Both observed and surrogate use the same baseline "
        "construction, so the comparison is fair. The shared-baseline ρ "
        "is artifactually inflated (numerator shares D_pre) but the "
        "inflation acts on both arms.",
        "",
        "## Provenance",
        "",
        f"- N_surrogates = {n_surr} per cell",
        f"- n_swaps_per_surrogate = SWAP_FACTOR ({n_swap_factor}) × N(N-1)/2",
        "- Cohort: " + ", ".join(COHORT),
        "- Bands: " + ", ".join(TARGET_BANDS),
        "- Phases: rest_pre, task_test, rest_post (shared-baseline)",
        "- FC method: imcoh_abs",
        f"- Wall-clock runtime: {runtime_s:.1f} s",
        "- Build script: `scripts/01_compute/audit/audit_62_matched_strength_surrogate.py`",
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
    ap.add_argument("--seed", type=int, default=20260510)
    args = ap.parse_args()

    bands = args.bands
    patients = args.patients
    n_surr = args.n_surrogates

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    rows: list[dict] = []
    for band in bands:
        for pat in patients:
            t_cell = time.time()
            r = per_cell(pat, band, n_surr, rng)
            if r is None:
                continue
            rows.append(r)
            dt = time.time() - t_cell
            print(f"[audit_62] {pat}/{band}: rho_obs={r['rho_observed_shared']:+.3f} "
                  f"surr_mean={r['surr_mean']:+.3f} z={r['z_observed_vs_surr']:+.2f} "
                  f"({dt:.1f}s)")

    runtime = time.time() - t0
    print(f"[audit_62] total: {runtime:.1f}s for {len(rows)} cells")

    # Per-cell CSV (strip surr arrays)
    per_pat = pd.DataFrame([{k: v for k, v in r.items()
                              if k != "_surr_array"} for r in rows])
    per_pat.to_csv(OUT / "per_patient_per_band.csv", index=False)

    # Cohort
    cohort = cohort_summary(rows)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)
    print(cohort.to_string(index=False))

    # Figures
    make_cohort_figure(rows, FIG / "cohort_distribution.pdf")
    make_per_patient_panel(rows, FIG / "per_patient_panel.pdf")

    write_readme(cohort, n_surr, SWAP_FACTOR, per_pat, runtime)
    print(f"[audit_62] outputs at {OUT}")


if __name__ == "__main__":
    main()
