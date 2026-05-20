#!/usr/bin/env python3
"""Audit 68 — α epi-exclusion sensitivity (matched-strength split-baseline).

The 2026-05-14 evidence review (§4 item 2) flagged that the per-patient
α trace correlates with implant epileptic-zone fraction (Spearman
ρ(frac_epi, α obs_ρ) = -0.41), suggesting α reorganization lives in
clinically healthy circuits. This script re-runs the audit_63
split-baseline matched-strength surrogate test at α band ONLY, after
removing epileptic-zone contacts from each patient's FC matrix.

Predicate: if the α trace is genuinely in non-epi circuits, removing
epi nodes should NOT weaken the cohort signal (and may sharpen it).
If the α trace is partly carried by epi-zone activity, removing epi
nodes should weaken the cohort signal.

Inputs (reused from audit_63):
    - cached half FCs at data/cache/imcoh_halves_fc/Pat_NN/...
    - load_fc_matrix(pat, phase, "alpha", "imcoh_abs") for task_test, rest_post
    - load_epileptic_nodes(pat) + channel_labels.csv for epi masking

Outputs:
    data/audit/alpha_epi_exclusion/
        per_patient.csv         — per-patient observed + surrogate stats
        cohort_summary.csv      — cohort verdict
        comparison.csv          — vs audit_63 α (full FC, no exclusion)
        README.md               — head + verdict
        figures/cohort.pdf      — overlay of per-patient observed vs surrogate
        figures/comparison.pdf  — full vs excluded ρ_split scatter
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_epileptic_nodes
from lrg_eegfc.workflow.fc import load_fc_matrix

# Reuse audit_63 helpers — strength-preserving rewiring + LRG ultrametric.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    HALVES_FC_CACHE,
    ensure_half_fcs,
    lrg_ultrametric_condensed,
    rho_split_from_phases,
    strength_preserving_shuffle,
    verify_strengths,
    _half_fc_path,
)


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BAND = "alpha"
PHASES = ("rest_pre_A", "rest_pre_B", "task_test", "rest_post")
N_SURROGATES = 200
SWAP_FACTOR = 20

OUT = ROOT / "data" / "audit" / "alpha_epi_exclusion"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)


def load_channels(pat: str) -> list[str]:
    """FC-index-aligned channel labels (matches load_ch in scripts/08_epileptic)."""
    p = SEEG_DATAPATH / pat / "channel_labels.csv"
    if not p.exists():
        p = SEEG_DATAPATH / pat / "channel_labels.txt"
    with open(p) as f:
        first = f.readline().strip()
    skip = 1 if first.lower() == "label" or first.strip('"').lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [
        str(l).strip('"').split(",")[0].strip().replace(" ", "")
        for l in df.iloc[:, 0]
    ]


def epi_mask_for(pat: str) -> np.ndarray:
    ch = load_channels(pat)
    epi = set(load_epileptic_nodes(pat))
    return np.array([l in epi for l in ch], dtype=bool)


def load_phase_fc(pat: str, phase: str, band: str) -> np.ndarray:
    """Load a phase FC matrix (full or half), normalize to (N,N) float64."""
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


def per_cell(pat: str, n_surr: int, swap_factor: int,
             rng: np.random.Generator, verbose: bool = False) -> dict | None:
    mask = epi_mask_for(pat)
    keep = ~mask
    N_full = mask.size
    N_keep = int(keep.sum())
    n_epi = int(mask.sum())

    Ws: dict[str, np.ndarray] = {}
    for phase in PHASES:
        try:
            W = load_phase_fc(pat, phase, BAND)
        except Exception as e:
            if verbose:
                print(f"[audit_68] SKIP {pat}/{phase}: {e}")
            return None
        if W.shape[0] != N_full:
            if verbose:
                print(f"[audit_68] SKIP {pat}/{phase}: "
                      f"FC shape {W.shape} vs mask {N_full}")
            return None
        Ws[phase] = W[keep][:, keep]

    # Observed
    D_obs = {phase: lrg_ultrametric_condensed(Ws[phase]) for phase in PHASES}
    Lc = D_obs[PHASES[0]].size
    if any(D.size != Lc for D in D_obs.values()):
        if verbose:
            print(f"[audit_68] SKIP {pat}: ultrametric size mismatch")
        return None
    rho_obs = rho_split_from_phases(
        D_obs["rest_pre_A"], D_obs["rest_pre_B"],
        D_obs["task_test"], D_obs["rest_post"])

    n_swaps = swap_factor * (N_keep * (N_keep - 1)) // 2
    rho_surr = np.empty(n_surr, dtype=float)
    for r in range(n_surr):
        ok = True
        D_surr = {}
        for phase in PHASES:
            W_s = strength_preserving_shuffle(Ws[phase], n_swaps, rng)
            if not verify_strengths(Ws[phase], W_s, tol=1e-4):
                ok = False
                break
            D_surr[phase] = lrg_ultrametric_condensed(W_s)
            if D_surr[phase].size != Lc:
                ok = False
                break
        if not ok:
            rho_surr[r] = np.nan
            continue
        rho_surr[r] = rho_split_from_phases(
            D_surr["rest_pre_A"], D_surr["rest_pre_B"],
            D_surr["task_test"], D_surr["rest_post"])

    surr = rho_surr[np.isfinite(rho_surr)]
    if surr.size == 0:
        return None
    surr_mean = float(np.mean(surr))
    surr_std = float(np.std(surr, ddof=1))
    z = (rho_obs - surr_mean) / surr_std if surr_std > 0 else float("nan")
    p_one = float(np.mean(surr >= rho_obs))

    return {
        "patient": pat,
        "band": BAND,
        "N_full": int(N_full),
        "N_keep": N_keep,
        "n_epi_excluded": n_epi,
        "n_pairs": int(Lc),
        "n_swaps_per_surrogate": int(n_swaps),
        "n_surrogates": int(surr.size),
        "obs_rho": rho_obs,
        "surr_mean_rho": surr_mean,
        "surr_std_rho": surr_std,
        "surr_p5": float(np.quantile(surr, 0.05)),
        "surr_p50": float(np.quantile(surr, 0.50)),
        "surr_p95": float(np.quantile(surr, 0.95)),
        "obs_z": z,
        "obs_p_one_sided": p_one,
        "_surr_array": surr,
    }


def cohort_summary(rows: list[dict]) -> pd.DataFrame:
    rho_obs = np.array([r["obs_rho"] for r in rows])
    surr_med_per_pat = np.array([r["surr_p50"] for r in rows])
    n_above = int(sum(1 for r in rows if r["obs_p_one_sided"] < 0.05))
    try:
        wz, wp = wilcoxon(rho_obs - surr_med_per_pat, alternative="greater")
        cohort_z = float(wz)
        cohort_p = float(wp)
    except Exception:
        cohort_z = float("nan")
        cohort_p = float("nan")
    med_obs = float(np.median(rho_obs))
    med_surr_med = float(np.median(surr_med_per_pat))
    ratio = (abs(med_obs) / abs(med_surr_med)
             if abs(med_surr_med) > 1e-9 else float("inf"))
    if (cohort_z > 2 and abs(med_surr_med) < 0.05 and n_above >= 8):
        verdict = "separated"
    elif med_surr_med > 0.5 * med_obs and med_obs > 0:
        verdict = "also_positive"
    else:
        verdict = "intermediate"
    return pd.DataFrame([{
        "band": BAND,
        "n_patients": len(rows),
        "obs_median_rho": med_obs,
        "surr_median_rho_median": med_surr_med,
        "ratio": ratio,
        "paired_wilcoxon_z": cohort_z,
        "paired_wilcoxon_p": cohort_p,
        "n_above_surrogate": f"{n_above}/{len(rows)}",
        "verdict": verdict,
    }])


def make_cohort_figure(rows: list[dict], cohort_row: dict,
                        out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    sub_sorted = sorted(rows, key=lambda r: r["obs_rho"])
    ys = np.arange(len(sub_sorted))
    for y, r in zip(ys, sub_sorted):
        ax.hlines(y, r["surr_p5"], r["surr_p95"],
                  color="#aaaaaa", lw=4, alpha=0.85)
        ax.plot(r["surr_p50"], y, marker="|", color="#444444",
                markersize=10, mew=1.2)
        ax.plot(r["obs_rho"], y, marker="o", color="#1f3d6e",
                markersize=7, mec="white", mew=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels([r["patient"].replace("Pat_", "P")
                        + f" (-{r['n_epi_excluded']})"
                        for r in sub_sorted], fontsize=8)
    ax.axvline(0, color="0.6", lw=0.7, ls="--", zorder=0)
    ax.set_xlabel(r"$\rho_{\mathrm{split}}$")
    ax.set_title(
        f"{BRAIN_BAND_TEX_DICT[BAND]} (epi-excluded) — "
        f"verdict: {cohort_row['verdict']}\n"
        rf"Wilcoxon $p={cohort_row['paired_wilcoxon_p']:.4f}$,  "
        rf"$n_{{>}} = {cohort_row['n_above_surrogate']}$",
        fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_comparison_figure(rows: list[dict], full_csv: Path,
                            out_path: Path) -> None:
    if not full_csv.exists():
        return
    full = pd.read_csv(full_csv)
    full = full[full.band == BAND].set_index("patient")
    xs = []
    ys = []
    labels = []
    for r in rows:
        if r["patient"] not in full.index:
            continue
        xs.append(float(full.loc[r["patient"], "obs_rho"]))
        ys.append(r["obs_rho"])
        labels.append(r["patient"].replace("Pat_", "P"))
    if not xs:
        return
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    ax.scatter(xs, ys, s=40, color="#1f3d6e", edgecolor="white", linewidth=0.6)
    for x, y, lbl in zip(xs, ys, labels):
        ax.annotate(lbl, xy=(x, y), xytext=(4, 4),
                    textcoords="offset points", fontsize=7)
    lo = min(min(xs), min(ys)) - 0.05
    hi = max(max(xs), max(ys)) + 0.05
    ax.plot([lo, hi], [lo, hi], color="0.7", lw=0.7, ls="--")
    ax.axhline(0, color="0.7", lw=0.5, ls=":")
    ax.axvline(0, color="0.7", lw=0.5, ls=":")
    ax.set_xlabel(r"$\rho_{\mathrm{split}}$ (full FC, audit_63 α)")
    ax.set_ylabel(r"$\rho_{\mathrm{split}}$ (epi-excluded)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def comparison_table(rows: list[dict], full_csv: Path) -> pd.DataFrame:
    if not full_csv.exists():
        return pd.DataFrame()
    full = pd.read_csv(full_csv)
    full = (full[full.band == BAND]
            [["patient", "obs_rho", "surr_mean_rho", "obs_z",
              "obs_p_one_sided"]]
            .rename(columns=lambda c: c if c == "patient"
                    else f"{c}_full"))
    excl = pd.DataFrame([{
        "patient": r["patient"],
        "obs_rho_excluded": r["obs_rho"],
        "surr_mean_rho_excluded": r["surr_mean_rho"],
        "obs_z_excluded": r["obs_z"],
        "obs_p_one_sided_excluded": r["obs_p_one_sided"],
        "n_epi_excluded": r["n_epi_excluded"],
        "N_keep": r["N_keep"],
    } for r in rows])
    return full.merge(excl, on="patient")


def write_readme(rows: list[dict], cohort: pd.DataFrame,
                 comp: pd.DataFrame, runtime_s: float) -> None:
    row = cohort.iloc[0]
    full_med = comp.obs_rho_full.median() if "obs_rho_full" in comp else float("nan")
    excl_med = row["obs_median_rho"]
    direction = ("strengthened" if excl_med >= full_med
                 else "weakened")
    delta = excl_med - full_med
    epi_counts = [r["n_epi_excluded"] for r in rows]

    lines = [
        "---",
        "name: alpha_epi_exclusion",
        "scope: section_5_3_alpha_trace_anatomy_sensitivity",
        f"date: {time.strftime('%Y-%m-%d')}",
        "status: complete",
        "---",
        "",
        "# α epi-exclusion sensitivity (matched-strength split-baseline)",
        "",
        "**Head.** α LRG ρ_split is re-run after removing epileptic-zone "
        "contacts from each patient's FC (per-patient epi-mask via "
        "`load_epileptic_nodes` × `channel_labels.csv`). Cohort median "
        f"obs_rho {direction} from "
        f"{full_med:+.3f} (full FC, audit_63) to "
        f"{excl_med:+.3f} (epi-excluded) — Δ = {delta:+.3f}. "
        f"Wilcoxon p = {row['paired_wilcoxon_p']:.4f}, "
        f"n above own surrogate {row['n_above_surrogate']}, "
        f"verdict **{row['verdict']}**.",
        "",
        "## Cohort summary (epi-excluded α only)",
        "",
        "| n | obs median ρ | surr median ρ | ratio | Wilcoxon z | p | n>surr | verdict |",
        "|---|---|---|---|---|---|---|---|",
        f"| {row['n_patients']} "
        f"| {row['obs_median_rho']:+.3f} "
        f"| {row['surr_median_rho_median']:+.3f} "
        f"| {row['ratio']:.2f}× "
        f"| {row['paired_wilcoxon_z']:+.2f} "
        f"| {row['paired_wilcoxon_p']:.4f} "
        f"| {row['n_above_surrogate']} "
        f"| **{row['verdict']}** |",
        "",
        "## Per-patient comparison vs full FC",
        "",
        "| patient | N_full | n_epi excluded | obs_ρ full | obs_ρ excluded | Δρ |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in comp.iterrows():
        lines.append(
            f"| {r['patient']} | {r.get('N_keep', '?') + r.get('n_epi_excluded', 0)} "
            f"| {int(r.get('n_epi_excluded', 0))} "
            f"| {r['obs_rho_full']:+.3f} "
            f"| {r['obs_rho_excluded']:+.3f} "
            f"| {r['obs_rho_excluded'] - r['obs_rho_full']:+.3f} |"
        )

    lines.extend([
        "",
        "## Interpretation",
        "",
        "The §4 anatomy correlate noted Spearman ρ(frac_epi, α obs_ρ) = "
        "-0.41 (uncorrected p≈0.24) — patients with fewer epi contacts "
        "had stronger α trace. The predicate is: if α reorganization "
        "lives in clinically healthy circuits, removing epi nodes from "
        "the FC should leave the cohort signal intact or strengthen it; "
        "if α reorganization is partially carried by epi-zone activity, "
        "removing epi nodes should weaken the cohort signal. The result "
        f"is that the cohort median {direction} ({delta:+.3f}), with "
        f"epi-node counts per patient ranging "
        f"{min(epi_counts)}–{max(epi_counts)}.",
        "",
        "## Provenance",
        "",
        f"- N_surrogates = {N_SURROGATES} per cell",
        f"- n_swaps_per_surrogate = SWAP_FACTOR ({SWAP_FACTOR}) × N(N-1)/2 "
        "(N = epi-excluded count)",
        "- Cohort: " + ", ".join(COHORT),
        "- Band: " + BAND,
        "- Phases: " + ", ".join(PHASES),
        "- FC method: imcoh_abs",
        "- Epi mask: load_epileptic_nodes(pat) ∩ channel_labels.csv "
        "(FC-index-aligned)",
        "- LRG: τ = 1/λ_max, ultrametric via average linkage on Trho.",
        f"- Wall-clock runtime: {runtime_s:.1f} s",
        "- Build script: "
        "`scripts/01_compute/audit/audit_68_alpha_epi_exclusion.py`",
        "",
        "## Files",
        "",
        "- `per_patient.csv` — per-patient observed + surrogate stats",
        "- `cohort_summary.csv` — cohort verdict",
        "- `comparison.csv` — per-patient full vs excluded ρ_split",
        "- `figures/cohort.pdf` — per-patient observed vs surrogate "
        "P5-P95 spans",
        "- `figures/comparison.pdf` — full vs excluded scatter",
        "",
    ])
    (OUT / "README.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--n-surrogates", type=int, default=N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=SWAP_FACTOR)
    ap.add_argument("--seed", type=int, default=20260514)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    patients = args.patients
    n_surr = args.n_surrogates
    swap_factor = args.swap_factor

    print("[audit_68] pre-flight: ensure half FCs cached for α")
    for pat in patients:
        ensure_half_fcs(pat, [BAND])

    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    rows: list[dict] = []
    for pat in patients:
        t_cell = time.time()
        r = per_cell(pat, n_surr, swap_factor, rng, verbose=args.verbose)
        if r is None:
            continue
        rows.append(r)
        dt = time.time() - t_cell
        print(f"[audit_68] {pat}: N_full={r['N_full']} keep={r['N_keep']} "
              f"epi_excl={r['n_epi_excluded']} obs={r['obs_rho']:+.3f} "
              f"surr_mean={r['surr_mean_rho']:+.3f} z={r['obs_z']:+.2f} "
              f"p={r['obs_p_one_sided']:.3f} ({dt:.1f}s)")

    runtime = time.time() - t0
    print(f"[audit_68] total: {runtime:.1f}s for {len(rows)} cells")

    per_pat = pd.DataFrame([{k: v for k, v in r.items()
                              if k != "_surr_array"} for r in rows])
    per_pat.to_csv(OUT / "per_patient.csv", index=False)

    cohort = cohort_summary(rows)
    cohort.to_csv(OUT / "cohort_summary.csv", index=False)

    full_csv = ROOT / "data" / "audit" / \
        "matched_strength_surrogate_split_baseline" \
        / "per_patient_per_band.csv"
    comp = comparison_table(rows, full_csv)
    if not comp.empty:
        comp.to_csv(OUT / "comparison.csv", index=False)

    make_cohort_figure(rows, cohort.iloc[0].to_dict(), FIG / "cohort.pdf")
    make_comparison_figure(rows, full_csv, FIG / "comparison.pdf")

    write_readme(rows, cohort, comp, runtime)
    print(f"[audit_68] outputs at {OUT}")


if __name__ == "__main__":
    main()
