#!/usr/bin/env python3
"""Compute multiscale VI profiles and hypothesis contrasts on ImCoh LRG results.

Post-reset canonical pipeline: FC = |ImCoh| (transform = "abs" applied per
frequency bin, then band-averaged).  Only `imcoh_abs` LRG results are read.

For each (patient, band, phase-pair), compute VI(k) across the full range of
dendrogram cut levels k = 2, 3, ..., N-1 (singleton leaves excluded). Then
evaluate the hypotheses at each scale.

Hypotheses (positive contrast = supported). H2 (H2a + H2b) is the central
hypothesis of the research:
  H1  (task stability):    VI(TL,TT) < mean(other pairs)
  H2a (task trace):        VI(rest_pre,rest_post) > VI(taskT,rest_post)
  H2b (task approach):     VI(rest_pre,taskT)  > VI(taskT,rest_post)
  H3  (within < cross):    mean(VI_within) < mean(VI_cross)

H4 (frequency gradient) is a cross-patient concordance analysis; it is
computed downstream in scripts/01_compute/report_h1h4_vi.py from the raw
VI profiles (Kendall's W + Spearman rho, following the method in
scripts/05_multiscale/analyze_h4_gradient.py).

Outputs (under data/reports/imcoh_vi/):
  vi_raw_profiles.csv          -- full VI(k) for every combination
  hypothesis_contrasts.csv     -- per-patient contrast at each (band, k, hyp)

Run: python scripts/01_compute/compute_imcoh_vi.py [--patients Pat_02] [-v]
     python scripts/01_compute/compute_imcoh_vi.py --dry-run   # coverage only
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    ALL_PHASE_PAIRS,
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PATIENTS_LIST,
    PHASE_LABELS,
    classify_pair,
)
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.lrg import LRGResult

PATIENTS = list(PATIENTS_LIST)  # all 6 patients, no exclusions
OUT = REPORTS_ROOT / "imcoh_vi"

BAND_ORDER = BRAIN_BANDS_NAMES
PHASES = list(PHASE_LABELS)
ALL_PAIRS = list(ALL_PHASE_PAIRS)

WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]

PAIR_SHORT = {
    ("rest_pre", "rest_post"): "Pre-Post",
    ("task_learn", "task_test"): "TL-TT",
    ("rest_pre", "task_learn"): "Pre-TL",
    ("rest_pre", "task_test"): "Pre-TT",
    ("task_learn", "rest_post"): "TL-Post",
    ("task_test", "rest_post"): "TT-Post",
}


# ── VI computation (canonical implementation) ──────────────────────────
from lrg_eegfc.utils.metrics import compute_vi


# ── Load LRG result from ImCoh cache ───────────────────────────────────
def load_imcoh_lrg(pat: str, phase: str, band: str) -> LRGResult | None:
    path = IMCOH_LRG_CACHE / pat / f"{band}_{phase}_lrg_imcoh-abs.npz"
    if not path.exists():
        return None
    data = np.load(path)
    return LRGResult(
        ultrametric_matrix=data["ultrametric_matrix"],
        linkage_matrix=data["linkage_matrix"],
        entropy_tau=data["entropy_tau"],
        entropy_1_minus_S=data["entropy_1_minus_S"],
        entropy_C=data["entropy_C"],
        optimal_threshold=float(data["optimal_threshold"]),
        patient=str(data["patient"]),
        phase=str(data["phase"]),
        band=str(data["band"]),
        fc_method=str(data["fc_method"]),
        n_nodes=int(data["n_nodes"]),
    )


# ── VI helper ──────────────────────────────────────────────────────────
def _get_vi(pdf: pd.DataFrame, pa: str, pb: str) -> float:
    """Get VI value for a specific pair from a patient-k slice."""
    row = pdf[
        ((pdf["phase_a"] == pa) & (pdf["phase_b"] == pb)) |
        ((pdf["phase_a"] == pb) & (pdf["phase_b"] == pa))
    ]
    if row.empty:
        return np.nan
    return row["vi"].iloc[0]


# ── Hypothesis contrasts ───────────────────────────────────────────────
def contrast_H1(pdf: pd.DataFrame) -> float:
    """H1: task stability. POSITIVE = supported.
    mean(other pairs VI) - VI(taskL, taskT)."""
    vi_tt = _get_vi(pdf, "task_learn", "task_test")
    other_vis = []
    for pa, pb in ALL_PAIRS:
        if (pa, pb) == ("task_learn", "task_test"):
            continue
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            other_vis.append(v)
    if np.isnan(vi_tt) or len(other_vis) == 0:
        return np.nan
    return np.mean(other_vis) - vi_tt


def contrast_H2a(pdf: pd.DataFrame) -> float:
    """H2a: task trace. POSITIVE = trace detected.
    VI(rest_pre,rest_post) - VI(taskT,rest_post)."""
    return _get_vi(pdf, "rest_pre", "rest_post") - _get_vi(pdf, "task_test", "rest_post")


def contrast_H2b(pdf: pd.DataFrame) -> float:
    """H2b: task approach. POSITIVE = approach detected.
    VI(rest_pre,taskT) - VI(taskT,rest_post)."""
    return _get_vi(pdf, "rest_pre", "task_test") - _get_vi(pdf, "task_test", "rest_post")


def contrast_H3(pdf: pd.DataFrame) -> float:
    """H3: within < cross. POSITIVE = supported.
    mean(cross VI) - mean(within VI)."""
    within_vis = []
    for pa, pb in WITHIN_PAIRS:
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            within_vis.append(v)
    cross_vis = []
    for pa, pb in CROSS_PAIRS:
        v = _get_vi(pdf, pa, pb)
        if not np.isnan(v):
            cross_vis.append(v)
    if not within_vis or not cross_vis:
        return np.nan
    return np.mean(cross_vis) - np.mean(within_vis)


HYPOTHESES = {
    "H1": contrast_H1,
    "H2a": contrast_H2a,
    "H2b": contrast_H2b,
    "H3": contrast_H3,
}


# ── Coverage guard ─────────────────────────────────────────────────────
def check_coverage(patients: list[str], bands: list[str]) -> tuple[list[str], list[str]]:
    """List expected vs missing LRG cache files for the post-reset grid.

    Returns (expected_paths, missing_paths).  Pat_06 is skipped because it
    lacks task phases (CLAUDE.md invariant 3) — excluded upstream via
    PATIENTS_4PHASE.
    """
    expected: list[str] = []
    missing: list[str] = []
    for pat in patients:
        for band in bands:
            for phase in PHASES:
                p = IMCOH_LRG_CACHE / pat / f"{band}_{phase}_lrg_imcoh-abs.npz"
                expected.append(str(p))
                if not p.exists():
                    missing.append(str(p))
    return expected, missing


# ── Step 1: Compute VI(k) profiles ────────────────────────────────────
def compute_all_vi_profiles(patients: list[str], bands: list[str],
                            verbose: bool) -> pd.DataFrame:
    """Compute VI at every k from 2 to N/2 for all combinations."""
    print("Computing VI profiles...")
    rows = []

    for pat in patients:
        for band in bands:
            # Load all 4 phase results
            results = {}
            for phase in PHASES:
                res = load_imcoh_lrg(pat, phase, band)
                if res is not None:
                    results[phase] = res

            if len(results) < 2:
                if verbose:
                    print(f"  SKIP {pat} {band}: only {len(results)} phases")
                continue

            # k range: 2 to N-1 (full dendrogram, excluding trivial k=1 and
            # the singleton-leaves cut k=N). Covers both coarse and fine scales.
            n_nodes = min(r.n_nodes for r in results.values())
            max_k = max(n_nodes - 1, 3)
            k_range = list(range(2, max_k + 1))

            for pa, pb in ALL_PAIRS:
                if pa not in results or pb not in results:
                    continue
                ra, rb = results[pa], results[pb]

                # Require same n_nodes for partition comparison
                if ra.n_nodes != rb.n_nodes:
                    if verbose:
                        print(f"  WARN {pat} {band} {pa}-{pb}: "
                              f"n_nodes mismatch ({ra.n_nodes} vs {rb.n_nodes})")
                    continue

                Za, Zb = ra.linkage_matrix, rb.linkage_matrix

                for k in k_range:
                    la = fcluster(Za, k, criterion="maxclust")
                    lb = fcluster(Zb, k, criterion="maxclust")
                    vi = compute_vi(la, lb)
                    nvi = vi / np.log(max(k, 2))

                    rows.append({
                        "patient": pat,
                        "band": band,
                        "phase_a": pa,
                        "phase_b": pb,
                        "pair": PAIR_SHORT.get((pa, pb), f"{pa}-{pb}"),
                        "pair_type": classify_pair(pa, pb),
                        "k": k,
                        "vi": vi,
                        "nvi": nvi,
                    })

        if verbose:
            print(f"  {pat} done")

    df = pd.DataFrame(rows)
    print(f"  Total: {len(df)} profile points")
    return df


# ── Step 2: Compute hypothesis contrasts ──────────────────────────────
def compute_hypothesis_contrasts(vi_df: pd.DataFrame,
                                 patients: list[str]) -> pd.DataFrame:
    """Compute per-patient contrast at each (band, k, hypothesis)."""
    print("Computing hypothesis contrasts...")
    rows = []

    for hyp_name, contrast_fn in HYPOTHESES.items():
        for band in BAND_ORDER:
            bdf = vi_df[vi_df["band"] == band]
            k_values = sorted(bdf["k"].unique())

            for k in k_values:
                kdf = bdf[bdf["k"] == k]

                for pat in patients:
                    pdf = kdf[kdf["patient"] == pat]
                    if pdf.empty:
                        continue
                    contrast = contrast_fn(pdf)
                    if contrast is not None and not np.isnan(contrast):
                        rows.append({
                            "hypothesis": hyp_name,
                            "band": band,
                            "k": k,
                            "patient": pat,
                            "contrast": contrast,
                            "sign": int(np.sign(contrast)),
                        })

    df = pd.DataFrame(rows)
    print(f"  H1-H3 total: {len(df)} contrast values")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Compute multiscale VI profiles on ImCoh LRG results."
    )
    parser.add_argument("--patients", nargs="+", default=PATIENTS)
    parser.add_argument("--bands", nargs="+", default=BAND_ORDER)
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--dry-run", action="store_true",
                        help="Only report cache coverage; do not compute.")
    args = parser.parse_args()

    # Coverage guard (fail-fast, but do not abort the run — report gaps)
    expected, missing = check_coverage(args.patients, args.bands)
    print(f"Coverage: {len(expected) - len(missing)}/{len(expected)} "
          f"imcoh-abs LRG files present")
    if missing:
        print(f"  MISSING ({len(missing)}):")
        for m in missing[:20]:
            print(f"    {m}")
        if len(missing) > 20:
            print(f"    ... and {len(missing) - 20} more")
    if args.dry_run:
        return

    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # Step 1: VI profiles
    vi_df = compute_all_vi_profiles(args.patients, args.bands, args.verbose)
    if vi_df.empty:
        print("ERROR: No VI profiles computed. Check that ImCoh LRG results exist.")
        sys.exit(1)

    csv_path = OUT / "vi_raw_profiles.csv"
    vi_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Step 2: Hypothesis contrasts (H1-H3 per-patient-per-(band,k))
    contrasts_df = compute_hypothesis_contrasts(vi_df, args.patients)
    contrasts_path = OUT / "hypothesis_contrasts.csv"
    contrasts_df.to_csv(contrasts_path, index=False)
    print(f"  Saved: {contrasts_path}")

    # Summary: count unanimous cells per hypothesis.
    # Unanimity = all patients PRESENT at that (band, k) cell agree on sign.
    # Pat_06 only has rest phases → absent from hypothesis cells requiring
    # task phases, which is fine — it simply doesn't contribute to those cells.
    print("\n--- Summary ---")
    for hyp in HYPOTHESES:
        hdf = contrasts_df[contrasts_df["hypothesis"] == hyp]
        if hdf.empty:
            print(f"  {hyp}: no data")
            continue

        n_unanimous_pos = 0
        n_unanimous_neg = 0
        n_cells = 0
        for (band, k), grp in hdf.groupby(["band", "k"]):
            n_cells += 1
            signs = grp["sign"].values
            if len(signs) >= 2:
                if (signs > 0).all():
                    n_unanimous_pos += 1
                elif (signs < 0).all():
                    n_unanimous_neg += 1

        print(f"  {hyp}: {n_cells} (band,k) cells, "
              f"{n_unanimous_pos} unanimous +, "
              f"{n_unanimous_neg} unanimous -")

    elapsed = time.time() - t0
    print(f"\nDone ({elapsed:.1f}s). Outputs in {OUT}")


if __name__ == "__main__":
    main()
