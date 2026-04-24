#!/usr/bin/env python3
"""Verify that ImCoh FC removes (or reduces) same-probe enrichment bias.

ImCoh (imaginary coherence) discards zero-lag interactions, which should
attenuate volume conduction between physically close contacts on the same
probe. This script computes same-probe enrichment in LRG communities at
multiple granularities and compares against the standard MSC enrichment.

Target: enrichment < 1.5x at all scales.

Run: python scripts/py/verify_imcoh_probe_bias.py [--patients Pat_02 Pat_05] [-v]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, nperseg_for_fs
from lrg_eegfc.config.paths import (
    IMCOH_CACHE, IMCOH_LRG_CACHE, LRG_CACHE as MSC_LRG_CACHE,
    MSC_CACHE, SEEG_DATAPATH,
)
from lrg_eegfc.workflow.lrg import LRGResult

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
DATASET_ROOT = SEEG_DATAPATH

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
N_COMMS = [3, 5, 10, 15, 20, 30]
from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical


def get_nperseg(pat: str) -> int:
    fs = FS_MAP.get(pat, 2048.0)
    return int(fs * 2.0)


def load_channel_labels(patient: str) -> list[str]:
    """Load monopolar channel labels (same channels used by ImCoh)."""
    # Try .txt first
    txt_file = DATASET_ROOT / patient / "channel_labels.txt"
    if txt_file.exists():
        labels = []
        with open(txt_file, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    if "," in line:
                        label = line.split(",")[0]
                    else:
                        label = line
                    label = label.replace(" ", "")
                    labels.append(label)
        if labels:
            return labels

    # Try .csv
    csv_file = DATASET_ROOT / patient / "channel_labels.csv"
    if csv_file.exists():
        labels = []
        with open(csv_file, "r") as f:
            for i, line in enumerate(f):
                if i == 0 and "label" in line.lower():
                    continue
                line = line.strip().strip('"')
                if line:
                    if "," in line:
                        label = line.split(",")[0]
                    else:
                        label = line
                    label = label.replace(" ", "")
                    labels.append(label)
        if labels:
            return labels

    raise FileNotFoundError(
        f"No channel_labels.txt or .csv found for {patient} in {DATASET_ROOT}"
    )


def probe_from_monopolar_label(label: str) -> str:
    """Extract probe name from monopolar sEEG label like 'A1' -> 'A'."""
    # Normalise prime characters
    label = label.replace("\u00ec", "'")
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


def compute_enrichment(
    labels_comm: np.ndarray,
    probe_labels: list[str],
    N: int,
) -> float:
    """Fraction of same-community pairs that are same-probe, normalized.

    Returns enrichment ratio: 1.0 = no bias, >1.0 = same-probe overrepresented.
    """
    same_comm = 0
    same_comm_and_probe = 0
    same_probe_total = 0
    for i in range(N):
        for j in range(i + 1, N):
            if labels_comm[i] == labels_comm[j]:
                same_comm += 1
                if probe_labels[i] == probe_labels[j]:
                    same_comm_and_probe += 1
            if probe_labels[i] == probe_labels[j]:
                same_probe_total += 1

    if same_comm == 0:
        return 0.0
    frac_sp_in_sc = same_comm_and_probe / same_comm
    expected = same_probe_total / (N * (N - 1) / 2)
    return frac_sp_in_sc / expected if expected > 0 else 0.0


def load_lrg_result(cache_dir: Path, pat: str, phase: str, band: str,
                    fc_tag: str) -> LRGResult | None:
    """Load an LRG result from a custom cache directory."""
    path = cache_dir / pat / f"{band}_{phase}_lrg_{fc_tag}.npz"
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


def main():
    parser = argparse.ArgumentParser(
        description="Verify ImCoh probe bias removal."
    )
    parser.add_argument("--patients", nargs="+", default=PATIENTS,
                        help="Patient IDs to verify")
    parser.add_argument("--bands", nargs="+", default=BANDS,
                        help="Frequency bands to check")
    parser.add_argument("--phase", default="rest_pre",
                        help="Phase to use for verification (default: rest_pre)")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    print("=" * 70)
    print("IMCOH PROBE BIAS VERIFICATION")
    print("Target: enrichment < 1.5x at all scales")
    print("=" * 70)

    all_pass = True
    results_table = []  # (pat, band, n_comm, imcoh_enrich, msc_enrich)

    for pat in args.patients:
        nperseg = get_nperseg(pat)

        # Load channel labels (monopolar -- same as ImCoh)
        try:
            ch_labels = load_channel_labels(pat)
        except FileNotFoundError as e:
            print(f"\n  {pat}: {e}")
            continue

        N_ch = len(ch_labels)
        probe_labels = [probe_from_monopolar_label(l) for l in ch_labels]

        print(f"\n{'='*60}")
        print(f"{pat}: {N_ch} channels, {len(set(probe_labels))} probes")
        print(f"{'='*60}")

        # Same-probe pair fraction
        n_sp = sum(
            1 for i in range(N_ch) for j in range(i + 1, N_ch)
            if probe_labels[i] == probe_labels[j]
        )
        n_total = N_ch * (N_ch - 1) // 2
        print(f"  Same-probe pairs: {n_sp}/{n_total} ({n_sp/n_total:.1%})")

        for band in args.bands:
            # Load ImCoh LRG result
            imcoh_res = load_lrg_result(
                IMCOH_LRG_CACHE, pat, args.phase, band, "imcoh-abs"
            )
            if imcoh_res is None:
                if args.verbose:
                    print(f"  {band}: no ImCoh LRG result, skipping")
                continue

            # Load MSC LRG result for comparison
            msc_res = load_lrg_result(
                MSC_LRG_CACHE, pat, args.phase, band, "msc"
            )

            # ImCoh preserves all channels; the LRG giant component may be
            # smaller than N_ch. Use probe labels for the giant component
            # nodes only.
            n_gc = imcoh_res.n_nodes
            gc_probes = probe_labels[:n_gc]  # giant component is first n nodes

            print(f"\n  {band} (n_gc={n_gc}):")
            print(f"  {'n':>4s} {'ImCoh':>10s} {'MSC':>10s} {'status':>8s}")

            Z = imcoh_res.linkage_matrix
            Z_msc = msc_res.linkage_matrix if msc_res is not None else None

            for nc in N_COMMS:
                if nc >= n_gc:
                    continue
                labels_imcoh = fcluster(Z, nc, criterion="maxclust")
                enrich_imcoh = compute_enrichment(labels_imcoh, gc_probes, n_gc)

                enrich_msc = np.nan
                if Z_msc is not None and msc_res.n_nodes == n_gc:
                    labels_msc = fcluster(Z_msc, nc, criterion="maxclust")
                    enrich_msc = compute_enrichment(labels_msc, gc_probes, n_gc)

                status = "PASS" if enrich_imcoh < 1.5 else "FAIL"
                if enrich_imcoh >= 1.5:
                    all_pass = False

                msc_str = f"{enrich_msc:.2f}x" if not np.isnan(enrich_msc) else "n/a"
                print(f"  {nc:4d} {enrich_imcoh:9.2f}x {msc_str:>10s} {status:>8s}")

                results_table.append({
                    "patient": pat, "band": band, "n_comm": nc,
                    "imcoh_enrichment": enrich_imcoh,
                    "msc_enrichment": enrich_msc,
                    "pass": enrich_imcoh < 1.5,
                })

        # Same-probe vs cross-probe ImCoh comparison (alpha band)
        imcoh_path = (
            IMCOH_CACHE / pat
            / f"alpha_{args.phase}_imcoh_sparsify-none_nperseg-{nperseg}.npy"
        )
        if imcoh_path.exists():
            A = np.load(imcoh_path)
            np.fill_diagonal(A, 0)
            N = A.shape[0]
            sp_mask = np.zeros((N, N), dtype=bool)
            for i in range(N):
                for j in range(N):
                    if probe_labels[i] == probe_labels[j]:
                        sp_mask[i, j] = True
            np.fill_diagonal(sp_mask, False)

            sp_fc = np.abs(A[sp_mask]).mean() if sp_mask.any() else 0
            cp_fc = np.abs(A[~sp_mask & ~np.eye(N, dtype=bool)]).mean()
            ratio = sp_fc / (cp_fc + 1e-30)
            print(f"\n  ImCoh ratio (same-probe/cross-probe): {ratio:.2f}x")
            print(f"    Same-probe |ImCoh|: {sp_fc:.4f}")
            print(f"    Cross-probe |ImCoh|: {cp_fc:.4f}")

        # Also load MSC for same comparison
        msc_path = (
            MSC_CACHE / pat
            / f"alpha_{args.phase}_msc_sparsify-none_nperseg-{nperseg}.npy"
        )
        if msc_path.exists():
            A_msc = np.load(msc_path)
            np.fill_diagonal(A_msc, 0)
            N_msc = A_msc.shape[0]
            sp_mask_msc = np.zeros((N_msc, N_msc), dtype=bool)
            for i in range(N_msc):
                for j in range(N_msc):
                    if probe_labels[i] == probe_labels[j]:
                        sp_mask_msc[i, j] = True
            np.fill_diagonal(sp_mask_msc, False)

            sp_msc = A_msc[sp_mask_msc].mean() if sp_mask_msc.any() else 0
            cp_msc = A_msc[~sp_mask_msc & ~np.eye(N_msc, dtype=bool)].mean()
            ratio_msc = sp_msc / (cp_msc + 1e-30)
            print(f"  MSC ratio (same-probe/cross-probe):   {ratio_msc:.2f}x")
            print(f"    Same-probe MSC: {sp_msc:.4f}")
            print(f"    Cross-probe MSC: {cp_msc:.4f}")

    # Save results table
    if results_table:
        import pandas as pd
        from lrg_eegfc.config.paths import REPORTS_ROOT
        out_dir = REPORTS_ROOT / "imcoh"
        out_dir.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(results_table)
        csv_path = out_dir / "probe_bias_verification.csv"
        df.to_csv(csv_path, index=False)
        print(f"\nResults saved: {csv_path}")

    print(f"\n{'='*70}")
    print(f"OVERALL: {'ALL PASS' if all_pass else 'SOME FAILED'}")
    if not all_pass:
        print("  Some enrichment values exceed 1.5x. ImCoh may not fully")
        print("  remove volume-conduction bias at all scales.")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
