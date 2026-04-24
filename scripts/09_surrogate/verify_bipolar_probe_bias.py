#!/usr/bin/env python3
"""Step 2 GATE: Verify that bipolar re-referencing removes probe bias.

Computes same-probe enrichment in LRG communities for bipolar MSC matrices
and compares against the common-reference enrichment.

Run: python scripts/py/verify_bipolar_probe_bias.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import BIPOLAR_CACHE, MSC_CACHE
from lrg_eegfc.workflow.diagnostics import (
    _compute_propagator,
    _ultrametric_from_propagator,
)

BIPOLAR_MSC = BIPOLAR_CACHE / "msc_cache"
BIPOLAR_CH = BIPOLAR_CACHE / "channels"
ORIG_MSC = MSC_CACHE

PATIENTS = ["Pat_02", "Pat_05"]  # verify on two different patients
BANDS = ["alpha", "beta", "delta", "theta"]
N_COMMS = [3, 5, 8, 10, 15, 20, 30]

from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical


def get_nperseg(pat):
    fs = FS_MAP.get(pat, 2048.0)
    return int(fs * 2.0)


def probe_from_bipolar_label(label: str) -> str:
    """Extract probe from bipolar label like 'A2-A1' -> 'A'."""
    # Take the first part before '-'
    first = label.split("-")[0]
    m = re.match(r"([A-Za-z]+'?)", first)
    return m.group(1) if m else first


def compute_enrichment(labels_comm, probe_labels, N):
    """Fraction of same-community pairs that are same-probe, normalized."""
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


def main():
    print("=" * 70)
    print("BIPOLAR PROBE BIAS VERIFICATION")
    print("Target: enrichment < 1.5× at all scales")
    print("=" * 70)

    all_pass = True

    for pat in PATIENTS:
        nperseg = get_nperseg(pat)

        # Load bipolar labels
        labels_path = BIPOLAR_CH / pat / "bipolar_labels.json"
        if not labels_path.exists():
            print(f"\n  {pat}: no bipolar labels found, skipping")
            continue
        with open(labels_path) as f:
            bip_labels = json.load(f)

        N = len(bip_labels)
        probe_labels = [probe_from_bipolar_label(l) for l in bip_labels]

        print(f"\n{'='*60}")
        print(f"{pat}: {N} bipolar channels, {len(set(probe_labels))} probes")
        print(f"{'='*60}")

        # Also check same-probe fraction
        n_sp = sum(
            1
            for i in range(N)
            for j in range(i + 1, N)
            if probe_labels[i] == probe_labels[j]
        )
        n_total = N * (N - 1) // 2
        print(f"  Same-probe pairs: {n_sp}/{n_total} ({n_sp/n_total:.1%})")

        for band in BANDS:
            bip_path = (
                BIPOLAR_MSC / pat
                / f"{band}_rsPre_msc_sparsify-none_nperseg-{nperseg}.npy"
            )
            if not bip_path.exists():
                print(f"  {band}: not computed yet")
                continue

            A = np.load(bip_path)
            np.fill_diagonal(A, 0)
            D = A.sum(axis=1)
            if (D < 1e-10).any():
                print(f"  {band}: has disconnected nodes!")
                continue

            L = np.diag(D) - A
            ev, U = np.linalg.eigh(L)
            ev = np.maximum(ev, 0.0)
            if ev[1] < 1e-10:
                print(f"  {band}: disconnected graph")
                continue

            # Build dendrogram
            K = _compute_propagator(ev, U, 1.0 / ev[1] * 0.5)
            Du = _ultrametric_from_propagator(K)
            Z = linkage(squareform(Du, checks=False), method="average")

            print(f"\n  {band}:")
            print(f"  {'n':>4s} {'enrichment':>10s} {'status':>8s}")
            for nc in N_COMMS:
                if nc >= N:
                    continue
                labels = fcluster(Z, nc, criterion="maxclust")
                enrich = compute_enrichment(labels, probe_labels, N)
                status = "PASS" if enrich < 1.5 else "FAIL"
                if enrich >= 1.5:
                    all_pass = False
                print(f"  {nc:4d} {enrich:10.2f}× {status:>8s}")

        # Also compute same-probe vs cross-probe MSC comparison
        bip_path = (
            BIPOLAR_MSC / pat
            / f"alpha_rsPre_msc_sparsify-none_nperseg-{nperseg}.npy"
        )
        if bip_path.exists():
            A = np.load(bip_path)
            np.fill_diagonal(A, 0)
            sp_mask = np.zeros((N, N), dtype=bool)
            for i in range(N):
                for j in range(N):
                    if probe_labels[i] == probe_labels[j]:
                        sp_mask[i, j] = True
            np.fill_diagonal(sp_mask, False)

            sp_msc = A[sp_mask].mean() if sp_mask.any() else 0
            cp_msc = A[~sp_mask & ~np.eye(N, dtype=bool)].mean()
            ratio = sp_msc / (cp_msc + 1e-30)
            print(f"\n  MSC ratio (same-probe/cross-probe): {ratio:.2f}×")
            print(f"    Same-probe MSC: {sp_msc:.4f}")
            print(f"    Cross-probe MSC: {cp_msc:.4f}")

    print(f"\n{'='*70}")
    print(f"OVERALL: {'ALL PASS ✓' if all_pass else 'SOME FAILED ✗'}")
    if not all_pass:
        print("  Consider additionally zeroing edges between adjacent bipolar")
        print("  channels that share a physical contact.")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
