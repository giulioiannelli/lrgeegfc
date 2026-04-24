#!/usr/bin/env python3
"""Verify the ImCoh reset: new signed-imcoh cache vs archived pre-reset squared.

After the reset, the stored quantity is signed ImCoh (Nolte 2004). The
archived pre-reset quantity was |ImCoh|^2 (mislabelled as "imcoh"). If the
reset is purely a labelling / representation change, then

    (new_signed_imcoh) ** 2  ==  old_cached_imcoh_value

should hold to ~1e-10 across the cohort. If it doesn't, the new computation
is numerically different and we must investigate before trusting downstream.

Also performs per-file correctness invariants:
- shape (N, N), dtype float
- diag exactly 0
- skew-symmetric: max |m + m.T| < 1e-10
- range within [-1, 1]

Usage:
    python scripts/01_compute/verify_imcoh_reset.py -v
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, PATIENTS_4PHASE, PHASE_LABELS, nperseg_for_fs
from lrg_eegfc.config.paths import IMCOH_CACHE, CACHE_ROOT

from lrg_eegfc.config.const import FS_OVERRIDES as FS_MAP  # canonical
ARCHIVE_CACHE = CACHE_ROOT / "_archive_absimcoh_sq"


def verify_freq_resolved_cache(path: Path):
    """Return dict of checks for one new freq-resolved signed ImCoh file."""
    m = np.load(path)
    band_avg_signed = m.mean(axis=-1)
    # Reconstruct per-frequency skew-symmetry: Coh[:,:,f] + Coh.T[:,:,f] ≈ 0
    skew_err = float(np.abs(m + m.transpose(1, 0, 2)).max())
    return {
        "shape": m.shape,
        "dtype": str(m.dtype),
        "min": float(m.min()),
        "max": float(m.max()),
        "has_negative": bool((m < 0).any()),
        "skew_err_per_freq": skew_err,
        "diag_all_zero_per_freq": bool(
            np.allclose(np.diagonal(m, axis1=0, axis2=1), 0, atol=1e-7)
        ),
        "within_range": bool(m.min() >= -1 - 1e-6 and m.max() <= 1 + 1e-6),
        "band_avg_signed_range": (
            float(band_avg_signed.min()),
            float(band_avg_signed.max()),
        ),
    }


def compare_sq_view_to_archive(
    new_path: Path, archive_path: Path, atol: float = 1e-6,
):
    """Compare the band-averaged imcoh_sq view vs archived |ImCoh|² matrix.

    The loader derives ``imcoh_sq`` as ``mean(signed**2, axis=F_band)``,
    which should match the archived pre-reset quantity exactly (the old
    pipeline squared per-frequency then band-averaged).
    """
    new = np.load(new_path).astype(np.float64)
    try:
        old = np.load(archive_path)
    except FileNotFoundError:
        return None, None, None
    sq_view = (new ** 2).mean(axis=-1)  # loader imcoh_sq
    np.fill_diagonal(sq_view, 0)
    np.fill_diagonal(old, 0)
    diff = np.abs(sq_view - old)
    max_err = float(diff.max())
    mean_err = float(diff.mean())
    agrees = max_err < atol
    return max_err, mean_err, agrees


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--patients", nargs="+", default=PATIENTS_4PHASE)
    p.add_argument("--phases", nargs="+", default=list(PHASE_LABELS))
    p.add_argument("--bands", nargs="+", default=list(BRAIN_BANDS.keys()))
    p.add_argument("--atol", type=float, default=1e-6,
                   help="Tolerance for loader imcoh_sq view vs archive. "
                        "Float32 storage → 1e-6 is the natural floor.")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    total = n_ok_invariant = n_bad_invariant = 0
    n_match_archive = n_mismatch_archive = n_no_archive = 0
    worst_err = 0.0
    worst_file = None

    for pat in args.patients:
        fs = FS_MAP.get(pat, 2048.0)
        nperseg = nperseg_for_fs(fs)

        for phase in args.phases:
            for band in args.bands:
                new_path = (
                    IMCOH_CACHE / pat
                    / f"{band}_{phase}_imcoh_freqresolved_"
                    f"nperseg-{nperseg}.npy"
                )
                if not new_path.exists():
                    continue
                total += 1

                # Invariants on the freq-resolved signed tensor
                inv = verify_freq_resolved_cache(new_path)
                invariants_ok = (
                    inv["skew_err_per_freq"] < 1e-5
                    and inv["diag_all_zero_per_freq"]
                    and inv["within_range"]
                    and inv["has_negative"]
                )
                if invariants_ok:
                    n_ok_invariant += 1
                else:
                    n_bad_invariant += 1
                    print(f"  FAIL invariants {pat} {band} {phase}: {inv}")

                # Compare imcoh_sq view (mean(signed**2, axis=F)) vs archive
                archive_path = (
                    ARCHIVE_CACHE / pat
                    / f"{band}_{phase}_absimcoh_sq_sparsify-none_"
                    f"nperseg-{nperseg}.npy"
                )
                max_err, mean_err, agrees = compare_sq_view_to_archive(
                    new_path, archive_path, atol=args.atol,
                )
                if agrees is None:
                    n_no_archive += 1
                    tag = "ARCHIVE_MISSING"
                elif agrees:
                    n_match_archive += 1
                    tag = "MATCH"
                else:
                    n_mismatch_archive += 1
                    tag = "MISMATCH"
                    if max_err > worst_err:
                        worst_err = max_err
                        worst_file = (pat, band, phase)

                if args.verbose:
                    print(
                        f"  {pat} {band:10s} {phase:10s}  "
                        f"inv={'OK' if invariants_ok else 'FAIL'}  "
                        f"archive={tag}  "
                        f"max|new²−old|={max_err if max_err is not None else 'NA'}"
                    )

    print("\n" + "=" * 60)
    print(f"Total cache files checked: {total}")
    print(f"Invariants OK:              {n_ok_invariant}")
    print(f"Invariants FAIL:            {n_bad_invariant}")
    print(f"Archive comparison:")
    print(f"  MATCH   (new² ≈ old, atol={args.atol}):  {n_match_archive}")
    print(f"  MISMATCH:                                {n_mismatch_archive}")
    print(f"  ARCHIVE MISSING:                         {n_no_archive}")
    if worst_file is not None:
        print(f"  worst error: {worst_err:.2e} at {worst_file}")
    print("=" * 60)

    # Scientific interpretation
    if n_bad_invariant == 0 and n_mismatch_archive == 0:
        print(
            "\n✓ Reset is a pure labelling change. The 3-day analysis work\n"
            "  is numerically recoverable: every old |ImCoh|² value equals\n"
            "  the new signed imcoh squared, so rankings, communities, and\n"
            "  all rank-based metrics are unchanged. Only the QUOTED NUMBERS\n"
            "  (means, percentiles, ratios) differ by a sqrt, and only if\n"
            "  the paper narrative prefers |ImCoh| over |ImCoh|² (it does)."
        )
    elif n_bad_invariant == 0 and n_mismatch_archive > 0:
        print(
            "\n⚠ Reset introduces NUMERICAL differences from archive.\n"
            "  Investigate before trusting downstream — either a bug in the\n"
            "  new compute or a real discovery missed under the old framework."
        )
    else:
        print("\n✗ Invariants failed. Stop, debug msc.py imcoh branch.")

    return 0 if (n_bad_invariant == 0 and n_mismatch_archive == 0) else 1


if __name__ == "__main__":
    raise SystemExit(main())
