#!/usr/bin/env python3
"""Batch visualization of dense vs validated MSC matrices."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, MSC_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.io import list_patients
from lrg_eegfc.visuals.msc import plot_msc_comparison_dense_vs_validated


def _resolve_list(value: List[str] | None, default: List[str]) -> List[str]:
    return value if value else default


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create dense vs validated MSC comparison figures."
    )
    parser.add_argument("--patients", nargs="+", default=None)
    parser.add_argument("--bands", nargs="+", default=None)
    parser.add_argument("--phases", nargs="+", default=None)
    parser.add_argument("--n-surrogates", type=int, default=200)
    parser.add_argument("--nperseg", type=int, default=1024)
    parser.add_argument("--cache-root", type=Path, default=MSC_CACHE)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=FIGURES_ROOT / "msc_validation",
    )
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    patients = _resolve_list(args.patients, list_patients(SEEG_DATAPATH))
    bands = _resolve_list(args.bands, list(BRAIN_BANDS.keys()))
    phases = _resolve_list(args.phases, list(PHASE_LABELS))

    for patient in patients:
        for band in bands:
            for phase in phases:
                output_path = (
                    args.output_dir
                    / patient
                    / f"{band}_{phase}_msc_validation_nsurr{args.n_surrogates}.png"
                )
                if args.skip_existing and output_path.exists():
                    if args.verbose:
                        print(f"• Skipping (exists): {output_path}")
                    continue
                try:
                    plot_msc_comparison_dense_vs_validated(
                        patient=patient,
                        phase=phase,
                        band=band,
                        n_surrogates=args.n_surrogates,
                        cache_root=args.cache_root,
                        nperseg=args.nperseg,
                        output_path=output_path,
                    )
                    if args.verbose:
                        print(f"✓ Saved: {output_path}")
                except FileNotFoundError as exc:
                    if args.verbose:
                        print(f"[skip] {exc}")
                except Exception as exc:  # noqa: BLE001
                    print(f"[error] {patient} {band} {phase}: {exc}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
