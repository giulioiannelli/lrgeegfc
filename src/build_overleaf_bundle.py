#!/usr/bin/env python3
"""Collect figures/tables into an Overleaf bundle with a manifest."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Dict, List


def _patterns(patient: str | None) -> Dict[str, List[str]]:
    pat = patient or "Pat_*"
    return {
        "A": [f"data/figures/timeseries/{pat}/*.png"],
        "B": [f"data/figures/msc/{pat}/*msc_*network.png"],
        "C": [f"data/figures/msc/{pat}/phase_band_grid.png"],
        "D": [f"data/figures/msc/all_patients_*.png"],
        "E": [f"data/figures/comparison/{pat}/*_fc_comparison.png"],
        "F": [f"data/figures/comparison/{pat}/*network_variants*.png"],
        "G": [f"data/figures/cleaning/{pat}/*.png"],
        "H": [f"data/figures/msc_validation/{pat}/*.png"],
        "I": [f"data/figures/lrg/{pat}/*.png"],
        "L": [
            f"data/figures/lrg_video/{pat}/*/*.png",
            f"data/figures/lrg_video/{pat}/*/*.gif",
        ],
        "M": [f"data/figures/lrg/{pat}/*entropy*.png"],
        "N": [f"data/figures/lrg/{pat}/*dendrogram*.png"],
        "O": [f"data/figures/reorganization/{pat}/*metrics_grid.png"],
        "P": [f"data/figures/reorganization/metric_correlation_*.png"],
        "Q": [f"data/figures/metastable/{pat}/*.png"],
        "R": [f"data/figures/summary/*.png"],
        "S": [f"data/figures/time_windows/{pat}/*.png"],
        "T": [f"data/figures/summary/*.png"],
        "U": [f"data/figures/spatial/{pat}/*.png"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Collect outputs into outputs/overleaf with a manifest."
    )
    parser.add_argument("--patient", default=None, help="Optional patient filter (Pat_XX)")
    parser.add_argument("--output-root", type=Path, default=Path("outputs/overleaf"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    output_root = args.output_root
    figures_dir = output_root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    manifest_lines: List[str] = ["# Overleaf Bundle Manifest", ""]

    for letter, pattern_list in _patterns(args.patient).items():
        manifest_lines.append(f"## {letter}")
        found_files: List[Path] = []
        for pattern in pattern_list:
            matches = [Path(p) for p in sorted(Path().glob(pattern))]
            found_files.extend(matches)
        if not found_files:
            manifest_lines.append("- status: missing")
            manifest_lines.append(f"- patterns: {', '.join(pattern_list)}")
            manifest_lines.append("")
            continue

        manifest_lines.append(f"- status: {len(found_files)} file(s)")
        manifest_lines.append(f"- patterns: {', '.join(pattern_list)}")
        manifest_lines.append("- files:")

        for src in found_files:
            dest_name = f"{letter}_{src.name}"
            dest = figures_dir / dest_name
            manifest_lines.append(f"  - {src} -> {dest}")
            if not args.dry_run:
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dest)

        manifest_lines.append("")

    manifest_path = output_root / "manifest.md"
    manifest_path.write_text("\n".join(manifest_lines) + "\n")

    if args.dry_run:
        print(f"Dry run: wrote manifest to {manifest_path}")
    else:
        print(f"Bundle ready: {output_root}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
