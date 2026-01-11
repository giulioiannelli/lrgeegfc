#!/usr/bin/env python3
"""Inspect patient data and generate comprehensive report.

This script inspects all patient .mat files and metadata, generating a detailed
report about data completeness, sampling rates, file formats, and issues.

Usage:
    python src/inspect_patient_data.py
    python src/inspect_patient_data.py --output data_report.txt
    python src/inspect_patient_data.py --patients Pat_02 Pat_03
"""

import argparse
from datetime import date
from pathlib import Path

from lrg_eegfc.utils.io import (
    generate_csv_rows,
    generate_report,
    inspect_all_patients,
    save_csv,
    save_report,
)


def main():
    parser = argparse.ArgumentParser(
        description="Inspect patient data and generate report"
    )

    parser.add_argument(
        "--root-path",
        type=Path,
        default=Path("data/stereoeeg_patients"),
        help="Root directory containing patient data (default: data/stereoeeg_patients)"
    )
    parser.add_argument(
        "--patients",
        nargs="+",
        default=None,
        help="Specific patients to inspect (default: auto-detect all)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output report file (legacy single-file output)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(".agents/plans/active"),
        help="Directory for dated Markdown/CSV outputs (default: .agents/plans/active)"
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="Prefix for dated output files (default: <date>_data_inventory)"
    )
    parser.add_argument(
        "--print-only",
        action="store_true",
        help="Print report to console only (don't save file)"
    )

    args = parser.parse_args()

    print("Inspecting patient data...")
    print(f"Root path: {args.root_path}")

    # Run inspection
    results = inspect_all_patients(args.root_path, args.patients)

    # Generate report
    report = generate_report(results)

    # Output
    print(report)

    if not args.print_only:
        prefix = args.output_prefix or f\"{date.today().isoformat()}_data_inventory\"
        args.output_dir.mkdir(parents=True, exist_ok=True)
        md_path = args.output_dir / f\"{prefix}.md\"
        csv_path = args.output_dir / f\"{prefix}.csv\"

        save_report(results, md_path)
        save_csv(generate_csv_rows(results), csv_path)

        if args.output is not None:
            save_report(results, args.output)

        print(f\"\\n✓ Report saved to: {md_path}\")\n        print(f\"✓ CSV saved to: {csv_path}\")\n        if args.output is not None:\n            print(f\"✓ Legacy report saved to: {args.output}\")\n        print(f\"\\nSummary:\")\n        print(f\"  Total patients inspected: {results['summary']['total_patients']}\")\n        print(f\"  Patients with issues: {results['summary']['patients_with_issues']}\")\n        print(f\"  Total issues found: {results['summary']['total_issues']}\")


if __name__ == "__main__":
    main()
