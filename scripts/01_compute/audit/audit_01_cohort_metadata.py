#!/usr/bin/env python3
"""Audit Step 1 — per-patient cohort metadata for the N=9 diagnostic.

Writes ``data/audit/cohort_metadata.csv`` with columns:

    patient_id, sampling_rate_Hz, n_channels, n_seconds_rsPre,
    n_seconds_taskLearn, n_seconds_taskTest, n_seconds_rsPost,
    is_core_N5, notes

Missing phases -> NaN. Pat_03 flagged as 1024 Hz in `notes`.
"""
from __future__ import annotations

import csv

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import FS_OVERRIDES, PATIENTS_4PHASE, PHASE_LABELS
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io import load_timeseries

CORE_N5 = {"Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"}
OUT = ROOT / "data" / "audit" / "cohort_metadata.csv"


def _duration_seconds(patient: str, phase: str, fs: float) -> float:
    try:
        ts = load_timeseries(patient, phase, SEEG_DATAPATH)
    except Exception as exc:  # noqa: BLE001
        print(f"  {patient} {phase}: load failed ({exc})")
        return np.nan
    # ts shape is (n_channels, n_samples) after the loader's transpose fix
    n_samples = ts.shape[1]
    return n_samples / fs


def main() -> None:
    rows = []
    for pat in PATIENTS_4PHASE:
        fs = FS_OVERRIDES.get(pat, 2048.0)
        print(f"{pat}  fs={fs:.0f} Hz")
        durations = {}
        n_channels = None
        for phase in PHASE_LABELS:
            d = _duration_seconds(pat, phase, fs)
            durations[phase] = d
            if n_channels is None and not np.isnan(d):
                try:
                    ts = load_timeseries(pat, phase, SEEG_DATAPATH)
                    n_channels = int(ts.shape[0])
                except Exception:  # noqa: BLE001
                    pass

        notes = []
        if pat == "Pat_03":
            notes.append("1024 Hz outlier")
        if pat == "Pat_10":
            notes.append("113-ch via PATIENT_CHANNEL_DROP")
        if pat == "Pat_14":
            notes.append("task_test vendor-corrupt -> excluded from cross-phase")
        if pat == "Pat_13":
            notes.append("rest_pre vendor-replaced 2026-04-23")

        rows.append({
            "patient_id": pat,
            "sampling_rate_Hz": fs,
            "n_channels": n_channels if n_channels is not None else np.nan,
            "n_seconds_rsPre": durations["rest_pre"],
            "n_seconds_taskLearn": durations["task_learn"],
            "n_seconds_taskTest": durations["task_test"],
            "n_seconds_rsPost": durations["rest_post"],
            "is_core_N5": pat in CORE_N5,
            "notes": "; ".join(notes),
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} rows -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
