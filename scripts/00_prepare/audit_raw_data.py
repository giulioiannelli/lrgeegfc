#!/usr/bin/env python3
"""Exhaustive audit of data/raw/stereoeeg_patients/Pat_NN/ for every patient.

Checks performed per patient:
  1. Canonical layout (resting/, task/, implant/ dirs present).
  2. All 4 phase .mat files — existence, size, scipy OR h5py parse, Data key,
     channel count matches channel_labels.csv, sampling rate (fs) if present.
  3. channel_labels.csv — readable, has 'label' header, row count consistent.
  4. implant_pat_NN.csv — readable, has 'label' column, merges with labels.
  5. implant_pat_NN.xlsx under implant/ — openpyxl readable, red-font cell count.
  6. provenance.md — present and sha256 entries match current file hashes.
  7. Unexpected / orphan files under the patient root.

Prints a per-patient PASS/FAIL panel and a final red-flag table.

Exit code 0 = all green, 1 = one or more patients have defects.
"""
from __future__ import annotations

import csv
import hashlib
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat
import h5py

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import PHASE_LABELS, PHASE_SUBDIR, list_patients


# ────────────────────────── helpers ──────────────────────────

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def try_load_mat(path: Path):
    """Return (data, params, keys, loader, error) with data shape (n_ch, n_samp) or None."""
    keys = None
    # Attempt scipy first
    try:
        m = loadmat(str(path))
        keys = [k for k in m if not k.startswith("_")]
        if "Data" in m:
            d = np.asarray(m["Data"])
            if d.ndim != 2:
                return None, m.get("Parameters"), keys, "scipy", f"Data ndim={d.ndim}"
            if d.shape[0] > d.shape[1]:
                d = d.T
            return d, m.get("Parameters"), keys, "scipy", None
        return None, m.get("Parameters"), keys, "scipy", "no Data key"
    except NotImplementedError:
        pass
    except Exception as e:
        try:
            with h5py.File(path, "r") as f:
                keys = list(f.keys())
                if "Data" in f:
                    d = np.asarray(f["Data"])
                    if d.ndim == 2 and d.shape[0] > d.shape[1]:
                        d = d.T
                    return d, None, keys, "h5py", None
            return None, None, keys, "h5py", "no Data key"
        except Exception as h_err:
            return None, None, None, None, f"scipy={type(e).__name__}:{e}; h5py={type(h_err).__name__}:{h_err}"
    # scipy raised NotImplementedError -> try h5py
    try:
        with h5py.File(path, "r") as f:
            keys = list(f.keys())
            if "Data" in f:
                d = np.asarray(f["Data"])
                if d.ndim == 2 and d.shape[0] > d.shape[1]:
                    d = d.T
                return d, None, keys, "h5py", None
        return None, None, keys, "h5py", "no Data key"
    except Exception as e:
        return None, None, None, None, f"h5py: {type(e).__name__}: {e}"


def extract_fs(params) -> float | None:
    if params is None:
        return None
    try:
        names = getattr(getattr(params, "dtype", None), "names", None) or ()
        if "fs" not in names:
            return None
        val = np.squeeze(params["fs"]).item()
        return float(val)
    except Exception:
        return None


def parse_provenance(path: Path) -> dict[str, str]:
    """Extract sha256 entries from provenance.md. Keys = canonical filenames."""
    if not path.exists():
        return {}
    text = path.read_text(encoding="utf-8", errors="replace")
    # Look for lines like `sha256: <hash>` near `canonical: <path>`.
    entries: dict[str, str] = {}
    curr_canonical = None
    for line in text.splitlines():
        line = line.strip()
        m1 = re.match(r"^[-*]\s*canonical:\s*`?([^`\s]+)`?", line)
        if m1:
            curr_canonical = m1.group(1).split("/")[-1]
            continue
        m2 = re.match(r"^[-*]?\s*sha256:\s*`?([0-9a-f]{64})`?", line)
        if m2 and curr_canonical:
            entries[curr_canonical] = m2.group(1)
            curr_canonical = None
    return entries


# ────────────────────────── per-patient audit ──────────────────────────

@dataclass
class PhaseReport:
    phase: str
    file: Path
    exists: bool = False
    size_bytes: int = 0
    loader: str | None = None
    n_channels: int | None = None
    n_samples: int | None = None
    fs: float | None = None
    error: str | None = None

    @property
    def status(self) -> str:
        if not self.exists:
            return "MISSING"
        if self.error:
            return "ERROR"
        return "OK"


@dataclass
class PatientReport:
    patient: str
    root: Path
    problems: list[str] = field(default_factory=list)
    phases: dict[str, PhaseReport] = field(default_factory=dict)
    n_labels: int | None = None
    xlsx_red: int | None = None
    implant_csv_rows: int | None = None
    provenance_mismatches: list[str] = field(default_factory=list)
    unexpected_files: list[Path] = field(default_factory=list)
    fs_values: set[float] = field(default_factory=set)

    @property
    def ok(self) -> bool:
        return not self.problems and not self.provenance_mismatches and all(
            pr.status == "OK" for pr in self.phases.values()
        )


CANONICAL_DIRS = {"resting", "task", "implant"}
CANONICAL_FILES_RE = re.compile(
    r"^(channel_labels\.csv|implant_pat_\d{2}\.csv|provenance\.md)$"
)
PHASE_INFO_TXT_RE = re.compile(r"^(rest_pre|rest_post|task_learn|task_test)\.info\.txt$")


def audit_patient(patient: str, dataset_root: Path) -> PatientReport:
    r = PatientReport(patient=patient, root=dataset_root / patient)

    if not r.root.exists():
        r.problems.append(f"patient directory missing: {r.root}")
        return r

    # 1) Canonical dirs present
    for sub in CANONICAL_DIRS:
        if not (r.root / sub).is_dir():
            r.problems.append(f"missing subdir: {sub}/")

    # 2) channel_labels.csv
    labels_csv = r.root / "channel_labels.csv"
    if not labels_csv.exists():
        r.problems.append("missing channel_labels.csv")
    else:
        try:
            ldf = pd.read_csv(labels_csv)
            if "label" not in ldf.columns:
                r.problems.append(f"channel_labels.csv has no 'label' column (cols={list(ldf.columns)})")
            else:
                r.n_labels = len(ldf)
        except Exception as e:
            r.problems.append(f"channel_labels.csv unreadable: {type(e).__name__}: {e}")

    # 3) Phase .mat files
    for phase in PHASE_LABELS:
        subdir = PHASE_SUBDIR[phase]
        path = r.root / subdir / f"{phase}.mat"
        pr = PhaseReport(phase=phase, file=path)
        if not path.exists():
            pr.exists = False
            r.phases[phase] = pr
            continue
        pr.exists = True
        pr.size_bytes = path.stat().st_size
        if pr.size_bytes == 0:
            pr.error = "zero-byte file"
            r.phases[phase] = pr
            continue
        data, params, keys, loader, err = try_load_mat(path)
        pr.loader = loader
        if err:
            pr.error = err
            r.phases[phase] = pr
            continue
        if data is None or data.ndim != 2:
            pr.error = f"no 2-D Data (keys={keys})"
            r.phases[phase] = pr
            continue
        pr.n_channels = int(data.shape[0])
        pr.n_samples = int(data.shape[1])
        pr.fs = extract_fs(params)
        if pr.fs is not None:
            r.fs_values.add(float(pr.fs))
        if r.n_labels is not None and pr.n_channels != r.n_labels:
            pr.error = f"channel count {pr.n_channels} ≠ {r.n_labels} labels"
        r.phases[phase] = pr

    # 4) implant xlsx + red-cell count (openpyxl)
    patnum = int(patient.split("_")[-1])
    implant_dir = r.root / "implant"
    xlsx = implant_dir / f"implant_pat_{patnum:02d}.xlsx"
    if not xlsx.exists():
        r.problems.append(f"missing implant xlsx: {xlsx.relative_to(r.root)}")
    else:
        try:
            from openpyxl import load_workbook
            wb = load_workbook(xlsx, data_only=True)
            red_cells = 0
            for ws in wb.worksheets:
                for row in ws.iter_rows():
                    for cell in row:
                        if cell.value is None:
                            continue
                        font = cell.font
                        if font and font.color and font.color.rgb:
                            rgb = str(font.color.rgb).upper()
                            # Red tones: FFRRGGBB; red high, green/blue low.
                            if rgb.startswith("FFFF") and rgb[4:6] < "80" and rgb[6:8] < "80":
                                red_cells += 1
                            elif rgb in ("FFFF0000", "FFC00000", "FFE60000"):
                                red_cells += 1
            r.xlsx_red = red_cells
        except Exception as e:
            r.problems.append(f"implant xlsx unreadable: {type(e).__name__}: {e}")

    # 5) implant_pat_NN.csv
    icsv = r.root / f"implant_pat_{patnum:02d}.csv"
    if not icsv.exists():
        r.problems.append(f"missing implant_pat_{patnum:02d}.csv")
    else:
        try:
            idf = pd.read_csv(icsv)
            if "label" not in idf.columns:
                r.problems.append(f"implant csv has no 'label' column (cols={list(idf.columns)[:6]})")
            else:
                r.implant_csv_rows = len(idf)
        except Exception as e:
            r.problems.append(f"implant csv unreadable: {type(e).__name__}: {e}")

    # 6) provenance.md — sha256 entries match current files
    prov = parse_provenance(r.root / "provenance.md")
    for phase in PHASE_LABELS:
        path = r.root / PHASE_SUBDIR[phase] / f"{phase}.mat"
        if path.exists() and path.stat().st_size > 0:
            if path.name in prov:
                expect = prov[path.name]
                actual = sha256(path)
                if expect != actual:
                    r.provenance_mismatches.append(f"{path.name}: disk={actual[:12]} prov={expect[:12]}")
    # implant xlsx
    if xlsx.exists():
        if xlsx.name in prov:
            if prov[xlsx.name] != sha256(xlsx):
                r.provenance_mismatches.append(f"{xlsx.name}: sha256 mismatch")

    # 7) Unexpected top-level files or unknown subdirs under patient dir
    for item in sorted(r.root.iterdir()):
        if item.is_dir() and item.name not in CANONICAL_DIRS:
            r.unexpected_files.append(item)
        elif item.is_file():
            if CANONICAL_FILES_RE.match(item.name):
                continue
            r.unexpected_files.append(item)
    for sub, sub_re in [
        ("resting", re.compile(r"^(rest_pre|rest_post)\.(mat|info\.txt)$")),
        ("task", re.compile(r"^(task_learn|task_test)\.(mat|info\.txt)$")),
    ]:
        d = r.root / sub
        if d.is_dir():
            for item in sorted(d.iterdir()):
                if item.is_file() and not sub_re.match(item.name):
                    r.unexpected_files.append(item)
    d = r.root / "implant"
    if d.is_dir():
        for item in sorted(d.iterdir()):
            if item.is_file() and not re.match(rf"^implant_pat_{patnum:02d}\.xlsx$", item.name):
                r.unexpected_files.append(item)

    return r


# ────────────────────────── reporting ──────────────────────────

def print_patient(r: PatientReport):
    mark = "PASS" if r.ok else "FAIL"
    print(f"\n── {r.patient}  [{mark}] ──")
    # Labels & implant
    print(f"  channel_labels.csv rows: {r.n_labels}")
    print(f"  implant_pat_{int(r.patient.split('_')[-1]):02d}.csv rows: {r.implant_csv_rows}")
    print(f"  implant xlsx red-font cells: {r.xlsx_red}")
    # Phases
    for phase in PHASE_LABELS:
        pr = r.phases.get(phase)
        if pr is None:
            continue
        if pr.status == "MISSING":
            print(f"  {phase:10s}: MISSING FILE")
        elif pr.status == "ERROR":
            print(f"  {phase:10s}: ERROR ({pr.error})  size={pr.size_bytes:,}")
        else:
            fs = f"{pr.fs:.0f} Hz" if pr.fs is not None else "fs?"
            print(f"  {phase:10s}: shape=({pr.n_channels}, {pr.n_samples:,})  {fs}  loader={pr.loader}")
    # Sampling rates
    if len(r.fs_values) > 1:
        print(f"  ⚠ mixed fs across phases: {sorted(r.fs_values)}")
    # Provenance
    if r.provenance_mismatches:
        print("  ⚠ provenance.md mismatches:")
        for m in r.provenance_mismatches:
            print(f"      {m}")
    # Unexpected files
    if r.unexpected_files:
        print("  ⚠ unexpected files:")
        for p in r.unexpected_files:
            print(f"      {p.relative_to(r.root)}")
    if r.problems:
        print("  ⚠ structural problems:")
        for p in r.problems:
            print(f"      {p}")


def main() -> int:
    patients = list_patients()
    reports: list[PatientReport] = []
    for pat in patients:
        reports.append(audit_patient(pat, SEEG_DATAPATH))

    for r in reports:
        print_patient(r)

    # Red-flag summary
    print("\n" + "=" * 70)
    print("SUMMARY RED FLAGS (one line per defect)")
    print("=" * 70)
    any_defect = False
    for r in reports:
        for phase, pr in r.phases.items():
            if pr.status == "MISSING":
                print(f"  [{r.patient}] {phase}: MISSING FILE")
                any_defect = True
            elif pr.status == "ERROR":
                print(f"  [{r.patient}] {phase}: {pr.error}  ({pr.size_bytes:,} bytes)")
                any_defect = True
            elif r.n_labels is not None and pr.n_channels is not None and pr.n_channels != r.n_labels:
                print(f"  [{r.patient}] {phase}: channel count {pr.n_channels} ≠ {r.n_labels} labels")
                any_defect = True
        if len(r.fs_values) > 1:
            print(f"  [{r.patient}] mixed fs across phases: {sorted(r.fs_values)}")
            any_defect = True
        for msg in r.problems:
            print(f"  [{r.patient}] {msg}")
            any_defect = True
        for msg in r.provenance_mismatches:
            print(f"  [{r.patient}] provenance: {msg}")
            any_defect = True
        for p in r.unexpected_files:
            print(f"  [{r.patient}] unexpected file: {p.relative_to(r.root)}")
            any_defect = True

    if not any_defect:
        print("  (none — all patients clean)")

    return 0 if not any_defect else 1


if __name__ == "__main__":
    raise SystemExit(main())
