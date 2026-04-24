from pathlib import Path

import numpy as np
from scipy.io import savemat

from lrg_eegfc.utils.io import (
    generate_csv_rows,
    inspect_all_patients,
    inspect_patient,
)
from lrg_eegfc.config.const import PHASE_LABELS, PHASE_SUBDIR


def _write_mat(path: Path) -> None:
    data = np.random.randn(2, 10)
    savemat(path, {"Data": data, "fs": 2048.0})


def _setup_patient(root: Path, patient: str) -> None:
    patient_dir = root / patient
    patient_dir.mkdir(parents=True, exist_ok=True)
    for phase in PHASE_LABELS:
        subdir = patient_dir / PHASE_SUBDIR[phase]
        subdir.mkdir(parents=True, exist_ok=True)
        _write_mat(subdir / f"{phase}.mat")

    (patient_dir / "channel_labels.csv").write_text("label\nA1\n", encoding="utf-8")

    patnum = int(patient.split("_")[-1])
    implant_path = patient_dir / f"implant_pat_{patnum:02d}.csv"
    implant_path.write_text("label,x,y,z\nA1,0,0,0\n", encoding="utf-8")


def test_inspect_patient_detects_metadata(tmp_path: Path) -> None:
    root = tmp_path / "data"
    _setup_patient(root, "Pat_01")

    result = inspect_patient("Pat_01", root)
    assert result["directory_exists"] is True
    assert result["has_channel_labels"] is True
    assert result["has_implant"] is True
    assert result["phases"]["rest_pre"]["data_variable"] == "Data"
    assert result["phases"]["rest_pre"]["fs"] == 2048.0


def test_generate_csv_rows(tmp_path: Path) -> None:
    root = tmp_path / "data"
    _setup_patient(root, "Pat_02")

    results = inspect_all_patients(root, ["Pat_02"])
    rows = generate_csv_rows(results)

    assert rows, "Expected CSV rows to be generated"
    row = rows[0]
    assert row["patient"] == "Pat_02"
    assert row["file_exists"] is True
    assert row["data_variable"] == "Data"
    assert row["has_implant"] is True
