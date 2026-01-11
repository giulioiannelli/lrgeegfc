from pathlib import Path

from lrg_eegfc.config.const import list_patients, PATIENTS_LIST


def test_list_patients_from_custom_root(tmp_path: Path) -> None:
    (tmp_path / "Pat_02").mkdir()
    (tmp_path / "Pat_01").mkdir()
    (tmp_path / "NotAPatient").mkdir()

    patients = list_patients(tmp_path)

    assert patients == ["Pat_01", "Pat_02"]


def test_patients_list_lazy_type() -> None:
    assert isinstance(PATIENTS_LIST, list)
