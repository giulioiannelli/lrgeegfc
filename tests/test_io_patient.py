"""Smoke tests for lrg_eegfc.utils.io.patient promotions.

Focused on the helpers added 2026-05-28 (``PatientMasks`` +
``build_epi_masks``). The earlier functions (``load_timeseries``,
``load_channel_labels``, etc.) have integration coverage via
``test_data_inspect.py``.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from lrg_eegfc.utils.io.patient import (
    PatientMasks,
    _probe_of,
    build_epi_masks,
)


@pytest.fixture
def fake_patient(tmp_path: Path) -> Path:
    """Create a minimal patient layout with the channel_labels.csv stub."""
    pat_dir = tmp_path / "Pat_TEST"
    pat_dir.mkdir()
    labels = [
        '"label"\n',  # header row tolerated by load_channel_labels
        '"A1,G2"\n',
        '"A2,G2"\n',
        '"H 1,G2"\n',   # whitespace inside label gets stripped
        '"H 2,G2"\n',
        '"P\'1,G2"\n',  # apostrophe-suffixed probe (matches re pattern)
    ]
    (pat_dir / "channel_labels.csv").write_text("".join(labels))
    return tmp_path


def test_probe_of_alphabetic_prefix():
    assert _probe_of("A1") == "A"
    assert _probe_of("H12") == "H"
    assert _probe_of("P'5") == "P'"
    assert _probe_of("") == ""


def test_build_epi_masks_with_empty_epi_set(fake_patient: Path):
    with patch(
        "lrg_eegfc.utils.io.patient.load_epileptic_nodes",
        return_value=[],
    ):
        m = build_epi_masks("Pat_TEST", root_path=fake_patient)

    assert isinstance(m, PatientMasks)
    assert m.patient == "Pat_TEST"
    assert m.channels == ["A1", "A2", "H1", "H2", "P'1"]
    assert m.epi_mask.dtype == bool
    assert m.epi_mask.tolist() == [False] * 5
    assert list(m.probes) == ["A", "A", "H", "H", "P'"]


def test_build_epi_masks_marks_red_labels(fake_patient: Path):
    with patch(
        "lrg_eegfc.utils.io.patient.load_epileptic_nodes",
        return_value=["A1", "H1"],
    ):
        m = build_epi_masks("Pat_TEST", root_path=fake_patient)

    assert m.epi_mask.tolist() == [True, False, True, False, False]


def test_build_epi_masks_string_membership_no_silent_all_false(fake_patient: Path):
    """Audit_71 bug regression: np.isin(int, str) silently returned all-False.

    The canonical helper must use string set membership directly.
    """
    with patch(
        "lrg_eegfc.utils.io.patient.load_epileptic_nodes",
        return_value=["A2", "P'1"],
    ):
        m = build_epi_masks("Pat_TEST", root_path=fake_patient)

    assert m.epi_mask.sum() == 2
    assert m.channels[int(np.where(m.epi_mask)[0][0])] == "A2"
    assert m.channels[int(np.where(m.epi_mask)[0][1])] == "P'1"


def test_build_epi_masks_unknown_red_label_does_not_error(fake_patient: Path):
    with patch(
        "lrg_eegfc.utils.io.patient.load_epileptic_nodes",
        return_value=["XX99"],  # label not present in channels
    ):
        m = build_epi_masks("Pat_TEST", root_path=fake_patient)

    assert m.epi_mask.sum() == 0
    assert len(m.channels) == 5
