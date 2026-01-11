import os
from pathlib import Path

import numpy as np
import pytest

from lrg_eegfc.workflow.msc import compute_msc_matrix
from lrg_eegfc.config.const import PHASE_LABELS, BRAIN_BANDS


if os.getenv("LRG_EEGFC_SURROGATE_TEST") != "1":
    pytest.skip(
        "Set LRG_EEGFC_SURROGATE_TEST=1 to run surrogate scaling tests",
        allow_module_level=True,
    )

DATA_ROOT = Path("data/stereoeeg_patients")
if not DATA_ROOT.exists():
    pytest.skip("Dataset not available", allow_module_level=True)


def _pick_patient() -> str:
    patients = sorted(p.name for p in DATA_ROOT.iterdir() if p.is_dir() and p.name.startswith("Pat_"))
    if not patients:
        pytest.skip("No patient folders found", allow_module_level=True)
    return patients[0]


@pytest.mark.slow
def test_msc_surrogate_scaling_smoke() -> None:
    patient = _pick_patient()
    phase = PHASE_LABELS[0]
    band = next(iter(BRAIN_BANDS.keys()))

    results = []
    for n_surrogates in (0, 25):
        result = compute_msc_matrix(
            patient,
            phase,
            band,
            sparsify="soft" if n_surrogates > 0 else "none",
            n_surrogates=n_surrogates,
            filter_time=5000,
            overwrite_cache=False,
            use_cache=True,
            verbose=False,
        )
        results.append(result.adjacency_matrix)

    assert results[0].shape == results[1].shape
    assert np.isfinite(results[0]).all()
    assert np.isfinite(results[1]).all()
