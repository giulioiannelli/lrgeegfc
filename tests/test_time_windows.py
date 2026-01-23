from lrg_eegfc.workflow.time_windows import (
    build_window_run_id,
    compute_window_params,
    generate_window_indices,
    suggest_window_sec,
)


def test_suggest_window_sec():
    assert suggest_window_sec(10.0, min_window_sec=10.0, cycles=10.0) == 10.0
    assert suggest_window_sec(0.5, min_window_sec=10.0, cycles=10.0) == 20.0


def test_compute_window_params():
    window_len, step = compute_window_params(1000, 100.0, 1.0, 0.25)
    assert window_len == 100
    assert step == 75


def test_generate_window_indices():
    starts, total = generate_window_indices(1000, 100, 75)
    assert total == 13
    assert starts[0] == 0
    assert starts[-1] == 900


def test_build_window_run_id_corr():
    run_id = build_window_run_id(
        fc_method="corr",
        window_sec=10.0,
        overlap=0.25,
        filter_type="abs",
        zero_diagonal=True,
        filter_order=4,
    )
    assert "winsec-10" in run_id
    assert "ov-0p25" in run_id
    assert "corr-ftype-abs" in run_id
    assert "zdiag-True" in run_id
    assert "forder-4" in run_id


def test_build_window_run_id_msc():
    run_id = build_window_run_id(
        fc_method="msc",
        window_sec=10.0,
        overlap=0.25,
        sparsify="soft",
        n_surrogates=200,
        nperseg=1024,
        noverlap=512,
    )
    assert "msc-sparsify-soft" in run_id
    assert "nsurr-200" in run_id
    assert "nperseg-1024" in run_id
    assert "noverlap-512" in run_id
