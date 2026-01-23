"""Time-window utilities for FC analysis."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

from lrg_eegfc.config.const import DEFAULT_SAMPLE_RATE

DEFAULT_WINDOW_OVERLAP = 0.25
DEFAULT_WINDOW_MIN_SEC = 10.0
DEFAULT_WINDOW_CYCLES = 10.0

DEFAULT_CORR_WINDOWS_CACHE_ROOT = Path("data/corr_cache_windows")
DEFAULT_MSC_WINDOWS_CACHE_ROOT = Path("data/msc_cache_windows")
DEFAULT_LRG_WINDOWS_CACHE_ROOT = Path("data/lrg_cache_windows")

DEFAULT_CORR_WINDOWS_DEV_CACHE_ROOT = Path("data/corr_cache_windows_dev")
DEFAULT_MSC_WINDOWS_DEV_CACHE_ROOT = Path("data/msc_cache_windows_dev")
DEFAULT_LRG_WINDOWS_DEV_CACHE_ROOT = Path("data/lrg_cache_windows_dev")

__all__ = [
    "DEFAULT_WINDOW_OVERLAP",
    "DEFAULT_WINDOW_MIN_SEC",
    "DEFAULT_WINDOW_CYCLES",
    "DEFAULT_CORR_WINDOWS_CACHE_ROOT",
    "DEFAULT_MSC_WINDOWS_CACHE_ROOT",
    "DEFAULT_LRG_WINDOWS_CACHE_ROOT",
    "DEFAULT_CORR_WINDOWS_DEV_CACHE_ROOT",
    "DEFAULT_MSC_WINDOWS_DEV_CACHE_ROOT",
    "DEFAULT_LRG_WINDOWS_DEV_CACHE_ROOT",
    "build_window_run_id",
    "get_window_cache_dir",
    "get_lrg_window_cache_dir",
    "suggest_window_sec",
    "compute_window_params",
    "generate_window_indices",
]


def _float_tag(value: float, digits: int = 3) -> str:
    tag = f"{value:.{digits}f}".rstrip("0").rstrip(".")
    return tag.replace(".", "p")


def _resolve_cache_root(cache_root: Optional[Path], default_root: Path, dev_root: Path, filter_time: Optional[int]) -> Path:
    root = cache_root or default_root
    if filter_time is not None and filter_time > 0 and root == default_root:
        return dev_root
    return root


def suggest_window_sec(
    low_freq: float,
    *,
    min_window_sec: float = DEFAULT_WINDOW_MIN_SEC,
    cycles: float = DEFAULT_WINDOW_CYCLES,
) -> float:
    """Suggest a window length that spans a minimum number of cycles."""
    if low_freq <= 0:
        return float(min_window_sec)
    return float(max(min_window_sec, cycles / low_freq))


def build_window_run_id(
    *,
    fc_method: str,
    window_sec: float,
    overlap: float,
    filter_time: Optional[int] = None,
    filter_type: str = "abs",
    zero_diagonal: bool = True,
    filter_order: int = 4,
    sparsify: str = "none",
    n_surrogates: int = 0,
    nperseg: int = 1024,
    noverlap: Optional[int] = None,
) -> str:
    parts = [
        f"winsec-{_float_tag(window_sec)}",
        f"ov-{_float_tag(overlap)}",
    ]

    if fc_method == "corr":
        parts.extend(
            [
                f"corr-ftype-{filter_type}",
                f"zdiag-{zero_diagonal}",
                f"forder-{filter_order}",
            ]
        )
    elif fc_method == "msc":
        actual_noverlap = nperseg // 2 if noverlap is None else int(noverlap)
        parts.extend(
            [
                f"msc-sparsify-{sparsify}",
                f"nsurr-{int(n_surrogates)}",
                f"nperseg-{int(nperseg)}",
                f"noverlap-{int(actual_noverlap)}",
            ]
        )
    else:
        raise ValueError(f"Unknown fc_method: {fc_method}")

    if filter_time is not None and filter_time > 0:
        parts.append(f"ftime-{int(filter_time)}")

    return "_".join(parts)


def get_window_cache_dir(
    *,
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    window_sec: float,
    overlap: float,
    filter_time: Optional[int] = None,
    cache_root: Optional[Path] = None,
    filter_type: str = "abs",
    zero_diagonal: bool = True,
    filter_order: int = 4,
    sparsify: str = "none",
    n_surrogates: int = 0,
    nperseg: int = 1024,
    noverlap: Optional[int] = None,
    ensure_dir: bool = True,
) -> Path:
    if fc_method == "corr":
        root = _resolve_cache_root(
            cache_root,
            DEFAULT_CORR_WINDOWS_CACHE_ROOT,
            DEFAULT_CORR_WINDOWS_DEV_CACHE_ROOT,
            filter_time,
        )
    else:
        root = _resolve_cache_root(
            cache_root,
            DEFAULT_MSC_WINDOWS_CACHE_ROOT,
            DEFAULT_MSC_WINDOWS_DEV_CACHE_ROOT,
            filter_time,
        )

    run_id = build_window_run_id(
        fc_method=fc_method,
        window_sec=window_sec,
        overlap=overlap,
        filter_time=filter_time,
        filter_type=filter_type,
        zero_diagonal=zero_diagonal,
        filter_order=filter_order,
        sparsify=sparsify,
        n_surrogates=n_surrogates,
        nperseg=nperseg,
        noverlap=noverlap,
    )
    run_dir = root / patient / phase / band / run_id
    if ensure_dir:
        run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def get_lrg_window_cache_dir(
    *,
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    window_sec: float,
    overlap: float,
    filter_time: Optional[int] = None,
    cache_root: Optional[Path] = None,
    filter_type: str = "abs",
    zero_diagonal: bool = True,
    filter_order: int = 4,
    sparsify: str = "none",
    n_surrogates: int = 0,
    nperseg: int = 1024,
    noverlap: Optional[int] = None,
    ensure_dir: bool = True,
) -> Path:
    root = _resolve_cache_root(
        cache_root,
        DEFAULT_LRG_WINDOWS_CACHE_ROOT,
        DEFAULT_LRG_WINDOWS_DEV_CACHE_ROOT,
        filter_time,
    )

    run_id = build_window_run_id(
        fc_method=fc_method,
        window_sec=window_sec,
        overlap=overlap,
        filter_time=filter_time,
        filter_type=filter_type,
        zero_diagonal=zero_diagonal,
        filter_order=filter_order,
        sparsify=sparsify,
        n_surrogates=n_surrogates,
        nperseg=nperseg,
        noverlap=noverlap,
    )
    run_dir = root / patient / phase / band / run_id
    if ensure_dir:
        run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def compute_window_params(
    n_samples: int,
    fs: float,
    window_sec: float,
    overlap: float,
) -> Tuple[int, int]:
    if fs <= 0:
        fs = DEFAULT_SAMPLE_RATE
    if window_sec <= 0:
        raise ValueError("window_sec must be > 0")
    if not (0 <= overlap < 1):
        raise ValueError("overlap must satisfy 0 <= overlap < 1")

    window_len = int(round(window_sec * fs))
    if window_len < 1:
        window_len = 1
    step = int(round(window_len * (1.0 - overlap)))
    if step < 1:
        step = 1
    if window_len > n_samples:
        raise ValueError(
            f"window_len {window_len} exceeds n_samples {n_samples}. "
            "Reduce window_sec or filter_time."
        )
    return window_len, step


def generate_window_indices(
    n_samples: int,
    window_len: int,
    step: int,
    *,
    max_windows: Optional[int] = None,
) -> Tuple[List[int], int]:
    starts = list(range(0, n_samples - window_len + 1, step))
    total = len(starts)
    if max_windows is not None:
        starts = starts[:max_windows]
    return starts, total
