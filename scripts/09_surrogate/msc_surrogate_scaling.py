#!/usr/bin/env python3
"""Run an MSC surrogate scaling sweep and report stability metrics."""

from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter
from typing import Iterable, List

import numpy as np
from scipy.stats import spearmanr

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.utils.io import list_patients
from lrg_eegfc.workflow.msc import compute_msc_matrix, get_msc_cache_path


def _parse_int_list(value: str) -> List[int]:
    return [int(v.strip()) for v in value.split(",") if v.strip()]


def _upper_triangle(mat: np.ndarray) -> np.ndarray:
    idx = np.triu_indices_from(mat, k=1)
    return mat[idx]


def _safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0 or b.size == 0:
        return float("nan")
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    try:
        return float(spearmanr(a, b).correlation)
    except Exception:
        return float("nan")


def _select_default_patient(dataset_root: Path) -> str:
    patients = list_patients(dataset_root)
    if not patients:
        raise FileNotFoundError(f"No patients found under {dataset_root}")
    return patients[0]


def run_sweep(
    patient: str,
    phase: str,
    band: str,
    surrogates: Iterable[int],
    *,
    dataset_root: Path,
    cache_root: Path,
    nperseg: int,
    filter_time: int | None,
    overwrite: bool,
    use_cache: bool,
    verbose: bool,
) -> List[dict]:
    results = []
    matrices = []
    surrogates = list(surrogates)

    for idx, ns in enumerate(surrogates, start=1):
        sparsify = "soft" if ns > 0 else "none"
        cache_path = get_msc_cache_path(
            patient,
            phase,
            band,
            cache_root=cache_root,
            sparsify=sparsify,
            n_surrogates=ns,
            nperseg=nperseg,
            filter_time=filter_time,
        )
        cache_hit = cache_path.exists() and use_cache and not overwrite

        if verbose:
            print(
                f"[{idx}/{len(surrogates)}] n_surrogates={ns} "
                f"sparsify={sparsify} nperseg={nperseg} cache_hit={cache_hit}"
            )
        start = perf_counter()
        res = compute_msc_matrix(
            patient,
            phase,
            band,
            dataset_root=dataset_root,
            cache_root=cache_root,
            sparsify=sparsify,
            n_surrogates=ns,
            nperseg=nperseg,
            filter_time=filter_time,
            overwrite_cache=overwrite,
            use_cache=use_cache,
            verbose=False,
        )
        elapsed = perf_counter() - start
        if verbose:
            print(f"    done in {elapsed:.2f}s; matrix shape={res.adjacency_matrix.shape}")

        matrices.append(res.adjacency_matrix)
        results.append(
            {
                "patient": patient,
                "phase": phase,
                "band": band,
                "nperseg": nperseg,
                "filter_time": filter_time,
                "n_surrogates": ns,
                "sparsify": sparsify,
                "time_s": elapsed,
                "cache_hit": cache_hit,
            }
        )

    ref_vec = _upper_triangle(matrices[-1])
    ref_norm = np.linalg.norm(ref_vec)

    for result, mat in zip(results, matrices):
        vec = _upper_triangle(mat)
        result["pearson_vs_max"] = _safe_pearson(ref_vec, vec)
        result["spearman_vs_max"] = _safe_spearman(ref_vec, vec)
        result["rel_frobenius"] = float(np.linalg.norm(ref_vec - vec) / ref_norm) if ref_norm else float("nan")

    return results


def _write_csv(rows: List[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    lines = [",".join(keys)]
    for row in rows:
        lines.append(",".join(str(row[k]) for k in keys))
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run MSC surrogate scaling sweep and report stability metrics."
    )
    parser.add_argument("--patient", default=None, help="Patient ID (default: first found)")
    parser.add_argument("--phase", default=PHASE_LABELS[0], help="Phase label")
    parser.add_argument("--band", default=next(iter(BRAIN_BANDS.keys())), help="Band name")
    parser.add_argument(
        "--surrogates",
        default="0,25,50,100,200,500",
        help="Comma-separated surrogate counts",
    )
    parser.add_argument("--nperseg", type=int, default=1024, help="Welch nperseg")
    parser.add_argument(
        "--filter-time",
        type=int,
        default=None,
        help="Limit samples (dev-only); omit for full-length data",
    )
    from lrg_eegfc.config.paths import MSC_CACHE, SEEG_DATAPATH
    parser.add_argument("--dataset-root", type=Path, default=SEEG_DATAPATH)
    parser.add_argument("--cache-root", type=Path, default=MSC_CACHE)
    parser.add_argument("--overwrite", action="store_true", help="Recompute even if cached")
    parser.add_argument("--no-cache", action="store_true", help="Disable cache reads")
    parser.add_argument("--output", type=Path, default=None, help="Optional CSV output path")
    parser.add_argument("--verbose", action="store_true", help="Print progress updates")

    args = parser.parse_args()

    patient = args.patient or _select_default_patient(args.dataset_root)
    surrogates = _parse_int_list(args.surrogates)

    rows = run_sweep(
        patient,
        args.phase,
        args.band,
        surrogates,
        dataset_root=args.dataset_root,
        cache_root=args.cache_root,
        nperseg=args.nperseg,
        filter_time=args.filter_time,
        overwrite=args.overwrite,
        use_cache=not args.no_cache,
        verbose=args.verbose,
    )

    if args.verbose:
        print("\nSummary:")
    print("n_surrogates,time_s,pearson_vs_max,spearman_vs_max,rel_frobenius,cache_hit")
    for row in rows:
        print(
            f"{row['n_surrogates']},{row['time_s']:.3f},{row['pearson_vs_max']:.4f},"
            f"{row['spearman_vs_max']:.4f},{row['rel_frobenius']:.4f},{row['cache_hit']}"
        )

    if args.output:
        _write_csv(rows, args.output)
        print(f"\nSaved CSV: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
