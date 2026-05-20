#!/usr/bin/env python3
"""Audit 60 — anchor anatomy baseline + mutual-exclusivity tie-breaker.

For each (patient, band) cell, compute the fraction of intra-module
*pairs* that share the same sEEG probe (probe = leading letters of the
channel label, e.g. ``A1`` and ``A2`` share probe ``A``) for the
anchor and trace classes, and compare against the cohort same-probe
pair fraction (the per-patient baseline given the implant geometry
alone). Anchor's question: are anchor modules driven by probe
geometry (anatomy-anchored) or by genuine functional anchoring across
phases?

Strict-intersection leaf sets follow the §5.6 definitions:
    trace  : ``leaves_tt ∩ leaves_post``       (matched pair, two-way)
    anchor : ``leaves_pre ∩ leaves_tt ∩ leaves_post`` (triple match)

Implements the **mutual-exclusivity tie-breaker** that the §5.6
manuscript needs: leaves classified into a tighter class are removed
from the candidate pool of looser classes before each pass runs. The
canonical class-ordering is::

    anchor (tightest gate, three-way Jaccard ≥ 0.6)
    trace  (two-way matched pair, fragmentation in rPre)
    reset  (two-way matched pair, fragmentation in tt)
    rearrange (residual rPost-only)
    diffuse  (leftover)

Outputs
-------
``data/audit/anchor_anatomy_baseline/per_module_same_probe_fraction.csv``
``data/audit/anchor_anatomy_baseline/cohort_aggregate.csv``

Cohort
------
n=10 — Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15. Pat_03 outlier
(1024 Hz) included with the rest.

Scope: ``.agents/guides/task-persistence-investigation/2026-05-08_taxonomy_mutual_exclusivity.md``
"""
from __future__ import annotations

import sys
import traceback
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.io.patient import load_channel_labels, parse_seeg_label
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

# Reuse the canonical detectors from the per-class audits — no
# private re-implementations.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_47_kc_trace_network_view import find_trace_pairs
from audit_48_kc_reset_module_view import find_reset_pairs
from audit_49_kc_rearrangement_module_view import find_rearrangement_modules
from audit_50_kc_anchor_module_view import find_anchor_triples


PATIENTS = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
PHASES = ("rest_pre", "task_test", "rest_post")

OUT = ROOT / "data" / "audit" / "anchor_anatomy_baseline"


# ---------------------------------------------------------------------------
# Probe geometry
# ---------------------------------------------------------------------------

def derive_probe_ids(patient: str) -> list[str]:
    """Return per-channel probe label (e.g. ``A``, ``H'``) using the
    canonical loader and parser. Unparseable labels (rare scalp-EEG
    contacts) keep the cleaned label itself as a singleton probe id."""
    labels = load_channel_labels(patient)
    probes: list[str] = []
    for lab in labels:
        probe, _ = parse_seeg_label(lab)
        probes.append(probe if probe is not None else lab)
    return probes


def cohort_same_probe_pair_fraction(probe_ids: list[str]) -> float:
    """Baseline fraction of same-probe pairs among all unordered pairs.

    Equivalent to ``Σ_p binom(n_p, 2) / binom(N, 2)`` where ``n_p`` is
    the number of contacts on probe ``p``.
    """
    n = len(probe_ids)
    if n < 2:
        return 0.0
    counts: dict[str, int] = {}
    for p in probe_ids:
        counts[p] = counts.get(p, 0) + 1
    same = sum(c * (c - 1) // 2 for c in counts.values())
    total = n * (n - 1) // 2
    return same / total if total else 0.0


def module_same_probe_pair_stats(
    leaves: set[int], probe_ids: list[str]
) -> tuple[int, int, int]:
    """Return ``(n_pairs, n_same_probe_pairs, n_cross_probe_pairs)`` for
    the unordered pairs of contacts inside ``leaves``."""
    if len(leaves) < 2:
        return 0, 0, 0
    leaves_list = sorted(leaves)
    same = 0
    total = 0
    for i, j in combinations(leaves_list, 2):
        total += 1
        if probe_ids[i] == probe_ids[j]:
            same += 1
    return total, same, total - same


# ---------------------------------------------------------------------------
# Strict-intersection leaf sets (§5.6 definitions)
# ---------------------------------------------------------------------------

def strict_anchor_leaves(m: dict) -> set[int]:
    return set(m["leaves_pre"]) & set(m["leaves_tt"]) & set(m["leaves_post"])


def strict_trace_leaves(m: dict) -> set[int]:
    return set(m["leaves_tt"]) & set(m["leaves_post"])


def strict_reset_leaves(m: dict) -> set[int]:
    return set(m["leaves_pre"]) & set(m["leaves_post"])


def strict_rearrange_leaves(m: dict) -> set[int]:
    # rearrange is rPost-anchored; in §5.6 the residual is computed
    # at the assemble step. Here we keep the raw rPost leaves and let
    # the tie-breaker subtract the higher-priority leaves outside.
    return set(m["leaves_post"])


# ---------------------------------------------------------------------------
# Mutual-exclusivity tie-breaker (anchor > trace > reset > rearrange)
# ---------------------------------------------------------------------------

def apply_tie_breaker(modules_by_class: dict, n: int) -> dict:
    """Return per-class strict leaf sets after subtracting higher-priority
    leaves. Class-ordering: anchor > trace > reset > rearrange.

    The output is a list of ``(module_id, residual_leaves)`` per class
    in the order modules are returned by the per-class detectors.
    Modules whose residual size drops below ``MIN_SIZE = 3`` are
    dropped; their leaves fall through to lower-priority classes.
    """
    MIN_RESIDUAL = 3

    used: set[int] = set()
    out: dict[str, list[tuple[int, set[int]]]] = {}

    for cls in ("anchor", "trace", "reset", "rearrange"):
        out[cls] = []
        for idx, m in enumerate(modules_by_class.get(cls, [])):
            if cls == "anchor":
                strict = strict_anchor_leaves(m)
            elif cls == "trace":
                strict = strict_trace_leaves(m)
            elif cls == "reset":
                strict = strict_reset_leaves(m)
            else:
                strict = strict_rearrange_leaves(m)
            residual = strict - used
            if len(residual) >= MIN_RESIDUAL:
                out[cls].append((idx, residual))
                used |= residual

    diffuse_leaves = set(range(n)) - used
    out["diffuse"] = (
        [(0, diffuse_leaves)] if len(diffuse_leaves) >= MIN_RESIDUAL else []
    )
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []

    for patient in PATIENTS:
        try:
            probe_ids = derive_probe_ids(patient)
        except Exception as e:
            print(f"[audit_60] WARN {patient}: cannot load channel labels: {e}")
            traceback.print_exc()
            continue

        for band in BANDS:
            try:
                Z_by_phase = {}
                for phase in PHASES:
                    res = load_lrg_result(
                        patient, phase, band, fc_method="imcoh_abs"
                    )
                    if res is None or res.linkage_matrix is None:
                        raise FileNotFoundError(
                            f"{patient} {phase} {band} LRG missing"
                        )
                    Z_by_phase[phase] = np.asarray(res.linkage_matrix)

                n_lrg = int(Z_by_phase["rest_pre"].shape[0]) + 1

                # Sanity-check: probe_ids length must match the FC matrix
                # used at LRG time. If channel_labels.csv has a different
                # length, truncate the longer side conservatively to the
                # smaller — same convention as audit_47/50.
                A = load_fc_matrix(patient, "rest_pre", band, fc_method="imcoh_abs")
                n_fc = int(np.asarray(A).shape[0])
                n = min(n_lrg, n_fc, len(probe_ids))
                pid = probe_ids[:n]

                cohort_baseline = cohort_same_probe_pair_fraction(pid)

                anchor_lambda = find_anchor_triples(
                    Z_by_phase["rest_pre"],
                    Z_by_phase["task_test"],
                    Z_by_phase["rest_post"],
                    n,
                )
                trace_lambda = find_trace_pairs(
                    Z_by_phase["task_test"],
                    Z_by_phase["rest_post"],
                    Z_by_phase["rest_pre"],
                    n,
                )
                reset_lambda = find_reset_pairs(
                    Z_by_phase["rest_pre"],
                    Z_by_phase["rest_post"],
                    Z_by_phase["task_test"],
                    n,
                )
                rearrange_modules = find_rearrangement_modules(
                    Z_by_phase["rest_pre"],
                    Z_by_phase["task_test"],
                    Z_by_phase["rest_post"],
                    n,
                )

                modules_by_class = {
                    "anchor": list(anchor_lambda[0.0]),
                    "trace": list(trace_lambda[0.0]),
                    "reset": list(reset_lambda[0.0]),
                    "rearrange": list(rearrange_modules),
                }

                tie_broken = apply_tie_breaker(modules_by_class, n)

                # Per-module same-probe stats (using strict-intersection
                # leaves, NOT the tie-broken residual — the question is
                # whether anchors at their natural strict definition are
                # anatomy-anchored. The tie-breaker is reported in a
                # separate column as a sensitivity check.)
                for cls in ("anchor", "trace", "reset", "rearrange"):
                    raw_modules = modules_by_class[cls]
                    residuals = {idx: leaves for idx, leaves in tie_broken[cls]}
                    for idx, m in enumerate(raw_modules):
                        if cls == "anchor":
                            strict = strict_anchor_leaves(m)
                        elif cls == "trace":
                            strict = strict_trace_leaves(m)
                        elif cls == "reset":
                            strict = strict_reset_leaves(m)
                        else:
                            strict = strict_rearrange_leaves(m)
                        n_pairs, n_same, n_cross = module_same_probe_pair_stats(
                            strict, pid
                        )
                        same_frac = n_same / n_pairs if n_pairs else float("nan")
                        cross_frac = n_cross / n_pairs if n_pairs else float("nan")

                        residual = residuals.get(idx)
                        if residual is None:
                            res_n_pairs = res_n_same = 0
                            res_same_frac = float("nan")
                        else:
                            res_n_pairs, res_n_same, _ = (
                                module_same_probe_pair_stats(residual, pid)
                            )
                            res_same_frac = (
                                res_n_same / res_n_pairs
                                if res_n_pairs else float("nan")
                            )

                        rows.append({
                            "patient": patient,
                            "band": band,
                            "class": cls,
                            "module_id": idx,
                            "n_leaves": len(strict),
                            "n_pairs": n_pairs,
                            "n_same_probe_pairs": n_same,
                            "n_cross_probe_pairs": n_cross,
                            "same_probe_frac": same_frac,
                            "cross_probe_frac": cross_frac,
                            "cohort_baseline_same_probe_pair_frac": cohort_baseline,
                            "n_residual_leaves_after_tiebreak": (
                                len(residual) if residual is not None else 0
                            ),
                            "residual_n_pairs": res_n_pairs,
                            "residual_same_probe_frac": res_same_frac,
                        })
                print(
                    f"[audit_60] {patient} {band}: cohort_baseline="
                    f"{cohort_baseline:.3f} | anchor={len(modules_by_class['anchor'])} "
                    f"trace={len(modules_by_class['trace'])} "
                    f"reset={len(modules_by_class['reset'])} "
                    f"rearr={len(modules_by_class['rearrange'])}"
                )
            except Exception as e:
                print(f"[audit_60] WARN {patient} {band}: {e}")
                traceback.print_exc()

    df = pd.DataFrame(rows)
    csv_path = OUT / "per_module_same_probe_fraction.csv"
    df.to_csv(csv_path, index=False)
    print(f"[audit_60] per-module CSV: {csv_path}")

    # Cohort aggregate per (band, class): median same-probe-pair fraction.
    # Drop rows with NaN (zero-leaf or zero-pair modules).
    df_valid = df.dropna(subset=["same_probe_frac"]).copy()
    agg = (
        df_valid.groupby(["band", "class"])
        .agg(
            n_modules=("same_probe_frac", "size"),
            median_same_probe_frac=("same_probe_frac", "median"),
            mean_same_probe_frac=("same_probe_frac", "mean"),
            q25=("same_probe_frac", lambda x: x.quantile(0.25)),
            q75=("same_probe_frac", lambda x: x.quantile(0.75)),
            median_leaves=("n_leaves", "median"),
            median_cohort_baseline=("cohort_baseline_same_probe_pair_frac", "median"),
        )
        .reset_index()
    )
    agg["enrichment_ratio_median"] = (
        agg["median_same_probe_frac"] / agg["median_cohort_baseline"]
    )
    agg_path = OUT / "cohort_aggregate.csv"
    agg.to_csv(agg_path, index=False)
    print(f"[audit_60] aggregate CSV: {agg_path}")

    # Print β headline numbers for the report head.
    beta_rows = agg[agg["band"] == "beta"].sort_values("class")
    print("\n[audit_60] BETA cohort-median same-probe-pair fractions:")
    for _, r in beta_rows.iterrows():
        print(
            f"  {r['class']:>10s}: n={int(r['n_modules']):3d}  "
            f"median_same_probe={r['median_same_probe_frac']:.3f}  "
            f"baseline={r['median_cohort_baseline']:.3f}  "
            f"enrichment={r['enrichment_ratio_median']:.2f}x"
        )


if __name__ == "__main__":
    main()
