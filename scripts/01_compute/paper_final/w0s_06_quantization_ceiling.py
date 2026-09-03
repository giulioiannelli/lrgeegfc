#!/usr/bin/env python3
r"""Why the cophenetic readout expresses so little, and why that makes it scale-flat.

Two numbers in this project have been reported as separate unexplained facts:

  * the cophenetic geometry expresses only **1-6%** of the raw per-pair task
    reorganisation (W0-C: squared Spearman 0.013-0.062, every band, every scale);
  * a Spearman over ~6 900 contact pairs sees only about **20 distinct
    cophenetic values**, at every scale.

They are the same fact, and this stage shows it with an exact ceiling rather than
an analogy.

The mechanism, derived
----------------------
A UPGMA tree on ``N`` leaves has exactly ``K = N-1`` merges. The cophenetic
distance of a pair is the height of the merge that joins it, so

    C = sum_{g=1..K} h_g * 1_{G_g}

where ``G_1..G_K`` **partition** the ``M = N(N-1)/2`` pairs (each pair is joined
at exactly one merge) and the indicators are therefore orthogonal. The cophenetic
representation is thus a *coordinate subspace* of dimension ``K`` inside
``R^M``: it can carry only the between-group component of any per-pair signal,
and it is blind to everything that varies within a group. Two exact consequences.

**Ceiling 1 -- expressible variance.** For any per-pair signal ``u``, the largest
fraction of its variance any cophenetic vector can reproduce is the one-way
between-group ``R^2`` of ``u`` on the merge partition. For a signal carrying no
special relation to the tree this is ``(K-1)/(M-1)``; at ``N = 118`` that is
``117/6903 = 1.7%``. **The observed 1.3-6.2% is that number.** The cophenetic
geometry is not a filter that keeps the robust part of the raw signal, and it is
not orthogonal to it -- it is a ``K``-cell coarse-graining of a ``M``-cell object,
and 1.7% is what a ``K``-cell coarse-graining of an arbitrary signal retains.

**Ceiling 2 -- attainable rank correlation.** With tie-group counts ``n_g``, the
largest Spearman a tied variable can have with any untied one is exactly

    rho_max = sqrt( 1 - sum_g (n_g^3 - n_g) / (M^3 - M) )

(the standard tie correction; it is the ratio of the standard deviation of the
average-rank vector to that of the untied rank vector, and it is attained by any
target that orders the groups correctly).

A prediction, made before the candidates were scored
----------------------------------------------------
If quantisation is the mechanism, a scale-local readout helps **only if it
escapes the quantisation**, not if it merely reweights the same quantised values.
That sorts this lane's own candidates in advance:

  * ``Tloc`` / ``Qloc`` restrict the correlation to pairs inside one tree-level
    stratum -- and a stratum is *defined* by merge level, so it contains FEWER
    distinct cophenetic values than the whole. The band-pass makes the
    quantisation **worse**. It should fail.
  * ``Ccon`` / ``Qcon`` reweight the same quantised values without adding any.
    Neutral at best.
  * ``Thei`` correlates the ``K = N-1`` merge heights themselves, which are
    generically all distinct. It is the **only** candidate that escapes the
    quantisation, and the only one with a mechanistic reason to work.

So this stage reports, for every construction, how many distinct values its input
actually takes and how tied it is. A candidate that shows scale variation while
still quantised to ~20 levels should be disbelieved, not celebrated.

Everything here is observed-only and exact -- ceilings and counts, not tests.
Nothing is compared to zero, and nothing needs a surrogate, because a ceiling
derived from the tree's own combinatorics is not an inference about the data.
"""
from __future__ import annotations

import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import rankdata

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PATIENTS_4PHASE
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, linkage_at_scale
from lrg_eegfc.utils.metrics.tree import pair_merge_level
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.substrate import CANONICAL, canonical_graph

ROOT = setup_script_env()

from w0s_01_scale_locality_grid import (                       # noqa: E402
    BACKBONE, FRACS, NQ, PHASES, SGRID, _cr, _sp, _shannon_neff, _strata,
)

WORKERS = int(os.environ.get("W0S_WORKERS", "8"))
OUT = Path(os.environ.get(
    "W0S_OUT_ANALYSIS", ROOT / "data" / "paper_final" / "lane_s_scale"))


def _tie_ceiling(x: np.ndarray) -> float:
    """Largest Spearman a variable with these ties can reach against an untied one.

    ``sqrt(1 - sum_g (n_g^3 - n_g) / (M^3 - M))`` -- exact, attained when the
    target orders the tie groups correctly. It is a property of the value
    multiplicities alone, so it is a ceiling on *any* use of this vector in a
    rank correlation, whatever the other side is.
    """
    M = x.size
    if M < 3:
        return np.nan
    _, cnt = np.unique(x, return_counts=True)
    return float(np.sqrt(max(0.0, 1.0 - np.sum(cnt ** 3 - cnt)
                             / (M ** 3 - M))))


def _between_group_r2(u: np.ndarray, lab: np.ndarray) -> float:
    """One-way between-group ``R^2`` of ``u`` on a partition -- the expressible share.

    The exact fraction of ``u``'s variance that any function constant on the
    partition's cells can reproduce. For the cophenetic representation the cells
    are the merge groups, so this is the ceiling on what the hierarchy can say
    about a per-pair signal, before any question of whether it says it well.
    """
    u = np.asarray(u, float)
    u = u - u.mean()
    tot = float(u @ u)
    if tot <= 0:
        return np.nan
    idx, inv = np.unique(lab, return_inverse=True)
    cnt = np.bincount(inv).astype(float)
    ssum = np.bincount(inv, weights=u)
    return float(np.sum(ssum ** 2 / cnt) / tot)


def _reorg(d: dict) -> np.ndarray:
    """Symmetric per-pair task reorganisation ``1/2[(t-A)+(t-B)]``."""
    return 0.5 * ((d["task_test"] - d["A"]) + (d["task_test"] - d["B"]))


def per_cell(job):
    pat, band = job
    t0 = time.time()
    try:
        Ws = {ph: canonical_graph(pat, ph, band, dense=True) for ph in PHASES}
    except Exception as exc:                                     # noqa: BLE001
        return [dict(patient=pat, band=band, error=str(exc))]
    n = Ws["A"].shape[0]
    iu = np.triu_indices(n, 1)
    M = iu[0].size
    K = n - 1
    Araw = {ph: Ws[ph][iu] for ph in PHASES}
    u_raw = _reorg(Araw)
    r_raw = _cr(u_raw)
    rng = np.random.default_rng(0)

    rows = []
    for fi, f in enumerate(FRACS):
        eig = {ph: laplacian_eig(select_backbone(Ws[ph], BACKBONE, frac=float(f)))
               for ph in PHASES}
        for j, s in enumerate(SGRID):
            try:
                Z = {ph: linkage_at_scale(*eig[ph], s) for ph in PHASES}
                C = {ph: cophenet(Z[ph]) for ph in PHASES}
            except Exception:                                    # noqa: BLE001
                continue
            if any((not np.all(np.isfinite(v))) or np.std(v) == 0
                   for v in C.values()):
                continue
            levA = pair_merge_level(Z["A"])
            levT = pair_merge_level(Z["task_test"])
            u_coph = _reorg(C)

            # --- what the representation can express at all ------------------ #
            r2_raw_on_A = _between_group_r2(r_raw, levA)
            refine = levA.astype(np.int64) * (K + 1) + levT
            d_refine = int(np.unique(refine).size)
            r2_raw_on_refine = _between_group_r2(r_raw, refine)
            r2_rand = float(np.mean([
                _between_group_r2(rng.normal(size=M), levA) for _ in range(3)]))

            # --- how quantised the readout's input actually is --------------- #
            row = dict(
                patient=pat, band=band, frac=float(f), s=float(s), N=n, M=M,
                K_merges=K,
                n_distinct_coph=int(np.unique(C["A"]).size),
                n_eff_levels=_shannon_neff(np.bincount(levA, minlength=K + 1)[1:]),
                tie_ceiling_coph=_tie_ceiling(C["A"]),
                tie_ceiling_reorg=_tie_ceiling(u_coph),
                n_distinct_reorg=int(np.unique(u_coph).size),
                # ceilings on expressible variance
                r2_ceiling_uniform=(K - 1) / (M - 1),
                r2_raw_on_merge_partition=r2_raw_on_A,
                r2_raw_on_refined_partition=r2_raw_on_refine,
                d_refine=d_refine,
                r2_ceiling_refine_uniform=(d_refine - 1) / (M - 1),
                r2_random_signal_on_merge_partition=r2_rand,
                # the quantity W0-C measured, recomputed fresh
                r2_observed_coph_vs_raw=float(_sp(u_raw, u_coph) ** 2),
                # merge heights: the one representation that is NOT quantised
                n_distinct_heights=int(np.unique(Z["A"][:, 2]).size),
                tie_ceiling_heights=_tie_ceiling(Z["A"][:, 2]),
            )
            # --- per-stratum quantisation of this lane's own candidates ------ #
            _, q = _strata(levA, n)
            nd, tc = [], []
            for qq in range(1, NQ + 1):
                m = q == qq
                if m.sum() >= 3:
                    nd.append(np.unique(C["A"][m]).size)
                    tc.append(_tie_ceiling(C["A"][m]))
            row["n_distinct_coph_within_quintile_med"] = (
                float(np.median(nd)) if nd else np.nan)
            row["tie_ceiling_within_quintile_med"] = (
                float(np.median(tc)) if tc else np.nan)
            rows.append(row)
    rows[0]["elapsed_s"] = time.time() - t0
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, b) for b in BRAIN_BANDS_NAMES for p in PATIENTS_4PHASE]
    print(f"[w0s-quant] {len(jobs)} cells | {CANONICAL.label()} "
          f"| observed only, exact ceilings | {WORKERS} workers", flush=True)
    t0 = time.time()
    rows = []
    with Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows.extend(res)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:6.0f}s "
                  f"ETA {el/i*(len(jobs)-i):6.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "quantization_cells.csv", index=False)

    cols = ["n_distinct_coph", "n_eff_levels", "tie_ceiling_coph",
            "n_distinct_reorg", "tie_ceiling_reorg", "r2_ceiling_uniform",
            "r2_raw_on_merge_partition", "d_refine",
            "r2_raw_on_refined_partition",
            "r2_random_signal_on_merge_partition",
            "r2_observed_coph_vs_raw", "n_distinct_heights",
            "tie_ceiling_heights", "n_distinct_coph_within_quintile_med",
            "tie_ceiling_within_quintile_med"]
    summ = df.groupby("band")[cols].median().reset_index()
    summ.to_csv(OUT / "quantization_summary.csv", index=False)
    by_s = df.groupby("s")[cols].median().reset_index()
    by_s.to_csv(OUT / "quantization_by_scale.csv", index=False)

    print("\n-- the ceiling, per band (cohort median over patients, "
          "fractions and scales) --")
    print(summ[["band", "n_distinct_coph", "n_eff_levels", "tie_ceiling_coph",
                "r2_ceiling_uniform", "r2_raw_on_merge_partition",
                "r2_observed_coph_vs_raw", "n_distinct_heights",
                "n_distinct_coph_within_quintile_med"]].to_string(index=False))
    print(f"\n[w0s-quant] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
