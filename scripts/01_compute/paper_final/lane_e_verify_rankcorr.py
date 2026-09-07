#!/usr/bin/env python3
"""Verify the rank-correlation tensor reproduces the direct functionals exactly.

The lane's whole re-derivation strategy rests on one claim: the ``(8, 8)``
Spearman matrix of the reorganisation vectors is a *sufficient statistic* for
the cross-phase family, so storing it lets any estimator -- including ones not
yet written -- be evaluated post hoc without recomputing an eigendecomposition.
That claim is checkable, and this checks it on random graphs:

1. ``cross_phase_functionals_from_corr(R)`` equals ``cross_phase_functionals(D)``
   to floating point, for the four canonical functionals.
2. ``cross_phase_functionals_from_corr(R, swap=True)`` equals the direct
   computation with the encode/probe roles exchanged -- i.e. the swap really is
   an exact signed relabelling and needs no second pass over the data.
3. The algebraic swap identities hold: ``T_test(swap) = T_learn``,
   ``T_learn(swap) = T_test``, ``T_infspec(swap) = -T_infspec``.
4. :func:`partial_from_corr_matrix` at first order equals the closed form.
5. The drift-controlled functionals are *exactly* the uncontrolled ones when
   ``d`` carries no information (constructed case), which is the sanity check
   that conditioning is not silently shifting the estimator on its own.
6. The vectorised stack derivation (recursive first-order partial identity,
   used to evaluate millions of surrogate matrices in one pass) agrees
   elementwise with the per-matrix version that goes through a matrix inverse.
"""
from __future__ import annotations

import numpy as np

from lrg_eegfc.utils.fc.heat_multiscale import (
    CROSS_PHASE_FUNCTIONALS,
    CROSS_PHASE_ROLES,
    CROSS_PHASE_VECTOR_NAMES,
    _partial_from_corr,
    cophenetic_at_scale,
    cross_phase_functionals,
    cross_phase_functionals_from_corr,
    cross_phase_functionals_from_corr_stack,
    cross_phase_rank_corr,
    laplacian_eig,
)

PHASES = ("A", "B", "task_learn", "task_test", "rest_post")
SWAP_ROLES = dict(CROSS_PHASE_ROLES, encode=CROSS_PHASE_ROLES["probe"],
                  probe=CROSS_PHASE_ROLES["encode"])


def random_arc(rng, n=40):
    """Five correlated random weighted graphs -- a stand-in for one arc."""
    base = rng.random((n, n))
    base = 0.5 * (base + base.T)
    out = {}
    for ph in PHASES:
        W = np.clip(base + 0.35 * rng.standard_normal((n, n)), 0.01, 1.0)
        W = 0.5 * (W + W.T)
        np.fill_diagonal(W, 0.0)
        out[ph] = W
    return out


def main():
    rng = np.random.default_rng(11)
    worst = dict(direct=0.0, swap=0.0, ident=0.0, partial=0.0)
    for trial in range(6):
        Ws = random_arc(rng)
        eig = {ph: laplacian_eig(Ws[ph]) for ph in PHASES}
        for s in (1.0, 5.6, 45.1, 180.0):
            D = {ph: cophenetic_at_scale(*eig[ph], s) for ph in PHASES}
            R = cross_phase_rank_corr(D)
            got = cross_phase_functionals_from_corr(R)
            ref = cross_phase_functionals(D)
            got_sw = cross_phase_functionals_from_corr(R, swap=True)
            ref_sw = cross_phase_functionals(D, SWAP_ROLES)
            for k in CROSS_PHASE_FUNCTIONALS:
                worst["direct"] = max(worst["direct"], abs(got[k] - ref[k]))
                worst["swap"] = max(worst["swap"], abs(got_sw[k] - ref_sw[k]))
            worst["ident"] = max(
                worst["ident"],
                abs(ref_sw["T_test"] - ref["T_learn"]),
                abs(ref_sw["T_learn"] - ref["T_test"]),
                abs(ref_sw["T_infspec"] + ref["T_infspec"]),
            )
            ix = {k: i for i, k in enumerate(CROSS_PHASE_VECTOR_NAMES)}
            from lrg_eegfc.utils.fc.heat_multiscale import partial_from_corr_matrix
            a, b, c = ix["f"], ix["p"], ix["e"]
            worst["partial"] = max(worst["partial"], abs(
                partial_from_corr_matrix(R, a, b, [c])
                - _partial_from_corr(R[a, b], R[a, c], R[b, c])))

    # drift conditioning on an uninformative d: build R with d decorrelated
    R = cross_phase_rank_corr(D)
    R0 = R.copy()
    k = CROSS_PHASE_VECTOR_NAMES.index("d")
    R0[k, :] = 0.0
    R0[:, k] = 0.0
    R0[k, k] = 1.0
    g0 = cross_phase_functionals_from_corr(R0)
    d_null = max(abs(g0[f] - g0[f + "_d"]) for f in ("T_test", "T_learn", "T_infspec"))

    # vectorised stack derivation vs the per-matrix one, on random correlations
    stack = np.empty((4, 5, 8, 8))
    for a in range(4):
        for b in range(5):
            stack[a, b] = np.corrcoef(rng.standard_normal((8, 60)))
    s_worst = 0.0
    for sw in (False, True):
        got = cross_phase_functionals_from_corr_stack(stack, swap=sw)
        for a in range(4):
            for b in range(5):
                ref = cross_phase_functionals_from_corr(stack[a, b], swap=sw)
                for k in ref:
                    s_worst = max(s_worst, abs(float(got[k][a, b]) - float(ref[k])))

    print("max |tensor - direct|              :", worst["direct"])
    print("max |tensor(swap) - direct(swap)|  :", worst["swap"])
    print("max |swap algebraic identities|    :", worst["ident"])
    print("max |partial matrix - closed form| :", worst["partial"])
    print("max |drift-ctrl - plain| at d=0    :", d_null)
    print("max |stack - per-matrix|           :", s_worst)
    tol = 1e-9
    bad = ([k for k, v in worst.items() if v > tol]
           + (["d_null"] if d_null > tol else [])
           + (["stack"] if s_worst > tol else []))
    print("\nVERDICT:", "PASS" if not bad else f"FAIL {bad}")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
