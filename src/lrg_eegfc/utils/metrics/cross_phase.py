"""Cross-phase decomposition of per-pair distance trajectories.

Given a per-pair (ultrametric / cophenetic) distance in three phases
(rest_pre, task, rest_post), decompose the cross-phase fluctuation of each pair
into two orthonormal contrasts — *persistence* (``phi1``) and *excursion*
(``phi2``) — plus the DC level, and derive the anchor / trace / reset /
reorganize taxonomy channel scores. General graph/statistics primitive: no
manuscript-local scope.

The split is **exact**: ``var(d_pre, d_task, d_post) = (phi1**2 + phi2**2) / 3``,
so every pair's cross-phase fluctuation energy decomposes with no residual into
persistence energy ``phi1**2`` and excursion energy ``phi2**2``.

Consistency with the validated split-baseline trace
``rho_split = Spearman(d_task - d_pre_A, d_post - d_pre_B)``: the persistence
contrast is the rest axis up to baseline noise —
``d_post - d_pre_A = sqrt(2) * phi1`` — and a perfect reset
(``d_post == d_pre``) has ``phi1 = 0``, i.e. it lies in the **null-space** of
``rho_split``. The excursion (reset) channel therefore cannot perturb the trace
statistic; computing it leaves ``rho_split`` unchanged.
"""
from __future__ import annotations

import numpy as np

_SQRT2 = np.sqrt(2.0)
_SQRT6 = np.sqrt(6.0)


def phase_contrasts(d_pre, d_task, d_post):
    """Orthonormal persistence / excursion contrasts of a 3-phase trajectory.

    Parameters
    ----------
    d_pre, d_task, d_post : ndarray
        Per-pair distance in each phase (condensed, elementwise aligned).

    Returns
    -------
    phi1 : ndarray
        Persistence contrast ``(d_post - d_pre) / sqrt(2)`` — RPost vs RPre.
    phi2 : ndarray
        Excursion contrast ``(2*d_task - d_pre - d_post) / sqrt(6)`` —
        Task vs the rest-average.
    var : ndarray
        ``(phi1**2 + phi2**2) / 3`` = cross-phase variance of the triple.
    """
    phi1 = (d_post - d_pre) / _SQRT2
    phi2 = (2.0 * d_task - d_pre - d_post) / _SQRT6
    var = (phi1 * phi1 + phi2 * phi2) / 3.0
    return phi1, phi2, var


def cross_phase_channels(d_pre_a, d_pre_b, d_task, d_post, kappa=1.0):
    """Per-pair taxonomy channels from split-baseline cophenetic distances.

    ``d_pre_a`` is the rest_pre representative for the contrast triple and for
    the task displacement; ``d_pre_b`` is the independent baseline for the rest
    displacement (the validated split-baseline construction). ``sigma`` is the
    within-rest_pre noise scale ``sqrt(mean (d_pre_a - d_pre_b)**2)`` (RMS).

    RMS, not ``median |·|``: cophenetic distances are quantized to <= N-1 merge
    heights, so ~98% of baseline displacements are exact ties (zeros) and the
    median is 0, which would clamp sigma and label every pair a mover. The RMS
    is robustly non-zero and is the natural amplitude scale to compare against
    ``var`` (an energy).

    Returns a dict of per-pair arrays plus the scalar ``sigma``:
    ``phi1, phi2, var, sigma, anchor_weight, anchor_score, reset_score,
    mover_mask, dD_task, dD_rest, dD_base``.

    - ``anchor_weight = exp(-var / (2 sigma**2))`` in ``(0, 1]`` (1 = rigid).
    - ``mover_mask = var > kappa * sigma**2``.
    - ``reset_score = (phi2 / sigma)`` for movers that *returned*
      (``|dD_task| > sigma`` and ``|dD_rest| <= sigma``), else 0 — sign-aware.
    """
    d_pre_a = np.asarray(d_pre_a, dtype=float)
    d_pre_b = np.asarray(d_pre_b, dtype=float)
    d_task = np.asarray(d_task, dtype=float)
    d_post = np.asarray(d_post, dtype=float)

    phi1, phi2, var = phase_contrasts(d_pre_a, d_task, d_post)
    dD_task = d_task - d_pre_a
    dD_rest = d_post - d_pre_b
    dD_base = d_pre_a - d_pre_b
    sigma = max(float(np.sqrt(np.mean(dD_base * dD_base))), 1e-12)
    sig2 = sigma * sigma

    anchor_weight = np.exp(-var / (2.0 * sig2))
    anchor_score = -var / sig2
    mover_mask = var > kappa * sig2
    moved_task = np.abs(dD_task) > sigma
    returned = np.abs(dD_rest) <= sigma
    reset_score = (phi2 / sigma) * (moved_task & returned).astype(float)

    return dict(
        phi1=phi1, phi2=phi2, var=var, sigma=sigma,
        anchor_weight=anchor_weight, anchor_score=anchor_score,
        reset_score=reset_score, mover_mask=mover_mask,
        dD_task=dD_task, dD_rest=dD_rest, dD_base=dD_base,
    )


def cohort_share_row(channels):
    """Aggregate per-pair ``channels`` (``cross_phase_channels`` output) to a
    cell-level composition.

    Returns dict with:
    ``anchor_mass`` (mean anchor weight), ``mover_frac`` (fraction of movers),
    ``E_T``/``E_R`` (mover persistence/excursion energy), ``share_T``/``share_R``
    (mover energy split, sum to 1), and the mass×energy composition
    ``comp_anchor = 1 - mover_frac``, ``comp_trace = mover_frac * share_T``,
    ``comp_reset = mover_frac * share_R`` (these three sum to 1; a ``reorganize``
    slice is later carved from comp_trace+comp_reset by the coupled null).
    """
    mov = channels["mover_mask"]
    phi1 = channels["phi1"]
    phi2 = channels["phi2"]
    E_T = float(np.sum(phi1[mov] ** 2))
    E_R = float(np.sum(phi2[mov] ** 2))
    tot = E_T + E_R
    mover_frac = float(np.mean(mov))
    share_T = E_T / tot if tot > 0 else float("nan")
    share_R = E_R / tot if tot > 0 else float("nan")
    return dict(
        anchor_mass=float(np.mean(channels["anchor_weight"])),
        mover_frac=mover_frac,
        E_T=E_T, E_R=E_R,
        share_T=share_T, share_R=share_R,
        comp_anchor=1.0 - mover_frac,
        comp_trace=(mover_frac * share_T if tot > 0 else 0.0),
        comp_reset=(mover_frac * share_R if tot > 0 else 0.0),
    )
