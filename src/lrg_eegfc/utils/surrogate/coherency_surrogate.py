"""On-manifold coherency surrogates for ImCoh-derived network statistics.

.. note::
   **ORPHANED (2026-06-12).** The analyses these helpers powered (CRC + the
   coordinated cross-phase null / SB-CRC) were archived as opaque, construction-
   dependent matrix-level nulls — see
   ``scripts/archive/2026-06_opaque-matrix-nulls/POSTMORTEM.md``. The functions are
   left in place because they are general (Haar sampling, complex-coherency
   recompute) and harmless, but nothing in the live pipeline imports them. The
   principled null for this project is matched-strength
   (:mod:`lrg_eegfc.utils.surrogate.matched_strength`). Do not rebuild
   coherency-rotation surrogates for the cross-phase trace.


Coupling-Randomized Coherency (CRC): a surrogate null that keeps every
realization a *genuine* (Hermitian positive-semidefinite, unit-diagonal)
coherency — hence a realizable ``|ImCoh|`` adjacency with non-negative edges,
valid graph Laplacian and well-defined LRG propagator — while randomizing
*which* node-sets carry the lagged coupling. It is the on-manifold complement
of the off-manifold matched-strength 4-cycle rewiring in
:mod:`lrg_eegfc.utils.surrogate.matched_strength`: matched-strength preserves
node strength but leaves the realizable set; CRC preserves realizability and
the per-frequency coupling content but lets strength float. Run the two as
brackets — a cross-phase trace that clears both is robust to the strongest
version of each objection.

Scope report (math, 5-point preamble, properties):
``.agents/guides/task-persistence-investigation/2026-06-08_on-manifold-coherency-surrogate.md``

Construction (``real_congruence``, the validated default)
---------------------------------------------------------
For each in-band frequency ``f`` the complex coherency
``C(f) = D^{-1/2} S(f) D^{-1/2}`` (``S`` = Welch cross-spectral density,
``D = diag(S)``) is conjugated by a single Haar-random **real orthogonal**
matrix ``O`` (shared across the band, fresh per surrogate, independent per
phase):

    C̃(f) = O · C(f) · Oᵀ                       # Hermitian PSD; Λ(f) and ‖Im C(f)‖_F preserved
    C̃(f) ← diag(C̃(f))^{-1/2} C̃(f) diag(C̃(f))^{-1/2}   # unit-diagonal congruence → valid coherency
    W_surr = mean_f |Im C̃(f)|                   # ≥ 0, [0,1], zero-diagonal

A real-orthogonal congruence preserves the Frobenius norm of the imaginary
part — the actual ImCoh energy — so the surrogate ``|ImCoh|`` magnitude matches
the data (calibration ratio ≈ 1). The earlier ``complex_eigbasis`` variant
(``C̃ = Q diag(Λ) Qᴴ`` with a complex Haar ``Q``) is kept only as a diagnostic:
it *generates* imaginary mass out of real coherence and inflates ``|ImCoh|``
several-fold (Pat_05 β: ratio 4.2 vs 1.1), so it is **not** a valid null.

References
----------
Mezzadri (2007) "How to generate random matrices from the classical compact
groups." Notices AMS 54(5):592-604 — Haar sampling of U(N)/O(N).
Nolte et al. (2004) — imaginary coherency. Ewald et al. (2012) — |ImCoh|.
"""
from __future__ import annotations

import numpy as np

from lrg_eegfc.utils.fc.coherence._common import welch_csd


__all__ = [
    "haar_unitary",
    "haar_orthogonal",
    "complex_coherency_band",
    "complex_coherency_bands",
    "coupling_randomized_coherency",
    "shared_backbone_deviation",
    "deviation_rotated_coherency",
]


def haar_unitary(N: int, rng: np.random.Generator) -> np.ndarray:
    """Haar-distributed complex unitary in U(N) (Mezzadri 2007).

    Used only by the diagnostic ``complex_eigbasis`` mode; the validated
    ``real_congruence`` null uses :func:`haar_orthogonal`.
    """
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    ph = np.diagonal(R).copy()
    ph /= np.abs(ph)
    return Q * ph                                   # Q @ diag(ph)


def haar_orthogonal(N: int, rng: np.random.Generator) -> np.ndarray:
    """Haar-distributed real orthogonal matrix in O(N) (Mezzadri 2007, real case)."""
    Z = rng.standard_normal((N, N))
    Q, R = np.linalg.qr(Z)
    return Q * np.sign(np.diagonal(R))


def complex_coherency_band(
    X: np.ndarray, fs: float, band: tuple[float, float], nperseg: int
) -> np.ndarray:
    """Full complex coherency stack ``C(f)`` for ``f`` in ``band``.

    ``C(f) = D^{-1/2} S(f) D^{-1/2}``, ``D = diag(S)``; Hermitian PSD, unit
    diagonal. Normalization is identical to
    :func:`lrg_eegfc.utils.fc.coherence.imcoh.compute_imcoh`, so the observed
    adjacency ``mean_f |Im C(f)|`` reproduces the cached ``imcoh_abs`` matrix
    bit-for-bit. The signed/derived caches store only ``Im C`` — the full
    complex ``C`` (real coherence + imaginary part) is recomputed here because
    the orthogonal congruence needs the whole Hermitian structure.

    Parameters
    ----------
    X : (N, L) ndarray
        Channel x time. (Caller ensures channels-first; split halves on axis 1.)
    fs : float
        Sampling rate (Hz).
    band : (fmin, fmax)
        Inclusive frequency band in Hz.
    nperseg : int
        Welch segment length (use ``nperseg_for_fs(fs)``; halve it for split-half
        baselines, matching ``compute_imcoh_abs_halves``).

    Returns
    -------
    (F_b, N, N) complex ndarray
        Coherency per in-band frequency bin.
    """
    freqs, CSD = welch_csd(X, fs, nperseg=nperseg)          # (N, N, F) complex Hermitian PSD
    PSD = np.real(np.diagonal(CSD, axis1=0, axis2=1).T)     # (N, F)
    denom = np.sqrt(PSD[:, None, :] * PSD[None, :, :])
    C = np.divide(CSD, denom, out=np.zeros_like(CSD), where=denom > 0)
    mask = (freqs >= band[0]) & (freqs <= band[1])
    return np.moveaxis(C[:, :, mask], 2, 0)                 # (F_b, N, N)


def complex_coherency_bands(
    X: np.ndarray, fs: float, bands: dict, nperseg: int
) -> dict:
    """Multi-band sibling of :func:`complex_coherency_band` — one Welch pass.

    Shares the single ``welch_csd`` FFT across all requested bands (a perf
    wrapper); the per-band normalization and frequency mask are identical to
    :func:`complex_coherency_band`, so ``out[b]`` is bit-for-bit what that function
    returns for band ``b``.

    Parameters
    ----------
    X : (N, L) ndarray
        Channels-first timeseries (split halves on axis 1 before calling).
    fs : float
        Sampling rate (Hz).
    bands : dict[str, (fmin, fmax)]
        Inclusive frequency bands in Hz.
    nperseg : int
        Welch segment length (``nperseg_for_fs(fs)``; for a coordinated cross-phase
        null all phases must share one ``nperseg`` so the grids align).

    Returns
    -------
    dict[str, (F_b, N, N) complex]
    """
    freqs, CSD = welch_csd(X, fs, nperseg=nperseg)
    PSD = np.real(np.diagonal(CSD, axis1=0, axis2=1).T)
    denom = np.sqrt(PSD[:, None, :] * PSD[None, :, :])
    C = np.divide(CSD, denom, out=np.zeros_like(CSD), where=denom > 0)
    out = {}
    for b, (flo, fhi) in bands.items():
        mask = (freqs >= flo) & (freqs <= fhi)
        out[b] = np.moveaxis(C[:, :, mask], 2, 0)
    return out


def coupling_randomized_coherency(
    C_band: np.ndarray,
    rng: np.random.Generator,
    mode: str = "real_congruence",
    eps: float = 1e-12,
) -> tuple[np.ndarray, float]:
    """One CRC surrogate adjacency ``W_surr = mean_f |Im C̃(f)|``.

    Parameters
    ----------
    C_band : (F_b, N, N) complex ndarray
        Per-frequency complex coherency (from :func:`complex_coherency_band`).
    rng : numpy Generator
        Randomness for the Haar rotation (one rotation shared across the band).
    mode : {"real_congruence", "complex_eigbasis"}
        ``real_congruence`` (default, validated): ``C̃ = O C(f) Oᵀ`` with real
        orthogonal ``O``; preserves ``Λ(f)`` and the imaginary energy
        ``‖Im C(f)‖_F`` (magnitude-faithful). ``complex_eigbasis`` (diagnostic
        only): ``C̃ = Q diag(Λ) Qᴴ`` with complex unitary ``Q``; preserves
        ``Λ`` but inflates ``|ImCoh|`` — NOT a valid null.
    eps : float
        Floor on the renormalization diagonal to avoid division blow-up.

    Returns
    -------
    (W_surr, renorm_residual) : ((N, N) ndarray, float)
        ``W_surr`` ≥ 0, symmetric, zero-diagonal — a valid Laplacian input.
        ``renorm_residual`` = mean over frequencies of ``max_i |diag(C̃_raw)_i − 1|``,
        the perturbation the unit-diagonal congruence applies to the preserved
        spectrum (report as calibration; large values mean the spectrum is held
        only approximately even though the magnitude is faithful).
    """
    Fb, N, _ = C_band.shape
    accum = np.zeros((N, N), dtype=float)
    resid = 0.0
    if mode == "real_congruence":
        O = haar_orthogonal(N, rng)
        # Batched over frequency: C_raw[f] = O @ C_band[f] @ Oᵀ
        C_raw = np.einsum("ij,fjk,lk->fil", O, C_band, O, optimize=True)
        d = np.real(np.diagonal(C_raw, axis1=1, axis2=2))           # (Fb, N)
        resid = float(np.mean(np.max(np.abs(d - 1.0), axis=1)))
        inv = 1.0 / np.sqrt(np.clip(d, eps, None))                  # (Fb, N)
        C_til = C_raw * inv[:, :, None] * inv[:, None, :]
        accum = np.abs(np.imag(C_til)).sum(axis=0)
        return accum / Fb, resid
    if mode == "complex_eigbasis":
        Q = haar_unitary(N, rng)
        for f in range(Fb):
            lam = np.clip(np.linalg.eigvalsh(C_band[f]), 0.0, None)
            C_raw = (Q * lam) @ Q.conj().T
            d = np.real(np.diag(C_raw))
            resid += float(np.max(np.abs(d - 1.0)))
            inv = 1.0 / np.sqrt(np.clip(d, eps, None))
            accum += np.abs(np.imag(C_raw * inv[:, None] * inv[None, :]))
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return accum / Fb, resid / Fb


# ---------------------------------------------------------------------------
# Coordinated cross-phase null (SB-CRC): shared backbone + deviation rotation
# ---------------------------------------------------------------------------
# Scope: .agents/guides/task-persistence-investigation/2026-06-11_coordinated-cross-phase-null.md
#
# CRC (above) randomizes each phase INDEPENDENTLY, so under the §5.3
# split-baseline protocol its cross-phase ρ_split floor collapses to ≈ 0 — it
# tests "structure vs four independent draws", not "task-specific vs a stable
# coupling backbone shared across all phases". SB-CRC closes that gap: it holds a
# cross-phase shared backbone ``B(f)`` fixed across all surrogate phases and
# rotates only each phase's deviation ``Δ_φ(f) = C_φ(f) − B(f)`` with an
# INDEPENDENT real-orthogonal ``O_φ``. Because the cophenetic map is nonlinear, a
# shared backbone does not cancel in the Δ_task / Δ_rest differences, so the floor
# lifts off zero to exactly the shared-backbone contribution; the observed trace
# is thereby decomposed into a shared-backbone part and a task-specific part.


def shared_backbone_deviation(
    C_by_phase: dict[str, np.ndarray],
    weights: dict[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Cross-phase shared backbone ``B(f)`` and per-phase deviations ``Δ_φ(f)``.

    ``B(f) = Σ_φ w_φ C_φ(f)`` (a convex combination of valid coherencies, hence
    Hermitian PSD with unit diagonal — itself a valid coherency); ``Δ_φ = C_φ − B``
    is Hermitian with zero diagonal. The weighting is a caller concern: the audit
    script knows the phase semantics and passes an explicit ``weights`` dict
    (e.g. equal-condition, counting the two ``rest_pre`` halves once between them);
    ``weights=None`` is the plain mean over the supplied phases. The library stays
    phase-name agnostic.

    Parameters
    ----------
    C_by_phase : dict[str, (F_b, N, N) complex]
        Per-phase complex coherency stacks on a COMMON frequency grid (all phases
        computed at the same ``nperseg`` — see scope §5.1). Mismatched ``F_b``
        across phases is an error (the per-frequency average is undefined).
    weights : dict[str, float] or None
        Convex weights ``Σ_φ w_φ = 1``. ``None`` → uniform ``1/len``.

    Returns
    -------
    (B_band, Delta_by_phase) : ((F_b, N, N) complex, dict[str, (F_b, N, N) complex])
    """
    phases = list(C_by_phase)
    if not phases:
        raise ValueError("C_by_phase is empty")
    shapes = {p: C_by_phase[p].shape for p in phases}
    if len({s for s in shapes.values()}) != 1:
        raise ValueError(f"phases must share a common (F_b, N, N) grid; got {shapes}")
    if weights is None:
        weights = {p: 1.0 / len(phases) for p in phases}
    wsum = float(sum(weights[p] for p in phases))
    if not np.isclose(wsum, 1.0):
        raise ValueError(f"weights must sum to 1 over the supplied phases; got {wsum}")
    B = sum(weights[p] * C_by_phase[p] for p in phases)
    Delta = {p: C_by_phase[p] - B for p in phases}
    return B, Delta


def deviation_rotated_coherency(
    B_band: np.ndarray,
    Delta_band: np.ndarray,
    rng: np.random.Generator,
    eps: float = 1e-12,
    psd_project: bool = True,
    return_psd_resid: bool = False,
) -> tuple[np.ndarray, float, float]:
    """One phase's SB-CRC surrogate adjacency ``W_surr = mean_f |Im C̃(f)|``.

    ``C̃(f) = renorm_unitdiag( PSD[ B(f) + O Δ(f) Oᵀ ] )`` with one Haar
    real-orthogonal ``O`` shared across the band (fresh per surrogate, INDEPENDENT
    per phase). The real-orthogonal congruence preserves ``‖Im Δ(f)‖_F`` exactly
    (the deviation's ImCoh energy) and leaves ``B``'s imaginary part untouched.

    **PSD projection (``psd_project=True``, default).** The recombination
    ``B + O Δ Oᵀ`` is NOT guaranteed positive-semidefinite — the cross term
    ``−O B Oᵀ`` can drive eigenvalues (and hence diagonal entries) negative. Without
    a projection the unit-diagonal renormalization then divides by a near-zero
    ``√diag`` and the surrogate magnitude **blows up** (cohort calibration: some
    patients hit a magnitude ratio of 1e5–1e8 — Pat_06 β). Projecting each
    ``B + O Δ Oᵀ`` onto the nearest PSD matrix (clip negative eigenvalues to 0,
    batched ``eigh`` over frequency) restores a genuine Hermitian-PSD coherency, so
    every surrogate is on-manifold, the renorm is stable, and the magnitude stays
    faithful. The single-patient pilot (Pat_05, well-conditioned) missed this; the
    cohort exposed it — so projection is the default, not an option. Set
    ``psd_project=False`` only to reproduce the raw (off-manifold) construction.

    Parameters
    ----------
    B_band, Delta_band : (F_b, N, N) complex
        Shared backbone and this phase's deviation, from
        :func:`shared_backbone_deviation`.
    rng : numpy Generator
        Randomness for the one real-orthogonal rotation.
    eps : float
        Floor on the renormalization diagonal.
    psd_project : bool
        Project ``B + O Δ Oᵀ`` onto the nearest PSD matrix before the unit-diagonal
        renorm (default True — required for cohort-stable magnitudes, see above).
    return_psd_resid : bool
        If True, also return the most-negative eigenvalue of the raw (pre-projection)
        recombination — the on-manifold-fidelity diagnostic (how far the projection
        had to move it). With ``psd_project=True`` this is free (the ``eigh`` is
        already computed); with both False the third return is ``nan``. ``W_surr ≥ 0``
        holds regardless, so the Laplacian / propagator are always valid.

    Returns
    -------
    (W_surr, renorm_resid, psd_resid) : ((N, N) float ≥ 0, float, float)
    """
    Fb, N, _ = B_band.shape
    O = haar_orthogonal(N, rng)
    Drot = np.einsum("ij,fjk,lk->fil", O, Delta_band, O, optimize=True)
    C_raw = B_band + Drot
    psd_resid = float("nan")
    if psd_project:
        w, V = np.linalg.eigh(C_raw)                            # batched (Fb,N), (Fb,N,N)
        psd_resid = float(w.min())
        w = np.clip(w, 0.0, None)
        C_raw = (V * w[:, None, :]) @ np.conj(np.swapaxes(V, 1, 2))
    elif return_psd_resid:
        psd_resid = float(min(np.linalg.eigvalsh(C_raw[f]).min() for f in range(Fb)))
    d = np.real(np.diagonal(C_raw, axis1=1, axis2=2))           # (Fb, N)
    resid = float(np.mean(np.max(np.abs(d - 1.0), axis=1)))
    inv = 1.0 / np.sqrt(np.clip(d, eps, None))                  # (Fb, N)
    C_til = C_raw * inv[:, :, None] * inv[:, None, :]
    W = np.abs(np.imag(C_til)).sum(axis=0) / Fb
    return W, resid, psd_resid
