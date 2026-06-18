"""audit_96 — Coupling-Randomized Coherency (CRC) surrogate, construction pilot.

Scope: .agents/guides/task-persistence-investigation/2026-06-08_on-manifold-coherency-surrogate.md
(5-point critical preamble + math live there; this is the construction-validation
pilot, NOT the cohort run).

Goal of this pilot
------------------
1. Confirm the CRC construction runs end-to-end and produces VALID surrogates:
   non-negative edges (valid Laplacian), genuine Hermitian-PSD unit-diagonal
   coherency (on-manifold).
2. Confirm the observed |ImCoh| recomputed from the timeseries matches the cached
   `load_fc_matrix(..., 'imcoh_abs')` (bit-faithful observed path).
3. Show the CRC null WIDTH for the shared-baseline cophenetic trace statistic on
   one patient x beta, and where the observed value sits.

This pilot uses the SHARED-baseline rho_split proxy
(Delta_task = coph_task - coph_pre ; Delta_rest = coph_post - coph_pre).
The cohort runner will switch to the split-baseline halves + add the
matched-strength bracket (audit_63 cache), per the scope.

The three CRC core functions below (`haar_unitary`,
`complex_coherency_band`, `coupling_randomized_coherency`) PROMOTE to
`src/lrg_eegfc/utils/surrogate/coherency_surrogate.py` on the cohort caller
(2nd use), per coding-rules.md. Kept local here for the pilot only.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.stats import spearmanr, ks_2samp

from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import BRAIN_BANDS, nperseg_for_fs, DEFAULT_SAMPLE_RATE
from lrg_eegfc.config.const import FS_OVERRIDES
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.fc.coherence._common import welch_csd
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.utils.surrogate import cophenetic_condensed_from_eigs


# --------------------------------------------------------------------------
# CRC core (promotes to lrg_eegfc.utils.surrogate.coherency_surrogate on 2nd caller)
# --------------------------------------------------------------------------
def haar_unitary(N: int, rng: np.random.Generator) -> np.ndarray:
    """Haar-distributed complex unitary in U(N) (Mezzadri 2007)."""
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    ph = np.diagonal(R).copy()
    ph /= np.abs(ph)
    return Q * ph  # Q @ diag(ph)


def complex_coherency_band(
    X: np.ndarray, fs: float, band: tuple[float, float], nperseg: int
) -> np.ndarray:
    """Full complex coherency stack C(f) for f in band -> (F_b, N, N) complex.

    C(f) = D^{-1/2} S(f) D^{-1/2}, D = diag(S). Hermitian PSD, unit diagonal.
    Matches `compute_imcoh`'s normalization; the observed |ImCoh| adjacency is
    mean_f |Im C(f)|.
    """
    freqs, CSD = welch_csd(X, fs, nperseg=nperseg)          # CSD: (N, N, F)
    PSD = np.real(np.diagonal(CSD, axis1=0, axis2=1).T)      # (N, F)
    denom = np.sqrt(PSD[:, None, :] * PSD[None, :, :])       # (N, N, F)
    C = np.divide(CSD, denom, out=np.zeros_like(CSD), where=denom > 0)
    fmin, fmax = band
    mask = (freqs >= fmin) & (freqs <= fmax)
    return np.moveaxis(C[:, :, mask], 2, 0)                  # (F_b, N, N)


def haar_orthogonal(N: int, rng: np.random.Generator) -> np.ndarray:
    """Haar-distributed real orthogonal matrix in O(N) (Mezzadri 2007, real case)."""
    Z = rng.standard_normal((N, N))
    Q, R = np.linalg.qr(Z)
    return Q * np.sign(np.diagonal(R))


def crc_surrogate(
    C_band: np.ndarray,
    Lam: np.ndarray,
    rng: np.random.Generator,
    mode: str = "real_congruence",
    eps: float = 1e-12,
) -> tuple[np.ndarray, float]:
    """One CRC surrogate adjacency W_surr = mean_f |Im C̃(f)|.

    mode='complex_eigbasis' (v1): C̃ = Q diag(Λ) Qᴴ, complex Haar Q on the
        eigenvalue matrix. Preserves Λ; GENERATES imaginary mass (inflates |ImCoh|).
    mode='real_congruence'  (v2): C̃ = O C(f) Oᵀ, real Haar orthogonal O on the
        FULL complex coherency. Preserves Λ AND ‖Im C(f)‖_F exactly (no inflation);
        rotates which nodes carry the lagged coupling.
    Both renormalize to unit diagonal (PSD-preserving congruence). Returns
    (W_surr >= 0, mean_f max|diag(C̃_raw) - 1|).
    """
    Fb, N = Lam.shape
    accum = np.zeros((N, N), dtype=float)
    resid = 0.0
    if mode == "complex_eigbasis":
        Q = haar_unitary(N, rng)
        for f in range(Fb):
            lam = np.clip(Lam[f], 0.0, None)
            C_raw = (Q * lam) @ Q.conj().T
            d = np.real(np.diag(C_raw))
            resid += float(np.max(np.abs(d - 1.0)))
            inv = 1.0 / np.sqrt(np.clip(d, eps, None))
            accum += np.abs(np.imag(C_raw * inv[:, None] * inv[None, :]))
    elif mode == "real_congruence":
        O = haar_orthogonal(N, rng)
        for f in range(Fb):
            C_raw = O @ C_band[f] @ O.T                      # real-orthogonal congruence
            d = np.real(np.diag(C_raw))
            resid += float(np.max(np.abs(d - 1.0)))
            inv = 1.0 / np.sqrt(np.clip(d, eps, None))
            accum += np.abs(np.imag(C_raw * inv[:, None] * inv[None, :]))
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return accum / Fb, resid / Fb


def _coph_from_adjacency(W: np.ndarray) -> np.ndarray:
    """W (>=0 symmetric) -> LRG cophenetic condensed vector at tau = 1/lam_max."""
    L = np.diag(W.sum(axis=1)) - W
    lam, V = np.linalg.eigh(L)
    return cophenetic_condensed_from_eigs(lam, V)


# --------------------------------------------------------------------------
# Pilot run
# --------------------------------------------------------------------------
def main() -> None:
    PATIENT = "Pat_05"
    BAND = "beta"
    PHASES = ("rest_pre", "task_test", "rest_post")
    R = 200
    SEED = 20260608

    fs = FS_OVERRIDES.get(PATIENT, DEFAULT_SAMPLE_RATE)
    nperseg = nperseg_for_fs(fs)
    band = BRAIN_BANDS[BAND]
    print(f"[audit_96] {PATIENT} x {BAND}  fs={fs} nperseg={nperseg} band={band} R={R}")

    # ---- observed: complex coherency per phase, |ImCoh| adjacency, cophenetic ----
    C_band = {}
    A_obs = {}
    coph_obs = {}
    for ph in PHASES:
        X = load_timeseries(PATIENT, ph, SEEG_DATAPATH)
        Cb = complex_coherency_band(X, fs, band, nperseg)   # (F_b, N, N)
        C_band[ph] = Cb
        A = np.mean(np.abs(np.imag(Cb)), axis=0)            # observed |ImCoh|
        np.fill_diagonal(A, 0.0)
        A_obs[ph] = A
        coph_obs[ph] = _coph_from_adjacency(A)
    N = A_obs[PHASES[0]].shape[0]
    Fb = C_band[PHASES[0]].shape[0]
    print(f"[audit_96] N={N} contacts, F_b={Fb} in-band freq bins")

    # ---- sanity: observed |ImCoh| vs cached load_fc_matrix ----
    A_cache = load_fc_matrix(PATIENT, "rest_pre", BAND, "imcoh_abs")
    if A_cache.shape == A_obs["rest_pre"].shape:
        rel = np.max(np.abs(A_cache - A_obs["rest_pre"])) / (np.max(A_cache) + 1e-12)
        print(f"[audit_96] observed-vs-cache max rel-diff (rest_pre): {rel:.3e}")
    else:
        print(f"[audit_96] WARN shape mismatch obs {A_obs['rest_pre'].shape} "
              f"vs cache {A_cache.shape} — skipping bit-faithful check")

    # ---- observed statistic: shared-baseline rho_split (cophenetic) ----
    d_task_obs = coph_obs["task_test"] - coph_obs["rest_pre"]
    d_rest_obs = coph_obs["rest_post"] - coph_obs["rest_pre"]
    T_obs = spearmanr(d_task_obs, d_rest_obs).statistic
    print(f"[audit_96] OBSERVED rho_split^coph (shared-baseline) = {T_obs:+.4f}")

    # ---- precompute per-frequency coherency eigenvalues per phase ----
    Lam = {ph: np.array([np.linalg.eigvalsh(C_band[ph][f]) for f in range(Fb)])
           for ph in PHASES}  # (F_b, N) each

    # ---- CRC null: run BOTH modes for the v1-vs-v2 calibration comparison ----
    mag_obs = float(np.mean(A_obs['rest_post'][np.triu_indices(N, 1)]))
    triu = np.triu_indices(N, 1)
    for mode in ("complex_eigbasis", "real_congruence"):
        rng = np.random.default_rng(SEED)
        T_null = np.empty(R)
        strength_drift_ks, mag_surr, resid_all = [], [], []
        valid_nonneg = True
        t0 = time.time()
        for r in range(R):
            coph_s = {}
            for ph in PHASES:
                W, resid = crc_surrogate(C_band[ph], Lam[ph], rng, mode=mode)
                resid_all.append(resid)
                if W.min() < -1e-12:
                    valid_nonneg = False
                coph_s[ph] = _coph_from_adjacency(W)
                if ph == "rest_post":
                    mag_surr.append(float(np.mean(W[triu])))
                    strength_drift_ks.append(
                        ks_2samp(W.sum(1), A_obs[ph].sum(1)).statistic)
            dt = coph_s["task_test"] - coph_s["rest_pre"]
            dr = coph_s["rest_post"] - coph_s["rest_pre"]
            T_null[r] = spearmanr(dt, dr).statistic
        elapsed = time.time() - t0

        finite = T_null[np.isfinite(T_null)]
        p_upper = float(np.mean(finite >= T_obs))
        z = (T_obs - finite.mean()) / (finite.std() + 1e-12)
        tag = "v1" if mode == "complex_eigbasis" else "v2"
        print(f"\n[audit_96] ===== CRC {tag} ({mode}) — {R} surr, {elapsed:.1f}s =====")
        print(f"    null mean +/- std : {finite.mean():+.4f} +/- {finite.std():.4f}")
        print(f"    null  5/50/95 pct : {np.percentile(finite,5):+.4f} / "
              f"{np.percentile(finite,50):+.4f} / {np.percentile(finite,95):+.4f}")
        print(f"    OBSERVED          : {T_obs:+.4f}")
        print(f"    upper-tail p      : {p_upper:.4f}   (positive = trace)")
        print(f"    observed z vs null: {z:+.2f}")
        print(f"    calib non-neg (valid Laplacian) : {valid_nonneg}")
        print(f"    calib mean|ImCoh| obs={mag_obs:.4f} surr={np.mean(mag_surr):.4f} "
              f"ratio={np.mean(mag_surr)/mag_obs:.2f}   <-- target ~1")
        print(f"    calib strength KS surr-vs-obs   : {np.mean(strength_drift_ks):.3f} "
              f"(>0 expected; strength floats = bracket vs matched-strength)")
        print(f"    calib Lambda renorm residual    : {np.mean(resid_all):.3e}")


if __name__ == "__main__":
    main()
