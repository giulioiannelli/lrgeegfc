#!/usr/bin/env python3
"""Audit 61 — Direction E: epi Grassmann embedding (E.align + E.resect).

Scope: ``.agents/guides/task-persistence-investigation/2026-05-08_epi-grassmann-embedding.md``.

Two scalars per (patient, band, phase, k):

- ``eps_align`` — chordal Grassmann distance from the top-k eigenspace
  ``V_k`` of the full FC graph's Laplacian to the epi coordinate subspace
  ``S_E = span{e_i : i ∈ E_p} ⊆ R^N``. Equivalent reading: average
  fraction ``f_E_k = ‖V_k[E_p, :]‖_F² / k`` of the top-k mode mass that
  sits on epi rows.
- ``delta_resect`` — chordal Grassmann distance between the full-graph
  eigenspace restricted to non-epi rows (orthonormalised via QR) and the
  eigenspace freshly computed on the resected FC graph (epi rows/cols
  deleted, Laplacian recomputed). A complementary scalar
  ``lambda_shift_k = (λ_k^full − λ_k^resect) / λ_max^full`` is reported
  per cell.

Strength-stratified null: per (patient, phase, band) FC strength ``s_i``
is binned into ``Q = 5`` quintiles; ``R = 100`` random draws sample, in
each quintile, the same number of nodes as ``E_p`` has from the
non-epi pool. For each draw the same two scalars are computed with
``S_r`` substituted for ``E_p``.

**No acceptance gate.** The script emits the per-patient observation,
the per-patient null mean/std, the per-patient z-score, and a
descriptive cohort summary. Signal-vs-noise judgment happens post hoc
from the figures (audit_61b).

Outputs
-------
``data/audit/epi_grassmann/E_per_patient.csv``
    One row per (patient, band, phase, k) with both variants' raw
    scalars, null statistics, z-scores, and `rank_block` /
    `lambda_shift` diagnostics.

``data/audit/epi_grassmann/E_cohort.csv``
    One row per (band, phase, k, variant) with cohort median / Q1 / Q3 /
    min / max for the observed scalar AND the z-score, plus footnote
    columns Wilcoxon two-sided p and BH-FDR q within band.

CLI
---
``--sanity-only``   Run the SBM sanity gate and exit.
``--align-only``    Skip E.resect (faster, ~30 s cohort).
``--resect-only``   Skip E.align.
``--n-null R``      Override R = 100 null draws.
``--n-quintile Q``  Override Q = 5 strength quintiles.
``--seed S``        Override RNG seed (default 20260510).

Run inside the ``lapbrain`` conda env.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES
from lrg_eegfc.utils.io import load_channel_labels, load_epileptic_nodes
from lrg_eegfc.utils.metrics import (
    bh_fdr,
    chordal_distance,
    chordal_full_vs_resect,
    grassmann_to_coord_subspace,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

# ---------------------------------------------------------------------------
# Cohort + grid (n=9 epi-annotated; Pat_15 has no red-coloured contacts)
# ---------------------------------------------------------------------------

COHORT = (
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14",
)
PHASES = ("rest_pre", "task_test", "rest_post")
K_GRID = (2, 3, 5, 8, 13, 21)
DEFAULT_N_NULL = 100
DEFAULT_Q_QUINTILE = 5
DEFAULT_SEED = 20260510

OUT_DIR = ROOT / "data" / "audit" / "epi_grassmann"
FC_METHOD = "imcoh_abs"


# ---------------------------------------------------------------------------
# Per-patient masks
# ---------------------------------------------------------------------------


def epi_indices(patient: str) -> np.ndarray:
    """Return integer indices of epi contacts in FC-matrix order."""
    labels = load_channel_labels(patient)
    epi_set = set(load_epileptic_nodes(patient))
    return np.array(
        [i for i, lbl in enumerate(labels) if lbl in epi_set], dtype=int
    )


# ---------------------------------------------------------------------------
# Strength-stratified null
# ---------------------------------------------------------------------------


def quintile_bin(strength: np.ndarray, Q: int) -> np.ndarray:
    """Return per-node quintile id in [0, Q-1]."""
    edges = np.quantile(strength, np.linspace(0.0, 1.0, Q + 1))
    edges[-1] = np.inf       # ensure max is included
    edges[0] = -np.inf
    return np.digitize(strength, edges[1:-1])  # length Q-1 internal cuts


def stratified_draw(
    rng: np.random.Generator,
    quint: np.ndarray,
    epi_idx: np.ndarray,
    pool: np.ndarray,
    Q: int,
) -> np.ndarray:
    """Draw a size-matched node subset from ``pool`` matching the per-quintile
    composition of ``epi_idx``. If a quintile contains too few pool nodes,
    fall back to sampling without replacement from the whole pool for the
    short quintile."""
    n_per_q = np.zeros(Q, dtype=int)
    for q in range(Q):
        n_per_q[q] = int(np.sum(quint[epi_idx] == q))
    drawn = []
    short = 0
    for q in range(Q):
        avail = pool[quint[pool] == q]
        if avail.size >= n_per_q[q]:
            drawn.append(rng.choice(avail, size=n_per_q[q], replace=False))
        else:
            drawn.append(avail)
            short += int(n_per_q[q] - avail.size)
    out = np.concatenate(drawn) if drawn else np.array([], dtype=int)
    if short > 0:
        rest_pool = np.setdiff1d(pool, out)
        # rare: backfill from outside-quintile pool nodes
        rest = rng.choice(
            rest_pool, size=min(short, rest_pool.size), replace=False
        )
        out = np.concatenate([out, rest])
    return out


# ---------------------------------------------------------------------------
# Resected eigendecomposition
# ---------------------------------------------------------------------------


def laplacian_eigh(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Symmetric Laplacian eigendecomposition (ascending eigenvalues)."""
    A = np.asarray(A, dtype=float)
    A = 0.5 * (A + A.T)              # enforce symmetry
    np.fill_diagonal(A, 0.0)
    deg = A.sum(axis=1)
    L = np.diag(deg) - A
    L = 0.5 * (L + L.T)
    lam, V = np.linalg.eigh(L)
    return lam, V


def topk(V: np.ndarray, k: int) -> np.ndarray:
    return np.ascontiguousarray(V[:, 1 : k + 1])


# ---------------------------------------------------------------------------
# Per-cell computation (one (p, b, phi))
# ---------------------------------------------------------------------------


def compute_cell(
    *,
    patient: str,
    phase: str,
    band: str,
    epi_idx: np.ndarray,
    rng: np.random.Generator,
    n_null: int,
    Q: int,
    do_align: bool,
    do_resect: bool,
) -> list[dict]:
    """Compute per-k rows for one (patient, phase, band) cell."""
    A = load_fc_matrix(patient, phase, band, FC_METHOD)
    if A is None:
        raise FileNotFoundError(f"FC missing: {patient} {phase} {band}")
    res = load_lrg_result(patient, phase, band, fc_method=FC_METHOD)
    if res is None or res.eigenvectors is None or res.eigenvalues is None:
        raise FileNotFoundError(f"LRG eigvecs missing: {patient} {phase} {band}")

    N = A.shape[0]
    if res.eigenvectors.shape != (N, N):
        raise ValueError(
            f"eigenvector shape {res.eigenvectors.shape} != (N, N)=({N},{N})"
        )

    V_full = res.eigenvectors
    lam_full = res.eigenvalues
    lam_max_full = float(lam_full.max())

    pool = np.setdiff1d(np.arange(N), epi_idx)
    strength = A.sum(axis=1)
    quint = quintile_bin(strength, Q)

    # Resected Laplacian (E.resect only)
    V_R = None
    lam_R = None
    if do_resect:
        A_R = A[np.ix_(pool, pool)]
        lam_R, V_R = laplacian_eigh(A_R)

    rows: list[dict] = []
    for k in K_GRID:
        if k + 1 >= N:
            continue
        Vk_full = topk(V_full, k)

        # ---- E.align observed ---------------------------------------
        if do_align:
            eps_obs, f_obs = grassmann_to_coord_subspace(Vk_full, epi_idx)
        else:
            eps_obs, f_obs = (np.nan, np.nan)

        # ---- E.resect observed --------------------------------------
        if do_resect and V_R is not None:
            if k + 1 >= V_R.shape[1]:
                delta_obs, rank_obs = (np.nan, 0)
                lam_shift = np.nan
            else:
                Vk_R = topk(V_R, k)
                delta_obs, rank_obs = chordal_full_vs_resect(
                    Vk_full, Vk_R, pool
                )
                # eigenvalue shift (relative): use λ_{k+1} (top-k after dropping v_1)
                lam_shift = float(
                    (lam_full[k] - lam_R[k]) / lam_max_full
                ) if lam_max_full > 0 else np.nan
        else:
            delta_obs, rank_obs, lam_shift = (np.nan, 0, np.nan)

        # ---- nulls ---------------------------------------------------
        eps_null = []
        delta_null = []
        for _ in range(n_null):
            S_r = stratified_draw(rng, quint, epi_idx, pool, Q)
            if S_r.size == 0:
                continue
            if do_align:
                eps_S, _ = grassmann_to_coord_subspace(Vk_full, S_r)
                eps_null.append(eps_S)
            if do_resect:
                pool_S = np.setdiff1d(np.arange(N), S_r)
                if pool_S.size <= k + 1:
                    continue
                A_RS = A[np.ix_(pool_S, pool_S)]
                lam_RS, V_RS = laplacian_eigh(A_RS)
                if k + 1 >= V_RS.shape[1]:
                    continue
                Vk_RS = topk(V_RS, k)
                delta_S, _ = chordal_full_vs_resect(Vk_full, Vk_RS, pool_S)
                delta_null.append(delta_S)

        eps_null_arr = np.asarray(eps_null) if eps_null else np.array([])
        delta_null_arr = np.asarray(delta_null) if delta_null else np.array([])
        eps_mu = float(eps_null_arr.mean()) if eps_null_arr.size else np.nan
        eps_sd = float(eps_null_arr.std(ddof=1)) if eps_null_arr.size > 1 else np.nan
        delta_mu = float(delta_null_arr.mean()) if delta_null_arr.size else np.nan
        delta_sd = float(delta_null_arr.std(ddof=1)) if delta_null_arr.size > 1 else np.nan
        z_align = (
            (eps_obs - eps_mu) / eps_sd
            if (do_align and not np.isnan(eps_obs) and eps_sd not in (0.0, np.nan))
            else np.nan
        )
        z_resect = (
            (delta_obs - delta_mu) / delta_sd
            if (do_resect and not np.isnan(delta_obs) and delta_sd not in (0.0, np.nan))
            else np.nan
        )

        rows.append(dict(
            patient=patient,
            phase=phase,
            band=band,
            k=k,
            n_nodes=N,
            n_epi=int(epi_idx.size),
            n_pool=int(pool.size),
            eps_align_obs=eps_obs,
            f_align_obs=f_obs,
            eps_align_null_mean=eps_mu,
            eps_align_null_std=eps_sd,
            z_align=z_align,
            delta_resect_obs=delta_obs,
            rank_block=rank_obs,
            lambda_shift=lam_shift,
            delta_resect_null_mean=delta_mu,
            delta_resect_null_std=delta_sd,
            z_resect=z_resect,
            n_null_align=int(eps_null_arr.size),
            n_null_resect=int(delta_null_arr.size),
        ))
    return rows


# ---------------------------------------------------------------------------
# Cohort aggregation
# ---------------------------------------------------------------------------


def cohort_summary(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (band, phase, k, variant) with descriptive cohort stats."""
    summary = []
    for variant, obs_col, z_col in (
        ("align", "eps_align_obs", "z_align"),
        ("resect", "delta_resect_obs", "z_resect"),
    ):
        for band in BRAIN_BANDS_NAMES:
            # collect per-(phase, k) p-values for BH-FDR within band
            cells = []
            for phase in PHASES:
                for k in K_GRID:
                    sub = df[
                        (df.band == band) & (df.phase == phase) & (df.k == k)
                    ]
                    obs = sub[obs_col].dropna().values
                    z = sub[z_col].dropna().values
                    if z.size >= 2 and np.any(z != 0):
                        try:
                            _, p_two = wilcoxon(z, alternative="two-sided")
                        except Exception:
                            p_two = np.nan
                    else:
                        p_two = np.nan
                    cells.append(dict(
                        band=band, phase=phase, k=k, variant=variant,
                        n_patients=int(z.size),
                        obs_median=float(np.median(obs)) if obs.size else np.nan,
                        obs_q1=float(np.quantile(obs, 0.25)) if obs.size else np.nan,
                        obs_q3=float(np.quantile(obs, 0.75)) if obs.size else np.nan,
                        obs_min=float(obs.min()) if obs.size else np.nan,
                        obs_max=float(obs.max()) if obs.size else np.nan,
                        z_median=float(np.median(z)) if z.size else np.nan,
                        z_q1=float(np.quantile(z, 0.25)) if z.size else np.nan,
                        z_q3=float(np.quantile(z, 0.75)) if z.size else np.nan,
                        n_pos_z_neg=int(np.sum(z < 0)),  # epi-aligned direction for align
                        n_pos_z_pos=int(np.sum(z > 0)),  # resection-disrupts direction
                        wilcoxon_p_two_sided=float(p_two),
                    ))
            # BH-FDR within band over all (phase × k) cells of this variant
            ps = np.array([c["wilcoxon_p_two_sided"] for c in cells], dtype=float)
            qs = np.full_like(ps, np.nan)
            mask = ~np.isnan(ps)
            if mask.any():
                qs[mask] = bh_fdr(ps[mask])
            for c, q in zip(cells, qs):
                c["bh_q_within_band"] = float(q)
                summary.append(c)
    return pd.DataFrame(summary)


# ---------------------------------------------------------------------------
# SBM sanity gate
# ---------------------------------------------------------------------------


def sbm_block_adjacency(
    rng: np.random.Generator,
    n_blocks: int,
    nodes_per_block: int,
    p_intra: float,
    p_inter: float,
    sigma: float,
) -> np.ndarray:
    N = n_blocks * nodes_per_block
    A = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            same = (i // nodes_per_block) == (j // nodes_per_block)
            mu = p_intra if same else p_inter
            w = max(0.0, mu + rng.normal(0, sigma))
            A[i, j] = A[j, i] = w
    return A


def sanity_gate(seed: int = 12345) -> bool:
    """Two-part SBM sanity construction.

    PART 1 (E.align distinguishability) — *asymmetric* graph: one dense
    cluster of size ``n_E`` embedded in a sparse background of size
    ``N - n_E``. The Fiedler vector v_2 separates the dense cluster
    from the rest with anomalously high mass on the cluster. E.align
    against ``E_p = cluster`` should be much smaller than E.align
    against a random size-matched subset of background nodes.

    A *symmetric* k-block SBM does NOT pass this part: by symmetry every
    block carries equal top-(n_blocks−1)-mode mass (= 1/n_blocks of
    total), so E.align cannot distinguish block1 from any size-matched
    random subset. This is why we use the asymmetric construction.

    PART 2 (E.resect distinguishability) — symmetric 4-block SBM:
    removing a whole block creates a structurally different graph (the
    block-indicator subspace shrinks from 3D to 2D); removing 25
    random nodes leaves a perturbed 4-block graph whose top-3 modes
    still resemble the full top-3.
    """
    print("\n=== SBM sanity gate ===")
    rng = np.random.default_rng(seed)

    # ---- PART 1: asymmetric "dense cluster + background" ----
    n_E, n_B = 30, 70
    p_E, p_B, p_inter, sigma1 = 0.9, 0.3, 0.05, 0.05
    N1 = n_E + n_B
    A1 = np.zeros((N1, N1))
    for i in range(N1):
        for j in range(i + 1, N1):
            same_E = i < n_E and j < n_E
            same_B = i >= n_E and j >= n_E
            mu = p_E if same_E else (p_B if same_B else p_inter)
            w = max(0.0, mu + rng.normal(0, sigma1))
            A1[i, j] = A1[j, i] = w

    _, V1 = laplacian_eigh(A1)
    cluster = np.arange(n_E)
    background = np.arange(n_E, N1)
    rand_match = rng.choice(background, size=n_E, replace=False)

    Vk1 = topk(V1, 1)                                # k = 1: just the Fiedler
    eps_cluster, f_cluster = grassmann_to_coord_subspace(Vk1, cluster)
    eps_rand, f_rand = grassmann_to_coord_subspace(Vk1, rand_match)
    reflex = chordal_distance(Vk1, Vk1)
    print(f"  PART-1 (asymmetric): eps_cluster={eps_cluster:.4f} (f={f_cluster:.3f}); "
          f"eps_rand={eps_rand:.4f} (f={f_rand:.3f}); reflex={reflex:.2e}")
    assert reflex < 1e-6, "reflexive check failed"

    # Cluster's mass on the Fiedler should be much higher than a random
    # background sample's mass; eps_cluster much smaller than eps_rand.
    pass1 = eps_cluster < eps_rand - 0.2
    print(f"  PASS-1 eps(cluster) << eps(rand) at k=1? {pass1}")

    # ---- PART 2: symmetric 4-block SBM, E.resect ----
    n_blocks, npb, p_in, p_out, sigma2 = 4, 25, 0.8, 0.1, 0.05
    A2 = sbm_block_adjacency(rng, n_blocks, npb, p_in, p_out, sigma2)
    _, V2 = laplacian_eigh(A2)
    block1 = np.arange(npb)
    rand25 = rng.choice(np.arange(npb, n_blocks * npb), size=npb, replace=False)
    pool_b = np.setdiff1d(np.arange(A2.shape[0]), block1)
    pool_r = np.setdiff1d(np.arange(A2.shape[0]), rand25)

    Vk2 = topk(V2, 3)
    A_resect_b = A2[np.ix_(pool_b, pool_b)]
    A_resect_r = A2[np.ix_(pool_r, pool_r)]
    _, V_Rb = laplacian_eigh(A_resect_b)
    _, V_Rr = laplacian_eigh(A_resect_r)

    delta_block1, rank_b1 = chordal_full_vs_resect(
        Vk2, topk(V_Rb, 3), pool_b
    )
    delta_rand, rank_r = chordal_full_vs_resect(
        Vk2, topk(V_Rr, 3), pool_r
    )
    print(f"  PART-2 (symmetric 4-block): delta_block1={delta_block1:.4f} "
          f"(rank={rank_b1}); delta_rand={delta_rand:.4f} (rank={rank_r})")
    pass2 = delta_block1 > delta_rand
    print(f"  PASS-2 delta(block1) > delta(rand25)? {pass2}")

    ok = pass1 and pass2
    print(f"=== SBM sanity gate: {'PASS' if ok else 'FAIL'} ===\n")
    return ok


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit 61 — epi Grassmann embedding")
    parser.add_argument("--sanity-only", action="store_true",
                        help="Run SBM sanity gate and exit")
    parser.add_argument("--align-only", action="store_true")
    parser.add_argument("--resect-only", action="store_true")
    parser.add_argument("--n-null", type=int, default=DEFAULT_N_NULL)
    parser.add_argument("--n-quintile", type=int, default=DEFAULT_Q_QUINTILE)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.sanity_only:
        ok = sanity_gate(args.seed)
        raise SystemExit(0 if ok else 1)

    if args.align_only and args.resect_only:
        raise SystemExit("--align-only and --resect-only are mutually exclusive")
    do_align = not args.resect_only
    do_resect = not args.align_only

    # Always run the sanity gate before the cohort
    if not sanity_gate(args.seed):
        raise SystemExit("SBM sanity gate failed; aborting cohort run")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows: list[dict] = []
    for patient in COHORT:
        try:
            E_p = epi_indices(patient)
        except Exception as e:
            print(f"[audit_61] WARN {patient}: cannot load epi indices ({e})")
            continue
        if E_p.size < 2:
            print(f"[audit_61] WARN {patient}: |E_p|={E_p.size} too small; skip")
            continue
        for band in BRAIN_BANDS_NAMES:
            for phase in PHASES:
                try:
                    out = compute_cell(
                        patient=patient,
                        phase=phase,
                        band=band,
                        epi_idx=E_p,
                        rng=rng,
                        n_null=args.n_null,
                        Q=args.n_quintile,
                        do_align=do_align,
                        do_resect=do_resect,
                    )
                    rows.extend(out)
                    if args.verbose:
                        print(f"  [{patient} {band} {phase}] {len(out)} rows")
                except Exception as e:
                    print(f"[audit_61] WARN {patient} {band} {phase}: {e}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "E_per_patient.csv", index=False)
    print(f"[audit_61] wrote {OUT_DIR / 'E_per_patient.csv'} ({len(df)} rows)")

    summary = cohort_summary(df)
    summary.to_csv(OUT_DIR / "E_cohort.csv", index=False)
    print(f"[audit_61] wrote {OUT_DIR / 'E_cohort.csv'} ({len(summary)} rows)")

    if args.verbose:
        print(summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
