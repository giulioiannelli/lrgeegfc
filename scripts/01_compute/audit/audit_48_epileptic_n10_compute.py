#!/usr/bin/env python3
"""Audit 48 — Epileptic-node revisit on the n=10 cohort.

Computes three measure batches into CSVs; figures live in
``audit_48b_epileptic_n10_figures.py``.

Cohort (n=9 with epi annotations; Pat_15 dropped — no red-coloured contacts).
Bands: BRAIN_BANDS_NAMES (six). Phases: rest_pre / task_test / rest_post.

M1 — Raw |ImCoh| epi enrichment
    R_x = mean(|ImCoh|_ij | i,j ∈ epi, cls=x) / mean(|ImCoh|_ij | non-non, cls=x)
    where cls ∈ {all, cp = different probe, sp = same probe}.

M2 — LRG ρ̂(τ) heat-kernel communicability enrichment
    K(τ) = V exp(-τ Λ) V^T from cached eigvals/eigvecs.
    ρ̂_ij(τ) = K_ij / sqrt(K_ii · K_jj).
    R_ρ_x = mean(ρ̂_ij | epi-epi, cls=x) / mean(ρ̂_ij | non-non, cls=x).
    Two τ regimes:
        τ_max = 1/λ_max (default LRG scale, mirrors Section 5).
        τ-sweep over a log grid from τ_max to argmax(C(τ)) ≈ τ*.
    Also: cophenetic ultrametric counterpart from r.ultrametric_matrix.

M3 — Hierarchical-position diagnostics
    (a) Per-pair MRCA merge height M[k] (= height of MRCA in the dendrogram)
        from kc_vectors(Z), pooled over the cohort, stratified by pair class
        {epi-epi, non-non, cross}.
    (b) KC distance restricted to the epi-epi pair subset across phase pairs
        (rest_pre↔task, task↔rest_post, rest_pre↔rest_post) at λ ∈ {0, 0.5, 1}.
        T_KC_epi = d_KC(rest_pre, task; λ) − d_KC(task, rest_post; λ).
        Positive = trace at the epi sub-network level (parallels Result 2).
        (Sign convention locked 2026-05-26: T_d > 0 = TRACE.)

All math via library helpers; no private forks.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io import load_epileptic_nodes
from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result

COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14",  # Pat_15 has no annotation
]
BANDS = list(BRAIN_BANDS_NAMES)
PHASES = ("rest_pre", "task_test", "rest_post")
TAU_SWEEP_NPTS = 10  # log-spaced points from τ_max to τ*
KC_LAMBDAS = (0.0, 0.5, 1.0)

OUT_DIR = ROOT / "data" / "audit" / "epileptic_n10_revisit"


# ---------------------------------------------------------------------------
# Channel + epi mask helpers (mirrors scripts/08_epileptic/_validation_one_patient.py)
# ---------------------------------------------------------------------------

def _load_channel_labels(patient: str) -> list[str]:
    """Return cleaned channel labels for a patient (matches FC-matrix ordering)."""
    p = SEEG_DATAPATH / patient / "channel_labels.csv"
    with open(p) as fh:
        first = fh.readline().strip().strip('"').strip("'")
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(p, header=None, skiprows=skip)
    return [str(x).strip('"').split(",")[0].strip().replace(" ", "")
            for x in df.iloc[:, 0]]


def _probe_of(label: str) -> str:
    m = re.match(r"([A-Za-z]+'?)", label)
    return m.group(1) if m else label


@dataclass
class PatientMasks:
    patient: str
    channels: list[str]
    epi_mask: np.ndarray         # (N,) bool
    probes: np.ndarray           # (N,) object (probe name)


def _build_masks(patient: str) -> PatientMasks:
    ch = _load_channel_labels(patient)
    epi_set = set(load_epileptic_nodes(patient))
    epi_mask = np.array([c in epi_set for c in ch])
    probes = np.array([_probe_of(c) for c in ch], dtype=object)
    return PatientMasks(patient, ch, epi_mask, probes)


# ---------------------------------------------------------------------------
# Pair stratification
# ---------------------------------------------------------------------------

def _pair_class_masks(epi_mask: np.ndarray, probes: np.ndarray
                      ) -> dict[str, np.ndarray]:
    """Return upper-triangle (i<j) boolean masks for the pair-classes we care about.

    Returns a dict with keys::
        ee_all, ee_cp, ee_sp,  nn_all, nn_cp, nn_sp,
        cross_all, cross_cp, cross_sp.
    Each value is a (npairs,) bool array aligned with the scipy
    upper-triangle row-major flattening (matches `kc_vectors` ordering).
    """
    n = len(epi_mask)
    iu, ju = np.triu_indices(n, k=1)
    same_probe = probes[iu] == probes[ju]
    e_i, e_j = epi_mask[iu], epi_mask[ju]

    ee = e_i & e_j
    nn = (~e_i) & (~e_j)
    cross = (e_i & ~e_j) | (~e_i & e_j)

    return {
        "ee_all": ee, "ee_cp": ee & ~same_probe, "ee_sp": ee & same_probe,
        "nn_all": nn, "nn_cp": nn & ~same_probe, "nn_sp": nn & same_probe,
        "cross_all": cross, "cross_cp": cross & ~same_probe, "cross_sp": cross & same_probe,
    }


def _ratio(num_vals: np.ndarray, den_vals: np.ndarray) -> float:
    """Mean(num) / Mean(den), guarded against empty / zero den."""
    if num_vals.size == 0 or den_vals.size == 0:
        return float("nan")
    den = float(den_vals.mean())
    if den <= 1e-30:
        return float("nan")
    return float(num_vals.mean()) / den


# ---------------------------------------------------------------------------
# Heat-kernel ρ̂(τ) from cached eigendecomposition
# ---------------------------------------------------------------------------

def heat_kernel_rho(eigvals: np.ndarray, eigvecs: np.ndarray, tau: float
                    ) -> np.ndarray:
    """Communicability ρ̂_ij(τ) = K_ij(τ) / sqrt(K_ii(τ) · K_jj(τ)).

    K(τ) = V · diag(exp(-τ λ)) · V^T (real, symmetric, PSD).
    """
    decay = np.exp(-tau * eigvals)
    K = (eigvecs * decay[None, :]) @ eigvecs.T
    diag = np.sqrt(np.maximum(np.diag(K), 1e-300))
    rho = K / np.outer(diag, diag)
    return rho


def _tau_grid(eigvals: np.ndarray, entropy_tau: np.ndarray,
              entropy_C: np.ndarray, npts: int) -> np.ndarray:
    """Log-spaced τ from 1/λ_max up to τ* = argmax C(τ)."""
    lam_max = float(eigvals[-1])
    tau_min = 1.0 / lam_max
    if entropy_C is None or entropy_tau is None or len(entropy_C) == 0:
        tau_star = tau_min * 100.0
    else:
        tau_star = float(entropy_tau[int(np.argmax(entropy_C))])
        tau_star = max(tau_star, tau_min * 1.5)
    return np.geomspace(tau_min, max(tau_star, tau_min * 1.5), npts)


# ---------------------------------------------------------------------------
# M1 + M2 driver
# ---------------------------------------------------------------------------

def compute_enrichment_for_patient(masks: PatientMasks
                                   ) -> tuple[list[dict], list[dict], list[dict]]:
    """Return (M1 rows, M2 rows at τ_max + cophenetic, M2 τ-sweep rows)."""
    pat = masks.patient
    pair_masks = _pair_class_masks(masks.epi_mask, masks.probes)

    m1_rows: list[dict] = []
    m2_rows: list[dict] = []
    sweep_rows: list[dict] = []

    iu, ju = np.triu_indices(len(masks.epi_mask), k=1)

    for band in BANDS:
        for phase in PHASES:
            # ---- M1: raw |ImCoh| ratios -----------------------------------
            A = load_fc_matrix(pat, phase, band, "imcoh_abs")
            if A is None:
                continue
            A_flat = A[iu, ju]
            for cls in ("all", "cp", "sp"):
                ee = A_flat[pair_masks[f"ee_{cls}"]]
                nn = A_flat[pair_masks[f"nn_{cls}"]]
                m1_rows.append(dict(
                    patient=pat, band=band, phase=phase, cls=cls,
                    n_ee=int(ee.size), n_nn=int(nn.size),
                    ee_mean=float(ee.mean()) if ee.size else float("nan"),
                    nn_mean=float(nn.mean()) if nn.size else float("nan"),
                    R_imcoh=_ratio(ee, nn),
                ))

            # ---- M2: ρ̂(τ_max) + cophenetic, plus τ-sweep ----------------
            r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs")
            if r is None or r.eigenvalues is None or r.eigenvectors is None:
                continue
            ev = np.asarray(r.eigenvalues)
            V = np.asarray(r.eigenvectors)
            tau_max = 1.0 / float(ev[-1])
            rho = heat_kernel_rho(ev, V, tau_max)
            rho_flat = rho[iu, ju]

            # cophenetic ultrametric distance — squareform(1D condensed)
            from scipy.spatial.distance import squareform
            um_full = squareform(np.asarray(r.ultrametric_matrix))
            um_flat = um_full[iu, ju]

            for cls in ("all", "cp", "sp"):
                ee_r = rho_flat[pair_masks[f"ee_{cls}"]]
                nn_r = rho_flat[pair_masks[f"nn_{cls}"]]
                ee_u = um_flat[pair_masks[f"ee_{cls}"]]
                nn_u = um_flat[pair_masks[f"nn_{cls}"]]
                m2_rows.append(dict(
                    patient=pat, band=band, phase=phase, cls=cls,
                    n_ee=int(ee_r.size), n_nn=int(nn_r.size),
                    rho_ee_mean=float(ee_r.mean()) if ee_r.size else float("nan"),
                    rho_nn_mean=float(nn_r.mean()) if nn_r.size else float("nan"),
                    R_rho=_ratio(ee_r, nn_r),
                    um_ee_mean=float(ee_u.mean()) if ee_u.size else float("nan"),
                    um_nn_mean=float(nn_u.mean()) if nn_u.size else float("nan"),
                    R_um=_ratio(ee_u, nn_u),
                    tau_max=tau_max,
                ))

            # τ-sweep (cross-probe only, fastest informative track).
            # Each patient/phase/band has its own (λ_max, τ*) so we record the
            # absolute τ AND the normalized τ × λ_max (= 1 at τ_min). The figure
            # aggregates over τ-index (relative position) to compare same scale.
            tau_grid = _tau_grid(ev, np.asarray(r.entropy_tau),
                                 np.asarray(r.entropy_1_minus_S), TAU_SWEEP_NPTS)
            lam_max = float(ev[-1])
            for tau_idx, tau in enumerate(tau_grid):
                rho_t = heat_kernel_rho(ev, V, float(tau))
                rho_t_flat = rho_t[iu, ju]
                ee = rho_t_flat[pair_masks["ee_cp"]]
                nn = rho_t_flat[pair_masks["nn_cp"]]
                sweep_rows.append(dict(
                    patient=pat, band=band, phase=phase,
                    tau_idx=int(tau_idx),
                    tau=float(tau),
                    tau_norm=float(tau * lam_max),
                    R_rho_cp=_ratio(ee, nn),
                ))

    return m1_rows, m2_rows, sweep_rows


# ---------------------------------------------------------------------------
# M3 driver: MRCA depth + KC on epi-induced subtree
# ---------------------------------------------------------------------------

def compute_mrca_and_kc(masks: PatientMasks) -> tuple[list[dict], list[dict]]:
    """Return (MRCA rows, KC epi-subtree rows)."""
    pat = masks.patient
    pair_masks = _pair_class_masks(masks.epi_mask, masks.probes)

    mrca_rows: list[dict] = []
    kc_rows: list[dict] = []

    # cache (m, M) per (band, phase)
    kc_cache: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}

    for band in BANDS:
        for phase in PHASES:
            r = load_lrg_result(pat, phase, band, fc_method="imcoh_abs")
            if r is None:
                continue
            Z = np.asarray(r.linkage_matrix)
            m, M = kc_vectors(Z)
            kc_cache[(band, phase)] = (m, M)

            for cls in ("ee", "nn", "cross"):
                msk = pair_masks[f"{cls}_all"]
                vals = M[msk]
                if vals.size == 0:
                    continue
                mrca_rows.append(dict(
                    patient=pat, band=band, phase=phase, pair_class=cls,
                    n_pairs=int(vals.size),
                    mrca_height_median=float(np.median(vals)),
                    mrca_height_q25=float(np.percentile(vals, 25)),
                    mrca_height_q75=float(np.percentile(vals, 75)),
                ))

        # KC on epi-induced subtree, three phase pairs × λ
        if masks.epi_mask.sum() < 2:
            continue
        ee_subset = pair_masks["ee_all"]
        for (a, b) in [("rest_pre", "task_test"),
                       ("task_test", "rest_post"),
                       ("rest_pre", "rest_post")]:
            if (band, a) not in kc_cache or (band, b) not in kc_cache:
                continue
            m_a, M_a = kc_cache[(band, a)]
            m_b, M_b = kc_cache[(band, b)]
            mmax = max(m_a.max(), m_b.max(), 1)
            Mmax = max(M_a.max(), M_b.max(), 1e-12)
            m_an, m_bn = m_a / mmax, m_b / mmax
            M_an, M_bn = M_a / Mmax, M_b / Mmax
            for lam in KC_LAMBDAS:
                v_a = (1 - lam) * m_an + lam * M_an
                v_b = (1 - lam) * m_bn + lam * M_bn
                d_full = float(np.linalg.norm(v_a - v_b))
                d_epi = float(np.linalg.norm(v_a[ee_subset] - v_b[ee_subset]))
                kc_rows.append(dict(
                    patient=pat, band=band, phase_a=a, phase_b=b,
                    lam=lam,
                    d_kc_full=d_full,
                    d_kc_epi=d_epi,
                    n_epi_pairs=int(ee_subset.sum()),
                ))

    return mrca_rows, kc_rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    m1_all, m2_all, sweep_all = [], [], []
    mrca_all, kc_all = [], []

    for pat in COHORT:
        try:
            masks = _build_masks(pat)
        except Exception as e:
            print(f"[audit_48] skip {pat}: cannot build masks ({e})")
            continue
        if masks.epi_mask.sum() == 0:
            print(f"[audit_48] skip {pat}: no epi annotation")
            continue
        print(f"[audit_48] {pat}: n_epi={int(masks.epi_mask.sum())}, "
              f"n_probes={len(set(masks.probes[masks.epi_mask]))}")
        m1_rows, m2_rows, sweep_rows = compute_enrichment_for_patient(masks)
        mrca_rows, kc_rows = compute_mrca_and_kc(masks)
        m1_all.extend(m1_rows)
        m2_all.extend(m2_rows)
        sweep_all.extend(sweep_rows)
        mrca_all.extend(mrca_rows)
        kc_all.extend(kc_rows)

    pd.DataFrame(m1_all).to_csv(OUT_DIR / "M1_raw_imcoh_enrichment.csv", index=False)
    pd.DataFrame(m2_all).to_csv(OUT_DIR / "M2_lrg_enrichment.csv", index=False)
    pd.DataFrame(sweep_all).to_csv(OUT_DIR / "M2b_tau_sweep.csv", index=False)
    pd.DataFrame(mrca_all).to_csv(OUT_DIR / "M3_mrca_height.csv", index=False)
    pd.DataFrame(kc_all).to_csv(OUT_DIR / "M3_kc_epi_subtree.csv", index=False)

    print(f"\n[audit_48] CSVs written to {OUT_DIR}")
    print("  M1: raw |ImCoh| enrichment per (pat, band, phase, cls)")
    print("  M2: LRG ρ̂(τ_max) + cophenetic-ultrametric enrichment")
    print("  M2b: τ-sweep of R_ρ_cp")
    print("  M3 MRCA: MRCA merge-height per (pat, band, phase, pair_class)")
    print("  M3 KC : epi-subtree KC distance per (pat, band, phase_pair, λ)")


if __name__ == "__main__":
    main()
