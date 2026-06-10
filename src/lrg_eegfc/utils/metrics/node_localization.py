"""Node-level signed localization contributions from LRG trace / anchor caches.

Maps two cohort LRG probes to per-node signed scores keyed to Desikan-Killiany
(DK) regions, for anatomical-localization testing. General graph/network
primitive: every builder returns ``(values, node_idx, node_region)`` consumable
by any region-aggregation / label-shuffle localization test (cohort-level or
per-patient). No manuscript-local scope in the names.

Two probes
----------
- **cophenet trace score** ``s_ij = dD_task * dD_rest`` — the signed,
  threshold-free continuation of the split-baseline cophenet trace flagging
  quantity (``s_ij > 0`` = sign-consistent trace direction). Each pair
  contributes its ``s`` to BOTH endpoint nodes (endpoint logic), so ``values``
  and ``node_idx`` have length ``2 * n_pairs`` and ``node_region`` is the full
  per-node DK label vector (the object a label-shuffle permutes).

- **Grassmann participation deviation** ``p_i - mean_i p_i`` — phase-AVERAGED
  eigenmode participation ``p_i(k) = sum_{j in U_k} phi_j(i)^2`` (zero mode
  skipped) averaged over (rest_pre, task_test, rest_post), then centred on the
  patient's own mean so above-baseline concentration is positive.

  WARNING: the Grassmann quantity averages over phases — it does NOT difference
  them. It is therefore an ANCHOR-flavoured "where do the leading modes live"
  quantity, not a cross-phase TRACE. Localization of this quantity answers a
  different question than the cophenet trace localization. Read accordingly.

Inputs (read-only; no FC/LRG recomputation)
-------------------------------------------
- cophenet: ``{root}/data/reports/imcoh_continuous_trace/per_pair_split/{pat}_{band}.npz``
  (keys ``dD_task, dD_rest, iu_i, iu_j``).
- grassmann: ``{root}/data/cache/imcoh_lrg/{pat}/{band}_{phase}_lrg_imcoh-abs.npz``
  (key ``eigenvectors``), phases rest_pre / task_test / rest_post.
- DK regions: ``lrg_eegfc.utils.io.regions.load_channel_regions``.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.io.patient import load_epileptic_nodes
from lrg_eegfc.utils.io.regions import load_channel_regions, _normalise_label

#: Phases averaged for the Grassmann participation (anchor) quantity.
GRASSMANN_PHASES = ("rest_pre", "task_test", "rest_post")

#: Non-anatomical region labels — white matter, unparsed/unknown contacts. NOT
#: Desikan-Killiany regions; downstream localization tests drop these from the
#: anatomical universe (they pollute it and, being sampled by every patient,
#: manufacture spurious cohort significance). Exposed here so callers share one
#: definition.
NON_ANATOMICAL = frozenset({"Wm", "Unk", "unknown"})


def cophenet_trace_contributions(pat: str, band: str, root: Path):
    """Signed per-endpoint cophenet trace contributions for one patient.

    Returns ``(values, node_idx, node_region)`` where each pair contributes its
    signed score ``s_ij = dD_task * dD_rest`` to both endpoints (length
    ``2 * n_pairs``); ``node_region`` is the full ``(N,)`` DK-label vector.
    """
    p = (Path(root) / "data/reports/imcoh_continuous_trace/per_pair_split"
         / f"{pat}_{band}.npz")
    d = np.load(p)
    dD_task = np.asarray(d["dD_task"], dtype=np.float64)
    dD_rest = np.asarray(d["dD_rest"], dtype=np.float64)
    iu_i = np.asarray(d["iu_i"]).astype(int)
    iu_j = np.asarray(d["iu_j"]).astype(int)
    s = dD_task * dD_rest                       # signed per-pair trace score
    node_region = load_channel_regions(pat)["region"].to_numpy()
    values = np.concatenate([s, s])             # each pair -> 2 endpoint contribs
    node_idx = np.concatenate([iu_i, iu_j])
    return values, node_idx, node_region


def participation_kset(eigvecs: np.ndarray, k_iter) -> np.ndarray:
    """Per-node participation ``sum_{j in U_k} phi_j(i)^2`` averaged over a k set.

    Zero mode (column 0) skipped; ``U_k`` = columns ``1..k+1``. Replicates
    ``audit_72.participation_kset`` exactly.
    """
    N = eigvecs.shape[0]
    weights = np.zeros(N, dtype=np.float64)
    n_k = 0
    for k in k_iter:
        if k + 1 > N:
            continue
        U = eigvecs[:, 1:k + 1]
        weights += (U * U).sum(axis=1)
        n_k += 1
    return weights / max(n_k, 1)


def epi_keep_mask(regions_df, pat: str) -> np.ndarray:
    """Label-based non-epileptic node keep-mask (audit_71 pattern).

    ``load_epileptic_nodes`` returns string labels, NOT integer indices, so
    ``np.isin(np.arange(N), epi)`` is an all-False silent no-op. Normalise each
    contact label and test membership in the epileptic-label set instead.
    """
    epi = set(load_epileptic_nodes(pat))
    labels = [_normalise_label(l) for l in regions_df["label_raw"]]
    return np.array([l not in epi for l in labels], dtype=bool)


def grassmann_trace_contributions(pat: str, band: str, root: Path,
                                  k_set, epi_x: bool = False):
    """Per-node CROSS-PHASE Grassmann participation TRACE for one patient.

    The trace analog of :func:`grassmann_participation_contributions` (which is a
    phase-averaged ANCHOR). Using the subspace-invariant participation
    ``p_i = sum_{j in U_k} phi_j(i)^2`` (projector diagonal onto the top-k LRG
    subspace, comparable across phases without eigenvector sign-matching),

        dP_task = p_i(task)      - p_i(rest_pre)      # task reorganization
        dP_rest = p_i(rest_post) - p_i(rest_pre)      # what persists into rest
        s_i     = dP_task * dP_rest                   # >0 = sign-consistent TRACE

    This is the exact node-level analog of the cophenet ``s_ij = dD_task*dD_rest``
    split-baseline trace score (audit_71 form) and respects the locked
    positive=trace convention. Returns ``(values, node_idx, node_region)`` with
    ``node_idx = arange(n)``.

    SCOPE NOTE: a single (unsplit) rest_pre baseline is used, so the shared
    ``-p_pre`` term gives ``s_i`` a positive shared-baseline bias. That bias is
    irrelevant to the LOCALIZATION question — the within-patient label-shuffle
    null operates on the identical ``s_i`` values, so it tests only SPATIAL
    concentration, not trace existence (trace existence is the matched-strength
    rung, already locked for beta Grassmann at §5.4). It does NOT re-verify the
    trace; it asks only whether, given the trace, it concentrates anatomically.
    """
    k_list = list(k_set)
    p = {}
    for ph in GRASSMANN_PHASES:
        fp = (Path(root) / "data/cache/imcoh_lrg" / pat
              / f"{band}_{ph}_lrg_imcoh-abs.npz")
        eig = np.asarray(np.load(fp)["eigenvectors"], dtype=np.float64)
        p[ph] = participation_kset(eig, k_list)
    dP_task = p["task_test"] - p["rest_pre"]
    dP_rest = p["rest_post"] - p["rest_pre"]
    s = dP_task * dP_rest                                   # signed node trace
    regions_df = load_channel_regions(pat)
    node_region = regions_df["region"].to_numpy()
    if epi_x:
        keep = epi_keep_mask(regions_df, pat)
        s = s[keep]
        node_region = node_region[keep]
    node_idx = np.arange(s.shape[0])
    return s, node_idx, node_region


def grassmann_participation_contributions(pat: str, band: str, root: Path,
                                          k_set, epi_x: bool = False):
    """Phase-averaged Grassmann participation deviation for one patient.

    Per-node value = (phase-averaged participation over ``k_set``) minus the
    patient's own mean participation — signed deviation from baseline, so
    above-baseline concentration is positive. If ``epi_x``, participation is
    computed in full-graph spectral space then masked to non-epileptic nodes
    (BUGFIXED label-based mask via :func:`epi_keep_mask`).

    Returns ``(values, node_idx, node_region)`` with ``node_idx = arange(n)``
    (no endpoint duplication for the node-level Grassmann quantity).
    """
    k_list = list(k_set)
    eig = {}
    for ph in GRASSMANN_PHASES:
        fp = (Path(root) / "data/cache/imcoh_lrg" / pat
              / f"{band}_{ph}_lrg_imcoh-abs.npz")
        eig[ph] = np.asarray(np.load(fp)["eigenvectors"], dtype=np.float64)
    parts = [participation_kset(eig[ph], k_list) for ph in GRASSMANN_PHASES]
    part = np.mean(np.stack(parts, axis=0), axis=0)        # phase-avg, (N,)
    regions_df = load_channel_regions(pat)
    node_region = regions_df["region"].to_numpy()
    if epi_x:
        keep = epi_keep_mask(regions_df, pat)
        part = part[keep]
        node_region = node_region[keep]
    dev = part - part.mean()                                # signed deviation
    node_idx = np.arange(part.shape[0])
    return dev, node_idx, node_region


# ---------------------------------------------------------------------------
# spatial-grouping primitives for localization / concentration tests
# (promoted from diag_per_patient_concentration_control.py — ≥2 callers:
#  the shaft control + the localization atlas). General names, no manuscript scope.
# ---------------------------------------------------------------------------

_SHAFT_RE = re.compile(r"\s*([A-Za-z]+'?)")


def shaft_of(label_raw: str) -> str:
    """Electrode-shaft id = leading letters (+ optional prime) of a clinical label.

    ``'A 1,G2'`` → ``'A'``; ``"B' 3,G2"`` → ``"B'"``. A shaft is a purely spatial
    unit (where the surgeon placed the depth electrode); contacts on one shaft are
    spatially adjacent and measure near-identical signal.
    """
    m = _SHAFT_RE.match(str(label_raw))
    return m.group(1) if m else str(label_raw)


def eta2(values: np.ndarray, group_codes: np.ndarray, n_codes: int) -> float:
    """Between-group variance fraction (size-weighted) of ``values``.

    ``eta^2 = Σ_g n_g (mean_g − grand_mean)^2 / Σ_i (x_i − grand_mean)^2`` ∈ [0,1].
    """
    grand = values.mean()
    sums = np.bincount(group_codes, weights=values, minlength=n_codes)
    counts = np.bincount(group_codes, minlength=n_codes).astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        means = np.where(counts > 0, sums / counts, grand)
    between = float(np.sum(counts * (means - grand) ** 2))
    total = float(np.sum((values - grand) ** 2))
    return between / total if total > 0 else 0.0


def perm_p_eta2(values, label_vec, node_idx, n_perm, rng):
    """Floor-free permutation p (greater) for ``eta2`` under a label shuffle.

    ``label_vec`` is the length-N node→group label vector (region / shaft /
    system); contributions are grouped by ``label_vec[node_idx]``. The shuffle
    permutes ``label_vec`` (preserves the group-size multiset), matching the
    localization null exactly. Returns ``(obs_eta2, perm_p, n_groups)``.
    """
    uniq, codes = np.unique(label_vec, return_inverse=True)
    R = uniq.size
    obs = eta2(values, codes[node_idx], R)
    ge = 0
    for _ in range(n_perm):
        cs = rng.permutation(codes)
        if eta2(values, cs[node_idx], R) >= obs:
            ge += 1
    return obs, (1 + ge) / (n_perm + 1), R


def nmi(a_labels: np.ndarray, b_labels: np.ndarray) -> float:
    """Normalized mutual information between two node partitions (sklearn-free).

    ``NMI = I(A;B) / sqrt(H(A) H(B))`` ∈ [0,1]; 1.0 = identical partitions
    (collinear, e.g. regions that exactly follow shafts).
    """
    a_u, a = np.unique(a_labels, return_inverse=True)
    b_u, b = np.unique(b_labels, return_inverse=True)
    n = a.size
    Pa = np.bincount(a) / n
    Pb = np.bincount(b) / n
    joint = np.zeros((a_u.size, b_u.size))
    for i, j in zip(a, b):
        joint[i, j] += 1
    joint /= n
    Ha = -np.sum(Pa[Pa > 0] * np.log(Pa[Pa > 0]))
    Hb = -np.sum(Pb[Pb > 0] * np.log(Pb[Pb > 0]))
    nz = joint > 0
    I = float(np.sum(joint[nz] * np.log(joint[nz] / np.outer(Pa, Pb)[nz])))
    return I / np.sqrt(Ha * Hb) if Ha > 0 and Hb > 0 else 0.0
