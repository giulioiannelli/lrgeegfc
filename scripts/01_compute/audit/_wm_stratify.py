"""Shared helpers for the white-matter-node-stratified trace audits.

Private sibling of ``_epi_stratify.py`` (cf. its docstring). Imported by
``audit_83_wm_stratified_cophenetic`` and ``audit_84_wm_stratified_grassmann``
so the WM configuration vocabulary, node/pair masking, surrogate-cache
namespacing and per-cell RNG live in exactly one place.

This module is the *white-matter analogue* of ``_epi_stratify``: the ONLY
substantive change is the node mask. Everything statistical — the cohort
paired-Wilcoxon ``cohort_verdict``, the persist/weaken/emerge/absent
``sensitivity_flag``, the condensed upper-triangle index order — is *imported*
from ``_epi_stratify`` (reuse, never a fork) because those helpers are
stratification-agnostic. The matched-strength eigendecomposition machinery
comes from ``lrg_eegfc.utils.surrogate``.

WM definition (locked by the user, 2026-06-08)
----------------------------------------------
A contact is **white matter** iff its *dominant* Desikan-Killiany tissue is
``Wm`` — i.e. ``load_channel_regions(pat)["region"] == "Wm"`` (argmax of the
per-tissue weights, PTD excluded; see ``utils/io/regions.py``). This is the
atlas-dominant criterion, consistent with ``regions._NON_ANATOMICAL_REGIONS``.
Cohort prevalence is large: 30–57 % of contacts per patient (Pat_14 lowest,
Pat_15 highest, cohort ≈ 40 %), so ``exclude_wm`` is a structural change, not a
fine perturbation. The mask is FC-index-ordered (``channel_labels.csv`` row
order == FC / dendrogram-leaf order) and length-checked against the FC ``N``,
the same alignment contract ``epi_mask_for`` enforces.

Two families of "WM view" on the cross-phase trace (mirroring epi)
------------------------------------------------------------------
- **Subgraph configs** rebuild the LRG on a node submatrix ``W[ix(keep,keep)]``:
  ``full`` (all nodes), ``exclude_wm`` (gray-only — the headline "remove WM
  from the timeseries" view), ``wm_only`` (white-matter contacts alone).
  Because ImCoh / corr / MSC are *pairwise* spectral estimators, the FC on a
  gray-only montage equals the gray-only submatrix of the cached full FC to
  numerical precision — so the submatrix path faithfully realizes "drop the
  WM channels before computing FC", and only the (node-coupled) LRG step
  genuinely changes. The matched-strength surrogate is regenerated on the
  submatrix in its own cache directory (so an ``N``-mismatched ensemble never
  clobbers the canonical full-graph one).

- **Pair-class configs** keep the *full-graph* LRG and restrict the *pair set*:
  ``gray_gray`` (both endpoints non-WM), ``cross`` (gray↔WM interface),
  ``wm_wm`` (both endpoints WM). The surrogate is the full-graph
  matched-strength ensemble with the per-pair statistic recomputed on the
  restricted pair set (the audit_71 ``pair_keep`` mechanism). Pair-class is
  cophenetic/raw only — the Grassmann subspace statistic has no
  pair-restriction analogue (so audit_84 is subgraph-only).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.surrogate import surrogate_cache_path

# Reuse the proven, stratification-agnostic machinery (import, no fork). These
# functions know nothing about epilepsy — they are cohort-statistics helpers
# that happen to live in the epi sibling module (first writer wins). audit_83/84
# put the audit dir on sys.path before importing this module, so `_epi_stratify`
# resolves here too.
import _epi_stratify as _es  # type: ignore

cohort_verdict = _es.cohort_verdict
sensitivity_flag = _es.sensitivity_flag
upper_tri_indices = _es.upper_tri_indices


# ---------------------------------------------------------------------------
# Vocabulary (single source of truth)
# ---------------------------------------------------------------------------
COHORT = list(_es.COHORT)
ALL_BANDS = list(_es.ALL_BANDS)
PHASES_4 = _es.PHASES_4
PHASES_3 = _es.PHASES_3

SUBGRAPH_CONFIGS = ("full", "exclude_wm", "wm_only")
PAIRCLASS_CONFIGS = ("gray_gray", "cross", "wm_wm")
ALL_CONFIGS = SUBGRAPH_CONFIGS + PAIRCLASS_CONFIGS

CONFIG_CLASS = {**{c: "subgraph" for c in SUBGRAPH_CONFIGS},
                **{c: "pairclass" for c in PAIRCLASS_CONFIGS}}

# Honest hard floors (surfaced downstream as defined=False, never silent).
MIN_NODES = _es.MIN_NODES   # an LRG on < 3 leaves is degenerate
MIN_PAIRS = _es.MIN_PAIRS   # a Spearman over < 10 pairs is uninterpretable

N_SURROGATES = _es.N_SURROGATES
SWAP_FACTOR = _es.SWAP_FACTOR

# full + pair-class reuse the canonical seed-20260511 full-graph ensemble
# (audit_63/66 cache → guaranteed cache hits). The two WM submatrix configs get
# a fresh ensemble keyed on a WM seed in their own directories.
SEED_CANONICAL = 20260511
SEED_WM = 20260608


# ---------------------------------------------------------------------------
# WM mask (FC-index-ordered; same alignment contract as epi_mask_for)
# ---------------------------------------------------------------------------
def wm_mask_for(pat: str, n_expected: int | None = None) -> np.ndarray:
    """FC-index-ordered boolean white-matter mask (dominant region == ``Wm``).

    Raises if its length disagrees with the FC matrix size — the same
    contract audit_67/_epi_stratify enforce, so a montage/atlas misalignment
    can never silently degrade into an all-False no-op.
    """
    reg = load_channel_regions(pat)
    mask = (reg["region"].astype(str) == "Wm").to_numpy()
    if n_expected is not None and mask.size != n_expected:
        raise ValueError(
            f"{pat}: wm_mask length {mask.size} != FC N {n_expected}"
        )
    return mask


# ---------------------------------------------------------------------------
# Node / pair masks per config
# ---------------------------------------------------------------------------
def node_mask_for_config(wm_mask: np.ndarray, config: str,
                         min_nodes: int = MIN_NODES) -> np.ndarray | None:
    """Boolean node-keep mask for a SUBGRAPH config. ``None`` if structurally
    undefined for this patient (below ``min_nodes``)."""
    if config == "full":
        return np.ones_like(wm_mask, dtype=bool)
    if config == "exclude_wm":
        keep = ~wm_mask
    elif config == "wm_only":
        keep = wm_mask.copy()
    else:
        raise ValueError(f"not a subgraph config: {config}")
    return keep if int(keep.sum()) >= min_nodes else None


def pair_keep_for_config(wm_mask: np.ndarray, config: str,
                         min_pairs: int = MIN_PAIRS) -> np.ndarray | None:
    """Boolean pair-keep mask (condensed upper-tri order) for a PAIRCLASS
    config. ``None`` if the class has fewer than ``min_pairs`` pairs."""
    n = wm_mask.size
    iu_i, iu_j = upper_tri_indices(n)
    wi, wj = wm_mask[iu_i], wm_mask[iu_j]
    if config == "wm_wm":
        pk = wi & wj
    elif config == "gray_gray":
        pk = (~wi) & (~wj)
    elif config == "cross":
        pk = wi ^ wj            # exactly one endpoint white matter
    else:
        raise ValueError(f"not a pair-class config: {config}")
    return pk if int(pk.sum()) >= min_pairs else None


# ---------------------------------------------------------------------------
# Surrogate-eig cache path per config (no clobber across configs / strata)
# ---------------------------------------------------------------------------
SURR_DIR_EXCLUDE_WM = CACHE_ROOT / "matched_strength_surrogate_wm_excluded_lrg"
SURR_DIR_WM_ONLY = CACHE_ROOT / "matched_strength_surrogate_wm_only_lrg"


def config_seed(config: str) -> int:
    return (SEED_CANONICAL if (config == "full" or config in PAIRCLASS_CONFIGS)
            else SEED_WM)


def surr_eig_path(config: str, pat: str, band: str, phase: str,
                  n_surr: int = N_SURROGATES, swap_factor: int = SWAP_FACTOR,
                  ) -> Path:
    """Cache path for the matched-strength eigendecomposition ensemble a config
    needs. ``full``/pair-class reuse the canonical full-graph cache; the two WM
    submatrix configs get their own ``wmX`` / ``wmONLY`` directories."""
    seed = config_seed(config)
    if config == "full" or config in PAIRCLASS_CONFIGS:
        return surrogate_cache_path(pat, band, phase, n_surr, swap_factor, seed)
    if config == "exclude_wm":
        return (SURR_DIR_EXCLUDE_WM / pat
                / f"{band}_{phase}_wmX_R{n_surr}_swap{swap_factor}"
                  f"_seed{seed}_imcoh_abs.npz")
    if config == "wm_only":
        return (SURR_DIR_WM_ONLY / pat
                / f"{band}_{phase}_wmONLY_R{n_surr}_swap{swap_factor}"
                  f"_seed{seed}_imcoh_abs.npz")
    raise ValueError(config)


# ---------------------------------------------------------------------------
# Deterministic per-cell RNG (cache reproducibility regardless of loop order)
# ---------------------------------------------------------------------------
_PHASE_IDX = {p: i for i, p in enumerate(PHASES_4)}


def cell_rng(pat: str, band: str, phase: str, config: str) -> np.random.Generator:
    """Reproducible generator seeded by ``(seed, config, patient, band, phase)``
    so a freshly-computed ensemble is identical on re-run regardless of
    iteration order. Phase index is taken over ``PHASES_4`` so audit_83 (4-phase)
    and audit_84 (3-phase Grassmann) agree on any shared cell. The ``+ 600``
    config offset keeps the WM seed-sequence path distinct from the epi one even
    where ``config_seed`` coincides (full / pair-class)."""
    pat_n = int(pat.split("_")[-1])
    ss = np.random.SeedSequence([
        config_seed(config), 600 + ALL_CONFIGS.index(config), pat_n,
        ALL_BANDS.index(band), _PHASE_IDX[phase],
    ])
    return np.random.default_rng(ss)
