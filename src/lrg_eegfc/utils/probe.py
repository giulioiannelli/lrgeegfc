"""Probe-level analysis utilities for sEEG electrode data.

Provides reusable functions for:

* Extracting probe identifiers from channel labels
* Building same-probe boolean masks
* Computing same-probe vs cross-probe weight ratios
* Computing community-level probe enrichment at arbitrary LRG scales

All functions are FC-method agnostic — they operate on generic adjacency
matrices or community partitions.
"""

from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.utils.io.patient import parse_seeg_label

__all__ = [
    "probe_from_label",
    "extract_probe_labels",
    "contact_labels",
    "split_label",
    "build_probe_mask",
    "compute_probe_weight_ratio",
    "compute_community_probe_enrichment",
    "compute_enrichment_vs_scale",
]


# ---------------------------------------------------------------------------
# Probe extraction
# ---------------------------------------------------------------------------

def probe_from_label(label: str) -> str:
    """Extract the probe identifier from an sEEG channel label.

    Delegates to :func:`~lrg_eegfc.utils.io.patient.parse_seeg_label` which
    handles spaces, primes, reference suffixes, and typos.  Falls back to a
    simple regex if the parser returns ``None``.

    Parameters
    ----------
    label : str
        Raw channel label (e.g. ``"A 1,G2"``, ``"G' 3"``, ``"A1"``).

    Returns
    -------
    str
        Probe identifier (e.g. ``"A"``, ``"G'"``).
    """
    import re

    probe, _ = parse_seeg_label(label)
    if probe is not None:
        return probe
    # Fallback for unusual labels
    cleaned = label.strip().replace("\u00ec", "'")
    # Handle bipolar labels like 'A2-A1'
    cleaned = cleaned.split("-")[0]
    m = re.match(r"([A-Za-z]+'?)", cleaned)
    return m.group(1) if m else cleaned


def extract_probe_labels(channel_labels: Sequence[str]) -> list[str]:
    """Convert a sequence of channel labels to their probe identifiers.

    Parameters
    ----------
    channel_labels : sequence of str
        Raw channel labels.

    Returns
    -------
    list of str
        Probe identifier for each channel.
    """
    return [probe_from_label(l) for l in channel_labels]


def contact_labels(probes: Sequence[str], n: int) -> List[str]:
    """Build sEEG contact names (shaft identifier + running index within shaft).

    Given the per-contact probe identifiers in FC/matrix order, assign each
    contact a name ``"{probe}{k}"`` where *k* is its 1-based position within
    that probe (e.g. ``["A", "A", "G"] -> ["A1", "A2", "G1"]``).

    Parameters
    ----------
    probes : sequence of str
        Probe identifier per contact (as returned by
        :func:`extract_probe_labels`), FC-aligned.
    n : int
        Number of contacts to label (leading ``n`` entries of *probes*).

    Returns
    -------
    list of str
        Contact name per contact.
    """
    counts: dict[str, int] = {}
    labels: List[str] = []
    for p in probes[:n]:
        counts[p] = counts.get(p, 0) + 1
        labels.append(f"{p}{counts[p]}")
    return labels


def split_label(label: str) -> tuple[str, str]:
    """Split a contact name into (shaft, index): ``"A10" -> ("A", "10")``.

    The shaft is the leading non-digit run (including a prime, e.g. ``"G'"``),
    the index the trailing digit run (rendered as a superscript by callers).
    Falls back to ``(label, "")`` if the pattern does not match.
    """
    import re

    m = re.match(r"^(\D*)(\d*)$", label)
    return (m.group(1), m.group(2)) if m else (label, "")


# ---------------------------------------------------------------------------
# Probe mask
# ---------------------------------------------------------------------------

def build_probe_mask(channel_labels: Sequence[str]) -> NDArray:
    """Build a boolean mask where ``True`` = both contacts on the same probe.

    The diagonal is set to ``False`` (a contact is not its own pair).

    Parameters
    ----------
    channel_labels : sequence of str
        Channel labels of length *N*.

    Returns
    -------
    NDArray, shape (N, N)
        Symmetric boolean mask.
    """
    probes = extract_probe_labels(channel_labels)
    N = len(probes)
    probe_arr = np.array(probes)
    mask = probe_arr[:, None] == probe_arr[None, :]
    np.fill_diagonal(mask, False)
    return mask


# ---------------------------------------------------------------------------
# Weight-ratio statistics
# ---------------------------------------------------------------------------

def compute_probe_weight_ratio(
    W: NDArray,
    channel_labels: Sequence[str],
) -> float:
    """Compute the same-probe / cross-probe mean absolute weight ratio.

    A ratio of 1.0 means no bias; values > 1 indicate same-probe inflation.

    Parameters
    ----------
    W : NDArray, shape (N, N)
        Symmetric adjacency matrix (any FC method).
    channel_labels : sequence of str
        Channel labels (length N).

    Returns
    -------
    float
        ``mean(|W[same-probe]|) / mean(|W[cross-probe]|)``.
        Returns 0.0 if either set is empty.
    """
    mask = build_probe_mask(channel_labels)
    A = np.abs(W.copy())
    np.fill_diagonal(A, 0.0)

    off_diag = ~np.eye(A.shape[0], dtype=bool)
    sp_vals = A[mask]
    cp_vals = A[off_diag & ~mask]

    if sp_vals.size == 0 or cp_vals.size == 0:
        return 0.0
    return float(sp_vals.mean() / (cp_vals.mean() + 1e-30))


def probe_weight_distributions(
    W: NDArray,
    channel_labels: Sequence[str],
) -> tuple[NDArray, NDArray]:
    """Return same-probe and cross-probe weight arrays (upper triangle only).

    Parameters
    ----------
    W : NDArray, shape (N, N)
        Symmetric adjacency matrix.
    channel_labels : sequence of str
        Channel labels.

    Returns
    -------
    same_probe : NDArray
        Absolute weights for same-probe pairs.
    cross_probe : NDArray
        Absolute weights for cross-probe pairs.
    """
    mask = build_probe_mask(channel_labels)
    A = np.abs(W.copy())
    np.fill_diagonal(A, 0.0)

    # Upper triangle only to avoid double-counting
    triu = np.triu_indices_from(A, k=1)
    triu_mask = mask[triu]
    triu_vals = A[triu]

    return triu_vals[triu_mask], triu_vals[~triu_mask]


# ---------------------------------------------------------------------------
# Community-level probe enrichment
# ---------------------------------------------------------------------------

def compute_community_probe_enrichment(
    community_labels: NDArray,
    probe_labels: Sequence[str],
) -> float:
    """Fraction of same-community pairs that are same-probe, normalised.

    Returns the enrichment ratio:

    .. math::

        \\text{enrichment} = \\frac{P(\\text{same-probe} \\mid \\text{same-community})}
                                    {P(\\text{same-probe})}

    A value of 1.0 means communities are unbiased by probe geometry.
    Values > 1 indicate that same-probe pairs are over-represented within
    communities.

    Parameters
    ----------
    community_labels : NDArray, shape (N,)
        Integer community assignment for each node.
    probe_labels : sequence of str
        Probe identifier for each node (length N).

    Returns
    -------
    float
        Enrichment ratio (1.0 = no bias, > 1.0 = biased).
    """
    N = len(community_labels)
    if N < 2:
        return 0.0

    comm = np.asarray(community_labels)
    probes = np.array(list(probe_labels))

    # Vectorised: build upper-triangle pair masks
    i_idx, j_idx = np.triu_indices(N, k=1)
    same_comm = comm[i_idx] == comm[j_idx]
    same_probe = probes[i_idx] == probes[j_idx]

    n_same_comm = same_comm.sum()
    if n_same_comm == 0:
        return 0.0

    frac_sp_in_sc = (same_comm & same_probe).sum() / n_same_comm
    expected = same_probe.sum() / len(i_idx)
    return float(frac_sp_in_sc / expected) if expected > 0 else 0.0


def compute_enrichment_vs_scale(
    linkage_matrix: NDArray,
    probe_labels: Sequence[str],
    n_communities: Sequence[int] = (3, 5, 10, 15, 20, 30),
    n_nodes: Optional[int] = None,
) -> dict[int, float]:
    """Compute community-probe enrichment at multiple LRG scales.

    Parameters
    ----------
    linkage_matrix : NDArray
        Scipy-format linkage matrix from LRG analysis.
    probe_labels : sequence of str
        Probe identifier per node (length >= *n_nodes*).
    n_communities : sequence of int
        Community counts (scales) to evaluate.
    n_nodes : int, optional
        Number of nodes in the giant component.  If ``None``, inferred from
        the linkage matrix as ``linkage.shape[0] + 1``.

    Returns
    -------
    dict[int, float]
        Mapping ``n_communities → enrichment_ratio``.
    """
    if n_nodes is None:
        n_nodes = linkage_matrix.shape[0] + 1
    pl = list(probe_labels[:n_nodes])

    result: dict[int, float] = {}
    for nc in n_communities:
        if nc >= n_nodes:
            continue
        labels = fcluster(linkage_matrix, nc, criterion="maxclust")
        result[nc] = compute_community_probe_enrichment(labels, pl)
    return result
