"""Spatial-scale summaries of flat partitions on physical coordinates.

Bridges the graph-theoretic partition (cluster labels) and the brain's
physical geometry (electrode coordinates). Asks: at this partition
resolution, what is the typical *physical* size / spacing of a module?
The output unit follows the input ``coords`` unit (millimetres,
micrometres, etc.). Used to translate dendrogram-cut scale into a
biologically interpretable spatial axis.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist, pdist

__all__ = ["SUMMARY_NAMES", "cluster_spatial_scale"]

SUMMARY_NAMES = (
    "mean_diameter",
    "mean_radius",
    "nn_centroid",
    "size_weighted_radius",
)


def cluster_spatial_scale(
    labels: np.ndarray,
    coords: np.ndarray,
    *,
    summary: str = "mean_diameter",
    drop_singletons: bool = False,
) -> dict:
    """Compute four spatial summaries of a flat partition in one pass.

    Parameters
    ----------
    labels : array (N,)
        Integer cluster labels for the ``N`` electrodes.
    coords : array (N, 3)
        Electrode coordinates in any consistent unit (mm preferred).
    summary : str
        Selects which summary populates the ``scale`` field of the result.
        One of ``mean_diameter`` (default), ``mean_radius``, ``nn_centroid``,
        ``size_weighted_radius``. All four are computed and returned in
        ``all`` regardless of this choice.
    drop_singletons : bool
        If True, exclude clusters of size 1 from the per-cluster aggregates
        (their diameter and radius are 0; their nearest-neighbour distance
        depends only on centroid placement). Default False.

    Returns
    -------
    dict with keys
        ``scale`` — the selected summary (float, in coords units).
        ``all`` — dict of the four summaries:
            * ``mean_diameter``: mean over clusters of max pairwise distance
              within the cluster (cluster ``size`` if cluster has > 1 leaf;
              0 for singletons).
            * ``mean_radius``: mean over clusters of mean Euclidean distance
              from each cluster member to its centroid.
            * ``nn_centroid``: mean over clusters of the distance to the
              nearest other cluster centroid (NaN if only one cluster).
            * ``size_weighted_radius``: ``Σ_l |C_l| · radius(C_l) / N``.
        ``n_clusters`` — number of distinct labels.
        ``n_singletons`` — count of clusters of size 1.
        ``n_leaves`` — total electrodes ``N``.
        ``per_cluster`` — dict with arrays
            ``label, size, centroid (n_clusters, 3), diameter, radius,
            nn_distance``.

    Notes
    -----
    All four summaries are computed in one pass to avoid recomputing the
    centroids and pairwise distances; the ``summary`` argument only controls
    which scalar is featured in ``scale``. Complexity: ``O(N^2 / k)``
    average per (``k`` clusters) due to the within-cluster pdist; for
    typical ``N ≈ 120`` this is sub-millisecond.
    """
    if summary not in SUMMARY_NAMES:
        raise ValueError(
            f"summary must be one of {SUMMARY_NAMES}; got {summary!r}"
        )
    labels = np.asarray(labels)
    coords = np.asarray(coords, dtype=float)
    if labels.ndim != 1:
        raise ValueError("labels must be 1-D")
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coords must be shape (N, 3)")
    if labels.size != coords.shape[0]:
        raise ValueError("labels and coords disagree on N")
    if not np.all(np.isfinite(coords)):
        raise ValueError("coords contains non-finite values")

    unique, counts = np.unique(labels, return_counts=True)
    n_clusters = len(unique)
    n_singletons = int((counts == 1).sum())
    n_leaves = int(labels.size)

    diameters = np.zeros(n_clusters)
    radii = np.zeros(n_clusters)
    centroids = np.zeros((n_clusters, 3))
    sizes = counts.astype(int)

    for i, lab in enumerate(unique):
        members = coords[labels == lab]
        centroids[i] = members.mean(axis=0)
        if members.shape[0] > 1:
            diameters[i] = pdist(members).max()
            radii[i] = np.linalg.norm(members - centroids[i], axis=1).mean()
        # singletons: diameter = radius = 0 by definition

    # nearest-neighbour centroid distance per cluster
    if n_clusters > 1:
        cd = cdist(centroids, centroids)
        np.fill_diagonal(cd, np.inf)
        nn_dists = cd.min(axis=1)
    else:
        nn_dists = np.full(n_clusters, np.nan)

    # apply singleton drop to the aggregates that include them
    if drop_singletons:
        mask = sizes > 1
        if not mask.any():
            mean_diameter = float("nan")
            mean_radius = float("nan")
            nn_centroid = float("nan")
        else:
            mean_diameter = float(diameters[mask].mean())
            mean_radius = float(radii[mask].mean())
            nn_centroid = (float(nn_dists[mask].mean())
                           if np.all(np.isfinite(nn_dists[mask])) else float("nan"))
    else:
        mean_diameter = float(diameters.mean())
        mean_radius = float(radii.mean())
        nn_centroid = (float(nn_dists.mean())
                       if np.all(np.isfinite(nn_dists)) else float("nan"))
    # size-weighted radius is well-defined regardless
    size_weighted_radius = (
        float((sizes * radii).sum() / n_leaves) if n_leaves > 0 else float("nan")
    )

    summaries = {
        "mean_diameter": mean_diameter,
        "mean_radius": mean_radius,
        "nn_centroid": nn_centroid,
        "size_weighted_radius": size_weighted_radius,
    }
    return {
        "scale": summaries[summary],
        "all": summaries,
        "n_clusters": n_clusters,
        "n_singletons": n_singletons,
        "n_leaves": n_leaves,
        "per_cluster": {
            "label": unique,
            "size": sizes,
            "centroid": centroids,
            "diameter": diameters,
            "radius": radii,
            "nn_distance": nn_dists,
        },
    }
