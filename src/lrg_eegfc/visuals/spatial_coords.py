"""Spatial coordinate loaders + MNI-space transforms.

Split out of `spatial.py` on 2026-05-29 (Phase 4-B split 6/7). Holds the
electrode-implant metadata loader (`load_spatial_metadata`), the
coordinate dispatcher / unit converter (`prepare_spatial_coordinates`),
and the MNI affine estimator (`_estimate_mni_transform`).

Re-exported from `spatial.py` for backwards compatibility with callers
that still do `from lrg_eegfc.visuals.spatial import load_spatial_metadata`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from lrg_eegfc.config.paths import SEEG_DATAPATH


__all__ = [
    "load_spatial_metadata",
    "prepare_spatial_coordinates",
]


def _normalize_label(label: str) -> str:
    """Normalize electrode label for matching.

    Handles differences like 'A1' vs 'A 1' vs 'A 1,G2' vs "G'1" vs "G' 1".
    """
    import re

    # Remove quotes and strip whitespace
    s = str(label).strip().strip('"').strip("'")
    # Remove suffix like ',G2'
    if "," in s:
        s = s.split(",")[0]
    # Remove ALL spaces (A 1 -> A1, G' 1 -> G'1)
    s = s.replace(" ", "")
    # Normalize prime/accent characters
    s = s.replace("ì", "'")  # Replace accented i with prime
    return s.strip()


def load_spatial_metadata(
    patient: str,
    dataset_root: Path,
) -> pd.DataFrame:
    """Load channel metadata with spatial coordinates.

    Handles label mismatches between channel_labels.csv and Implant_pat_*.csv
    by normalizing labels before merging.

    Parameters
    ----------
    patient : str
        Patient identifier (e.g., "Pat_02").
    dataset_root : Path
        Root directory containing patient folders.

    Returns
    -------
    metadata : pd.DataFrame
        DataFrame with columns: label, x, y, z, and optionally atlas info.

    Raises
    ------
    FileNotFoundError
        If required files are not found.
    """
    patient_path = Path(dataset_root) / patient

    # Find implant file
    try:
        patnum = int(patient.split("_")[-1])
    except ValueError:
        patnum = 0
    implant_csv = patient_path / f"implant_pat_{patnum:02d}.csv"
    if not implant_csv.exists():
        implant_csv = patient_path / f"Implant_pat_{patnum:02d}.csv"

    channel_labels_csv = patient_path / "channel_labels.csv"

    if not implant_csv.exists():
        raise FileNotFoundError(f"Implant file not found: {implant_csv}")
    if not channel_labels_csv.exists():
        raise FileNotFoundError(f"Channel labels not found: {channel_labels_csv}")

    # Load files
    # Implant file may have quoted Desikan-Killany field with internal commas
    # Only read the columns we need — try both column name variants
    try:
        implant = pd.read_csv(implant_csv, usecols=["label", "x", "y", "z", "Desikan-Killany"])
    except ValueError:
        # Some patients use "Desikan-K" instead of "Desikan-Killany"
        implant = pd.read_csv(implant_csv, usecols=["label", "x", "y", "z", "Desikan-K"])
        implant = implant.rename(columns={"Desikan-K": "Desikan-Killany"})

    # Channel labels file may or may not have a header
    # Also may have extra columns (e.g., "F 1,G2" format)
    labels = pd.read_csv(channel_labels_csv)
    if "label" not in labels.columns:
        # No header - read again, use first column as label
        labels = pd.read_csv(channel_labels_csv, header=None)
        labels = labels.rename(columns={labels.columns[0]: "label"})
        # Keep only the label column
        labels = labels[["label"]]

    # Normalize labels for matching
    implant["_norm_label"] = implant["label"].apply(_normalize_label)
    labels["_norm_label"] = labels["label"].apply(_normalize_label)

    # Merge on normalized labels
    metadata = labels.merge(
        implant[["_norm_label", "x", "y", "z", "Desikan-Killany"]],
        on="_norm_label",
        how="left",
    )

    # Clean up coordinate columns (handle comma decimals)
    for col in ("x", "y", "z"):
        if col in metadata.columns:
            metadata[col] = (
                metadata[col]
                .astype(str)
                .str.replace(",", ".", regex=False)
                .replace("nan", np.nan)
                .astype(float)
            )

    # Drop helper column
    metadata = metadata.drop(columns=["_norm_label"])

    # Check how many matched
    n_matched = metadata["x"].notna().sum()
    n_total = len(metadata)
    if n_matched == 0:
        raise ValueError(
            f"No coordinates matched for {patient}. "
            f"Check label formats in implant and channel_labels files."
        )
    if n_matched < n_total:
        import warnings
        warnings.warn(
            f"Only {n_matched}/{n_total} channels matched coordinates for {patient}."
        )

    return metadata


def prepare_spatial_coordinates(
    metadata: pd.DataFrame,
    scale: str = "mm",
    center: bool = True,
    to_mni: bool = False,
) -> np.ndarray:
    """Extract and transform electrode coordinates from metadata.

    Parameters
    ----------
    metadata : pd.DataFrame
        Channel metadata with 'x', 'y', 'z' columns (in micrometers).
    scale : str, optional
        Output scale: "mm" (divide by 1000) or "um" (raw micrometers).
        Default is "mm".
    center : bool, optional
        If True, center coordinates around origin. Default is True.
        Ignored if to_mni=True.
    to_mni : bool, optional
        If True, apply estimated transform to MNI-like space for brain
        visualization. This shifts coordinates to match typical brain
        anatomy based on atlas labels. Default is False.

    Returns
    -------
    coords : np.ndarray
        Array of shape (n_channels, 3) with coordinates.

    Raises
    ------
    ValueError
        If required coordinate columns are missing or contain NaN values.
    """
    required_cols = ["x", "y", "z"]
    missing = [c for c in required_cols if c not in metadata.columns]
    if missing:
        raise ValueError(f"Missing coordinate columns: {missing}")

    coords = metadata[required_cols].to_numpy(dtype=float)

    if np.isnan(coords).any():
        n_nan = np.isnan(coords).any(axis=1).sum()
        raise ValueError(
            f"Found {n_nan} channels with missing coordinate values. "
            "Ensure implant metadata is complete."
        )

    if scale == "mm":
        coords = coords / 1000.0
    elif scale != "um":
        raise ValueError(f"Unknown scale '{scale}'. Use 'mm' or 'um'.")

    if to_mni:
        coords = _estimate_mni_transform(coords, metadata)
    elif center:
        coords = coords - coords.mean(axis=0)

    return coords


def _estimate_mni_transform(coords: np.ndarray, metadata: pd.DataFrame) -> np.ndarray:
    """Estimate transform from native space to MNI-like coordinates.

    Uses Desikan-Killiany atlas labels and their known MNI coordinates to
    estimate a translation that aligns native coordinates to MNI space.

    This is an approximation - for accurate MNI coordinates, use proper
    registration with patient MRI data.

    Parameters
    ----------
    coords : np.ndarray
        Coordinates in native space (mm).
    metadata : pd.DataFrame
        Metadata with Desikan-Killany atlas labels.

    Returns
    -------
    mni_coords : np.ndarray
        Coordinates in approximate MNI space.
    """
    # Known MNI coordinates for Desikan-Killiany regions (approximate centroids)
    # Values from standard MNI atlases
    MNI_REGION_COORDS = {
        # Left hemisphere cortical regions
        "ctx-lh-fusiform": (-35, -50, -18),
        "ctx-lh-parahippocampal": (-25, -25, -20),
        "ctx-lh-middletemporal": (-55, -25, -10),
        "ctx-lh-inferiortemporal": (-50, -30, -25),
        "ctx-lh-superiortemporal": (-55, -15, 0),
        "ctx-lh-transversetemporal": (-45, -20, 10),
        "ctx-lh-temporalpole": (-35, 15, -30),
        "ctx-lh-bankssts": (-55, -45, 5),
        "ctx-lh-entorhinal": (-25, -10, -30),
        "ctx-lh-insula": (-38, 0, 5),
        "ctx-lh-lateralorbitofrontal": (-30, 35, -15),
        "ctx-lh-medialorbitofrontal": (-8, 45, -15),
        "ctx-lh-parsorbitalis": (-40, 40, -10),
        "ctx-lh-parstriangularis": (-48, 30, 5),
        "ctx-lh-parsopercularis": (-50, 15, 10),
        "ctx-lh-rostralmiddlefrontal": (-35, 45, 20),
        "ctx-lh-caudalmiddlefrontal": (-40, 15, 45),
        "ctx-lh-superiorfrontal": (-15, 35, 45),
        "ctx-lh-frontalpole": (-10, 65, -5),
        "ctx-lh-precentral": (-40, -10, 55),
        "ctx-lh-postcentral": (-45, -25, 55),
        "ctx-lh-paracentral": (-10, -30, 65),
        "ctx-lh-supramarginal": (-55, -40, 35),
        "ctx-lh-inferiorparietal": (-45, -60, 40),
        "ctx-lh-superiorparietal": (-25, -60, 55),
        "ctx-lh-precuneus": (-10, -60, 40),
        "ctx-lh-cuneus": (-10, -85, 20),
        "ctx-lh-lingual": (-15, -70, -5),
        "ctx-lh-pericalcarine": (-10, -85, 5),
        "ctx-lh-lateraloccipital": (-40, -80, 10),
        "ctx-lh-rostralanteriorcingulate": (-8, 35, 10),
        "ctx-lh-caudalanteriorcingulate": (-8, 15, 30),
        "ctx-lh-posteriorcingulate": (-8, -40, 30),
        "ctx-lh-isthmuscingulate": (-10, -45, 10),
        "ctx-lh-hippocampus": (-28, -20, -15),
        "ctx-lh-amygdala": (-25, -5, -20),
        # Right hemisphere (mirror of left)
        "ctx-rh-fusiform": (35, -50, -18),
        "ctx-rh-parahippocampal": (25, -25, -20),
        "ctx-rh-middletemporal": (55, -25, -10),
        "ctx-rh-inferiortemporal": (50, -30, -25),
        "ctx-rh-superiortemporal": (55, -15, 0),
        "ctx-rh-transversetemporal": (45, -20, 10),
        "ctx-rh-temporalpole": (35, 15, -30),
        "ctx-rh-bankssts": (55, -45, 5),
        "ctx-rh-entorhinal": (25, -10, -30),
        "ctx-rh-insula": (38, 0, 5),
        "ctx-rh-lateralorbitofrontal": (30, 35, -15),
        "ctx-rh-medialorbitofrontal": (8, 45, -15),
        "ctx-rh-parsorbitalis": (40, 40, -10),
        "ctx-rh-parstriangularis": (48, 30, 5),
        "ctx-rh-parsopercularis": (50, 15, 10),
        "ctx-rh-rostralmiddlefrontal": (35, 45, 20),
        "ctx-rh-caudalmiddlefrontal": (40, 15, 45),
        "ctx-rh-superiorfrontal": (15, 35, 45),
        "ctx-rh-frontalpole": (10, 65, -5),
        "ctx-rh-precentral": (40, -10, 55),
        "ctx-rh-postcentral": (45, -25, 55),
        "ctx-rh-paracentral": (10, -30, 65),
        "ctx-rh-supramarginal": (55, -40, 35),
        "ctx-rh-inferiorparietal": (45, -60, 40),
        "ctx-rh-superiorparietal": (25, -60, 55),
        "ctx-rh-precuneus": (10, -60, 40),
        "ctx-rh-rostralanteriorcingulate": (8, 35, 10),
        "ctx-rh-caudalanteriorcingulate": (8, 15, 30),
        "ctx-rh-posteriorcingulate": (8, -40, 30),
        "ctx-rh-isthmuscingulate": (10, -45, 10),
        "ctx-rh-hippocampus": (28, -20, -15),
        "ctx-rh-amygdala": (25, -5, -20),
        # Subcortical
        "Left-Hippocampus": (-28, -20, -15),
        "Right-Hippocampus": (28, -20, -15),
        "Left-Amygdala": (-25, -5, -20),
        "Right-Amygdala": (25, -5, -20),
        "Left-Thalamus": (-12, -18, 8),
        "Right-Thalamus": (12, -18, 8),
        "Left-Caudate": (-12, 12, 10),
        "Right-Caudate": (12, 12, 10),
        "Left-Putamen": (-25, 5, 2),
        "Right-Putamen": (25, 5, 2),
        "Left-Accumbens-area": (-10, 10, -8),
        "Right-Accumbens-area": (10, 10, -8),
        "Left-Pallidum": (-18, 0, 0),
        "Right-Pallidum": (18, 0, 0),
        "Cerebellum-Cortex": (0, -55, -35),
        "Brain-Stem": (0, -30, -30),
    }

    def get_primary_region(dk_str):
        """Extract primary region from Desikan-Killiany string."""
        if pd.isna(dk_str):
            return None
        parts = str(dk_str).strip().split(",")
        if parts:
            return parts[0].strip().strip('"').strip()
        return None

    # Detect coordinate system orientation using atlas labels
    # In MNI: inferior structures (cerebellum, temporal pole) have negative Z
    #         superior structures (superior frontal) have positive Z
    # Some patients have inverted Z (more negative = more inferior)
    # Others have standard Z (more positive = more superior)

    z_inverted = False
    if "Desikan-Killany" in metadata.columns:
        # Find inferior regions (cerebellum, temporal pole, fusiform)
        inferior_z = []
        superior_z = []
        for idx, row in metadata.iterrows():
            region = get_primary_region(row.get("Desikan-Killany"))
            if region:
                region_lower = region.lower()
                if any(r in region_lower for r in ["cerebellum", "temporalpole", "fusiform", "entorhinal"]):
                    inferior_z.append(coords[idx, 2])
                elif any(r in region_lower for r in ["superiorfrontal", "precentral", "postcentral", "superiorparietal"]):
                    superior_z.append(coords[idx, 2])

        if inferior_z and superior_z:
            # If inferior regions have MORE negative Z than superior, it's inverted
            z_inverted = np.mean(inferior_z) < np.mean(superior_z)
        elif inferior_z:
            # If only inferior regions exist and Z is very negative, likely inverted
            z_inverted = np.mean(inferior_z) < -20
        elif superior_z:
            # If only superior regions exist and Z is very positive, likely standard
            z_inverted = np.mean(superior_z) < 0

    # Work with potentially flipped coordinates
    working_coords = coords.copy()

    # Collect matched region coordinates for calibration
    native_coords_matched = []
    mni_coords_matched = []

    if "Desikan-Killany" in metadata.columns:
        for idx, row in metadata.iterrows():
            region = get_primary_region(row.get("Desikan-Killany"))
            if region and region in MNI_REGION_COORDS:
                native_coords_matched.append(working_coords[idx])
                mni_coords_matched.append(MNI_REGION_COORDS[region])

    # Affine least-squares fit when enough anchors exist (per-axis scale + shift
    # corrects for native frames that differ from MNI by more than a translation).
    # Falls back to median translation, then to a centroid-heuristic.
    mni_coords = None
    if len(native_coords_matched) >= 6:
        native_arr = np.array(native_coords_matched, dtype=float)
        mni_arr = np.array(mni_coords_matched, dtype=float)
        # Per-axis diagonal affine: mni_k = s_k * native_k + b_k.
        # Less expressive than a full 3×4 affine (no rotation/shear), but robust
        # because the anchors per region are coincident in MNI — a 12-parameter
        # fit would over-rotate. Diagonal is the right model for sEEG.
        scales = np.empty(3)
        shifts = np.empty(3)
        for k in range(3):
            A = np.column_stack([native_arr[:, k], np.ones(len(native_arr))])
            sol, *_ = np.linalg.lstsq(A, mni_arr[:, k], rcond=None)
            scales[k], shifts[k] = sol
        # Sanity: scale should be ≈ ±1 within ~30%; reject if degenerate.
        if np.all(np.abs(scales) > 0.4) and np.all(np.abs(scales) < 2.5):
            mni_coords = working_coords * scales + shifts

    if mni_coords is None:
        if len(native_coords_matched) >= 3:
            native_arr = np.array(native_coords_matched)
            mni_arr = np.array(mni_coords_matched)
            translation = np.median(mni_arr - native_arr, axis=0)
        else:
            if "Desikan-Killany" in metadata.columns:
                atlas = metadata["Desikan-Killany"].dropna().astype(str)
                lh_count = sum("lh" in a or "Left" in a for a in atlas)
                rh_count = sum("rh" in a or "Right" in a for a in atlas)
                is_bilateral = lh_count > 10 and rh_count > 10
                is_predominantly_left = lh_count > rh_count
            else:
                is_bilateral = False
                is_predominantly_left = working_coords[:, 0].mean() < 0
            centroid = working_coords.mean(axis=0)
            if is_bilateral:
                target_centroid = np.array([0.0, -20.0, 10.0])
            elif is_predominantly_left:
                target_centroid = np.array([-35.0, -20.0, 10.0])
            else:
                target_centroid = np.array([35.0, -20.0, 10.0])
            translation = target_centroid - centroid
        mni_coords = working_coords + translation

    # Clamp to brain-hull MNI bounds (slightly tighter than before).
    mni_coords[:, 0] = np.clip(mni_coords[:, 0], -72, 72)
    mni_coords[:, 1] = np.clip(mni_coords[:, 1], -105, 70)
    mni_coords[:, 2] = np.clip(mni_coords[:, 2], -50, 80)

    # Ellipsoidal cap: pull any contact outside a brain-shaped envelope back to
    # its surface along the radial direction from the brain centre.  Box-clip
    # alone leaves contacts at the cube corners that the MNI152 mesh does not
    # cover; the ellipsoid is the right brain approximation for 3-D glass-brain
    # rendering.
    centre = np.array([0.0, -15.0, 10.0])
    semi = np.array([66.0, 84.0, 56.0])  # M-L, A-P, S-I half-extents
    inset = 0.95  # pull projected contacts slightly inside the surface
    rel = (mni_coords - centre) / semi
    r = np.linalg.norm(rel, axis=1)
    outside = r > 1.0
    if outside.any():
        mni_coords[outside] = (
            centre + (mni_coords[outside] - centre) * (inset / r[outside, None])
        )

    return mni_coords
