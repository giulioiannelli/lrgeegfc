"""High-level patient dataset loaders (moved from :mod:`lrg_eegfc.io`)."""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

from ...config.const import PARAMETER_KEYS, PATIENT_CHANNEL_DROP, PHASE_LABELS
from ...config.paths import SEEG_DATAPATH
from .loaders import load_mat_pat_data


__all__ = [
    "PatientRecording",
    "load_timeseries",
    "load_patient_metadata",
    "load_patient_dataset",
    "load_dataset",
    "load_channel_labels",
    "load_epileptic_nodes",
    "bipolar_rereference",
    "parse_seeg_label",
]


def load_channel_labels(
    patient: str,
    root_path: Path | str = SEEG_DATAPATH,
) -> list[str]:
    """Return cleaned channel labels for ``patient`` (matches FC-matrix ordering).

    Reads ``<root>/<patient>/channel_labels.csv``; tolerates an optional
    quoted ``"label"`` header, ``,G2`` reference suffixes, surrounding
    quotes, and intra-label whitespace ("``A 1``" → "``A1``"). Returns
    one entry per recorded channel in the same order used everywhere
    else in the pipeline (FC matrices, LRG eigenvectors, anatomy CSVs).

    Promoted to library 2026-05-08 — there were ≥ 9 private copies
    across scripts/ before this; new callers should import this
    instead of re-implementing.
    """
    csv_path = Path(root_path) / patient / "channel_labels.csv"
    with open(csv_path) as fh:
        first = fh.readline().strip().strip('"').strip("'")
    skip = 1 if first.lower() == "label" else 0
    df = pd.read_csv(csv_path, header=None, skiprows=skip)
    return [
        str(x).strip('"').split(",")[0].strip().replace(" ", "")
        for x in df.iloc[:, 0]
    ]


@dataclass
class PatientRecording:
    """Container holding one patient's time series and optional metadata."""

    timeseries: np.ndarray
    parameters: Mapping[str, float]
    channel_metadata: Optional[pd.DataFrame] = None


def _apply_channel_drop(data: np.ndarray, patient: str, phase: str) -> np.ndarray:
    """Apply patient-specific row-index drops to a (n_channels, n_samples) matrix.

    Per-patient drops are declared in :data:`PATIENT_CHANNEL_DROP`; see the
    docstring there for the rationale (e.g. Pat_10 task phases ship 3 extra
    rows that are absent from the resting recordings — dropping them at
    load time makes every Pat_10 phase 113-channel, enabling cross-phase
    analyses). Drops are applied against the *raw* row indexing, i.e. the
    indices in the vendor's ``Data`` matrix before any reshuffling.
    """
    drops = PATIENT_CHANNEL_DROP.get(patient, {}).get(phase, [])
    if not drops:
        return data
    keep = [i for i in range(data.shape[0]) if i not in set(drops)]
    return data[keep, :]


def _extract_timeseries(mat: Mapping[str, object], patient: str, phase: str) -> np.ndarray:
    if "Data" not in mat:
        raise KeyError(f"'Data' entry not found in MAT file for {patient} {phase}.")

    data = np.asarray(mat["Data"])
    if data.ndim != 2:
        raise ValueError(
            f"Expected 2-D SEEG matrix for {patient} {phase}; received shape {data.shape}."
        )
    if data.shape[0] > data.shape[1]:
        data = data.T
    return _apply_channel_drop(data, patient, phase)


def load_timeseries(patient: str, phase: str, root_path: Path) -> np.ndarray:
    mat = load_mat_pat_data(patient, phase, root_path)
    return _extract_timeseries(mat, patient, phase)


def _load_parameters(mat: Mapping[str, object]) -> Mapping[str, float]:
    parameters: Dict[str, float] = {}
    params = mat.get("Parameters")
    if params is None:
        return parameters

    def _unwrap_scalar(value: object) -> object:
        while isinstance(value, np.ndarray) and value.size == 1:
            value = value.item()
        return value

    names = getattr(getattr(params, "dtype", None), "names", None) or ()
    for name in names:
        if name not in PARAMETER_KEYS:
            continue
        try:
            raw = _unwrap_scalar(params[name])
            parameters[name] = float(raw)
        except Exception:  # noqa: BLE001
            logging.debug("Could not parse parameter %s", name, exc_info=True)
    return parameters


def load_patient_metadata(patient: str, root_path: Path) -> Optional[pd.DataFrame]:
    patient_path = root_path / patient
    try:
        patnum = int(patient.split("_")[-1])
    except ValueError:
        patnum = 0

    # Canonical lowercase first; legacy mixed-case kept as fallback.
    for cand in (
        patient_path / f"implant_pat_{patnum:02d}.csv",
        patient_path / f"Implant_pat_{patnum:02d}.csv",
    ):
        if cand.exists():
            implant_csv = cand
            break
    else:
        return None
    channel_labels_csv = patient_path / "channel_labels.csv"
    if not channel_labels_csv.exists():
        return None

    implant = pd.read_csv(implant_csv)
    labels = pd.read_csv(channel_labels_csv)

    # Apply per-patient label drops so metadata rows match the canonical
    # post-drop channel count used everywhere else in the pipeline.
    label_drops = PATIENT_CHANNEL_DROP.get(patient, {}).get("__labels__", [])
    if label_drops:
        labels = labels.drop(index=list(label_drops)).reset_index(drop=True)

    metadata = labels.merge(implant, on="label", how="left")

    for column in ("x", "y", "z"):
        if column in metadata.columns:
            metadata[column] = (
                metadata[column]
                .astype(str)
                .str.replace(",", ".", regex=False)
                .astype(float)
            )
    return metadata


def load_patient_dataset(
    patient: str,
    root_path: Path,
    phases: Iterable[str] = PHASE_LABELS,
) -> Dict[str, PatientRecording]:
    dataset: Dict[str, PatientRecording] = {}
    metadata = load_patient_metadata(patient, root_path)

    for phase in phases:
        mat = load_mat_pat_data(patient, phase, root_path)
        dataset[phase] = PatientRecording(
            timeseries=_extract_timeseries(mat, patient, phase),
            parameters=_load_parameters(mat),
            channel_metadata=metadata,
        )
    return dataset


def load_dataset(
    root_path: Path,
    patients: Iterable[str],
    phases: Iterable[str] = PHASE_LABELS,
) -> Dict[str, Dict[str, PatientRecording]]:
    return {patient: load_patient_dataset(patient, root_path, phases) for patient in patients}


def parse_seeg_label(label: str) -> Tuple[Optional[str], Optional[int]]:
    """Parse an sEEG channel label into (probe, contact_number).

    Handles the varied formatting found across patients:
    - Spaces between probe and number: ``'A 1'``, ``'A1'``, ``'H  3'``
    - Reference suffix: ``'A 1,G2'``
    - Primed (contralateral) probes: ``"G' 1"``, ``"H' 2"``
    - Typos: ``'Pì 12'`` (ì -> ``'``)
    - Lowercase scalp channels: ``'f3'``, ``'c4'``, ``'o2'``

    Parameters
    ----------
    label : str
        Raw channel label string.

    Returns
    -------
    probe : str or None
        Probe identifier (e.g. ``'A'``, ``"G'"``) or *None* if unparseable.
    contact : int or None
        Contact number or *None* if unparseable.
    """
    # Strip whitespace and surrounding quotes
    s = label.strip().strip('"')
    # Remove ',G2' reference suffix
    s = re.sub(r",G2$", "", s).strip()
    # Normalise the ì typo for prime to a proper apostrophe
    s = s.replace("\u00ec", "'")
    # Collapse spaces around the prime so "X '4" -> "X'4" and "G' 1" -> "G'1"
    s = re.sub(r"\s*'\s*", "'", s)
    # Match probe (letters + optional prime) then optional whitespace then digits
    m = re.match(r"([A-Za-z]+[']?)\s*(\d+)$", s)
    if m is None:
        return None, None
    probe = m.group(1)
    contact = int(m.group(2))
    return probe, contact


def bipolar_rereference(
    timeseries: np.ndarray,
    channel_labels: List[str],
    coordinates: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, List[str], Optional[np.ndarray]]:
    """Apply bipolar re-referencing to sEEG timeseries.

    For each probe, subtracts consecutive contacts (sorted by contact
    number).  Only pairs whose contact numbers differ by exactly 1 are
    used -- gaps are skipped so that distant tissue sites are never
    mixed.

    The bipolar signal for contact pair *(i, i+1)* is defined as::

        bipolar_k = V_{i+1} - V_i

    This cancels the common reference (typically G2) and shared far-field
    activity, leaving only the local voltage gradient between adjacent
    tissue sites.

    Parameters
    ----------
    timeseries : np.ndarray
        Shape ``(n_channels, n_samples)``.
    channel_labels : list of str
        Channel labels in the same order as rows of *timeseries*.
        Labels are parsed with :func:`parse_seeg_label`.
    coordinates : np.ndarray, optional
        Shape ``(n_channels, 3)`` MNI coordinates.  When provided, the
        returned coordinates are midpoints of each consecutive pair.

    Returns
    -------
    bipolar_ts : np.ndarray
        Shape ``(n_bipolar, n_samples)`` where *n_bipolar* is the number
        of valid consecutive pairs found across all probes.
    bipolar_labels : list of str
        Labels of the form ``'{high}-{low}'``, e.g. ``'A2-A1'``.
    bipolar_coords : np.ndarray or None
        Shape ``(n_bipolar, 3)`` midpoint coordinates, or *None* if
        *coordinates* was not provided.

    Raises
    ------
    ValueError
        If *timeseries* and *channel_labels* have incompatible lengths,
        or if no bipolar pairs can be formed.

    Notes
    -----
    Channels whose labels cannot be parsed (e.g. scalp EEG contacts
    ``'f3'``, ``'c4'``) are silently skipped when they form single-
    contact "probes" with no consecutive partner.  If a probe has only
    one contact, it contributes no bipolar channels.
    """
    n_channels, n_samples = timeseries.shape
    if n_channels != len(channel_labels):
        raise ValueError(
            f"timeseries has {n_channels} channels but "
            f"{len(channel_labels)} labels were provided."
        )

    # ---- 1. Parse labels and group by probe ----
    # Map (probe, contact_number) -> row index in timeseries
    probe_contacts: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    n_unparsed = 0
    for idx, raw_label in enumerate(channel_labels):
        probe, contact = parse_seeg_label(raw_label)
        if probe is None:
            n_unparsed += 1
            logging.debug("bipolar_rereference: could not parse label %r", raw_label)
            continue
        probe_contacts[probe].append((contact, idx))

    if n_unparsed > 0:
        logging.info(
            "bipolar_rereference: %d/%d labels could not be parsed and were skipped.",
            n_unparsed,
            n_channels,
        )

    # ---- 2. Build bipolar pairs ----
    bipolar_rows: List[np.ndarray] = []
    bipolar_labels: List[str] = []
    bipolar_coords_list: List[np.ndarray] = []

    for probe_name in sorted(probe_contacts.keys()):
        contacts = probe_contacts[probe_name]
        # Sort by contact number
        contacts.sort(key=lambda x: x[0])

        for k in range(len(contacts) - 1):
            num_lo, idx_lo = contacts[k]
            num_hi, idx_hi = contacts[k + 1]

            # Only consecutive contact numbers
            if num_hi - num_lo != 1:
                continue

            # bipolar = V_{hi} - V_{lo}
            bipolar_signal = timeseries[idx_hi] - timeseries[idx_lo]
            bipolar_rows.append(bipolar_signal)

            # Label: "{probe}{hi}-{probe}{lo}"
            bipolar_labels.append(f"{probe_name}{num_hi}-{probe_name}{num_lo}")

            # Midpoint coordinates
            if coordinates is not None:
                bipolar_coords_list.append(
                    (coordinates[idx_hi] + coordinates[idx_lo]) / 2.0
                )

    if len(bipolar_rows) == 0:
        raise ValueError(
            "No consecutive bipolar pairs found. "
            "Check that channel_labels contain valid sEEG contact labels."
        )

    bipolar_ts = np.stack(bipolar_rows, axis=0)
    bipolar_coords: Optional[np.ndarray] = None
    if coordinates is not None and bipolar_coords_list:
        bipolar_coords = np.stack(bipolar_coords_list, axis=0)

    logging.info(
        "bipolar_rereference: %d monopolar -> %d bipolar channels "
        "(%d probes, %d unparsed labels skipped).",
        n_channels,
        len(bipolar_labels),
        len(probe_contacts),
        n_unparsed,
    )

    return bipolar_ts, bipolar_labels, bipolar_coords


def _find_implant_xlsx(patient: str, patient_dir: Path) -> Optional[Path]:
    """Locate the implant Excel file for a patient.

    The canonical location (post-``lrg-eegfc data normalize``) is
    ``<patient_dir>/implant/implant_pat_{NN}.xlsx``. Legacy paths
    are still searched so code remains functional during transitions.
    See ``.agents/guides/03_implementation/DATA_LAYOUT.md``.
    """
    try:
        patnum = int(patient.split("_")[-1])
    except ValueError:
        patnum = 0

    candidates = [
        # Canonical (post-normalizer)
        patient_dir / "implant" / f"implant_pat_{patnum:02d}.xlsx",
        # Legacy paths, kept for backward compatibility during the transition
        patient_dir / f"Implant_pat_{patnum}.xlsx",
        patient_dir / f"Implant_pat_{patnum:02d}.xlsx",
        patient_dir / "Implant_locations" / f"Implant_pat_{patnum}.xlsx",
        patient_dir / "Implant_locations" / f"Implant_pat_{patnum:02d}.xlsx",
    ]
    return next((c for c in candidates if c.exists()), None)


def load_epileptic_nodes(
    patient: str,
    root_path: Path | str = SEEG_DATAPATH,
) -> list[str]:
    """Extract epileptic node labels (red-coloured) from implant Excel file.

    Parameters
    ----------
    patient : str
        Patient identifier (e.g., ``'Pat_02'``).
    root_path : Path or str
        Root directory containing patient folders.

    Returns
    -------
    list of str
        Electrode labels marked as epileptic (red font) in the implant file.
        Empty list if no implant file is found or *openpyxl* is not installed.
    """
    try:
        import openpyxl  # lazy import — optional dependency
    except ImportError:
        logging.warning("openpyxl is required to read epileptic nodes from Excel")
        return []

    patient_dir = Path(root_path) / patient
    xlsx = _find_implant_xlsx(patient, patient_dir)
    if xlsx is None:
        return []

    wb = openpyxl.load_workbook(xlsx, data_only=True)

    # Prefer a sheet named 'all_leads' / 'ALL_LEADS'
    ws = wb.active
    for name in ("all_leads", "ALL_LEADS"):
        if name in wb.sheetnames:
            ws = wb[name]
            break

    red_labels: list[str] = []
    for row in range(2, ws.max_row + 1):
        cell = ws.cell(row, 1)
        if cell.value is None:
            continue
        font = cell.font
        if font.color is None or font.color.type != "rgb":
            continue
        try:
            rgb = font.color.rgb
        except (TypeError, ValueError):
            continue
        if isinstance(rgb, str) and len(rgb) == 8:
            r, g, b = int(rgb[2:4], 16), int(rgb[4:6], 16), int(rgb[6:8], 16)
            if r > 200 and g < 50 and b < 50:
                red_labels.append(str(cell.value).strip())

    return red_labels
