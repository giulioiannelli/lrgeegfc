"""Data loading and inspection utilities."""

from lrg_eegfc.config.const import list_patients
from .loaders import load_data_dict, load_mat_pat_data
from .patient import (
    PatientRecording,
    bipolar_rereference,
    load_channel_labels,
    load_dataset,
    load_epileptic_nodes,
    load_patient_dataset,
    load_patient_metadata,
    load_timeseries,
    parse_seeg_label,
)
from .patient_robust import load_patient_dataset_robust
from .inspect import (
    inspect_mat_file,
    inspect_patient,
    inspect_all_patients,
    generate_report,
    save_report,
    generate_csv_rows,
    save_csv,
)

__all__ = [
    "load_data_dict",
    "load_mat_pat_data",
    "PatientRecording",
    "bipolar_rereference",
    "load_channel_labels",
    "load_dataset",
    "load_epileptic_nodes",
    "load_patient_dataset",
    "load_patient_metadata",
    "load_timeseries",
    "load_patient_dataset_robust",
    "parse_seeg_label",
    "inspect_mat_file",
    "inspect_patient",
    "inspect_all_patients",
    "generate_report",
    "save_report",
    "generate_csv_rows",
    "save_csv",
    "list_patients",
]
