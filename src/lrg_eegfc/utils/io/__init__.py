"""Data loading and inspection utilities."""

from lrg_eegfc.config.const import list_patients
from .loaders import load_data_dict, load_mat_pat_data
from .patient import (
    PatientRecording,
    load_dataset,
    load_patient_dataset,
    load_patient_metadata,
    load_timeseries,
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
    "load_dataset",
    "load_patient_dataset",
    "load_patient_metadata",
    "load_timeseries",
    "load_patient_dataset_robust",
    "inspect_mat_file",
    "inspect_patient",
    "inspect_all_patients",
    "generate_report",
    "save_report",
    "generate_csv_rows",
    "save_csv",
    "list_patients",
]
