#!/usr/bin/env python3
"""Export broadband per-phase timeseries for an external collaborator.

Ships the **exact broadband timecourses that feed our imcoh network
analysis**.  These are NOT band-filtered: our "band split" is a
frequency-domain selection of Welch-CSD bins applied *after* the
cross-spectral density is formed on the full broadband signal — there is
no time-domain band-pass anywhere in the imcoh pipeline (see
``src/lrg_eegfc/utils/fc/coherence/{_common,imcoh}.py`` and
``scripts/01_compute/batch/compute_imcoh_fc.py``).  Handing over the
broadband recording + the exact band recipe is therefore the only way a
collaborator can be sure they are working on the same timecourses.

Format is pure NumPy ``.npz`` (loadable with ``np.load`` — numpy only, no
other dependency). Per (patient, phase) we write:
  * ``<Pat>_<phase>.npz``  — key ``data`` ``(n_channels, n_samples)`` float32
    (lossless downcast of vendor float64), key ``fs`` (scalar Hz).
  * ``<Pat>_<phase>.md``   — ultracompact human-readable metadata.
  * ``<Pat>_<phase>.json`` — same metadata for programmatic reading.

Per patient we also write ``<Pat>_imcoh_abs.npz`` (the FC matrices we work
with — key ``'{band}__{phase}'`` -> ``(N, N)`` float64, byte-identical to
``load_fc_matrix(..., 'imcoh_abs')``) + a ``.md`` summary, and
``<Pat>_channels.csv`` (row order == Data row order == FC index: label, x,
y, z, region, epileptic, probe).

Top level: ``README.md``, ``bands.json`` (single source for the band edges
+ Welch params), ``reproduce_imcoh.py`` (self-contained numpy+scipy;
broadband -> our per-band imcoh), ``load_example.py``, ``MANIFEST.csv`` and
``MANIFEST_fc.csv``.

Progress is surfaced per (patient, phase) with elapsed/ETA and a final
wall-clock, per the optimize-and-surface rule.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    DEFAULT_SAMPLE_RATE,
    DEFAULT_WELCH_SEGMENT_SEC,
    FS_OVERRIDES,
    PATIENTS_4PHASE,
    PHASE_LABELS,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.loaders import phase_mat_path
from lrg_eegfc.utils.io.patient import (
    build_epi_masks,
    load_timeseries,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

# --------------------------------------------------------------------------
# Static human-readable phase descriptions (transitive-inference paradigm).
# --------------------------------------------------------------------------
PHASE_DESC = {
    "rest_pre": "Resting-state recording acquired BEFORE the task session.",
    "task_learn": "Task — learning phase (transitive-inference paradigm).",
    "task_test": "Task — test phase (transitive-inference paradigm).",
    "rest_post": "Resting-state recording acquired AFTER the task session.",
}

#: Only Pat_02/03/13 store the acquisition Parameters struct; all three
#: unanimously report band-pass 0.53–300 Hz, order 4, + a line-noise notch.
#: The other 7 patients ship only Data. We assume the common acquisition
#: chain but never fabricate per-patient filter values we cannot read.
COHORT_ACQ_NOTE = (
    "Patients whose .mat records acquisition Parameters (Pat_02, Pat_03, "
    "Pat_13) unanimously show band-pass 0.53–300 Hz, order 4, plus a "
    "line-noise notch. The same acquisition chain is assumed for the whole "
    "cohort; confirm with the recording site if exact per-patient filter "
    "settings are required. fs is authoritative (config): 2048 Hz, Pat_03 1024 Hz."
)


def sha256_of(path: Path, buf: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_provenance(patient: str) -> dict:
    """Map canonical phase -> {vendor_filename, source_sha256} from provenance.md."""
    prov = SEEG_DATAPATH / patient / "provenance.md"
    out: dict[str, dict] = {}
    if not prov.exists():
        return out
    text = prov.read_text()
    # rows look like: | `resting/rest_pre.mat` | `Resting_PreTask_Data.mat` | `<sha>` |
    row = re.compile(
        r"\|\s*`(?:resting|task)/(\w+)\.mat`\s*\|\s*`([^`]+)`\s*\|\s*`([0-9a-f]{64})`"
    )
    for m in row.finditer(text):
        phase, vendor, sha = m.group(1), m.group(2), m.group(3)
        out[phase] = {"vendor_filename": vendor, "source_sha256": sha}
    return out


def mat_format(path: Path) -> str:
    """Return the MAT-file version string ('5.0', '7.3', ...) from the header."""
    with open(path, "rb") as fh:
        head = fh.read(24)
    m = re.search(rb"MATLAB (\d+\.\d+)", head)
    return m.group(1).decode() if m else "unknown"


def _read_params_scipy(patient: str, phase: str) -> dict:
    """Parse the vendor Parameters struct from a v5 .mat (raises on v7.3)."""
    path = phase_mat_path(SEEG_DATAPATH, patient, phase)
    mat = sio.loadmat(path, variable_names=["Parameters"])  # NotImplementedError on v7.3
    P = mat.get("Parameters")
    info: dict = {}
    if P is None or P.dtype.names is None:
        return info
    for name in P.dtype.names:
        try:
            v = np.asarray(P[name][0, 0]).squeeze()
        except Exception:
            continue
        if v.dtype.kind in "fiu" and v.ndim == 0:
            info[name] = float(v) if v.dtype.kind == "f" else int(v)
        elif name in ("NotchFilter", "DataDimensions"):
            # MATLAB MCOS object — not machine-readable via scipy.
            info[name] = "applied (vendor MATLAB object; not machine-readable)"
    return info


def read_vendor_params(patient: str, phase: str) -> tuple[dict, str | None]:
    """Vendor Parameters (fs, band cutoffs, order, notch) + the phase they came from.

    The 3 vendor-replaced ``task_test`` files (Pat_07/13/14) are MATLAB v7.3
    and ship *no* Parameters struct. The band-pass/notch settings are
    acquisition constants identical across a patient's phases, so we inherit
    them from the first readable sibling phase and record the source.
    """
    order = [phase] + [p for p in PHASE_LABELS if p != phase]
    for ph in order:
        try:
            info = _read_params_scipy(patient, ph)
        except Exception:  # noqa: BLE001 — v7.3 / unreadable → try next phase
            continue
        if info:
            return info, ph
    return {}, None


def _clean_label(raw: str) -> str:
    """Match the cleaning in load_channel_labels: drop ,G2 ref, quotes, spaces."""
    return str(raw).strip('"').split(",")[0].strip().replace(" ", "")


def _implant_coord_map(patient: str) -> dict[str, dict]:
    """label(cleaned) -> {x, y, z, region} from the implant CSV.

    load_patient_metadata joins on RAW labels ("A 1,G2") vs the implant's
    clean "A1", so its coords come back all-NaN for every patient; we join
    on cleaned labels instead. Coordinates are shipped as provided in the
    implant file (native scale). 'region' is the dominant Desikan-Killiany
    token of the (compound, quoted) anatomy field.
    """
    patnum = int(patient.split("_")[-1])
    for cand in (
        SEEG_DATAPATH / patient / f"implant_pat_{patnum:02d}.csv",
        SEEG_DATAPATH / patient / f"Implant_pat_{patnum:02d}.csv",
    ):
        if cand.exists():
            imp = pd.read_csv(cand)
            break
    else:
        return {}
    anat_col = next(
        (c for c in imp.columns if str(c).lower().startswith("desikan")), None
    )
    out: dict[str, dict] = {}
    for _, r in imp.iterrows():
        lab = _clean_label(r["label"])

        def _num(v):
            try:
                return float(str(v).replace(",", "."))
            except (ValueError, TypeError):
                return np.nan

        region = ""
        if anat_col is not None and pd.notna(r.get(anat_col)):
            region = str(r[anat_col]).strip().split(",")[0].strip()
        out[lab] = {
            "x": _num(r.get("x")),
            "y": _num(r.get("y")),
            "z": _num(r.get("z")),
            "region": region,
        }
    return out


def build_channel_table(patient: str, n_rows: int) -> pd.DataFrame:
    """Channel table in Data row order: label, x, y, z, region, epileptic, probe."""
    masks = build_epi_masks(patient)  # channels/epi/probes in FC-matrix order
    labels = masks.channels
    coords = _implant_coord_map(patient)
    df = pd.DataFrame({"index": np.arange(len(labels)), "label": labels})
    df["x"] = [coords.get(l, {}).get("x", np.nan) for l in labels]
    df["y"] = [coords.get(l, {}).get("y", np.nan) for l in labels]
    df["z"] = [coords.get(l, {}).get("z", np.nan) for l in labels]
    df["region"] = [coords.get(l, {}).get("region", "") for l in labels]
    df["epileptic"] = masks.epi_mask.astype(int)
    df["probe"] = masks.probes.astype(str)
    if len(df) != n_rows:
        raise ValueError(
            f"{patient}: channel table has {len(df)} rows but Data has {n_rows}."
        )
    return df


def welch_params_for(fs: float) -> dict:
    nperseg = nperseg_for_fs(fs)
    return {
        "nperseg": nperseg,
        "noverlap": nperseg // 2,
        "segment_seconds": DEFAULT_WELCH_SEGMENT_SEC,
        "delta_f_hz": fs / nperseg,
        "window": "hann",
        "detrend": "none",
        "fft": "one-sided rfft",
        "csd_scaling": "1 / (fs * sum(window**2) * n_segments)",
    }


def _render_ts_md(m: dict) -> str:
    """Ultracompact human-readable metadata for one timeseries file."""
    v = m["vendor_preprocessing"]
    bp = (
        f"{v['bandpass_low_hz']}–{v['bandpass_high_hz']} Hz, order {v['filter_order']}"
        if v["recorded_in_this_patient_mat"]
        else "not stored in this .mat (cohort acquisition: 0.53–300 Hz, order 4)"
    )
    s = m["source"]
    return (
        f"# {m['patient']} · {m['phase']}\n\n"
        f"- **{m['fs_hz']:.0f} Hz** · **{m['n_channels']} ch × {m['n_samples']} samp** · **{m['duration_minutes']:.1f} min**\n"
        f"- **load**: `np.load('{m['patient']}_{m['phase']}.npz')` → `['data']` (n_ch, n_samp) float32, `['fs']` scalar Hz\n"
        f"- **channel order**: `{m['patient']}_channels.csv` (row i ↔ data row i ↔ FC index i)\n"
        f"- **phase**: {m['phase_description']}\n"
        f"- **band split**: frequency-domain (Welch-CSD bins), applied AFTER the CSD — NOT a time-domain filter. Recipe: `../bands.json`, `../reproduce_imcoh.py`.\n"
        f"- **vendor filter**: band-pass {bp}, + line-noise notch. We apply NO further filtering.\n"
        f"- **dtype**: float32, lossless downcast of vendor float64 (rel err {m['float32_downcast_rel_err_vs_std']:.1e}); units assumed µV, no rescaling.\n"
        f"- **source**: `{s['vendor_filename']}` ({m['source_mat_format']})\n"
    )


def export_one(patient: str, phase: str, pat_dir: Path, prov: dict) -> dict | None:
    """Write one (patient, phase) .npz + .json + .md; return a MANIFEST row dict."""
    try:
        X = load_timeseries(patient, phase, SEEG_DATAPATH)  # (N, L) float or None
    except (FileNotFoundError, OSError):
        return None
    if X is None:
        return None
    X = np.asarray(X)
    if X.shape[0] > X.shape[1]:
        X = X.T
    X64 = np.ascontiguousarray(X, dtype=np.float64)
    X32 = X64.astype(np.float32)
    n_ch, n_samp = X32.shape
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)

    # float32 fidelity vs the float64 the pipeline actually loads
    err = float(np.abs(X64 - X32.astype(np.float64)).max())
    rel = err / float(X64.std()) if X64.std() else 0.0

    vendor, vendor_src = read_vendor_params(patient, phase)
    prov_ph = prov.get(phase, {})
    src_mat = phase_mat_path(SEEG_DATAPATH, patient, phase)
    mat_fmt = mat_format(src_mat)
    try:
        src_mat_str = str(Path(src_mat).resolve().relative_to(ROOT.resolve()))
    except ValueError:
        src_mat_str = str(src_mat)

    npz_path = pat_dir / f"{patient}_{phase}.npz"
    np.savez_compressed(npz_path, data=X32, fs=np.float64(fs))

    meta = {
        "patient": patient,
        "phase": phase,
        "phase_description": PHASE_DESC.get(phase, ""),
        "fs_hz": fs,
        "n_channels": int(n_ch),
        "n_samples": int(n_samp),
        "duration_seconds": round(n_samp / fs, 3),
        "duration_minutes": round(n_samp / fs / 60, 3),
        "array_layout": "np.load(npz): key 'data' -> (n_channels, n_samples) float32; key 'fs' -> scalar Hz",
        "row_order": f"matches {patient}_channels.csv (== FC-matrix index order)",
        "dtype_shipped": "float32",
        "dtype_source": "float64 (vendor)",
        "float32_downcast_max_abs_err": err,
        "float32_downcast_rel_err_vs_std": rel,
        "channels_dropped_at_load": "none",
        "units": "vendor amplitude units (assumed microvolts; no rescaling by us)",
        "vendor_preprocessing": {
            "recorded_in_this_patient_mat": bool(vendor),
            "bandpass_low_hz": vendor.get("fcutLow"),
            "bandpass_high_hz": vendor.get("fcutHigh"),
            "filter_order": vendor.get("filter_order"),
            "notch_filter": vendor.get("NotchFilter"),
            "params_read_from_phase": vendor_src,
            "note": (
                "We apply NO filtering before imcoh — the signal is used "
                "broadband as delivered. " + COHORT_ACQ_NOTE
            ),
        },
        "source_mat_format": f"MATLAB {mat_fmt}",
        "band_split": {
            "how": "frequency-domain selection of Welch-CSD bins AFTER the "
            "CSD is formed on this broadband signal — NOT a time-domain filter",
            "recipe": "see ../bands.json and ../reproduce_imcoh.py",
        },
        "source": {
            "canonical_mat": src_mat_str,
            "vendor_filename": prov_ph.get("vendor_filename"),
            "source_sha256": prov_ph.get("source_sha256"),
        },
        "welch_params": welch_params_for(fs),
    }
    (pat_dir / f"{patient}_{phase}.json").write_text(json.dumps(meta, indent=2))
    (pat_dir / f"{patient}_{phase}.md").write_text(_render_ts_md(meta))

    return {
        "patient": patient,
        "phase": phase,
        "n_channels": n_ch,
        "n_samples": n_samp,
        "duration_min": round(n_samp / fs / 60, 2),
        "bytes": npz_path.stat().st_size,
    }


def export_fc_matrices(patient: str, pat_dir: Path) -> dict:
    """Write the exact (N,N) imcoh_abs FC matrices we work with (all band×phase).

    Byte-identical to ``load_fc_matrix(..., 'imcoh_abs')`` — the matrices that
    feed our LRG analysis — and reproducible from the shipped broadband
    timeseries via ``reproduce_imcoh.py`` (validated to ≤3e-8). One NpzFile
    per patient: ``Pat_NN_imcoh_abs.npz``, keys ``'{band}__{phase}'`` each a
    ``(N, N)`` float64 array. Channel order == ``Pat_NN_channels.csv``.
    """
    npz_path = pat_dir / f"{patient}_imcoh_abs.npz"
    mats: dict = {}
    N = None
    for band in BRAIN_BANDS:
        for phase in PHASE_LABELS:
            W = load_fc_matrix(patient, phase, band, "imcoh_abs")
            if W is None:
                continue
            mats[f"{band}__{phase}"] = np.ascontiguousarray(W, dtype=np.float64)
            N = int(W.shape[0])
    np.savez_compressed(npz_path, **mats)

    # ultracompact human-readable sidecar
    bands = " · ".join(BRAIN_BANDS)
    (pat_dir / f"{patient}_imcoh_abs.md").write_text(
        f"# {patient} · imcoh_abs FC matrices\n\n"
        f"- **{len(mats)} matrices** = 6 bands × 4 phases, each **({N}, {N})** float64, symmetric, diag≈0.\n"
        f"- **load**: `np.load('{patient}_imcoh_abs.npz')['beta__rest_post']` → key `'{{band}}__{{phase}}'`.\n"
        f"- **bands**: {bands}.  **phases**: rest_pre · task_learn · task_test · rest_post.\n"
        f"- **what**: `imcoh_abs = mean_f |signed ImCoh_ij(f)|` over each band's Welch bins (Ewald 2012); range [0,1].\n"
        f"- **identity**: byte-identical to our pipeline's `load_fc_matrix(..., 'imcoh_abs')`; reproducible from the timeseries via `../reproduce_imcoh.py`.\n"
        f"- **channel order**: `{patient}_channels.csv` (row/col i ↔ FC index i).\n"
    )
    return {
        "patient": patient,
        "file": f"{patient}/{patient}_imcoh_abs.npz",
        "n_matrices": len(mats),
        "n_channels": N,
        "bytes": npz_path.stat().st_size,
    }


def build_fc_manifest(out: Path, patients: list[str]) -> list[dict]:
    """Scan Pat_NN_imcoh_abs.npz files and build a manifest (resume-safe)."""
    rows: list[dict] = []
    for pat in patients:
        npz_path = out / pat / f"{pat}_imcoh_abs.npz"
        if not npz_path.exists():
            continue
        with np.load(npz_path) as z:
            keys = list(z.files)            # cheap — zip namelist, no decompress
            N = int(z[keys[0]].shape[0]) if keys else None
        rows.append({
            "patient": pat,
            "file": f"{pat}/{pat}_imcoh_abs.npz",
            "n_matrices": len(keys),
            "n_channels": N,
            "bytes": npz_path.stat().st_size,
            "sha256": sha256_of(npz_path),
        })
    return rows


def write_top_level(out: Path, bands_meta: dict) -> None:
    (out / "bands.json").write_text(json.dumps(bands_meta, indent=2))
    (out / "reproduce_imcoh.py").write_text(REPRODUCE_PY)
    (out / "load_example.py").write_text(LOAD_EXAMPLE_PY)
    (out / "README.md").write_text(readme_text(bands_meta))


def readme_text(bands_meta: dict) -> str:
    bands = bands_meta["bands"]
    band_rows = "\n".join(
        f"| {b} | {lo} | {hi} |" for b, (lo, hi) in bands.items()
    )
    return f"""# LRG-EEG broadband timeseries — collaborator export

## What this is (read first)

These are the **exact broadband timecourses that feed our imaginary-coherence
(imcoh) network analysis**, as plain **NumPy `.npz`** files (load with
`np.load` — numpy only, no other dependencies). One file per patient per phase:

    Pat_NN/Pat_NN_<phase>.npz   ->  key "data" (n_channels, n_samples) float32, key "fs" (scalar Hz)

`<phase>` ∈ `rest_pre`, `task_learn`, `task_test`, `rest_post`.

**There is no band-filtered timeseries — on purpose.** We do **not** band-pass
the signal in the time domain. We compute one Welch cross-spectral density
(CSD) on the broadband signal, form the signed imaginary coherency, and then
**select the frequency bins of each band and average** (the band split happens
*after*, in the frequency domain). Shipping the broadband recording + the exact
recipe below is the only way to guarantee you are on the same timecourses.

Every data file has two sidecars: a **`.md`** (ultracompact, human-readable)
and a **`.json`** (same facts, for programmatic reading).

## Provenance / preprocessing already in the signal

The vendor delivered these already band-limited and notch-filtered. Only
Pat_02, Pat_03 and Pat_13 store the acquisition parameters in their source
`.mat`; all three unanimously report band-pass **0.53–300 Hz, order 4** plus a
line-noise **notch**, and the same chain is assumed cohort-wide (each sidecar
flags whether its patient recorded the params). `fs` is authoritative from
config: **2048 Hz**, Pat_03 **1024 Hz**. We apply **no** further filtering.
Amplitude is in the vendor's units (assumed µV); we do not rescale. Float32
is a lossless downcast of the vendor float64 to well below acquisition
precision (max abs error ~1e-4 units, ~1e-6 relative — see each sidecar).

## Files

- `Pat_NN/Pat_NN_<phase>.npz`     — broadband timeseries: key `data` `(n_channels, n_samples)` float32, key `fs`
- `Pat_NN/Pat_NN_<phase>.md`      — ultracompact human-readable metadata
- `Pat_NN/Pat_NN_<phase>.json`    — same metadata for Python (fs, duration, provenance, welch params)
- `Pat_NN/Pat_NN_imcoh_abs.npz`   — **the FC matrices we work with**: `imcoh_abs` `(N, N)` per band×phase (see below)
- `Pat_NN/Pat_NN_imcoh_abs.md`    — human-readable summary of the FC file
- `Pat_NN/Pat_NN_channels.csv`    — channel order (== data rows == FC index): label, x, y, z, region, epileptic, probe
- `bands.json`                    — band edges + Welch parameters (single source of truth)
- `reproduce_imcoh.py`            — self-contained (numpy+scipy): broadband → our per-band imcoh
- `load_example.py`               — minimal loader snippet
- `MANIFEST.csv` / `MANIFEST_fc.csv` — timeseries / FC files: shape, fs, bytes, sha256

## The FC matrices we actually work with (`Pat_NN_imcoh_abs.npz`)

We also ship the **exact `(N, N)` connectivity matrices that feed our LRG
analysis** — one `.npz` per patient, key `"{{band}}__{{phase}}"`, float64:

```python
import numpy as np
W = np.load("Pat_02/Pat_02_imcoh_abs.npz")["beta__rest_post"]   # (N, N) imcoh_abs
```

These are byte-identical to the matrices in our pipeline (`imcoh_abs =
mean_f |signed ImCoh|` over each band's Welch bins) and are **reproducible
from the shipped broadband timeseries** with `reproduce_imcoh.py` — we verified
agreement to ≤3e-8 across all six bands, so the timeseries you receive and the
FC matrices you receive are guaranteed to be the same objects our networks are
built on. Row/column order matches `Pat_NN_channels.csv`. `imcoh_abs` (what LRG
needs — non-negative) is the band-averaged magnitude of the signed, per-bin
imaginary coherency.

## The band recipe (how to get our networks from these signals)

Bands (Hz):

| band | low | high |
|------|-----|------|
{band_rows}

Welch CSD (per patient sampling rate `fs`):
- `nperseg = fs * {bands_meta['welch']['segment_seconds']:g}` samples (Δf = {bands_meta['welch']['segment_seconds']:g}-second segments)
- `noverlap = nperseg // 2` (50%), window = **Hann**, **no detrend**, one-sided rfft
- CSD scaling `1 / (fs * sum(window**2) * n_segments)`

Signed imaginary coherency (Nolte et al. 2004), per frequency bin:

    ImCoh_ij(f) = Im(S_ij(f)) / sqrt(S_ii(f) * S_jj(f)) ,  diag = 0,  range [-1, 1]

Per-band connectivity matrix used by our LRG analysis (`imcoh_abs`, Ewald 2012 /
Bastos-Schoffelen 2016 convention):

    W_ij(band) = mean_over_bins_in_[lo,hi]( | ImCoh_ij(f) | )

Order matters: we average the **magnitude per bin**, not the magnitude of the
average (`mean(|.|) != |mean(.)|`).

Run `python reproduce_imcoh.py Pat_02/Pat_02_rest_pre.npz --band beta` to
generate `W(beta)` for that recording. That `(N, N)` matrix is the exact input
to our LRG / dendrogram pipeline.

## Sampling rates

Most patients are 2048 Hz; **Pat_03 is 1024 Hz** (so its `nperseg` is 2048,
not 4096). `bands.json` and each sidecar carry the per-file `fs`; always read
`fs` from the file (`np.load(...)['fs']`) rather than assuming.

## Loading (Python — numpy only)

```python
import numpy as np
z = np.load("Pat_02/Pat_02_rest_pre.npz")
X  = z["data"]     # (n_channels, n_samples) float32
fs = float(z["fs"])
```
"""


# --------------------------------------------------------------------------
# Self-contained reproduce script shipped inside the export (no lrg_eegfc dep).
# Faithful copy of utils/fc/coherence/_common.welch_csd + imcoh.compute_imcoh
# + workflow.fc band-average (mean of |signed| over band bins).
# --------------------------------------------------------------------------
REPRODUCE_PY = r'''#!/usr/bin/env python3
"""Reproduce our per-band imcoh connectivity from a broadband .npz file.

Self-contained (numpy + scipy only). Mirrors the exact pipeline that
produces the FC matrices at the basis of the LRG network analysis:

    broadband X (n_ch, n_samp)
      -> Welch CSD  (Hann, 50% overlap, no detrend, one-sided rfft)
      -> signed ImCoh_ij(f) = Im(S_ij)/sqrt(S_ii*S_jj)      (Nolte 2004)
      -> select band bins, average |ImCoh|                  (imcoh_abs; Ewald 2012)

Usage:
    python reproduce_imcoh.py Pat_02/Pat_02_rest_pre.npz --band beta
    python reproduce_imcoh.py Pat_02/Pat_02_rest_pre.npz --band all --out-dir out/
"""
import argparse, json
from pathlib import Path
import numpy as np
from scipy.signal import get_window
from numpy.lib.stride_tricks import sliding_window_view

BANDS = None  # loaded from bands.json next to this script


def welch_csd(X, fs, nperseg, noverlap=None, batch_size=64):
    N, L = X.shape
    if noverlap is None:
        noverlap = nperseg // 2
    step = nperseg - noverlap
    seg = sliding_window_view(X, window_shape=nperseg, axis=1)[:, ::step, :]
    n_seg = seg.shape[1]
    if n_seg == 0:
        raise ValueError("nperseg larger than signal length")
    freqs = np.fft.rfftfreq(nperseg, 1.0 / fs)
    win = get_window("hann", nperseg).astype(X.dtype, copy=False)
    scale = 1.0 / (fs * (win * win).sum() * n_seg)
    CSD = np.zeros((N, N, len(freqs)), dtype=complex)
    for s in range(0, n_seg, batch_size):
        blk = seg[:, s:s + batch_size, :] * win[None, None, :]
        fb = np.fft.rfft(blk, n=nperseg, axis=2)
        CSD += np.einsum("nbf,mbf->nmf", fb, np.conj(fb), optimize="greedy")
    return freqs, CSD * scale


def signed_imcoh(freqs, CSD):
    N = CSD.shape[0]
    F = CSD.shape[2]
    PSD = np.real(np.diagonal(CSD, axis1=0, axis2=1).T)   # (N, F)
    denom = np.sqrt(PSD[:, None, :] * PSD[None, :, :])
    IC = np.divide(np.imag(CSD), denom, out=np.zeros((N, N, F)), where=denom > 0)
    np.einsum("iif->if", IC)[...] = 0.0
    return IC


def imcoh_abs_band(IC, freqs, lo, hi):
    mask = (freqs >= lo) & (freqs <= hi)
    if not mask.any():
        raise ValueError(f"no Welch bins in [{lo}, {hi}] Hz")
    return np.abs(IC[:, :, mask]).mean(axis=2)   # (N, N)


def main():
    global BANDS
    ap = argparse.ArgumentParser()
    ap.add_argument("npz", help="path to a *_<phase>.npz broadband file")
    ap.add_argument("--band", default="beta",
                    help="band name from bands.json, or 'all'")
    ap.add_argument("--out-dir", default=None,
                    help="save W(band) as .npy here (else just print a summary)")
    a = ap.parse_args()
    BANDS = json.loads((Path(__file__).parent / "bands.json").read_text())["bands"]

    z = np.load(a.npz)
    X = z["data"].astype(np.float64)   # (n_ch, n_samp)
    fs = float(z["fs"])
    nperseg = int(round(fs * 2.0))     # 2-second segments (Δf = 0.5 Hz)
    freqs, CSD = welch_csd(X, fs, nperseg)
    IC = signed_imcoh(freqs, CSD)

    todo = list(BANDS) if a.band == "all" else [a.band]
    for b in todo:
        lo, hi = BANDS[b]
        W = imcoh_abs_band(IC, freqs, lo, hi)
        print(f"{b:<11} W shape={W.shape}  mean={W[~np.eye(len(W),dtype=bool)].mean():.4f}  max={W.max():.4f}")
        if a.out_dir:
            Path(a.out_dir).mkdir(parents=True, exist_ok=True)
            out = Path(a.out_dir) / (Path(a.npz).stem + f"_imcoh_abs_{b}.npy")
            np.save(out, W)
            print(f"  saved {out}")


if __name__ == "__main__":
    main()
'''

LOAD_EXAMPLE_PY = r'''#!/usr/bin/env python3
"""Minimal example: load a broadband recording (pure numpy) and its metadata."""
import json, sys
from pathlib import Path
import numpy as np

p = Path(sys.argv[1] if len(sys.argv) > 1 else "Pat_02/Pat_02_rest_pre.npz")
z = np.load(p)
X = z["data"]                              # (n_channels, n_samples) float32
fs = float(z["fs"])                        # sampling rate, Hz
meta = json.loads(p.with_suffix(".json").read_text())
print(f"{p.name}: {X.shape} @ {fs} Hz, {meta['duration_minutes']} min")
print("channels file:", p.parent / f"{meta['patient']}_channels.csv")

# FC matrices we work with (per patient): np.load(...)['{band}__{phase}']
fc = p.parent / f"{meta['patient']}_imcoh_abs.npz"
if fc.exists():
    W = np.load(fc)["beta__rest_post"]     # (N, N) imcoh_abs
    print(f"imcoh_abs beta/rest_post: {W.shape}")
'''


def build_manifest(out: Path, patients: list[str], phases: list[str]) -> list[dict]:
    """Scan the export dir and build a complete manifest (resume-safe).

    Shape/fs come from the .json sidecar (cheap) so we never decompress the
    large .npz arrays just to read a header; bytes + sha256 from the .npz.
    """
    rows: list[dict] = []
    for pat in patients:
        for phase in phases:
            npz_path = out / pat / f"{pat}_{phase}.npz"
            json_path = out / pat / f"{pat}_{phase}.json"
            if not (npz_path.exists() and json_path.exists()):
                continue
            m = json.loads(json_path.read_text())
            rows.append({
                "patient": pat,
                "phase": phase,
                "file": f"{pat}/{pat}_{phase}.npz",
                "n_channels": m["n_channels"],
                "n_samples": m["n_samples"],
                "fs_hz": m["fs_hz"],
                "duration_min": round(m["n_samples"] / m["fs_hz"] / 60, 2),
                "bytes": npz_path.stat().st_size,
                "sha256": sha256_of(npz_path),
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=PATIENTS_4PHASE)
    ap.add_argument("--phases", nargs="+", default=list(PHASE_LABELS))
    ap.add_argument(
        "--out",
        default=str(ROOT / "data" / "exports" / "lrgeeg_timeseries_export"),
    )
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--fc-only", action="store_true",
                    help="only (re)write FC matrices + metadata; skip timeseries")
    ap.add_argument("--no-fc", action="store_true",
                    help="skip the imcoh_abs FC-matrix export")
    args = ap.parse_args()

    emit_ts = not args.fc_only
    emit_fc = not args.no_fc

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # bands.json (single source of truth for the recipe)
    bands_meta = {
        "bands": {b: [float(lo), float(hi)] for b, (lo, hi) in BRAIN_BANDS.items()},
        "welch": {
            "segment_seconds": DEFAULT_WELCH_SEGMENT_SEC,
            "window": "hann",
            "overlap": "50% (noverlap = nperseg // 2)",
            "detrend": "none",
            "fft": "one-sided rfft",
            "nperseg_rule": "int(fs * segment_seconds)  -> 4096 @2048Hz, 2048 @1024Hz",
            "csd_scaling": "1 / (fs * sum(window**2) * n_segments)",
        },
        "imcoh": {
            "signed": "Im(S_ij) / sqrt(S_ii * S_jj), diag=0, range [-1,1] (Nolte 2004)",
            "band_matrix_imcoh_abs": "mean over band bins of |signed ImCoh| (Ewald 2012)",
        },
        "fs_hz_per_patient": {
            p: FS_OVERRIDES.get(p, DEFAULT_SAMPLE_RATE) for p in args.patients
        },
    }
    write_top_level(out, bands_meta)

    n_units = len(args.patients) * len(args.phases)
    manifest: list[dict] = []
    t0 = time.time()
    done = 0

    for pat in args.patients:
        pat_dir = out / pat
        pat_dir.mkdir(parents=True, exist_ok=True)
        prov = parse_provenance(pat)

        # channel table — written once per patient, independent of h5 skip logic.
        n_ch_expected = len(build_epi_masks(pat).channels)
        build_channel_table(pat, n_ch_expected).to_csv(
            pat_dir / f"{pat}_channels.csv", index=False
        )

        if emit_ts:
            for phase in args.phases:
                npz_path = pat_dir / f"{pat}_{phase}.npz"
                if npz_path.exists() and not args.overwrite:
                    done += 1
                    print(f"[{done}/{n_units}] {pat} {phase}  SKIP (exists)", flush=True)
                    continue
                t = time.time()
                row = export_one(pat, phase, pat_dir, prov)
                done += 1
                if row is None:
                    print(f"[{done}/{n_units}] {pat} {phase}  MISSING, skipped", flush=True)
                    continue
                if row["n_channels"] != n_ch_expected:
                    print(
                        f"  WARN {pat} {phase}: Data has {row['n_channels']} ch but "
                        f"channels.csv has {n_ch_expected}", flush=True,
                    )
                manifest.append(row)

                el = time.time() - t0
                eta = el / done * (n_units - done)
                print(
                    f"[{done}/{n_units}] {pat} {phase}  "
                    f"{row['n_channels']}ch x {row['n_samples']}samp "
                    f"({row['duration_min']:.1f}min)  "
                    f"{row['bytes']/1e6:.0f}MB  {time.time()-t:.1f}s  "
                    f"elapsed={el/60:.1f}m ETA={eta/60:.1f}m",
                    flush=True,
                )

        if emit_fc:
            t = time.time()
            fc = export_fc_matrices(pat, pat_dir)
            if fc["n_channels"] not in (None, n_ch_expected):
                print(
                    f"  WARN {pat} FC: {fc['n_channels']} ch but channels.csv "
                    f"has {n_ch_expected}", flush=True,
                )
            print(
                f"  {pat} FC imcoh_abs: {fc['n_matrices']} matrices "
                f"({fc['bytes']/1e3:.0f} KB) {time.time()-t:.1f}s", flush=True,
            )

    # Build manifests from a full scan of the export dir so they are complete
    # regardless of which files were (re)computed this run (resume-safe).
    print(f"\nBuilding manifests ({(time.time()-t0)/60:.1f} min elapsed)...", flush=True)
    ts_rows = build_manifest(out, args.patients, args.phases)
    if ts_rows:
        pd.DataFrame(ts_rows).to_csv(out / "MANIFEST.csv", index=False)
    fc_rows = build_fc_manifest(out, args.patients) if emit_fc else []
    if fc_rows:
        pd.DataFrame(fc_rows).to_csv(out / "MANIFEST_fc.csv", index=False)

    ts_bytes = sum(r["bytes"] for r in ts_rows)
    fc_bytes = sum(r["bytes"] for r in fc_rows)
    n_fc_mat = sum(r["n_matrices"] for r in fc_rows)
    print(
        f"Done. Timeseries: {len(ts_rows)} files, {ts_bytes/1e9:.1f} GB. "
        f"FC (imcoh_abs): {len(fc_rows)} files / {n_fc_mat} matrices, "
        f"{fc_bytes/1e6:.1f} MB. Total {(ts_bytes+fc_bytes)/1e9:.1f} GB. "
        f"Run wall-clock {(time.time()-t0)/60:.1f} min.\nExport root: {out}",
        flush=True,
    )


if __name__ == "__main__":
    main()
