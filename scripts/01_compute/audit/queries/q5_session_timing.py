"""Q5: session timing per patient × phase."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    PATIENTS_4PHASE, PHASE_LABELS, FS_OVERRIDES, DEFAULT_SAMPLE_RATE
)
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries

rows = []
for p in PATIENTS_4PHASE:
    fs = FS_OVERRIDES.get(p, DEFAULT_SAMPLE_RATE)
    durs = {}
    for ph in PHASE_LABELS:
        try:
            ts = load_timeseries(p, ph, SEEG_DATAPATH)
            n_samples = ts.shape[1]
            dur_s = n_samples / fs
            durs[ph] = dur_s
            durs[f"n_{ph}"] = n_samples
            durs[f"shape_{ph}"] = ts.shape
        except Exception as e:
            durs[ph] = float("nan")
            durs[f"n_{ph}"] = -1
            print(f"[{p}/{ph}] FAIL: {e}", file=sys.stderr)
    inter_rest_gap = durs.get("task_learn", 0) + durs.get("task_test", 0)
    rows.append({
        "patient": p,
        "fs_Hz": fs,
        "rest_pre_s": durs.get("rest_pre", float("nan")),
        "task_learn_s": durs.get("task_learn", float("nan")),
        "task_test_s": durs.get("task_test", float("nan")),
        "rest_post_s": durs.get("rest_post", float("nan")),
        "inter_rest_gap_s": inter_rest_gap,
        "rest_pre_n_samples": durs.get("n_rest_pre", -1),
        "rest_post_n_samples": durs.get("n_rest_post", -1),
        "task_test_n_samples": durs.get("n_task_test", -1),
    })

df = pd.DataFrame(rows)
out = Path("data/audit/raw_fc_phase_distance/session_timing.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
print(df.to_string(index=False))
print(f"\nwrote {out}")
