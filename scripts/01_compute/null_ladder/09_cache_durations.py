#!/usr/bin/env python3
"""Cache per-(patient, phase) sample counts so the calibration workers do not
reload multi-GB recordings just to read a shape.

The .info.txt sidecars carry the sampling rate but not the sample count, and the
five-phase sham construction needs all four real durations to size its pseudo-
phases in proportion. Loading four recordings per cell purely for their shapes
dominated the calibration runtime (ETA 9.5 h); this reduces it to one cheap pass.

Writes data/paper_final/w0b_nulls/durations.json
"""
from __future__ import annotations
import json
import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.config.const import FS_OVERRIDES, DEFAULT_SAMPLE_RATE
from lrg_eegfc.utils.io.patient import load_timeseries

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
OUT = ROOT / "data" / "paper_final" / "w0b_nulls" / "durations.json"


def main():
    out = {}
    for pat in COHORT:
        fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
        rec = {"fs": float(fs)}
        for ph in PHASES:
            X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
            if X.shape[0] > X.shape[1]:
                X = X.T
            rec[ph] = int(X.shape[1])
            rec["n_ch"] = int(X.shape[0])
            del X
        out[pat] = rec
        print(f"{pat}: n_ch={rec['n_ch']} " +
              " ".join(f"{p}={rec[p]/fs:.0f}s" for p in PHASES), flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
