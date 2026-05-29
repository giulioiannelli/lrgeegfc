#!/usr/bin/env python3
"""Polished bundle figures for the writing-agent handoff -- orchestrator.

Phase 4-B split (2026-05-29): the 9 bundle figures previously lived in
this single 1317-LOC file. Each figure is now its own script under the
same folder:

- q_fig1_trace_scatter.py
- q_fig2_phase_geometry.py
- q_fig3_dS_dP_convergence.py
- q_fig4_drift.py
- q_fig5_distance_class.py
- q_fig6_chord_arc.py
- q_figS1_z_inflation.py
- q_figS2_swarm.py
- q_figS3_stoplight.py

This orchestrator runs all 9 in sequence, preserving the original
bundle behaviour for callers who want everything in one go. Pass
`imcoh_abs` (default) or `imcoh_sq` as the first CLI arg to pick
the substrate.

Also replaces the previous `sys.path.insert(scripts/archive/...)`
injection -- network layout/drawing helpers now import directly from
`lrg_eegfc.visuals.network_layouts` / `network_drawing` (promoted
2026-05-29 in Phase 4-B split 1/7).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from q_fig1_trace_scatter import main as run_fig1
from q_fig2_phase_geometry import main as run_fig2
from q_fig3_dS_dP_convergence import main as run_fig3
from q_fig4_drift import main as run_fig4
from q_fig5_distance_class import main as run_fig5
from q_fig6_chord_arc import main as run_fig6
from q_figS1_z_inflation import main as run_figS1
from q_figS2_swarm import main as run_figS2
from q_figS3_stoplight import main as run_figS3


def main() -> None:
    print("=== q_bundle figure bundle ===")
    run_fig1()
    run_fig2()
    run_fig3()
    run_fig4()
    run_fig5()
    run_fig6()
    run_figS1()
    run_figS2()
    run_figS3()
    print("=== q_bundle done ===")


if __name__ == "__main__":
    main()
