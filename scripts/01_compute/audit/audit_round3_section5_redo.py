#!/usr/bin/env python3
"""Round-3 Section 5 figure redo bundle — orchestrator.

Phase 4-B split (2026-05-29): the 10 redos previously lived in this
single 2236-LOC file. Each redo family is now its own script under
the same folder:

- audit_round3_section5_psi.py         (redo1 + redo2)
- audit_round3_section5_kc.py          (redo3 + redo6..redo10)
- audit_round3_section5_ctm.py         (redo4)
- audit_round3_section5_cross_probe.py (redo5)

This orchestrator runs all four in sequence, preserving the original
bundle behaviour for callers who want everything in one go.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Folder-local imports
sys.path.insert(0, str(Path(__file__).parent))

from audit_round3_section5_psi import main as run_psi
from audit_round3_section5_kc import main as run_kc
from audit_round3_section5_ctm import main as run_ctm
from audit_round3_section5_cross_probe import main as run_cp


def main() -> None:
    print("=== Round-3 Section-5 redo bundle ===")
    run_psi()
    run_kc()
    run_ctm()
    run_cp()
    print("=== all done ===")


if __name__ == "__main__":
    main()
