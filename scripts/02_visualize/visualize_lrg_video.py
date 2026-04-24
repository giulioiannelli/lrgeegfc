#!/usr/bin/env python3
"""Generate ultrametric threshold frames for LRG tau values."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, LRG_CACHE
from lrg_eegfc.workflow.lrg import load_lrg_result


def _resolve_list(value: List[str] | None, default: List[str]) -> List[str]:
    return value if value else default


def _save_frame(
    matrix: np.ndarray,
    tau: float,
    output_path: Path,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(matrix, cmap="viridis", vmin=0, vmax=1, interpolation="none")
    ax.set_title(f"{title}\n tau={tau:.4e}", fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _maybe_make_gif(frame_paths: List[Path], gif_path: Path, fps: int) -> None:
    try:
        import imageio.v2 as imageio
    except Exception:
        return
    images = [imageio.imread(path) for path in frame_paths]
    gif_path.parent.mkdir(parents=True, exist_ok=True)
    imageio.mimsave(gif_path, images, fps=fps)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create ultrametric threshold frames across tau."
    )
    parser.add_argument("--patient", required=True)
    parser.add_argument("--band", required=True, choices=list(BRAIN_BANDS.keys()))
    parser.add_argument("--phase", required=True, choices=list(PHASE_LABELS))
    parser.add_argument("--fc-method", required=True, choices=["msc", "corr"])
    parser.add_argument("--cache-root", type=Path, default=LRG_CACHE)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--make-gif", action="store_true")
    parser.add_argument("--fps", type=int, default=6)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    result = load_lrg_result(
        args.patient,
        args.phase,
        args.band,
        args.fc_method,
        cache_root=args.cache_root,
    )
    if result is None:
        print("Missing LRG cache.")
        return 1

    ultrametric_square = squareform(result.ultrametric_matrix)
    tau_values = result.entropy_tau

    output_dir = (
        FIGURES_ROOT / "lrg_video"
        / args.patient
        / f"{args.band}_{args.phase}_{args.fc_method}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    frame_paths: List[Path] = []
    total = 0
    for idx, tau in enumerate(tau_values[:: args.stride]):
        if args.max_frames is not None and total >= args.max_frames:
            break
        mask = (ultrametric_square <= tau).astype(float)
        output_path = output_dir / f"frame_{idx:04d}.png"
        _save_frame(mask, float(tau), output_path, f"{args.patient} {args.band} {args.phase}")
        frame_paths.append(output_path)
        total += 1

    if args.make_gif:
        gif_path = output_dir / "lrg_ultrametric.gif"
        _maybe_make_gif(frame_paths, gif_path, args.fps)
        if args.verbose and gif_path.exists():
            print(f"Saved gif: {gif_path}")

    if args.verbose:
        print(f"Frames: {len(frame_paths)} in {output_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
