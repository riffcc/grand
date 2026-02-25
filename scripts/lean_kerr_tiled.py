#!/usr/bin/env python3
"""
Process-level tiled renderer for Lean Kerr reference frames.

Runs `kerr_ref_frame` on independent tiles in parallel, then stitches PNG.
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from PIL import Image


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def tile_size(total: int, n: int, i: int) -> tuple[int, int]:
    base = total // n
    rem = total % n
    w = base + (1 if i < rem else 0)
    x0 = i * base + min(i, rem)
    return x0, w


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True, help="Final output PNG path")
    ap.add_argument("--width", type=int, default=768)
    ap.add_argument("--height", type=int, default=768)
    ap.add_argument("--tiles-x", type=int, default=4)
    ap.add_argument("--tiles-y", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--fov", type=float, default=60.0)
    ap.add_argument("--r-obs", type=float, default=40.0)
    ap.add_argument("--r-s", type=float, default=1.0)
    ap.add_argument("--astar", type=float, default=0.9)
    ap.add_argument("--theta-obs", type=float, default=0.296705972839036)
    ap.add_argument("--exposure", type=float, default=1.15)
    ap.add_argument("--gamma", type=float, default=1.55)
    ap.add_argument("--black-level", type=float, default=0.10)
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    exe = repo / "lean" / ".lake" / "build" / "bin" / "kerr_ref_frame"
    out = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    tiles_dir = out.parent / (out.stem + "_tiles")
    tiles_dir.mkdir(parents=True, exist_ok=True)

    jobs: list[tuple[int, int, int, int, Path]] = []
    for ty in range(args.tiles_y):
        y0, h = tile_size(args.height, args.tiles_y, ty)
        for tx in range(args.tiles_x):
            x0, w = tile_size(args.width, args.tiles_x, tx)
            path = tiles_dir / f"tile_{tx}_{ty}.ppm"
            jobs.append((x0, y0, w, h, path))

    def render_tile(job: tuple[int, int, int, int, Path]) -> tuple[int, int, Path]:
        x0, y0, w, h, path = job
        cmd = [
            str(exe),
            str(path),
            str(w),
            str(h),
            f"{args.fov}",
            f"{args.r_obs}",
            f"{args.r_s}",
            f"{args.astar}",
            f"{args.theta_obs}",
            str(args.width),
            str(args.height),
            str(x0),
            str(y0),
            f"{args.exposure}",
            f"{args.gamma}",
            f"{args.black_level}",
        ]
        run(cmd)
        return x0, y0, path

    print(f"Rendering {len(jobs)} Lean tiles with jobs={args.jobs} ...")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = [pool.submit(render_tile, job) for job in jobs]
        for fut in as_completed(futs):
            x0, y0, p = fut.result()
            print(f"  done tile x={x0} y={y0}: {p.name}")

    canvas = Image.new("RGB", (args.width, args.height))
    for x0, y0, path in sorted(((j[0], j[1], j[4]) for j in jobs), key=lambda t: (t[1], t[0])):
        tile_img = Image.open(path).convert("RGB")
        canvas.paste(tile_img, (x0, y0))
    canvas.save(out)
    print(f"wrote stitched image: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
