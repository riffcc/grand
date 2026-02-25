#!/usr/bin/env python3
"""
Generate Lean-vs-Rust Kerr parity artifacts:
  - lean_reference.csv / lean_reference.png
  - rust_render.png / rust_luma.png
  - abs_diff.png (heatmap)

This is a visual parity scaffold (shape-level), not a claim of exact physical
equivalence between the current Lean reference proxy and full bh_render output.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from pathlib import Path

from PIL import Image


def run(cmd: list[str], cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def load_csv_grid(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    with path.open("r", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            rows.append([float(x) for x in row])
    if not rows:
        raise RuntimeError(f"empty CSV: {path}")
    w = len(rows[0])
    if any(len(r) != w for r in rows):
        raise RuntimeError(f"ragged CSV: {path}")
    return rows


def save_norm_gray(grid: list[list[float]], out_path: Path) -> Image.Image:
    h = len(grid)
    w = len(grid[0])
    vals = [v for row in grid for v in row]
    vmin = min(vals)
    vmax = max(vals)
    scale = 0.0 if abs(vmax - vmin) < 1e-12 else 255.0 / (vmax - vmin)
    img = Image.new("L", (w, h))
    px = img.load()
    for y, row in enumerate(grid):
        for x, v in enumerate(row):
            u = 0 if scale == 0.0 else int(round((v - vmin) * scale))
            px[x, y] = max(0, min(255, u))
    img.save(out_path)
    return img


def save_luma_norm(img_rgb: Image.Image, out_path: Path) -> Image.Image:
    img = img_rgb.convert("RGB")
    w, h = img.size
    src = img.load()
    vals: list[float] = []
    for y in range(h):
        for x in range(w):
            r, g, b = src[x, y]
            vals.append(0.2126 * r + 0.7152 * g + 0.0722 * b)
    vmin = min(vals)
    vmax = max(vals)
    scale = 0.0 if abs(vmax - vmin) < 1e-12 else 255.0 / (vmax - vmin)
    out = Image.new("L", (w, h))
    dst = out.load()
    i = 0
    for y in range(h):
        for x in range(w):
            v = vals[i]
            i += 1
            u = 0 if scale == 0.0 else int(round((v - vmin) * scale))
            dst[x, y] = max(0, min(255, u))
    out.save(out_path)
    return out


def heat_color(u: float) -> tuple[int, int, int]:
    # blue -> cyan -> yellow -> red
    u = max(0.0, min(1.0, u))
    r = max(0, min(255, int(255 * max(0.0, min(1.0, 1.5 - abs(4 * u - 3))))))
    g = max(0, min(255, int(255 * max(0.0, min(1.0, 1.5 - abs(4 * u - 2))))))
    b = max(0, min(255, int(255 * max(0.0, min(1.0, 1.5 - abs(4 * u - 1))))))
    return r, g, b


def save_abs_diff_heatmap(a: Image.Image, b: Image.Image, out_path: Path) -> None:
    if a.size != b.size:
        raise RuntimeError(f"size mismatch: {a.size} vs {b.size}")
    w, h = a.size
    pa = a.convert("L").load()
    pb = b.convert("L").load()
    out = Image.new("RGB", (w, h))
    dst = out.load()
    for y in range(h):
        for x in range(w):
            d = abs(pa[x, y] - pb[x, y]) / 255.0
            dst[x, y] = heat_color(d)
    out.save(out_path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=128)
    ap.add_argument("--height", type=int, default=128)
    ap.add_argument("--astar", type=float, default=0.9)
    ap.add_argument("--fov", type=float, default=14.0)
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/kerr_parity"))
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    lean_dir = repo / "lean"
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    lean_csv = outdir / "lean_reference.csv"
    lean_png = outdir / "lean_reference.png"
    rust_png = outdir / "rust_render.png"
    rust_luma_png = outdir / "rust_luma.png"
    diff_png = outdir / "abs_diff.png"

    run(["lake", "build", "kerr_ref_frame"], cwd=lean_dir)
    run(
        [
            str(lean_dir / ".lake" / "build" / "bin" / "kerr_ref_frame"),
            str(lean_csv),
            str(args.width),
            str(args.height),
        ],
        cwd=lean_dir,
    )

    grid = load_csv_grid(lean_csv)
    lean_norm = save_norm_gray(grid, lean_png)

    run(["cargo", "build", "--release", "-p", "gutoe-gpu", "--bin", "bh_render"], cwd=repo)
    env = os.environ.copy()
    env.update(
        {
            "BH_FORCE_CPU": "1",
            "BH_KERR_ASTAR": str(args.astar),
            "BH_DISK_MODEL": "riaf",
            "BH_PLASMA_MODEL": "grmhd",
            "BH_USE_TRANSFER": "1",
            "BH_TAU_SCALE": "0.08",
            "BH_SPECTRUM": "millimeter",
            "BH_FOV_OVERRIDE": str(args.fov),
        }
    )
    run(
        [
            str(repo / "target" / "release" / "bh_render"),
            "m87star",
            f"{args.width}x{args.height}",
        ],
        cwd=repo,
        env=env,
    )
    shutil.copyfile("/tmp/bh_renders/m87star.png", rust_png)
    rust_luma = save_luma_norm(Image.open(rust_png), rust_luma_png)
    save_abs_diff_heatmap(lean_norm, rust_luma, diff_png)

    print(f"wrote: {lean_csv}")
    print(f"wrote: {lean_png}")
    print(f"wrote: {rust_png}")
    print(f"wrote: {rust_luma_png}")
    print(f"wrote: {diff_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

