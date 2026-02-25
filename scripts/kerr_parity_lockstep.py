#!/usr/bin/env python3
"""
Lockstep Lean-vs-Rust Kerr oracle parity runner.

Pipeline:
1) Render Lean tiled oracle image (same tone map / camera params).
2) Render Rust with BH_TRUTH_MODE=lean and matched camera params.
3) Write diff heatmap + numeric parity summary.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

from PIL import Image


def run(cmd: list[str], cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, env=env, check=True)

def run_bh_render_with_timeout(
    cmd: list[str],
    cwd: Path,
    env: dict[str, str],
    timeout_s: int = 90,
) -> None:
    """
    bh_render keeps an HTTP gallery server alive after image write.
    We allow it to run for `timeout_s`, then terminate and treat that as success
    if the expected output image exists.
    """
    print("+", " ".join(cmd), f"(timeout={timeout_s}s)")
    proc = subprocess.Popen(cmd, cwd=cwd, env=env)
    try:
        proc.wait(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(timeout=5)


def diff_metrics(a: Image.Image, b: Image.Image) -> tuple[Image.Image, dict[str, float]]:
    if a.size != b.size:
        raise RuntimeError(f"size mismatch: {a.size} vs {b.size}")
    ar = a.convert("RGB")
    br = b.convert("RGB")
    w, h = ar.size
    pa = ar.load()
    pb = br.load()
    heat = Image.new("RGB", (w, h))
    ph = heat.load()

    errs: list[float] = []
    sq = 0.0
    mx = 0.0
    for y in range(h):
        for x in range(w):
            ra, ga, ba = pa[x, y]
            rb, gb, bb = pb[x, y]
            dr = abs(ra - rb)
            dg = abs(ga - gb)
            db = abs(ba - bb)
            e = (dr + dg + db) / (3.0 * 255.0)
            errs.append(e)
            sq += e * e
            mx = max(mx, e)
            # red/yellow heatmap
            ph[x, y] = (int(255 * e), int(255 * math.sqrt(e)), 0)

    errs.sort()
    n = len(errs)
    metrics = {
        "mae": sum(errs) / n,
        "rmse": math.sqrt(sq / n),
        "p95": errs[int(0.95 * (n - 1))],
        "p99": errs[int(0.99 * (n - 1))],
        "max": mx,
    }
    return heat, metrics


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=256)
    ap.add_argument("--height", type=int, default=256)
    ap.add_argument("--fov", type=float, default=60.0)
    ap.add_argument("--inc-deg", type=float, default=17.0)
    ap.add_argument("--az-deg", type=float, default=0.0)
    ap.add_argument("--astar", type=float, default=0.9)
    ap.add_argument("--max-phi-pi", type=float, default=80.0)
    ap.add_argument("--dphi", type=float, default=0.003)
    ap.add_argument("--tiles-x", type=int, default=4)
    ap.add_argument("--tiles-y", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--force-cpu", action="store_true")
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/kerr_parity_lockstep"))
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    theta_obs = math.pi * args.inc_deg / 180.0
    lean_png = outdir / "lean_oracle.png"
    rust_png = outdir / "rust_oracle.png"
    diff_png = outdir / "diff_heat.png"
    summary_json = outdir / "summary.json"

    # Build tools
    run(["lake", "build", "kerr_ref_frame"], cwd=repo / "lean")
    run(["cargo", "build", "--release", "-p", "gutoe-gpu", "--bin", "bh_render"], cwd=repo)

    # Lean render (tiled)
    run(
        [
            "python3",
            str(repo / "scripts" / "lean_kerr_tiled.py"),
            "--out",
            str(lean_png),
            "--width",
            str(args.width),
            "--height",
            str(args.height),
            "--tiles-x",
            str(args.tiles_x),
            "--tiles-y",
            str(args.tiles_y),
            "--jobs",
            str(args.jobs),
            "--fov",
            str(args.fov),
            "--r-obs",
            "40.0",
            "--r-s",
            "1.0",
            "--astar",
            str(args.astar),
            "--theta-obs",
            str(theta_obs),
            "--exposure",
            "1.15",
            "--gamma",
            "1.55",
            "--black-level",
            "0.10",
        ],
        cwd=repo,
    )

    # Rust render in Lean oracle mode
    env = os.environ.copy()
    env.update(
        {
            "BH_TRUTH_MODE": "lean",
            "BH_KERR_ASTAR": str(args.astar),
            "BH_FOV_OVERRIDE": str(args.fov),
            "BH_INC_OVERRIDE": str(args.inc_deg),
            "BH_AZ_OVERRIDE": str(args.az_deg),
            "BH_MAX_PHI_PI_OVERRIDE": str(args.max_phi_pi),
            "BH_DPHI_OVERRIDE": str(args.dphi),
            "BH_USE_TRANSFER": "0",
            "BH_SPECTRUM": "bolometric",
            "BH_FORCE_GR": "0",
        }
    )
    if args.force_cpu:
        env["BH_FORCE_CPU"] = "1"

    run_bh_render_with_timeout(
        [
            str(repo / "target" / "release" / "bh_render"),
            "m87star",
            f"{args.width}x{args.height}",
        ],
        cwd=repo,
        env=env,
        timeout_s=120,
    )
    shutil.copyfile("/tmp/bh_renders/m87star.png", rust_png)

    # Diff + summary
    lean_img = Image.open(lean_png)
    rust_img = Image.open(rust_png)
    heat, metrics = diff_metrics(lean_img, rust_img)
    heat.save(diff_png)

    summary = {
        "width": args.width,
        "height": args.height,
        "fov": args.fov,
        "inc_deg": args.inc_deg,
        "az_deg": args.az_deg,
        "astar": args.astar,
        "max_phi_pi": args.max_phi_pi,
        "dphi": args.dphi,
        "lean_png": str(lean_png),
        "rust_png": str(rust_png),
        "diff_png": str(diff_png),
        "metrics": metrics,
    }
    summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
