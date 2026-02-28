#!/usr/bin/env python3
"""
Rust-native wrapper for first-second multispectral movie rendering.
"""

from __future__ import annotations

import argparse
import subprocess
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/tmp/bh_renders/first_second_multispectral")
    ap.add_argument("--work-dir", default="/tmp")
    ap.add_argument("--clip-seconds", type=float, default=120.0)
    ap.add_argument("--first-seconds", type=float, default=1.0)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--size", type=int, default=420)
    ap.add_argument("--skip-gif", action="store_true")
    args = ap.parse_args()

    cmd = [
        "cargo",
        "run",
        "-q",
        "-p",
        "gutoe-physics",
        "--bin",
        "universe_first_second_multispectral_movie",
        "--",
        "--out-dir",
        args.out_dir,
        "--work-dir",
        args.work_dir,
        "--clip-seconds",
        str(args.clip_seconds),
        "--first-seconds",
        str(args.first_seconds),
        "--fps",
        str(args.fps),
        "--size",
        str(args.size),
    ]
    if args.skip_gif:
        cmd.append("--skip-gif")

    return subprocess.run(cmd, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
