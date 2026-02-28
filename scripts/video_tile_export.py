#!/usr/bin/env python3
"""
Export AI-friendly tiled/frame bundle from a video.

Outputs:
- tile_contact_sheet.png
- frames/frame_XXXX.png
- tile_manifest.json
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import shutil
import subprocess


def run(cmd: list[str]) -> str:
    p = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return p.stdout.strip()


def ffprobe_duration(video: pathlib.Path) -> float:
    out = run([
        "ffprobe",
        "-v",
        "error",
        "-show_entries",
        "format=duration",
        "-of",
        "default=noprint_wrappers=1:nokey=1",
        str(video),
    ])
    return float(out)


def ffprobe_fps(video: pathlib.Path) -> float:
    out = run([
        "ffprobe",
        "-v",
        "error",
        "-select_streams",
        "v:0",
        "-show_entries",
        "stream=r_frame_rate",
        "-of",
        "default=noprint_wrappers=1:nokey=1",
        str(video),
    ])
    # format like 24/1
    if "/" in out:
        a, b = out.split("/", 1)
        return float(a) / float(b)
    return float(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--video", default="/tmp/bh_renders/first_second/universe_first_second_120s_1080p.mp4")
    ap.add_argument("--out-dir", default="/tmp/bh_renders/first_second/ai_tiles")
    ap.add_argument("--samples", type=int, default=256)
    ap.add_argument("--cols", type=int, default=16)
    ap.add_argument("--tile-width", type=int, default=240)
    ap.add_argument("--frame-width", type=int, default=512)
    args = ap.parse_args()

    ffmpeg = shutil.which("ffmpeg")
    ffprobe = shutil.which("ffprobe")
    if not ffmpeg or not ffprobe:
        raise RuntimeError("ffmpeg/ffprobe not found in PATH")

    video = pathlib.Path(args.video)
    if not video.exists():
        raise FileNotFoundError(video)

    out_dir = pathlib.Path(args.out_dir)
    frames_dir = out_dir / "frames"
    out_dir.mkdir(parents=True, exist_ok=True)
    frames_dir.mkdir(parents=True, exist_ok=True)

    duration = ffprobe_duration(video)
    fps_src = ffprobe_fps(video)

    samples = max(1, int(args.samples))
    cols = max(1, int(args.cols))
    rows = int(math.ceil(samples / cols))
    fps_sample = samples / max(duration, 1e-9)

    # Extract sampled frames for AI inspection.
    subprocess.run(
        [
            ffmpeg,
            "-y",
            "-i",
            str(video),
            "-vf",
            f"fps={fps_sample:.8f},scale={args.frame_width}:-1",
            str(frames_dir / "frame_%04d.png"),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Build a dense contact sheet.
    contact = out_dir / "tile_contact_sheet.png"
    subprocess.run(
        [
            ffmpeg,
            "-y",
            "-i",
            str(video),
            "-vf",
            f"fps={fps_sample:.8f},scale={args.tile_width}:-1,tile={cols}x{rows}",
            "-frames:v",
            "1",
            str(contact),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Count emitted frames and compute timestamps mapping.
    frame_files = sorted(frames_dir.glob("frame_*.png"))
    n = len(frame_files)
    step = duration / max(n, 1)

    manifest = {
        "video": str(video),
        "duration_seconds": duration,
        "source_fps": fps_src,
        "samples_requested": samples,
        "samples_exported": n,
        "sample_step_seconds": step,
        "grid": {
            "cols": cols,
            "rows": rows,
            "tile_width": args.tile_width,
        },
        "frames": {
            "width": args.frame_width,
            "path": str(frames_dir),
            "pattern": "frame_%04d.png",
        },
        "artifacts": {
            "contact_sheet": str(contact),
            "frames_dir": str(frames_dir),
        },
    }
    manifest_path = out_dir / "tile_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"wrote {contact}")
    print(f"wrote {manifest_path}")
    print(f"wrote frames: {n} in {frames_dir}")


if __name__ == "__main__":
    main()
