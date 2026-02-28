#!/usr/bin/env python3
"""
Render CTC real-matter paradox transformations as an animation.

Reads the JSON report from ctc_real_matter_paradox_sim and produces:
- ctc_real_matter_transform.mp4
- ctc_real_matter_transform_keyframe.png
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path
from typing import Tuple

from PIL import Image, ImageDraw, ImageFont


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def ease_in_out(t: float) -> float:
    t = clamp(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def load_fonts() -> Tuple[ImageFont.FreeTypeFont, ImageFont.FreeTypeFont, ImageFont.FreeTypeFont]:
    try:
        title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 44)
        body = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", 24)
        small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", 20)
        return title, body, small
    except Exception:
        default = ImageFont.load_default()
        return default, default, default


def draw_packet(
    draw: ImageDraw.ImageDraw,
    center: Tuple[float, float],
    radius: float,
    fill: Tuple[int, int, int],
    outline: Tuple[int, int, int],
    alpha_scale: float = 1.0,
) -> None:
    x, y = center
    r = radius
    a = clamp(alpha_scale, 0.0, 1.0)
    col = tuple(int(c * a) for c in fill)
    out = tuple(int(c * a) for c in outline)
    draw.ellipse((x - r, y - r, x + r, y + r), fill=col, outline=out, width=3)


def lane_points(w: int, h: int):
    # Origin branch (top), target branch (mid), loop channel (bottom)
    return {
        "origin_left": (int(0.10 * w), int(0.26 * h)),
        "origin_mid": (int(0.45 * w), int(0.26 * h)),
        "origin_right": (int(0.80 * w), int(0.26 * h)),
        "target_left": (int(0.10 * w), int(0.53 * h)),
        "target_mid": (int(0.45 * w), int(0.53 * h)),
        "target_right": (int(0.80 * w), int(0.53 * h)),
        "channel_left": (int(0.12 * w), int(0.79 * h)),
        "channel_right": (int(0.78 * w), int(0.79 * h)),
    }


def draw_background(draw: ImageDraw.ImageDraw, w: int, h: int) -> None:
    draw.rectangle((0, 0, w, h), fill=(15, 20, 34))
    draw.rectangle((0, 0, w, 96), fill=(24, 34, 54))
    draw.rounded_rectangle((40, 140, w - 40, int(0.38 * h)), radius=16, outline=(100, 135, 210), width=2)
    draw.rounded_rectangle((40, int(0.44 * h), w - 40, int(0.64 * h)), radius=16, outline=(210, 170, 90), width=2)
    draw.rounded_rectangle((40, int(0.70 * h), w - 40, h - 36), radius=16, outline=(120, 170, 120), width=2)


def draw_lanes(draw: ImageDraw.ImageDraw, p, small_font) -> None:
    draw.line((*p["origin_left"], *p["origin_right"]), fill=(90, 140, 230), width=5)
    draw.line((*p["target_left"], *p["target_right"]), fill=(230, 180, 90), width=5)
    draw.line((*p["channel_left"], *p["channel_right"]), fill=(120, 180, 120), width=5)
    draw.text((p["origin_left"][0], p["origin_left"][1] - 36), "Origin branch O", fill=(210, 225, 255), font=small_font)
    draw.text((p["target_left"][0], p["target_left"][1] - 36), "Target branch T", fill=(245, 220, 180), font=small_font)
    draw.text((p["channel_left"][0], p["channel_left"][1] - 36), "Loop channel", fill=(210, 240, 210), font=small_font)


def stage_text(t: float) -> str:
    if t < 0.18:
        return "Stage 1: setup (single history seed)"
    if t < 0.36:
        return "Stage 2: branch split (O and T both instantiated)"
    if t < 0.58:
        return "Stage 3: intervention branch (ancestor state transforms)"
    if t < 0.80:
        return "Stage 4: traveler packet in loop channel complement"
    return "Stage 5: Deutsch fixed-point ensemble (p*=0.5)"


def main() -> None:
    ap = argparse.ArgumentParser(description="Render CTC real-matter transformation movie")
    ap.add_argument(
        "--report-json",
        default="/tmp/bh_renders/ctc_real_matter_paradox/ctc_real_matter_paradox_report.json",
        help="Input simulation JSON",
    )
    ap.add_argument(
        "--out-dir",
        default="/tmp/bh_renders/ctc_real_matter_paradox",
        help="Output directory",
    )
    ap.add_argument("--work-dir", default="", help="Temporary frame directory (defaults under out-dir)")
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--seconds", type=float, default=14.0)
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--height", type=int, default=1080)
    ap.add_argument("--keep-frames", action="store_true")
    args = ap.parse_args()

    in_path = Path(args.report_json)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(args.work_dir) if args.work_dir else out_dir / "frames_ctc_real_matter"
    work_dir.mkdir(parents=True, exist_ok=True)

    data = json.loads(in_path.read_text())
    m_a = float(data["inputs"]["m_ancestor_kg"])
    m_t = float(data["inputs"]["m_traveler_kg"])
    p_star = float(data["test_3_deutsch_fixed_point"]["p_star"])
    freq = float(data["test_4_monte_carlo_real_matter"]["freq_traveler_present"])

    title_font, body_font, small_font = load_fonts()
    w, h = args.width, args.height
    p = lane_points(w, h)

    n_frames = max(1, int(args.fps * args.seconds))
    key_idx = n_frames - 1

    for i in range(n_frames):
        tn = i / max(1, n_frames - 1)
        img = Image.new("RGB", (w, h))
        draw = ImageDraw.Draw(img)

        draw_background(draw, w, h)
        draw_lanes(draw, p, small_font)

        draw.text((42, 24), "CTC Matter Transform: Branch + Loop-Channel Bookkeeping", fill=(235, 243, 255), font=title_font)
        draw.text((44, 96), stage_text(tn), fill=(180, 204, 240), font=body_font)

        # Ancestor packet in origin branch: persists in O.
        a_x = lerp(p["origin_left"][0], p["origin_right"][0], tn)
        a_y = p["origin_left"][1]
        draw_packet(draw, (a_x, a_y), 22, (90, 215, 150), (230, 255, 240), 1.0)
        draw.text((a_x - 52, a_y + 30), f"A(O) {m_a:.0f}kg", fill=(200, 245, 225), font=small_font)

        # Target ancestor state transforms from alive -> dead around stage 3.
        t_x = lerp(p["target_left"][0], p["target_right"][0], tn)
        t_y = p["target_left"][1]
        kill_progress = ease_in_out((tn - 0.40) / 0.20)
        alive_col = (
            int(lerp(80, 145, 1 - kill_progress)),
            int(lerp(215, 80, kill_progress)),
            int(lerp(145, 90, kill_progress)),
        )
        draw_packet(draw, (t_x, t_y), 22, alive_col, (245, 230, 200), 1.0)
        state_label = "alive" if kill_progress < 0.5 else "dead"
        draw.text((t_x - 64, t_y + 30), f"A(T) {state_label}", fill=(245, 220, 200), font=small_font)

        # Traveler packet branches:
        # - local on target branch with weight p*=0.5
        # - channel packet with weight 1-p*=0.5
        traveler_alpha = 1.0 if tn > 0.25 else ease_in_out(tn / 0.25)
        local_x = lerp(p["target_left"][0] + 80, p["target_mid"][0] + 40, tn)
        local_y = p["target_left"][1] - 46
        draw_packet(draw, (local_x, local_y), 18, (95, 170, 255), (230, 245, 255), traveler_alpha)
        draw.text((local_x - 100, local_y - 44), f"T local {p_star:.1f} x {m_t:.0f}kg", fill=(200, 225, 255), font=small_font)

        channel_x = lerp(p["channel_left"][0] + 60, p["channel_right"][0] - 90, tn)
        channel_y = p["channel_left"][1]
        draw_packet(draw, (channel_x, channel_y), 18, (110, 220, 120), (240, 255, 240), traveler_alpha)
        draw.text((channel_x - 126, channel_y + 30), f"T channel {(1-p_star):.1f} x {m_t:.0f}kg", fill=(205, 245, 205), font=small_font)

        # Flow arrows for transformations
        if tn >= 0.30:
            draw.line((local_x - 26, local_y + 14, t_x - 20, t_y - 6), fill=(150, 200, 255), width=3)
            draw.polygon([(t_x - 20, t_y - 6), (t_x - 32, t_y - 10), (t_x - 28, t_y + 3)], fill=(150, 200, 255))
        if tn >= 0.55:
            draw.line((local_x + 20, local_y + 20, channel_x - 16, channel_y - 16), fill=(145, 220, 145), width=3)
            draw.polygon([(channel_x - 16, channel_y - 16), (channel_x - 28, channel_y - 18), (channel_x - 20, channel_y - 6)], fill=(145, 220, 145))

        # Mass ledger panel
        panel_x0, panel_y0 = int(0.72 * w), int(0.72 * h)
        panel_x1, panel_y1 = w - 50, h - 56
        draw.rounded_rectangle((panel_x0, panel_y0, panel_x1, panel_y1), radius=12, fill=(28, 38, 56), outline=(110, 145, 200), width=2)
        m_local = m_a + p_star * m_t
        m_channel = (1.0 - p_star) * m_t
        m_total = m_local + m_channel
        draw.text((panel_x0 + 14, panel_y0 + 10), "Mass ledger", fill=(225, 238, 255), font=small_font)
        draw.text((panel_x0 + 14, panel_y0 + 42), f"local   = {m_local:.1f} kg", fill=(170, 205, 255), font=small_font)
        draw.text((panel_x0 + 14, panel_y0 + 70), f"channel = {m_channel:.1f} kg", fill=(180, 235, 180), font=small_font)
        draw.text((panel_x0 + 14, panel_y0 + 98), f"total   = {m_total:.1f} kg", fill=(220, 248, 220), font=small_font)
        draw.text((panel_x0 + 14, panel_y0 + 126), f"MC freq = {freq:.6f}", fill=(230, 220, 180), font=small_font)

        # footer theorem hook
        draw.text(
            (44, h - 36),
            "notMap fixed point: p = 1 - p => p = 1/2  |  traveler packet conserved across local+channel",
            fill=(164, 183, 214),
            font=small_font,
        )

        frame_path = work_dir / f"frame_{i:05d}.png"
        img.save(frame_path)

        if i == key_idx:
            img.save(out_dir / "ctc_real_matter_transform_keyframe.png")

    out_mp4 = out_dir / "ctc_real_matter_transform.mp4"
    ffmpeg_cmd = [
        "ffmpeg",
        "-y",
        "-framerate",
        str(args.fps),
        "-i",
        str(work_dir / "frame_%05d.png"),
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        "18",
        str(out_mp4),
    ]
    subprocess.run(ffmpeg_cmd, check=True)

    if not args.keep_frames:
        for pth in work_dir.glob("frame_*.png"):
            pth.unlink()
        work_dir.rmdir()

    print(out_mp4)
    print((out_dir / "ctc_real_matter_transform_keyframe.png"))


if __name__ == "__main__":
    main()
