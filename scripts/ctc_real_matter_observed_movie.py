#!/usr/bin/env python3
"""
Observed-style CTC matter renderer.

This renders matter as textured, glowing packets with visible transformations:
- ancestor packet (origin branch) remains coherent,
- target ancestor packet fragments during paradox intervention,
- traveler packet splits between local branch and loop channel at p*=0.5.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import subprocess
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def smoothstep(t: float) -> float:
    t = clamp(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def mix(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def mix_color(c0, c1, t: float):
    return tuple(int(mix(c0[i], c1[i], t)) for i in range(3))


def load_font(size: int, mono: bool = False):
    try:
        if mono:
            return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", size)
        return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", size)
    except Exception:
        return ImageFont.load_default()


def draw_starfield(img: Image.Image, seed: int = 137) -> None:
    rng = random.Random(seed)
    d = ImageDraw.Draw(img)
    w, h = img.size
    for _ in range(600):
        x = rng.randint(0, w - 1)
        y = rng.randint(0, h - 1)
        b = rng.randint(70, 180)
        d.point((x, y), fill=(b, b, min(255, b + 25)))
    for _ in range(120):
        x = rng.randint(0, w - 1)
        y = rng.randint(0, h - 1)
        b = rng.randint(170, 255)
        r = rng.randint(1, 2)
        d.ellipse((x - r, y - r, x + r, y + r), fill=(b, b, b))


def draw_lane_glow(img: Image.Image, y: int, color: tuple[int, int, int], width: int = 4) -> None:
    w, h = img.size
    layer = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    d.line((int(0.05 * w), y, int(0.95 * w), y), fill=(*color, 170), width=width)
    layer = layer.filter(ImageFilter.GaussianBlur(radius=7))
    img.alpha_composite(layer)


def draw_portal(img: Image.Image, x: float, y: float, r: float, phase: float) -> None:
    layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    for i in range(8):
        rr = r + i * 5.0
        a = int(100 - i * 10)
        col = (80, 220, 130) if i % 2 == 0 else (120, 180, 255)
        d.ellipse((x - rr, y - rr, x + rr, y + rr), outline=(*col, max(0, a)), width=2)
    # swirl spokes
    for k in range(10):
        ang = phase + k * (2.0 * math.pi / 10.0)
        x1 = x + math.cos(ang) * (r - 6)
        y1 = y + math.sin(ang) * (r - 6)
        x2 = x + math.cos(ang + 0.4) * (r + 18)
        y2 = y + math.sin(ang + 0.4) * (r + 18)
        d.line((x1, y1, x2, y2), fill=(130, 240, 160, 120), width=2)
    layer = layer.filter(ImageFilter.GaussianBlur(radius=1.7))
    img.alpha_composite(layer)


def draw_packet_core(
    img: Image.Image,
    x: float,
    y: float,
    radius: float,
    color: tuple[int, int, int],
    alpha: float = 1.0,
) -> None:
    alpha = clamp(alpha, 0.0, 1.0)
    layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    # Halo
    for i in range(8, 0, -1):
        rr = radius * (1.0 + i * 0.25)
        a = int(14 * i * alpha)
        d.ellipse((x - rr, y - rr, x + rr, y + rr), fill=(*color, a))
    # Core
    d.ellipse((x - radius, y - radius, x + radius, y + radius), fill=(*color, int(180 * alpha)))
    # Specular highlight
    d.ellipse(
        (x - radius * 0.55, y - radius * 0.60, x - radius * 0.10, y - radius * 0.15),
        fill=(245, 250, 255, int(90 * alpha)),
    )
    layer = layer.filter(ImageFilter.GaussianBlur(radius=2.2))
    img.alpha_composite(layer)


def draw_packet_texture(
    img: Image.Image,
    x: float,
    y: float,
    radius: float,
    t: float,
    seed: int,
    tone: tuple[int, int, int],
    alpha: float = 1.0,
) -> None:
    rng = random.Random(seed)
    layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    n = int(140 * alpha)
    for _ in range(n):
        a = rng.random() * 2.0 * math.pi
        rr = radius * math.sqrt(rng.random()) * 0.95
        px = x + math.cos(a + 1.4 * t) * rr
        py = y + math.sin(a - 1.1 * t) * rr
        s = rng.uniform(1.0, 2.4)
        c = mix_color((tone[0], tone[1], tone[2]), (255, 255, 255), rng.random() * 0.5)
        d.ellipse((px - s, py - s, px + s, py + s), fill=(*c, int(90 * alpha)))
    layer = layer.filter(ImageFilter.GaussianBlur(radius=0.6))
    img.alpha_composite(layer)


def draw_fragment_cloud(
    img: Image.Image,
    x: float,
    y: float,
    radius: float,
    progress: float,
    t: float,
    seed: int,
    color: tuple[int, int, int],
) -> None:
    progress = clamp(progress, 0.0, 1.0)
    if progress <= 0.0:
        return
    rng = random.Random(seed)
    layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    count = 260
    for _ in range(count):
        ang = rng.random() * 2.0 * math.pi
        base = rng.random()
        speed = mix(8.0, 95.0, progress) * (0.55 + base)
        drift = 8.0 * math.sin(t * 3.1 + ang)
        rr = radius * 0.5 + progress * speed + drift
        px = x + math.cos(ang) * rr
        py = y + math.sin(ang) * rr * 0.65
        s = mix(2.2, 1.0, progress) * (0.7 + 0.6 * rng.random())
        a = int(mix(120, 25, progress))
        c = mix_color(color, (255, 220, 210), rng.random() * 0.45)
        d.ellipse((px - s, py - s, px + s, py + s), fill=(*c, a))
    layer = layer.filter(ImageFilter.GaussianBlur(radius=1.0 + 1.4 * progress))
    img.alpha_composite(layer)


def draw_trajectory_trail(
    img: Image.Image,
    points: list[tuple[float, float]],
    color: tuple[int, int, int],
) -> None:
    layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    if len(points) >= 2:
        d.line(points, fill=(*color, 130), width=3)
    layer = layer.filter(ImageFilter.GaussianBlur(radius=2.5))
    img.alpha_composite(layer)


def main() -> None:
    ap = argparse.ArgumentParser(description="Render observed-style CTC matter movie")
    ap.add_argument(
        "--report-json",
        default="/tmp/bh_renders/ctc_real_matter_paradox/ctc_real_matter_paradox_report.json",
    )
    ap.add_argument("--out-dir", default="/tmp/bh_renders/ctc_real_matter_paradox")
    ap.add_argument("--work-dir", default="", help="Frame work dir (override /tmp usage)")
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--seconds", type=float, default=16.0)
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--height", type=int, default=1080)
    ap.add_argument("--keep-frames", action="store_true")
    args = ap.parse_args()

    in_path = Path(args.report_json)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(args.work_dir) if args.work_dir else out_dir / "frames_ctc_real_matter_observed"
    work_dir.mkdir(parents=True, exist_ok=True)

    data = json.loads(in_path.read_text())
    m_a = float(data["inputs"]["m_ancestor_kg"])
    m_t = float(data["inputs"]["m_traveler_kg"])
    p_star = float(data["test_3_deutsch_fixed_point"]["p_star"])
    freq = float(data["test_4_monte_carlo_real_matter"]["freq_traveler_present"])

    title_font = load_font(42, mono=False)
    body_font = load_font(22, mono=True)
    small_font = load_font(18, mono=True)

    w, h = args.width, args.height
    y_origin = int(0.26 * h)
    y_target = int(0.53 * h)
    y_channel = int(0.79 * h)
    x_l = int(0.09 * w)
    x_r = int(0.90 * w)

    n_frames = max(1, int(args.seconds * args.fps))
    keyframe_idx = int(0.73 * (n_frames - 1))

    for i in range(n_frames):
        tn = i / max(1, n_frames - 1)
        t_real = tn * args.seconds
        base = Image.new("RGBA", (w, h), (10, 14, 24, 255))
        draw_starfield(base, seed=137)

        # Gradient wash
        wash = Image.new("RGBA", (w, h), (0, 0, 0, 0))
        wd = ImageDraw.Draw(wash)
        wd.ellipse((int(-0.2 * w), int(-0.3 * h), int(0.55 * w), int(0.55 * h)), fill=(40, 55, 95, 65))
        wd.ellipse((int(0.35 * w), int(0.35 * h), int(1.25 * w), int(1.2 * h)), fill=(25, 60, 45, 60))
        wash = wash.filter(ImageFilter.GaussianBlur(radius=45))
        base.alpha_composite(wash)

        draw_lane_glow(base, y_origin, (95, 145, 255), width=4)
        draw_lane_glow(base, y_target, (255, 180, 95), width=4)
        draw_lane_glow(base, y_channel, (120, 220, 135), width=4)

        d = ImageDraw.Draw(base)
        d.rectangle((0, 0, w, 86), fill=(18, 26, 43, 220))
        d.text((32, 22), "Observed Matter View: CTC Branch and Packet Transformations", fill=(234, 242, 255), font=title_font)

        # Position curves
        x_prog = smoothstep(tn)
        xo = mix(x_l, x_r, x_prog)
        xt = mix(x_l, x_r, x_prog)
        xl = mix(int(0.2 * w), int(0.62 * w), smoothstep((tn - 0.15) / 0.60))
        xc = mix(int(0.28 * w), int(0.80 * w), smoothstep((tn - 0.50) / 0.45))

        # Trails
        draw_trajectory_trail(base, [(x_l, y_origin), (xo, y_origin)], (115, 160, 255))
        draw_trajectory_trail(base, [(x_l, y_target), (xt, y_target)], (255, 190, 120))
        draw_trajectory_trail(base, [(int(0.2 * w), y_target - 45), (xl, y_target - 45)], (130, 190, 255))
        draw_trajectory_trail(base, [(int(0.28 * w), y_channel), (xc, y_channel)], (135, 230, 150))

        # Origin ancestor: coherent packet
        draw_packet_core(base, xo, y_origin, 25, (95, 220, 155), alpha=1.0)
        draw_packet_texture(base, xo, y_origin, 21, t_real, seed=11, tone=(85, 190, 130), alpha=1.0)

        # Target ancestor: transforms (alive -> fragmented)
        kill_p = smoothstep((tn - 0.43) / 0.20)
        intact_alpha = 1.0 - 0.82 * kill_p
        target_col = mix_color((100, 215, 145), (215, 95, 90), kill_p)
        draw_packet_core(base, xt, y_target, 24, target_col, alpha=intact_alpha)
        draw_packet_texture(base, xt, y_target, 20, t_real * 1.2, seed=21, tone=target_col, alpha=intact_alpha)
        draw_fragment_cloud(base, xt, y_target, 22, progress=kill_p, t=t_real, seed=33, color=(230, 120, 100))

        # Traveler local packet
        traveler_alpha = smoothstep((tn - 0.18) / 0.18)
        draw_packet_core(base, xl, y_target - 45, 19, (105, 175, 255), alpha=traveler_alpha)
        draw_packet_texture(base, xl, y_target - 45, 16, t_real * 1.5, seed=44, tone=(95, 160, 245), alpha=traveler_alpha)

        # Loop-channel packet + portal
        portal_phase = 1.8 * t_real
        draw_portal(base, int(0.23 * w), y_channel, 36, portal_phase)
        channel_alpha = smoothstep((tn - 0.54) / 0.22)
        draw_packet_core(base, xc, y_channel, 18, (125, 230, 140), alpha=channel_alpha)
        draw_packet_texture(base, xc, y_channel, 15, t_real * 1.4, seed=55, tone=(110, 220, 130), alpha=channel_alpha)

        # Transfer arc from local to channel
        if tn > 0.53:
            arc = Image.new("RGBA", (w, h), (0, 0, 0, 0))
            ad = ImageDraw.Draw(arc)
            x1, y1 = xl + 16, y_target - 30
            x2, y2 = xc - 10, y_channel - 8
            ctrlx, ctrly = mix(x1, x2, 0.5), y1 - 130
            pts = []
            for k in range(80):
                u = k / 79.0
                px = (1 - u) ** 2 * x1 + 2 * (1 - u) * u * ctrlx + u * u * x2
                py = (1 - u) ** 2 * y1 + 2 * (1 - u) * u * ctrly + u * u * y2
                pts.append((px, py))
            ad.line(pts, fill=(145, 240, 170, 170), width=3)
            arc = arc.filter(ImageFilter.GaussianBlur(radius=2.0))
            base.alpha_composite(arc)

        # Labels + values
        d = ImageDraw.Draw(base)
        d.text((48, y_origin - 38), f"Origin ancestor packet: {m_a:.0f} kg (coherent)", fill=(190, 225, 255), font=body_font)
        d.text((48, y_target - 38), "Target ancestor packet: intervention transform", fill=(245, 210, 170), font=body_font)
        d.text((48, y_channel - 38), "Loop channel packet: traveler complement branch", fill=(190, 240, 190), font=body_font)

        local_m = m_a + p_star * m_t
        channel_m = (1.0 - p_star) * m_t
        total_m = local_m + channel_m
        d.rounded_rectangle((int(0.69 * w), int(0.72 * h), int(0.97 * w), int(0.95 * h)), radius=14, fill=(22, 34, 50, 220), outline=(95, 140, 205), width=2)
        d.text((int(0.705 * w), int(0.745 * h)), f"p* = {p_star:.3f}", fill=(220, 235, 255), font=body_font)
        d.text((int(0.705 * w), int(0.785 * h)), f"local mass   {local_m:6.1f} kg", fill=(165, 210, 255), font=small_font)
        d.text((int(0.705 * w), int(0.815 * h)), f"channel mass {channel_m:6.1f} kg", fill=(175, 235, 175), font=small_font)
        d.text((int(0.705 * w), int(0.845 * h)), f"total mass   {total_m:6.1f} kg", fill=(225, 245, 225), font=small_font)
        d.text((int(0.705 * w), int(0.875 * h)), f"MC freq      {freq:.6f}", fill=(240, 220, 180), font=small_font)

        d.text((34, h - 32), "Visual lane: matter packets + fragmentation + branch/channel transfer (CTC fixed-point ensemble)", fill=(155, 176, 208), font=small_font)

        frame_path = work_dir / f"frame_{i:05d}.png"
        base.convert("RGB").save(frame_path)
        if i == keyframe_idx:
            base.convert("RGB").save(out_dir / "ctc_real_matter_observed_keyframe.png")

    out_mp4 = out_dir / "ctc_real_matter_observed.mp4"
    cmd = [
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
        "17",
        str(out_mp4),
    ]
    subprocess.run(cmd, check=True)

    if not args.keep_frames:
        for p in work_dir.glob("frame_*.png"):
            p.unlink()
        work_dir.rmdir()

    print(out_mp4)
    print(out_dir / "ctc_real_matter_observed_keyframe.png")


if __name__ == "__main__":
    main()
