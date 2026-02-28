#!/usr/bin/env python3
"""
Render a cinematic 3D CMB movie with transparent tomographic shells.

Default outputs:
- /tmp/bh_renders/cmb_3d_transparent_movie.mp4
- /tmp/bh_renders/cmb_3d_transparent_movie.gif
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import shutil
import subprocess
import tempfile

import matplotlib
import numpy as np
from scipy.special import sph_harm

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm  # noqa: E402


def read_tau_derived(path: pathlib.Path) -> float:
    if path.exists():
        try:
            d = json.loads(path.read_text(encoding="utf-8"))
            return float(d["reionization"]["tau_reio_derived"])
        except Exception:
            pass
    return 0.067531477


def run_class(class_bin: str, params: dict[str, float], lmax: int = 2500) -> pathlib.Path:
    td = tempfile.mkdtemp(prefix="cmb_movie_class_")
    run_dir = pathlib.Path(td)
    ini = run_dir / "run.ini"
    root = run_dir / "g_"
    ini.write_text(
        "\n".join(
            [
                f"h = {params['h']}",
                f"omega_b = {params['omega_b']}",
                f"omega_cdm = {params['omega_cdm']}",
                f"Omega_k = {params['Omega_k']}",
                f"Omega_Lambda = {params['Omega_Lambda']}",
                f"A_s = {params['A_s']}",
                f"n_s = {params['n_s']}",
                f"tau_reio = {params['tau_reio']}",
                "output=tCl,lCl,pCl",
                "lensing=yes",
                f"l_max_scalars={lmax}",
                "format=camb",
                f"root = {root}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    subprocess.run([class_bin, str(ini)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    cands = sorted(run_dir.glob("*cl*.dat"))
    cands.sort(key=lambda p: (0 if "lensedcls" in p.name.lower() else 1))
    if not cands:
        raise RuntimeError("No CLASS C_l output found")
    return cands[0]


def parse_class_dl(path: pathlib.Path, col_idx_1based: int, lmin: int = 2, lmax: int = 2500):
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        if len(f) < col_idx_1based:
            continue
        try:
            ell = int(float(f[0]))
            val = float(f[col_idx_1based - 1])
        except ValueError:
            continue
        if lmin <= ell <= lmax:
            out.append((ell, val))
    out.sort(key=lambda t: t[0])
    return out


def cl_from_dl(curve: list[tuple[int, float]], lmax: int):
    cl = np.zeros(lmax + 1, dtype=float)
    for ell, dl in curve:
        if 2 <= ell <= lmax:
            cl[ell] = dl * (2.0 * math.pi) / (ell * (ell + 1))
    return cl


def synth_map_from_cl(cl_tt: np.ndarray, lmax: int, nlat: int, nlon: int, seed: int):
    rng = np.random.default_rng(seed)
    theta = np.linspace(0.0, math.pi, nlat)
    phi = np.linspace(0.0, 2.0 * math.pi, nlon, endpoint=False)
    PHI, THETA = np.meshgrid(phi, theta)
    sky = np.zeros_like(PHI, dtype=float)

    for ell in range(2, lmax + 1):
        c = float(max(cl_tt[ell], 0.0))
        if c <= 0.0:
            continue
        a_l0 = rng.normal(scale=math.sqrt(c))
        y0 = sph_harm(0, ell, PHI, THETA).real
        sky += a_l0 * y0
        s = math.sqrt(c / 2.0)
        for m in range(1, ell + 1):
            a = rng.normal(scale=s) + 1j * rng.normal(scale=s)
            y = sph_harm(m, ell, PHI, THETA)
            sky += 2.0 * np.real(a * y)

    sky -= np.mean(sky)
    return sky, theta, phi


def render_frame(frame_path: pathlib.Path, norm: np.ndarray, lat: np.ndarray, lon: np.ndarray, frame_idx: int, frames: int):
    LON, LAT = np.meshgrid(lon, lat)
    phase = 2.0 * math.pi * frame_idx / max(frames, 1)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    ax = fig.add_subplot(111, projection="3d")

    shells = np.linspace(0.88, 1.12, 9)
    for i, base_r in enumerate(shells):
        shell_weight = 1.0 - abs((i - (len(shells) - 1) / 2.0) / ((len(shells) - 1) / 2.0))
        breathe = 1.0 + 0.08 * math.sin(phase + i * 0.55)
        r = base_r + 0.020 * norm * breathe

        X = r * np.cos(LAT) * np.cos(LON)
        Y = r * np.cos(LAT) * np.sin(LON)
        Z = r * np.sin(LAT)

        rgba = cm.coolwarm((norm + 1.0) / 2.0)
        alpha = 0.06 + 0.20 * shell_weight + 0.07 * np.abs(norm)
        rgba[..., 3] = np.clip(alpha, 0.05, 0.40)

        ax.plot_surface(
            X,
            Y,
            Z,
            facecolors=rgba,
            rstride=3,
            cstride=3,
            linewidth=0,
            antialiased=False,
            shade=False,
        )

    for ref_r in [0.9, 1.0, 1.1]:
        Xw = ref_r * np.cos(LAT) * np.cos(LON)
        Yw = ref_r * np.cos(LAT) * np.sin(LON)
        Zw = ref_r * np.sin(LAT)
        ax.plot_wireframe(Xw, Yw, Zw, rstride=14, cstride=18, color=(0.10, 0.10, 0.10, 0.12), linewidth=0.4)

    elev = 20.0 + 5.5 * math.sin(phase * 0.5)
    azim = 35.0 + 360.0 * frame_idx / max(frames, 1)
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.set_title("CMB 3D Transparent Tomography — Derived TT Lane", pad=12)

    mappable = matplotlib.cm.ScalarMappable(cmap="coolwarm")
    mappable.set_array(norm)
    cbar = fig.colorbar(mappable, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("delta T (arb)")

    fig.tight_layout()
    fig.savefig(frame_path)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="/tmp/bh_renders", help="output directory")
    p.add_argument("--class-bin", default="/tmp/class_public/class", help="CLASS binary")
    p.add_argument("--frames", type=int, default=120, help="number of frames")
    p.add_argument("--fps", type=int, default=24, help="frames per second")
    p.add_argument("--lmax-map", type=int, default=36, help="lmax for synthetic map")
    p.add_argument("--nlat", type=int, default=120)
    p.add_argument("--nlon", type=int, default=240)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--skip-gif", action="store_true")
    args = p.parse_args()

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tau = read_tau_derived(out_dir / "cmb_tau_derived_report.json")
    params = {
        "h": 0.680163311753,
        "omega_b": 0.022807271041,
        "omega_cdm": 0.124403296589,
        "Omega_k": 0.0,
        "Omega_Lambda": 0.681700909091,
        "n_s": 0.965416666667,
        "A_s": 2.219284522994e-9,
        "tau_reio": tau,
    }

    class_out = run_class(args.class_bin, params)
    pred_tt = parse_class_dl(class_out, 2)
    cl_tt = cl_from_dl(pred_tt, lmax=max(args.lmax_map, 2))
    sky, theta, phi = synth_map_from_cl(
        cl_tt,
        lmax=args.lmax_map,
        nlat=args.nlat,
        nlon=args.nlon,
        seed=args.seed,
    )

    sky_plot = np.roll(sky, sky.shape[1] // 2, axis=1)
    v = np.percentile(np.abs(sky_plot), 99.0)
    norm = np.clip(sky_plot / max(v, 1e-9), -1.0, 1.0)
    lat = (math.pi / 2.0) - theta
    lon = phi - math.pi

    frame_dir = pathlib.Path(tempfile.mkdtemp(prefix="cmb_movie_frames_"))
    try:
        for i in range(args.frames):
            frame_path = frame_dir / f"frame_{i:04d}.png"
            render_frame(frame_path, norm, lat, lon, i, args.frames)

        mp4_path = out_dir / "cmb_3d_transparent_movie.mp4"
        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise RuntimeError("ffmpeg not found in PATH")

        subprocess.run(
            [
                ffmpeg,
                "-y",
                "-framerate",
                str(args.fps),
                "-i",
                str(frame_dir / "frame_%04d.png"),
                "-c:v",
                "libx264",
                "-pix_fmt",
                "yuv420p",
                str(mp4_path),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        gif_path = out_dir / "cmb_3d_transparent_movie.gif"
        if not args.skip_gif:
            palette = frame_dir / "palette.png"
            subprocess.run(
                [
                    ffmpeg,
                    "-y",
                    "-i",
                    str(frame_dir / "frame_%04d.png"),
                    "-vf",
                    "palettegen",
                    str(palette),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                [
                    ffmpeg,
                    "-y",
                    "-framerate",
                    str(args.fps),
                    "-i",
                    str(frame_dir / "frame_%04d.png"),
                    "-i",
                    str(palette),
                    "-lavfi",
                    "paletteuse",
                    str(gif_path),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

        summary = {
            "tau_reio": tau,
            "class_output": str(class_out),
            "frames": args.frames,
            "fps": args.fps,
            "lmax_map": args.lmax_map,
            "nlat": args.nlat,
            "nlon": args.nlon,
            "seed": args.seed,
            "artifacts": {
                "mp4": str(mp4_path),
                "gif": str(gif_path) if not args.skip_gif else None,
            },
        }
        (out_dir / "cmb_3d_transparent_movie_summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )

        print("wrote", mp4_path)
        if not args.skip_gif:
            print("wrote", gif_path)
        print("wrote", out_dir / "cmb_3d_transparent_movie_summary.json")
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
