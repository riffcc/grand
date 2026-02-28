#!/usr/bin/env python3
"""
CMB holy-shit gallery + 3D render (derived tau lane).

Outputs (default: /tmp/bh_renders):
- cmb_holyshit_spectra_triptych.png
- cmb_holyshit_pull_triptych.png
- cmb_holyshit_projection_gallery.png
- cmb_holyshit_3d_render.png
- cmb_holyshit_summary.json
"""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import tempfile
from dataclasses import dataclass

import matplotlib
import numpy as np
from scipy.special import sph_harm

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm  # noqa: E402


@dataclass(frozen=True)
class SpectrumPoint:
    ell: int
    d_ell: float
    sigma: float


def read_tau_derived(path: pathlib.Path) -> float:
    if path.exists():
        try:
            d = json.loads(path.read_text(encoding="utf-8"))
            return float(d["reionization"]["tau_reio_derived"])
        except Exception:
            pass
    return 0.067531477


def run_class(class_bin: str, params: dict[str, float], lmax: int = 2500) -> pathlib.Path:
    td = tempfile.mkdtemp(prefix="cmb_holyshit_")
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


def load_planck(path: pathlib.Path, lmin: int = 2, lmax: int = 2500):
    pts = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        if len(f) < 4:
            continue
        ell = int(round(float(f[0])))
        if not (lmin <= ell <= lmax):
            continue
        dl = float(f[1])
        sigma = 0.5 * (abs(float(f[2])) + abs(float(f[3])))
        pts.append(SpectrumPoint(ell=ell, d_ell=dl, sigma=sigma))
    pts.sort(key=lambda p: p.ell)
    return pts


def interp(curve: list[tuple[int, float]], ell: int) -> float | None:
    if ell < curve[0][0] or ell > curve[-1][0]:
        return None
    for (xa, ya), (xb, yb) in zip(curve, curve[1:]):
        if xa <= ell <= xb:
            t = (ell - xa) / (xb - xa) if xb > xa else 0.0
            return ya * (1.0 - t) + yb * t
    return curve[-1][1]


def channel_stats(pred: list[tuple[int, float]], obs: list[SpectrumPoint]):
    chi2 = 0.0
    n = 0
    pulls = []
    for p in obs:
        y = interp(pred, p.ell)
        if y is None:
            continue
        pull = (y - p.d_ell) / p.sigma
        chi2 += pull * pull
        pulls.append((p.ell, pull))
        n += 1
    red = chi2 / max(1, n - 1)
    return chi2, red, pulls


def cl_from_dl(curve: list[tuple[int, float]], lmax: int):
    cl = np.zeros(lmax + 1, dtype=float)
    for ell, dl in curve:
        if ell <= lmax and ell >= 2:
            cl[ell] = dl * (2.0 * math.pi) / (ell * (ell + 1))
    return cl


def synth_map_from_cl(cl_tt: np.ndarray, lmax: int = 36, nlat: int = 180, nlon: int = 360, seed: int = 42):
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


def plot_spectra_triptych(out: pathlib.Path, pred_tt, pred_te, pred_ee, obs_tt_b, obs_te_b, obs_ee_b, obs_tt_f, obs_te_f, obs_ee_f):
    fig, axs = plt.subplots(3, 1, figsize=(13, 12), sharex=True)
    channels = [
        ("TT", pred_tt, obs_tt_b, obs_tt_f, axs[0], (0, 6800)),
        ("TE", pred_te, obs_te_b, obs_te_f, axs[1], (-220, 220)),
        ("EE", pred_ee, obs_ee_b, obs_ee_f, axs[2], (0, 65)),
    ]
    for name, pred, ob, of, ax, ylim in channels:
        ax.scatter([p.ell for p in of], [p.d_ell for p in of], s=4, alpha=0.18, color="#7f8c8d", label=f"Planck {name} full")
        ax.errorbar([p.ell for p in ob], [p.d_ell for p in ob], yerr=[p.sigma for p in ob], fmt="o", markersize=3, linewidth=0.7,
                    color="#f39c12", ecolor="#f39c12", alpha=0.9, label=f"Planck {name} binned")
        ax.plot([e for e, _ in pred], [v for _, v in pred], lw=2.0, color="#2e86de", label=f"GUTOE CLASS {name}")
        ax.set_ylabel(rf"$D_\ell^{{{name}}}\ [\mu K^2]$")
        ax.set_ylim(*ylim)
        ax.grid(alpha=0.2)
        ax.legend(loc="upper right", fontsize=8)
    axs[-1].set_xlabel(r"Multipole $\ell$")
    axs[0].set_title("CMB TT/TE/EE vs Planck (derived tau lane)")
    fig.tight_layout()
    fig.savefig(out / "cmb_holyshit_spectra_triptych.png", dpi=180)


def plot_pull_triptych(out: pathlib.Path, pulls_tt, pulls_te, pulls_ee):
    fig, axs = plt.subplots(3, 1, figsize=(13, 8), sharex=True)
    for (name, pulls, ax) in [("TT", pulls_tt, axs[0]), ("TE", pulls_te, axs[1]), ("EE", pulls_ee, axs[2])]:
        ax.axhline(0, color="black", lw=1)
        ax.plot([e for e, _ in pulls], [p for _, p in pulls], color="#2e86de", lw=1.2)
        ax.set_ylabel(f"{name} pull")
        ax.grid(alpha=0.2)
    axs[-1].set_xlabel(r"Multipole $\ell$")
    axs[0].set_title("Residual pulls (binned) — derived tau lane")
    fig.tight_layout()
    fig.savefig(out / "cmb_holyshit_pull_triptych.png", dpi=180)


def plot_projection_gallery(out: pathlib.Path, sky: np.ndarray, theta: np.ndarray, phi: np.ndarray):
    lat = (math.pi / 2.0) - theta
    lon = phi - math.pi
    LON, LAT = np.meshgrid(lon, lat)

    sky_plot = np.roll(sky, sky.shape[1] // 2, axis=1)
    v = np.percentile(np.abs(sky_plot), 99.0)

    fig = plt.figure(figsize=(16, 12))

    ax1 = fig.add_subplot(2, 3, 1)
    im = ax1.imshow(sky_plot, origin="lower", cmap="coolwarm", vmin=-v, vmax=v, aspect="auto")
    ax1.set_title("Equirectangular")
    ax1.set_xlabel("lon")
    ax1.set_ylabel("lat")

    for i, proj in enumerate(["mollweide", "aitoff", "hammer", "lambert"], start=2):
        ax = fig.add_subplot(2, 3, i, projection=proj)
        ax.pcolormesh(LON, LAT, sky_plot, shading="auto", cmap="coolwarm", vmin=-v, vmax=v)
        ax.set_title(proj.capitalize())
        ax.grid(True, alpha=0.25)

    ax6 = fig.add_subplot(2, 3, 6)
    ax6.hist(sky_plot.ravel(), bins=120, color="#2e86de", alpha=0.8)
    ax6.set_title("Temperature anisotropy histogram")
    ax6.set_xlabel("delta T [muK] (arb scale)")
    ax6.set_ylabel("count")

    cbar = fig.colorbar(im, ax=fig.axes, fraction=0.02, pad=0.02)
    cbar.set_label("delta T")
    fig.suptitle("CMB projection gallery (single realization from derived TT spectrum)", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out / "cmb_holyshit_projection_gallery.png", dpi=180)


def plot_3d_render(out: pathlib.Path, sky: np.ndarray, theta: np.ndarray, phi: np.ndarray):
    lat = (math.pi / 2.0) - theta
    lon = phi - math.pi
    LON, LAT = np.meshgrid(lon, lat)

    sky_plot = np.roll(sky, sky.shape[1] // 2, axis=1)
    v = np.percentile(np.abs(sky_plot), 99.0)
    norm = np.clip(sky_plot / max(v, 1e-9), -1.0, 1.0)

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Multi-shell transparent render: lets you see through successive slices
    # while preserving 3D depth cues.
    shells = np.linspace(0.88, 1.12, 9)
    for i, base_r in enumerate(shells):
        shell_weight = 1.0 - abs((i - (len(shells) - 1) / 2.0) / ((len(shells) - 1) / 2.0))
        # Keep anisotropy deformation on each shell so structure is coherent.
        r = base_r + 0.020 * norm
        X = r * np.cos(LAT) * np.cos(LON)
        Y = r * np.cos(LAT) * np.sin(LON)
        Z = r * np.sin(LAT)

        rgba = cm.coolwarm((norm + 1.0) / 2.0)
        # Transparency scale: deeper shells are faint, outer/mid shells clearer.
        alpha = 0.07 + 0.18 * shell_weight + 0.08 * np.abs(norm)
        rgba[..., 3] = np.clip(alpha, 0.05, 0.38)

        ax.plot_surface(
            X,
            Y,
            Z,
            facecolors=rgba,
            rstride=2,
            cstride=2,
            linewidth=0,
            antialiased=False,
            shade=False,
        )

    # Subtle wireframe references for stronger depth perception.
    for ref_r in [0.9, 1.0, 1.1]:
        Xw = ref_r * np.cos(LAT) * np.cos(LON)
        Yw = ref_r * np.cos(LAT) * np.sin(LON)
        Zw = ref_r * np.sin(LAT)
        ax.plot_wireframe(Xw, Yw, Zw, rstride=18, cstride=24, color=(0.08, 0.08, 0.08, 0.11), linewidth=0.4)

    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=38)
    ax.set_axis_off()
    ax.set_title("3D CMB Render — Transparent Tomographic Slices")

    mappable = matplotlib.cm.ScalarMappable(cmap="coolwarm")
    mappable.set_array(norm)
    cbar = fig.colorbar(mappable, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("delta T")

    fig.tight_layout()
    fig.savefig(out / "cmb_holyshit_3d_render.png", dpi=220)


def main():
    out = pathlib.Path("/tmp/bh_renders")
    out.mkdir(parents=True, exist_ok=True)

    root = pathlib.Path("/mnt/riffcastle/castle/garage/grand-2026")
    class_bin = "/tmp/class_public/class"

    tau = read_tau_derived(out / "cmb_tau_derived_report.json")
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

    class_out = run_class(class_bin, params)
    pred_tt = parse_class_dl(class_out, 2)
    pred_ee = parse_class_dl(class_out, 3)
    pred_te = parse_class_dl(class_out, 5)

    obs_tt_b = load_planck(root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt")
    obs_tt_f = load_planck(root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt")
    obs_te_b = load_planck(root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-binned_R3.02.txt")
    obs_te_f = load_planck(root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-full_R3.01.txt")
    obs_ee_b = load_planck(root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-binned_R3.02.txt")
    obs_ee_f = load_planck(root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-full_R3.01.txt")

    chi2_tt_f, red_tt_f, pulls_tt = channel_stats(pred_tt, obs_tt_f)
    chi2_te_f, red_te_f, pulls_te = channel_stats(pred_te, obs_te_f)
    chi2_ee_f, red_ee_f, pulls_ee = channel_stats(pred_ee, obs_ee_f)

    _, red_tt_b, pulls_tt_b = channel_stats(pred_tt, obs_tt_b)
    _, red_te_b, pulls_te_b = channel_stats(pred_te, obs_te_b)
    _, red_ee_b, pulls_ee_b = channel_stats(pred_ee, obs_ee_b)

    # Use binned pulls for clearer plot readability
    plot_spectra_triptych(out, pred_tt, pred_te, pred_ee, obs_tt_b, obs_te_b, obs_ee_b, obs_tt_f, obs_te_f, obs_ee_f)
    plot_pull_triptych(out, pulls_tt_b, pulls_te_b, pulls_ee_b)

    cl_tt = cl_from_dl(pred_tt, lmax=64)
    sky, theta, phi = synth_map_from_cl(cl_tt, lmax=36, nlat=180, nlon=360, seed=7)
    plot_projection_gallery(out, sky, theta, phi)
    plot_3d_render(out, sky, theta, phi)

    summary = {
        "tau_reio": tau,
        "class_output": str(class_out),
        "full_reduced_chi2": {"TT": red_tt_f, "TE": red_te_f, "EE": red_ee_f},
        "binned_reduced_chi2": {"TT": red_tt_b, "TE": red_te_b, "EE": red_ee_b},
        "full_chi2": {"TT": chi2_tt_f, "TE": chi2_te_f, "EE": chi2_ee_f},
        "artifacts": [
            str(out / "cmb_holyshit_spectra_triptych.png"),
            str(out / "cmb_holyshit_pull_triptych.png"),
            str(out / "cmb_holyshit_projection_gallery.png"),
            str(out / "cmb_holyshit_3d_render.png"),
        ],
    }
    (out / "cmb_holyshit_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("wrote", out / "cmb_holyshit_spectra_triptych.png")
    print("wrote", out / "cmb_holyshit_pull_triptych.png")
    print("wrote", out / "cmb_holyshit_projection_gallery.png")
    print("wrote", out / "cmb_holyshit_3d_render.png")
    print("wrote", out / "cmb_holyshit_summary.json")
    print("full reduced chi2:", summary["full_reduced_chi2"])


if __name__ == "__main__":
    main()
