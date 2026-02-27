#!/usr/bin/env python3
"""
CMB pull-structure diagnostics for GRAND-355.

Produces:
- /tmp/bh_renders/cmb_vs_planck.png
- /tmp/bh_renders/cmb_vs_planck_zoom.png
- /tmp/bh_renders/cmb_pull_ratio_denoise.png
- /tmp/bh_renders/cmb_pull_ratio_denoise.json
- /tmp/bh_renders/cmb_full_delta_denoise.png
- /tmp/bh_renders/cmb_full_delta_denoise.json

This script is diagnostic-only: it does not alter model parameters.
"""

from __future__ import annotations

import json
import pathlib
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import matplotlib
import numpy as np
from scipy.signal import savgol_filter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


@dataclass(frozen=True)
class DlPoint:
    ell: int
    dl: float
    sigma: float


def run_class(class_bin: str, params: Dict[str, float]) -> List[Tuple[int, float]]:
    with tempfile.TemporaryDirectory(prefix="cmb_pull_") as tmp:
        td = pathlib.Path(tmp)
        ini = td / "run.ini"
        root = td / "g_"
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
                    "l_max_scalars=2500",
                    "format=camb",
                    f"root = {root}",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        subprocess.run([class_bin, str(ini)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        cands = sorted(td.glob("*cl*.dat"))
        cands.sort(key=lambda p: (0 if "lensedcls" in p.name.lower() else 1))
        out: List[Tuple[int, float]] = []
        for line in cands[0].read_text(encoding="utf-8").splitlines():
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            f = s.split()
            if len(f) < 2:
                continue
            try:
                ell = int(float(f[0]))
                dl = float(f[1])
            except ValueError:
                continue
            if 2 <= ell <= 2500:
                out.append((ell, dl))
        out.sort()
        return out


def interp(curve: Sequence[Tuple[int, float]], ell: int) -> float | None:
    if ell < curve[0][0] or ell > curve[-1][0]:
        return None
    for (xa, ya), (xb, yb) in zip(curve, curve[1:]):
        if xa <= ell <= xb:
            t = (ell - xa) / (xb - xa) if xb > xa else 0.0
            return ya * (1.0 - t) + yb * t
    return curve[-1][1]


def load_planck(path: pathlib.Path) -> List[DlPoint]:
    out: List[DlPoint] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        if len(f) < 4:
            continue
        ell = int(round(float(f[0])))
        if not (2 <= ell <= 2500):
            continue
        dl = float(f[1])
        sigma = 0.5 * (abs(float(f[2])) + abs(float(f[3])))
        out.append(DlPoint(ell=ell, dl=dl, sigma=sigma))
    out.sort(key=lambda p: p.ell)
    return out


def main() -> None:
    class_bin = "/tmp/class_public/class"
    root = pathlib.Path("/mnt/riffcastle/castle/garage/grand-2026")
    planck_binned = root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt"
    planck_full = root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt"
    out_dir = pathlib.Path("/tmp/bh_renders")
    out_dir.mkdir(parents=True, exist_ok=True)

    baseline = {
        "h": 0.680163311753,
        "omega_b": 0.022807271041,
        "omega_cdm": 0.124403296589,
        "Omega_k": 0.0,
        "Omega_Lambda": 0.681700909091,
        "n_s": 0.965416666667,
        "A_s": 2.219284522994e-9,
        "tau_reio": 0.054,
    }
    tuned = dict(baseline)
    tuned["A_s"] = 2.116643e-9
    tuned["tau_reio"] = 0.051

    cb = run_class(class_bin, baseline)
    ct = run_class(class_bin, tuned)
    obs_b = load_planck(planck_binned)
    obs_f = load_planck(planck_full)

    # Plot: ours vs theirs
    ells_base = [x for x, _ in cb]
    dl_base = [y for _, y in cb]
    ells_tuned = [x for x, _ in ct]
    dl_tuned = [y for _, y in ct]
    ells_full = [p.ell for p in obs_f]
    dl_full = [p.dl for p in obs_f]
    ells_bin = [p.ell for p in obs_b]
    dl_bin = [p.dl for p in obs_b]
    sig_bin = [p.sigma for p in obs_b]

    pull_base = np.array([(interp(cb, p.ell) - p.dl) / p.sigma for p in obs_b], dtype=float)
    pull_tuned = np.array([(interp(ct, p.ell) - p.dl) / p.sigma for p in obs_b], dtype=float)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax1.scatter(ells_full, dl_full, s=5, alpha=0.18, color="#7f8c8d", label="Planck TT full (unbinned)")
    ax1.errorbar(
        ells_bin,
        dl_bin,
        yerr=sig_bin,
        fmt="o",
        markersize=3.5,
        linewidth=0.8,
        color="#f39c12",
        ecolor="#f39c12",
        alpha=0.9,
        label="Planck TT binned",
    )
    ax1.plot(ells_base, dl_base, color="#2e86de", lw=2.0, label="GUTOE CLASS baseline")
    ax1.plot(ells_tuned, dl_tuned, color="#27ae60", lw=2.0, ls="--", label="GUTOE CLASS tuned A_s/tau")
    ax1.set_ylabel(r"$D_\ell^{TT}\ [\mu K^2]$")
    ax1.set_xlim(2, 2500)
    ax1.set_title("CMB TT: GUTOE prediction vs Planck (full + binned)")
    ax1.grid(alpha=0.2)
    ax1.legend(loc="upper right", fontsize=9)

    ax2.axhline(0.0, color="black", lw=1)
    ax2.plot(ells_bin, pull_base, color="#2e86de", lw=1.5, label="baseline pull vs binned")
    ax2.plot(ells_bin, pull_tuned, color="#27ae60", lw=1.5, ls="--", label="tuned pull vs binned")
    ax2.set_xlabel(r"Multipole $\ell$")
    ax2.set_ylabel("pull")
    ax2.set_ylim(-15, 15)
    ax2.grid(alpha=0.2)
    ax2.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "cmb_vs_planck.png", dpi=180)

    fig2, ax = plt.subplots(figsize=(13, 5.5))
    ax.errorbar(
        ells_bin,
        dl_bin,
        yerr=sig_bin,
        fmt="o",
        markersize=3.5,
        linewidth=0.8,
        color="#f39c12",
        ecolor="#f39c12",
        alpha=0.9,
        label="Planck binned",
    )
    ax.plot(ells_base, dl_base, color="#2e86de", lw=2.0, label="Baseline")
    ax.plot(ells_tuned, dl_tuned, color="#27ae60", lw=2.0, ls="--", label="Tuned A_s/tau")
    ax.set_xlim(30, 2500)
    ax.set_ylim(0, 6500)
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"$D_\ell^{TT}\ [\mu K^2]$")
    ax.set_title("CMB TT peak structure and damping tail")
    ax.grid(alpha=0.2)
    ax.legend(loc="upper right")
    fig2.tight_layout()
    fig2.savefig(out_dir / "cmb_vs_planck_zoom.png", dpi=180)

    # Binned pull transform fit
    ells_b = np.array(ells_bin, dtype=float)
    pb = pull_base
    pt = pull_tuned
    A = np.column_stack([pb, np.ones_like(pb)])
    a_glob, b_glob = np.linalg.lstsq(A, pt, rcond=None)[0]
    pt_glob = a_glob * pb + b_glob
    ss_res = np.sum((pt - pt_glob) ** 2)
    ss_tot = np.sum((pt - pt.mean()) ** 2)
    r2_glob = 1.0 - ss_res / ss_tot

    w = 9
    a_loc = np.zeros_like(pb)
    b_loc = np.zeros_like(pb)
    for i in range(len(pb)):
        lo = max(0, i - w)
        hi = min(len(pb), i + w + 1)
        X = np.column_stack([pb[lo:hi], np.ones(hi - lo)])
        y = pt[lo:hi]
        aa, bb = np.linalg.lstsq(X, y, rcond=None)[0]
        a_loc[i] = aa
        b_loc[i] = bb
    pt_hat = a_loc * pb + b_loc
    rmse_before = float(np.sqrt(np.mean((pt - pb) ** 2)))
    rmse_after = float(np.sqrt(np.mean((pt - pt_hat) ** 2)))

    fig3, axs = plt.subplots(2, 2, figsize=(13, 8.5))
    ax = axs[0, 0]
    ax.scatter(pb, pt, s=22, alpha=0.75, color="#2e86de", label="binned multipoles")
    xx = np.linspace(pb.min(), pb.max(), 200)
    ax.plot(xx, a_glob * xx + b_glob, color="#e74c3c", lw=2, label=f"global fit: y={a_glob:.3f}x+{b_glob:.3f}")
    ax.set_xlabel("baseline pull")
    ax.set_ylabel("tuned pull")
    ax.set_title(f"Scatter fit (R²={r2_glob:.3f})")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axs[0, 1]
    ax.plot(ells_b, a_loc, color="#8e44ad", lw=1.8, label="a(ℓ) local")
    ax.plot(ells_b, b_loc, color="#16a085", lw=1.4, label="b(ℓ) local")
    ax.axhline(1, color="gray", ls="--", lw=1)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("ℓ")
    ax.set_ylabel("coeff")
    ax.set_title("Local affine coefficients: p_t ≈ a(ℓ) p_b + b(ℓ)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axs[1, 0]
    ax.plot(ells_b, pb, color="#2e86de", lw=1.2, label="baseline pull")
    ax.plot(ells_b, pt, color="#27ae60", lw=1.2, label="tuned pull")
    ax.plot(ells_b, pt_hat, color="#d35400", lw=1.4, ls="--", label="affine denoised fit")
    ax.set_xlabel("ℓ")
    ax.set_ylabel("pull")
    ax.set_title("Pull transform and denoised model")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axs[1, 1]
    ax.plot(ells_b, pt - pb, color="#7f8c8d", lw=1.2, label="Δ pull (tuned-baseline)")
    ax.plot(ells_b, pt - pt_hat, color="#c0392b", lw=1.2, label="residual after local affine")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("ℓ")
    ax.set_ylabel("pull residual")
    ax.set_title(f"Residual compression RMSE: before={rmse_before:.3f}, after={rmse_after:.3f}")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig3.tight_layout()
    fig3.savefig(out_dir / "cmb_pull_ratio_denoise.png", dpi=180)

    # Full-spectrum smooth component extraction from delta pull
    pull_b_full = np.array([(interp(cb, p.ell) - p.dl) / p.sigma for p in obs_f], dtype=float)
    pull_t_full = np.array([(interp(ct, p.ell) - p.dl) / p.sigma for p in obs_f], dtype=float)
    delta = pull_t_full - pull_b_full
    ell_full = np.array([p.ell for p in obs_f], dtype=float)

    win = 301
    if win >= len(delta):
        win = len(delta) - 1 if len(delta) % 2 == 0 else len(delta)
    if win % 2 == 0:
        win += 1
    smooth = savgol_filter(delta, window_length=win, polyorder=3, mode="interp")
    resid = delta - smooth
    rmse_delta = float(np.sqrt(np.mean(delta**2)))
    rmse_resid = float(np.sqrt(np.mean(resid**2)))
    compress = rmse_delta / rmse_resid if rmse_resid > 0 else float("inf")

    fig4, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13.5, 7.8), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
    )
    ax1.plot(ell_full, delta, color="#34495e", lw=1.0, label="Δ pull (tuned - baseline)")
    ax1.plot(ell_full, smooth, color="#e67e22", lw=2.0, label="smooth structured component")
    ax1.axhline(0, color="black", lw=0.8)
    ax1.set_ylabel("Δ pull")
    ax1.set_title("Full-spectrum structured component in pull delta")
    ax1.grid(alpha=0.2)
    ax1.legend(fontsize=9)

    ax2.plot(ell_full, resid, color="#16a085", lw=1.0, label="residual after subtracting smooth component")
    ax2.axhline(0, color="black", lw=0.8)
    ax2.set_xlabel(r"Multipole $\ell$")
    ax2.set_ylabel("residual")
    ax2.set_title(f"Residual RMSE: raw={rmse_delta:.4f}, after={rmse_resid:.4f} (compression={compress:.2f}x)")
    ax2.grid(alpha=0.2)
    ax2.legend(fontsize=9)
    fig4.tight_layout()
    fig4.savefig(out_dir / "cmb_full_delta_denoise.png", dpi=200)

    (out_dir / "cmb_pull_ratio_denoise.json").write_text(
        json.dumps(
            {
                "n_points": int(len(pb)),
                "global_affine": {"a": float(a_glob), "b": float(b_glob), "r2": float(r2_glob)},
                "rmse_before_vs_baseline": rmse_before,
                "rmse_after_local_affine": rmse_after,
                "compression_factor": float(rmse_before / rmse_after if rmse_after > 0 else np.inf),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (out_dir / "cmb_full_delta_denoise.json").write_text(
        json.dumps(
            {
                "rmse_delta": rmse_delta,
                "rmse_resid": rmse_resid,
                "compression_factor": compress,
                "window": int(win),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print("wrote:")
    for name in [
        "cmb_vs_planck.png",
        "cmb_vs_planck_zoom.png",
        "cmb_pull_ratio_denoise.png",
        "cmb_pull_ratio_denoise.json",
        "cmb_full_delta_denoise.png",
        "cmb_full_delta_denoise.json",
    ]:
        print(" ", out_dir / name)


if __name__ == "__main__":
    main()
