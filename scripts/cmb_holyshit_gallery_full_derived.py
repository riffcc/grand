#!/usr/bin/env python3
from __future__ import annotations

import json
import pathlib
import subprocess
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def run_class(class_bin: str, params: dict[str, float], lmax: int = 2500) -> pathlib.Path:
    td = tempfile.mkdtemp(prefix="cmb_full_derived_gallery_")
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
        ell = int(float(f[0]))
        val = float(f[col_idx_1based - 1])
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
        pts.append((ell, dl, sigma))
    pts.sort(key=lambda p: p[0])
    return pts


def interp(curve, ell: int):
    if ell < curve[0][0] or ell > curve[-1][0]:
        return None
    for (xa, ya), (xb, yb) in zip(curve, curve[1:]):
        if xa <= ell <= xb:
            t = (ell - xa) / (xb - xa) if xb > xa else 0.0
            return ya * (1.0 - t) + yb * t
    return curve[-1][1]


def pulls(pred, obs):
    out = []
    for ell, dl, sig in obs:
        y = interp(pred, ell)
        if y is not None and sig > 0:
            out.append((ell, (y - dl) / sig))
    return out


def main():
    root = pathlib.Path("/mnt/riffcastle/castle/garage/grand-2026")
    out = pathlib.Path("/tmp/bh_renders/cmb_full_derived")
    out.mkdir(parents=True, exist_ok=True)
    class_bin = "/tmp/class_public/class"

    report_path = out / "cmb_full_derived_report.json"
    rep = json.loads(report_path.read_text(encoding="utf-8"))
    inp = rep["inputs"]
    params = {
        "h": float(inp["h"]),
        "omega_b": float(inp["omega_b"]),
        "omega_cdm": float(inp["omega_cdm"]),
        "Omega_k": 0.0,
        "Omega_Lambda": float(1.0 - (inp["omega_b"] + inp["omega_cdm"]) / (inp["h"] * inp["h"]) - 9.0e-5),
        "A_s": float(inp["A_s"]),
        "n_s": float(inp["n_s"]),
        "tau_reio": float(inp["tau_reio"]),
    }

    class_out = run_class(class_bin, params)
    pred_ee = parse_class_dl(class_out, 3)
    obs_ee_b = load_planck(root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-binned_R3.02.txt")
    obs_ee_f = load_planck(root / "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-full_R3.01.txt")
    p_b = pulls(pred_ee, obs_ee_b)
    p_f = pulls(pred_ee, obs_ee_f)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
    ax1.scatter([x for x, y, _ in obs_ee_f], [y for _, y, _ in obs_ee_f], s=4, alpha=0.18, color="#7f8c8d", label="Planck EE full")
    ax1.errorbar(
        [x for x, y, _ in obs_ee_b],
        [y for _, y, _ in obs_ee_b],
        yerr=[s for _, _, s in obs_ee_b],
        fmt="o",
        markersize=3,
        linewidth=0.7,
        color="#f39c12",
        ecolor="#f39c12",
        alpha=0.9,
        label="Planck EE binned",
    )
    ax1.plot([e for e, _ in pred_ee], [v for _, v in pred_ee], lw=2.0, color="#2e86de", label="GUTOE EE full-derived")
    ax1.set_ylabel(r"$D_\ell^{EE}\ [\mu K^2]$")
    ax1.set_ylim(0, 65)
    ax1.grid(alpha=0.2)
    ax1.legend(loc="upper right", fontsize=8)
    ax1.set_title("CMB EE — Full-derived lane (delta=5/2, C_inf=1+1/66)")

    ax2.axhline(0, color="black", lw=1)
    ax2.plot([e for e, p in p_b], [p for _, p in p_b], color="#2e86de", lw=1.2, label="binned pull")
    ax2.plot([e for e, p in p_f], [p for _, p in p_f], color="#95a5a6", lw=0.8, alpha=0.6, label="full pull")
    ax2.set_xlabel(r"Multipole $\ell$")
    ax2.set_ylabel("EE pull")
    ax2.grid(alpha=0.2)
    ax2.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    out_png = out / "cmb_full_derived_ee_gallery.png"
    fig.savefig(out_png, dpi=200)

    summary = {
        "report": str(report_path),
        "class_output": str(class_out),
        "artifact": str(out_png),
        "inputs": inp,
    }
    (out / "cmb_full_derived_ee_gallery.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("wrote", out_png)
    print("wrote", out / "cmb_full_derived_ee_gallery.json")


if __name__ == "__main__":
    main()
