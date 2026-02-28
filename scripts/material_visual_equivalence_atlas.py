#!/usr/bin/env python3
"""
Render a normalized visual-equivalence atlas for predicted stable elements.

This is a first-pass physically motivated proxy renderer:
- Uses predicted stable element set from periodic_table_scoreboard.csv
- Assigns a crystal-family motif from chemistry category heuristics
- Uses normalized lighting/size so all materials are visually comparable

Outputs:
- <out-dir>/material_visual_equivalence_atlas.png
- <out-dir>/material_visual_equivalence_atlas.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle


SYMBOLS = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]


ALKALI = {3, 11, 19, 37, 55, 87}
ALKALINE = {4, 12, 20, 38, 56, 88}
HALOGEN = {9, 17, 35, 53, 85, 117}
NOBLE = {2, 10, 18, 36, 54, 86, 118}
NONMETAL = {1, 6, 7, 8, 15, 16, 34}
METALLOID = {5, 14, 32, 33, 51, 52}
LANTHANIDE = set(range(57, 72))
ACTINIDE = set(range(89, 104))
POST_TRANSITION = {13, 31, 49, 50, 81, 82, 83}


@dataclass
class ElementRow:
    z: int
    symbol: str
    stable_like: int
    family: str
    crystal: str
    color: tuple[float, float, float]


def family_of_z(z: int) -> str:
    if z in ALKALI:
        return "alkali"
    if z in ALKALINE:
        return "alkaline"
    if z in HALOGEN:
        return "halogen"
    if z in NOBLE:
        return "noble"
    if z in NONMETAL:
        return "nonmetal"
    if z in METALLOID:
        return "metalloid"
    if z in LANTHANIDE:
        return "lanthanide"
    if z in ACTINIDE:
        return "actinide"
    if z in POST_TRANSITION:
        return "post-transition"
    return "transition"


def crystal_proxy(family: str, z: int, stable_like: int) -> str:
    if family in {"alkali", "transition"}:
        return ["bcc", "fcc", "hcp"][z % 3]
    if family in {"alkaline", "post-transition"}:
        return ["hcp", "fcc", "bcc"][stable_like % 3]
    if family == "metalloid":
        return "diamond"
    if family in {"nonmetal", "halogen", "noble"}:
        return "molecular"
    if family in {"lanthanide", "actinide"}:
        return "hcp"
    return "bcc"


def color_proxy(family: str, z: int, stable_like: int) -> tuple[float, float, float]:
    # Family base hues chosen for strong visual separation in the atlas.
    family_hue = {
        "alkali": 0.02,
        "alkaline": 0.09,
        "transition": 0.60,
        "post-transition": 0.55,
        "metalloid": 0.35,
        "nonmetal": 0.72,
        "halogen": 0.78,
        "noble": 0.82,
        "lanthanide": 0.14,
        "actinide": 0.00,
    }.get(family, 0.60)

    h = (family_hue + 0.003 * (z % 11)) % 1.0
    s = min(0.95, 0.35 + 0.08 * stable_like)
    v = min(0.98, 0.62 + 0.03 * (z % 7))
    return hsv_to_rgb(h, s, v)


def hsv_to_rgb(h: float, s: float, v: float) -> tuple[float, float, float]:
    i = int(h * 6.0)
    f = h * 6.0 - i
    p = v * (1.0 - s)
    q = v * (1.0 - f * s)
    t = v * (1.0 - (1.0 - f) * s)
    i %= 6
    if i == 0:
        return (v, t, p)
    if i == 1:
        return (q, v, p)
    if i == 2:
        return (p, v, t)
    if i == 3:
        return (p, q, v)
    if i == 4:
        return (t, p, v)
    return (v, p, q)


def draw_lattice(ax: plt.Axes, crystal: str, color: tuple[float, float, float]) -> None:
    c = np.array(color)
    dot = np.clip(c * 0.85 + 0.15, 0, 1)
    if crystal == "bcc":
        pts = [(0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75), (0.5, 0.5)]
    elif crystal == "fcc":
        pts = [
            (0.2, 0.2),
            (0.8, 0.2),
            (0.2, 0.8),
            (0.8, 0.8),
            (0.5, 0.2),
            (0.2, 0.5),
            (0.8, 0.5),
            (0.5, 0.8),
        ]
    elif crystal == "hcp":
        pts = []
        for j in range(4):
            y = 0.18 + j * 0.2
            offset = 0.08 if j % 2 else 0.0
            for i in range(4):
                x = 0.18 + i * 0.2 + offset
                if x < 0.9:
                    pts.append((x, y))
    elif crystal == "diamond":
        pts = [(0.25, 0.25), (0.75, 0.25), (0.5, 0.5), (0.25, 0.75), (0.75, 0.75)]
    else:  # molecular
        pts = [(0.3, 0.3), (0.7, 0.3), (0.5, 0.5), (0.3, 0.7), (0.7, 0.7)]
        for p1, p2 in [((0.3, 0.3), (0.5, 0.5)), ((0.5, 0.5), (0.7, 0.7)), ((0.7, 0.3), (0.5, 0.5)), ((0.5, 0.5), (0.3, 0.7))]:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=dot * 0.8, lw=0.8, alpha=0.5)

    for x, y in pts:
        ax.add_patch(Circle((x, y), 0.018, facecolor=dot, edgecolor="none", alpha=0.55))


def draw_sphere(ax: plt.Axes, color: tuple[float, float, float]) -> None:
    base = np.array(color)
    ax.add_patch(Circle((0.5, 0.54), 0.22, facecolor=base, edgecolor=(0, 0, 0), lw=0.8, alpha=0.95))
    ax.add_patch(Circle((0.43, 0.63), 0.075, facecolor=(1, 1, 1), edgecolor="none", alpha=0.35))
    ax.add_patch(Circle((0.58, 0.45), 0.16, facecolor=(0, 0, 0), edgecolor="none", alpha=0.12))


def load_predicted_stable(scoreboard: pathlib.Path) -> list[ElementRow]:
    rows: list[ElementRow] = []
    with scoreboard.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            z = int(r["Z"])
            pred_has_stable = r["predicted_has_stable"].strip().lower() == "true"
            if not pred_has_stable:
                continue
            stable_like = int(r["predicted_stable_like_isotopes"])
            symbol = SYMBOLS[z] if 0 <= z < len(SYMBOLS) else f"E{z}"
            family = family_of_z(z)
            crystal = crystal_proxy(family, z, stable_like)
            color = color_proxy(family, z, stable_like)
            rows.append(
                ElementRow(
                    z=z,
                    symbol=symbol,
                    stable_like=stable_like,
                    family=family,
                    crystal=crystal,
                    color=color,
                )
            )
    rows.sort(key=lambda e: e.z)
    return rows


def render_atlas(elements: list[ElementRow], out_png: pathlib.Path, cols: int) -> None:
    n = len(elements)
    rows = int(math.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.15, rows * 2.45), dpi=220)
    axes = np.array(axes).reshape(rows, cols)

    bg = np.array([0.06, 0.07, 0.09])
    for i in range(rows * cols):
        r = i // cols
        c = i % cols
        ax = axes[r, c]
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.add_patch(Rectangle((0, 0), 1, 1, facecolor=bg, edgecolor=(0.2, 0.22, 0.27), lw=0.8))
        if i >= n:
            continue

        e = elements[i]
        draw_lattice(ax, e.crystal, e.color)
        draw_sphere(ax, e.color)
        ax.text(0.08, 0.92, f"{e.symbol}", color="white", fontsize=11.5, fontweight="bold", ha="left", va="top")
        ax.text(0.92, 0.92, f"Z={e.z}", color="#D6D9DF", fontsize=7.5, ha="right", va="top")
        ax.text(0.5, 0.17, e.crystal.upper(), color="#E2E7F2", fontsize=7.5, ha="center", va="center")
        ax.text(0.5, 0.08, f"{e.family}", color="#98A4B8", fontsize=6.9, ha="center", va="center")

    fig.suptitle(
        "GUTOE Material Visual Equivalence Atlas (Predicted Stable Elements)\n"
        "Normalized Lighting + Crystal-Class Motif Proxy",
        fontsize=14,
        color="#EAEFF7",
        y=0.995,
    )
    fig.patch.set_facecolor("#0B0D10")
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.96])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Render visual-equivalence atlas for predicted stable elements.")
    ap.add_argument(
        "--scoreboard",
        default="/tmp/nuclear_chart/periodic_table_scoreboard.csv",
        help="Path to periodic_table_scoreboard.csv",
    )
    ap.add_argument("--out-dir", default="/tmp/bh_renders/material_atlas", help="Output directory")
    ap.add_argument("--cols", type=int, default=10, help="Tile columns")
    args = ap.parse_args()

    scoreboard = pathlib.Path(args.scoreboard)
    if not scoreboard.exists():
        raise SystemExit(f"missing scoreboard: {scoreboard}")

    elements = load_predicted_stable(scoreboard)
    if not elements:
        raise SystemExit("no predicted stable elements in scoreboard")

    out_dir = pathlib.Path(args.out_dir)
    out_png = out_dir / "material_visual_equivalence_atlas.png"
    out_json = out_dir / "material_visual_equivalence_atlas.json"

    render_atlas(elements, out_png, cols=max(1, args.cols))

    payload = {
        "scoreboard": str(scoreboard),
        "elements_rendered": len(elements),
        "z_min": elements[0].z,
        "z_max": elements[-1].z,
        "output_png": str(out_png),
        "elements": [
            {
                "z": e.z,
                "symbol": e.symbol,
                "stable_like_isotopes": e.stable_like,
                "family": e.family,
                "crystal_proxy": e.crystal,
                "rgb": [float(e.color[0]), float(e.color[1]), float(e.color[2])],
            }
            for e in elements
        ],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("wrote", out_png)
    print("wrote", out_json)


if __name__ == "__main__":
    main()

