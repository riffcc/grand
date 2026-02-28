#!/usr/bin/env python3
"""
Material Raymarch Atlas
-----------------------
Physically-inspired SDF raymarch renderer for predicted-stable elements.

Outputs:
- <out-dir>/material_raymarch_atlas.png
- <out-dir>/material_raymarch_atlas.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import time
from dataclasses import dataclass

import numpy as np
from PIL import Image, ImageDraw, ImageFont


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
class ElementMat:
    z: int
    symbol: str
    stable_like: int
    family: str
    crystal: str
    albedo: np.ndarray
    metallic: float
    roughness: float
    ior: float


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


def srgb(hex_rgb: str) -> np.ndarray:
    hex_rgb = hex_rgb.lstrip("#")
    return np.array(
        [int(hex_rgb[0:2], 16), int(hex_rgb[2:4], 16), int(hex_rgb[4:6], 16)],
        dtype=np.float64,
    ) / 255.0


def material_params(z: int, family: str, stable_like: int) -> tuple[np.ndarray, float, float, float]:
    # A few anchor colors for recognizable appearance.
    color_overrides = {
        6: "#5A5E66",   # C (graphite-like)
        13: "#C6CDD6",  # Al
        26: "#8F959A",  # Fe
        28: "#A0A8AF",  # Ni
        29: "#C77E52",  # Cu
        47: "#D2D8DE",  # Ag
        50: "#AEB5BF",  # Sn
        74: "#8A9098",  # W
        78: "#BFC6CF",  # Pt
        79: "#D0B13F",  # Au
        82: "#7F858C",  # Pb
    }
    if z in color_overrides:
        albedo = srgb(color_overrides[z])
    else:
        family_color = {
            "alkali": "#D6CFC3",
            "alkaline": "#D7D5C8",
            "transition": "#AAB4C0",
            "post-transition": "#B9B5AF",
            "metalloid": "#8095A6",
            "nonmetal": "#5F7388",
            "halogen": "#7B8AA9",
            "noble": "#96A7C3",
            "lanthanide": "#BFC2C9",
            "actinide": "#B8B0A3",
        }.get(family, "#AAB4C0")
        albedo = srgb(family_color)
        # subtle deterministic hue drift by Z
        shift = 0.015 * math.sin(0.37 * z)
        albedo = np.clip(albedo + shift, 0.06, 0.95)

    if family in {"transition", "post-transition", "lanthanide", "actinide", "alkali", "alkaline"}:
        metallic = 0.93
        roughness = 0.12 + 0.045 * (stable_like % 5)
        ior = 2.0
    elif family == "metalloid":
        metallic = 0.20
        roughness = 0.15 + 0.03 * (z % 4)
        ior = 2.6
    else:
        metallic = 0.04
        roughness = 0.06 + 0.03 * ((z + stable_like) % 5)
        ior = 1.45 + 0.08 * ((z // 10) % 4)
    return np.clip(albedo, 0.03, 0.98), float(np.clip(metallic, 0.0, 1.0)), float(np.clip(roughness, 0.03, 0.92)), ior


def load_elements(scoreboard_csv: pathlib.Path, max_elements: int | None = None) -> list[ElementMat]:
    out: list[ElementMat] = []
    with scoreboard_csv.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["predicted_has_stable"].strip().lower() != "true":
                continue
            z = int(row["Z"])
            stable_like = int(row["predicted_stable_like_isotopes"])
            symbol = SYMBOLS[z] if z < len(SYMBOLS) else f"E{z}"
            family = family_of_z(z)
            crystal = crystal_proxy(family, z, stable_like)
            albedo, metallic, roughness, ior = material_params(z, family, stable_like)
            out.append(
                ElementMat(
                    z=z,
                    symbol=symbol,
                    stable_like=stable_like,
                    family=family,
                    crystal=crystal,
                    albedo=albedo,
                    metallic=metallic,
                    roughness=roughness,
                    ior=ior,
                )
            )
    out.sort(key=lambda e: e.z)
    if max_elements is not None:
        out = out[: max_elements]
    return out


def normalize(v: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, eps)


def rotate_x(p: np.ndarray, ang: float) -> np.ndarray:
    c = math.cos(ang)
    s = math.sin(ang)
    y = c * p[..., 1] - s * p[..., 2]
    z = s * p[..., 1] + c * p[..., 2]
    out = p.copy()
    out[..., 1] = y
    out[..., 2] = z
    return out


def rotate_y(p: np.ndarray, ang: float) -> np.ndarray:
    c = math.cos(ang)
    s = math.sin(ang)
    x = c * p[..., 0] + s * p[..., 2]
    z = -s * p[..., 0] + c * p[..., 2]
    out = p.copy()
    out[..., 0] = x
    out[..., 2] = z
    return out


def sdf_sphere(p: np.ndarray, r: float) -> np.ndarray:
    return np.linalg.norm(p, axis=-1) - r


def sdf_box(p: np.ndarray, b: np.ndarray) -> np.ndarray:
    q = np.abs(p) - b
    outside = np.linalg.norm(np.maximum(q, 0.0), axis=-1)
    inside = np.minimum(np.max(q, axis=-1), 0.0)
    return outside + inside


def sdf_octahedron(p: np.ndarray, s: float) -> np.ndarray:
    return (np.abs(p[..., 0]) + np.abs(p[..., 1]) + np.abs(p[..., 2]) - s) * 0.57735026919


def sdf_hex_prism(p: np.ndarray, h: np.ndarray) -> np.ndarray:
    # https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
    k = np.array([-0.8660254, 0.5, 0.57735], dtype=np.float64)
    p = np.abs(p)
    d = np.maximum(p[..., 2] - h[1], np.maximum(p[..., 0] * k[0] + p[..., 1] * k[1], p[..., 1]) - h[0])
    return d


def sdf_crystal(p: np.ndarray, crystal: str) -> np.ndarray:
    if crystal == "bcc":
        d1 = sdf_octahedron(p, 0.55)
        d2 = sdf_box(p, np.array([0.44, 0.44, 0.44]))
        return np.maximum(d1, d2)  # faceted intersection
    if crystal == "fcc":
        d1 = sdf_box(p, np.array([0.48, 0.38, 0.48]))
        d2 = sdf_octahedron(p, 0.67)
        return np.maximum(d1, d2)
    if crystal == "hcp":
        return sdf_hex_prism(p, np.array([0.36, 0.40]))
    if crystal == "diamond":
        # tetrahedral-ish faceted form
        p2 = rotate_y(p, math.radians(35))
        d1 = sdf_octahedron(p2, 0.48)
        d2 = sdf_box(p2, np.array([0.34, 0.50, 0.34]))
        return np.maximum(d1, d2)
    # molecular cluster
    centers = np.array(
        [[-0.18, 0.02, -0.12], [0.18, 0.02, -0.12], [-0.07, 0.16, 0.17], [0.10, -0.15, 0.13]],
        dtype=np.float64,
    )
    d = np.full(p.shape[:-1], 1e9, dtype=np.float64)
    for c in centers:
        d = np.minimum(d, sdf_sphere(p - c, 0.20))
    return d


def scene_sdf(p: np.ndarray, crystal: str, ang_x: float, ang_y: float) -> tuple[np.ndarray, np.ndarray]:
    q = rotate_y(rotate_x(p, ang_x), ang_y)
    d_obj = sdf_crystal(q, crystal)
    d_floor = p[..., 1] + 0.54
    obj = d_obj < d_floor
    d = np.where(obj, d_obj, d_floor)
    mid = np.where(obj, 1, 2)  # 1=material, 2=floor
    return d, mid


def estimate_normal(pos: np.ndarray, crystal: str, ang_x: float, ang_y: float) -> np.ndarray:
    e = 1.2e-3
    ex = np.array([e, 0.0, 0.0], dtype=np.float64)
    ey = np.array([0.0, e, 0.0], dtype=np.float64)
    ez = np.array([0.0, 0.0, e], dtype=np.float64)
    dx = scene_sdf(pos + ex, crystal, ang_x, ang_y)[0] - scene_sdf(pos - ex, crystal, ang_x, ang_y)[0]
    dy = scene_sdf(pos + ey, crystal, ang_x, ang_y)[0] - scene_sdf(pos - ey, crystal, ang_x, ang_y)[0]
    dz = scene_sdf(pos + ez, crystal, ang_x, ang_y)[0] - scene_sdf(pos - ez, crystal, ang_x, ang_y)[0]
    n = np.stack([dx, dy, dz], axis=-1)
    return normalize(n)


def env_color(rd: np.ndarray) -> np.ndarray:
    t = np.clip(0.5 * (rd[..., 1] + 1.0), 0.0, 1.0)
    sky = np.array([0.12, 0.17, 0.26], dtype=np.float64)
    horizon = np.array([0.28, 0.26, 0.30], dtype=np.float64)
    return sky[None, :] * (t[..., None]) + horizon[None, :] * (1.0 - t[..., None])


def fresnel_schlick(cos_theta: np.ndarray, f0: np.ndarray) -> np.ndarray:
    return f0 + (1.0 - f0) * ((1.0 - cos_theta[..., None]) ** 5)


def render_material_tile(elem: ElementMat, size: int, max_steps: int) -> np.ndarray:
    w = size
    h = size
    xs = (np.arange(w, dtype=np.float64) + 0.5) / w * 2.0 - 1.0
    ys = (np.arange(h, dtype=np.float64) + 0.5) / h * 2.0 - 1.0
    xx, yy = np.meshgrid(xs, ys)

    ro = np.array([0.0, 0.07, 2.35], dtype=np.float64)
    rd = np.stack([xx * 0.95, -yy * 0.95, -np.ones_like(xx)], axis=-1)
    rd = normalize(rd)

    t = np.zeros((h, w), dtype=np.float64)
    active = np.ones((h, w), dtype=bool)
    hit = np.zeros((h, w), dtype=bool)
    mid = np.zeros((h, w), dtype=np.int32)

    # deterministic orientation from Z
    ang_x = math.radians(11.0 + (elem.z % 13))
    ang_y = math.radians(22.0 + (elem.z % 29))

    for _ in range(max_steps):
        if not active.any():
            break
        p = ro[None, None, :] + rd * t[..., None]
        d, m = scene_sdf(p, elem.crystal, ang_x, ang_y)
        just_hit = (d < 9e-4) & active
        hit |= just_hit
        mid = np.where(just_hit, m, mid)
        step = np.clip(d, 1e-4, 0.5)
        t = t + np.where(active, step, 0.0)
        active = (~hit) & (t < 8.0)

    img = env_color(rd.reshape(-1, 3)).reshape(h, w, 3)
    if not hit.any():
        return np.clip(img, 0, 1)

    hp = ro[None, None, :] + rd * t[..., None]
    hit_idx = np.where(hit)
    p_hit = hp[hit_idx]
    rd_hit = rd[hit_idx]
    mid_hit = mid[hit_idx]

    n = estimate_normal(p_hit, elem.crystal, ang_x, ang_y)
    v = normalize(-rd_hit)
    l = normalize(np.array([0.42, 0.86, 0.33], dtype=np.float64)[None, :]).repeat(len(p_hit), axis=0)
    hvec = normalize(v + l)

    ndotl = np.clip(np.sum(n * l, axis=-1), 0.0, 1.0)
    ndotv = np.clip(np.sum(n * v, axis=-1), 0.0, 1.0)
    ndoth = np.clip(np.sum(n * hvec, axis=-1), 0.0, 1.0)
    vdoth = np.clip(np.sum(v * hvec, axis=-1), 0.0, 1.0)

    # Material BRDF params
    albedo = elem.albedo[None, :].repeat(len(p_hit), axis=0)
    metallic = float(elem.metallic)
    roughness = float(elem.roughness)
    a = max(0.035, roughness * roughness)
    a2 = a * a

    denom = (ndoth * ndoth * (a2 - 1.0) + 1.0)
    D = a2 / np.maximum(np.pi * denom * denom, 1e-7)
    k = (roughness + 1.0) ** 2 / 8.0
    Gv = ndotv / np.maximum(ndotv * (1.0 - k) + k, 1e-7)
    Gl = ndotl / np.maximum(ndotl * (1.0 - k) + k, 1e-7)
    G = Gv * Gl

    f0_dielectric = ((elem.ior - 1.0) / (elem.ior + 1.0)) ** 2
    f0 = (1.0 - metallic) * f0_dielectric + metallic * albedo
    F = fresnel_schlick(vdoth, np.asarray(f0))
    spec = (D * G)[:, None] * F / np.maximum(4.0 * ndotv * ndotl, 1e-6)[:, None]

    diffuse = (1.0 - metallic) * albedo / np.pi
    refl = v - 2.0 * np.sum(v * n, axis=-1, keepdims=True) * n
    env_refl = env_color(refl)
    ambient = 0.18 * albedo + 0.22 * env_refl * (1.0 - roughness)
    direct = (diffuse + spec) * ndotl[:, None] * np.array([1.5, 1.45, 1.35], dtype=np.float64)[None, :]

    # Floor shading override
    floor_mask = (mid_hit == 2)
    if floor_mask.any():
        p_f = p_hit[floor_mask]
        checker = (((np.floor((p_f[:, 0] + 4.0) * 4.0) + np.floor((p_f[:, 2] + 4.0) * 4.0)) % 2.0) * 0.06) + 0.15
        floor_col = np.stack([checker, checker, checker + 0.01], axis=-1)
        ambient[floor_mask] = floor_col * 0.9
        direct[floor_mask] = floor_col * ndotl[floor_mask, None] * 0.6

    col = np.clip(ambient + direct, 0.0, 1.0)
    # simple distance fog
    fog = np.exp(-0.06 * np.linalg.norm(p_hit - ro[None, :], axis=-1))
    sky = env_color(rd_hit)
    col = col * fog[:, None] + sky * (1.0 - fog[:, None])

    # write hit pixels
    img_hit = img.reshape(-1, 3)
    flat_idx = hit_idx[0] * w + hit_idx[1]
    img_hit[flat_idx] = col
    img = img_hit.reshape(h, w, 3)

    # gamma
    img = np.clip(img, 0.0, 1.0) ** (1.0 / 2.2)
    return img


def render_atlas(elements: list[ElementMat], out_png: pathlib.Path, out_json: pathlib.Path, cols: int, tile_size: int, max_steps: int) -> None:
    rows = int(math.ceil(len(elements) / cols))
    pad = 4
    atlas_w = cols * tile_size + (cols + 1) * pad
    atlas_h = rows * tile_size + (rows + 1) * pad + 42
    canvas = Image.new("RGB", (atlas_w, atlas_h), (9, 11, 14))
    draw = ImageDraw.Draw(canvas)
    font_title = ImageFont.load_default()
    font_small = ImageFont.load_default()

    draw.text((10, 8), "Material Raymarch Atlas (Normalized Camera/Light)", fill=(226, 233, 242), font=font_title)

    t0 = time.time()
    metadata = []
    for i, e in enumerate(elements):
        r = i // cols
        c = i % cols
        x0 = pad + c * (tile_size + pad)
        y0 = 42 + pad + r * (tile_size + pad)
        tile = render_material_tile(e, size=tile_size, max_steps=max_steps)
        tile_u8 = (np.clip(tile, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
        im = Image.fromarray(tile_u8, mode="RGB")
        canvas.paste(im, (x0, y0))

        # overlay label band
        draw.rectangle((x0 + 3, y0 + 3, x0 + 66, y0 + 28), fill=(0, 0, 0, 120), outline=(100, 108, 120))
        draw.text((x0 + 8, y0 + 7), f"{e.symbol}  Z={e.z}", fill=(240, 244, 250), font=font_small)
        draw.text((x0 + 8, y0 + tile_size - 18), f"{e.crystal} | {e.family}", fill=(220, 226, 238), font=font_small)

        metadata.append(
            {
                "z": e.z,
                "symbol": e.symbol,
                "stable_like_isotopes": e.stable_like,
                "family": e.family,
                "crystal": e.crystal,
                "metallic": e.metallic,
                "roughness": e.roughness,
                "ior": e.ior,
                "albedo": [float(x) for x in e.albedo],
            }
        )

    dt = time.time() - t0
    draw.text((atlas_w - 230, 8), f"tiles={len(elements)}  time={dt:.1f}s", fill=(180, 189, 204), font=font_small)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_png)

    payload = {
        "elements_rendered": len(elements),
        "cols": cols,
        "rows": rows,
        "tile_size": tile_size,
        "raymarch_steps": max_steps,
        "render_seconds": dt,
        "output_png": str(out_png),
        "elements": metadata,
    }
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Render physically-inspired material atlas via SDF raymarching.")
    ap.add_argument("--scoreboard", default="/tmp/nuclear_chart/periodic_table_scoreboard.csv")
    ap.add_argument("--out-dir", default="/tmp/bh_renders/material_atlas_raymarch")
    ap.add_argument("--cols", type=int, default=10)
    ap.add_argument("--tile-size", type=int, default=132)
    ap.add_argument("--steps", type=int, default=70)
    ap.add_argument("--max-elements", type=int, default=0, help="0 means all predicted stable")
    args = ap.parse_args()

    scoreboard = pathlib.Path(args.scoreboard)
    if not scoreboard.exists():
        raise SystemExit(f"missing scoreboard: {scoreboard}")

    max_e = None if args.max_elements <= 0 else args.max_elements
    elements = load_elements(scoreboard, max_elements=max_e)
    if not elements:
        raise SystemExit("no predicted stable elements")

    out_dir = pathlib.Path(args.out_dir)
    out_png = out_dir / "material_raymarch_atlas.png"
    out_json = out_dir / "material_raymarch_atlas.json"
    render_atlas(elements, out_png, out_json, cols=max(1, args.cols), tile_size=max(80, args.tile_size), max_steps=max(24, args.steps))
    print("wrote", out_png)
    print("wrote", out_json)


if __name__ == "__main__":
    main()

