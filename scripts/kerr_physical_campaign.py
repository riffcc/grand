#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


def run_render(repo: Path, env: dict[str, str], size: str, timeout_s: int) -> None:
    cmd = [str(repo / "target" / "release" / "bh_render"), "m87star", size]
    try:
        subprocess.run(cmd, cwd=repo, env=env, timeout=timeout_s, check=True)
    except subprocess.TimeoutExpired:
        # bh_render may linger in gallery mode; image is usually already written.
        pass


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", default="256x256")
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/kerr_campaign"))
    ap.add_argument("--incs", default="17,45,85")
    ap.add_argument("--fovs", default="18,36,60")
    ap.add_argument("--spectra", default="radio,millimeter,optical")
    ap.add_argument("--force-cpu", action="store_true")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    incs = [x.strip() for x in args.incs.split(",") if x.strip()]
    fovs = [x.strip() for x in args.fovs.split(",") if x.strip()]
    spectra = [x.strip() for x in args.spectra.split(",") if x.strip()]

    records = []
    for inc in incs:
        for fov in fovs:
            for spec in spectra:
                slug = f"m87_i{inc}_f{fov}_s{spec}"
                env = os.environ.copy()
                env.update(
                    {
                        "BH_KERR_ASTAR": "0.9",
                        "BH_USE_TRANSFER": "1",
                        "BH_TAU_SCALE": "0.2",
                        "BH_SPECTRUM": spec,
                        "BH_FOV_OVERRIDE": fov,
                        "BH_INC_OVERRIDE": inc,
                        "BH_AZ_OVERRIDE": "0",
                        "BH_MAX_PHI_PI_OVERRIDE": "80",
                        "BH_DPHI_OVERRIDE": "0.003",
                    }
                )
                if args.force_cpu:
                    env["BH_FORCE_CPU"] = "1"

                run_render(repo, env, args.size, args.timeout)
                src_png = Path("/tmp/bh_renders/m87star.png")
                src_json = Path("/tmp/bh_renders/m87star.json")
                dst_png = out / f"{slug}.png"
                dst_json = out / f"{slug}.json"
                if src_png.exists():
                    shutil.copyfile(src_png, dst_png)
                if src_json.exists():
                    shutil.copyfile(src_json, dst_json)
                records.append(
                    {
                        "slug": slug,
                        "inc": float(inc),
                        "fov": float(fov),
                        "spectrum": spec,
                        "png": str(dst_png),
                        "json": str(dst_json),
                    }
                )
                print(f"saved {slug}")

    (out / "campaign.json").write_text(json.dumps(records, indent=2) + "\n")
    print(f"wrote {out / 'campaign.json'} ({len(records)} frames)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
