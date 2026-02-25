#!/usr/bin/env python3
"""
Plot AME2020 residual structure from `ame2020_residuals.csv`.

Input CSV columns:
  Z,N,A,pred_binding_mev,obs_binding_mev,obs_unc_mev,residual_mev,abs_residual_mev
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from statistics import pstdev

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot AME2020 residual structure")
    p.add_argument(
        "--csv",
        default="/tmp/nuclear_chart/ame2020_residuals.csv",
        help="Path to AME residual CSV",
    )
    p.add_argument(
        "--out",
        default="/tmp/nuclear_chart/ame2020_residual_structure.png",
        help="Output PNG path",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    rows = []
    with open(args.csv, newline="") as f:
        for r in csv.DictReader(f):
            rows.append(
                {
                    "Z": int(r["Z"]),
                    "N": int(r["N"]),
                    "A": int(r["A"]),
                    "residual_mev": float(r["residual_mev"]),
                }
            )

    by_n = defaultdict(list)
    by_a = defaultdict(list)
    for row in rows:
        by_n[row["N"]].append(row["residual_mev"])
        by_a[row["A"]].append(row["residual_mev"])

    ns = sorted(by_n.keys())
    as_ = sorted(by_a.keys())
    n_mean = [sum(by_n[n]) / len(by_n[n]) for n in ns]
    n_std = [pstdev(by_n[n]) if len(by_n[n]) > 1 else 0.0 for n in ns]
    a_mean = [sum(by_a[a]) / len(by_a[a]) for a in as_]
    a_std = [pstdev(by_a[a]) if len(by_a[a]) > 1 else 0.0 for a in as_]

    fig, ax = plt.subplots(2, 1, figsize=(12, 9))
    ax[0].plot(ns, n_mean, color="tab:blue", lw=1.5, label="mean residual")
    ax[0].fill_between(
        ns,
        [m - s for m, s in zip(n_mean, n_std)],
        [m + s for m, s in zip(n_mean, n_std)],
        color="tab:blue",
        alpha=0.25,
        label="±1σ across Z",
    )
    for magic_n in (50, 82, 126, 184):
        ax[0].axvline(magic_n, color="gray", ls="--", lw=0.8)
    ax[0].axhline(0.0, color="black", lw=0.8)
    ax[0].set_title("AME2020 Residual Structure vs Neutron Number (N)")
    ax[0].set_ylabel("Residual (MeV): pred - obs")
    ax[0].legend(loc="upper right", fontsize=9)

    ax[1].plot(as_, a_mean, color="tab:orange", lw=1.5, label="mean residual")
    ax[1].fill_between(
        as_,
        [m - s for m, s in zip(a_mean, a_std)],
        [m + s for m, s in zip(a_mean, a_std)],
        color="tab:orange",
        alpha=0.25,
        label="±1σ across isobars",
    )
    ax[1].axhline(0.0, color="black", lw=0.8)
    ax[1].set_title("AME2020 Residual Trend vs Mass Number (A)")
    ax[1].set_xlabel("A")
    ax[1].set_ylabel("Residual (MeV): pred - obs")
    ax[1].legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

