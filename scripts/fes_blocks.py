#!/usr/bin/env python3
"""Reweight OPES COLVAR files into block free energies.

Each file has columns time, collective variable, bias, and optionally
the OPES rct. Pass the walker files in any order:

    fes_blocks.py w0.COLVAR w1.COLVAR ...

The last 20% of the run is split into five equal blocks. A bin needs
at least 10 samples. Each block is shifted so its minimum is zero.
The script writes the table to stdout.
"""
from __future__ import annotations

import math
import sys

import numpy as np

KT = 8.314462618 * 210.0 / 1000.0


def load(paths):
    frames = []
    for path in paths:
        rows = []
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split()
                rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
        frames.append(np.asarray(rows))
    return np.vstack(frames)


def block_fes(sample, bins):
    centers = 0.5 * (bins[:-1] + bins[1:])
    free = np.full(len(centers), np.nan)
    for i in range(len(centers)):
        mask = (sample[:, 1] >= bins[i]) & (sample[:, 1] < bins[i + 1])
        if int(mask.sum()) < 10:
            continue
        bias = sample[mask, 2] / KT
        shift = bias.max()
        logw = shift + math.log(np.exp(np.clip(bias - shift, -80.0, 0.0)).sum())
        free[i] = -KT * logw
    finite = np.isfinite(free)
    if finite.any():
        free = free - np.nanmin(free)
    return centers, free


def main(argv):
    if len(argv) < 2:
        print("usage: fes_blocks.py COLVAR [COLVAR ...]", file=sys.stderr)
        return 2
    data = load(argv[1:])
    t0 = 0.8 * data[:, 0].max()
    late = data[data[:, 0] >= t0]
    edges = np.linspace(t0, data[:, 0].max(), 6)
    bins = np.linspace(0.10, 0.98, 45)
    columns = []
    centers = None
    for left, right in zip(edges[:-1], edges[1:]):
        if right < data[:, 0].max():
            sample = late[(late[:, 0] >= left) & (late[:, 0] < right)]
        else:
            sample = late[late[:, 0] >= left]
        centers, free = block_fes(sample, bins)
        columns.append(free)
    stacked = np.vstack(columns)
    span = np.nanmax(stacked, axis=0) - np.nanmin(stacked, axis=0)
    print("# s  F_b0 F_b1 F_b2 F_b3 F_b4  span")
    print("# Last 20% of the run, five blocks, each zeroed at its minimum.")
    for i, score in enumerate(centers):
        values = " ".join(
            "nan" if not np.isfinite(stacked[k, i]) else f"{stacked[k, i]:.3f}"
            for k in range(stacked.shape[0])
        )
        width = "nan" if not np.isfinite(span[i]) else f"{span[i]:.3f}"
        print(f"{score:.3f} {values} {width}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
