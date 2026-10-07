#!/usr/bin/env python3
"""Generate a deterministic synthetic test case for GGR.

This example exists only to demonstrate and test the open computational
workflow. It is not part of the field dataset and does not reproduce the
Petrel interpolation used for the field results.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import griddata

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import ggr_workflow as ggr  # noqa: E402


def build_synthetic_surface():
    x = np.arange(0.0, 6000.0 + ggr.GRID_CELL_M, ggr.GRID_CELL_M)
    y = np.arange(0.0, 6000.0 + ggr.GRID_CELL_M, ggr.GRID_CELL_M)
    X, Y = np.meshgrid(x, y)

    broad_structure = (
        1450.0
        + 0.020 * X
        + 0.010 * Y
        + 38.0 * np.exp(-((X - 3650.0) ** 2 + (Y - 3100.0) ** 2) / (2.0 * 1050.0 ** 2))
    )

    vertical_lines = np.array([1500.0, 3000.0, 4500.0])
    distance_vertical = np.min(np.abs(X[..., None] - vertical_lines), axis=2)
    distance_horizontal = np.abs(Y - 3000.0)
    distance_to_lines = np.minimum(distance_vertical, distance_horizontal)

    imprint = (
        7.0
        * np.exp(-(distance_to_lines / 140.0) ** 2)
        * np.sin(2.0 * np.pi * (Y / 850.0 + X / 1700.0))
    )

    observed = broad_structure + imprint
    return x, y, X, Y, observed, distance_to_lines


def reconstruct(points, values, X, Y):
    surface = griddata(points, values, (X, Y), method="linear")
    if np.any(~np.isfinite(surface)):
        nearest = griddata(points, values, (X, Y), method="nearest")
        surface = np.where(np.isfinite(surface), surface, nearest)
    return surface


def write_xyz(path: Path, X, Y, Z):
    pd.DataFrame(
        {"X": X.ravel(), "Y": Y.ravel(), "Z": Z.ravel()}
    ).to_csv(path, sep=" ", index=False, header=False, float_format="%.6f")


def run_demo(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    x, y, X, Y, original, distance_grid = build_synthetic_surface()
    original_grid = ggr.SurfaceGrid(x=x, y=y, z=original, name="Synthetic original")
    ggr.validate_grid(original_grid)

    sampled = ggr.sample_surface(original_grid, ggr.SAMPLE_SPACING_M)

    # Distance to the synthetic line network can be evaluated analytically.
    vertical_lines = np.array([1500.0, 3000.0, 4500.0])
    dv = np.min(np.abs(sampled["X"].to_numpy()[:, None] - vertical_lines), axis=1)
    dh = np.abs(sampled["Y"].to_numpy() - 3000.0)
    sampled["distance_to_line_m"] = np.minimum(dv, dh)

    points = sampled[["X", "Y"]].to_numpy(float)
    values = sampled["Z"].to_numpy(float)
    control = reconstruct(points, values, X, Y)

    control_lii = ggr.lii(
        control,
        distance_grid,
        sigma_m=ggr.BASE_SIGMA_M,
        near_m=ggr.BASE_NEAR_M,
    )

    rows = []
    for D in ggr.D_VALUES_M:
        retained = sampled.loc[sampled["distance_to_line_m"] > D]
        surface = reconstruct(
            retained[["X", "Y"]].to_numpy(float),
            retained["Z"].to_numpy(float),
            X,
            Y,
        )

        value = ggr.lii(
            surface,
            distance_grid,
            sigma_m=ggr.BASE_SIGMA_M,
            near_m=ggr.BASE_NEAR_M,
        )
        reduction = 100.0 * (1.0 - value / control_lii)
        _, rmse = ggr.mae_rmse(surface - control)

        rows.append(
            {
                "D_m": int(D),
                "Retained_nodes": int(len(retained)),
                "LII": float(value),
                "LII_Reduction_pct": float(reduction),
                "RMSE_vs_synthetic_control_ms": float(rmse),
            }
        )

    metrics = pd.DataFrame(rows)
    write_xyz(output_dir / "synthetic_original_surface.xyz", X, Y, original)
    sampled.to_csv(output_dir / "synthetic_sampled_nodes.csv", index=False)
    metrics.to_csv(output_dir / "demo_metrics.csv", index=False)

    print("Synthetic demonstration completed.")
    print("These values are for software testing only and are not field-study results.")
    print(metrics.to_string(index=False))
    return metrics


def main():
    run_demo(Path(__file__).resolve().parent / "output")


if __name__ == "__main__":
    main()
