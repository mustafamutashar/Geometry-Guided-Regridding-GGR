#!/usr/bin/env python3
"""Geometry-Guided Regridding (GGR) computational workflow.

This module implements the open computational stages used in the GGR study:
regular surface sampling, distance-based node exclusion, quantitative
surface diagnostics, and sensitivity analysis. Field-surface reconstruction
remains an external Petrel step because the study used one fixed Petrel
gridding configuration for the control and all GGR cases.

Coordinates are expected in metres in a projected CRS. TWT values are in ms.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

try:
    import geopandas as gpd
    import shapely
    from shapely.geometry import Point
except Exception:
    gpd = None
    shapely = None
    Point = None

CODE_VERSION = "1.0.0"
GRID_CELL_M = 50.0
SAMPLE_SPACING_M = 100.0
D_VALUES_M = np.arange(100, 1001, 100, dtype=int)
BASE_SIGMA_M = 250.0
BASE_NEAR_M = 250.0
FAR_MIN_M = 1500.0
FAR_MAX_M = 3000.0
TREND_SIGMA_M = 1000.0
SENSITIVITY_VALUES_M = np.array([150, 200, 250, 300, 400, 500], dtype=float)


@dataclass
class SurfaceGrid:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    name: str = ""

    @property
    def dx(self) -> float:
        return float(np.median(np.diff(self.x)))

    @property
    def dy(self) -> float:
        return float(np.median(np.diff(self.y)))


def _read_xyz(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    rows: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = re.split(r"[,;\s]+", line.strip())
            if len(parts) < 3:
                continue
            try:
                rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No numeric XYZ rows found in {path}")
    return pd.DataFrame(rows, columns=["X", "Y", "Z"])


def load_xyz_grid(path: str | Path, name: str = "") -> SurfaceGrid:
    df = _read_xyz(path)
    if df.duplicated(["X", "Y"]).any():
        raise ValueError(f"Duplicate XY coordinates found in {path}")

    x = np.sort(df["X"].unique())
    y = np.sort(df["Y"].unique())
    xi = pd.Series(np.arange(x.size), index=x)
    yi = pd.Series(np.arange(y.size), index=y)

    z = np.full((y.size, x.size), np.nan, dtype=float)
    ix = df["X"].map(xi).to_numpy(dtype=int)
    iy = df["Y"].map(yi).to_numpy(dtype=int)
    z[iy, ix] = df["Z"].to_numpy(float)

    surface = SurfaceGrid(x=x, y=y, z=z, name=name or Path(path).stem)
    validate_grid(surface)
    return surface


def validate_grid(surface: SurfaceGrid, expected_cell_m: float = GRID_CELL_M) -> None:
    if surface.x.size < 2 or surface.y.size < 2:
        raise ValueError(f"{surface.name}: grid must contain at least two X and Y coordinates")
    dx = np.diff(surface.x)
    dy = np.diff(surface.y)
    if np.any(dx <= 0) or np.any(dy <= 0):
        raise ValueError(f"{surface.name}: X and Y coordinates must increase monotonically")
    if not np.allclose(dx, dx[0], atol=1e-6) or not np.allclose(dy, dy[0], atol=1e-6):
        raise ValueError(f"{surface.name}: grid spacing must be uniform")
    if not np.isclose(dx[0], dy[0], atol=1e-6):
        raise ValueError(f"{surface.name}: square grid cells are required")
    if expected_cell_m is not None and not np.isclose(dx[0], expected_cell_m, atol=1e-6):
        raise ValueError(f"{surface.name}: expected {expected_cell_m:g} m grid spacing, found {dx[0]:g} m")


def assert_same_grid(a: SurfaceGrid, b: SurfaceGrid) -> None:
    if a.z.shape != b.z.shape:
        raise ValueError(f"Grid-shape mismatch: {a.name} vs {b.name}")
    if not np.allclose(a.x, b.x) or not np.allclose(a.y, b.y):
        raise ValueError(f"Grid-coordinate mismatch: {a.name} vs {b.name}")
    if not np.array_equal(np.isfinite(a.z), np.isfinite(b.z)):
        raise ValueError(f"Valid-cell mask mismatch: {a.name} vs {b.name}")


def sample_surface(surface: SurfaceGrid, spacing_m: float = SAMPLE_SPACING_M) -> pd.DataFrame:
    sx = spacing_m / surface.dx
    sy = spacing_m / surface.dy
    if not np.isclose(sx, round(sx)) or not np.isclose(sy, round(sy)):
        raise ValueError("Sampling spacing must be an integer multiple of grid spacing")
    sx, sy = int(round(sx)), int(round(sy))
    X, Y = np.meshgrid(surface.x[::sx], surface.y[::sy])
    Z = surface.z[::sy, ::sx]
    mask = np.isfinite(Z)
    return pd.DataFrame({"X": X[mask], "Y": Y[mask], "Z": Z[mask]})


def load_line_union(path: str | Path, crs: str):
    if gpd is None:
        raise ImportError("geopandas and shapely are required for seismic-line distances")
    gdf = gpd.read_file(path)
    if gdf.empty or gdf.crs is None:
        raise ValueError("Seismic-line file is empty or has no CRS")
    gdf = gdf.loc[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if not set(gdf.geom_type).issubset({"LineString", "MultiLineString"}):
        raise ValueError("Seismic-line input must contain LineString or MultiLineString geometries")
    gdf = gdf.to_crs(crs)
    return gdf.geometry.union_all() if hasattr(gdf.geometry, "union_all") else gdf.unary_union


def nearest_line_distance(x: np.ndarray, y: np.ndarray, line_union) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if shapely is not None and hasattr(shapely, "points"):
        return np.asarray(shapely.distance(shapely.points(x, y), line_union), dtype=float)
    return np.array([Point(float(xx), float(yy)).distance(line_union) for xx, yy in zip(x, y)])


def save_xyz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df[["X", "Y", "Z"]].to_csv(path, sep=" ", index=False, header=False, float_format="%.7f")


def prepare_nodes(original_xyz: str | Path, lines_path: str | Path, output_dir: str | Path, crs: str) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    original = load_xyz_grid(original_xyz, "Original")
    nodes = sample_surface(original)
    line_union = load_line_union(lines_path, crs)
    nodes["distance_to_line_m"] = nearest_line_distance(nodes["X"].to_numpy(), nodes["Y"].to_numpy(), line_union)

    save_xyz(nodes, out / "CONTROL_0m_100x100m.xyz")
    summary = []
    n0 = len(nodes)
    for D in D_VALUES_M:
        retained = nodes.loc[nodes["distance_to_line_m"] > D, ["X", "Y", "Z"]].copy()
        save_xyz(retained, out / f"GGR_D{D:04d}m_RETAINED.xyz")
        summary.append({
            "D_m": int(D),
            "Removed_points": int(n0 - len(retained)),
            "Retained_points": int(len(retained)),
            "Removed_pct": 100.0 * (n0 - len(retained)) / n0,
        })

    nodes.to_csv(out / "CONTROL_nodes_with_distance.csv", index=False)
    pd.DataFrame(summary).to_csv(out / "GGR_node_exclusion_summary.csv", index=False)
    print(f"Regularly sampled nodes: {n0:,}")
    print("Regrid the control and all retained-node files using the same Petrel settings.")


def normalized_gaussian(z: np.ndarray, sigma_cells: float) -> np.ndarray:
    valid = np.isfinite(z)
    values = np.where(valid, z, 0.0)
    weights = valid.astype(float)
    smooth_values = gaussian_filter(values, sigma=sigma_cells, mode="nearest")
    smooth_weights = gaussian_filter(weights, sigma=sigma_cells, mode="nearest")
    out = np.full_like(z, np.nan, dtype=float)
    good = valid & (smooth_weights > 1e-12)
    out[good] = smooth_values[good] / smooth_weights[good]
    return out


def lii(z: np.ndarray, distances: np.ndarray, sigma_m: float, near_m: float) -> float:
    highpass = z - normalized_gaussian(z, sigma_m / GRID_CELL_M)
    near = np.isfinite(highpass) & (distances <= near_m)
    far = np.isfinite(highpass) & (distances >= FAR_MIN_M) & (distances < FAR_MAX_M)
    if not np.any(near) or not np.any(far):
        return np.nan
    denominator = np.median(np.abs(highpass[far]))
    return float(np.median(np.abs(highpass[near])) / denominator) if denominator > 0 else np.nan


def surface_distance_grid(surface: SurfaceGrid, line_union) -> np.ndarray:
    X, Y = np.meshgrid(surface.x, surface.y)
    return nearest_line_distance(X.ravel(), Y.ravel(), line_union).reshape(surface.z.shape)


def interpolate_at_picks(surface: SurfaceGrid, picks: pd.DataFrame) -> np.ndarray:
    fn = RegularGridInterpolator((surface.y, surface.x), surface.z, bounds_error=False, fill_value=np.nan)
    return fn(np.column_stack([picks["Y"].to_numpy(float), picks["X"].to_numpy(float)]))


def mae_rmse(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan
    return float(np.mean(np.abs(values))), float(np.sqrt(np.mean(values ** 2)))


def trend_metrics(processed: SurfaceGrid, original: SurfaceGrid) -> tuple[float, float, float]:
    sigma_cells = TREND_SIGMA_M / GRID_CELL_M
    a = normalized_gaussian(processed.z, sigma_cells)
    b = normalized_gaussian(original.z, sigma_cells)
    valid = np.isfinite(a) & np.isfinite(b)
    trend_rmse = float(np.sqrt(np.mean((a[valid] - b[valid]) ** 2)))

    ay, ax = np.gradient(a, processed.dy, processed.dx)
    by, bx = np.gradient(b, original.dy, original.dx)
    angle_a = np.degrees(np.arctan2(ay, ax))
    angle_b = np.degrees(np.arctan2(by, bx))
    delta = np.abs(angle_a - angle_b) % 360.0
    delta = np.minimum(delta, 360.0 - delta)
    direction_diff = float(np.nanmedian(delta[valid]))

    mag_a = np.hypot(ax, ay)
    mag_b = np.hypot(bx, by)
    corr = float(np.corrcoef(mag_a[valid], mag_b[valid])[0, 1]) if np.sum(valid) > 2 else np.nan
    return trend_rmse, direction_diff, corr


def find_surface(directory: Path, D: int) -> Path:
    candidates = [
        directory / f"GGR_D{D:04d}m.xyz",
        directory / f"GGR_{D}m.xyz",
    ]
    for path in candidates:
        if path.exists():
            return path
    matches = sorted(directory.glob(f"*{D}m*.xyz"))
    if len(matches) == 1:
        return matches[0]
    raise FileNotFoundError(f"Could not uniquely identify the D={D} m GGR surface in {directory}")


def analyze(original_xyz: str | Path, control_xyz: str | Path, ggr_dir: str | Path,
            lines_path: str | Path, picks_xyz: str | Path, output_dir: str | Path, crs: str) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    original = load_xyz_grid(original_xyz, "Original")
    control = load_xyz_grid(control_xyz, "Control")
    assert_same_grid(original, control)

    line_union = load_line_union(lines_path, crs)
    distances = surface_distance_grid(control, line_union)
    picks = _read_xyz(picks_xyz).rename(columns={"Z": "TWT"})

    control_lii = lii(control.z, distances, BASE_SIGMA_M, BASE_NEAR_M)
    if not np.isfinite(control_lii) or control_lii <= 0:
        raise ValueError("Control LII is invalid; check near-line and far-field coverage")

    rows = []
    surfaces: dict[int, SurfaceGrid] = {}
    for D in D_VALUES_M:
        surface = load_xyz_grid(find_surface(Path(ggr_dir), int(D)), f"GGR {D} m")
        assert_same_grid(surface, control)
        surfaces[int(D)] = surface

        diff = surface.z - control.z
        surf_mae, surf_rmse = mae_rmse(diff)
        abs_diff = np.abs(diff[np.isfinite(diff)])
        p95 = float(np.percentile(abs_diff, 95)) if abs_diff.size else np.nan
        far = np.isfinite(diff) & (distances >= FAR_MIN_M)
        far_rmse = float(np.sqrt(np.mean(diff[far] ** 2))) if np.any(far) else np.nan

        lii_value = lii(surface.z, distances, BASE_SIGMA_M, BASE_NEAR_M)
        reduction = 100.0 * (1.0 - lii_value / control_lii)

        pick_residual = interpolate_at_picks(surface, picks) - picks["TWT"].to_numpy(float)
        pick_mae, pick_rmse = mae_rmse(pick_residual)
        trend_rmse, direction_diff, gradient_corr = trend_metrics(surface, original)

        rows.append({
            "D_m": int(D),
            "LII": lii_value,
            "LII_Reduction_pct": reduction,
            "Surface_MAE_vs_Control_ms": surf_mae,
            "Surface_RMSE_vs_Control_ms": surf_rmse,
            "Surface_P95_abs_change_ms": p95,
            "Far_field_RMSE_ge1500m_ms": far_rmse,
            "Pick_MAE_ms": pick_mae,
            "Pick_RMSE_ms": pick_rmse,
            "LowPassTrend_RMSE_ms": trend_rmse,
            "Median_TrendDirectionDiff_deg": direction_diff,
            "TrendGradientMagnitudeCorr": gradient_corr,
        })

    metrics = pd.DataFrame(rows)
    benefit = metrics["LII_Reduction_pct"].to_numpy(float)
    cost = metrics["Surface_RMSE_vs_Control_ms"].to_numpy(float)
    B = benefit / np.nanmax(benefit)
    C = cost / np.nanmax(cost)
    metrics["J_D"] = np.sqrt((1.0 - B) ** 2 + C ** 2)
    metrics.to_csv(out / "GGR_full_metrics.csv", index=False)

    control_diff = control.z - original.z
    valid0 = np.isfinite(control_diff)
    qa = pd.DataFrame([{
        "Common_valid_cells": int(np.sum(valid0)),
        "Mean_difference_ms": float(np.mean(control_diff[valid0])),
        "RMSE_ms": float(np.sqrt(np.mean(control_diff[valid0] ** 2))),
        "R2": float(np.corrcoef(control.z[valid0], original.z[valid0])[0, 1] ** 2),
    }])
    qa.to_csv(out / "Control_vs_Original_QA.csv", index=False)

    sensitivity_rows = []
    selected = []
    cost_norm = cost / np.nanmax(cost)
    for sigma_m in SENSITIVITY_VALUES_M:
        for near_m in SENSITIVITY_VALUES_M:
            base = lii(control.z, distances, sigma_m, near_m)
            reductions = np.array([
                100.0 * (1.0 - lii(surfaces[int(D)].z, distances, sigma_m, near_m) / base)
                for D in D_VALUES_M
            ])
            benefit_norm = reductions / np.nanmax(reductions)
            J = np.sqrt((1.0 - benefit_norm) ** 2 + cost_norm ** 2)
            best_index = int(np.nanargmin(J))
            best_D = int(D_VALUES_M[best_index])
            selected.append(best_D)
            for D, reduction, score in zip(D_VALUES_M, reductions, J):
                sensitivity_rows.append({
                    "Gaussian_sigma_m": float(sigma_m),
                    "Near_line_threshold_m": float(near_m),
                    "D_m": int(D),
                    "LII_Reduction_pct": float(reduction),
                    "Surface_RMSE_vs_Control_ms": float(cost[np.where(D_VALUES_M == D)[0][0]]),
                    "J_D": float(score),
                    "Selected_for_configuration": int(D) == best_D,
                })

    pd.DataFrame(sensitivity_rows).to_csv(out / "GGR_sensitivity_all_36x10.csv", index=False)
    counts = pd.Series(selected).value_counts().reindex(D_VALUES_M, fill_value=0)
    counts.rename_axis("D_m").reset_index(name="Selection_count").to_csv(
        out / "GGR_sensitivity_selection_counts.csv", index=False
    )
    print(metrics.to_string(index=False))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Geometry-Guided Regridding (GGR) computational workflow")
    parser.add_argument("--version", action="version", version=f"GGR workflow {CODE_VERSION}")
    sub = parser.add_subparsers(dest="command", required=True)

    prepare = sub.add_parser("prepare", help="Create control and GGR retained-node files")
    prepare.add_argument("--original", required=True, help="Original 50 m structural surface XYZ")
    prepare.add_argument("--lines", required=True, help="2D seismic-line vector file")
    prepare.add_argument("--out", required=True, help="Output directory")
    prepare.add_argument("--crs", default="EPSG:32638", help="Projected CRS of the XYZ coordinates")

    analyse = sub.add_parser("analyze", help="Analyze externally regridded control and GGR surfaces")
    analyse.add_argument("--original", required=True)
    analyse.add_argument("--control", required=True)
    analyse.add_argument("--ggr-dir", required=True)
    analyse.add_argument("--lines", required=True)
    analyse.add_argument("--picks", required=True)
    analyse.add_argument("--out", required=True)
    analyse.add_argument("--crs", default="EPSG:32638")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "prepare":
        prepare_nodes(args.original, args.lines, args.out, args.crs)
    elif args.command == "analyze":
        analyze(args.original, args.control, args.ggr_dir, args.lines, args.picks, args.out, args.crs)


if __name__ == "__main__":
    main()
