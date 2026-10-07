#!/usr/bin/env python3
"""Fast executable check for the open GGR workflow."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
import importlib.util

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import ggr_workflow as ggr  # noqa: E402


def load_demo_module():
    demo_path = ROOT / "examples" / "demonstration_dataset" / "generate_demo.py"
    spec = importlib.util.spec_from_file_location("ggr_demo", demo_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load demonstration script")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    x = np.arange(0.0, 6000.0 + 50.0, 50.0)
    y = np.arange(0.0, 6000.0 + 50.0, 50.0)
    X, Y = np.meshgrid(x, y)
    distance = np.abs(X - 3000.0)

    z = 1500.0 + 0.01 * X + 0.005 * Y
    z += 5.0 * np.exp(-(distance / 150.0) ** 2) * np.sin(2.0 * np.pi * Y / 900.0)

    surface = ggr.SurfaceGrid(x=x, y=y, z=z, name="Smoke-test surface")
    ggr.validate_grid(surface)

    sampled = ggr.sample_surface(surface)
    assert len(sampled) > 0

    value = ggr.lii(z, distance, ggr.BASE_SIGMA_M, ggr.BASE_NEAR_M)
    assert np.isfinite(value) and value > 0.0

    mae, rmse = ggr.mae_rmse(np.array([-2.0, 0.0, 2.0]))
    assert np.isclose(mae, 4.0 / 3.0)
    assert np.isclose(rmse, np.sqrt(8.0 / 3.0))

    demo = load_demo_module()
    with tempfile.TemporaryDirectory() as temp_dir:
        metrics = demo.run_demo(Path(temp_dir))
        assert len(metrics) == len(ggr.D_VALUES_M)
        assert metrics["D_m"].tolist() == ggr.D_VALUES_M.tolist()
        assert np.all(np.isfinite(metrics["LII"]))
        assert np.all(np.isfinite(metrics["RMSE_vs_synthetic_control_ms"]))

    print("GGR smoke test passed.")


if __name__ == "__main__":
    main()
