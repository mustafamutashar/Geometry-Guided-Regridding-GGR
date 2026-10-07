# User guide

This guide describes how to use the open computational components of Geometry-Guided Regridding (GGR).

## 1. Scope

The repository supports two reproducible stages:

1. preparation of geometry-guided retained-node datasets from an original gridded structural surface and a 2D seismic-line network;
2. quantitative analysis of externally reconstructed control and GGR surfaces.

For the field study associated with the manuscript, surface reconstruction itself was performed in Petrel using one fixed gridding configuration. This repository does not reproduce Petrel interpolation internally.

A deterministic synthetic demonstration is included only for software testing and method illustration. It is not part of the field study.

## 2. Software requirements

Python 3.10 or later is recommended.

Install dependencies from the repository root:

```bash
pip install -r requirements.txt
```

Main packages:

- NumPy
- pandas
- SciPy
- Matplotlib
- GeoPandas
- Shapely

## 3. Coordinate system and units

The workflow assumes projected planar coordinates in metres.

For the associated field study, the working coordinate system was:

`WGS 84 / UTM Zone 38N (EPSG:32638)`

Structural-surface values are interpreted as two-way travel time in milliseconds when reproducing the manuscript workflow.

Do not use geographic longitude/latitude coordinates directly for distance calculations.

## 4. Input formats

### 4.1 Structural surface

The original and reconstructed surfaces are supplied as ASCII XYZ files with three numeric columns:

```text
X Y Z
```

Example:

```text
702500.000 3512000.000 1532.418
702550.000 3512000.000 1532.511
702600.000 3512000.000 1532.607
```

The field-study workflow used a regular 50 × 50 m grid.

### 4.2 Seismic lines

The 2D seismic-line network must be provided as a vector file readable by GeoPandas, such as:

- GeoPackage (`.gpkg`)
- ESRI Shapefile (`.shp`)
- GeoJSON (`.geojson`)

Line geometries must be `LineString` or `MultiLineString`.

### 4.3 Interpretation picks

For pick-based validation, provide an ASCII XYZ file:

```text
X Y TWT
```

The script reads the third column as the observed structural value.

## 5. Stage A — prepare retained-node datasets

Run:

```bash
python ggr_workflow.py prepare \
  --original path/to/original_surface.xyz \
  --lines path/to/seismic_lines.gpkg \
  --out path/to/prepared_nodes \
  --crs EPSG:32638
```

The command:

1. reads the original 50 m structural surface;
2. samples exact surface nodes at 100 × 100 m spacing;
3. calculates the shortest planar Euclidean distance from each sampled node to the nearest seismic line;
4. writes the `D = 0` control node set;
5. removes nodes satisfying `d <= D` for exclusion distances from 100 to 1000 m;
6. writes the retained node sets without modifying their structural values.

Key outputs include:

```text
CONTROL_0m_100x100m.xyz
CONTROL_nodes_with_distance.csv
GGR_D0100m_RETAINED.xyz
GGR_D0200m_RETAINED.xyz
...
GGR_D1000m_RETAINED.xyz
GGR_node_exclusion_summary.csv
```

## 6. External surface reconstruction

The prepared retained-node files are intended to be reconstructed using one fixed interpolation/gridding configuration.

For the field study, this stage was performed in Petrel.

To preserve the experimental design, the following must remain identical for the control and every GGR case:

- gridding algorithm/configuration;
- grid extent;
- cell size;
- output geometry;
- mask;
- interpolation parameters.

The exclusion distance `D` should be the only systematic GGR variable.

This external reconstruction step is not replaced by SciPy in the field workflow.

## 7. Stage B — analyze reconstructed surfaces

After the control and all GGR surfaces have been reconstructed externally, export them as XYZ files on the same 50 m grid.

Place the GGR surfaces in one directory and use names containing the exclusion distance, for example:

```text
GGR_D0100m.xyz
GGR_D0200m.xyz
...
GGR_D1000m.xyz
```

Then run:

```bash
python ggr_workflow.py analyze \
  --original path/to/original_surface.xyz \
  --control path/to/control_surface.xyz \
  --ggr-dir path/to/ggr_surfaces \
  --lines path/to/seismic_lines.gpkg \
  --picks path/to/interpretation_picks.xyz \
  --out path/to/results \
  --crs EPSG:32638
```

## 8. Quantitative outputs

The analysis command reports:

- Local Imprint Index (LII);
- LII reduction relative to the control;
- five-point Laplacian near-line/far-field ratio and reduction relative to the control;
- surface MAE relative to the control;
- surface RMSE relative to the control;
- 95th percentile of absolute surface change;
- far-field RMSE;
- pick MAE and RMSE;
- broad-scale trend RMSE;
- median structural-direction difference;
- gradient-magnitude correlation;
- distance-to-ideal score `J_D`;
- sensitivity results across 36 analytical parameter combinations.

Output tables include:

```text
GGR_full_metrics.csv
GGR_distance_band_metrics.csv
Control_vs_Original_QA.csv
GGR_sensitivity_all_36x10.csv
GGR_sensitivity_selection_counts.csv
```

## 9. Base analytical parameters

The manuscript field case used:

- grid cell size: 50 m;
- sampled-node spacing: 100 m;
- exclusion distances: 100–1000 m at 100 m increments;
- Gaussian high-pass scale: 250 m;
- near-line threshold: 250 m;
- far-field interval: 1500–3000 m;
- broad-scale trend Gaussian scale: 1000 m;
- sensitivity values: 150, 200, 250, 300, 400, and 500 m.

These are dataset-specific analytical settings, not universal recommendations.

## 10. Demonstration dataset

Run the open demonstration with:

```bash
python examples/demonstration_dataset/generate_demo.py
```

The demonstration generates a deterministic synthetic surface and test outputs.

It is provided only to verify the software path and illustrate the method. It is not used for any numerical result reported in the manuscript.

## 11. Smoke test

Run:

```bash
python tests/smoke_test.py
```

A successful run ends with:

```text
GGR smoke test passed.
```

The smoke test checks basic numerical functions and executes the demonstration workflow in a temporary directory.

## 12. Reproducibility notes

Before comparing reconstructed field surfaces, verify that:

- all surfaces share the same X and Y coordinates;
- all surfaces have the same valid-cell mask;
- units are consistent;
- the seismic-line geometry is in the same projected CRS;
- the same external gridding settings were used for every case.

The code intentionally raises errors when incompatible grids are detected.

## 13. Field-data availability

The original field seismic interpretation data and Petrel-derived structural surfaces are not distributed in this repository. Any applicable ownership or access conditions should be documented in the manuscript Data Availability statement.

The demonstration dataset provides an open test case for the computational procedures independently of the field data.

## 14. Citation

Citation metadata are provided in `CITATION.cff`.

The repository citation should be updated with the final publication DOI after the associated article is published.
