# Geometry-Guided Regridding (GGR)

Geometry-Guided Regridding (GGR) is a post-interpretation computational workflow developed to reduce seismic-line imprint in structural surfaces derived from 2D seismic interpretation.

The method modifies the spatial support used for surface reconstruction rather than smoothing the completed surface. It identifies grid nodes according to their distance from the original 2D seismic survey lines, excludes nodes within a specified distance threshold, and reconstructs the surface from the retained spatial support.

The workflow is intended for structural surfaces affected by anisotropic sampling associated with 2D seismic survey geometry.

## Method overview

The main GGR workflow consists of the following steps:

1. Sample the original interpreted structural surface on a regular spatial grid.
2. Calculate the shortest planar distance from each sampled node to the nearest 2D seismic line.
3. Apply an exclusion distance `D`.
4. Remove sampled nodes satisfying `d <= D`.
5. Retain the original values of nodes satisfying `d > D`.
6. Reconstruct the structural surface using the retained nodes.
7. Quantify seismic-line imprint, surface modification, structural preservation, and sensitivity to the exclusion distance.

The exclusion distance is treated as a dataset-specific parameter and should not be interpreted as a universal value.

## Repository structure

```text
Geometry-Guided-Regridding-GGR/
│
├── README.md
├── LICENSE
├── requirements.txt
├── CITATION.cff
├── ALGORITHM_S1.md
│
├── src/
│   └── ggr_workflow.py
│
├── examples/
│   └── synthetic_example/
│
├── data/
│   └── synthetic/
│
├── docs/
│   └── user_guide.md
│
└── tests/
    └── test_workflow.py
```

## Computational components

The repository provides reproducible implementations for:

- geometric distance calculation between sampled nodes and 2D seismic lines;
- distance-based node exclusion;
- seismic-line imprint diagnostics;
- surface difference metrics;
- distance-band analysis;
- broad-scale structural preservation diagnostics;
- exclusion-distance sensitivity analysis;
- benefit-cost evaluation for selecting a preferred exclusion distance.

## Seismic-line imprint metric

The local imprint index is based on the ratio between short-wavelength surface variability close to seismic lines and the corresponding variability in a far-field reference zone.

The high-pass residual is defined as:

`H_S = S - G_sigma[S]`

where `G_sigma` represents Gaussian smoothing.

The Local Imprint Index (LII) is calculated as:

`LII(S) = median(|H_S|, d <= d_near) / median(|H_S|, d_far_min <= d < d_far_max)`

For the study described in the associated manuscript, the base diagnostic parameters were:

- Gaussian scale: 250 m
- near-line distance: 250 m
- far-field reference interval: 1500–3000 m

These values are dataset-specific analytical parameters.

## Sensitivity analysis

The workflow evaluates the stability of the preferred exclusion distance using multiple Gaussian scales and near-line thresholds.

For each parameter combination, imprint-reduction benefit and surface-modification cost are normalized and combined using a distance-to-ideal criterion.

The preferred exclusion distance is therefore selected as a trade-off between reducing line-related imprint and limiting unnecessary modification of the structural surface.

## Synthetic example

A synthetic dataset is included to demonstrate the open computational stages of the GGR workflow.

The synthetic example is provided for:

- software testing;
- method demonstration;
- reproducibility checks;
- sensitivity-analysis verification.

Synthetic results are not intended to reproduce the numerical results of the field case study.

## Field-data limitations

The field structural surfaces and seismic interpretation data used in the associated study are subject to data-access and ownership restrictions and are therefore not distributed in this repository.

The field surfaces were reconstructed using a fixed Petrel gridding configuration. The open-source implementation provided here reproduces the geometry processing, quantitative diagnostics, sensitivity analysis, and synthetic workflow, but does not claim to reproduce proprietary Petrel interpolation internally.

## Requirements

Python 3.10 or later is recommended.

Main dependencies are listed in `requirements.txt`.

Installation:

```bash
pip install -r requirements.txt
```

## Basic usage

Example:

```bash
python src/ggr_workflow.py
```

Detailed instructions, expected inputs, outputs, and example workflows are provided in `docs/user_guide.md`.

## Reproducibility

The repository is designed to allow independent evaluation of the open computational stages of Geometry-Guided Regridding.

The synthetic example can be used to verify:

- geometry processing;
- node exclusion;
- diagnostic calculations;
- sensitivity analysis;
- expected program behavior.

## License

This project is distributed under the MIT License.

See `LICENSE`.

## Citation

If this repository is used in academic work, please cite the associated publication.

Citation metadata will be provided in `CITATION.cff`.

## Contact

For scientific questions related to the method, please contact the corresponding author of the associated manuscript.
