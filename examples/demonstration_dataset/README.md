# Demonstration dataset

This directory contains a deterministic synthetic demonstration used only to test the open computational stages of Geometry-Guided Regridding (GGR).

## Important distinction from the field study

The demonstration dataset is not part of the field case study and is not used to generate, validate, replace, or support any numerical result reported for the field data in the manuscript.

All numerical results reported in the manuscript are derived from the field case study. The purpose of this synthetic example is limited to software testing, method demonstration, and reproducibility of the open computational procedures when the original field data cannot be distributed.

The field surfaces in the study were reconstructed in Petrel using one fixed gridding configuration. The demonstration below uses SciPy interpolation only to provide an open test case. It must not be interpreted as a reproduction of the Petrel interpolation algorithm.

## What the example does

The script `generate_demo.py`:

1. creates a smooth synthetic structural surface on a 50 m grid;
2. adds a controlled short-wavelength imprint concentrated near a simple synthetic 2D line network;
3. samples the surface at 100 m spacing;
4. excludes sampled nodes according to GGR distance thresholds;
5. reconstructs demonstration surfaces using SciPy linear interpolation;
6. calculates LII reduction and RMSE relative to the synthetic control reconstruction;
7. writes the demonstration metrics and generated test files to an `output/` directory.

The example is deterministic and does not use random numbers.

## Run

From the repository root:

```bash
python examples/demonstration_dataset/generate_demo.py
```

Expected outputs include:

- `synthetic_original_surface.xyz`
- `synthetic_sampled_nodes.csv`
- `demo_metrics.csv`

The exact numerical values of this demonstration are not expected to match the field-study results.
