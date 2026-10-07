# Algorithm S1 — Geometry-Guided Regridding (GGR)

## Purpose

Geometry-Guided Regridding (GGR) is a post-interpretation workflow designed to reduce seismic-line imprint in structural surfaces derived from 2D seismic interpretation.

The method changes the spatial support used for reconstruction rather than applying smoothing directly to the completed structural surface.

## Inputs

The workflow requires:

1. an interpreted structural surface;
2. the geometry of the original 2D seismic lines;
3. a regular set of sampled surface nodes;
4. one or more exclusion distances, `D`;
5. a fixed gridding configuration for all reconstructed surfaces.

For the field application described in the associated manuscript, the original structural surface used a 50 × 50 m grid, and surface nodes were sampled at 100 × 100 m spacing.

## Step 1 — Sample the original surface

Sample the original interpreted surface on a regular spatial lattice.

No interpolation is applied during this sampling step. Each sampled node retains the value of the original structural surface at that location.

Let the sampled nodes be:

`P_i = (x_i, y_i, S_i)`

where `S_i` is the structural-surface value at node `i`.

## Step 2 — Calculate distance to the seismic-line network

For every sampled node, calculate the shortest planar Euclidean distance to the nearest 2D seismic line.

`d_i = min distance(P_i, seismic-line network)`

The resulting distance field describes the spatial relationship between the sampled structural surface and the survey-line geometry.

## Step 3 — Define the control case

The control reconstruction uses all sampled nodes.

This corresponds to:

`D = 0`

The control surface is reconstructed using the same gridding configuration later used for all GGR cases.

## Step 4 — Apply geometry-based exclusion

For a selected exclusion distance `D`, remove nodes satisfying:

`d_i <= D`

Retain nodes satisfying:

`d_i > D`

The retained node values are not modified.

Therefore, GGR changes the spatial support of the reconstruction but does not alter the retained structural values.

## Step 5 — Reconstruct the surface

Reconstruct a structural surface from the retained nodes.

For valid comparison, all GGR cases and the control case must use the same:

- gridding algorithm or software configuration;
- grid extent;
- grid-cell size;
- output geometry;
- mask;
- interpolation settings.

The exclusion distance `D` should be the only systematic GGR variable.

In the associated field study, reconstruction was performed in Petrel using one fixed gridding configuration for the control and all GGR cases.

## Step 6 — Quantify seismic-line imprint

For a reconstructed surface `S`, calculate a short-wavelength residual:

`H_S = S - G_sigma[S]`

where `G_sigma` denotes Gaussian smoothing with scale `sigma`.

The Local Imprint Index (LII) is:

`LII(S) = median(|H_S|, d <= d_near) / median(|H_S|, d_far_min <= d < d_far_max)`

The base field-study parameters were:

- `sigma = 250 m`
- `d_near = 250 m`
- far-field interval = 1500–3000 m

These are analytical parameters for the dataset and are not universal thresholds.

## Step 7 — Calculate imprint reduction

Relative imprint reduction for exclusion distance `D` is calculated against the control surface:

`R_D = 100 × (1 - LII_D / LII_Control)`

Higher `R_D` indicates stronger suppression of seismic-line imprint.

## Step 8 — Quantify surface modification

Compare every GGR surface with the control surface using metrics such as:

- mean absolute error (MAE);
- root mean square error (RMSE);
- 95th percentile of absolute change;
- far-field RMSE;
- distance-band statistics.

For example:

`MAE_D = (1/N) × sum(|S_D(i) - S_Control(i)|)`

`RMSE_D = sqrt((1/N) × sum((S_D(i) - S_Control(i))^2))`

## Step 9 — Evaluate structural preservation

Broad-scale structural preservation can be evaluated after Gaussian smoothing at a larger scale.

The associated study used a 1000 m broad-scale trend and evaluated:

- trend RMSE;
- median structural-direction difference;
- gradient-magnitude correlation.

These diagnostics are complementary and are not used as substitutes for the imprint metric.

## Step 10 — Sensitivity analysis

Repeat the diagnostic calculations across multiple combinations of:

- Gaussian scale;
- near-line threshold;
- exclusion distance.

The associated study tested Gaussian scales and near-line thresholds of:

`150, 200, 250, 300, 400, 500 m`

for a total of 36 parameter combinations.

## Step 11 — Select the preferred exclusion distance

Normalize imprint-reduction benefit and surface-modification cost:

`B_D = R_D / max(R_D)`

`C_D = RMSE_D / max(RMSE_D)`

Calculate the distance-to-ideal score:

`J_D = sqrt((1 - B_D)^2 + C_D^2)`

The preferred distance is:

`D* = argmin(J_D)`

This criterion balances imprint suppression against unnecessary surface modification.

The preferred `D` is dataset-specific and must be recalibrated for each survey.

## Reproducibility boundary

The open implementation in this repository reproduces the geometry processing, diagnostics, sensitivity analysis, and synthetic demonstration.

It does not claim to reproduce proprietary Petrel interpolation internally.

Where field seismic interpretation data cannot be distributed because of ownership or access restrictions, the synthetic example provides an open test case for evaluating the computational workflow.
