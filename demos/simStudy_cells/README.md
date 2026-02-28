# Cell-restricted simulation studies

In these simulation studies, the number of participants is fixed to 160, and we only considered two levels of trials per condition (40 and 160).

There are four simulation studies in this folder, 

**2 levels of Drift:** The population intercept for the drift rate `mu_drift` could be extracted from a uniform defined within the `0 to 1` range (a.k.a. `lowDrift`) or within a `2 to 3` range (a.k.a. `highDrift`)

**2 levels of Bound:** The hierarchical mean for the boundary parameter `mu_bound` could be extracted from a uniform defined within the `2 to 2.5` range (a.k.a. `lowBound`) or within a `3.5 to 4` range (a.k.a. `highBound`)

All four possible combinations defined by the resulting 2x2 factorial design are explored through different simulation study conditions:

1. `lowDrift-lowBound` 

2. `lowDrift-highBound`

3. `highDrift-lowBound`

4. `highDrift-highBound`

The ranges of values used to define the "low" and "high" levels across these two parameters were selected based on the results found from the full simulation study conducted with a wide range of values (see `/output/figures/summary-statistics/full/`). For each parameter, these were the ranges of values that either maximized or minimized the difference between the mean and median RT (`lowDrift` and `highBound`, and `highDrift` and `lowBound`, respectively).

All remaining simulation settings are done in accordance with the `simulation_settings.R` configuration:

- The simulation design has a 2 x 4 factorial design with 2 trial levels and 4 fixed true beta levels (0.0, 0.1, 0.2, 0.4)
- For each design cell, we generate 1,000 datasets from a Wiener diffusion process
- Each dataset is generated using a unique set of true parameter values following a within-subjects t-test design, where we start by selecting random values for the population intercept for the drift-rate parameter, and the hierarchical mean of the boundary and nondecision time parameters from three uniform distributions. All cell-restricted simulation studies share the same range of nondecision time parameters, but the ranges of the population intercept for the drift rate, and the population mean boundary separation vary by study.
- For each simulated data-set, we created a second version where some percentage of the data for each participant was replaced by contaminated data (see the README file under the /demos/ folder for full information on how contamination was handled)
- For both clean and contaminated datasets, we computed the following summary statistics: the total number of correct responses, the mean and the median of all RTs observed, and the variance and IQR-derived variance of the observed RTs.
- The conditions of primary interest for comparison follow a 2x2 design for Contamination level (Clean vs Contaminated) and Summary statistic fed into the model (Robust vs EZ-standard)
