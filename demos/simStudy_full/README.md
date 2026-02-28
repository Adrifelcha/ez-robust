# Full simulation studies

In these simulation studies, we generate datasets with four participant size levels (20, 40, 80, 160), five levels of trial-per-condition (20, 40, 80, 160, 320).

The full simulation design has a 4 x 5 x 4 factorial design with 4 participant levels, 5 trial-per-condition levels and 4 fixed true beta levels (0.0, 0.1, 0.2, 0.4)

All remaining simulation settings are done in accordance with the `simulation_settings.R` configuration:

- For each design cell, we generate 1,000 datasets from a Wiener diffusion process
- Each dataset is generated using a unique set of true parameter values following a within-subjects t-test design, where we start by selecting random values for the population intercept for the drift-rate parameter, and the hierarchical mean of the boundary and nondecision time parameters from three uniform distributions. All cell-restricted simulation studies share the same range of nondecision time parameters, but the ranges of the population intercept for the drift rate, and the population mean boundary separation vary by study.
    - The `wide-parameters` study
    - The `restricted-parameters`
- For each simulated data-set, we created a second version where some percentage of the data for each participant was replaced by contaminated data (see the README file under the /demos/ folder for full information on how contamination was handled)
- For both clean and contaminated datasets, we computed the following summary statistics: the total number of correct responses, the mean and the median of all RTs observed, and the variance and IQR-derived variance of the observed RTs.
- The conditions of primary interest for comparison follow a 2x2 design for Contamination level (Clean vs Contaminated) and Summary statistic fed into the model (Robust vs EZ-standard)


### Two simulation studies


# About the simulation design

We evaluated the performance of the robust implementation relative to the standard hierarchical Bayesian EZ-DDM in hypothesis testing. To this end, we simulated trial-level data from a Wiener diffusion process with a within-subject $t$ test design on the drift rate parameter $\delta_{p,k}$, specified through a hierarchical meta-regression model with the following structure:

$$
\begin{align}
    \delta^{\text{pred}}_{p,k} &= \mu_\delta + \beta X_k \\
    \delta^{\text{obs}}_{p,k} &\sim \mathcal{N}\left(\delta^{\text{pred}}_{p,k},\ \sigma^2_\delta\right) \\
    \alpha_{p} &\sim \mathcal{N}\left(\mu_\alpha,\ \sigma^2_\alpha\right) \\
    \tau_{p} &\sim \mathcal{N}\left(\mu_\tau,\ \sigma^2_\tau\right)
\end{align}
$$

where:
- $p$ is the participant index
- $k$ is the condition index
- $X_k \in \{0,1\}$ is a binary condition indicator
- $\beta$ captures the population-level regression coefficient for the condition effect
- $\mu_\delta$ is the population-level intercept capturing the grand mean drift rate

As such, the predicted drift rate per participant per condition $\delta_{p,k}$ is given by $\mu_\delta + \beta X_k$, with variance $\sigma^2_\delta$. Hierarchical parameters $\mu_\alpha$, $\mu_\tau$, and $\sigma^2_\alpha$, $\sigma^2_\tau$ capture the mean and variance of the group-level distributions of the boundary separation $\alpha_p$ and nondecision time $\tau_p$ parameters, respectively.

### Simulation Design

In our simulation, we generated 1,000 independent datasets per design cell in a $5 \times 4 \times 4$ factorial design defined by:

1. **Number of participants**: $P \in \{20, 40, 80, 160\}$
2. **Number of trials per condition**: $T \in \{20, 40, 80, 160, 320\}$
3. **True fixed effect size**: $\beta \in \{0.0, 0.1, 0.2, 0.4\}$

Each dataset was generated from its own set of true parameters, with fixed values for $\beta$ (according to the factorial design cell) and all variance parameters:
- $\sigma^2_\delta = 0.75$
- $\sigma^2_\alpha = 0.5$
- $\sigma^2_\tau = 0.1$

The rest of the true parameters were generated as follows. First, we drew values for $\mu_\delta$, $\mu_\alpha$, and $\mu_\tau$ from uniform distributions:

$$
\begin{align}
    \mu_\delta &\sim \mathcal{U}(-3, 3) \\
    \mu_\alpha &\sim \mathcal{U}(2, 4) \\
    \mu_\tau &\sim \mathcal{U}(0.2, 0.4)
\end{align}
$$

Then, we generated participant-specific boundary separation ($\alpha_p$) and nondecision time ($\tau_p$) parameters from the corresponding group-level normal distributions. Similarly, we drew participant-specific drift rates $\delta_{p,k}$ from the regression structure described above.

## Comparison Conditions

Each simulated dataset consisted of clean trial-level data from a Wiener diffusion process. For each dataset, we constructed a contaminated counterpart by replacing 5% of trials per participant with contaminant data, while leaving the remaining 95% of trials identical to the clean data. On both versions (clean and contaminated), we computed standard and robust summary statistics, which yielded a $2 \times 2$ factorial comparison design:

- **Data integrity**: Clean vs. Contaminated (identical except for the first 5% of trials per participant)
- **Summary statistics**: Standard EZ-DDM (mean and variance) vs. robust EZ-DDM (median and IQR-based variance)

## Contamination Procedure

We replaced the first 5% trials per participant using a two-stage procedure. For each contaminant trial, we generated a binary indicator to determine the type of contamination:

- **RT noise** (50% probability): Added uniform noise between 2s and 3s to the observed RT
- **Decision noise** (50% probability): Set the drift rate to 0 and redrew the observation
