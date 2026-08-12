# Bayesian Spatial Heterogeneity Index (SpHI)

This repository contains R code for the simulation study accompanying the Bayesian spatial heterogeneity index (SpHI) framework. The framework is designed to quantify and rank regional heterogeneity in health outcomes across population subgroups while accounting for patient-level covariates and lower-level spatial dependence.

## Simulation design

The simulation study considers three levels of lower-level spatial dependence:

| Scenario | Spatial dependence | Purpose |
|---|---|---|
| 1 | None | Evaluate both models when no lower-level spatial effect is present |
| 2 | Low | Evaluate performance under modest lower-level spatial dependence |
| 3 | High | Evaluate performance under strong lower-level spatial dependence |

Each scenario contains 100 independently generated datasets. Two models are fitted to every dataset:

- **HI model:** the nonspatial heterogeneity-index model.
- **SpHI model:** the spatial heterogeneity-index model, which accounts for lower-level spatial dependence.

This design permits assessment of parameter recovery, estimation of the heterogeneity index, and recovery of the regional heterogeneity rankings under correct model specification and spatial misspecification.

## Repository structure

```text
sphi/
├── README.md
├── gen_nrep_datasets.R
├── heter_data_gen_fn.R
├── spheter_data_gen_fn.R
├── mcmc_fn_heter.R
├── mcmc_fn_spheter.R
├── fit_nonspheter_sce1.R
├── fit_nonspheter_sce2.R
├── fit_nonspheter_sce3.R
├── fit_spheter_sce1.R
├── fit_spheter_sce2.R
└── fit_spheter_sce3.R
```

### Data-generation scripts

- `gen_nrep_datasets.R`: generates the replicate datasets for all simulation scenarios.
- `heter_data_gen_fn.R`: contains functions for generating data without lower-level spatial dependence.
- `spheter_data_gen_fn.R`: contains functions for generating data with lower-level spatial dependence.

### Model-fitting functions

- `mcmc_fn_heter.R`: implements the MCMC sampler for the nonspatial HI model.
- `mcmc_fn_spheter.R`: implements the MCMC sampler for the spatial SpHI model.

### Model-fitting scripts

| Script | Data scenario | Fitted model |
|---|---|---|
| `fit_nonspheter_sce1.R` | No spatial dependence | HI |
| `fit_nonspheter_sce2.R` | Low spatial dependence | HI |
| `fit_nonspheter_sce3.R` | High spatial dependence | HI |
| `fit_spheter_sce1.R` | No spatial dependence | SpHI |
| `fit_spheter_sce2.R` | Low spatial dependence | SpHI |
| `fit_spheter_sce3.R` | High spatial dependence | SpHI |

## Running the simulation study

1. Clone or download this repository.
2. Open the project directory in R or RStudio.
3. Install the R packages imported by the scripts if they are not already installed.
4. Replace any computer-specific absolute paths in the scripts with the location of the project on your computer, or use project-relative paths.
5. Run `gen_nrep_datasets.R` to generate the simulated datasets.
6. Run the six model-fitting scripts to fit the HI and SpHI models under the three scenarios.

The fitting scripts save one RDS file per simulation replicate. They are written to skip completed fits, which allows an interrupted simulation study to be restarted without overwriting existing results.

## Expected output directories

```text
simdata_scenario1/
simdata_scenario2/
simdata_scenario3/
mcmc_nonsp_scenario1/
mcmc_nonsp_scenario2/
mcmc_nonsp_scenario3/
mcmc_sp_scenario1/
mcmc_sp_scenario2/
mcmc_sp_scenario3/
```

The generated datasets and MCMC output files are not stored in this repository because they may be large and can be reproduced from the supplied scripts.

## Reproducibility

The data-generation and model-fitting scripts use replicate-specific random seeds. Before running the full study, users should verify that:

- the required packages load successfully;
- the spatial correlation function uses the same range-parameter convention in data generation and model fitting;
- the MCMC chains have adequate convergence and effective sample sizes; and
- all 100 output files are present for each model–scenario combination.

## Data availability

This repository contains simulation code only. It does not include patient-level Healthcare Cost and Utilization Project (HCUP) data or other restricted health data.


