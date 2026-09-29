# Nonlinear SCAR study

Reconstruction and refitting of the nonlinear SCAR simulation study with all 18
methods.

## Run

`scripts/04_simulation_full.R` rebuilds the tables of the paper from the 100
archived replications in `archive/simulations/nonlinear_SCAR_all18/`. It fits no
model and runs no MCMC. Output goes to
`results/simulations/nonlinear_SCAR/all18_archived_full/`.

`scripts/03_simulation_quick.R` fits all 18 methods on two new datasets of 400
units with short runs. It needs the refitting packages, and its estimates only
check that the code runs. `scripts/04b_nonlinear_SCAR_refit_OPTIONAL.R` refits the
whole study and is disabled by default.

## Bayesian linear priors

The coefficient, Pólya-Gamma, latent label, eta and tau updates are imported
unchanged from `reference/sources/Bayesian_linear.R`. The models with
covariates only use covariate scale 0.83. The models with covariates and
network use covariate scale 0.59 and a network scale calibrated from the
standardized network matrix of each dataset, rounded to two decimals. The
intercept prior is N(-3.40, 0.61^2), and eta in the PU models has a
Beta(2, 6.39) prior.

## Archive and tables

The Bayesian linear folder of the archive contains 100 production
checkpoints, 20 pilot checkpoints, five pilot geometries, the global
calibration and the run configuration. The archive reader rejects extra
Bayesian linear files and inconsistent prior or calibration signatures. The
Bayesian linear production datasets were regenerated from the recorded seeds,
and the checkpoints store the seeds, labels and summaries but not every
covariate matrix and adjacency matrix.

The competitor results are read from `replication_results` in the archived
`SIMULATION_100_COMPETING_METHODS_COMPLETE_RESULTS.RData`.

The tables are recomputed from the results of each replication and compared
with `reference/EXPECTED_ARCHIVED_SUMMARIES.csv` and
`reference/EXPECTED_ETA_SUMMARIES.csv`, which hold full precision values
computed from the archive.
