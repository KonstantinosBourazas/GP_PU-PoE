# GP-PU-PoE

This repository contains the R code for the illustrative example, the prior
sensitivity analysis and the three simulation studies of the paper

> Bourazas, K., Alexopoulos, A., Dellaportas, P. and Kalogeropoulos, K.
> Tax Fraud Detection via Bayesian Positive–Unlabeled Learning with Gaussian
> Processes and Multilayer Networks.

The model is a Bayesian positive unlabeled (PU) classifier that combines
covariates and network information through a product of Gaussian process (GP)
experts (PoE). The code and data of the application are not included.

## Requirements

R 4.3 or later. The core packages are Matrix, igraph, RSpectra, coda,
ggplot2 (3.4.0 or later) and patchwork. They are enough for the illustrative
example, the prior sensitivity analysis and the reconstruction of all
simulation tables from the archived results. Refitting the simulation methods
also needs e1071, AdaSampling and BayesLogit from CRAN and INLA from the INLA
repository.

```r
source("scripts/00_install_dependencies.R")
```

The script installs the core packages. Set `INSTALL_REFIT <- TRUE` at its top
to install the refitting packages as well. It installs only missing packages
and never upgrades installed ones.

INLA is used only by the INLA and INLA+Net competitors. A specific INLA
version can be installed with

```r
remotes::install_version("INLA", version = "24.12.11",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/testing"),
  dependencies = TRUE)
```

## Package versions

`renv.lock` records the R version and the package versions of the environment
that produced the results, together with the CRAN and INLA repositories. To
install exactly these versions in a private library for this project:

```r
install.packages("renv")
renv::activate()   # then restart R in the project folder
renv::restore()
```

R uses the private library when the project is opened in RStudio or when R is
started in the project folder. `Rscript --vanilla` skips it, so with renv run
the commands of this README without `--vanilla`. If the recorded INLA version
is no longer in the INLA repository, install it with the
`remotes::install_version()` call above.

The archived results were produced with three R versions on Windows 11, and
each archived run folder contains the `sessionInfo` of the session that
produced it. `renv.lock` records the R 4.4.0 environment.

| Results | R | Main package versions |
|---|---|---|
| Illustrative example, competitors and nnPU | 4.4.0 | igraph 2.1.4, RSpectra 0.16-1, INLA 24.12.11, AdaSampling 1.3, e1071 1.7-14 |
| GP models of the linear SCAR and nonlinear SAR studies | 4.5.2 | igraph 2.3.3, RSpectra 0.16-2 |
| GP models of the nonlinear SCAR study and all Bayesian linear models | 4.5.3 | igraph 2.3.3, RSpectra 0.16-2, BayesLogit 2.4 |

## Getting started

Open `GP-PU-PoE.Rproj` in RStudio, or set the working directory to the
project folder. Every script finds the project folder by itself and writes its
output under `results/`. The following commands take about a
minute and check that the installation works.

```sh
Rscript --vanilla tests/test_static.R
Rscript --vanilla tests/test_inputs.R
Rscript --vanilla scripts/01_illustrative_example.R
Rscript --vanilla scripts/01_illustrative_example.R --mode=quick
```

The third command draws the figure and the tables of the illustrative example
from the saved fits in a few seconds. The last one runs short MCMC chains and
checks that the sampler works.

In RStudio, open a script, set the options at its top and click Source.

## Scripts

| Script | Output | MCMC |
|---|---|---|
| `01_illustrative_example.R` | Four GP models on one synthetic dataset | none by default (saved fits), quick or full |
| `01b_expert_conditioning.R` | Covariate and network expert scores of the GP-PU-PoE fit | none, reads the archived full fit |
| `02_prior_sensitivity.R` | 13 prior settings on the dataset of the illustrative example | quick or full |
| `02b_prior_sensitivity_MC100.R` | 13 prior settings over 100 datasets, from the archive | none |
| `02c_prior_sensitivity_MC100_refit_OPTIONAL.R` | Refits the 1,300 fits read by script 02b | several days, disabled by default |
| `03_simulation_quick.R` | Nonlinear SCAR, 18 methods, two datasets | short runs |
| `04_simulation_full.R` | Nonlinear SCAR tables, 100 replications, from the archive | none |
| `04b_nonlinear_SCAR_refit_OPTIONAL.R` | Refits the nonlinear SCAR study | long, disabled by default |
| `05_nonlinear_SAR_quick.R` | Nonlinear SAR, 18 methods, two datasets | short runs |
| `06_nonlinear_SAR_full.R` | Nonlinear SAR tables, 100 replications, from the archive | none |
| `06b_nonlinear_SAR_refit_OPTIONAL.R` | Refits the nonlinear SAR study | long, disabled by default |
| `07_linear_SCAR_quick.R` | Linear SCAR, 18 methods, two datasets | short runs |
| `08_linear_SCAR_full.R` | Linear SCAR tables, 100 replications, from the archive | none |
| `08b_linear_SCAR_refit_OPTIONAL.R` | Refits the linear SCAR study | long, disabled by default |

The scripts 03, 05 and 07 and the refit scripts 04b, 06b and 08b need the
refitting packages. All other scripts need only the core packages.

## Saved, quick and full modes

Script 01 runs in three modes. The default, `saved`, draws the figure and the
tables of the paper from the saved full fits in `archive/illustration/` in a
few seconds and runs no MCMC. The modes `quick` and `full` fit the four models
again. Script 02 runs in the modes `quick` and `full`. The two MCMC modes use
the same code, data, priors and seeds, and differ only in the length of the
MCMC run.

| MCMC schedule per model | quick | full |
|---|---:|---:|
| Blocked warm up | 200 | 10,000 |
| Joint warm up | 300 | 20,000 |
| Sampling iterations | 600 | 150,000 |
| Thinning | 3 | 30 |
| Retained draws | 200 | 5,000 |

Full mode is the schedule used in the paper. The four full fits of script 01
take about 45 minutes on a laptop, and one full GP-PU-PoE fit about 13
minutes. Quick mode only checks that the code runs. Its output goes to a
separate folder, and its figures are marked as quick mode.

In full mode, script 01 can also run the MCMC with another seed or another
sampling length. Set `MCMC_SEED`, `ITERATIONS` and `THIN` at the top of the
script, or use `--mcmc-seed`, `--iterations` and `--thin`. The data, the priors
and the warm up stay those of the paper, `MCMC_SEED = 53` gives the seeds of
the paper, and the output goes to its own folder, for example
`results/illustration/custom_seed54_iter150000_thin30/`. Script 01b reads such
a run with `--input=PATH`.

## Archived results

`archive/` holds the completed fits behind the tables and figures of the
paper: the full schedule fit of the illustrative example
(`archive/illustration/`), 100 replications of the 18 methods for each
simulation study, and 1,300 fits (13 prior settings times 100 datasets) for
the prior sensitivity analysis over 100 datasets. Scripts 01 and 01b read the
fits of the illustrative example, and scripts 02b, 04, 06 and 08 rebuild the
tables and figures from the other files, all without MCMC. These four scripts
check every file against a manifest of MD5 checksums before reading it. The
archive is never modified, and new fits are always written to `results/`.

## Simulation studies

Each study compares 18 methods on 100 synthetic datasets of 400 units. The
methods are Bayesian linear models with and without the PU likelihood and the
network (Lin-Cov, Lin-PU-Cov, Lin-Cov+Net, Lin-PU-Cov+Net), GP models (GP-Cov,
GP-PU-Cov, GP-PoE, GP-PU-PoE, GP-PU-PoE_beta0) and the competitors Ada,
Ada+Net, nnPU, nnPU+Net, WLR, WLR+Net, INLA, INLA+Net and GCN, the graph
convolutional network of Kipf and Welling (2017). The original competitor code
and the archived files call the GCN method GNN. The three
studies differ in the data generating process (linear or nonlinear) and in
how the hidden positives are selected (SCAR or SAR).

Each study folder contains `reference/sources/`, the original simulation
scripts. The study code imports the functions and settings it needs from
these scripts after checking their MD5 checksums, and never runs them as a
whole.

## Reproducibility

* Every model has a fixed seed, and the random number generator is set
  explicitly with `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`.
* The synthetic data and the prior calibration of the illustrative example
  are identical on Windows (R 4.4.0, igraph 2.1.4) and Linux (R 4.3.3,
  igraph 1.6.0). The network modularity features differ by less than 1e-13.
* A full run of script 01 with the code of this repository reproduced the
  posterior draws of the original illustration script exactly (R 4.4.0,
  Windows). The expert scores of script 01b are Monte Carlo summaries of
  Gaussian draws built from an eigen decomposition. They are identical on the
  same platform and differ slightly across platforms (posterior means by less
  than 0.002).
* Completed fits are saved as checkpoints and reused only when the code, the
  settings and the R environment match. An interrupted fit restarts from the
  beginning.
* A `.run_lock` folder prevents two runs from writing to the same output
  folder. After a crash, delete only that folder, once no run is active.
* The Bayesian linear models draw the Pólya-Gamma variables with BayesLogit, as
  in the archived fits (BayesLogit 2.4). The refit scripts stop if BayesLogit is
  not installed, since pgdraw gives different draws for the same seed.
* The GP samplers keep an option for a spike and slab prior on the expert
  scales. It is switched off in every analysis (`PI_SPIKE = 0`), so the prior
  of each scale is its lognormal slab alone, and the code returns the slab
  density without extra computation or random draws.
* `.gitattributes` keeps every file byte for byte. The reference sources are
  checked by MD5 checksum, so any edit to them stops the workflows.

## Tests

Run a test from the project folder with `Rscript --vanilla tests/<name>.R`.

| Tests | What they check | Requirements |
|---|---|---|
| `test_static`, `test_sensitivity_static`, `test_*_all18_static`, `test_*_all18_function_lookup`, `test_bayesian_linear_priors`, `test_linear_SCAR_root_locator`, `test_locking` | Parsing, reference checksums, unchanged original functions, settings, priors, output locks | core packages, no MCMC |
| `test_inputs` | Synthetic data, network features, prior calibration, likelihood identities | core packages, no MCMC |
| `test_reproducibility` | Two fresh quick runs and a shortened run of the original functions give identical draws | core packages, about one minute |
| `test_sensitivity_reproducibility` | The same for the prior sensitivity code | core packages, about four minutes |
| `test_*_all18_archive`, `test_linear_SCAR_all18_failure_fields` | Archived results against stored reference values | `archive/`, no MCMC |
| `test_*_all18_imports`, `test_*_all18_quick_outputs`, `test_*_all18_reproducibility` | Import of the 18 methods and short refits | refitting packages |
| `test_expert_conditioning` | Expert scores on eight draws of the archived full fit | `archive/illustration/`, no MCMC |
| `test_illustration_saved` | Saved mode of script 01 against the archived tables, and the settings of a new run | `archive/illustration/`, no MCMC |

The workflow `.github/workflows/check.yaml` runs on GitHub Actions the tests
of the first three rows, the illustrative example and the expert scores in
quick mode, `test_illustration_saved`, `test_expert_conditioning`, and script
02b with its test. All of these pass with R 4.3.3 on Linux.

## Repository layout

```
scripts/                           entry points (see the table above)
R/illustration/, config/           illustrative example
reference/                         original script of the illustrative example
illustration_expert_conditioning/  expert scores of the illustrative example
sensitivity/                       prior sensitivity on one dataset
sensitivity_mc100/                 prior sensitivity over 100 datasets
simulation_linear_SCAR_all18/      linear SCAR study
simulation_nonlinear_SCAR_all18/   nonlinear SCAR study
simulation_nonlinear_SAR_all18/    nonlinear SAR study
archive/                           completed fits behind the tables and figures, read only
docs/                              checksum and source map of the original illustration script
tests/                             tests
results/                           all new output, not tracked by Git
```

## Funding

Alexopoulos and Bourazas gratefully acknowledge funding from the Hellenic Foundation for Research and Innovation (H.F.R.I.) under the National Recovery and Resilience Plan “Greece 2.0”, funded by the European Union—NextGenerationEU (H.F.R.I. Project No. 15973).



## Citation and license

Please cite the paper above when you use this code. `CITATION.cff` gives the
citation details. The code is released under the MIT license (`LICENSE`).

## References

Kipf, T. N. and Welling, M. (2017). Semi-Supervised Classification with Graph
Convolutional Networks. In International Conference on Learning
Representations (ICLR).
