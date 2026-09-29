# Prior sensitivity on one dataset

This module refits GP-PU-PoE on the dataset of the illustrative example under
13 prior specifications. The dataset is the synthetic dataset with seed 53:
100 units, 70 true zeros, 10 hidden positives, 20 observed positives, ten
covariates and one network. The modularity features use seed 53043, and the
MCMC seed 530401 is reset before each fit.

## Run

```sh
Rscript --vanilla scripts/02_prior_sensitivity.R --mode=quick --action=all
Rscript --vanilla scripts/02_prior_sensitivity.R --mode=full --action=all
Rscript --vanilla scripts/02_prior_sensitivity.R --mode=full --action=plot
```

or open `scripts/02_prior_sensitivity.R`, set `MODE` and click Source. The
packages are those of the illustrative example. `--out=PATH` selects another
output folder. With `--action=fit`, `--fit-ids=DEFAULT,ETA_DIFFUSE` fits a
subset, and a later run in the same folder completes the rest. The final figure
needs all 13 fits.

| MCMC schedule per fit | quick | full |
|---|---:|---:|
| Blocked warm up | 200 | 10,000 |
| Joint warm up | 300 | 20,000 |
| Sampling iterations | 600 | 150,000 |
| Thinning | 3 | 30 |
| Retained draws | 200 | 5,000 |

A full fit takes about 13 minutes on a laptop, so the 13 full fits take about
three hours. Quick mode only checks that the code runs, and its figures are
marked as quick mode.

## Prior specifications

The 13 specifications are the default, ten alternatives that each change one
prior group, and two joint stress tests.

| Group | Conservative | Default | Diffuse |
|---|---|---|---|
| Total residual GP SD target | 1 | 2 | 4 |
| Fractions of the attainable correlation range | 0.05, 0.25 | 0.10, 0.50 | 0.10, 0.90 |
| Mean block RMS multiplier | 2.08 | 4.16 | 8.32 |
| Central 95% interval of the baseline risk | 0.005, 0.05 | 0.01, 0.10 | 0.02, 0.20 |
| Eta prior | Beta(1, 5) | Beta(2, 6.39) | Beta(1, 1) |
| Joint | all conservative | default | all diffuse |

The mean block multiplier is the calibration parameter `MEAN_RMS_95_MULT`, not a
Half-t scale. When only the amplitude changes, the mean scales stay tied to the
default GP SD budget. When the lengthscale targets change, the amplitude
calibration is recomputed at the new lengthscale to keep the effective SD
target. The diffuse amplitude and the joint diffuse specifications use caps of
8 instead of 5 on the covariate amplitude and on the mean marginal network SD,
as in the original code.

The figure has six rows and three columns. The default fit appears in the
middle column of every row, so the 18 panels come from 13 fits. All panels use
the ordering of the default fit within groups and show posterior means and 80%
HPD intervals of `sigma(f_i)` from `coda::HPDinterval`.

## Code and outputs

The sampler and the other numerical functions come from `R/illustration/`. The
prior functions of this analysis are in `sensitivity/R/source_functions.R`. They
replace the calibration helper of the illustrative example only within this
analysis. The original script is kept in `sensitivity/reference/`, and
`SOURCE_FUNCTION_MAP.csv` there gives the lines of each copied function. The
script is never run as a whole. The checks before fitting and the tests read
the original specifications and functions from it.

Output goes to `results/prior_sensitivity/quick/` or
`results/prior_sensitivity/full/`. The main files are

- `fit_checkpoints/<SPECIFICATION>_complete.rds`
- `sensitivity_summaries.rds` and `completed_run_metadata.rds`
- `toy_seed53_PU_PoE_prior_sensitivity_13fits_fit_summary.csv`
- `toy_seed53_PU_PoE_prior_sensitivity_13fits_eta_HPD80.csv`
- `toy_seed53_PU_PoE_prior_sensitivity_13fits_ESS_core_parameters.csv`
- `toy_seed53_PU_PoE_prior_sensitivity_13fits_ESS_latent_probabilities_summary.csv`
- `PU_PoE_PRIOR_SENSITIVITY_6x3_JOINT_STRESS_HPD80.pdf` (full mode only)
- `default_vs_illustration.csv`

`plot` reads saved summaries and never calls the sampler. `diagnostics`
recomputes the summaries from the saved checkpoints. `fit` completes the fits
without the figure. `check` and `prepare` run no MCMC. Completed fits are reused
only when the code, the settings, the priors, the data and the R environment
match.

The default specification is always fitted here. If the GP-PU-PoE fit of the
illustrative example exists in `results/illustration/<mode>/`, the two fits are
compared, and different draws under identical inputs, settings, environment
and seed stop the run. Without that fit, the comparison is recorded as not
available.

## Test

`tests/test_sensitivity_reproducibility.R` runs the static checks, two fresh
quick runs of all 13 fits and shortened runs of the original functions for five
specifications (DEFAULT, ELL_DIFFUSE, ETA_CONSERVATIVE, JOINT_CONSERVATIVE and
JOINT_DIFFUSE), and checks plotting, reuse and locks.
