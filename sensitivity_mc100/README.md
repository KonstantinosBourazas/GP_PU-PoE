# Prior sensitivity over 100 datasets

This module rebuilds the prior sensitivity results over 100 datasets from the
1,300 archived fits (13 prior specifications times 100 datasets). It reads saved
summaries only and runs no MCMC. The fits can be regenerated with
`scripts/02c_prior_sensitivity_MC100_refit_OPTIONAL.R`, which takes several
days.

## Run

Open `scripts/02b_prior_sensitivity_MC100.R` and click Source with
`ACTION <- "all"`, or run

```sh
Rscript --vanilla scripts/02b_prior_sensitivity_MC100.R --action=all
```

The input is `archive/prior_sensitivity_mc100/`, and the output goes to
`results/prior_sensitivity/mc100/`. After a first run with `all`,
`ACTION <- "plot"` redraws the two figures only.

## Figures

* The prediction figure keeps the 6 by 3 layout of the analysis on one dataset.
  The 100 units stay at their original indices in every dataset, without
  sorting by posterior score, and the value at each index is averaged over the
  100 datasets. A unit index is a position in the simulation design, not a unit
  shared across datasets. AUC and rank summaries are computed within each
  dataset, not from the averaged probabilities.
* The eta figure has six groups of three capped bars. Dots are averaged
  posterior means, and the caps are averaged 80% HPD limits within fits. The
  default estimate appears in every group, so the 18 bars show 13 settings.
* Eta is the false negative rate `Pr(Y=0 | T=1)`. In this PU model the false
  positive rate `Pr(Y=1 | T=0)` is zero.
* The figures have no embedded captions. Caption texts are written to
  `figures/Figure_captions.txt`.

The averaged HPD limits are averages of the limits within fits. They are not a
pooled posterior interval or a confidence interval. Coverage is computed from
the 100 intervals within fits.

## Tables in `combined/`

- `MC100_AUC_summary_and_paired_differences.csv`: overall and hidden AUC by
  setting (mean, SD, minimum, maximum, 5th and 95th percentiles), paired
  differences from DEFAULT with their Monte Carlo standard errors, and AUC
  correlations with DEFAULT.
- `MC100_AUC_pairwise_correlations.csv` and four
  `MC100_AUC_correlation_matrix_*.csv` files: Pearson and Spearman correlations
  across datasets for every pair of settings, for overall and hidden AUC. A
  constant series has an undefined correlation, reported as NA.
- `MC100_AUC_ranges_by_replication.csv` and `MC100_AUC_range_summary.csv`: the
  range across settings within each dataset, next to the range of the mean
  AUCs across settings.
- `MC100_ranking_stability_vs_DEFAULT.csv`: Spearman correlations of the unit
  scores with DEFAULT within each dataset, for all 100 units and for the 80
  unlabeled units.
- `MC100_posterior_mean_variation_by_original_node.csv`: SD, Monte Carlo
  standard error and 5th and 95th percentiles of the posterior means across
  datasets at each original index. The HPD bands of the figure do not show this
  variation.
- Coverage, bias, RMSE and Monte Carlo standard errors for eta.
- `MC100_ESS_by_replication.csv` and `MC100_ESS_summary.csv`: the effective
  sample sizes of the 100 fitted probabilities saved in each fit (median, mean,
  minimum, quartiles and maximum), and for each setting their mean, SD, median
  and minimum over the 100 datasets.

The paired differences and the ranges are reported next to the correlations,
since a high correlation of AUCs can hide a constant shift. The 5th and 95th
percentiles describe the spread over the 100 datasets.

## Checks

Before combining, the module checks the archive manifest, the checkpoint
identities, the seeds, the metadata, the common datasets, the full run
settings, the probability summaries, the AUCs and the eta summaries. The new
aggregates are compared with the reference tables in `reference/`. After a run
with `all`,

```r
source("sensitivity_mc100/tests/test_postprocessing.R")
```

checks the output data and labels. The HPD intervals and ESS values within fits
come from the saved checkpoints, since the posterior draws are not archived.
