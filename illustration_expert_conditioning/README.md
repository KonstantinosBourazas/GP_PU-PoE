# Expert scores of the illustrative example

This module decomposes the GP-PU-PoE fit of the illustrative example into its
two experts. It reports the posterior expert scores `sigma(x_0)` for the
covariates and `sigma(x_1)` for the network. These are the experts inside the
joint model, not separately fitted GP-Cov or GP-Net models. The module reads a
saved fit and runs no MCMC.

## Run

```sh
Rscript --vanilla scripts/01b_expert_conditioning.R --mode=full --action=all
```

or open `scripts/01b_expert_conditioning.R` in RStudio and click Source.

`all` computes the expert scores and writes the figure and summaries, `plot`
redraws the figure from saved summaries, and `check` validates the inputs
only. None of them refits the model. In full mode the input is the archived
full fit in `archive/illustration/`, so no MCMC run is needed, and the output
goes to `results/illustration/expert_conditioning/`. `INPUT_DIR` (or
`--input=PATH`) points to another folder with the same files, for example
`results/illustration/full/` after a full run of
`scripts/01_illustrative_example.R`, or the folder of a run with another MCMC
seed or length. The output then goes to a subfolder of
`results/illustration/expert_conditioning/` named after that folder. The mode
`quick` reads `results/illustration/quick/`, the output of a quick run of
script 01, and writes its figures, marked as quick mode, in the subfolder
`quick/`.

## Calculation

The saved fit contains draws of `plogis(f)`, `beta0`, `tau` and the kernel
parameters, but not the chains of the structured mean coefficients. These
coefficients are integrated out analytically. With

    x_j = Z_j gamma_j + g_j,
    G_j = C_j + tau_j^2 Z_j Z_j',
    H   = G_0 + G_1 + jitter I,
    r   = f - beta0 1,

the conditional distribution of expert j given a draw of `(f, beta0, tau, theta)`
has

    mean_j       = G_j H^{-1} r,
    covariance_j = G_j - G_j H^{-1} G_j.

Each expert includes its structured mean and its GP residual. The sigmoid is
applied to each conditional draw before the posterior means and the 80% HPD
intervals are computed with `coda::HPDinterval`. The Gaussian draws use seed
456 and the same RNG kind as the illustrative example, and the RNG kind and
state of the R session are restored afterwards. The draws are built from an
eigen decomposition of each conditional covariance, so their values depend on
the linear algebra library. They are identical on the same platform and agree
within Monte Carlo error across platforms.

## Interpretation

* A score above 0.5 means positive latent evidence from that expert. The two
  scores do not add up to the probability of the joint model.
* Each expert is drawn from its marginal conditional distribution, so the
  draws of the two experts are not a joint sample and should not be combined.
* The conditional means satisfy `m_0 + m_1 + jitter H^{-1} r = r`.
* `f` is recovered from `plogis(f)` with `qlogis`. Stored values of exactly 0 or
  1 are replaced by finite limits, and their counts are reported in the log and
  in `expert_conditioning_results.rds`.
* An HPD interval need not contain the posterior mean.

## Figure and files

The covariate panel uses `dodgerblue1` and the network panel `firebrick2`. Both
panels keep the ordering of the illustrative example within true zeros, hidden
positives and observed positives, based on the GP-PU-PoE posterior mean. The
figure is 14 by 5 inches, and the PNG has 300 dpi. A run writes

    HPD_1x2_EXPERTS_GP_PU_POE_TOY_SEED53_80.pdf
    HPD_1x2_EXPERTS_GP_PU_POE_TOY_SEED53_80.png
    toy_seed53_expert_probability_summaries_HPD80.csv
    expert_conditioning_results.rds
    console_all.log

`R/compute.R` contains the conditioning, `R/plot.R` the figure and
`R/workflow.R` the input checks and the output handling. The test
`tests/test_expert_conditioning.R` runs the calculation on eight draws of the
archived full fit. It checks that the results are reproducible and the inputs
unchanged, the order and labels of the figure, the RNG state and the identity
of the conditional means.
