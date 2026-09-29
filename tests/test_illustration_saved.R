## Run from the project root: Rscript --vanilla tests/test_illustration_saved.R
## Checks the saved mode of script 01 and the settings of a new run.
stopifnot(file.exists(".here"))
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
wf <- new.env(parent = as.environment("package:stats"))
sys.source(file.path(root, "R", "illustration", "workflow.R"), envir = wf)

## Saved mode: figure and tables from archive/illustration/, identical values.
out <- file.path(tempdir(), "illustration_saved_test")
unlink(out, recursive = TRUE)
wf$run_illustration(root, mode = "saved", action = "all", out_dir = out)
archive <- file.path(root, "archive", "illustration")

stopifnot(file.exists(file.path(
  out,
  "HPD_2x2_FOUR_GP_MODELS_TOY_SEED53_SORTED_BY_PUPOE_LATENT_MEAN_80.pdf"
)))

for (f in c(
  "toy_seed53_common_plot_order_by_PU_PoE_latent_mean.csv",
  "toy_seed53_ESS_core_parameters.csv",
  "toy_seed53_ESS_latent_probabilities_by_node.csv",
  "toy_seed53_ESS_latent_probabilities_summary.csv"
)) {
  a <- utils::read.csv(file.path(out, f), stringsAsFactors = FALSE)
  b <- utils::read.csv(file.path(archive, f), stringsAsFactors = FALSE)
  stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-10)))
}

stopifnot(!dir.exists(file.path(out, "model_checkpoints")))

## Saved mode fits nothing and does not write to archive/.
refused <- function(expr) {
  inherits(tryCatch(expr, error = function(e) e), "error")
}

stopifnot(
  refused(wf$run_illustration(
    root,
    mode = "saved",
    action = "fit",
    out_dir = out
  )),
  refused(wf$run_illustration(
    root,
    mode = "saved",
    action = "all",
    out_dir = file.path(archive, "new")
  )),
  refused(wf$run_illustration(
    root,
    mode = "quick",
    action = "check",
    out_dir = out,
    iterations = 300
  ))
)

## A new run keeps data, priors and warm up, and changes only seeds and length.
full <- wf$new_illustration_environment(root, "full", out_dir = out)

same <- wf$new_illustration_environment(
  root,
  "custom",
  out_dir = out,
  settings = wf$paper_run_settings()
)

other <- wf$new_illustration_environment(
  root,
  "custom",
  out_dir = out,
  settings = list(mcmc_seed = 54, iterations = 300, thin = 3)
)

seeds <- function(e) vapply(e$MODEL_SPECS, function(s) s$mcmc_seed, integer(1))

stopifnot(
  identical(seeds(same), seeds(full)),
  identical(seeds(other), c(540101L, 540201L, 540301L, 540401L)),
  identical(other$POSTERIOR_ITERATIONS, 300L),
  identical(other$THIN, 3L),
  identical(other$N_RETAINED, 100L),
  identical(other$BLOCKED_BURNIN, full$BLOCKED_BURNIN),
  identical(other$JOINT_BURNIN, full$JOINT_BURNIN),
  identical(other$TOY_SEED, full$TOY_SEED)
)

changed <- c("POSTERIOR_ITERATIONS", "THIN", "N_RETAINED", "MODEL_SPECS")

for (key in setdiff(full$CONFIGURATION_KEYS, changed)) {
  stopifnot(identical(get(key, envir = other), get(key, envir = full)))
}

stopifnot(identical(
  basename(
    wf$new_illustration_environment(
      root,
      "custom",
      settings = list(mcmc_seed = 54, iterations = 300, thin = 3)
    )$OUT_DIR
  ),
  "custom_seed54_iter300_thin3"
))

unlink(out, recursive = TRUE)

cat(
  "PASS: saved mode reproduces the archived tables without MCMC; new runs change only the MCMC seed and length.\n"
)
