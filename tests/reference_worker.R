## Reference runs of the original functions with the short schedule of quick
## mode, used by test_reproducibility.R.
stopifnot(file.exists(".here"))
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop("Supply one output directory to reference_worker.R")
}

wf <- new.env(parent = as.environment("package:stats"))
sys.source("R/illustration/workflow.R", envir = wf)
wf$check_illustration_dependencies()
e <- wf$new_illustration_environment(getwd(), "quick", args[1L])
dir.create(e$OUT_DIR, recursive = TRUE, showWarnings = FALSE)
e$MODEL_CHECKPOINT_DIR <- file.path(e$OUT_DIR, "model_checkpoints")
dir.create(e$MODEL_CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

old <- paste(
  readLines("reference/GP_four_model_illustrtation.R", warn = FALSE),
  collapse = "\n"
)

guard <- paste(
  readLines("tests/reference_schedule_guard.txt", warn = FALSE),
  collapse = "\n"
)

replacement <- paste(
  readLines("tests/reference_schedule_replacement.txt", warn = FALSE),
  collapse = "\n"
)

if (!grepl(guard, old, fixed = TRUE)) {
  stop("Schedule check not found in the original script.")
}

old <- sub(guard, replacement, old, fixed = TRUE)

for (expr in parse(text = old, keep.source = FALSE)) {
  if (
    is.call(expr) &&
      identical(expr[[1L]], as.name("<-")) &&
      is.call(expr[[3L]]) &&
      identical(expr[[3L]][[1L]], as.name("function"))
  ) {
    eval(expr, envir = e)
  }
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
## prepare.R holds the data preparation of the original script.
sys.source("R/illustration/prepare.R", envir = e)

for (spec in e$MODEL_SPECS) {
  p <- e$model_prior_bundles[[spec$model_id]]
  set.seed(as.integer(spec$mcmc_seed))
  fit <- e$run_exact_toy_model_chain(
    model_spec = spec,
    prior_bundle = p,
    Y = e$Y,
    X_cov = e$X_cov,
    A_network = e$A_network,
    Z_network = e$Z_network,
    network_geometry = e$network_geometry,
    m_beta = e$beta0_prior_for_sampler$mean,
    s_beta = e$beta0_prior_for_sampler$sd,
    a_eta = e$eta_prior_for_sampler$a_eta,
    b_eta = e$eta_prior_for_sampler$b_eta,
    n_blocked_burnin = e$BLOCKED_BURNIN,
    n_joint_burnin = e$JOINT_BURNIN,
    n_sample = e$POSTERIOR_ITERATIONS,
    thin = e$THIN,
    verbose = TRUE
  )
  summaries <- e$summarize_completed_model(fit, spec, p)
  result <- c(
    list(complete = TRUE, model_spec = spec, prior_bundle = p, fit = fit),
    summaries
  )
  saveRDS(
    result,
    file.path(e$MODEL_CHECKPOINT_DIR, paste0(spec$model_id, "_complete.rds"))
  )
}

cat("Reference runs completed (short schedule).\n")
