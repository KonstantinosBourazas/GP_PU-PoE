## Reference runs of the original functions with the short schedule of quick
## mode, used by test_sensitivity_reproducibility.R.
stopifnot(file.exists(".here"))
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
  stop("Supply an output directory and comma-separated fit IDs.")
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
wf <- new.env(parent = as.environment("package:stats"))
sys.source(file.path(root, "sensitivity", "R", "workflow.R"), envir = wf)
e <- wf$new_sensitivity_environment(root, "quick", args[1L])
dir.create(e$OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(e$FIT_CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

source <- paste(
  readLines(
    file.path(
      root,
      "sensitivity",
      "reference",
      "PU_PoE_prior_sensitivity_final.R"
    ),
    warn = FALSE
  ),
  collapse = "\n"
)

guard <- paste(
  readLines(
    file.path(root, "tests", "reference_schedule_guard.txt"),
    warn = FALSE
  ),
  collapse = "\n"
)

replacement <- paste(
  readLines(
    file.path(root, "tests", "reference_schedule_replacement.txt"),
    warn = FALSE
  ),
  collapse = "\n"
)

if (!grepl(guard, source, fixed = TRUE)) {
  stop("Schedule check not found in the original script.")
}

source <- sub(guard, replacement, source, fixed = TRUE)

for (expr in parse(text = source, keep.source = FALSE)) {
  if (
    is.call(expr) &&
      length(expr) == 3L &&
      identical(expr[[1L]], as.name("<-")) &&
      is.call(expr[[3L]]) &&
      identical(expr[[3L]][[1L]], as.name("function"))
  ) {
    eval(expr, e)
  }
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
sys.source(file.path(root, "sensitivity", "R", "prepare.R"), envir = e)
ids <- strsplit(args[2L], ",", fixed = TRUE)[[1L]]

for (id in ids) {
  spec <- e$SENSITIVITY_SPECS[[id]]
  p <- e$sensitivity_prior_bundles[[id]]
  ms <- e$make_sensitivity_model_spec(spec)
  set.seed(as.integer(spec$mcmc_seed))
  fit <- e$run_exact_toy_model_chain(
    model_spec = ms,
    prior_bundle = p,
    Y = e$Y,
    X_cov = e$X_cov,
    A_network = e$A_network,
    Z_network = e$Z_network,
    network_geometry = e$network_geometry,
    m_beta = p$beta0_prior$mean,
    s_beta = p$beta0_prior$sd,
    a_eta = p$eta_prior$a_eta,
    b_eta = p$eta_prior$b_eta,
    n_blocked_burnin = e$BLOCKED_BURNIN,
    n_joint_burnin = e$JOINT_BURNIN,
    n_sample = e$POSTERIOR_ITERATIONS,
    thin = e$THIN,
    max_log_sigma_cov = log(spec$max_sigma_cov),
    max_mean_sd_network = spec$max_mean_sd_network,
    verbose = TRUE
  )
  sums <- e$summarize_completed_model(fit, ms, p)
  result <- list(
    complete = TRUE,
    fit_spec = spec,
    model_spec = ms,
    prior_bundle = p,
    fit = fit,
    node_summary = e$add_fit_metadata(sums$node_summary, spec),
    eta_summary = e$add_fit_metadata(sums$eta_summary, spec)
  )
  saveRDS(result, file.path(e$FIT_CHECKPOINT_DIR, paste0(id, "_complete.rds")))
}

cat("Reference runs completed (short schedule).\n")
