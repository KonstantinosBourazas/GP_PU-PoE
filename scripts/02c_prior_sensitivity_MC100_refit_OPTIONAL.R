## Refits the 1,300 fits of the prior sensitivity analysis over 100 datasets.
## This takes several days and is not needed for the tables and figures, which
## script 02b rebuilds from the archive. The driver of the original runs is kept
## as a function, and the sampler is read from the original script.
RUN_LONG_MCMC <- FALSE
MODE <- "run" # "run", "combine", "plot"
SESSION_ID <- 1L # 1..5 as in the original driver
REPLICATIONS <- 1:100
KEEP_FULL_DRAWS <- TRUE
VERBOSE_MCMC <- TRUE

local({
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  ref <- file.path(root, "sensitivity_mc100", "reference")
  if (!file.exists(file.path(ref, "replications_driver_functions.R"))) {
    stop("Open GP-PU-PoE.Rproj first.")
  }
  if (identical(MODE, "run") && !isTRUE(RUN_LONG_MCMC)) {
    stop(
      "The refit is switched off. Set RUN_LONG_MCMC <- TRUE to refit. ",
      "Script 02b rebuilds the figures from the archive."
    )
  }
  e <- new.env(parent = as.environment("package:stats"))
  e$ORIGINAL_SCRIPT <- file.path(ref, "prior_sensitivity_single_run_original.R")
  ## New fits are written to results/, not to archive/.
  e$OUT_DIR <- file.path(root, "results", "prior_sensitivity", "mc100_new_runs")
  for (key in c(
    "MODE",
    "SESSION_ID",
    "REPLICATIONS",
    "KEEP_FULL_DRAWS",
    "VERBOSE_MCMC"
  )) {
    assign(key, get(key, inherits = TRUE), envir = e)
  }
  sys.source(file.path(ref, "replications_driver_functions.R"), envir = e)
  e$run_prior_mc()
})
