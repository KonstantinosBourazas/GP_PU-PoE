## Checks the priors of the Bayesian linear models in the three studies.
local({
  root <- getwd()
  for (i in seq_len(50L)) {
    if (
      file.exists(file.path(root, "simulation_linear_SCAR_all18/R/workflow.R"))
    ) {
      break
    }
    up <- dirname(root)
    if (identical(up, root)) {
      stop("Open GP-PU-PoE.Rproj before running this test.")
    }
    root <- up
  }
  e <- new.env(parent = parent.env(.GlobalEnv))
  for (n in c(
    "utilities",
    "registry",
    "archive",
    "tables",
    "engines",
    "checks"
  )) {
    sys.source(
      file.path(root, "simulation_nonlinear_SCAR_all18/R", paste0(n, ".R")),
      envir = e
    )
  }
  sys.source(
    file.path(root, "simulation_nonlinear_SCAR_all18/config/linear_priors.R"),
    envir = e
  )
  e$scar_static_checks(root)
  e <- new.env(parent = parent.env(.GlobalEnv))
  for (n in c(
    "utilities",
    "registry",
    "archive",
    "tables",
    "engines",
    "checks"
  )) {
    sys.source(
      file.path(root, "simulation_nonlinear_SAR_all18/R", paste0(n, ".R")),
      envir = e
    )
  }
  sys.source(
    file.path(root, "simulation_nonlinear_SAR_all18/config/linear_priors.R"),
    envir = e
  )
  e$sar_static_checks(root)
  e <- new.env(parent = parent.env(.GlobalEnv))
  for (n in c(
    "utilities",
    "registry",
    "archive",
    "tables",
    "engines",
    "checks"
  )) {
    sys.source(
      file.path(root, "simulation_linear_SCAR_all18/R", paste0(n, ".R")),
      envir = e
    )
  }
  sys.source(
    file.path(root, "simulation_linear_SCAR_all18/config/linear_priors.R"),
    envir = e
  )
  e$lscar_static_checks(root)
  cat("PASS: priors of the Bayesian linear models in the three studies.\n")
})
