local({
  root <- getwd()
  for (i in 1:20) {
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
  wf <- new.env(parent = parent.env(.GlobalEnv))
  sys.source(
    file.path(root, "simulation_linear_SCAR_all18/R/workflow.R"),
    envir = wf
  )
  e <- new.env(parent = parent.env(.GlobalEnv))
  wf$lscar_load(root, e)

  wf$run_linear_SCAR_all18(root, mode = "archived_full", action = "check")
  cat("PASS: archived results of the linear SCAR study.\n")
})
