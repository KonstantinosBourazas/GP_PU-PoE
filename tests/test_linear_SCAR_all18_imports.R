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

  e$lscar_static_checks(root)
  versions <- e$lscar_preflight(attach = TRUE)
  for (group in e$lscar_groups()) {
    x <- e$lscar_import(root, group, "quick", versions)
    e$lscar_assert(
      x$N_NODES == 400L &&
        identical(x$HIDDEN_POS_IDX, 356:370) &&
        x$N_HIDDEN == 15L,
      paste("linear SCAR constants were changed", group)
    )
    if (!group %in% c("Competitors", "nnPU")) {
      e$lscar_assert(
        x$N_PILOT == 5L &&
          x$PILOT_BURNIN == 120L &&
          x$PRODUCTION_SAMPLING == 200L &&
          x$PRODUCTION_THIN == 2L,
        paste("Wrong quick schedule", group)
      )
    }
  }
  cat(
    "PASS: imports of the eight method groups, constants and quick schedules.\n"
  )
})
