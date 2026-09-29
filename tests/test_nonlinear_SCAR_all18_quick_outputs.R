local({
  root <- getwd()
  for (i in 1:20) {
    if (
      file.exists(file.path(
        root,
        "simulation_nonlinear_SCAR_all18/R/workflow.R"
      ))
    ) {
      break
    }
    up <- dirname(root)
    if (identical(up, root)) {
      stop("Open GP-PU-PoE.Rproj before running this test.")
    }
    root <- up
  }
  wf <- new.env(parent = .GlobalEnv)
  sys.source(
    file.path(root, "simulation_nonlinear_SCAR_all18/R/workflow.R"),
    envir = wf
  )
  e <- new.env(parent = .GlobalEnv)
  wf$scar_load(root, e)

  dir <- file.path(
    root,
    "results/simulations/nonlinear_SCAR/all18_quick/groups"
  )
  e$scar_assert(dir.exists(dir), "Run scripts/03_simulation_quick.R first.")
  groups <- lapply(e$scar_groups(), function(g) {
    e$scar_read_group(dir, g, full = FALSE, reps = 1:2)
  })
  m <- e$scar_bind(lapply(groups, `[[`, "metrics"))
  d <- e$scar_bind(lapply(groups, `[[`, "scores"))
  e$scar_assert(
    nrow(m) == 36L && nrow(d) == 14400L,
    "Incomplete quick results."
  )
  e$scar_summarize(m, 1:2)
  cat("PASS: quick results of the 18 methods on two datasets.\n")
})
