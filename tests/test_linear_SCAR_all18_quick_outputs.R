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

  dir <- file.path(root, "results/simulations/linear_SCAR/all18_quick/groups")
  e$lscar_assert(dir.exists(dir), "Run scripts/07_linear_SCAR_quick.R first.")
  specs <- e$lscar_group_spec()
  pieces <- list()
  i <- 0L
  for (r in 1:2) {
    geom <- e$lscar_read(file.path(
      dir,
      "Competitors/generated_geometries",
      sprintf("replicate_%03d_geometry.rds", r)
    ))
    for (group in e$lscar_groups()) {
      sp <- specs[specs$group == group, , drop = FALSE]
      z <- e$lscar_read(file.path(
        dir,
        group,
        sp$checkpoint_dir,
        sprintf(sp$pattern, r)
      ))
      e$lscar_assert(isTRUE(z$complete), paste(group, r, "incomplete"))
      for (v in e$lscar_unpack_methods(z)) {
        i <- i + 1L
        pieces[[i]] <- e$lscar_validate_method(v, r, group, expected_y = geom$Y)
      }
    }
  }
  m <- e$lscar_bind(lapply(pieces, `[[`, "metrics"))
  d <- e$lscar_bind(lapply(pieces, `[[`, "scores"))
  e$lscar_assert(
    nrow(m) == 36L && nrow(d) == 14400L,
    "Incomplete quick results."
  )
  e$lscar_summarize(m, 1:2)
  cat("PASS: quick results of the 18 methods on two datasets.\n")
})
