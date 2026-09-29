local({
  root <- getwd()
  for (i in 1:20) {
    if (
      file.exists(file.path(
        root,
        "simulation_nonlinear_SAR_all18/R/workflow.R"
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
  wf <- new.env(parent = parent.env(.GlobalEnv))
  sys.source(
    file.path(root, "simulation_nonlinear_SAR_all18/R/workflow.R"),
    envir = wf
  )
  e <- new.env(parent = parent.env(.GlobalEnv))
  wf$sar_load(root, e)

  dir <- file.path(root, "results/simulations/nonlinear_SAR/all18_quick/groups")
  e$sar_assert(dir.exists(dir), "Run scripts/05_nonlinear_SAR_quick.R first.")
  specs <- e$sar_group_spec()
  pieces <- list()
  i <- 0L
  for (r in 1:2) {
    geom <- e$sar_read(file.path(
      dir,
      "Competitors/generated_geometries",
      sprintf("replicate_%03d_geometry.rds", r)
    ))
    for (group in e$sar_groups()) {
      sp <- specs[specs$group == group, , drop = FALSE]
      z <- e$sar_read(file.path(
        dir,
        group,
        sp$checkpoint_dir,
        sprintf(sp$pattern, r)
      ))
      e$sar_assert(isTRUE(z$complete), paste(group, r, "incomplete"))
      e$sar_validate_selection(z, r, as.integer(geom$Y))
      for (v in e$sar_unpack_methods(z)) {
        i <- i + 1L
        pieces[[i]] <- e$sar_validate_method(v, r, group, expected_y = geom$Y)
      }
    }
  }
  m <- e$sar_bind(lapply(pieces, `[[`, "metrics"))
  d <- e$sar_bind(lapply(pieces, `[[`, "scores"))
  e$sar_assert(nrow(m) == 36L && nrow(d) == 14400L, "Incomplete quick results.")
  e$sar_summarize(m, 1:2)
  cat("PASS: quick results of the 18 methods on two datasets.\n")
})
