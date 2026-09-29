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
  wf <- new.env(parent = .GlobalEnv)
  sys.source(
    file.path(root, "simulation_nonlinear_SAR_all18/R/workflow.R"),
    envir = wf
  )
  e <- new.env(parent = .GlobalEnv)
  wf$sar_load(root, e)

  e$sar_static_checks(root)
  # The output lock is released after a normal run and after an error.
  out <- tempfile("sar_all18_lock_")
  dir.create(out)
  lock <- e$sar_lock(out)
  blocked <- inherits(try(e$sar_lock(out), silent = TRUE), "try-error")
  e$sar_assert(blocked, "A second run was not blocked.")
  e$sar_unlock(lock)
  e$sar_assert(!dir.exists(lock), "Lock was not released.")
  f <- function() {
    l <- e$sar_lock(out)
    on.exit(e$sar_unlock(l))
    stop("test error")
  }
  invisible(try(f(), silent = TRUE))
  e$sar_assert(
    !dir.exists(file.path(out, ".run_lock")),
    "Error cleanup did not release lock."
  )
  unlink(out, recursive = TRUE, force = TRUE)
  cat("PASS: static checks and output lock.\n")
})
