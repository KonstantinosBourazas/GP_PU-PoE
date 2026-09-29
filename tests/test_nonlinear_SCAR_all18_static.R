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

  e$scar_static_checks(root)
  # The output lock is released after a normal run and after an error.
  out <- tempfile("scar_all18_lock_")
  dir.create(out)
  lock <- e$scar_lock(out)
  blocked <- inherits(try(e$scar_lock(out), silent = TRUE), "try-error")
  e$scar_assert(blocked, "A second run was not blocked.")
  e$scar_unlock(lock)
  e$scar_assert(!dir.exists(lock), "Lock was not released.")
  f <- function() {
    l <- e$scar_lock(out)
    on.exit(e$scar_unlock(l))
    stop("test error")
  }
  invisible(try(f(), silent = TRUE))
  e$scar_assert(
    !dir.exists(file.path(out, ".run_lock")),
    "Error cleanup did not release lock."
  )
  unlink(out, recursive = TRUE, force = TRUE)
  cat("PASS: static checks and output lock.\n")
})
