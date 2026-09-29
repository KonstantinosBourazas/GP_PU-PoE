scar_load <- function(root, envir) {
  mod <- file.path(root, "simulation_nonlinear_SCAR_all18")
  for (f in c(
    "utilities.R",
    "registry.R",
    "archive.R",
    "tables.R",
    "engines.R",
    "checks.R"
  )) {
    sys.source(file.path(mod, "R", f), envir = envir, keep.source = FALSE)
  }
  sys.source(
    file.path(mod, "config/linear_priors.R"),
    envir = envir,
    keep.source = FALSE
  )
}

run_nonlinear_SCAR_all18 <- function(
  root,
  mode = c("archived_full", "quick", "new_full"),
  action = c("all", "check"),
  out = NULL,
  archive = NULL,
  allow_long_MCMC = FALSE
) {
  mode <- match.arg(mode)
  action <- match.arg(action)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  e <- new.env(parent = .GlobalEnv)
  scar_load(root, e)
  if (mode == "new_full" && !isTRUE(allow_long_MCMC)) {
    stop(
      "The refit is switched off. Script 04 rebuilds the tables from the archive.",
      call. = FALSE
    )
  }
  if (is.null(out)) {
    out <- file.path(
      root,
      "results/simulations/nonlinear_SCAR",
      paste0("all18_", mode)
    )
  }
  if (is.null(archive)) {
    archive <- file.path(root, "archive/simulations/nonlinear_SCAR_all18")
  }
  e$scar_static_checks(root)
  e$scar_safe_output(root, out)
  # The output folder cannot be inside the archive folder, wherever it is.
  norm <- function(x) {
    tolower(gsub(
      "\\\\",
      "/",
      normalizePath(x, winslash = "/", mustWork = FALSE)
    ))
  }
  a <- norm(archive)
  o <- norm(out)
  e$scar_assert(
    o != a && !startsWith(paste0(o, "/"), paste0(a, "/")),
    "Output cannot be inside the selected input archive."
  )
  lock <- e$scar_lock(out)
  on.exit(e$scar_unlock(lock), add = TRUE)
  log <- e$scar_log_start(out, paste0("console_", mode, "_", action, ".log"))
  on.exit(e$scar_log_end(log), add = TRUE, after = FALSE)
  start <- Sys.time()
  cat("\nNONLINEAR SCAR STUDY, 18 METHODS\nMode: ", mode, "\n", sep = "")
  if (mode == "archived_full") {
    cat("Tables rebuilt from the archived results, without MCMC.\n")
    bundle <- e$scar_collect_archive(root, archive)
    if (action == "all") {
      e$scar_write_reports(root, bundle, out)
    } else {
      s <- e$scar_summarize(bundle$metrics, 1:100)
      eta <- e$scar_eta_table(bundle$metrics)
      c <- e$scar_check_reference(root, s, eta)
      e$scar_csv(c, file.path(out, "verification/comparison_results.csv"))
      e$scar_assert(all(c$pass), "Final archive comparison failed.")
    }
  } else if (action == "check") {
    versions <- e$scar_preflight(attach = TRUE)
    for (g in e$scar_groups()) {
      engine <- e$scar_import(
        root,
        g,
        if (mode == "quick") "quick" else "full",
        versions
      )
      cat(
        "PASS: imported ",
        g,
        "; functions, priors and packages checked.\n",
        sep = ""
      )
    }
  } else {
    if (mode == "quick") {
      cat(
        "Quick mode: 18 methods on 2 datasets with short runs, not the results of the paper.\n"
      )
    }
    bundle <- e$scar_fit_all(
      root,
      out,
      if (mode == "quick") "quick" else "full"
    )
    e$scar_write_reports(root, bundle, out)
    cat("\nAll 18 methods fitted and reported.\n")
  }
  cat(
    "\nCompleted in ",
    as.numeric(difftime(Sys.time(), start, units = "secs")),
    " seconds.\nOutput: ",
    normalizePath(out, winslash = "/"),
    "\n",
    sep = ""
  )
  invisible(out)
}
