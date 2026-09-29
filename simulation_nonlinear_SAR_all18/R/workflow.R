sar_load <- function(root, envir) {
  mod <- file.path(root, "simulation_nonlinear_SAR_all18")
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

run_nonlinear_SAR_all18 <- function(
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
  e <- new.env(parent = parent.env(.GlobalEnv))
  sar_load(root, e)
  if (mode == "new_full" && !isTRUE(allow_long_MCMC)) {
    stop(
      "The refit is switched off. Script 06 rebuilds the tables from the archive.",
      call. = FALSE
    )
  }
  if (is.null(out)) {
    out <- file.path(
      root,
      "results/simulations/nonlinear_SAR",
      paste0("all18_", mode)
    )
  }
  if (is.null(archive)) {
    archive <- file.path(root, "archive/simulations/nonlinear_SAR_all18")
  }
  e$sar_static_checks(root)
  e$sar_safe_output(root, out)
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
  e$sar_assert(
    o != a && !startsWith(paste0(o, "/"), paste0(a, "/")),
    "Output cannot be inside the selected input archive."
  )
  lock <- e$sar_lock(out)
  on.exit(e$sar_unlock(lock), add = TRUE)
  log <- e$sar_log_start(out, paste0("console_", mode, "_", action, ".log"))
  on.exit(e$sar_log_end(log), add = TRUE, after = FALSE)
  start <- Sys.time()
  cat("\nNONLINEAR SAR STUDY, 18 METHODS\nMode: ", mode, "\n", sep = "")
  if (mode == "archived_full") {
    cat("Tables rebuilt from the archived results, without MCMC.\n")
    bundle <- e$sar_collect_archive(root, archive)
    if (action == "all") {
      e$sar_write_reports(root, bundle, out)
    } else {
      s <- e$sar_summarize(bundle$metrics, 1:100)
      eta <- e$sar_eta_table(bundle$metrics)
      c <- e$sar_check_reference(root, s, eta)
      e$sar_csv(c, file.path(out, "verification/comparison_results.csv"))
      e$sar_assert(
        all(c$numeric_check_pass),
        "Archived numerical comparison failed."
      )
    }
  } else if (action == "check") {
    versions <- e$sar_preflight(attach = TRUE)
    for (g in e$sar_groups()) {
      engine <- e$sar_import(
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
    bundle <- e$sar_fit_all(root, out, if (mode == "quick") "quick" else "full")
    e$sar_write_reports(root, bundle, out)
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
