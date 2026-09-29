## Entry function of the prior sensitivity analysis over 100 datasets.

mc100_normal_path <- function(x, must_exist = FALSE) {
  normalizePath(path.expand(x), winslash = "/", mustWork = must_exist)
}

mc100_child_path <- function(path, parent) {
  path <- tolower(mc100_normal_path(path))
  parent <- tolower(mc100_normal_path(parent))
  identical(path, parent) || startsWith(path, paste0(parent, "/"))
}

mc100_release_lock <- function(lock) {
  for (i in seq_len(10L)) {
    if (!file.exists(lock)) {
      return(invisible(TRUE))
    }
    suppressWarnings(unlink(
      lock,
      recursive = TRUE,
      force = TRUE,
      expand = FALSE
    ))
    if (!file.exists(lock)) {
      return(invisible(TRUE))
    }
    Sys.sleep(0.15)
  }
  warning(
    "Could not remove the lock: ",
    lock,
    ". Remove it by hand once no run is active.",
    call. = FALSE
  )
  invisible(FALSE)
}

run_prior_sensitivity_mc100 <- function(
  root,
  action = "all",
  archive = NULL,
  out = NULL
) {
  root <- mc100_normal_path(root, TRUE)
  action <- match.arg(action, c("all", "plot", "check", "combine"))
  code_root <- file.path(root, "sensitivity_mc100")
  archive <- if (is.null(archive)) {
    file.path(root, "archive", "prior_sensitivity_mc100")
  } else {
    archive
  }
  archive <- mc100_normal_path(archive, TRUE)
  out <- if (is.null(out)) {
    file.path(root, "results", "prior_sensitivity", "mc100")
  } else {
    out
  }
  out <- mc100_normal_path(out)
  protected <- c(
    archive,
    code_root,
    file.path(root, "R"),
    file.path(root, "config"),
    file.path(root, "sensitivity"),
    file.path(root, "scripts"),
    file.path(root, "tests"),
    file.path(root, "reference"),
    file.path(root, "results", "illustration"),
    file.path(root, "results", "prior_sensitivity", "quick"),
    file.path(root, "results", "prior_sensitivity", "full")
  )
  if (
    identical(tolower(out), tolower(root)) ||
      any(vapply(
        protected,
        function(p) {
          mc100_child_path(out, p) || mc100_child_path(p, out)
        },
        logical(1)
      ))
  ) {
    stop(
      "The output folder must be outside the archive, the code folders and the other results."
    )
  }
  e <- new.env(parent = as.environment("package:stats"))
  for (name in c(
    "combine.R",
    "metrics.R",
    "original_plot_functions.R",
    "plot.R"
  )) {
    sys.source(file.path(code_root, "R", name), envir = e, keep.source = FALSE)
  }
  if (action %in% c("all", "plot")) {
    missing <- c("ggplot2", "patchwork")[
      !vapply(
        c("ggplot2", "patchwork"),
        requireNamespace,
        logical(1),
        quietly = TRUE
      )
    ]
    if (length(missing)) {
      stop(
        "Missing plotting packages: ",
        paste(missing, collapse = ", "),
        ". Run scripts/00_install_dependencies.R first."
      )
    }
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(out, ".postprocess_lock")
  if (file.exists(lock)) {
    stop(
      "The output folder is locked: ",
      lock,
      ". Remove the lock by hand once no run is active."
    )
  }
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("Cannot create the output lock: ", out)
  }
  on.exit(mc100_release_lock(lock), add = TRUE)
  start <- Sys.time()
  cat(
    "\nPRIOR SENSITIVITY OVER 100 DATASETS\nAction: ",
    action,
    "\nArchive: ",
    archive,
    "\nOutput: ",
    out,
    "\n\n",
    sep = ""
  )
  bundle_file <- file.path(out, "combined", "MC100_combined_summaries.rds")
  if (action == "check") {
    e$mc100_manifest_check(root, archive)
    for (p in list.files(
      file.path(code_root, "R"),
      "[.]R$",
      full.names = TRUE
    )) {
      parse(p)
    }
    cat(
      "Checked: 1,300 checkpoints and 39 tables against the manifest; the code parses.\n"
    )
    cat(
      "ACTION 'all' checks the contents, rebuilds the tables and draws the figures.\n"
    )
    return(invisible(TRUE))
  }
  if (action %in% c("all", "combine")) {
    bundle <- e$mc100_combine(root, archive, out)
  } else {
    if (!file.exists(bundle_file)) {
      stop("No combined bundle yet. Run the script once with ACTION 'all'.")
    }
    e$mc100_manifest_check(root, archive)
    bundle <- readRDS(bundle_file)
    e$mc100_plot_input_check(bundle)
    expected <- utils::read.csv(
      file.path(code_root, "reference", "ARCHIVE_MANIFEST.csv"),
      stringsAsFactors = FALSE
    )
    if (!identical(bundle$input_manifest, expected)) {
      stop("The combined bundle refers to a different archive.")
    }
  }
  if (action %in% c("all", "combine")) {
    ess <- e$mc100_ess_summaries(archive)
    e$mc100_table(
      ess$by_replication,
      file.path(out, "combined", "MC100_ESS_by_replication.csv")
    )
    e$mc100_table(
      ess$summary,
      file.path(out, "combined", "MC100_ESS_summary.csv")
    )
    cat("\nESS OF THE FITTED PROBABILITIES OVER THE 100 DATASETS\n")
    print(
      ess$summary[c(
        "Setting",
        "Replications",
        "Mean_of_median_ESS",
        "SD_of_median_ESS",
        "Mean_of_mean_ESS",
        "SD_of_mean_ESS"
      )],
      row.names = FALSE,
      digits = 7
    )
  }
  if (action %in% c("all", "plot")) {
    e$mc100_make_figures(bundle, out)
  }
  ## The archive is checked again after the output is written.
  e$mc100_manifest_check(root, archive)
  code_files <- list.files(
    file.path(code_root, "R"),
    "[.]R$",
    full.names = TRUE
  )
  e$mc100_table(
    data.frame(
      file = basename(code_files),
      md5 = unname(tools::md5sum(code_files))
    ),
    file.path(out, "verification", "postprocessing_code_md5.csv")
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(out, "verification", "postprocessing_sessionInfo.txt")
  )
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  writeLines(
    c(
      "Prior sensitivity over 100 datasets: tables and figures rebuilt from the archive.",
      paste0("Action: ", action),
      paste0("Elapsed seconds: ", elapsed),
      "Inputs: 1,300 archived full-fit summaries, 100 replications per prior setting.",
      bundle$interval_note,
      "The HPD intervals and ESS of each fit are the archived values; the posterior draws are not archived."
    ),
    file.path(out, "verification", "RUN_REPORT.txt")
  )
  cat("\nCompleted.\n")
  cat(
    "Elapsed: ",
    round(elapsed, 2),
    " seconds.\nResults: ",
    out,
    "\n",
    sep = ""
  )
  invisible(bundle)
}
