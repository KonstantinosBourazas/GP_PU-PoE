## Entry functions of the prior sensitivity analysis on the dataset of the
## illustrative example (13 fits).

sensitivity_core <- function(root) {
  e <- new.env(parent = as.environment("package:stats"))
  sys.source(file.path(root, "R", "illustration", "workflow.R"), envir = e)
  e
}

check_sensitivity_core <- function(root) {
  manifest <- utils::read.csv(
    file.path(root, "sensitivity", "reference", "ILLUSTRATION_CORE_MD5.csv"),
    stringsAsFactors = FALSE
  )
  paths <- file.path(root, manifest$file)
  bad <- !file.exists(paths)
  if (any(bad)) {
    stop(
      "Files of the illustrative example are missing: ",
      paste(manifest$file[bad], collapse = ", ")
    )
  }
  actual <- unname(tools::md5sum(paths))
  bad <- is.na(actual) | actual != manifest$md5
  if (any(bad)) {
    stop(
      "Files of the illustrative example have changed: ",
      paste(manifest$file[bad], collapse = ", "),
      ". See sensitivity/reference/ILLUSTRATION_CORE_MD5.csv."
    )
  }
  invisible(paths)
}

sensitivity_files <- function(root) {
  ## The fits of the illustrative example record the checksums of all files in
  ## R/ and config/, so the sensitivity code is kept in sensitivity/.
  old <- check_sensitivity_core(root)
  new <- c(
    file.path(root, "sensitivity", "config.R"),
    list.files(file.path(root, "sensitivity", "R"), "[.]R$", full.names = TRUE),
    file.path(
      root,
      "sensitivity",
      "reference",
      "PU_PoE_prior_sensitivity_final.R"
    )
  )
  sort(normalizePath(c(old, new), winslash = "/", mustWork = TRUE))
}

new_sensitivity_environment <- function(root, mode = "quick", out_dir = NULL) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  mode <- match.arg(mode, c("quick", "full"))
  check_sensitivity_core(root)
  core <- sensitivity_core(root)
  core$check_illustration_dependencies()
  out <- if (is.null(out_dir)) {
    file.path(root, "results", "prior_sensitivity", mode)
  } else {
    path.expand(out_dir)
  }
  out <- normalizePath(out, winslash = "/", mustWork = FALSE)
  protected <- c(
    root,
    file.path(root, "R"),
    file.path(root, "config"),
    file.path(root, "scripts"),
    file.path(root, "tests"),
    file.path(root, "sensitivity"),
    file.path(root, "reference"),
    file.path(root, "results", "illustration")
  )
  if (
    identical(out, root) ||
      any(vapply(
        protected[-1L],
        function(x) {
          identical(tolower(out), tolower(x)) ||
            startsWith(tolower(out), paste0(tolower(x), "/"))
        },
        logical(1)
      ))
  ) {
    stop(
      "The output folder must be outside the code folders and results/illustration/."
    )
  }
  e <- core$new_illustration_environment(root, mode = mode, out_dir = out)
  e$SENSITIVITY_MODULE_VERSION <- "prior_sensitivity_v1"
  for (name in c("source_functions.R", "checkpoints.R", "preflight.R")) {
    sys.source(
      file.path(root, "sensitivity", "R", name),
      envir = e,
      keep.source = FALSE
    )
  }
  sys.source(
    file.path(root, "sensitivity", "config.R"),
    envir = e,
    keep.source = FALSE
  )
  e$FIT_CHECKPOINT_DIR <- file.path(e$OUT_DIR, "fit_checkpoints")
  e$ENVIRONMENT_FINGERPRINT <- core$illustration_environment_fingerprint()
  files <- sensitivity_files(root)
  e$SENSITIVITY_SOURCE_MD5 <- setNames(
    unname(tools::md5sum(files)),
    substring(files, nchar(root) + 2L)
  )
  e$SENSITIVITY_CONFIGURATION <- list(
    mode = mode,
    common = mget(
      setdiff(e$CONFIGURATION_KEYS, "MODEL_SPECS"),
      envir = e,
      inherits = FALSE
    ),
    common_mcmc_seed = e$COMMON_MCMC_SEED,
    top_k = e$TOP_K_VALUES,
    specifications = e$SENSITIVITY_SPECS
  )
  e
}

sensitivity_save_summaries <- function(e) {
  keys <- c(
    "toy",
    "X_cov",
    "X_standardization",
    "A_network",
    "modularity",
    "Z_network",
    "network_geometry",
    "SENSITIVITY_SPECS",
    "SENSITIVITY_RUN_ORDER",
    "sensitivity_prior_bundles",
    "node_summaries_long",
    "node_summaries_long_all",
    "eta_summary",
    "fit_diagnostics",
    "prior_calibration_completed",
    "sensitivity_fit_summary",
    "topk_overlap_summary",
    "parameter_summary",
    "core_parameter_ess",
    "node_probability_ess",
    "node_ess_summary",
    "common_plot_order",
    "joint_eta_comparison",
    "joint_stress_fit_comparison"
  )
  bundle <- mget(keys, envir = e, inherits = FALSE)
  bundle$run_metadata <- list(
    module = e$SENSITIVITY_MODULE_VERSION,
    mode = e$RUN_MODE,
    created_at = as.character(Sys.time()),
    configuration = e$SENSITIVITY_CONFIGURATION,
    source_files_md5 = e$SENSITIVITY_SOURCE_MD5,
    environment = e$ENVIRONMENT_FINGERPRINT,
    distinct_fits = 13L,
    figure_panels = 18L,
    probability = "latent sigma(f_i)",
    interval = "80% HPD from coda::HPDinterval",
    fit_runtime_seconds = vapply(
      e$sensitivity_results,
      function(x) x$fit$runtime_sec,
      numeric(1)
    )
  )
  e$sensitivity_save_rds(
    bundle,
    file.path(e$OUT_DIR, "sensitivity_summaries.rds")
  )
  e$sensitivity_save_rds(
    bundle$run_metadata,
    file.path(e$OUT_DIR, "completed_run_metadata.rds")
  )
}

sensitivity_render <- function(root, e) {
  existing <- grDevices::dev.list()
  on.exit(
    {
      now <- setdiff(grDevices::dev.list(), existing)
      for (d in now) {
        try(grDevices::dev.off(d), silent = TRUE)
      }
    },
    add = TRUE
  )
  sys.source(file.path(root, "sensitivity", "R", "plot.R"), envir = e)
  utils::write.csv(
    data.frame(mode = e$RUN_MODE, distinct_fits = 13L, panels = 18L),
    file.path(e$OUT_DIR, "figure_run_mode.csv"),
    row.names = FALSE
  )
  cat("Saved sensitivity figure:\n", e$out_pdf_6x3, "\n", sep = "")
}

run_prior_sensitivity <- function(
  root,
  mode = "quick",
  action = "all",
  out_dir = NULL,
  fit_ids = NULL
) {
  mode <- match.arg(mode, c("quick", "full"))
  action <- match.arg(
    action,
    c("check", "prepare", "fit", "plot", "diagnostics", "all")
  )
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  e <- new_sensitivity_environment(root, mode, out_dir)
  core <- sensitivity_core(root)
  parse_files <- c(
    sensitivity_files(root),
    file.path(root, "scripts", "02_prior_sensitivity.R")
  )
  invisible(lapply(parse_files, function(p) {
    parse(file = p, keep.source = FALSE)
  }))
  if (is.null(fit_ids)) {
    fit_ids <- e$SENSITIVITY_RUN_ORDER
  }
  if (
    !is.character(fit_ids) ||
      !length(fit_ids) ||
      anyNA(fit_ids) ||
      anyDuplicated(fit_ids) ||
      any(!fit_ids %in% e$SENSITIVITY_RUN_ORDER)
  ) {
    stop("Unknown fit_ids. See the names in sensitivity/config.R.")
  }
  selected <- e$SENSITIVITY_RUN_ORDER[e$SENSITIVITY_RUN_ORDER %in% fit_ids]
  if (
    action %in%
      c("all", "plot", "diagnostics") &&
      !identical(selected, e$SENSITIVITY_RUN_ORDER)
  ) {
    stop(
      "all/plot/diagnostics require all 13 specifications. Use action=fit for a selected subset."
    )
  }
  old_options <- options(stringsAsFactors = FALSE)
  old_rng <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(
    {
      options(old_options)
      do.call(RNGkind, as.list(old_rng))
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  dir.create(e$OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(e$OUT_DIR)) {
    stop("Could not create output directory: ", e$OUT_DIR)
  }
  e$OUT_DIR <- normalizePath(e$OUT_DIR, winslash = "/", mustWork = TRUE)
  e$FIT_CHECKPOINT_DIR <- file.path(e$OUT_DIR, "fit_checkpoints")
  lock <- file.path(e$OUT_DIR, ".run_lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop(
      "Output directory is locked or not writable: ",
      lock,
      ". Remove .run_lock by hand once no run is active in this folder.",
      call. = FALSE
    )
  }
  on.exit(core$release_illustration_lock(lock), add = TRUE)
  writeLines(
    c(paste("pid", Sys.getpid()), as.character(Sys.time())),
    file.path(lock, "owner.txt")
  )
  log <- file.path(e$OUT_DIR, paste0("console_", action, ".log"))
  sink_depth <- sink.number()
  sink(log, split = TRUE)
  on.exit(
    {
      while (sink.number() > sink_depth) {
        sink()
      }
    },
    add = TRUE,
    after = FALSE
  )
  writeLines(
    c(paste("mode", mode), paste("action", action), "status running"),
    file.path(e$OUT_DIR, "LAST_ACTION.txt")
  )
  start <- Sys.time()
  tryCatch(
    {
      cat(
        "\nPRIOR SENSITIVITY ANALYSIS\nMode:",
        mode,
        "\nAction:",
        action,
        "\nOutput:",
        e$OUT_DIR,
        "\n"
      )
      if (mode == "quick") {
        cat(
          "Quick mode: short runs that check the code, not the results of the paper.\n"
        )
      }
      print(data.frame(
        blocked = e$BLOCKED_BURNIN,
        joint = e$JOINT_BURNIN,
        sampling = e$POSTERIOR_ITERATIONS,
        thinning = e$THIN,
        retained = e$N_RETAINED
      ))
      cat(
        "13 prior specifications on the dataset of the illustrative example.\n"
      )
      writeLines(
        capture.output(sessionInfo()),
        file.path(e$OUT_DIR, "sessionInfo_current.txt")
      )
      saveRDS(
        e$ENVIRONMENT_FINGERPRINT,
        file.path(e$OUT_DIR, "environment_current.rds")
      )
      saveRDS(
        e$SENSITIVITY_CONFIGURATION,
        file.path(e$OUT_DIR, "configuration_current.rds")
      )
      utils::write.csv(
        data.frame(
          package = names(e$ENVIRONMENT_FINGERPRINT$packages),
          version = unname(e$ENVIRONMENT_FINGERPRINT$packages)
        ),
        file.path(e$OUT_DIR, "package_versions_current.csv"),
        row.names = FALSE
      )
      if (action != "plot") {
        sys.source(file.path(root, "sensitivity", "R", "prepare.R"), envir = e)
        e$check_sensitivity_inputs()
        saveRDS(
          mget(
            c(
              "toy",
              "X_cov",
              "X_standardization",
              "A_network",
              "modularity",
              "Z_network",
              "network_geometry",
              "sensitivity_prior_bundles"
            ),
            envir = e,
            inherits = FALSE
          ),
          file.path(e$OUT_DIR, "prepared_inputs.rds")
        )
      }
      if (action %in% c("fit", "diagnostics", "all")) {
        dir.create(e$FIT_CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
        e$sensitivity_results <- setNames(
          vector("list", length(selected)),
          selected
        )
        for (id in selected) {
          e$sensitivity_results[[id]] <- e$run_sensitivity_fit(
            id,
            load_only = action == "diagnostics"
          )
          if (id == "DEFAULT") {
            e$compare_sensitivity_default(e$sensitivity_results[[id]])
          }
        }
        if (identical(selected, e$SENSITIVITY_RUN_ORDER)) {
          sys.source(
            file.path(root, "sensitivity", "R", "summarize.R"),
            envir = e
          )
          sensitivity_save_summaries(e)
          cat(
            "\n================ 13-FIT PRIOR SENSITIVITY COMPLETE ================\n"
          )
          cat("Retained draws per fit:", e$N_RETAINED, "\n")
        } else {
          cat(
            "Selected fits completed. The summaries and the figure need all 13 fits.\n"
          )
        }
      }
      if (action == "plot") {
        path <- file.path(e$OUT_DIR, "sensitivity_summaries.rds")
        if (!file.exists(path)) {
          stop("Completed summaries not found: ", path)
        }
        bundle <- readRDS(path)
        if (
          !is.list(bundle) ||
            !identical(bundle$run_metadata$mode, mode) ||
            !identical(bundle$run_metadata$distinct_fits, 13L)
        ) {
          stop("Summary bundle mode/schema mismatch.")
        }
        list2env(bundle[setdiff(names(bundle), "run_metadata")], envir = e)
      }
      if (action %in% c("plot", "all")) {
        sensitivity_render(root, e)
      }
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      writeLines(
        c(
          paste("mode", mode),
          paste("action", action),
          "status completed",
          paste("wall_seconds", format(elapsed, digits = 12))
        ),
        file.path(e$OUT_DIR, "LAST_ACTION.txt")
      )
      cat(
        "\nCompleted action",
        action,
        "in",
        elapsed,
        "seconds.\nLog:",
        log,
        "\n"
      )
      invisible(list(
        mode = mode,
        action = action,
        output_directory = e$OUT_DIR,
        elapsed_seconds = elapsed
      ))
    },
    error = function(err) {
      writeLines(
        c(
          paste("mode", mode),
          paste("action", action),
          "status failed",
          conditionMessage(err)
        ),
        file.path(e$OUT_DIR, "LAST_ACTION.txt")
      )
      cat("\nERROR:", conditionMessage(err), "\n")
      stop(err)
    }
  )
}
