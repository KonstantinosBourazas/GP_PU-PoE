## Entry functions of the illustrative example.

validate_illustration_schedule <- function(
  blocked,
  joint,
  sample,
  thin,
  paper_schedule = FALSE
) {
  x <- c(blocked, joint, sample, thin)
  if (
    length(x) != 4L ||
      any(!is.finite(x)) ||
      any(x < 1) ||
      any(x != floor(x)) ||
      sample %% thin != 0 ||
      sample / thin < 2
  ) {
    stop(
      "MCMC lengths must be positive integers, with >=2 retained draws and sample divisible by thin."
    )
  }
  if (isTRUE(paper_schedule) && !all(x == c(10000L, 20000L, 150000L, 30L))) {
    stop("Full mode uses the schedule of the paper: 10k/20k/150k, thin 30.")
  }
  invisible(TRUE)
}

illustration_dependencies <- function() {
  c("Matrix", "igraph", "RSpectra", "coda", "ggplot2", "patchwork")
}

check_illustration_dependencies <- function() {
  if (getRversion() < "4.0.0") {
    stop("R 4.0 or later is required.")
  }
  pkgs <- illustration_dependencies()
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing packages: ",
      paste(missing, collapse = ", "),
      ". Run scripts/00_install_dependencies.R first."
    )
  }
  if (utils::packageVersion("ggplot2") < "3.4.0") {
    stop("ggplot2 3.4.0 or later is required.")
  }
  invisible(setNames(
    vapply(
      pkgs,
      function(p) as.character(utils::packageVersion(p)),
      character(1)
    ),
    pkgs
  ))
}

illustration_module_files <- function(root) {
  files <- c(
    list.files(
      file.path(root, "R"),
      "[.]R$",
      recursive = TRUE,
      full.names = TRUE
    ),
    list.files(
      file.path(root, "config"),
      "[.]R$",
      recursive = TRUE,
      full.names = TRUE
    )
  )
  sort(normalizePath(files, winslash = "/", mustWork = TRUE))
}

illustration_environment_fingerprint <- function() {
  list(
    R = R.version.string,
    platform = R.version$platform,
    packages = check_illustration_dependencies(),
    external_libraries = extSoftVersion(),
    numerical_libraries = list(
      BLAS = sessionInfo()$BLAS,
      LAPACK = sessionInfo()$LAPACK
    ),
    RNG = c("Mersenne-Twister", "Inversion", "Rejection"),
    thread_environment = Sys.getenv(c(
      "OMP_NUM_THREADS",
      "OPENBLAS_NUM_THREADS",
      "MKL_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS"
    ))
  )
}

## Settings of a new run. The defaults are those of the paper.
paper_run_settings <- function() {
  list(mcmc_seed = 53L, iterations = 150000L, thin = 30L)
}

check_run_settings <- function(settings) {
  whole <- function(x, upper) {
    length(x) == 1L &&
      is.numeric(x) &&
      is.finite(x) &&
      x >= 1 &&
      x <= upper &&
      x == floor(x)
  }
  if (
    !is.list(settings) ||
      !whole(settings$mcmc_seed, 200000) ||
      !whole(settings$iterations, 1e8) ||
      !whole(settings$thin, 1e6) ||
      settings$iterations %% settings$thin != 0 ||
      settings$iterations %/% settings$thin < 2
  ) {
    stop(
      "MCMC_SEED, ITERATIONS and THIN must be positive whole numbers, ",
      "ITERATIONS must be a multiple of THIN and keep at least two draws, ",
      "and MCMC_SEED must be at most 200000.",
      call. = FALSE
    )
  }
  lapply(settings[c("mcmc_seed", "iterations", "thin")], as.integer)
}

custom_run_folder <- function(settings) {
  sprintf(
    "custom_seed%d_iter%d_thin%d",
    settings$mcmc_seed,
    settings$iterations,
    settings$thin
  )
}

new_illustration_environment <- function(
  root,
  mode = "quick",
  out_dir = NULL,
  fresh = FALSE,
  settings = NULL
) {
  mode <- match.arg(mode, c("quick", "full", "custom"))
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, ".here"))) {
    stop("Not the project folder: ", root)
  }
  ## Avoid inheriting unrelated objects from the user's global workspace.
  e <- new.env(parent = as.environment("package:stats"))
  e$RUN_MODE <- mode
  e$MODULE_VERSION <- "illustration_v1"
  e$PROJECT_ROOT <- root
  e$validate_illustration_schedule <- validate_illustration_schedule
  if (identical(mode, "custom")) {
    settings <- check_run_settings(settings)
    e$CUSTOM_MCMC_SEED <- settings$mcmc_seed
    e$CUSTOM_ITERATIONS <- settings$iterations
    e$CUSTOM_THIN <- settings$thin
  }
  e$OUT_DIR <- if (is.null(out_dir)) {
    file.path(
      root,
      "results",
      "illustration",
      if (identical(mode, "custom")) custom_run_folder(settings) else mode
    )
  } else {
    path.expand(out_dir)
  }
  e$OUT_DIR <- normalizePath(e$OUT_DIR, winslash = "/", mustWork = FALSE)
  e$MODEL_CHECKPOINT_DIR <- file.path(e$OUT_DIR, "model_checkpoints")
  e$RESUME_COMPLETED_MODELS <- TRUE
  e$FORCE_FRESH_RUN <- isTRUE(fresh)
  for (f in c(
    "helpers.R",
    "geometry.R",
    "priors.R",
    "sampler.R",
    "summaries.R"
  )) {
    sys.source(
      file.path(root, "R", "illustration", f),
      envir = e,
      keep.source = FALSE
    )
  }
  sys.source(
    file.path(root, "config", "illustration.R"),
    envir = e,
    keep.source = FALSE
  )
  e$CONFIGURATION_SNAPSHOT <- mget(
    e$CONFIGURATION_KEYS,
    envir = e,
    inherits = FALSE
  )
  source_files <- illustration_module_files(root)
  source_names <- substring(source_files, nchar(root) + 2L)
  e$SOURCE_FILE_MD5 <- setNames(
    unname(tools::md5sum(source_files)),
    source_names
  )
  e$ORIGINAL_SOURCE_SHA256 <- readLines(
    file.path(root, "docs", "ORIGINAL_SOURCE_SHA256.txt"),
    warn = FALSE
  )[1L]
  sys.source(
    file.path(root, "R", "illustration", "cache.R"),
    envir = e,
    keep.source = FALSE
  )
  e
}

## Removes the lock of this run. It is called only after this run has created
## the lock, so the lock of another run is never removed.
release_illustration_lock <- function(lock, attempts = 10L, delay = 0.2) {
  if (
    !is.character(lock) ||
      length(lock) != 1L ||
      is.na(lock) ||
      !identical(basename(lock), ".run_lock")
  ) {
    stop("Lock cleanup is restricted to a path named .run_lock.")
  }
  if (
    length(attempts) != 1L ||
      !is.finite(attempts) ||
      attempts < 1 ||
      attempts != floor(attempts) ||
      length(delay) != 1L ||
      !is.finite(delay) ||
      delay < 0
  ) {
    stop("Invalid lock-cleanup retry settings.")
  }
  lock_exists <- function() dir.exists(lock) || file.exists(lock)
  last_problem <- NULL
  for (i in seq_len(as.integer(attempts))) {
    if (!lock_exists()) {
      return(invisible(TRUE))
    }
    ## Do not interpret '[' or other path characters as wildcard syntax.
    ## force = TRUE also removes a read-only lock.
    status <- tryCatch(
      withCallingHandlers(
        unlink(lock, recursive = TRUE, force = TRUE, expand = FALSE),
        warning = function(w) {
          last_problem <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      ),
      error = function(err) {
        last_problem <<- conditionMessage(err)
        1L
      }
    )
    if (!lock_exists()) {
      return(invisible(TRUE))
    }
    if (i < attempts && delay > 0) Sys.sleep(delay)
  }
  warning(
    "Could not remove this run's lock after ",
    attempts,
    " attempts: ",
    lock,
    if (!is.null(last_problem)) {
      paste0(" (", last_problem, ")")
    } else {
      paste0(" (unlink status ", status, ")")
    },
    ". Remove .run_lock by hand once no run is active in this folder.",
    call. = FALSE
  )
  invisible(FALSE)
}

## The saved results of the paper: the full fits behind the figure and tables.
SAVED_BUNDLE_FILE <- "TOY_SEED53_FOUR_GP_MODELS_COMPLETE_RESULTS.RData"

saved_illustration_dir <- function(root) {
  file.path(root, "archive", "illustration")
}

## Loads a completed summary bundle and the settings that belong to it.
load_illustration_bundle <- function(e, bundle, expected_mode) {
  if (!file.exists(bundle)) {
    stop(
      "No completed summary bundle: ",
      bundle,
      ". Run action 'fit' or 'all' first.",
      call. = FALSE
    )
  }
  load(bundle, envir = e)
  if (!identical(e$run_config$execution_mode, expected_mode)) {
    stop("Mode does not match the saved results.", call. = FALSE)
  }
  ## Settings and nodes of the saved result.
  e$N_RETAINED <- e$run_config$MCMC$retained_draws
  e$N_NODES <- e$run_config$n_nodes
  e$Y <- as.integer(e$toy$Y)
  e$T_TRUE <- as.integer(e$toy$T)
  e$NODE_GROUP <- e$node_summaries_wide$node_group
  invisible(e)
}

## Checks that the saved fits are complete and belong to the saved data.
check_saved_illustration <- function(e, input_dir) {
  toy_file <- file.path(input_dir, "toy_seed53.rds")
  if (
    !file.exists(toy_file) ||
      !identical(
        unname(tools::md5sum(toy_file)),
        e$run_config$exact_toy_data_hash
      )
  ) {
    stop(
      "The saved data do not match the saved results in ",
      input_dir,
      call. = FALSE
    )
  }
  for (spec in e$MODEL_SPECS) {
    checkpoint <- file.path(
      input_dir,
      "model_checkpoints",
      paste0(spec$model_id, "_complete.rds")
    )
    if (!file.exists(checkpoint)) {
      stop("Missing saved fit: ", checkpoint, call. = FALSE)
    }
  }
  invisible(TRUE)
}

run_illustration <- function(
  root,
  mode = "quick",
  action = "all",
  out_dir = NULL,
  fresh = FALSE,
  mcmc_seed = 53L,
  iterations = 150000L,
  thin = 30L
) {
  mode <- match.arg(mode, c("saved", "quick", "full"))
  action <- match.arg(
    action,
    c("check", "prepare", "fit", "plot", "diagnostics", "all")
  )
  settings <- check_run_settings(list(
    mcmc_seed = mcmc_seed,
    iterations = iterations,
    thin = thin
  ))
  paper_settings <- identical(settings, paper_run_settings())
  saved <- identical(mode, "saved")
  if (!paper_settings && !identical(mode, "full")) {
    stop(
      "MCMC_SEED, ITERATIONS and THIN apply only to MODE 'full'. ",
      "MODE 'saved' shows the saved results of the paper and MODE 'quick' ",
      "uses its own short schedule.",
      call. = FALSE
    )
  }
  if (saved && action %in% c("prepare", "fit")) {
    stop(
      "MODE 'saved' reads the saved fits and runs no MCMC. ",
      "Use MODE 'full' or 'quick' to fit the models.",
      call. = FALSE
    )
  }
  run_mode <- if (saved) {
    "full"
  } else if (paper_settings) {
    mode
  } else {
    "custom"
  }
  check_illustration_dependencies() # includes the plotting packages
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  parse_files <- c(
    illustration_module_files(root),
    file.path(root, "scripts", "01_illustrative_example.R")
  )
  invisible(lapply(parse_files, function(p) {
    parse(file = p, keep.source = FALSE)
  }))

  if (saved) {
    input_dir <- normalizePath(
      saved_illustration_dir(root),
      winslash = "/",
      mustWork = FALSE
    )
    if (!dir.exists(input_dir)) {
      stop("The saved results are missing: ", input_dir, call. = FALSE)
    }
    if (is.null(out_dir)) {
      out_dir <- file.path(root, "results", "illustration", "saved")
    }
    inside <- function(path, folder) {
      path <- tolower(normalizePath(path, winslash = "/", mustWork = FALSE))
      folder <- tolower(folder)
      identical(path, folder) || startsWith(path, paste0(folder, "/"))
    }
    if (
      inside(
        out_dir,
        tolower(normalizePath(
          file.path(root, "archive"),
          winslash = "/",
          mustWork = FALSE
        ))
      )
    ) {
      stop("The output folder cannot be inside archive/.", call. = FALSE)
    }
  }
  e <- new_illustration_environment(
    root,
    run_mode,
    out_dir,
    fresh,
    settings = if (identical(run_mode, "custom")) settings
  )
  e$ENVIRONMENT_FINGERPRINT <- illustration_environment_fingerprint()
  if (saved) {
    e$MODEL_CHECKPOINT_DIR <- file.path(input_dir, "model_checkpoints")
  }
  old_options <- options(stringsAsFactors = FALSE)
  old_rng <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
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
  if (!saved) {
    e$MODEL_CHECKPOINT_DIR <- file.path(e$OUT_DIR, "model_checkpoints")
  }
  lock <- file.path(e$OUT_DIR, ".run_lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    if (dir.exists(lock) || file.exists(lock)) {
      owner_file <- file.path(lock, "owner.txt")
      owner <- if (file.exists(owner_file)) {
        tryCatch(
          paste(readLines(owner_file, warn = FALSE), collapse = "; "),
          error = function(err) "owner details unavailable"
        )
      } else {
        "owner details unavailable"
      }
      stop(
        "Output directory is locked: ",
        lock,
        " (",
        owner,
        "). ",
        "Remove .run_lock by hand once no run is active in this folder.",
        call. = FALSE
      )
    }
    stop(
      "Could not create the output-directory lock: ",
      lock,
      ". Check directory access and permissions.",
      call. = FALSE
    )
  }
  on.exit(release_illustration_lock(lock), add = TRUE)
  writeLines(
    c(paste("pid", Sys.getpid()), as.character(Sys.time())),
    file.path(lock, "owner.txt")
  )

  log_file <- file.path(e$OUT_DIR, paste0("console_", action, ".log"))
  sink_depth <- sink.number()
  sink(log_file, split = TRUE)
  ## Close this run's log before the other exit handlers release the lock.
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
  cat(
    "\nILLUSTRATIVE EXAMPLE\nMode:",
    mode,
    "\nAction:",
    action,
    "\nOutput:",
    e$OUT_DIR,
    "\n"
  )
  if (saved) {
    cat(
      "Saved mode: the full fits of the paper are read from ",
      input_dir,
      ". No MCMC is run.\n",
      sep = ""
    )
  } else {
    if (mode == "quick") {
      cat(
        "Quick mode: short runs that check the code, not the results of the paper.\n"
      )
    }
    if (identical(run_mode, "custom")) {
      cat(
        "New run with settings other than those of the paper: MCMC seed ",
        settings$mcmc_seed,
        ", ",
        settings$iterations,
        " iterations, thin ",
        settings$thin,
        ".\n",
        sep = ""
      )
    }
    print(data.frame(
      blocked = e$BLOCKED_BURNIN,
      joint = e$JOINT_BURNIN,
      sampling = e$POSTERIOR_ITERATIONS,
      thin = e$THIN,
      retained = e$N_RETAINED
    ))
  }
  writeLines(
    capture.output(sessionInfo()),
    file.path(e$OUT_DIR, "sessionInfo_current.txt")
  )
  saveRDS(
    e$ENVIRONMENT_FINGERPRINT,
    file.path(e$OUT_DIR, "environment_current.rds")
  )
  write.csv(
    data.frame(
      package = names(e$ENVIRONMENT_FINGERPRINT$packages),
      version = unname(e$ENVIRONMENT_FINGERPRINT$packages)
    ),
    file.path(e$OUT_DIR, "package_versions_current.csv"),
    row.names = FALSE
  )
  saveRDS(
    e$CONFIGURATION_SNAPSHOT,
    file.path(e$OUT_DIR, "configuration_current.rds")
  )

  start <- Sys.time()
  tryCatch(
    {
      if (saved) {
        load_illustration_bundle(
          e,
          file.path(input_dir, SAVED_BUNDLE_FILE),
          "full"
        )
        check_saved_illustration(e, input_dir)
        if (action %in% c("diagnostics", "all")) {
          sys.source(
            file.path(root, "R", "illustration", "diagnostics.R"),
            envir = e
          )
        }
        if (action %in% c("plot", "all")) {
          sys.source(file.path(root, "R", "illustration", "plot.R"), envir = e)
        }
        if (action == "check") cat("The saved results are complete.\n")
      } else {
        if (action %in% c("prepare", "fit", "all")) {
          dir.create(
            e$MODEL_CHECKPOINT_DIR,
            recursive = TRUE,
            showWarnings = FALSE
          )
          sys.source(
            file.path(root, "R", "illustration", "prepare.R"),
            envir = e
          )
          prep_keys <- c(
            "toy",
            "X_cov",
            "X_standardization",
            "A_network",
            "Y",
            "T_TRUE",
            "NODE_GROUP",
            "modularity",
            "Z_network",
            "network_geometry",
            "model_prior_bundles",
            "beta0_prior_for_sampler",
            "eta_prior_for_sampler"
          )
          saveRDS(
            mget(prep_keys, envir = e, inherits = FALSE),
            file.path(e$OUT_DIR, "prepared_inputs.rds")
          )
        }
        if (action %in% c("fit", "all")) {
          sys.source(
            file.path(root, "R", "illustration", "fit_and_summarize.R"),
            envir = e
          )
          sys.source(
            file.path(root, "R", "illustration", "diagnostics.R"),
            envir = e
          )
        }
        if (action %in% c("plot", "diagnostics")) {
          load_illustration_bundle(
            e,
            file.path(e$OUT_DIR, SAVED_BUNDLE_FILE),
            run_mode
          )
          if (action == "diagnostics") {
            sys.source(
              file.path(root, "R", "illustration", "diagnostics.R"),
              envir = e
            )
          }
        }
        if (action %in% c("plot", "all")) {
          sys.source(file.path(root, "R", "illustration", "plot.R"), envir = e)
        }
        if (action == "check") cat("Parsing, packages and settings checked.\n")
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
        log_file,
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
