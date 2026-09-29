## Imports functions and a list of settings from the original scripts. The
## scripts are never run as a whole, since their top-level code runs the study.
sar_assignment_name <- function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]])) {
    as.character(x[[2L]])
  } else {
    ""
  }
}

## Find function definitions at script scope, including definitions inside bare
## braces and if/else branches. This only inspects parsed language objects.
## Conditions and blocks are not evaluated, and function bodies, loops, local(),
## source(), quote() and other calls are not entered.
sar_script_function_index <- function(expressions) {
  index <- list()
  walk <- function(x) {
    if (!is.call(x)) {
      return(invisible(NULL))
    }
    head <- x[[1L]]
    if (identical(head, as.name("<-"))) {
      if (length(x) == 3L && is.symbol(x[[2L]])) {
        rhs <- x[[3L]]
        if (is.call(rhs) && identical(rhs[[1L]], as.name("function"))) {
          name <- as.character(x[[2L]])
          index[[name]] <<- c(index[[name]], list(x))
        }
      }
      # Functions defined inside functions are not read.
      return(invisible(NULL))
    }
    if (identical(head, as.name("{"))) {
      if (length(x) > 1L) {
        for (j in seq.int(2L, length(x))) {
          walk(x[[j]])
        }
      }
    } else if (identical(head, as.name("if"))) {
      # Both branches are read; the condition x[[2L]] is skipped.
      if (length(x) >= 3L) {
        walk(x[[3L]])
      }
      if (length(x) == 4L) walk(x[[4L]])
    }
    invisible(NULL)
  }
  for (x in expressions) {
    walk(x)
  }
  index
}

sar_original_function <- function(index, name, group) {
  definitions <- index[[name]]
  sar_assert(
    length(definitions) == 1L,
    paste0(
      "Expected exactly one original script-scope function ",
      group,
      " ",
      name,
      "; found ",
      length(definitions),
      ". Definitions in script-level if/else blocks are allowed; function-local ",
      "definitions and definitions inside other calls are not imported."
    )
  )
  x <- definitions[[1L]]
  sar_assert(
    identical(sar_assignment_name(x), name) &&
      is.call(x[[3L]]) &&
      identical(x[[3L]][[1L]], as.name("function")),
    paste("Bad original function", group, name)
  )
  x
}

sar_lhs_root <- function(x) {
  if (is.symbol(x)) {
    return(as.character(x))
  }
  if (is.call(x) && as.character(x[[1L]]) %in% c("[", "[[", "$")) {
    return(sar_lhs_root(x[[2L]]))
  }
  ""
}

sar_dependencies <- function(groups = sar_groups()) {
  packages <- c("Matrix", "igraph", "RSpectra")
  if ("Competitors" %in% groups) {
    packages <- c(packages, "AdaSampling", "e1071", "INLA")
  }
  if ("Bayesian_linear" %in% groups) {
    packages <- c(packages, "BayesLogit")
  }
  unique(packages)
}

sar_preflight <- function(groups = sar_groups(), attach = TRUE) {
  p <- sar_dependencies(groups)
  # The archived Bayesian linear fits drew the Polya-Gamma variables with BayesLogit.
  # pgdraw gives different draws for the same seed, so it is not used as a fallback.
  absent <- p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) {
    stop(
      "Refitting the 18 methods needs these packages: ",
      paste(absent, collapse = ", "),
      ".\nInstall them with scripts/00_install_dependencies.R and INSTALL_REFIT <- TRUE.\nScript 06, which rebuilds the tables from the archive, does not need them.",
      call. = FALSE
    )
  }
  if (attach) {
    for (pkg in p) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
  if ("INLA" %in% p) {
    INLA::inla.setOption(num.threads = "1:1")
  }
  stats::setNames(
    vapply(
      p,
      function(pkg) as.character(utils::packageVersion(pkg)),
      character(1)
    ),
    p
  )
}

sar_import <- function(
  root,
  group,
  mode = c("quick", "full"),
  package_versions = NULL
) {
  mode <- match.arg(mode)
  mod <- sar_source_root(root)
  manifest <- utils::read.csv(
    file.path(mod, "reference/SOURCE_MANIFEST.csv"),
    stringsAsFactors = FALSE
  )
  row <- manifest[manifest$group == group, , drop = FALSE]
  sar_assert(nrow(row) == 1L, "Unknown source group.")
  src <- file.path(mod, row$reference_file)
  sar_assert(
    identical(sar_md5(src), as.character(row$md5)),
    paste("Original reference source changed:", group)
  )
  expr <- parse(src, keep.source = FALSE)
  nm <- vapply(expr, sar_assignment_name, character(1))
  function_index <- sar_script_function_index(expr)
  get_expr <- function(name) {
    if (name %in% names(function_index)) {
      return(sar_original_function(function_index, name, group))
    }
    k <- which(nm == name)
    sar_assert(
      length(k) == 1L,
      paste("Expected a unique definition:", group, name)
    )
    expr[[k]]
  }
  e <- new.env(parent = parent.env(.GlobalEnv))
  fns <- readLines(
    file.path(mod, "reference", paste0(group, "_functions.txt")),
    warn = FALSE
  )
  settings <- readLines(
    file.path(mod, "reference", paste0(group, "_settings.txt")),
    warn = FALSE
  )
  for (name in fns) {
    # Only the function assignment is evaluated.
    x <- sar_original_function(function_index, name, group)
    eval(x, e)
  }
  # The original linear code read these numbers from a file. The fixed values
  # below are those recorded in the archived runs.
  ## The priors of the Bayesian linear models are evaluated from the original script.

  for (k in seq_along(expr)) {
    x <- expr[[k]]
    if (nm[k] %in% settings) {
      eval(x, e)
    }
    if (
      is.call(x) &&
        identical(x[[1L]], as.name("<-")) &&
        !is.symbol(x[[2L]]) &&
        sar_lhs_root(x[[2L]]) %in% c("TRUE_T", "POSITIVE_SIGNAL_GROUP")
    ) {
      eval(x, e)
    }
  }
  sar_assert(
    e$N_NODES == 400L &&
      e$N_REP == 100L &&
      identical(which(e$TRUE_T == 1L), 356:400) &&
      e$N_HIDDEN == 15L &&
      e$N_OBSERVED_POS == 30L &&
      e$LAMBDA_SAR == 2 &&
      e$SAR_Z_CLIP == 3,
    paste(group, "nonlinear SAR data design changed")
  )
  if (group == "Bayesian_linear") {
    sar_validate_linear_environment(e)
  }
  if (mode == "quick") {
    overrides <- list(
      PILOT_BURNIN = 120L,
      PILOT_TAIL_LENGTH = 40L,
      PRODUCTION_BURNIN = 120L,
      PRODUCTION_SAMPLING = 200L,
      PRODUCTION_THIN = 2L,
      MIN_BLOCKED_WARMUP = 40L,
      COR_ADAPT_START = 20L,
      COR_ADAPT_INTERVAL = 10L,
      COR_ADAPT_WINDOW = 40L,
      EMP_COV_TAIL = 40L,
      ADAPT_START = 20L,
      ADAPT_INTERVAL = 10L,
      PILOT_PROGRESS_EVERY = 40L,
      PRODUCTION_BURNIN_PROGRESS_EVERY = 40L,
      PRODUCTION_SAMPLE_PROGRESS_EVERY = 100L,
      N_PREDICTIVE_DRAWS = 100L,
      ADA_C = 3L,
      GNN_EPOCHS = 10L,
      NNPU_EPOCHS = 30L
    )
    for (name in names(overrides)) {
      if (exists(name, e, inherits = FALSE)) assign(name, overrides[[name]], e)
    }
  }
  e$RESUME_REPLICATES <- TRUE
  e$RESUME_PILOTS <- TRUE
  e$RESUME_PRODUCTION <- TRUE
  e$FORCE_FRESH_RUN <- FALSE
  e$FORCE_FRESH_PILOTS <- FALSE
  e$FORCE_FRESH_PRODUCTION <- FALSE
  e$FAIL_IF_ANY_METHOD_FAILS <- TRUE
  e$FAIL_IF_ANY_MODEL_FAILS <- TRUE
  e$FAIL_IF_PRODUCTION_FAILS <- TRUE
  e$SAVE_GENERATED_GEOMETRIES <- TRUE
  e$SAVE_PILOT_GEOMETRIES <- TRUE
  e$SAVE_TRAINING_TRACES <- TRUE
  e$USE_SAVED_COMPETING_GEOMETRIES <- TRUE
  e$REQUIRE_SAVED_COMPETING_GEOMETRIES <- TRUE
  # Storing these expressions runs nothing.
  e$SAR_GET_EXPR <- get_expr
  e$SAR_GROUP <- group
  e$SAR_MODE <- mode
  schedule_names <- c(
    "N_NODES",
    "N_PILOT",
    "PILOT_BURNIN",
    "PILOT_TAIL_LENGTH",
    "PRODUCTION_BURNIN",
    "PRODUCTION_SAMPLING",
    "PRODUCTION_THIN",
    "N_PREDICTIVE_DRAWS",
    "ADA_C",
    "GNN_EPOCHS",
    "NNPU_EPOCHS",
    "ADAPT_START",
    "ADAPT_INTERVAL",
    "COR_ADAPT_START",
    "COR_ADAPT_WINDOW",
    "EMP_COV_TAIL"
  )
  schedule_names <- schedule_names[vapply(
    schedule_names,
    exists,
    logical(1),
    envir = e,
    inherits = FALSE
  )]
  e$SAR_SIGNATURE <- list(
    version = "nonlinear_SAR_all18_v1",
    group = group,
    mode = mode,
    source_files = sar_source_signature(root),
    R = R.version.string,
    platform = R.version$platform,
    numerical_libraries = extSoftVersion(),
    packages = package_versions,
    computational_settings = mget(schedule_names, envir = e, inherits = FALSE),
    linear_priors = if (group == "Bayesian_linear") {
      sar_linear_priors()
    } else {
      NULL
    },
    RNG = c("Mersenne-Twister", "Inversion", "Rejection")
  )
  e
}

sar_configure_engine <- function(e, out, shared_geometry_dir) {
  group <- e$SAR_GROUP
  gs <- sar_group_spec()
  prefix <- gs$prefix[gs$group == group]
  e$OUT_DIR <- out
  e$CHECKPOINT_DIR <- file.path(out, "replicate_checkpoints")
  e$GEOMETRY_DIR <- file.path(out, "generated_geometries")
  e$PILOT_CHECKPOINT_DIR <- file.path(out, "pilot_checkpoints")
  e$PRODUCTION_CHECKPOINT_DIR <- file.path(out, "production_checkpoints")
  e$PILOT_GEOMETRY_DIR <- file.path(out, "pilot_geometries")
  e$PRIOR_CALIBRATION_DIR <- file.path(out, "prior_calibrations")
  e$COMPETING_GEOMETRY_DIR <- shared_geometry_dir
  e$COMPETING_OUT_DIR <- dirname(shared_geometry_dir)
  global_stem <- if (group == "Bayesian_linear") {
    "global_pilot_calibration_5_datasets"
  } else {
    paste0(prefix, "_global_pilot_calibration_5_datasets")
  }
  e$GLOBAL_CALIBRATION_FILE <- file.path(out, paste0(global_stem, ".rds"))
  e$GLOBAL_CALIBRATION_RDATA <- file.path(out, paste0(global_stem, ".RData"))
  for (d in c(
    out,
    e$CHECKPOINT_DIR,
    e$GEOMETRY_DIR,
    e$PILOT_CHECKPOINT_DIR,
    e$PRODUCTION_CHECKPOINT_DIR,
    e$PILOT_GEOMETRY_DIR,
    e$PRIOR_CALIBRATION_DIR
  )) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  p <- file.path(out, "ALL18_RUN_SIGNATURE.rds")
  if (file.exists(p)) {
    sar_assert(
      identical(sar_read(p), e$SAR_SIGNATURE),
      paste(
        "Existing",
        group,
        "output uses different code, settings or environment. Use a new output folder."
      )
    )
  } else {
    sar_assert(
      !length(list.files(out, pattern = "[.]rds$", recursive = TRUE)),
      paste("Unknown results in", out)
    )
    sar_atomic_rds(e$SAR_SIGNATURE, p)
  }
  if (group %in% c("Competitors", "nnPU")) {
    eval(e$SAR_GET_EXPR("SETTINGS_SIGNATURE"), e)
    e$SETTINGS_SIGNATURE$all18_execution <- e$SAR_SIGNATURE
  } else {
    eval(e$SAR_GET_EXPR("PILOT_SETTINGS_SIGNATURE"), e)
    e$PILOT_SETTINGS_SIGNATURE$all18_execution <- e$SAR_SIGNATURE
  }
  invisible(e)
}

sar_run_group <- function(e, reps) {
  g <- e$SAR_GROUP
  cat("\n============== ", g, " / ", e$SAR_MODE, " ==============\n", sep = "")
  print(e$SAR_SIGNATURE$computational_settings)
  if (!g %in% c("Competitors", "nnPU")) {
    e$global_pilot_calibration <- e$load_or_build_global_calibration()
    e$GLOBAL_CALIBRATION_HASH <- sar_md5(e$GLOBAL_CALIBRATION_FILE)
    eval(e$SAR_GET_EXPR("PRODUCTION_SETTINGS_SIGNATURE"), e)
    e$PRODUCTION_SETTINGS_SIGNATURE$all18_execution <- e$SAR_SIGNATURE
  }
  out <- lapply(reps, function(r) {
    z <- if (g %in% c("Competitors", "nnPU")) {
      e$run_one_replication(as.integer(r))
    } else {
      e$run_one_production_replication(as.integer(r))
    }
    sar_assert(isTRUE(z$complete), paste(g, r, "did not complete"))
    z
  })
  # One replication at a time keeps the method groups aligned.
  pp <- list()
  i <- 0L
  for (j in seq_along(reps)) {
    for (v in sar_unpack_methods(out[[j]])) {
      i <- i + 1L
      pp[[i]] <- sar_validate_method(
        v,
        reps[j],
        g,
        expected_y = sar_read(file.path(
          e$COMPETING_GEOMETRY_DIR,
          sprintf("replicate_%03d_geometry.rds", reps[j])
        ))$Y
      )
    }
  }
  selections <- lapply(seq_along(reps), function(j) {
    y <- sar_read(file.path(
      e$COMPETING_GEOMETRY_DIR,
      sprintf("replicate_%03d_geometry.rds", reps[j])
    ))$Y
    z <- sar_validate_selection(out[[j]], reps[j], as.integer(y))
    z$group <- g
    z
  })
  list(
    metrics = sar_bind(lapply(pp, `[[`, "metrics")),
    scores = sar_bind(lapply(pp, `[[`, "scores")),
    selection = sar_bind(selections)
  )
}

sar_fit_all <- function(root, out, mode = c("quick", "full")) {
  mode <- match.arg(mode)
  versions <- sar_preflight(attach = TRUE)
  reps <- if (mode == "quick") 1:2 else 1:100
  oldrng <- RNGkind()
  had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  seed <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(
    {
      do.call(RNGkind, as.list(oldrng))
      if (had) {
        assign(".Random.seed", seed, .GlobalEnv)
      } else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  shared <- file.path(out, "groups", "Competitors", "generated_geometries")
  engines <- setNames(
    lapply(sar_groups(), function(g) sar_import(root, g, mode, versions)),
    sar_groups()
  )
  results <- list()
  for (g in sar_groups()) {
    e <- engines[[g]]
    sar_configure_engine(e, file.path(out, "groups", g), shared)
    results[[g]] <- sar_run_group(e, reps)
    rm(e)
    invisible(gc())
  }
  m <- sar_bind(lapply(results, `[[`, "metrics"))
  d <- sar_bind(lapply(results, `[[`, "scores"))
  sar_assert(
    nrow(m) == 18L * length(reps) && nrow(d) == 18L * length(reps) * 400L,
    "Some results of the 18 methods are missing."
  )
  cat("\nAll 18 methods completed on ", length(reps), " datasets.\n", sep = "")
  list(
    metrics = m,
    scores = d,
    selection = sar_bind(lapply(results, `[[`, "selection")),
    mode = if (mode == "quick") "quick" else "new_full",
    replications = length(reps)
  )
}
