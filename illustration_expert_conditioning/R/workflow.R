## Expert scores from the saved GP-PU-PoE fit of the illustrative example.
ec_assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}

ec_normal <- function(p, must_work = TRUE) {
  normalizePath(path.expand(p), winslash = "/", mustWork = must_work)
}

ec_inside <- function(path, directory) {
  path <- tolower(ec_normal(path, FALSE))
  directory <- tolower(ec_normal(directory, FALSE))
  identical(path, directory) || startsWith(path, paste0(directory, "/"))
}

ec_fingerprints <- function(paths, labels = names(paths)) {
  setNames(unname(tools::md5sum(paths)), labels)
}

ec_dependencies <- function() {
  pkgs <- c("coda", "ggplot2", "patchwork")
  absent <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  ec_assert(
    !length(absent),
    paste("Missing packages:", paste(absent, collapse = ", "))
  )
  ec_assert(
    utils::packageVersion("ggplot2") >= "3.4.0",
    "ggplot2 3.4.0 or later is required."
  )
}

ec_read_inputs <- function(root, mode = "full", input_dir = NULL) {
  mode <- match.arg(mode, c("full", "quick"))
  if (is.null(input_dir)) {
    input_dir <- if (mode == "full") {
      file.path(root, "archive", "illustration")
    } else {
      file.path(root, "results", "illustration", mode)
    }
  }
  input_dir <- ec_normal(input_dir, FALSE)
  relative <- c(
    "prepared_inputs.rds",
    "model_checkpoints/GP_PU_POE_complete.rds",
    "toy_seed53.rds"
  )
  paths <- setNames(file.path(input_dir, relative), relative)
  absent <- paths[!file.exists(paths)]
  ec_assert(
    !length(absent),
    paste0(
      "Input files are missing:\n",
      paste(absent, collapse = "\n"),
      "\nSet INPUT_DIR to the folder of a ",
      mode,
      " run of script 01."
    )
  )
  ec_assert(
    !file.exists(file.path(input_dir, ".run_lock")),
    paste0("The input folder is locked by a running fit: ", input_dir)
  )
  order_file <- file.path(
    input_dir,
    "toy_seed53_common_plot_order_by_PU_PoE_latent_mean.csv"
  )
  if (file.exists(order_file)) {
    paths <- c(paths, setNames(order_file, basename(order_file)))
  }
  before <- ec_fingerprints(paths)
  prep <- readRDS(paths[["prepared_inputs.rds"]])
  cp <- readRDS(paths[["model_checkpoints/GP_PU_POE_complete.rds"]])
  toy <- readRDS(paths[["toy_seed53.rds"]])
  ec_assert(
    isTRUE(cp$complete) &&
      identical(cp$model_spec$model_id, "GP_PU_POE") &&
      isTRUE(cp$model_spec$use_pu) &&
      identical(
        as.character(cp$model_spec$active_experts),
        c("cov", "network")
      ),
    "Expected a completed GP-PU-PoE fit with covariate and network experts."
  )
  sig <- cp$model_signature
  cfg <- sig$configuration
  ## MODE "full" also reads a new run of script 01 with another MCMC seed or length.
  allowed <- if (mode == "full") c("full", "custom") else mode
  ec_assert(
    length(sig$execution_mode) == 1L &&
      sig$execution_mode %in% allowed &&
      is.list(cfg),
    "The checkpoint mode does not match MODE, or its saved configuration is missing."
  )
  ec_assert(
    isTRUE(all.equal(prep$toy, toy, tolerance = 0)) &&
      identical(
        unname(sig$original$toy_data_hash),
        unname(before["toy_seed53.rds"])
      ),
    "The prepared data, saved toy dataset and GP-PU-PoE checkpoint do not match."
  )
  ec_assert(
    identical(cp$prior_bundle, prep$model_prior_bundles$GP_PU_POE),
    "The prepared prior bundle does not match the completed fit."
  )
  d <- cp$fit$draws
  P <- d$latent_probability
  ec_assert(
    is.matrix(P) &&
      identical(dim(P), as.integer(c(cfg$N_RETAINED, cfg$N_NODES))) &&
      nrow(P) >= 2L &&
      ncol(P) == 100L &&
      cfg$TOY_SEED == 53L &&
      all(is.finite(P)) &&
      all(P >= 0 & P <= 1),
    "Invalid probability draws, draw count or toy design in the saved fit."
  )
  ec_assert(
    is.matrix(d$theta) &&
      is.matrix(d$tau) &&
      nrow(d$theta) == nrow(P) &&
      nrow(d$tau) == nrow(P) &&
      length(d$beta0) == nrow(P) &&
      all(
        c("ell_cov", "sigma_cov", "ell_network", "sigma_network") %in%
          colnames(d$theta)
      ) &&
      all(c("tau_cov", "tau_network") %in% colnames(d$tau)) &&
      all(is.finite(d$theta)) &&
      all(d$theta > 0) &&
      all(is.finite(d$tau)) &&
      all(d$tau > 0) &&
      all(is.finite(d$beta0)),
    "The saved kernel, tau or intercept draws are incomplete or invalid."
  )
  ec_assert(
    is.finite(cfg$JITTER) &&
      cfg$JITTER >= 0 &&
      cfg$HPD_PROBABILITY == 0.80 &&
      nrow(prep$X_cov) == ncol(P) &&
      nrow(prep$Z_network) == ncol(P) &&
      all(is.finite(prep$X_cov)) &&
      all(is.finite(prep$Z_network)) &&
      identical(as.integer(prep$Y), as.integer(toy$Y)) &&
      identical(as.integer(prep$T_TRUE), as.integer(toy$T)) &&
      isTRUE(all.equal(prep$network_geometry$A, prep$A_network, tolerance = 0)),
    "Saved geometry, labels or numerical settings are inconsistent."
  )
  ns <- cp$node_summary
  ns <- ns[order(ns$node), , drop = FALSE]
  ec_assert(
    identical(as.integer(ns$node), seq_len(100L)) &&
      identical(as.integer(ns$T), as.integer(prep$T_TRUE)) &&
      identical(as.integer(ns$Y), as.integer(prep$Y)) &&
      identical(as.character(ns$node_group), as.character(prep$NODE_GROUP)) &&
      isTRUE(all.equal(
        as.numeric(ns$latent_probability_mean),
        as.numeric(colMeans(P)),
        tolerance = 1e-12
      )),
    "Node summaries do not correspond to the saved GP-PU-PoE draws."
  )
  groups <- match(
    ns$node_group,
    c("true zero", "hidden positive", "observed positive")
  )
  ec_assert(
    !anyNA(groups) &&
      identical(as.integer(tabulate(groups, 3L)), c(70L, 10L, 20L)),
    "The toy plot must have the original 70/10/20 groups."
  )
  order_nodes <- ns$node[order(groups, ns$latent_probability_mean, ns$node)]
  if (file.exists(order_file)) {
    saved_order <- utils::read.csv(order_file, stringsAsFactors = FALSE)
    ec_assert(
      identical(as.integer(saved_order$original_node), as.integer(order_nodes)),
      "The recorded illustration ordering differs from the checkpoint ordering."
    )
  }
  helpers <- c(
    "R/illustration/helpers.R",
    "R/illustration/geometry.R",
    "R/illustration/summaries.R"
  )
  helper_paths <- setNames(file.path(root, helpers), helpers)
  expected <- sig$source_files_md5[helpers]
  actual <- ec_fingerprints(helper_paths)
  ec_assert(
    !anyNA(expected) &&
      !anyNA(actual) &&
      identical(unname(expected), unname(actual)),
    "The illustration helper/kernel/HPD source differs from that used by the saved fit."
  )
  ec_assert(
    identical(before, ec_fingerprints(paths)),
    "An illustration input changed while being read."
  )
  list(
    root = root,
    mode = mode,
    input_dir = input_dir,
    paths = paths,
    fingerprints = before,
    helper_paths = helper_paths,
    helper_fingerprints = actual,
    prepared = prep,
    checkpoint = cp,
    config = cfg,
    order_nodes = as.integer(order_nodes),
    tail_counts = c(zero = sum(P == 0), one = sum(P == 1), entries = length(P))
  )
}

ec_environment <- function(inputs) {
  e <- new.env(parent = as.environment("package:stats"))
  for (p in inputs$helper_paths) {
    sys.source(p, envir = e)
  }
  e$JITTER <- inputs$config$JITTER
  e$HPD_PROBABILITY <- inputs$config$HPD_PROBABILITY
  e$N_NODES <- inputs$config$N_NODES
  e$NODE_GROUP <- inputs$prepared$NODE_GROUP
  e$T_TRUE <- inputs$prepared$T_TRUE
  e$Y <- inputs$prepared$Y
  e$toy_common_plot_order <- inputs$order_nodes
  e$cols_toy_HPD80 <- c(
    "Covariates model" = "dodgerblue1",
    "Network model" = "firebrick2"
  )
  e$toy_block_lines_HPD80 <- data.frame(x = c(70.5, 80.5))
  for (file in c("compute.R", "plot.R")) {
    sys.source(
      file.path(inputs$root, "illustration_expert_conditioning/R", file),
      envir = e
    )
  }
  e
}

ec_with_rng <- function(code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit(
    {
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(456L)
  force(code)
}

ec_compute <- function(inputs, e) {
  latents <- ec_with_rng({
    ## draw_toy_experts() reads .Random.seed from its own environment.
    e$.Random.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    e$draw_toy_experts(
      inputs$checkpoint$fit,
      inputs$prepared$X_cov,
      inputs$prepared$Z_network,
      inputs$prepared$network_geometry
    )
  })
  probabilities <- lapply(latents, stats::plogis)
  summary <- do.call(
    rbind,
    lapply(names(probabilities), function(nm) {
      P <- probabilities[[nm]]
      h <- e$hpd_matrix_coda(P, probability = e$HPD_PROBABILITY)
      data.frame(
        expert = nm,
        node = seq_len(ncol(P)),
        node_group = e$NODE_GROUP,
        T = e$T_TRUE,
        Y = e$Y,
        mean = colMeans(P),
        lower = h[, "lower"],
        upper = h[, "upper"],
        plot_id = match(seq_len(ncol(P)), e$toy_common_plot_order),
        stringsAsFactors = FALSE
      )
    })
  )
  row.names(summary) <- NULL
  summary
}

ec_release_lock <- function(lock) {
  for (i in seq_len(10L)) {
    if (!file.exists(lock)) {
      return(invisible(TRUE))
    }
    unlink(lock, recursive = TRUE, force = TRUE, expand = FALSE)
    if (!file.exists(lock)) {
      return(invisible(TRUE))
    }
    Sys.sleep(0.2)
  }
  warning("Could not remove this run's own lock: ", lock, call. = FALSE)
  invisible(FALSE)
}

run_expert_conditioning <- function(
  root,
  mode = "full",
  action = "all",
  input_dir = NULL,
  out_dir = NULL
) {
  mode <- match.arg(mode, c("full", "quick"))
  action <- match.arg(action, c("all", "plot", "check"))
  root <- ec_normal(root)
  ec_dependencies()
  inputs <- ec_read_inputs(root, mode, input_dir)
  e <- ec_environment(inputs)
  cat(
    "\nILLUSTRATION EXPERT CONDITIONING\nMode:",
    mode,
    "\nAction:",
    action,
    "\nGP-PU-PoE draws:",
    nrow(inputs$checkpoint$fit$draws$latent_probability),
    "\n"
  )
  cat(
    "Inputs checked: saved fit, data, kernels, HPD function and plot order.\n"
  )
  if (action == "check") {
    return(invisible(TRUE))
  }
  dedicated <- file.path(root, "results/illustration/expert_conditioning")
  if (is.null(out_dir)) {
    archived <- ec_inside(
      inputs$input_dir,
      file.path(root, "archive", "illustration")
    )
    out_dir <- if (mode == "quick") {
      file.path(dedicated, "quick")
    } else if (archived) {
      dedicated
    } else {
      file.path(dedicated, basename(inputs$input_dir))
    }
  }
  out_dir <- ec_normal(out_dir, FALSE)
  ec_assert(
    !ec_inside(out_dir, inputs$input_dir) &&
      !ec_inside(inputs$input_dir, out_dir),
    "The output folder must differ from the input folder."
  )
  ec_assert(
    !ec_inside(out_dir, root) || ec_inside(out_dir, dedicated),
    "Inside the project, the output must go to results/illustration/expert_conditioning/."
  )
  compute_file <- file.path(
    root,
    "illustration_expert_conditioning/R/compute.R"
  )
  binding <- list(
    version = "expert_conditioning_1",
    mode = mode,
    input_md5 = inputs$fingerprints,
    helper_md5 = inputs$helper_fingerprints,
    compute_md5 = unname(tools::md5sum(compute_file)),
    conditioning_seed = 456L,
    interval_probability = 0.80,
    jitter = inputs$config$JITTER
  )
  cache <- file.path(out_dir, "expert_conditioning_results.rds")
  if (action == "plot") {
    ec_assert(
      file.exists(cache),
      "No saved expert scores in this folder. Run ACTION = 'all' first."
    )
    result <- readRDS(cache)
    ec_assert(
      identical(result$binding, binding),
      "The saved expert scores do not match these inputs. Run ACTION = 'all'."
    )
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ec_assert(dir.exists(out_dir), paste("Cannot create output folder:", out_dir))
  lock <- file.path(out_dir, ".run_lock")
  ec_assert(
    dir.create(lock, showWarnings = FALSE),
    paste("The output folder is locked:", lock)
  )
  on.exit(ec_release_lock(lock), add = TRUE)
  writeLines(
    c(paste("pid", Sys.getpid()), as.character(Sys.time())),
    file.path(lock, "owner.txt")
  )
  depth <- sink.number()
  sink(file.path(out_dir, paste0("console_", action, ".log")), split = TRUE)
  on.exit(
    {
      while (sink.number() > depth) {
        sink()
      }
    },
    add = TRUE,
    after = FALSE
  )
  cat(
    "ILLUSTRATION EXPERT CONDITIONING\nMode:",
    mode,
    "\nAction:",
    action,
    "\nInput:",
    inputs$input_dir,
    "\nOutput:",
    out_dir,
    "\n"
  )
  if (mode == "quick") {
    cat("Quick mode: these scores only check the code.\n")
  }
  cat(
    "Stored probabilities at 0 / 1:",
    inputs$tail_counts["zero"],
    "/",
    inputs$tail_counts["one"],
    "out of",
    inputs$tail_counts["entries"],
    "entries.\n"
  )
  if (action == "all") {
    summary <- withCallingHandlers(
      ec_compute(inputs, e),
      message = function(m) {
        cat(conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    result <- list(
      binding = binding,
      summary = summary,
      plot_order = inputs$order_nodes,
      tail_counts = inputs$tail_counts,
      retained_draws = nrow(inputs$checkpoint$fit$draws$latent_probability),
      quantities = c("sigma(x_0)", "sigma(x_1)"),
      interpretation = "Marginal posterior expert scores. The draws of the two experts are not a joint sample.",
      numerical_note = "f recovered by qlogis; exact stored 0/1 entries use finite tail limits."
    )
  }
  summary <- result$summary
  ec_assert(
    nrow(summary) == 200L &&
      all(is.finite(as.matrix(summary[c("mean", "lower", "upper")]))) &&
      all(
        summary$mean >= 0 &
          summary$mean <= 1 &
          summary$lower >= 0 &
          summary$upper <= 1 &
          summary$lower <= summary$upper
      ),
    "Invalid expert posterior summaries."
  )
  ec_assert(
    identical(inputs$fingerprints, ec_fingerprints(inputs$paths)) &&
      identical(
        inputs$helper_fingerprints,
        ec_fingerprints(inputs$helper_paths)
      ),
    "An input file changed while the scores were computed."
  )
  e$expert_summary_toy <- summary
  fig <- e$make_expert_figure()
  stem <- if (mode == "full") {
    "HPD_1x2_EXPERTS_GP_PU_POE_TOY_SEED53_80"
  } else {
    "HPD_1x2_EXPERTS_GP_PU_POE_TOY_SEED53_80_QUICK"
  }
  ggplot2::ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot = fig,
    width = 14,
    height = 5,
    units = "in",
    useDingbats = FALSE
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    plot = fig,
    width = 14,
    height = 5,
    units = "in",
    dpi = 300
  )
  utils::write.csv(
    summary,
    file.path(out_dir, "toy_seed53_expert_probability_summaries_HPD80.csv"),
    row.names = FALSE
  )
  if (action == "all") {
    saveRDS(result, cache)
  }
  if (interactive()) {
    print(fig)
  } # Rscript would write Rplots.pdf
  cat(
    "Done: two experts, 100 units, original order.\n",
    "================ EXPERT CONDITIONING COMPLETE ================\n",
    "Saved in:\n",
    out_dir,
    "\n",
    sep = ""
  )
  invisible(list(output_directory = out_dir, summary = summary))
}
