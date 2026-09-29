## Driver of the original runs over 100 datasets, without the machine-specific
## top-level paths and the automatic run_prior_mc() call. The rest is unchanged.
## Used by script 02c only; script 02b does not source it.

run_prior_mc <- function() {
  
  ## ------------------------------------------------------------------------
  ## 1. CHECKS AND IMPORT OF THE ORIGINAL FUNCTIONS
  ## ------------------------------------------------------------------------
  
  if (
    length(MODE) != 1L ||
    !MODE %in% c("run", "combine", "plot")
  ) {
    stop("Invalid MODE.")
  }
  
  old_options <- options(stringsAsFactors = FALSE)
  on.exit(options(old_options), add = TRUE)
  
  if (!file.exists(ORIGINAL_SCRIPT)) {
    stop("Original script not found: ", ORIGINAL_SCRIPT)
  }
  
  packages <- c(
    "Matrix", "igraph", "RSpectra",
    "coda", "ggplot2", "patchwork"
  )
  
  absent <- packages[
    !vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(absent)) {
    stop("Install packages: ", paste(absent, collapse = ", "))
  }
  
  suppressPackageStartupMessages({
    library(Matrix)
    library(igraph)
  })
  
  ## parse() rather than source(), so the original simulation loop is not run.
  expressions <- parse(
    file = ORIGINAL_SCRIPT,
    keep.source = FALSE
  )
  
  assigned_name <- function(x) {
    if (
      is.call(x) &&
      identical(x[[1L]], as.name("<-")) &&
      is.symbol(x[[2L]])
    ) {
      as.character(x[[2L]])
    } else {
      ""
    }
  }
  
  expression_names <- vapply(
    expressions,
    assigned_name,
    character(1)
  )
  
  original_expression <- function(name) {
    k <- which(expression_names == name)
    
    if (length(k) != 1L) {
      stop("Expected one original definition of ", name)
    }
    
    expressions[[k]]
  }
  
  E <- new.env(parent = globalenv())
  
  settings_end <- match("OUT_DIR", expression_names)
  
  if (is.na(settings_end)) {
    stop("The original settings block was not found.")
  }
  
  ## Original constants before its output-directory block.
  for (k in seq_len(settings_end - 1L)) {
    if (grepl("^[A-Z][A-Z0-9_]*$", expression_names[k])) {
      eval(expressions[[k]], E)
    }
  }
  
  ## Original function definitions only.
  for (k in seq_along(expressions)) {
    x <- expressions[[k]]
    
    if (
      nzchar(expression_names[k]) &&
      is.call(x[[3L]]) &&
      identical(x[[3L]][[1L]], as.name("function"))
    ) {
      eval(x, E)
    }
  }
  
  ## Original eleven settings plus the two original joint settings.
  for (
    name in c(
      "SENSITIVITY_SPECS",
      "BASE_PUPOE_MODEL_SPEC",
      "JOINT_STRESS_SPECS"
    )
  ) {
    eval(original_expression(name), E)
  }
  
  specs <- c(
    E$SENSITIVITY_SPECS,
    E$JOINT_STRESS_SPECS
  )
  
  fit_ids <- names(specs)
  
  if (length(fit_ids) != 13L || anyDuplicated(fit_ids)) {
    stop("Expected 13 unique prior settings.")
  }
  
  ## ------------------------------------------------------------------------
  ## 2. FIVE SESSIONS AND COMMON REPLICATION SEEDS
  ## ------------------------------------------------------------------------
  
  batches <- list(
    "DEFAULT",
    
    c(
      "SIGMA_CONSERVATIVE",
      "ELL_CONSERVATIVE",
      "GAMMA_CONSERVATIVE"
    ),
    
    c(
      "BETA0_CONSERVATIVE",
      "ETA_CONSERVATIVE",
      "JOINT_CONSERVATIVE"
    ),
    
    c(
      "SIGMA_DIFFUSE",
      "ELL_DIFFUSE",
      "GAMMA_DIFFUSE"
    ),
    
    c(
      "BETA0_DIFFUSE",
      "ETA_DIFFUSE",
      "JOINT_DIFFUSE"
    )
  )
  
  if (
    !setequal(unlist(batches), fit_ids) ||
    anyDuplicated(unlist(batches))
  ) {
    stop("The five batches do not partition the 13 settings.")
  }
  
  n_rep <- 100L
  
  ## Different seeds across replications.
  ## Within a replication, all settings use the same data and MCMC seeds.
  ## Replication 1 retains the three original seeds.
  seed_plan <- data.frame(
    replicate = seq_len(n_rep),
    
    data_seed = as.integer(
      E$TOY_SEED + 100003L * (0:99)
    ),
    
    modularity_seed = as.integer(
      E$MODULARITY_SEED + 100003L * (0:99)
    ),
    
    mcmc_seed = as.integer(
      E$COMMON_MCMC_SEED + 100003L * (0:99)
    )
  )
  
  ## SESSION_ID and REPLICATIONS are left out, so changing which batch is
  ## run does not invalidate completed fits.
  experiment <- list(
    driver_version = "MC100_HPD80_v1",
    
    source_md5 = unname(
      tools::md5sum(ORIGINAL_SCRIPT)
    ),
    
    package_versions = setNames(
      vapply(
        packages,
        function(p) as.character(packageVersion(p)),
        character(1)
      ),
      packages
    ),
    
    seed_plan = seed_plan,
    specifications = specs
  )
  
  group_levels <- c(
    "true zero",
    "hidden positive",
    "observed positive"
  )
  
  true_eta <- length(E$HIDDEN_POS_IDX) / length(E$TRUE_POS_IDX)
  
  RNGkind(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  
  dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ## ------------------------------------------------------------------------
  ## 3. FILE AND CHECKPOINT HELPERS
  ## ------------------------------------------------------------------------
  
  atomic_write <- function(object, path, csv = FALSE) {
    dir.create(
      dirname(path),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    tmp <- tempfile(
      ".writing_",
      tmpdir = dirname(path)
    )
    
    on.exit(unlink(tmp), add = TRUE)
    
    if (csv) {
      write.csv(object, tmp, row.names = FALSE)
    } else {
      saveRDS(object, tmp)
    }
    
    if (file.exists(path) && !file.remove(path)) {
      stop("Cannot replace ", path)
    }
    
    if (!file.rename(tmp, path)) {
      stop("Cannot finish writing ", path)
    }
    
    invisible(path)
  }
  
  write_table <- function(x, path) {
    atomic_write(x, path, csv = TRUE)
  }
  
  object_hash <- function(x) {
    tmp <- tempfile()
    on.exit(unlink(tmp), add = TRUE)
    
    saveRDS(
      x,
      tmp,
      compress = FALSE,
      version = 2
    )
    
    unname(tools::md5sum(tmp))
  }
  
  checkpoint_path <- function(id, r, raw = FALSE) {
    file.path(
      OUT_DIR,
      if (raw) "draws" else "checkpoints",
      id,
      sprintf("rep_%03d.rds", r)
    )
  }
  
  read_checked <- function(path, id, r) {
    z <- tryCatch(
      readRDS(path),
      error = function(e) {
        stop(
          "Unreadable checkpoint: ", path,
          "\n", conditionMessage(e)
        )
      }
    )
    
    if (
      !isTRUE(z$complete) ||
      !identical(z$experiment, experiment) ||
      !identical(z$fit_id, id) ||
      !identical(z$replicate, as.integer(r))
    ) {
      stop(
        "Incompatible checkpoint: ", path,
        "\nUse the same source/settings or a NEW output folder."
      )
    }
    
    z
  }
  
  bind_rows <- function(rows) {
    z <- do.call(rbind, rows)
    rownames(z) <- NULL
    z
  }
  
  ## ------------------------------------------------------------------------
  ## 4. SAME DATA GENERATOR, WITH A REPLICATION-SPECIFIC SEED
  ## ------------------------------------------------------------------------
  
  generate_toy <- function(seed) {
    n <- 100L
    p <- 10L
    pos <- 71:100
    
    set.seed(seed)
    
    X <- matrix(
      rnorm(n * p, 0, 1),
      n,
      p
    )
    
    X[pos, ] <- rnorm(
      length(pos) * p,
      0.5,
      1
    )
    
    colnames(X) <- paste0("X", seq_len(p))
    
    Pmat <- matrix(0.10, n, n)
    Pmat[pos, ] <- 0.05
    Pmat[, pos] <- 0.05
    Pmat[pos, pos] <- 0.16
    
    A <- matrix(0L, n, n)
    ut <- upper.tri(A)
    
    A[ut] <- as.integer(
      runif(sum(ut)) < Pmat[ut]
    )
    
    A <- A + t(A)
    diag(A) <- 0L
    
    T <- integer(n)
    T[pos] <- 1L
    
    Y <- integer(n)
    Y[81:100] <- 1L
    
    list(
      seed = as.integer(seed),
      N = n,
      X = X,
      A = A,
      T = T,
      Y = Y,
      pos_idx = pos,
      observed_idx = 81:100,
      hidden_idx = 71:80,
      zero_idx = 1:70,
      probability_matrix = Pmat
    )
  }
  
  ## Verify that the first dataset reproduces the original exactly.
  if (
    MODE == "run" &&
    !identical(
      generate_toy(53L),
      E$generate_exact_toy_seed53()
    )
  ) {
    stop(
      "The copied generator does not reproduce the original seed-53 data."
    )
  }
  
  prepare_replication <- function(r) {
    E$TOY_SEED <- seed_plan$data_seed[r]
    E$MODULARITY_SEED <- seed_plan$modularity_seed[r]
    
    E$toy <- generate_toy(E$TOY_SEED)
    E$TOY_DATA_HASH <- object_hash(E$toy)
    
    E$X_standardization <- E$standardize_columns(E$toy$X)
    E$X_cov <- E$X_standardization$scaled
    colnames(E$X_cov) <- colnames(E$toy$X)
    
    E$A_network <- E$sanitize_network_adjacency(E$toy$A)
    E$Y <- as.integer(E$toy$Y)
    E$T_TRUE <- as.integer(E$toy$T)
    
    E$NODE_GROUP <- rep("true zero", E$N_NODES)
    E$NODE_GROUP[E$HIDDEN_POS_IDX] <- "hidden positive"
    E$NODE_GROUP[E$OBSERVED_POS_IDX] <- "observed positive"
    
    E$modularity <- E$modularity_features_single_network(
      A = E$A_network,
      seed = E$MODULARITY_SEED,
      eig_tol = E$MOD_EIG_TOL,
      max_q = E$MOD_MAX_Q,
      k_max = E$MOD_K_MAX
    )
    
    E$Z_network <- as.matrix(E$modularity$Z)
    
    E$network_geometry <- E$prepare_network_components_for_blocks(
      E$A_network,
      network_name = "toy_network"
    )
  }
  
  metadata <- function(df, id, r) {
    df$fit_id <- id
    df$replicate <- as.integer(r)
    df$varying_block <- specs[[id]]$varying_block
    df$prior_role <- specs[[id]]$prior_role
    rownames(df) <- NULL
    df
  }
  
  ## ------------------------------------------------------------------------
  ## 5. RUN OR RESUME ONE CONFIGURATION ON ONE DATASET
  ## ------------------------------------------------------------------------
  
  run_one <- function(id, r) {
    path <- checkpoint_path(id, r)
    
    if (file.exists(path)) {
      z <- read_checked(path, id, r)
      
      if (!identical(z$toy_hash, E$TOY_DATA_HASH)) {
        stop("Dataset mismatch: ", path)
      }
      
      cat(sprintf(
        "Loaded %s, replication %03d\n",
        id,
        r
      ))
      
      return(invisible("loaded"))
    }
    
    dir.create(
      dirname(path),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    ## Prevent two processes from running the same fit simultaneously.
    lock <- paste0(path, ".lock")
    
    if (!dir.create(lock, showWarnings = FALSE)) {
      stop(
        "Fit already locked: ", lock,
        "\nIf the previous R process has stopped, ",
        "remove ONLY this .lock folder and retry."
      )
    }
    
    on.exit(
      unlink(lock, recursive = TRUE),
      add = TRUE
    )
    
    spec <- specs[[id]]
    spec$mcmc_seed <- seed_plan$mcmc_seed[r]
    
    ## Recalibrate from this dataset with the original calibration function.
    prior <- E$calibrate_sensitivity_fit_priors(
      fit_spec = spec,
      X_cov = E$X_cov,
      Z_network = E$Z_network,
      network_geometry = E$network_geometry
    )
    
    model <- E$make_sensitivity_model_spec(spec)
    
    signature <- E$make_sensitivity_signature(
      spec,
      model,
      prior
    )
    
    raw_path <- checkpoint_path(id, r, raw = TRUE)
    
    if (file.exists(raw_path)) {
      raw <- read_checked(raw_path, id, r)
      
      if (!identical(raw$fit_signature, signature)) {
        stop("Raw-fit signature mismatch.")
      }
      
      fit <- raw$fit
      status <- "raw_loaded"
      
    } else {
      cat(sprintf(
        "\nRUN %s | replication %03d/100 | data seed %d\n",
        id,
        r,
        E$TOY_SEED
      ))
      
      set.seed(spec$mcmc_seed)
      
      ## Original sampler, unchanged.
      fit <- E$run_exact_toy_model_chain(
        model_spec = model,
        prior_bundle = prior,
        Y = E$Y,
        X_cov = E$X_cov,
        A_network = E$A_network,
        Z_network = E$Z_network,
        network_geometry = E$network_geometry,
        
        m_beta = prior$beta0_prior$mean,
        s_beta = prior$beta0_prior$sd,
        a_eta = prior$eta_prior$a_eta,
        b_eta = prior$eta_prior$b_eta,
        
        n_blocked_burnin = E$BLOCKED_BURNIN,
        n_joint_burnin = E$JOINT_BURNIN,
        n_sample = E$POSTERIOR_ITERATIONS,
        thin = E$THIN,
        
        max_log_sigma_cov = log(spec$max_sigma_cov),
        max_mean_sd_network = spec$max_mean_sd_network,
        
        verbose = VERBOSE_MCMC
      )
      
      raw <- list(
        complete = TRUE,
        experiment = experiment,
        fit_id = id,
        replicate = as.integer(r),
        toy_hash = E$TOY_DATA_HASH,
        fit_signature = signature,
        fit_spec = spec,
        model_spec = model,
        prior_bundle = prior,
        fit = fit,
        toy = E$toy,
        X_standardization = E$X_standardization,
        modularity = E$modularity,
        session_info = capture.output(sessionInfo())
      )
      
      ## Save the finished chain before the HPD and ESS summaries.
      ## If postprocessing is interrupted, the chain need not run again.
      atomic_write(raw, raw_path)
      
      status <- "computed"
    }
    
    ## Original posterior summaries, including the three probability targets.
    s <- E$summarize_completed_model(
      fit,
      model,
      prior
    )
    
    s$parameter_summary <- E$summarize_parameter_draws(raw)
    
    ## Preserve the original ESS diagnostics.
    core <- cbind(
      beta0 = fit$draws$beta0,
      eta = fit$draws$eta,
      fit$draws$tau,
      fit$draws$theta
    )
    
    s$core_parameter_ess <- data.frame(
      parameter = colnames(core),
      retained_draws = nrow(core),
      ESS = E$effective_size_matrix(core)
    )
    
    s$core_parameter_ess$ESS_fraction_of_retained <-
      s$core_parameter_ess$ESS / nrow(core)
    
    ne <- E$effective_size_matrix(
      fit$draws$latent_probability
    )
    
    s$node_probability_ess <- data.frame(
      node = seq_len(E$N_NODES),
      node_group = E$NODE_GROUP,
      retained_draws = E$N_RETAINED,
      ESS = ne,
      ESS_fraction_of_retained = ne / E$N_RETAINED
    )
    
    s$node_ess_summary <- data.frame(
      retained_draws = E$N_RETAINED,
      min_ESS = min(ne),
      q25_ESS = unname(quantile(ne, 0.25)),
      median_ESS = median(ne),
      mean_ESS = mean(ne),
      q75_ESS = unname(quantile(ne, 0.75)),
      max_ESS = max(ne)
    )
    
    ## Actual replication-level coverage, separate from averaged endpoints.
    s$eta_summary$HPD80_covers_true <- as.integer(
      s$eta_summary$eta_HPD80_lower <= true_eta &
        true_eta <= s$eta_summary$eta_HPD80_upper
    )
    
    s <- lapply(
      s,
      metadata,
      id = id,
      r = r
    )
    
    result <- c(
      list(
        complete = TRUE,
        experiment = experiment,
        fit_id = id,
        replicate = as.integer(r),
        toy_hash = E$TOY_DATA_HASH,
        fit_signature = signature,
        fit_spec = spec,
        eta_prior = prior$eta_prior,
        beta0_prior = prior$beta0_prior,
        completed_at = E$timestamp_now()
      ),
      s
    )
    
    atomic_write(result, path)
    
    if (!KEEP_FULL_DRAWS) {
      unlink(raw_path)
    }
    
    cat(sprintf(
      "Saved %s, replication %03d\n",
      id,
      r
    ))
    
    invisible(status)
  }
  
  ## ------------------------------------------------------------------------
  ## 6. ETA TABLES
  ## ------------------------------------------------------------------------
  
  eta_mean_row <- function(df, id) {
    data.frame(
      Setting = id,
      Replications = nrow(df),
      True_eta = true_eta,
      
      Mean_posterior_mean = mean(
        df$eta_posterior_mean
      ),
      
      Mean_HPD80_lower = mean(
        df$eta_HPD80_lower
      ),
      
      Mean_HPD80_upper = mean(
        df$eta_HPD80_upper
      ),
      
      Mean_HPD80_width = mean(
        df$eta_HPD80_width
      ),
      
      Coverage_HPD80 = mean(
        df$HPD80_covers_true
      )
    )
  }
  
  eta_display <- function(df) {
    ans <- df[
      ,
      c(
        "Setting",
        "Replications",
        "Mean_posterior_mean",
        "Mean_HPD80_lower",
        "Mean_HPD80_upper"
      ),
      drop = FALSE
    ]
    
    ## Formatting only. Calculations and other output files remain unrounded.
    for (name in names(ans)[-(1:2)]) {
      ans[[name]] <- sprintf("%.2f", ans[[name]])
    }
    
    ans
  }
  
  setting_summary <- function(id, required_reps) {
    results <- lapply(
      required_reps,
      function(r) {
        read_checked(checkpoint_path(id, r), id, r)
      }
    )
    
    eta <- bind_rows(
      lapply(results, `[[`, "eta_summary")
    )
    
    folder <- file.path(OUT_DIR, "tables", id)
    
    write_table(
      eta,
      file.path(folder, "eta_by_replication.csv")
    )
    
    exact <- eta_mean_row(eta, id)
    
    write_table(
      exact,
      file.path(folder, "eta_means_full_precision.csv")
    )
    
    write_table(
      eta_display(exact),
      file.path(folder, "eta_means_2dp.csv")
    )
    
    results
  }
  
  ## ------------------------------------------------------------------------
  ## 7. WORKER MODE
  ## ------------------------------------------------------------------------
  
  if (MODE == "run") {
    if (
      length(SESSION_ID) != 1L ||
      !SESSION_ID %in% 1:5
    ) {
      stop("SESSION_ID must be 1, 2, 3, 4, or 5.")
    }
    
    reps <- as.integer(REPLICATIONS)
    
    if (
      !length(reps) ||
      anyNA(reps) ||
      any(REPLICATIONS != reps) ||
      anyDuplicated(reps) ||
      any(!reps %in% 1:100)
    ) {
      stop(
        "REPLICATIONS must contain unique integers from 1 to 100."
      )
    }
    
    selected <- batches[[SESSION_ID]]
    
    counts <- c(
      computed = 0L,
      loaded = 0L,
      raw_loaded = 0L
    )
    
    for (r in reps) {
      prepare_replication(r)
      
      for (id in selected) {
        status <- run_one(id, r)
        counts[status] <- counts[status] + 1L
        invisible(gc())
      }
    }
    
    ## Each session writes only tables for its own configurations.
    for (id in selected) {
      present <- which(
        vapply(
          1:100,
          function(r) file.exists(checkpoint_path(id, r)),
          logical(1)
        )
      )
      
      results <- setting_summary(id, present)
      
      eta_now <- bind_rows(
        lapply(results, `[[`, "eta_summary")
      )
      
      cat("\n", id, "\n", sep = "")
      
      print(
        eta_display(eta_mean_row(eta_now, id)),
        row.names = FALSE
      )
    }
    
    cat("\nSession complete.\n")
    print(counts)
    
    cat(
      "\nUse MODE = 'combine' after all five sessions finish.\n"
    )
    
    ## Return normally. Never quit the R session.
    return(invisible(NULL))
  }
  
  ## ------------------------------------------------------------------------
  ## 8. COMBINE MODE -- NO MCMC
  ## ------------------------------------------------------------------------
  
  combined_file <- file.path(
    OUT_DIR,
    "combined",
    "MC100_combined_summaries.rds"
  )
  
  if (MODE == "combine") {
    tasks <- expand.grid(
      replicate = 1:100,
      fit_id = fit_ids,
      stringsAsFactors = FALSE
    )
    
    paths <- mapply(
      checkpoint_path,
      tasks$fit_id,
      tasks$replicate,
      USE.NAMES = FALSE
    )
    
    missing <- which(!file.exists(paths))
    
    if (length(missing)) {
      print(
        head(tasks[missing, , drop = FALSE], 20L),
        row.names = FALSE
      )
      
      stop(
        length(missing),
        " of the 1,300 fit summaries are missing. ",
        "No MCMC was started."
      )
    }
    
    defaults <- lapply(
      1:100,
      function(r) {
        read_checked(
          checkpoint_path("DEFAULT", r),
          "DEFAULT",
          r
        )
      }
    )
    
    ## One DEFAULT-based ordering per replication.
    ## Every configuration uses that same ordering within that replication.
    orders <- lapply(
      defaults,
      function(z) {
        d <- z$node_summary
        
        d$node[
          order(
            match(d$node_group, group_levels),
            d$latent_probability_mean,
            d$node
          )
        ]
      }
    )
    
    common_orders <- bind_rows(
      lapply(
        1:100,
        function(r) {
          d <- defaults[[r]]$node_summary
          
          data.frame(
            replicate = r,
            plot_position = 1:100,
            original_node = orders[[r]],
            node_group = d$node_group[
              match(orders[[r]], d$node)
            ]
          )
        }
      )
    )
    
    plot_tables <- list()
    eta_rows <- list()
    robustness_rows <- list()
    topk_rows <- list()
    
    for (id in fit_ids) {
      results <- setting_summary(id, 1:100)
      
      ## Verify common datasets across all configurations.
      for (r in 1:100) {
        if (
          !identical(
            results[[r]]$toy_hash,
            defaults[[r]]$toy_hash
          )
        ) {
          stop(
            "Configurations did not use identical data in replication ",
            r
          )
        }
      }
      
      folder <- file.path(OUT_DIR, "tables", id)
      
      ## Preserve the original summaries and diagnostics by replication.
      table_names <- c(
        "node_summary",
        "eta_summary",
        "diagnostics",
        "prior_calibration",
        "parameter_summary",
        "core_parameter_ess",
        "node_probability_ess",
        "node_ess_summary"
      )
      
      for (name in table_names) {
        write_table(
          bind_rows(lapply(results, `[[`, name)),
          file.path(
            folder,
            paste0(name, "_by_replication.csv")
          )
        )
      }
      
      make_average <- function(use_default_order) {
        aligned <- lapply(
          1:100,
          function(r) {
            d <- results[[r]]$node_summary
            
            idx <- if (use_default_order) {
              orders[[r]]
            } else {
              1:100
            }
            
            d[
              match(idx, d$node),
              c(
                "latent_probability_mean",
                "latent_HPD80_lower",
                "latent_HPD80_upper"
              ),
              drop = FALSE
            ]
          }
        )
        
        means <- Reduce(
          `+`,
          lapply(aligned, as.matrix)
        ) / n_rep
        
        data.frame(
          fit_id = id,
          position = 1:100,
          
          node_group = c(
            rep("true zero", 70),
            rep("hidden positive", 10),
            rep("observed positive", 20)
          ),
          
          n_replications = n_rep,
          Mean_posterior_mean = means[, 1L],
          Mean_HPD80_lower = means[, 2L],
          Mean_HPD80_upper = means[, 3L]
        )
      }
      
      ## Primary plotting table: average at matched DEFAULT-rank positions.
      plot_tables[[id]] <- make_average(TRUE)
      
      write_table(
        plot_tables[[id]],
        file.path(
          folder,
          "predictive_means_by_default_rank.csv"
        )
      )
      
      ## Also retain the literal averages for each original node index.
      write_table(
        make_average(FALSE),
        file.path(
          folder,
          "predictive_means_by_original_node.csv"
        )
      )
      
      eta_rows[[id]] <- eta_mean_row(
        bind_rows(lapply(results, `[[`, "eta_summary")),
        id
      )
      
      ## Original robustness comparisons, calculated within each replication.
      rob <- list()
      top <- list()
      
      for (r in 1:100) {
        d <- results[[r]]$node_summary
        b <- defaults[[r]]$node_summary
        
        p <- d$latent_probability_mean[
          match(1:100, d$node)
        ]
        
        p0 <- b$latent_probability_mean[
          match(1:100, b$node)
        ]
        
        dg <- results[[r]]$diagnostics
        
        rob[[r]] <- data.frame(
          fit_id = id,
          replicate = r,
          
          spearman_with_default = suppressWarnings(
            cor(p, p0, method = "spearman")
          ),
          
          mean_absolute_probability_difference = mean(
            abs(p - p0)
          ),
          
          max_absolute_probability_difference = max(
            abs(p - p0)
          ),
          
          mean_HPD80_width = mean(
            d$latent_HPD80_width
          ),
          
          overall_AUC = dg$overall_AUC_latent_mean,
          hidden_AUC = dg$hidden_AUC_latent_mean,
          runtime_sec = dg$runtime_sec
        )
        
        top[[r]] <- bind_rows(
          lapply(
            E$TOP_K_VALUES,
            function(k) {
              overlap <- length(
                intersect(
                  E$ordered_top_nodes(p, k),
                  E$ordered_top_nodes(p0, k)
                )
              )
              
              data.frame(
                fit_id = id,
                replicate = r,
                k = k,
                overlap_count = overlap,
                overlap_fraction = overlap / k
              )
            }
          )
        )
      }
      
      robustness_rows[[id]] <- bind_rows(rob)
      topk_rows[[id]] <- bind_rows(top)
      
      write_table(
        robustness_rows[[id]],
        file.path(folder, "robustness_by_replication.csv")
      )
      
      write_table(
        topk_rows[[id]],
        file.path(folder, "topk_overlap_by_replication.csv")
      )
      
      metric_names <- setdiff(
        names(robustness_rows[[id]]),
        c("fit_id", "replicate")
      )
      
      write_table(
        data.frame(
          Metric = metric_names,
          
          Mean = vapply(
            robustness_rows[[id]][metric_names],
            mean,
            numeric(1)
          ),
          
          SD = vapply(
            robustness_rows[[id]][metric_names],
            sd,
            numeric(1)
          )
        ),
        file.path(folder, "robustness_MC_mean_sd.csv")
      )
    }
    
    eta_all <- bind_rows(eta_rows)
    
    combined <- list(
      experiment = experiment,
      plot_means = bind_rows(plot_tables),
      eta_means = eta_all,
      eta_table_2dp = eta_display(eta_all),
      common_orders = common_orders,
      robustness = bind_rows(robustness_rows),
      topk_overlap = bind_rows(topk_rows),
      
      interval_note = paste(
        "HPD endpoints were calculated separately in each fit and then averaged.",
        "They are NOT an HPD interval from pooled posterior draws."
      )
    )
    
    atomic_write(combined, combined_file)
    
    write_table(
      seed_plan,
      file.path(OUT_DIR, "combined", "seed_plan.csv")
    )
    
    write_table(
      common_orders,
      file.path(
        OUT_DIR,
        "combined",
        "default_order_by_replication.csv"
      )
    )
    
    write_table(
      combined$plot_means,
      file.path(
        OUT_DIR,
        "combined",
        "ALL_predictive_mean_HPD80.csv"
      )
    )
    
    write_table(
      eta_all,
      file.path(
        OUT_DIR,
        "combined",
        "ALL_eta_means_full_precision.csv"
      )
    )
    
    write_table(
      combined$eta_table_2dp,
      file.path(
        OUT_DIR,
        "combined",
        "ALL_eta_means_2dp.csv"
      )
    )
    
    print(
      combined$eta_table_2dp,
      row.names = FALSE
    )
    
    cat(
      "\nCombined 100 replications for all 13 settings.",
      "\nNo MCMC was run.\n"
    )
  }
  
  ## ------------------------------------------------------------------------
  ## 9. PLOT SAVED AVERAGES -- ORIGINAL PANEL FUNCTIONS AND STYLE
  ## ------------------------------------------------------------------------
  
  if (!file.exists(combined_file)) {
    stop("Run MODE = 'combine' first.")
  }
  
  combined <- readRDS(combined_file)
  
  if (!identical(combined$experiment, experiment)) {
    stop("Combined bundle settings do not match.")
  }
  
  d <- combined$plot_means
  
  ## Here 'node' is the common plotting position, not an individual
  ## with the same identity across independent datasets.
  E$node_summaries_long_all <- data.frame(
    fit_id = d$fit_id,
    node = d$position,
    latent_probability_mean = d$Mean_posterior_mean,
    latent_HPD80_lower = d$Mean_HPD80_lower,
    latent_HPD80_upper = d$Mean_HPD80_upper
  )
  
  E$common_plot_order <- 1:100
  E$N_NODES <- 100L
  
  E$baseline_eta_beta_display <-
    E$make_eta_prior_from_spec(specs$DEFAULT)$b_eta
  
  ## Reuse the colours and panel specifications of the original script.
  for (
    name in c(
      "conservative_base",
      "diffused_base",
      "column_colours",
      "column_order",
      "column_headers",
      "panel_specs",
      "row_order"
    )
  ) {
    eval(original_expression(name), E)
  }
  
  panels <- E$panel_specs
  
  panels$joint <- list(
    row_label = "Joint priors",
    
    fit_ids = c(
      "JOINT_CONSERVATIVE",
      "DEFAULT",
      "JOINT_DIFFUSE"
    ),
    
    titles = list(
      "All conservative",
      "Default",
      "All diffuse"
    )
  )
  
  rows <- c(E$row_order, "joint")
  
  plots <- list(
    E$make_text_strip("", 1),
    E$make_text_strip(E$column_headers["conservative"], 12),
    E$make_text_strip(E$column_headers["default"], 12),
    E$make_text_strip(E$column_headers["diffuse"], 12)
  )
  
  for (j in seq_along(rows)) {
    ps <- panels[[rows[j]]]
    
    plots[[length(plots) + 1L]] <- if (rows[j] == "joint") {
      E$make_vertical_text_strip(
        "Joint priors",
        fontsize = 11
      )
    } else {
      E$make_text_strip(ps$row_label, 15)
    }
    
    for (k in seq_along(E$column_order)) {
      plots[[length(plots) + 1L]] <-
        E$make_sensitivity_panel_6x3(
          fit_id_now = ps$fit_ids[k],
          column_key_now = E$column_order[k],
          panel_title_now = ps$titles[[k]],
          show_x_title = j == length(rows),
          show_y_title = k == 1L
        )
    }
  }
  
  figure <- patchwork::wrap_plots(
    plots,
    ncol = 4,
    widths = c(0.08, 1, 1, 1),
    heights = c(0.08, rep(1, 6)),
    byrow = TRUE
  ) +
    patchwork::plot_annotation(
      title = "Prior sensitivity analysis for the PU-PoE illustration",
      
      subtitle = paste(
        "100 replications:",
        "mean posterior probabilities and averaged 80% HPD endpoints"
      ),
      
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 14,
          hjust = 0.5,
          margin = ggplot2::margin(b = 4)
        ),
        
        plot.subtitle = ggplot2::element_text(
          size = 8,
          hjust = 0.5
        )
      )
    )
  
  dir.create(
    file.path(OUT_DIR, "figures"),
    showWarnings = FALSE
  )
  
  pdf_path <- file.path(
    OUT_DIR,
    "figures",
    "PU_PoE_PRIOR_SENSITIVITY_MC100_6x3.pdf"
  )
  
  ggplot2::ggsave(
    filename = pdf_path,
    plot = figure,
    device = "pdf",
    useDingbats = FALSE,
    width = 8.27,
    height = 11.69,
    units = "in",
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    filename = sub("\\.pdf$", ".png", pdf_path),
    plot = figure,
    width = 8.27,
    height = 11.69,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  cat(
    "\nSaved figure:\n",
    pdf_path,
    "\n",
    sep = ""
  )
  
  invisible(combined)
}
