## Tests the output of script 02b. Run after script 02b with ACTION "all".
local({
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  code <- file.path(root, "sensitivity_mc100")
  out <- file.path(root, "results", "prior_sensitivity", "mc100")
  archive <- file.path(root, "archive", "prior_sensitivity_mc100")
  if (!file.exists(file.path(code, "R", "workflow.R"))) {
    stop("Open GP-PU-PoE.Rproj before running this test.")
  }
  for (p in c(
    list.files(file.path(code, "R"), "[.]R$", full.names = TRUE),
    file.path(root, "scripts", "02b_prior_sensitivity_MC100.R")
  )) {
    parse(p)
  }
  e <- new.env(parent = as.environment("package:stats"))
  for (name in c(
    "combine.R",
    "metrics.R",
    "original_plot_functions.R",
    "plot.R"
  )) {
    sys.source(file.path(code, "R", name), envir = e, keep.source = FALSE)
  }
  before <- e$mc100_manifest_check(root, archive)
  p <- file.path(out, "combined", "MC100_combined_summaries.rds")
  if (!file.exists(p)) {
    stop("Run script 02b with ACTION 'all' first.")
  }
  bundle <- readRDS(p)
  e$mc100_plot_input_check(bundle)
  for (item in list(
    c("plot_means", "expected_prediction_means.csv"),
    c("eta_means", "expected_eta_means.csv")
  )) {
    expected <- utils::read.csv(
      file.path(code, "reference", item[2L]),
      stringsAsFactors = FALSE
    )
    e$mc100_compare_table(bundle[[item[1L]]], expected, item[1L])
  }
  for (item in list(
    c("summary", "expected_auc_summary.csv"),
    c("pairwise_correlations", "expected_auc_pairwise.csv"),
    c("range_summary", "expected_auc_range_summary.csv"),
    c("ranking_stability", "expected_ranking_stability.csv")
  )) {
    expected <- utils::read.csv(
      file.path(code, "reference", item[2L]),
      stringsAsFactors = FALSE
    )
    e$mc100_compare_table(
      bundle$auc_reports[[item[1L]]],
      expected,
      item[1L],
      tolerance = 1e-9
    )
  }
  for (r in 1:100) {
    d <- bundle$node_positions[bundle$node_positions$replicate == r, ]
    stopifnot(
      identical(as.integer(d$original_node), 1:100),
      identical(as.integer(d$plot_position), 1:100)
    )
  }
  d <- e$mc100_eta_plot_data(bundle)
  stopifnot(
    nrow(d) == 18L,
    length(unique(d$fit_id)) == 13L,
    length(unique(d$Mean_posterior_mean[d$fit_id == "DEFAULT"])) == 1L
  )
  # The four plotting functions are those of the original script.
  original <- parse(
    file.path(code, "reference", "prior_sensitivity_single_run_original.R"),
    keep.source = FALSE
  )
  for (key in c(
    "blend_colour",
    "make_sensitivity_panel_6x3",
    "make_text_strip",
    "make_vertical_text_strip"
  )) {
    hits <- which(vapply(
      original,
      function(x) {
        is.call(x) &&
          length(x) == 3L &&
          identical(x[[1L]], as.name("<-")) &&
          identical(x[[2L]], as.name(key))
      },
      logical(1)
    ))
    stopifnot(length(hits) == 1L)
    ref_env <- new.env(parent = as.environment("package:stats"))
    eval(original[[hits]], ref_env)
    stopifnot(
      identical(body(get(key, e)), body(get(key, ref_env))),
      identical(formals(get(key, e)), formals(get(key, ref_env)))
    )
  }
  eta_plot <- e$mc100_eta_plot(bundle)
  pred_plot <- e$mc100_prediction_plot(bundle)
  no_text <- function(x) is.null(x) || inherits(x, "waiver")
  stopifnot(
    no_text(eta_plot$labels$subtitle),
    no_text(eta_plot$labels$caption),
    no_text(pred_plot$patches$annotation$subtitle),
    no_text(pred_plot$patches$annotation$caption)
  )
  for (stem in c(
    "PU_PoE_PRIOR_SENSITIVITY_MC100_predictions_6x3",
    "PU_PoE_PRIOR_SENSITIVITY_MC100_eta_wide"
  )) {
    for (ext in c("pdf", "png")) {
      pp <- file.path(out, "figures", paste0(stem, ".", ext))
      stopifnot(file.exists(pp), file.info(pp)$size > 1000)
    }
  }
  after <- e$mc100_manifest_check(root, archive)
  stopifnot(identical(before, after))
  cat(
    "PASS: means at the original unit indices, eta values, AUC and rank tables, figures and archive.\n"
  )
})
