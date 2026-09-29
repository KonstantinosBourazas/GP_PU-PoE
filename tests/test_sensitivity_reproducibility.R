## Runs the 13 quick fits twice and the original functions for five
## specifications, and checks plotting, reuse, cleanup after errors and locks.
stopifnot(file.exists(".here"))
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
local_env <- new.env(parent = as.environment("package:stats"))

sys.source(
  file.path(root, "tests", "test_sensitivity_static.R"),
  envir = local_env
)

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

if (!file.exists(rscript)) {
  stop("Rscript not found in the current R installation.")
}

base <- tempfile(
  "sensitivity_reproducibility_",
  tmpdir = file.path(root, "results")
)

dir.create(base, recursive = TRUE)
base <- normalizePath(base, winslash = "/", mustWork = TRUE)
a <- file.path(base, "new_A")
b <- file.path(base, "new_B")
ref <- file.path(base, "reference")
entry <- file.path(root, "scripts", "02_prior_sensitivity.R")

run_child <- function(script, arguments, log, expect_error = FALSE) {
  status <- suppressWarnings(system2(
    rscript,
    args = c("--vanilla", shQuote(script), shQuote(arguments)),
    stdout = log,
    stderr = log
  ))
  ok <- identical(as.integer(status), 0L)
  if (ok == isTRUE(expect_error)) {
    tail_log <- tryCatch(
      tail(readLines(log, warn = FALSE), 30L),
      error = function(e) conditionMessage(e)
    )
    stop(
      "Unexpected exit status. See: ",
      log,
      "\n",
      paste(tail_log, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(status)
}

assert_unlocked <- function(dir) {
  stopifnot(
    !dir.exists(file.path(dir, ".run_lock")),
    !file.exists(file.path(dir, ".run_lock"))
  )
}

snapshot_files <- function(dir) {
  files <- sort(list.files(
    dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  ))
  files <- files[!dir.exists(files)]
  if (!length(files)) {
    return(character(0))
  }
  tools::md5sum(files)
}

old_illustration <- snapshot_files(file.path(root, "results", "illustration"))
# Cleanup after a normal run and after a failed plot; an existing lock is kept.
check <- file.path(base, "check")

run_child(
  entry,
  c("--mode=quick", "--action=check", paste0("--out=", check)),
  file.path(base, "check.log")
)

assert_unlocked(check)
failed <- file.path(base, "missing_plot")

run_child(
  entry,
  c("--mode=quick", "--action=plot", paste0("--out=", failed)),
  file.path(base, "expected_plot_failure.log"),
  TRUE
)

assert_unlocked(failed)
locked <- file.path(base, "preexisting_lock")
dir.create(file.path(locked, ".run_lock"), recursive = TRUE)
sentinel <- file.path(locked, ".run_lock", "sentinel.txt")
writeLines("existing lock", sentinel)

run_child(
  entry,
  c("--mode=quick", "--action=check", paste0("--out=", locked)),
  file.path(base, "expected_lock_refusal.log"),
  TRUE
)

stopifnot(file.exists(sentinel))
# This lock was created by the test and stays in the report folder.
cat("PASS: output lock and cleanup after errors.\n")

run_child(
  entry,
  c("--mode=quick", "--action=all", paste0("--out=", a)),
  file.path(base, "new_A.log")
)

assert_unlocked(a)

run_child(
  entry,
  c("--mode=quick", "--action=fit", paste0("--out=", b)),
  file.path(base, "new_B.log")
)

assert_unlocked(b)

ids <- c(
  "DEFAULT",
  "SIGMA_CONSERVATIVE",
  "SIGMA_DIFFUSE",
  "ELL_CONSERVATIVE",
  "ELL_DIFFUSE",
  "GAMMA_CONSERVATIVE",
  "GAMMA_DIFFUSE",
  "BETA0_CONSERVATIVE",
  "BETA0_DIFFUSE",
  "ETA_CONSERVATIVE",
  "ETA_DIFFUSE",
  "JOINT_CONSERVATIVE",
  "JOINT_DIFFUSE"
)

reference_ids <- c(
  "DEFAULT",
  "ELL_DIFFUSE",
  "ETA_CONSERVATIVE",
  "JOINT_CONSERVATIVE",
  "JOINT_DIFFUSE"
)

run_child(
  file.path(root, "sensitivity", "tests", "reference_worker.R"),
  c(ref, paste(reference_ids, collapse = ",")),
  file.path(base, "reference.log")
)

read_fit <- function(dir, id) {
  readRDS(file.path(dir, "fit_checkpoints", paste0(id, "_complete.rds")))
}

rows <- list()

for (id in ids) {
  x <- read_fit(a, id)
  others <- list(new_A_vs_new_B = read_fit(b, id))
  if (id %in% reference_ids) {
    others$new_A_vs_original_functions <- read_fit(ref, id)
  }
  for (comparison in names(others)) {
    y <- others[[comparison]]
    tests <- c(
      setNames(
        vapply(
          c("draws", "state", "proposal", "acceptance", "phases"),
          function(k) identical(x$fit[[k]], y$fit[[k]]),
          logical(1)
        ),
        c("draws", "state", "proposal", "acceptance", "phases")
      ),
      prior_bundle = identical(x$prior_bundle, y$prior_bundle),
      node_summary = identical(x$node_summary, y$node_summary),
      eta_summary = identical(x$eta_summary, y$eta_summary)
    )
    rows[[length(rows) + 1L]] <- data.frame(
      fit = id,
      comparison = comparison,
      check = names(tests),
      identical = unname(tests)
    )
  }
}

table <- do.call(rbind, rows)
write.csv(table, file.path(base, "comparison_results.csv"), row.names = FALSE)

if (!all(table$identical)) {
  stop("The runs differ. See ", file.path(base, "comparison_results.csv"))
}

cp <- sort(list.files(
  file.path(a, "fit_checkpoints"),
  "_complete[.]rds$",
  full.names = TRUE
))

before <- tools::md5sum(cp)

run_child(
  entry,
  c("--mode=quick", "--action=plot", paste0("--out=", a)),
  file.path(base, "plot_only.log")
)

assert_unlocked(a)
stopifnot(identical(before, tools::md5sum(cp)))

run_child(
  entry,
  c("--mode=quick", "--action=all", paste0("--out=", a)),
  file.path(base, "reuse.log")
)

assert_unlocked(a)

stopifnot(
  identical(before, tools::md5sum(cp)),
  identical(
    old_illustration,
    snapshot_files(file.path(root, "results", "illustration"))
  )
)

map <- read.csv(
  file.path(a, "toy_seed53_PU_PoE_prior_sensitivity_6x3_plot_grid_map.csv"),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(map) == 18L,
  sum(map$fit_id == "DEFAULT") == 6L,
  length(unique(map$fit_id)) == 13L
)

writeLines(
  capture.output(sessionInfo()),
  file.path(base, "sessionInfo_test.txt")
)

cat(
  "PASS: two quick runs of the 13 fits agree exactly and agree with the original functions for five specifications; 18 panels from 13 fits; plotting and reuse leave the checkpoints unchanged.\n"
)

cat("Reports: ", base, "\n", sep = "")
