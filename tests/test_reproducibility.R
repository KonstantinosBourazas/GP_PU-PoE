## Run from project root: Rscript --vanilla tests/test_reproducibility.R
## Runs all four quick models twice, and the original functions once, in three
## new R sessions. Also checks plotting and the reuse of completed fits.
stopifnot(file.exists(".here"))
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
## The output lock is tested first.
lock_test_env <- new.env(parent = as.environment("package:stats"))
sys.source(file.path(root, "tests", "test_locking.R"), envir = lock_test_env)

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

if (!file.exists(rscript)) {
  stop("Rscript executable was not found in the current R installation.")
}

base <- tempfile("reproducibility_", tmpdir = file.path(root, "results"))
dir.create(base, recursive = TRUE)
base <- normalizePath(base, winslash = "/", mustWork = TRUE)
a <- file.path(base, "new_A")
b <- file.path(base, "new_B")
ref <- file.path(base, "reference")

run_child <- function(script, arguments, logfile) {
  status <- system2(
    rscript,
    args = c("--vanilla", shQuote(script), shQuote(arguments)),
    stdout = logfile,
    stderr = logfile
  )
  if (!identical(as.integer(status), 0L)) {
    detail <- tryCatch(
      tail(readLines(logfile, warn = FALSE), 25L),
      error = function(err) "Could not read the log."
    )
    stop(
      "A run failed. See: ",
      logfile,
      "\n",
      paste(detail, collapse = "\n"),
      call. = FALSE
    )
  }
}

assert_unlocked <- function(output) {
  lock <- file.path(output, ".run_lock")
  if (dir.exists(lock) || file.exists(lock)) {
    stop(
      "A completed run left an output lock: ",
      lock,
      ". See its log.",
      call. = FALSE
    )
  }
}

entry <- file.path(root, "scripts", "01_illustrative_example.R")

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

run_child(
  file.path(root, "tests", "reference_worker.R"),
  ref,
  file.path(base, "reference.log")
)

ids <- c("GP_COV", "GP_NET", "GP_POE", "GP_PU_POE")
comparisons <- list()

read_fit <- function(dir, id) {
  readRDS(file.path(dir, "model_checkpoints", paste0(id, "_complete.rds")))
}

for (id in ids) {
  x <- read_fit(a, id)
  y <- read_fit(b, id)
  z <- read_fit(ref, id)
  for (pair_name in c("new_A_vs_new_B", "new_A_vs_reference")) {
    other <- if (pair_name == "new_A_vs_new_B") y else z
    fields <- list(
      draws = identical(x$fit$draws, other$fit$draws),
      state = identical(x$fit$state, other$fit$state),
      proposal = identical(x$fit$proposal, other$fit$proposal),
      acceptance = identical(x$fit$acceptance, other$fit$acceptance),
      phases = identical(x$fit$phases, other$fit$phases),
      prior_bundle = identical(x$prior_bundle, other$prior_bundle),
      node_summary = identical(x$node_summary, other$node_summary),
      eta_summary = identical(x$eta_summary, other$eta_summary)
    )
    comparisons[[length(comparisons) + 1L]] <- data.frame(
      model = id,
      comparison = pair_name,
      check = names(fields),
      identical = unlist(fields, use.names = FALSE)
    )
  }
}

table <- do.call(rbind, comparisons)
write.csv(table, file.path(base, "comparison_results.csv"), row.names = FALSE)

if (!all(table$identical)) {
  stop("The runs differ. See comparison_results.csv in ", base)
}

checkpoints <- list.files(
  file.path(a, "model_checkpoints"),
  "_complete[.]rds$",
  full.names = TRUE
)

before <- tools::md5sum(checkpoints)

run_child(
  entry,
  c("--mode=quick", "--action=plot", paste0("--out=", a)),
  file.path(base, "plot_only.log")
)

assert_unlocked(a)
stopifnot(identical(before, tools::md5sum(checkpoints)))

run_child(
  entry,
  c("--mode=quick", "--action=all", paste0("--out=", a)),
  file.path(base, "reuse.log")
)

assert_unlocked(a)
stopifnot(identical(before, tools::md5sum(checkpoints)))

writeLines(
  capture.output(sessionInfo()),
  file.path(base, "test_sessionInfo.txt")
)

cat(
  "PASS: two quick runs agree exactly and agree with the original functions; plotting and reuse leave the checkpoints unchanged.\n"
)

cat("Reports: ", base, "\n", sep = "")
