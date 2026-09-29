## Run from the project root: Rscript --vanilla tests/test_static.R
## Checks parsing, the functions of the original script and the settings.
stopifnot(file.exists(".here"))

files <- c(
  list.files("R", "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files("config", "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files("scripts", "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files("tests", "[.]R$", recursive = TRUE, full.names = TRUE)
)

for (f in files) {
  parse(file = f, keep.source = FALSE)
}

wf <- new.env(parent = as.environment("package:stats"))
sys.source("R/illustration/workflow.R", envir = wf)
quick <- wf$new_illustration_environment(getwd(), "quick")
full <- wf$new_illustration_environment(getwd(), "full")

stopifnot(
  identical(full$BLOCKED_BURNIN, 10000L),
  identical(full$JOINT_BURNIN, 20000L),
  identical(full$POSTERIOR_ITERATIONS, 150000L),
  identical(full$THIN, 30L),
  identical(full$N_RETAINED, 5000L),
  identical(quick$N_RETAINED, 200L)
)

stopifnot(
  unname(tools::md5sum("reference/GP_four_model_illustrtation.R")) ==
    "9ed9a4a9896420fd883b089b17061966"
)

reference <- parse(
  "reference/GP_four_model_illustrtation.R",
  keep.source = FALSE
)

reference_functions <- new.env(parent = as.environment("package:stats"))

is_function_assignment <- function(x) {
  is.call(x) &&
    identical(x[[1L]], as.name("<-")) &&
    is.call(x[[3L]]) &&
    identical(x[[3L]][[1L]], as.name("function"))
}

for (x in reference) {
  if (is_function_assignment(x)) eval(x, envir = reference_functions)
}

map <- read.csv(
  "docs/SOURCE_MAP.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

verbatim <- map[["function"]][map$status == "verbatim"]

for (name in verbatim) {
  old <- get(name, envir = reference_functions, inherits = FALSE)
  new <- get(name, envir = full, inherits = FALSE)
  stopifnot(
    identical(formals(old), formals(new)),
    identical(body(old), body(new))
  )
}

## The sampler differs from the original only in the schedule check.
old_source <- paste(
  readLines("reference/GP_four_model_illustrtation.R", warn = FALSE),
  collapse = "\n"
)

guard <- paste(
  readLines("tests/reference_schedule_guard.txt", warn = FALSE),
  collapse = "\n"
)

replacement <- paste(
  readLines("tests/reference_schedule_replacement.txt", warn = FALSE),
  collapse = "\n"
)

stopifnot(
  length(gregexpr(guard, old_source, fixed = TRUE)[[1L]]) == 1L,
  gregexpr(guard, old_source, fixed = TRUE)[[1L]][1L] > 0L
)

adjusted <- sub(guard, replacement, old_source, fixed = TRUE)
adjusted_functions <- new.env(parent = as.environment("package:stats"))

for (x in parse(text = adjusted, keep.source = FALSE)) {
  if (is_function_assignment(x)) eval(x, envir = adjusted_functions)
}

stopifnot(
  identical(
    formals(adjusted_functions$run_exact_toy_model_chain),
    formals(full$run_exact_toy_model_chain)
  ),
  identical(
    body(adjusted_functions$run_exact_toy_model_chain),
    body(full$run_exact_toy_model_chain)
  )
)

## This function differs from the original only in one error message.
old_summary <- reference_functions$summarize_completed_model

expected_summary <- gsub(
  " did not return the expected 5,000 x 100 probability matrices.",
  " did not return the configured probability-matrix dimensions.",
  paste(deparse(body(old_summary)), collapse = "\n"),
  fixed = TRUE
)

stopifnot(identical(
  parse(text = expected_summary, keep.source = FALSE)[[1L]],
  body(full$summarize_completed_model)
))

compute_keys <- c(
  "BLOCKED_BURNIN",
  "JOINT_BURNIN",
  "POSTERIOR_ITERATIONS",
  "THIN",
  "N_RETAINED",
  "ADAPT_START",
  "ADAPT_INTERVAL",
  "COR_ADAPT_START",
  "COR_ADAPT_INTERVAL",
  "COR_ADAPT_WINDOW",
  "EMP_COV_TAIL",
  "BLOCKED_PROGRESS_EVERY",
  "JOINT_PROGRESS_EVERY",
  "SAMPLE_PROGRESS_EVERY"
)

for (key in setdiff(full$CONFIGURATION_KEYS, compute_keys)) {
  stopifnot(identical(get(key, envir = full), get(key, envir = quick)))
}

cat("PASS: parsing, the functions of the original script and the settings.\n")
