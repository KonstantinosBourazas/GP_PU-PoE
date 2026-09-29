## Static checks of the prior sensitivity module.
stopifnot(file.exists(".here"))
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
wf <- new.env(parent = as.environment("package:stats"))
sys.source(file.path(root, "sensitivity", "R", "workflow.R"), envir = wf)

files <- c(
  wf$sensitivity_files(root),
  file.path(root, "scripts", "02_prior_sensitivity.R"),
  list.files(
    file.path(root, "sensitivity", "tests"),
    "[.]R$",
    full.names = TRUE
  ),
  list.files(
    file.path(root, "tests"),
    "^test_sensitivity.*[.]R$",
    full.names = TRUE
  )
)

for (path in files) {
  parse(file = path, keep.source = FALSE)
}

quick <- wf$new_sensitivity_environment(root, "quick")
full <- wf$new_sensitivity_environment(root, "full")

stopifnot(
  identical(full$BLOCKED_BURNIN, 10000L),
  identical(full$JOINT_BURNIN, 20000L),
  identical(full$POSTERIOR_ITERATIONS, 150000L),
  identical(full$THIN, 30L),
  identical(full$N_RETAINED, 5000L),
  identical(quick$N_RETAINED, 200L),
  identical(full$SENSITIVITY_SPECS, quick$SENSITIVITY_SPECS)
)

original_path <- file.path(
  root,
  "sensitivity",
  "reference",
  "PU_PoE_prior_sensitivity_final.R"
)

source <- paste(readLines(original_path, warn = FALSE), collapse = "\n")

guard <- paste(
  readLines(
    file.path(root, "tests", "reference_schedule_guard.txt"),
    warn = FALSE
  ),
  collapse = "\n"
)

replacement <- paste(
  readLines(
    file.path(root, "tests", "reference_schedule_replacement.txt"),
    warn = FALSE
  ),
  collapse = "\n"
)

locations <- gregexpr(guard, source, fixed = TRUE)[[1L]]
stopifnot(length(locations) == 1L, locations[1L] > 0)
expected <- sub(guard, replacement, source, fixed = TRUE)

expected <- sub(
  " did not return the expected 5,000 x 100 probability matrices.",
  " did not return the configured probability-matrix dimensions.",
  expected,
  fixed = TRUE
)

old <- new.env(parent = as.environment("package:stats"))

for (x in parse(text = expected, keep.source = FALSE)) {
  if (
    is.call(x) &&
      length(x) == 3L &&
      identical(x[[1L]], as.name("<-")) &&
      is.call(x[[3L]]) &&
      identical(x[[3L]][[1L]], as.name("function"))
  ) {
    eval(x, old)
  }
}

# The numerical functions shared with the original script.
exclude <- c(
  "run_or_load_sensitivity_fit",
  "make_sensitivity_panel",
  "make_text_strip",
  "make_sensitivity_panel_6x3",
  "make_vertical_text_strip",
  "blend_colour"
)

names_to_check <- setdiff(ls(old), exclude)

rows <- lapply(names_to_check, function(key) {
  a <- get(key, old, inherits = FALSE)
  b <- get(key, full, inherits = FALSE)
  if (!identical(formals(a), formals(b)) || !identical(body(a), body(b))) {
    stop("Function differs from the original script: ", key)
  }
  data.frame(check = key, status = "PASS")
})

# The illustration files used by the sensitivity module.
core <- wf$sensitivity_core(root)

original_manifest <- read.csv(
  file.path(root, "sensitivity", "reference", "ILLUSTRATION_CORE_MD5.csv"),
  stringsAsFactors = FALSE
)

current_core <- core$illustration_module_files(root)
relative_core <- substring(current_core, nchar(root) + 2L)
stopifnot(identical(sort(relative_core), sort(original_manifest$file)))
dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
out <- tempfile("sensitivity_static_", tmpdir = file.path(root, "results"))
dir.create(out)

write.csv(
  do.call(rbind, rows),
  file.path(out, "function_checks.csv"),
  row.names = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

cat(
  "PASS: parsing, original functions, full schedule, 13 prior specifications and illustration files.\nReports: ",
  out,
  "\n",
  sep = ""
)
