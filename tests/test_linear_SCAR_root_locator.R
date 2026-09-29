## Tests how scripts 07 and 08 find the project folder. Run from the project
## folder with source("tests/test_linear_SCAR_root_locator.R").
local({
  entry_files <- file.path(
    "scripts",
    c("07_linear_SCAR_quick.R", "08_linear_SCAR_full.R")
  )
  if (!all(file.exists(entry_files))) {
    stop("Open GP-PU-PoE.Rproj before running this test.")
  }
  extract_locator <- function(path) {
    exprs <- parse(path, keep.source = FALSE)
    is_local <- vapply(
      exprs,
      function(x) {
        is.call(x) && identical(x[[1L]], as.name("local"))
      },
      logical(1)
    )
    stopifnot(sum(is_local) == 1L)
    block <- exprs[[which(is_local)]][[2L]]
    assignments <- as.list(block)[-1L]
    hit <- vapply(
      assignments,
      function(x) {
        is.call(x) &&
          identical(x[[1L]], as.name("<-")) &&
          identical(x[[2L]], as.name("locate_root"))
      },
      logical(1)
    )
    stopifnot(sum(hit) == 1L)
    e <- new.env(parent = baseenv())
    eval(assignments[[which(hit)]], envir = e)
    e$locate_root
  }
  f <- extract_locator(entry_files[[1L]])
  g <- extract_locator(entry_files[[2L]])
  stopifnot(identical(formals(f), formals(g)), identical(body(f), body(g)))

  initial_wd <- getwd()
  temp <- tempfile("linear_scar_root_checks_")
  dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE, force = TRUE), add = TRUE)
  make_project <- function(name) {
    p <- file.path(temp, name)
    dir.create(
      file.path(p, "simulation_linear_SCAR_all18", "R"),
      recursive = TRUE
    )
    dir.create(file.path(p, "scripts"))
    writeLines(
      "# Test file.",
      file.path(p, "simulation_linear_SCAR_all18", "R", "workflow.R")
    )
    writeLines("Version: 1.0", file.path(p, "GP-PU-PoE.Rproj"))
    normalizePath(p, winslash = "/", mustWork = TRUE)
  }
  a <- make_project("project with spaces")
  b <- make_project("second project")
  elsewhere <- file.path(temp, "elsewhere")
  dir.create(elsewhere)
  expected_error <- function(expr, text) {
    err <- tryCatch(
      {
        force(expr)
        NULL
      },
      error = identity
    )
    stopifnot(
      inherits(err, "error"),
      grepl(text, conditionMessage(err), fixed = TRUE)
    )
    invisible(TRUE)
  }
  stopifnot(identical(f(a, starts = elsewhere, interactive_ok = FALSE), a))
  stopifnot(identical(
    f(
      file.path(a, "GP-PU-PoE.Rproj"),
      starts = elsewhere,
      interactive_ok = FALSE
    ),
    a
  ))
  stopifnot(identical(f(starts = a, interactive_ok = FALSE), a))
  stopifnot(identical(
    f(starts = file.path(a, "scripts"), interactive_ok = FALSE),
    a
  ))
  stopifnot(identical(
    f(starts = c(a, file.path(a, "scripts")), interactive_ok = FALSE),
    a
  ))
  expected_error(
    f(elsewhere, starts = a, interactive_ok = FALSE),
    "Missing file:"
  )
  expected_error(
    f(starts = elsewhere, interactive_ok = FALSE),
    "Set PROJECT_DIR"
  )
  expected_error(f(starts = c(a, b), interactive_ok = FALSE), "Set PROJECT_DIR")
  stopifnot(identical(
    f(starts = elsewhere, interactive_ok = TRUE, chooser = function() {
      file.path(a, "GP-PU-PoE.Rproj")
    }),
    a
  ))
  expected_error(
    f(starts = elsewhere, interactive_ok = TRUE, chooser = function() {
      file.path(a, "scripts", "08_linear_SCAR_full.R")
    }),
    "Select the .Rproj file"
  )
  expected_error(
    f(starts = elsewhere, interactive_ok = TRUE, chooser = function() {
      stop("cancelled")
    }),
    "No project was selected"
  )
  stopifnot(identical(getwd(), initial_wd))
  cat("PASS: scripts 07 and 08 find the project folder.\n")
})
