## Tests how functions are read from the original scripts, including functions
## defined inside if/else blocks. Needs no refit packages and no archive.
local({
  root <- getwd()
  for (i in 1:20) {
    if (
      file.exists(file.path(
        root,
        "simulation_nonlinear_SCAR_all18/R/workflow.R"
      ))
    ) {
      break
    }
    up <- dirname(root)
    if (identical(up, root)) {
      stop("Open GP-PU-PoE.Rproj before running this test.")
    }
    root <- up
  }
  wf <- new.env(parent = .GlobalEnv)
  sys.source(
    file.path(root, "simulation_nonlinear_SCAR_all18/R/workflow.R"),
    envir = wf
  )
  e <- new.env(parent = .GlobalEnv)
  wf$scar_load(root, e)

  ## The condition and both branches stop if they are evaluated. Only the
  ## function definitions may be read.
  fixture <- parse(
    text = c(
      'plain <- function(x=3L) x + 1L',
      'if (stop("The condition must not be evaluated")) {',
      ' stop("The then block must not be executed")',
      ' in_if <- function(x) x * 2L',
      '} else {',
      ' stop("The else block must not be executed")',
      ' in_else <- function(x=5L) { local_only <- function(y) y; local_only(x) }',
      '}',
      'holder <- function() { in_else <- function(x) -x; nested <- function() 1L }',
      'local({ concealed <- function() 1L; stop("Do not execute local") })',
      'for (j in stop("Do not execute a loop")) { loop_function <- function() j }',
      'quote(quoted_function <- function() 1L)',
      '{ in_block <- function() 9L }'
    ),
    keep.source = FALSE
  )
  index <- e$scar_script_function_index(fixture)
  expected <- c("plain", "in_if", "in_else", "holder", "in_block")
  e$scar_assert(
    setequal(names(index), expected),
    "Function lookup entered a function body/loop/other call or lost an allowed definition."
  )
  env <- new.env(parent = baseenv())
  for (name in expected) {
    x <- e$scar_original_function(index, name, "fixture")
    eval(x, env)
    f <- get(name, envir = env, inherits = FALSE)
    e$scar_assert(
      is.function(f) &&
        identical(formals(f), x[[3L]][[2L]]) &&
        identical(body(f), x[[3L]][[3L]]),
      paste("Extracted definition changed:", name)
    )
  }
  e$scar_assert(
    identical(env$plain(), 4L) &&
      identical(env$in_if(3L), 6L) &&
      identical(env$in_else(), 5L) &&
      identical(env$in_block(), 9L),
    "Extracted fixture functions did not behave as specified."
  )
  absent <- inherits(
    try(e$scar_original_function(index, "missing", "fixture"), silent = TRUE),
    "try-error"
  )
  e$scar_assert(absent, "A missing definition was not reported.")
  ambiguous <- e$scar_script_function_index(parse(
    text = 'if (TRUE) dup <- function() 1L else dup <- function() 2L',
    keep.source = FALSE
  ))
  rejected <- inherits(
    try(e$scar_original_function(ambiguous, "dup", "fixture"), silent = TRUE),
    "try-error"
  )
  e$scar_assert(rejected, "An ambiguous definition was not reported.")

  ## The eight original scripts: static checks, then the function definitions.
  e$scar_static_checks(root)
  mod <- e$scar_source_root(root)
  manifest <- utils::read.csv(
    file.path(mod, "reference/SOURCE_MANIFEST.csv"),
    stringsAsFactors = FALSE
  )
  total <- 0L
  for (i in seq_len(nrow(manifest))) {
    group <- manifest$group[i]
    expressions <- parse(
      file.path(mod, manifest$reference_file[i]),
      keep.source = FALSE
    )
    index <- e$scar_script_function_index(expressions)
    fns <- readLines(
      file.path(mod, "reference", paste0(group, "_functions.txt")),
      warn = FALSE
    )
    env <- new.env(parent = baseenv())
    for (name in fns) {
      x <- e$scar_original_function(index, name, group)
      eval(x, env)
      f <- get(name, envir = env, inherits = FALSE)
      e$scar_assert(
        is.function(f) &&
          identical(formals(f), x[[3L]][[2L]]) &&
          identical(body(f), x[[3L]][[3L]]),
        paste("Original function changed:", group, name)
      )
      total <- total + 1L
    }
    if (group %in% c("GP_PU_PoE", "GP_PU_PoE_beta0")) {
      top_names <- vapply(expressions, e$scar_assignment_name, character(1))
      e$scar_assert(
        !("format_mean_sd" %in% top_names),
        paste(
          group,
          "no longer reproduces the original nested-definition case."
        )
      )
      e$scar_assert(
        identical(env$format_mean_sd(0.125, 0.25), "0.125 (0.250)") &&
          identical(env$format_mean_sd(NA_real_, 0.25), ""),
        paste(group, "format_mean_sd did not match the original definition.")
      )
    }
  }
  e$scar_assert(
    total == 480L,
    "Expected the 480 functions of the original scripts."
  )
  cat("PASS: function lookup and the 480 functions of the original scripts.\n")
})
