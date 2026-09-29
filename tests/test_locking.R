## Tests the output lock of the illustration workflow.
## Run from the project root; also called by test_reproducibility.R.
test_illustration_locking <- function(root = getwd()) {
  stopifnot(file.exists(file.path(root, ".here")))
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  wf <- new.env(parent = as.environment("package:stats"))
  sys.source(file.path(root, "R", "illustration", "workflow.R"), envir = wf)
  wf$check_illustration_dependencies()
  dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
  test_base <- tempfile("lock_check_", tmpdir = file.path(root, "results"))
  ## A path with spaces and brackets, which must be treated literally.
  out <- file.path(test_base, "path with spaces [literal]")
  dir.create(out, recursive = TRUE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)
  lock <- file.path(out, ".run_lock")
  sentinel <- file.path(out, "result.txt")
  writeLines("result file", sentinel)
  sentinel_hash <- unname(tools::md5sum(sentinel))
  report <- character(0)
  check_unlocked <- function() {
    if (dir.exists(lock) || file.exists(lock)) {
      stop("Lock was not released: ", lock, call. = FALSE)
    }
    stopifnot(identical(unname(tools::md5sum(sentinel)), sentinel_hash))
  }
  run_check <- function() {
    invisible(capture.output(
      wf$run_illustration(root, mode = "quick", action = "check", out_dir = out)
    ))
  }
  ## Success cleanup, including a second action in the same directory.
  run_check()
  check_unlocked()
  run_check()
  check_unlocked()
  report <- c(report, "two_check_actions_release_lock")

  ## Cleanup after an error: this folder has no fitted summary.
  failure <- tryCatch(
    {
      invisible(capture.output(
        wf$run_illustration(
          root,
          mode = "quick",
          action = "plot",
          out_dir = out
        )
      ))
      NULL
    },
    error = function(err) err
  )
  stopifnot(
    inherits(failure, "error"),
    grepl(
      "No completed summary bundle",
      conditionMessage(failure),
      fixed = TRUE
    )
  )
  check_unlocked()
  report <- c(report, "missing_summary_error_releases_lock")

  ## An existing lock is neither removed nor bypassed.
  stopifnot(dir.create(lock))
  owner_file <- file.path(lock, "owner.txt")
  writeLines("test fixture representing another active run", owner_file)
  owner_hash <- unname(tools::md5sum(owner_file))
  failure <- tryCatch(
    {
      run_check()
      NULL
    },
    error = function(err) err
  )
  stopifnot(
    inherits(failure, "error"),
    grepl(
      "Output directory is locked",
      conditionMessage(failure),
      fixed = TRUE
    ),
    dir.exists(lock),
    identical(unname(tools::md5sum(owner_file)), owner_hash)
  )
  ## Only remove the lock just created by this test.
  stopifnot(isTRUE(wf$release_illustration_lock(lock)))
  check_unlocked()
  report <- c(report, "existing_lock_is_not_bypassed")

  ## Simulate two temporary unlink failures, then allow actual deletion.
  unlink_calls <- 0L
  wf$unlink <- function(x, recursive, force, expand) {
    unlink_calls <<- unlink_calls + 1L
    stopifnot(
      isTRUE(recursive),
      isTRUE(force),
      identical(expand, FALSE),
      identical(x, lock)
    )
    if (unlink_calls <= 2L) {
      return(1L)
    }
    base::unlink(x, recursive = recursive, force = force, expand = expand)
  }
  run_check()
  stopifnot(identical(unlink_calls, 3L))
  rm("unlink", envir = wf)
  check_unlocked()
  report <- c(report, "temporary_cleanup_failures_are_retried")

  ## A persistent failure is reported.
  stopifnot(dir.create(lock))
  writeLines("persistent-failure test fixture", file.path(lock, "owner.txt"))
  wf$unlink <- function(...) 1L
  warning_text <- character(0)
  released <- withCallingHandlers(
    wf$release_illustration_lock(lock, attempts = 2L, delay = 0),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(
    identical(released, FALSE),
    length(warning_text) == 1L,
    grepl("Could not remove this run's lock", warning_text, fixed = TRUE),
    dir.exists(lock)
  )
  rm("unlink", envir = wf)
  stopifnot(isTRUE(wf$release_illustration_lock(lock)))
  check_unlocked()
  report <- c(
    report,
    "persistent_cleanup_failure_is_reported",
    "result_file_is_unchanged"
  )

  write.csv(
    data.frame(check = report, passed = TRUE),
    file.path(test_base, "lock_checks.csv"),
    row.names = FALSE
  )
  cat("PASS: output lock. Reports: ", test_base, "\n", sep = "")
  invisible(test_base)
}

test_illustration_locking(getwd())
