## Tests the reader of the failure records in the linear SCAR archive.
local({
  args0 <- commandArgs(trailingOnly = FALSE)
  cli <- sub("^--file=", "", args0[grepl("^--file=", args0)])
  files <- unlist(
    lapply(sys.frames(), function(f) {
      x <- f$ofile
      if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
        x
      } else {
        NULL
      }
    }),
    use.names = FALSE
  )
  root <- NULL
  for (start in unique(c(getwd(), dirname(c(files, cli))))) {
    if (!dir.exists(start)) {
      next
    }
    p <- normalizePath(start, winslash = "/", mustWork = TRUE)
    for (i in seq_len(50L)) {
      if (
        file.exists(file.path(p, "simulation_linear_SCAR_all18/R/archive.R"))
      ) {
        root <- p
        break
      }
      up <- dirname(p)
      if (identical(up, p)) {
        break
      }
      p <- up
    }
    if (!is.null(root)) break
  }
  if (is.null(root)) {
    stop("Open GP-PU-PoE.Rproj before running this test.", call. = FALSE)
  }
  e <- new.env(parent = parent.env(.GlobalEnv))
  for (f in c("utilities.R", "registry.R", "archive.R")) {
    sys.source(file.path(root, "simulation_linear_SCAR_all18/R", f), envir = e)
  }
  must_error <- function(expr, expected_text) {
    error <- tryCatch(
      {
        force(expr)
        NULL
      },
      error = function(err) err
    )
    if (
      !inherits(error, "error") ||
        !grepl(expected_text, conditionMessage(error), fixed = TRUE)
    ) {
      stop("An invalid failure record was not rejected.", call. = FALSE)
    }
    invisible(TRUE)
  }
  # nrow() of an empty list is NULL, so a plain assertion fails on it.
  must_error(e$lscar_assert(nrow(list()) == 0L, "EMPTY_LIST"), "EMPTY_LIST")
  accepted <- list(
    NULL,
    list(),
    data.frame(),
    data.frame(error_message = character())
  )
  for (x in accepted) {
    stopifnot(isTRUE(e$lscar_assert_empty_method_failures(x, "empty test")))
  }
  # A nonempty list is a failure record, even when its elements are empty.
  rejected <- list(
    list(WLR = data.frame(error_message = "synthetic failure")),
    data.frame(error_message = "synthetic failure"),
    list(WLR = NULL),
    list(WLR = data.frame())
  )
  for (x in rejected) {
    must_error(
      e$lscar_assert_empty_method_failures(x, "negative test"),
      "saved method failures"
    )
  }
  malformed <- list(
    character(),
    numeric(),
    matrix(character(), nrow = 0L, ncol = 1L)
  )
  for (x in malformed) {
    must_error(
      e$lscar_assert_empty_method_failures(x, "malformed test"),
      "invalid method_failures format"
    )
  }
  cat(
    "PASS: empty NULL/list/data-frame cases accepted; real and malformed failure cases rejected.\n"
  )

  archive <- file.path(root, "archive/simulations/linear_SCAR_all18")
  specs <- e$lscar_group_spec()
  checked <- 0L
  for (k in seq_len(nrow(specs))) {
    sp <- specs[k, , drop = FALSE]
    for (r in seq_len(100L)) {
      path <- file.path(
        archive,
        sp$group,
        sp$checkpoint_dir,
        sprintf(sp$pattern, r)
      )
      x <- e$lscar_read(path)
      e$lscar_assert(
        isTRUE(x$complete) &&
          identical(as.integer(x$replicate), as.integer(r)),
        paste(sp$group, r, "incomplete checkpoint")
      )
      e$lscar_assert_empty_method_failures(
        x[["method_failures", exact = TRUE]],
        paste(sp$group, r)
      )
      e$lscar_assert(
        is.null(x[["failure", exact = TRUE]]),
        paste(sp$group, r, "saved failure")
      )
      checked <- checked + 1L
    }
  }
  stopifnot(checked == 800L)
  cat("PASS: failure fields in all 800 archived production checkpoints.\n")
})
