sar_assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

sar_near <- function(x, y, tolerance = 1e-9) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  length(x) == length(y) &&
    identical(is.na(x), is.na(y)) &&
    all(is.finite(x[!is.na(x)])) &&
    all(is.finite(y[!is.na(y)])) &&
    all(abs(x[!is.na(x)] - y[!is.na(y)]) <= tolerance)
}

sar_bind <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z) > 0L, x)
  sar_assert(length(x) > 0L, "No result rows to combine.")
  nm <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(d) {
    for (k in setdiff(nm, names(d))) {
      d[[k]] <- NA
    }
    d[, nm, drop = FALSE]
  })
  d <- do.call(rbind, x)
  rownames(d) <- NULL
  d
}

sar_read <- function(p) {
  if (!file.exists(p)) {
    stop("Missing archived result: ", p, call. = FALSE)
  }
  tryCatch(readRDS(p), error = function(e) {
    stop("Cannot read ", p, "\n", conditionMessage(e), call. = FALSE)
  })
}

sar_md5 <- function(p) unname(tools::md5sum(p))

sar_atomic_rds <- function(x, p) {
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".writing_", tmpdir = dirname(p))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (file.exists(p)) {
    sar_assert(file.remove(p), paste("Cannot replace", p))
  }
  sar_assert(file.rename(tmp, p), paste("Cannot finish writing", p))
  invisible(p)
}

sar_csv <- function(d, p) {
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(d, p, row.names = FALSE, na = "NA")
}

sar_source_root <- function(root) {
  file.path(root, "simulation_nonlinear_SAR_all18")
}

sar_source_signature <- function(root) {
  m <- sar_source_root(root)
  f <- sort(unlist(
    lapply(
      file.path(m, c("R", "config", "reference")),
      list.files,
      full.names = TRUE,
      recursive = TRUE
    ),
    use.names = FALSE
  ))
  stats::setNames(sar_md5(f), substring(f, nchar(m) + 2L))
}

sar_lock <- function(out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  sar_assert(dir.exists(out), paste("Cannot create output directory", out))
  lock <- file.path(out, ".run_lock")
  if (dir.exists(lock)) {
    stop(
      "Output directory is locked: ",
      out,
      "\nRemove .run_lock by hand once no run is active.",
      call. = FALSE
    )
  }
  sar_assert(
    dir.create(lock, showWarnings = FALSE),
    paste("Unable to create lock", lock)
  )
  writeLines(
    c(paste("PID", Sys.getpid()), format(Sys.time())),
    file.path(lock, "owner.txt")
  )
  lock
}

sar_unlock <- function(lock) {
  for (i in 1:8) {
    unlink(lock, recursive = TRUE, force = TRUE, expand = FALSE)
    if (!dir.exists(lock)) {
      return(invisible(TRUE))
    }
    Sys.sleep(.1 * i)
  }
  warning(
    "Cannot release output lock: ",
    lock,
    ". Remove it by hand once no run is active."
  )
  invisible(FALSE)
}

sar_safe_output <- function(root, out) {
  if (is.null(out) || !nzchar(out)) {
    stop("An output directory is required.")
  }
  clean <- function(x) {
    tolower(gsub(
      "\\\\",
      "/",
      normalizePath(x, winslash = "/", mustWork = FALSE)
    ))
  }
  a <- clean(file.path(root, "archive"))
  o <- clean(out)
  sar_assert(
    o != a && !startsWith(paste0(o, "/"), paste0(a, "/")),
    "The output folder cannot be inside archive/."
  )
  invisible(out)
}

sar_log_start <- function(out, name) {
  con <- file(file.path(out, name), open = "wt")
  sink(con, split = TRUE)
  list(connection = con, depth = sink.number())
}

sar_log_end <- function(x) {
  while (sink.number() >= x$depth) {
    sink()
  }
  close(x$connection)
}

sar_auc <- function(p, t) {
  n1 <- sum(t == 1L)
  n0 <- sum(t == 0L)
  sar_assert(n1 > 0 && n0 > 0, "AUC needs both classes.")
  (sum(rank(p, ties.method = "average")[t == 1L]) - n1 * (n1 + 1) / 2) /
    (n1 * n0)
}

sar_format <- function(mean, sd) {
  ifelse(is.na(mean), "--", sprintf("%.3f (%.3f)", mean, sd))
}
