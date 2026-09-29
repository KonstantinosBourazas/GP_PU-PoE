lscar_static_checks <- function(root) {
  mod <- lscar_source_root(root)
  for (f in list.files(
    mod,
    pattern = "[.]R$",
    full.names = TRUE,
    recursive = TRUE
  )) {
    parse(file = f, keep.source = FALSE)
  }
  manifest <- utils::read.csv(
    file.path(mod, "reference/SOURCE_MANIFEST.csv"),
    stringsAsFactors = FALSE
  )
  lscar_assert(
    setequal(manifest$group, lscar_groups()),
    "The source registry does not cover all eight groups."
  )
  forbidden <- c(
    "load",
    "source",
    "sys.source",
    "readRDS",
    "saveRDS",
    "write.csv",
    "unlink",
    "dir.create",
    "system",
    "system2"
  )
  calls <- function(x) {
    if (!is.call(x)) {
      return(character())
    }
    head <- if (is.symbol(x[[1L]])) as.character(x[[1L]]) else ""
    unique(c(head, unlist(lapply(as.list(x)[-1L], calls), use.names = FALSE)))
  }
  total <- 0L
  for (i in seq_len(nrow(manifest))) {
    g <- manifest$group[i]
    p <- file.path(mod, manifest$reference_file[i])
    lscar_assert(
      identical(lscar_md5(p), as.character(manifest$md5[i])),
      paste("Reference changed", g)
    )
    expr <- parse(p, keep.source = FALSE)
    nm <- vapply(expr, lscar_assignment_name, character(1))
    function_index <- lscar_script_function_index(expr)
    fns <- readLines(
      file.path(mod, "reference", paste0(g, "_functions.txt")),
      warn = FALSE
    )
    settings <- readLines(
      file.path(mod, "reference", paste0(g, "_settings.txt")),
      warn = FALSE
    )
    for (n in fns) {
      # Same function lookup as in the import.
      invisible(lscar_original_function(function_index, n, g))
    }
    lscar_assert(all(settings %in% nm), paste("Setting not found", g))
    for (k in which(nm %in% settings)) {
      lscar_assert(
        !length(intersect(calls(expr[[k]][[3L]]), forbidden)),
        paste("Imported setting calls a file or system function", g, nm[k])
      )
    }
    total <- total + length(fns)
  }
  lscar_assert(
    nrow(lscar_registry()) == 18L && !anyDuplicated(lscar_registry()$Method),
    "Expected 18 distinct methods."
  )
  cat(
    "PASS: parsing, ",
    total,
    " functions of the original scripts, imported settings, linear priors and the 18 methods.\n",
    sep = ""
  )
  lscar_linear_prior_selftest(root)
  invisible(TRUE)
}
