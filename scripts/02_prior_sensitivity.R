## Prior sensitivity analysis of the illustrative example: 13 prior settings.
## RStudio: open GP-PU-PoE.Rproj, set MODE and ACTION, then Source this file.
MODE <- "quick" # "quick" checks that the code runs; "full" is the schedule of the paper
ACTION <- "all" # "all", "plot", "check", "prepare", "fit", "diagnostics"

local({
  locate_root <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    command_file <- sub("^--file=", "", args[grepl("^--file=", args)])
    source_files <- unlist(
      lapply(sys.frames(), function(fr) {
        x <- fr$ofile
        if (is.character(x) && length(x) == 1L) x else NULL
      }),
      use.names = FALSE
    )
    candidates <- unique(c(getwd(), dirname(c(source_files, command_file))))
    for (candidate in candidates) {
      path <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
      for (i in seq_len(20L)) {
        if (
          file.exists(file.path(path, ".here")) &&
            file.exists(file.path(path, "sensitivity", "R", "workflow.R"))
        ) {
          return(path)
        }
        parent <- dirname(path)
        if (identical(parent, path)) {
          break
        }
        path <- parent
      }
    }
    stop(
      "Cannot find the project folder. Open GP-PU-PoE.Rproj or set the working directory to the project folder."
    )
  }
  root <- locate_root()
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat(
      "Rscript --vanilla scripts/02_prior_sensitivity.R --mode=quick --action=all\n",
      "Options: --mode=quick|full; --action=all|plot|check|prepare|fit|diagnostics;\n",
      "         --out=PATH; --fit-ids=DEFAULT,ETA_DIFFUSE (fit action only).\n",
      "Quick mode only checks that the code runs.\n",
      sep = ""
    )
  } else {
    unknown <- args[!grepl("^--(mode|action|out|fit-ids)=", args)]
    if (length(unknown)) {
      stop("Unknown arguments: ", paste(unknown, collapse = " "))
    }
    option_value <- function(name, default) {
      hit <- args[startsWith(args, paste0("--", name, "="))]
      if (!length(hit)) {
        return(default)
      }
      if (length(hit) != 1L) {
        stop("Repeated option --", name)
      }
      value <- substring(hit, nchar(name) + 4L)
      if (!nzchar(value)) {
        stop("Empty value for --", name)
      }
      value
    }
    fit_string <- option_value("fit-ids", NULL)
    fit_ids <- if (is.null(fit_string)) {
      NULL
    } else {
      strsplit(fit_string, ",", fixed = TRUE)[[1L]]
    }
    wf <- new.env(parent = as.environment("package:stats"))
    sys.source(file.path(root, "sensitivity", "R", "workflow.R"), envir = wf)
    wf$run_prior_sensitivity(
      root = root,
      mode = option_value("mode", MODE),
      action = option_value("action", ACTION),
      out_dir = option_value("out", NULL),
      fit_ids = fit_ids
    )
  }
})
