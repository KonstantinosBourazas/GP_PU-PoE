## Prior sensitivity analysis over 100 datasets. Checks the 1,300 archived fits
## in archive/prior_sensitivity_mc100/ and rebuilds the tables and the two
## figures, without MCMC. The figures keep the original order of the units. The
## script also writes paired AUC and rank summaries.
## RStudio: open GP-PU-PoE.Rproj, then Source this file. After a first run with
## ACTION "all", ACTION "plot" only redraws the figures.
ACTION <- "all" # "all", "plot", "check", "combine"
ARCHIVE_DIR <- NULL # NULL uses archive/prior_sensitivity_mc100
OUTPUT_DIR <- NULL # NULL uses results/prior_sensitivity/mc100

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
    for (candidate in unique(c(
      getwd(),
      dirname(c(source_files, command_file))
    ))) {
      path <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
      for (i in seq_len(20L)) {
        if (
          file.exists(file.path(path, "sensitivity_mc100", "R", "workflow.R"))
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
      "Rscript --vanilla scripts/02b_prior_sensitivity_MC100.R --action=all\n",
      "--action=all|plot|check|combine; --archive=PATH; --out=PATH\n",
      "The script reads the archived fits and runs no MCMC.\n",
      sep = ""
    )
  } else {
    unknown <- args[!grepl("^--(action|archive|out)=", args)]
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
    wf <- new.env(parent = as.environment("package:stats"))
    sys.source(
      file.path(root, "sensitivity_mc100", "R", "workflow.R"),
      envir = wf
    )
    wf$run_prior_sensitivity_mc100(
      root = root,
      action = option_value("action", ACTION),
      archive = option_value("archive", ARCHIVE_DIR),
      out = option_value("out", OUTPUT_DIR)
    )
  }
})
