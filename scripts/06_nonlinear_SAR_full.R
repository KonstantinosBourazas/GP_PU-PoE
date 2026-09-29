## Nonlinear SAR study: rebuilds the tables of the paper from the 100
## archived replications of the 18 methods, without MCMC.
ACTION <- "all" # "all" or "check"
OUTPUT_DIR <- NULL
ARCHIVE_DIR <- NULL

local({
  locate_root <- function() {
    args0 <- commandArgs(trailingOnly = FALSE)
    cli <- sub("^--file=", "", args0[grepl("^--file=", args0)])
    files <- unlist(
      lapply(sys.frames(), function(f) {
        if (is.character(f$ofile) && length(f$ofile) == 1L) f$ofile else NULL
      }),
      use.names = FALSE
    )
    for (start in unique(c(getwd(), dirname(c(files, cli))))) {
      p <- normalizePath(start, winslash = "/", mustWork = FALSE)
      for (i in 1:20) {
        if (
          file.exists(file.path(
            p,
            "simulation_nonlinear_SAR_all18/R/workflow.R"
          ))
        ) {
          return(p)
        }
        up <- dirname(p)
        if (identical(up, p)) {
          break
        }
        p <- up
      }
    }
    stop(
      "Cannot find the project folder. Open GP-PU-PoE.Rproj or set the working directory to the project folder."
    )
  }
  root <- locate_root()
  args <- commandArgs(trailingOnly = TRUE)
  bad <- args[!grepl("^--(action|out|archive)=", args)]
  if (length(bad)) {
    stop("Unknown command-line arguments: ", paste(bad, collapse = ", "))
  }
  opt <- function(key, default) {
    hit <- args[startsWith(args, paste0("--", key, "="))]
    if (!length(hit)) {
      return(default)
    }
    if (length(hit) != 1L) {
      stop("Duplicate option: ", key)
    }
    value <- substring(hit, nchar(key) + 4L)
    if (!nzchar(value)) {
      stop("Empty option: ", key)
    }
    value
  }
  wf <- new.env(parent = parent.env(.GlobalEnv))
  sys.source(
    file.path(root, "simulation_nonlinear_SAR_all18/R/workflow.R"),
    envir = wf
  )
  wf$run_nonlinear_SAR_all18(
    root,
    mode = "archived_full",
    action = opt("action", ACTION),
    out = opt("out", OUTPUT_DIR),
    archive = opt("archive", ARCHIVE_DIR),
    allow_long_MCMC = FALSE
  )
})
