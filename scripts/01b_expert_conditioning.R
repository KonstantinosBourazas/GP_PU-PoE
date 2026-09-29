## Expert decomposition of the GP-PU-PoE fit of the illustrative example.
## Reads a saved fit and runs no MCMC.
PROJECT_DIR <- NULL
MODE <- "full" # "full" or "quick": which illustration fit to read
ACTION <- "all" # "all": scores and figure; "plot": figure from saved scores; "check": inputs
INPUT_DIR <- NULL # Default: archive/illustration/ (full) or results/illustration/quick/ (quick)
OUTPUT_DIR <- NULL # Default: results/illustration/expert_conditioning/

local({
  args <- commandArgs(trailingOnly = TRUE)
  allowed <- "^--(project|mode|action|input|out)="
  if (any(!grepl(allowed, args))) {
    stop("Unknown command-line option.", call. = FALSE)
  }
  opt <- function(key, fallback) {
    hit <- args[startsWith(args, paste0("--", key, "="))]
    if (!length(hit)) {
      return(fallback)
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
  marker <- "illustration_expert_conditioning/R/workflow.R"
  is_root <- function(p) {
    file.exists(file.path(p, ".here")) &&
      file.exists(file.path(p, marker))
  }
  selected_root <- function(p) {
    if (
      file.exists(p) &&
        !dir.exists(p) &&
        grepl("[.]Rproj$", p, ignore.case = TRUE)
    ) {
      p <- dirname(p)
    }
    p <- normalizePath(path.expand(p), winslash = "/", mustWork = TRUE)
    if (!is_root(p)) {
      stop(
        "The selected folder does not contain the expert conditioning module:\n",
        p,
        call. = FALSE
      )
    }
    p
  }
  project <- opt("project", PROJECT_DIR)
  if (!is.null(project)) {
    root <- selected_root(project)
  } else {
    cli <- commandArgs(trailingOnly = FALSE)
    cli <- sub("^--file=", "", cli[grepl("^--file=", cli)])
    sourced <- unlist(
      lapply(sys.frames(), function(f) {
        if (is.character(f$ofile) && length(f$ofile) == 1L) f$ofile else NULL
      }),
      use.names = FALSE
    )
    roots <- character()
    for (start in unique(c(getwd(), dirname(c(sourced, cli))))) {
      if (!dir.exists(start)) {
        next
      }
      p <- normalizePath(start, winslash = "/", mustWork = TRUE)
      for (i in seq_len(50L)) {
        if (is_root(p)) {
          roots <- c(roots, p)
          break
        }
        up <- dirname(p)
        if (identical(up, p)) {
          break
        }
        p <- up
      }
    }
    roots <- unique(roots)
    if (length(roots) == 1L) {
      root <- roots[[1L]]
    } else {
      if (!interactive()) {
        stop("Set PROJECT_DIR or pass --project=PATH.", call. = FALSE)
      }
      cat("Select GP-PU-PoE.Rproj.\n")
      chosen <- file.choose()
      if (!grepl("[.]Rproj$", chosen, ignore.case = TRUE)) {
        stop("Select the .Rproj file of the project.")
      }
      root <- selected_root(chosen)
    }
  }
  wf <- new.env(parent = as.environment("package:stats"))
  sys.source(file.path(root, marker), envir = wf)
  wf$run_expert_conditioning(
    root,
    mode = opt("mode", MODE),
    action = opt("action", ACTION),
    input_dir = opt("input", INPUT_DIR),
    out_dir = opt("out", OUTPUT_DIR)
  )
})
