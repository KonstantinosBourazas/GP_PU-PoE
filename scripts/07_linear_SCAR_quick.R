## Linear SCAR study: fits all 18 methods on two new datasets with short runs.
## This checks that the code runs; the estimates are not those of the paper.

PROJECT_DIR <- NULL # Optional: the project folder or its .Rproj file
ACTION <- "all" # "all" or "check"
OUTPUT_DIR <- NULL
ARCHIVE_DIR <- NULL

## Leave PROJECT_DIR as NULL when the project is open in RStudio. If the
## project folder is not found, a dialog asks for its .Rproj file.

local({
  ## Finds the project folder.
  locate_root <- function(
    project_dir = NULL,
    starts = getwd(),
    interactive_ok = interactive(),
    chooser = file.choose
  ) {
    marker <- "simulation_linear_SCAR_all18/R/workflow.R"
    normalize <- function(p) {
      normalizePath(path.expand(p), winslash = "/", mustWork = TRUE)
    }
    is_root <- function(p) {
      dir.exists(p) &&
        file.exists(file.path(p, marker)) &&
        !dir.exists(file.path(p, marker))
    }
    exact_root <- function(selection, project_file_only = FALSE) {
      if (
        !is.character(selection) ||
          length(selection) != 1L ||
          is.na(selection) ||
          !nzchar(selection)
      ) {
        stop(
          "Select an existing project folder or its .Rproj file.",
          call. = FALSE
        )
      }
      selection <- path.expand(selection)
      if (
        project_file_only &&
          (!file.exists(selection) ||
            dir.exists(selection) ||
            !grepl("[.]Rproj$", selection, ignore.case = TRUE))
      ) {
        stop("Select the .Rproj file, not an .R script.", call. = FALSE)
      }
      if (dir.exists(selection)) {
        p <- normalize(selection)
      } else if (
        file.exists(selection) &&
          grepl("[.]Rproj$", selection, ignore.case = TRUE)
      ) {
        p <- dirname(normalize(selection))
      } else {
        stop(
          "PROJECT_DIR is not an existing folder or .Rproj file:\n",
          selection,
          call. = FALSE
        )
      }
      if (!is_root(p)) {
        stop(
          "The selected folder does not contain the linear SCAR module.\n",
          "Selected folder: ",
          p,
          "\nMissing file: ",
          file.path(p, marker),
          "\nSelect the project folder that directly contains scripts/,",
          " simulation_linear_SCAR_all18/ and archive/.",
          call. = FALSE
        )
      }
      p
    }

    ## A folder given explicitly is used as it is.
    if (!is.null(project_dir)) {
      return(exact_root(project_dir))
    }

    starts <- unique(as.character(starts))
    starts <- starts[!is.na(starts) & nzchar(starts)]
    hits <- character()
    for (start in starts) {
      if (!dir.exists(start)) {
        next
      }
      p <- normalize(start)
      for (i in seq_len(50L)) {
        if (is_root(p)) {
          hits <- c(hits, p)
          break
        }
        up <- dirname(p)
        if (identical(up, p)) {
          break
        }
        p <- up
      }
    }
    hits <- unique(hits)
    if (length(hits) == 1L) {
      return(hits[[1L]])
    }

    cat("\nThe project folder was not found automatically.\n")
    cat("Current working directory: ", getwd(), "\n", sep = "")
    if (length(hits) > 1L) {
      cat(
        "Several project folders were found:\n",
        paste(hits, collapse = "\n"),
        "\n",
        sep = ""
      )
    } else {
      cat("Expected relative path: ", marker, "\n", sep = "")
    }
    if (!isTRUE(interactive_ok)) {
      stop(
        "Set PROJECT_DIR at the top of this script or pass --project=PATH.",
        call. = FALSE
      )
    }
    cat("Select GP-PU-PoE.Rproj in the file dialog.\n")
    selected <- tryCatch(chooser(), error = function(e) {
      stop(
        "No project was selected. Set PROJECT_DIR and run the script again.",
        call. = FALSE
      )
    })
    exact_root(selected, project_file_only = TRUE)
  }

  ## Command line options, including --project.
  args <- commandArgs(trailingOnly = TRUE)
  bad <- args[!grepl("^--(action|out|archive|project)=", args)]
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
  action <- match.arg(opt("action", ACTION), c("all", "check"))

  ## The path of this script is known when it runs with Source or Rscript.
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
  starts <- unique(c(dirname(c(files, cli)), getwd()))
  root <- locate_root(opt("project", PROJECT_DIR), starts = starts)

  ## Check that the files of the module are present.
  needed <- c(
    paste0(
      "simulation_linear_SCAR_all18/R/",
      c(
        "workflow.R",
        "utilities.R",
        "registry.R",
        "archive.R",
        "tables.R",
        "engines.R",
        "checks.R"
      )
    ),
    "simulation_linear_SCAR_all18/config/linear_priors.R",
    "simulation_linear_SCAR_all18/reference/SOURCE_MANIFEST.csv",
    "simulation_linear_SCAR_all18/reference/ARCHIVE_MANIFEST.csv"
  )
  absent <- needed[!file.exists(file.path(root, needed))]
  if (length(absent)) {
    stop(
      "Some files of the linear SCAR module are missing:\n",
      paste(file.path(root, absent), collapse = "\n"),
      call. = FALSE
    )
  }
  archive <- opt("archive", ARCHIVE_DIR)

  cat("\nUsing linear SCAR project:\n", root, "\n", sep = "")
  wf <- new.env(parent = parent.env(.GlobalEnv))
  sys.source(
    file.path(root, "simulation_linear_SCAR_all18/R/workflow.R"),
    envir = wf
  )
  wf$run_linear_SCAR_all18(
    root,
    mode = "quick",
    action = action,
    out = opt("out", OUTPUT_DIR),
    archive = archive,
    allow_long_MCMC = FALSE
  )
})
