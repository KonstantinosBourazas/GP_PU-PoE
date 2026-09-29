## RStudio: open GP-PU-PoE.Rproj, then Source this file.
## Command line options override the settings below (see --help).
##
## MODE "saved": draws the figure and the tables of the paper from the saved
##               full fits in archive/illustration/, in a few seconds, no MCMC.
## MODE "full":  runs the MCMC of the four models again, with the settings
##               below. With the values of the paper it takes about 45 minutes
##               and reproduces the saved fits exactly.
## MODE "quick": short MCMC runs that only check that the code works.
MODE <- "saved"
ACTION <- "all" # "check", "prepare", "fit", "plot", "diagnostics", "all"

## Settings of a new run (MODE "full"). The values of the paper are shown.
## Other values write to a separate folder under results/illustration/.
MCMC_SEED <- 53L # another value gives new chains for the same data
ITERATIONS <- 150000L # sampling iterations after the warm up
THIN <- 30L # ITERATIONS / THIN draws are kept

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
        if (file.exists(file.path(path, ".here"))) {
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
      "Rscript --vanilla scripts/01_illustrative_example.R --mode=saved\n",
      "Options: --mode=saved|quick|full; --action=check|prepare|fit|plot|diagnostics|all;\n",
      "         --mcmc-seed=N; --iterations=N; --thin=N; --out=PATH; --fresh\n",
      "Saved draws the results of the paper from archive/illustration/ without MCMC.\n",
      "Quick only checks that the code runs. Full runs the MCMC with the settings\n",
      "of the paper, unless --mcmc-seed, --iterations or --thin change them.\n",
      sep = ""
    )
  } else {
    unknown <- args[
      !grepl("^--(mode|action|out|mcmc-seed|iterations|thin)=", args) &
        args != "--fresh"
    ]
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
    sys.source(file.path(root, "R", "illustration", "workflow.R"), envir = wf)
    number <- function(name, default) {
      value <- option_value(name, NULL)
      if (is.null(value)) {
        return(default)
      }
      x <- suppressWarnings(as.numeric(value))
      if (length(x) != 1L || is.na(x)) {
        stop("--", name, " must be a number.")
      }
      x
    }
    wf$run_illustration(
      root = root,
      mode = option_value("mode", MODE),
      action = option_value("action", ACTION),
      out_dir = option_value("out", NULL),
      fresh = "--fresh" %in% args,
      mcmc_seed = number("mcmc-seed", MCMC_SEED),
      iterations = number("iterations", ITERATIONS),
      thin = number("thin", THIN)
    )
  }
})
