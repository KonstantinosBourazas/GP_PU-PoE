## Runs the quick script twice in new R sessions and compares the results.
## This fits all 18 methods twice.
local({
  root <- getwd()
  for (i in 1:20) {
    if (
      file.exists(file.path(root, "simulation_linear_SCAR_all18/R/workflow.R"))
    ) {
      break
    }
    up <- dirname(root)
    if (identical(up, root)) {
      stop("Open GP-PU-PoE.Rproj before running this test.")
    }
    root <- up
  }
  wf <- new.env(parent = parent.env(.GlobalEnv))
  sys.source(
    file.path(root, "simulation_linear_SCAR_all18/R/workflow.R"),
    envir = wf
  )
  e <- new.env(parent = parent.env(.GlobalEnv))
  wf$lscar_load(root, e)

  e$lscar_preflight(attach = TRUE)
  rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) {
    rscript <- file.path(
      R.home("bin"),
      if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    )
  }
  e$lscar_assert(file.exists(rscript), "Rscript was not found.")
  testdir <- tempfile("lscar_repro_", tmpdir = file.path(root, "results"))
  dir.create(testdir, recursive = TRUE)
  entry <- file.path(root, "scripts/07_linear_SCAR_quick.R")
  modfiles <- list.files(
    file.path(root, "simulation_linear_SCAR_all18"),
    recursive = TRUE,
    full.names = TRUE
  )
  before <- e$lscar_md5(modfiles)
  for (k in 1:2) {
    output <- file.path(testdir, paste0("run", k))
    log <- file.path(testdir, paste0("run", k, ".log"))
    status <- system2(
      rscript,
      c("--vanilla", shQuote(entry), shQuote(paste0("--out=", output))),
      stdout = log,
      stderr = log
    )
    if (status != 0L) {
      stop(
        "A quick run failed. See ",
        log,
        "\n",
        paste(tail(readLines(log, warn = FALSE), 30), collapse = "\n")
      )
    }
  }
  read <- function(k, file) {
    utils::read.csv(
      file.path(testdir, paste0("run", k), "tables", file),
      stringsAsFactors = FALSE
    )
  }
  files <- c(
    "metrics_by_replication_all18_QUICK.csv",
    "node_scores_all18_QUICK.csv"
  )
  for (file in files) {
    a <- read(1, file)
    b <- read(2, file)
    keep <- names(a)[
      !grepl("runtime|elapsed|time_", names(a), ignore.case = TRUE)
    ]
    e$lscar_assert(
      isTRUE(all.equal(
        a[keep],
        b[keep],
        tolerance = 1e-10,
        check.attributes = TRUE
      )),
      paste("The two runs differ:", file)
    )
  }
  e$lscar_assert(
    identical(before, e$lscar_md5(modfiles)),
    "Module files were modified during fitting."
  )
  cat(
    "PASS: two quick runs in new R sessions agree within 1e-10.\nReports: ",
    testdir,
    "\n",
    sep = ""
  )
})
