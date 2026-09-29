## Installs the R packages used in this repository. Packages that are already
## installed are left unchanged, so no installed version is upgraded.
##
## The core packages are enough for the illustrative example, the prior
## sensitivity analysis and the tables of the simulation studies, which are
## rebuilt from the archived results. The refit packages are needed only to
## refit the simulation methods (scripts 03, 05 and 07, and the optional 04b,
## 06b and 08b). INLA is not on CRAN and is installed from its own repository.

INSTALL_REFIT <- FALSE # TRUE also installs the packages needed for refitting

core <- c("Matrix", "igraph", "RSpectra", "coda", "ggplot2", "patchwork")
refit <- c("e1071", "AdaSampling", "BayesLogit")
cran <- "https://cloud.r-project.org"
inla <- "https://inla.r-inla-download.org/R/stable"

absent <- function(p) {
  p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
}

wanted <- if (INSTALL_REFIT) c(core, refit) else core
todo <- absent(wanted)

if (length(todo)) {
  utils::install.packages(todo, repos = cran)
}

if (INSTALL_REFIT && length(absent("INLA"))) {
  utils::install.packages(
    "INLA",
    repos = c(CRAN = cran, INLA = inla),
    dependencies = TRUE
  )
}

needed <- if (INSTALL_REFIT) c(wanted, "INLA") else wanted
left <- absent(needed)

if (length(left)) {
  stop("Installation incomplete: ", paste(left, collapse = ", "), call. = FALSE)
}

if (utils::packageVersion("ggplot2") < "3.4.0") {
  stop(
    "ggplot2 3.4.0 or later is required. Update it before running the scripts.",
    call. = FALSE
  )
}

print(data.frame(
  package = needed,
  version = vapply(
    needed,
    function(p) as.character(utils::packageVersion(p)),
    character(1)
  ),
  row.names = NULL
))
