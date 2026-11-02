# Reproducible environment with renv. Run from revision/.
#   Rscript R/setup_env.R                            restore renv.lock
#   Rscript --no-init-file R/setup_env.R --from-library
#       write renv.lock from the library that produced the results (R_LIBS),
#       without installing anything

repos <- c(CRAN = "https://cloud.r-project.org")
options(repos = repos)

source(file.path("R", "registry.R"))
pkgs <- all_packages()

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

if ("--from-library" %in% commandArgs(TRUE)) {
  renv::snapshot(project = getwd(), library = .libPaths(), lockfile = "renv.lock",
                 packages = c("renv", pkgs), repos = repos, prompt = FALSE, force = TRUE)
  message("renv.lock written: ", length(pkgs), " packages and their dependencies.")
} else {
  renv::restore(project = getwd(), prompt = FALSE)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))
  message("Environment ready: ", length(pkgs), " packages.")
}
