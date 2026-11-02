# Reproducible environment with renv. Run from the repository root.
#   Rscript src/00_setup_env.R                            restore renv.lock
#   Rscript --no-init-file src/00_setup_env.R --from-library
#       write renv.lock from the library that produced the results (R_LIBS),
#       without installing anything

repos <- c(CRAN = "https://cloud.r-project.org")
options(repos = repos)

e <- new.env()
sys.source(file.path("src", "03_functions.R"), envir = e)
pkgs <- e$required_packages

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
