#!/usr/bin/env Rscript
# ---- Runner for the reviewer-response modules ----

args    <- commandArgs(trailingOnly = TRUE)
flags   <- args[startsWith(args, "--")]
targets <- setdiff(args, flags)

jobs_flag <- grep("^--jobs=[0-9]+$", flags, value = TRUE)
jobs  <- if (length(jobs_flag)) as.integer(sub("^--jobs=", "", jobs_flag[1])) else 1L
flags <- setdiff(flags, jobs_flag)
known_flags <- c("--preflight", "--clean", "--stop-on-error")
if (length(bad <- setdiff(flags, known_flags)))
  stop("unknown flag(s): ", paste(bad, collapse = ", "), call. = FALSE)

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- normalizePath(if (length(self)) dirname(self[1]) else getwd())
Sys.setenv(AMREV_ROOT = root)
setwd(root)
source(file.path(root, "R", "setup.R"))
P <- load_params(root)

# ---- listing ----
if (!length(targets)) {
  cat("Modules (key | questions | output folder):\n")
  for (k in names(MODULES))
    cat(sprintf("  %-13s %-12s %s\n", k,
                paste(MODULES[[k]]$questions, collapse = ","), module_id(k)))
  cat("\nResults root: ", P$roots$results, "\n",
      "Usage: Rscript run.R [--preflight] [--clean] [--stop-on-error] [--jobs=N] <all | key | tag> ...\n",
      sep = "")
  quit(save = "no", status = 0)
}

# ---- resolve targets ----
resolve <- function(t) {
  if (t == "all") return(names(MODULES))
  if (t %in% names(MODULES)) return(t)
  by_tag <- names(MODULES)[vapply(MODULES, function(m) t %in% m$questions, logical(1))]
  if (length(by_tag)) return(by_tag)
  by_id <- names(MODULES)[vapply(names(MODULES), module_id, "") == sub("\\.R$", "", t)]
  if (length(by_id)) return(by_id)
  stop("no module matches '", t, "'. Run without arguments to list modules.",
       call. = FALSE)
}
keys <- unique(unlist(lapply(targets, resolve)))

preflight <- "--preflight" %in% flags
rscript   <- file.path(R.home("bin"), "Rscript")
env <- c(paste0("AMREV_ROOT=", root),
         paste0("AMREV_PREFLIGHT=", as.integer(preflight)),
         paste0("AMREV_CLEAN=", as.integer("--clean" %in% flags)))

# ---- run ----
history <- file.path(P$roots$results, "run_history.tsv")

run_one <- function(k) {
  script <- file.path(root, "modules", MODULES[[k]]$script)
  t0 <- Sys.time()
  if (preflight)
    return(data.frame(module = k, exit_status = system2(rscript, shQuote(script), env = env)))

  log_dir <- file.path(P$roots$results, module_id(k), "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(log_dir, paste0(format(t0, "%Y%m%d_%H%M%S"), ".log"))
  cat(sprintf("[run] %-13s -> %s\n", k, log))

  st   <- system2(rscript, shQuote(script), stdout = log, stderr = log, env = env)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[end] %-13s %s (%.1f min)\n", k, if (st == 0) "ok" else "FAILED", mins))

  info_f <- file.path(P$roots$results, module_id(k), "provenance", "run_info.tsv")
  if (st != 0 && file.exists(info_f)) {
    info <- utils::read.delim(info_f, colClasses = "character")
    info$status <- "failed"
    utils::write.table(info, info_f, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  data.frame(started = format(t0, "%Y-%m-%d %H:%M:%S"), module = k,
             module_id = module_id(k), exit_status = st, minutes = round(mins, 2),
             log = log, git_commit = git_state(root)[["commit"]], jobs = jobs)
}

if (jobs > 1 && !preflight && .Platform$OS.type == "unix") {
  rows <- parallel::mclapply(keys, run_one, mc.cores = jobs, mc.preschedule = FALSE)
} else {
  rows <- list()
  for (k in keys) {
    rows[[k]] <- run_one(k)
    if (rows[[k]]$exit_status != 0 && "--stop-on-error" %in% flags) break
  }
}
res <- do.call(rbind, rows)
status <- setNames(res$exit_status, res$module)

if (!preflight) {
  dir.create(dirname(history), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(res, history, sep = "\t", quote = FALSE, row.names = FALSE,
                     append = file.exists(history), col.names = !file.exists(history))
}

# ---- summary ----
failed <- names(status)[status != 0]
cat(sprintf("\n%s: %d of %d module(s) run, %d failed%s\n",
            if (preflight) "Preflight" else "Run", length(status), length(keys), length(failed),
            if (length(failed)) paste0(" (", paste(failed, collapse = ", "), ")") else ""))
quit(save = "no", status = as.integer(length(failed) > 0))
