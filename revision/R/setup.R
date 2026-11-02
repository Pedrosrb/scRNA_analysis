# ---- Shared bootstrap, sourced by run.R and by every module. ----

amrev_root <- function() {
  r <- Sys.getenv("AMREV_ROOT")
  if (!nzchar(r))
    stop("AMREV_ROOT is not set. Run through run.R, run the module with ",
         "Rscript, or source it from the project root.", call. = FALSE)
  normalizePath(r, mustWork = TRUE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

source(file.path(amrev_root(), "R", "registry.R"))
for (.f in sort(list.files(file.path(amrev_root(), "R", "lib"), "\\.R$", full.names = TRUE)))
  source(.f)
rm(.f)

load_params <- function(root = amrev_root()) {
  f <- Sys.getenv("AMREV_PARAMS", unset = file.path(root, "config", "params.R"))
  if (!file.exists(f)) stop("parameter file not found: ", f, call. = FALSE)
  e <- new.env(parent = baseenv())
  sys.source(f, envir = e)
  if (!exists("PARAMS", envir = e, inherits = FALSE))
    stop(f, " does not define PARAMS", call. = FALSE)
  structure(e$PARAMS, params_file = normalizePath(f))
}

# ---- inputs ----
check_inputs <- function(inputs, optional = character(0), hash = FALSE) {
  if (!length(inputs))
    return(data.frame(input = character(0), path = character(0),
                      required = logical(0), exists = logical(0)))
  flat <- unlist(inputs)
  info <- file.info(flat)
  out <- data.frame(
    input    = names(flat),
    path     = unname(flat),
    required = !names(flat) %in% optional,
    exists   = file.exists(flat),
    size_mb  = round(info$size / 2^20, 2),
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE)
  if (hash)
    out$md5 <- ifelse(out$exists & !info$isdir,
                      unname(tools::md5sum(flat)), NA_character_)
  out
}

print_inputs <- function(chk, id) {
  cat(sprintf("\n[preflight] %s\n", id))
  for (i in seq_len(nrow(chk))) {
    tag <- if (chk$exists[i]) "ok     " else if (chk$required[i]) "MISSING" else "absent "
    cat(sprintf("  %s %-18s %s\n", tag, chk$input[i], chk$path[i]))
  }
  invisible(chk)
}

# ---- provenance ----
git_state <- function(root) {
  run <- function(...) tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(root), ...),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  commit <- run("rev-parse", "--short", "HEAD")
  if (!length(commit) || !is.null(attr(commit, "status")))
    return(c(commit = NA_character_, dirty = NA_character_))
  c(commit = commit[1], dirty = as.character(length(run("status", "--porcelain")) > 0))
}

write_run_info <- function(M, status) {
  g <- git_state(M$root)
  info <- data.frame(
    module      = M$key,
    module_id   = M$id,
    questions   = paste(MODULES[[M$key]]$questions, collapse = ","),
    status      = status,
    started     = format(M$started, "%Y-%m-%d %H:%M:%S"),
    finished    = if (status == "completed") format(Sys.time(), "%Y-%m-%d %H:%M:%S") else NA,
    minutes     = if (status == "completed")
      round(as.numeric(difftime(Sys.time(), M$started, units = "mins")), 2) else NA,
    params_file = M$params_file,
    git_commit  = g[["commit"]],
    git_dirty   = g[["dirty"]],
    r_version   = R.version.string,
    host        = Sys.info()[["nodename"]],
    user        = Sys.info()[["user"]],
    stringsAsFactors = FALSE)
  utils::write.table(info, file.path(M$provenance, "run_info.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- module lifecycle ----
module_init <- function(key) {
  if (!key %in% names(MODULES)) stop("unknown module: ", key, call. = FALSE)
  root <- amrev_root()
  P    <- load_params(root)
  mp   <- P$modules[[key]]
  if (is.null(mp)) stop("no block for '", key, "' in PARAMS$modules", call. = FALSE)

  id  <- module_id(key)
  dir <- file.path(P$roots$results, id)
  M <- list(key = key, id = id, root = root, dir = dir,
            tables     = file.path(dir, "tables"),
            figures    = file.path(dir, "figures"),
            reports    = file.path(dir, "reports"),
            provenance = file.path(dir, "provenance"),
            logs       = file.path(dir, "logs"),
            params = mp, shared = P, inputs = mp$inputs,
            params_file = attr(P, "params_file"),
            started = Sys.time())

  chk <- check_inputs(mp$inputs, mp$optional, hash = isTRUE(P$run$hash_inputs))
  print_inputs(chk, id)
  missing <- chk$input[chk$required & !chk$exists]

  pkgs <- MODULES[[key]]$packages
  no_pkg <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(no_pkg)) {
    cat(sprintf("  MISSING package(s): %s\n", paste(no_pkg, collapse = ", ")))
    missing <- c(missing, paste0("package:", no_pkg))
  }

  if (Sys.getenv("AMREV_PREFLIGHT") == "1")
    quit(save = "no", status = if (length(missing)) 1L else 0L)
  if (length(missing))
    stop("required input(s) or package(s) missing: ", paste(missing, collapse = ", "),
         call. = FALSE)

  if (Sys.getenv("AMREV_CLEAN") == "1")
    unlink(c(M$tables, M$figures, M$reports), recursive = TRUE)
  for (d in c(M$tables, M$figures, M$reports, M$provenance, M$logs))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  utils::write.table(chk, file.path(M$provenance, "inputs.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  shared <- P[setdiff(names(P), "modules")]
  writeLines(c(sprintf("# Resolved parameters for %s, %s", id, format(M$started)),
               sprintf("# Source: %s", M$params_file),
               "module <- ", deparse(mp), "", "shared <- ", deparse(shared)),
             file.path(M$provenance, "params_used.R"))
  write_run_info(M, "running")

  flow_reset()
  set.seed(mp$seed %||% 1L)
  cat(sprintf("[module] %s | seed %s | output %s\n\n", id, mp$seed %||% 1L, dir))
  M
}

module_finish <- function(M) {
  fl <- flow_table()
  if (!is.null(fl))
    utils::write.table(fl, file.path(M$tables, "sample_flow.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)

  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(M$provenance, "sessionInfo.txt"))
  write_run_info(M, "completed")
  cat(sprintf("\n[module] %s completed in %.1f min -> %s\n", M$id,
              as.numeric(difftime(Sys.time(), M$started, units = "mins")), M$dir))
  invisible(M)
}

