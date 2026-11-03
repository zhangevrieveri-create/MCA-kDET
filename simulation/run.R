# =============================================================================
# MCA-kDET simulation study -- command-line entry point
#
#   Rscript run.R --exp=smoke --cores=4
#   Rscript run.R --exp=E2 --profile=quick --cores=30 [--reps=..] [--B=..] [--filter=E2_A_]
#
# E1-E7 (revised design) need --profile=full (precise version, for the paper)
# or --profile=quick (fewer replications). Their default output directory is
# results_v2_full/ or results_v2_quick/; other experiments write to results/.
#
# Results:  <out>/<exp>/<cell>__r<from>-<to>.csv   (one file per chunk)
# Log:      <out>/<exp>.log
# Re-running the same command resumes; finished chunks are skipped.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1L]) else default
}
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
SRC_DIR <- file.path(script_dir, "R")
setwd(script_dir)

# Packages are never installed from here: on this machine a failed load is
# usually a sporadic lazy-load/decompression error, not a missing package, and
# an automatic reinstall once half-overwrote glmnet (21 Oct 2026). Install new
# packages separately (CLAUDE.md); here a failed load is retried, then stops.
load_pkg <- function(pkg, tries = 3L) {
  for (a in seq_len(tries))
    if (isTRUE(tryCatch(requireNamespace(pkg, quietly = TRUE), error = function(e) FALSE)))
      return(invisible(TRUE))
  stop("package '", pkg, "' could not be loaded after ", tries,
       " tries; if it is installed, simply run the same command again")
}
load_pkg("glmnet")

PROFILE <- get_arg("profile", NA_character_)   # read by experiments.R
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "garch_ext.R",
            "competitors2.R", "rotation.R", "learners_lib2.R", "models_s.R", "experiments.R"))
  source(file.path(SRC_DIR, f))

exp_name <- get_arg("exp", "smoke")
cores <- as.integer(get_arg("cores", max(1L, parallel::detectCores(logical = FALSE) - 1L)))
reps_over <- get_arg("reps"); B_over <- get_arg("B"); filt <- get_arg("filter")
chunk <- as.integer(get_arg("chunk", 50L))
is_v2 <- exp_name %in% setdiff(V2_EXPS, "calib2")   # E1-E8
if (!is.na(PROFILE) && !PROFILE %in% c("full", "quick")) stop("--profile must be full or quick")
# new library / new scenarios (third review round, task of 1 Nov 2026): never mixed
# with results_v2_*; the pilot has no profile
is_v3 <- exp_name %in% c("E11", "E12", "E13")
outdir <- get_arg("out", if (is_v2) paste0("results_v2_", PROFILE)
                         else if (is_v3) paste0("results_v3_", PROFILE)
                         else if (exp_name == "libpilot3") "results_v3_pilot" else "results")

grid <- make_grid(exp_name)
if (!is.null(filt)) grid <- grid[grepl(filt, grid$cell), , drop = FALSE]
if (!is.null(reps_over)) grid$reps <- as.integer(reps_over)
if (!is.null(B_over)) grid$B <- as.integer(B_over)
# validate the libraries before starting the cluster; ranger is never
# installed automatically (install it separately first, see CLAUDE.md)
for (r in seq_len(nrow(grid))) parse_libsets(grid$libsets[r], parse_lib(grid$lib[r]))
if (any(grepl("rf", grid$lib))) load_pkg("ranger")
if (any(grepl("songchen", grid$groups))) load_pkg("gTestsMulti")

dir.create(outdir, showWarnings = FALSE)
logf <- file.path(outdir, paste0(exp_name, ".log"))
log_message <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = logf, append = TRUE)
}
log_message(sprintf("=== %s%s | %d cells | %d reps | cores = %d | R %s | glmnet %s | started %s",
            exp_name, if (is.na(PROFILE)) "" else paste0(" [profile ", PROFILE, "]"),
            nrow(grid), sum(grid$reps), cores, getRversion(),
            as.character(utils::packageVersion("glmnet")), format(Sys.time())))
write.csv(grid, file.path(outdir, paste0(exp_name, "_grid.csv")), row.names = FALSE)

log_message(sprintf("Starting %d worker(s); loading R code and packages...", cores))
flush.console()
cl <- if (cores > 1L) make_cluster(cores, SRC_DIR) else NULL
log_message("Workers ready; starting cells.")
t0 <- proc.time()[["elapsed"]]

# overall progress bar, printed after every chunk; the remaining time is
# extrapolated from the chunks run in this session (chunks found on disk are
# counted as done but not used for the rate)
total_reps <- sum(grid$reps); done_reps <- 0; new_reps <- 0; cur_cell <- 0L
progress <- function(n, new) {
  done_reps <<- done_reps + n
  if (new) new_reps <<- new_reps + n
  el <- (proc.time()[["elapsed"]] - t0) / 60
  frac <- done_reps / total_reps
  eta <- if (new_reps > 0) el / new_reps * (total_reps - done_reps) else NA
  w <- 30L; nb <- round(w * frac)
  log_message(sprintf("  [%s%s] %3.0f%% | %d/%d reps | cell %d/%d | elapsed %.1f min | remaining ~%s min",
              strrep("#", nb), strrep("-", w - nb), 100 * frac, as.integer(done_reps),
              as.integer(total_reps), cur_cell, nrow(grid), el,
              if (is.na(eta)) "?" else sprintf("%.1f", eta)))
  flush.console()
}
for (r in seq_len(nrow(grid))) {
  cell <- as.list(grid[r, ])
  cur_cell <- r
  log_message(sprintf("[%s] START cell %d/%d: %s | %d reps | B=%d",
              format(Sys.time(), "%H:%M:%S"), r, nrow(grid), cell$cell, cell$reps, cell$B))
  flush.console()
  run_cell(cell, cell$reps, cl = cl, outdir = outdir, chunk = min(chunk, cell$reps), log = log_message,
           progress = progress)
}
if (!is.null(cl)) parallel::stopCluster(cl)
log_message(sprintf("=== %s finished in %.1f min", exp_name, (proc.time()[["elapsed"]] - t0) / 60))

source(file.path(SRC_DIR, "summarize.R"))
summ <- summarize_exp(exp_name, outdir)
if (!identical(get_arg("print-summary", "yes"), "no"))
  log_message(paste(capture.output(print(summ, row.names = FALSE, digits = 3)), collapse = "\n"))
