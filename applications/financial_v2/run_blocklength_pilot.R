# Null-only diagnostic pilot after the predeclared short-block candidate failed.
# Reuses exactly the first pilot's simulated panels/seeds; holds out later seeds.
args <- commandArgs(TRUE)
arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit[1L]) else default
}
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
src <- normalizePath(file.path(root, "../../R"), winslash = "/")
dll <- file.path(root, "native/panel_auc.dll")
config <- list(d = 5L, B = 999L, lengths = c(8L, 20L, 40L, 60L, 80L), seed0 = 270000000L,
               outdir = arg("out", file.path(root, "validation_longblock_pilot_d5")))
reps <- as.integer(arg("reps", "100")); cores <- as.integer(arg("cores", "30"))
first_rep <- as.integer(arg("first", "1"))
scenarios <- c("null_iid", "null_ar90", "null_garch99")
load_runtime <- function(root, src, dll) {
  for (pkg in c("glmnet", "ranger")) {
    okay <- FALSE
    for (a in 1:3) if (isTRUE(tryCatch(requireNamespace(pkg, quietly = TRUE), error = function(e) FALSE))) {okay <- TRUE; break}
    if (!okay) stop("Package unavailable: ", pkg)
  }
  for (f in c("learners.R", "learners_lib.R", "tests.R")) sys.source(file.path(src, f), envir = .GlobalEnv)
  dyn.load(dll)
  sys.source(file.path(root, "panel_calibration.R"), envir = .GlobalEnv)
  install_panel_knn()
  invisible(TRUE)
}
load_runtime(root, src, dll)
cal <- read.csv(file.path(root, "design_two_tests/block_calendar.csv"))
cal <- cal[cal$dimension == 5L, ]
valid <- tolower(as.character(cal$complete_all_main_stocks)) == "true"
masks <- list(fit = valid[cal$split == "fit"], test = valid[cal$split == "test"])
dir.create(config$outdir, recursive = TRUE, showWarnings = FALSE)
for (name in c("chunks", "audits", "score_cache")) dir.create(file.path(config$outdir, name), showWarnings = FALSE)
dir.create(file.path(config$outdir, "circular/chunks"), recursive = TRUE, showWarnings = FALSE)
fingerprint <- list(config = config, scenarios = scenarios,
                    hashes = tools::md5sum(c(file.path(root, c("run_blocklength_pilot.R", "panel_calibration.R")),
                                             dll, file.path(root, "design_two_tests/block_calendar.csv"),
                                             file.path(src, c("learners.R", "learners_lib.R", "tests.R")))))
fp <- file.path(config$outdir, "fingerprint.rds")
if (file.exists(fp)) stopifnot(identical(fingerprint, readRDS(fp))) else saveRDS(fingerprint, fp)
writeLines(capture.output(sessionInfo()), file.path(config$outdir, "sessionInfo.txt"))
panel_counts_circular <- function(Tn, length, B, seed) {
  set.seed(seed)
  out <- matrix(0L, Tn, B+1L); out[,1L] <- 1L
  for (b in seq_len(B)) {
    starts <- sample.int(Tn, ceiling(Tn/length), replace = TRUE) - 1L
    index <- (as.vector(t(outer(starts, seq.int(0L,length-1L), "+")))[seq_len(Tn)] %% Tn) + 1L
    out[,b+1L] <- tabulate(index, Tn)
  }
  out
}
# Mechanical circular-index check, separate from statistical calibration.
check <- panel_counts_circular(23L, 8L, 31L, 918L)
stopifnot(all(colSums(check) == 23L), all(check >= 0), all(check[,1L] == 1L))
block_job <- function(job, config, masks) {
  dest <- file.path(config$outdir, "chunks", sprintf("%s__r%05d-%05d.csv", job$scenario, job$start, job$end))
  circular_dest <- file.path(config$outdir, "circular/chunks", basename(dest))
  if (file.exists(dest) && file.exists(circular_dest)) return(invisible(TRUE))
  rows <- circular_rows <- list()
  for (r in seq.int(job$start, job$end)) {
    value <- NULL
    for (attempt in 1:3) {
      value <- tryCatch(panel_rep(job$scenario, r, config$d, config$B, config$lengths, masks, config$seed0, keep = TRUE), error = function(e) e)
      if (!inherits(value, "error")) break
    }
    if (inherits(value, "error")) stop(job$scenario, " r=", r, ": ", conditionMessage(value))
    # Retain all fitted scores for further calibration diagnostics without refits.
    path <- file.path(config$outdir, "score_cache", sprintf("%s_r%05d.rds", job$scenario, r))
    if (!file.exists(path)) saveRDS(list(questions = value$audit$questions, scenario = job$scenario, rep = r, seed = value$audit$seed), path)
    if (r == 1L) {
      ap <- file.path(config$outdir, "audits", sprintf("%s_r%05d.rds", job$scenario, r))
      if (!file.exists(ap)) saveRDS(value$audit, ap)
    }
    rows[[length(rows)+1L]] <- value$rows
    circular <- value$rows
    for (len in config$lengths) {
      counts <- panel_counts_circular(length(masks$test), len, config$B, value$audit$seed + 50L + len)
      calibrated <- lapply(value$audit$questions, panel_calibrate, counts = counts)
      ids <- which(circular$block_length == len)
      circular$p[ids] <- vapply(calibrated, `[[`, numeric(1), "p")
      circular$p_holm[ids] <- p.adjust(circular$p[ids], "holm")
      circular$critical95[ids] <- vapply(calibrated, `[[`, numeric(1), "critical95")
      if (r == 1L) {
        ap <- file.path(config$outdir, "audits", sprintf("circular_%s_L%d.rds", job$scenario, len))
        if (!file.exists(ap)) saveRDS(list(counts = counts, calibrated = calibrated), ap)
      }
    }
    circular_rows[[length(circular_rows)+1L]] <- circular
  }
  for (target in c("moving", "circular")) {
    path <- if (target == "moving") dest else circular_dest
    tab <- do.call(rbind, if (target == "moving") rows else circular_rows)
    if (file.exists(path)) {
      old <- read.csv(path)
      cols <- setdiff(names(tab), "elapsed_seconds")
      stopifnot(isTRUE(all.equal(old[,cols],tab[,cols],check.attributes=FALSE,tolerance=1e-12)))
    } else {
      write.csv(tab,paste0(path,".part"),row.names=FALSE)
      stopifnot(file.rename(paste0(path,".part"),path))
    }
  }
  cat(sprintf("Completed long-block pilot %s reps %d-%d\n", job$scenario, job$start, job$end)); flush.console()
  invisible(TRUE)
}
jobs <- list()
for (case in scenarios) for (first in seq.int(first_rep, first_rep+reps-1L, by = 10L)) {
  job <- list(scenario = case, start = first, end = min(first+9L, first_rep+reps-1L))
  path <- file.path(config$outdir, "chunks", sprintf("%s__r%05d-%05d.csv", case, job$start, job$end))
  if (!file.exists(path) || !file.exists(file.path(config$outdir,"circular/chunks",basename(path)))) jobs[[length(jobs)+1L]] <- job
}
cat(sprintf("Null-only long-block diagnostic pilot: %d pending chunks, %d cores, B=999, L=8/20/40/60/80.\n", length(jobs), cores)); flush.console()
started <- proc.time()[["elapsed"]]
if (length(jobs)) {
  cl <- parallel::makePSOCKcluster(min(cores, length(jobs)), outfile = file.path(config$outdir, "workers.log"))
  parallel::clusterCall(cl, load_runtime, root, src, dll)
  parallel::clusterExport(cl, c("block_job", "panel_counts_circular"))
  tryCatch(invisible(parallel::parLapplyLB(cl, jobs, function(job, config, masks) block_job(job, config, masks), config, masks)),
           finally = parallel::stopCluster(cl))
}
cat(sprintf("Long-block pilot completed in %.1f minutes.\n", (proc.time()[["elapsed"]]-started)/60))
