# Resumable end-to-end calibration validation, using the original learner library.
args <- commandArgs(TRUE)
arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit[1L]) else default
}
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
src <- normalizePath(file.path(root, "..", "..", "R"), winslash = "/")
dll <- file.path(root, "native", "panel_auc.dll")
mode <- arg("mode", "run")
registry <- c("null_iid", "null_ar50", "null_ar90", "null_garch95", "null_garch99", "null_t5_ar",
              "alt_sector_mean", "alt_sector_scale", "alt_sector_dependence", "alt_one_bank")
scenario <- strsplit(arg("scenarios", paste(registry[1:6], collapse = ",")), ",", fixed = TRUE)[[1L]]
if (!all(scenario %in% registry)) stop("Unknown validation scenario")
config <- list(d = as.integer(arg("d", "5")), B = as.integer(arg("B", if (mode == "smoke") "99" else "999")),
               lengths = as.integer(strsplit(arg("lengths", "4,8,12"), ",", fixed = TRUE)[[1L]]),
               seed0 = 270000000L, audit_reps = 1L,
               outdir = arg("out", file.path(root, if (mode == "smoke") "validation_smoke" else "validation_v1_d5")))
cores <- as.integer(arg("cores", "30")); reps <- as.integer(arg("reps", if (mode == "smoke") "1" else "100"))
chunk <- as.integer(arg("chunk", "25"))
dir.create(config$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$outdir, "chunks"), showWarnings = FALSE)
dir.create(file.path(config$outdir, "audits"), showWarnings = FALSE)
load_runtime <- function(root, src, dll) {
  for (pkg in c("glmnet", "ranger")) {
    available <- FALSE
    for (a in seq_len(3L)) if (isTRUE(tryCatch(requireNamespace(pkg, quietly = TRUE), error = function(e) FALSE))) { available <- TRUE; break }
    if (!available) stop("Package unavailable: ", pkg)
  }
  for (f in c("learners.R", "learners_lib.R", "tests.R")) sys.source(file.path(src, f), envir = .GlobalEnv)
  dyn.load(dll)
  sys.source(file.path(root, "panel_calibration.R"), envir = .GlobalEnv)
  install_panel_knn()
  invisible(TRUE)
}
load_runtime(root, src, dll)
calendar <- read.csv(file.path(root, "design_two_tests", "block_calendar.csv"))
calendar <- calendar[calendar$dimension == config$d, ]
mask_values <- tolower(as.character(calendar$complete_all_main_stocks))
stopifnot(all(mask_values %in% c("true", "false")))
mask_values <- mask_values == "true"
masks <- list(fit = mask_values[calendar$split == "fit"],
              test = mask_values[calendar$split == "test"])
stopifnot(length(masks$fit) > 0, all(config$lengths <= length(masks$test)))
fingerprint <- list(config = config[setdiff(names(config), "outdir")], scenarios = registry,
                    hashes = tools::md5sum(c(file.path(root, c("run_validation.R", "panel_calibration.R", "native/panel_auc.cpp")),
                                             dll, file.path(root, "design_two_tests", "block_calendar.csv"),
                                             file.path(src, c("learners.R", "learners_lib.R", "tests.R")))))
fingerprint_file <- file.path(config$outdir, "fingerprint.rds")
if (file.exists(fingerprint_file)) {
  if (!identical(readRDS(fingerprint_file), fingerprint)) stop("Validation fingerprint changed; use a new output directory")
} else saveRDS(fingerprint, fingerprint_file)
writeLines(capture.output(sessionInfo()), file.path(config$outdir, "sessionInfo.txt"))

if (mode == "smoke") {
  # Independent explicit resampling, with integer ties and missing calendar rows.
  set.seed(73)
  valid <- c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
  time <- rep(which(valid), 3L); labels <- rep(1:3, each = sum(valid)); targets <- rep(1:3, 2L)
  S <- matrix(sample(c(-1, 0, 1), length(time) * length(targets), replace = TRUE), ncol = length(targets))
  counts <- panel_counts(length(valid), 3L, 19L, 71L)
  weighted <- panel_auc(S, labels, targets, time, counts)
  direct <- matrix(NA_real_, nrow(weighted), ncol(weighted))
  for (b in seq_len(ncol(counts))) {
    ix <- rep(seq_len(nrow(S)), counts[time, b])
    for (c in seq_len(ncol(S))) direct[c, b] <- auc_vec(S[ix, c], labels[ix] == targets[c])
  }
  stopifnot(max(abs(weighted - direct)) < 1e-13)
  tie <- S * 0
  stopifnot(all(panel_auc(tie, labels, targets, time, counts) == .5))
  cat("Native weighted AUC matches explicit duplicated-sample ranks, including ties.\n")
  started <- proc.time()[["elapsed"]]
  result <- panel_rep(scenario[1L], 1L, config$d, config$B, config$lengths, masks, config$seed0, keep = TRUE)
  saveRDS(result$audit, file.path(config$outdir, "audits", "end_to_end.rds"))
  write.csv(result$rows, file.path(config$outdir, "smoke_rows.csv"), row.names = FALSE)
  print(result$rows[, c("question", "block_length", "p", "statistic", "learner_warning_count")], row.names = FALSE)
  cat(sprintf("End-to-end smoke completed in %.2f seconds\n", proc.time()[["elapsed"]] - started))
  quit(status = 0)
}

jobs <- list()
for (case in scenario) for (first in seq.int(1L, reps, by = chunk)) {
  job <- list(scenario = case, start = first, end = min(first + chunk - 1L, reps))
  path <- file.path(config$outdir, "chunks", sprintf("%s__r%05d-%05d.csv", case, job$start, job$end))
  if (!file.exists(path)) jobs[[length(jobs) + 1L]] <- job
}
cat(sprintf("Validation d=%d | %d scenarios | %d reps each | B=%d | %d pending chunks | cores=%d\n",
            config$d, length(scenario), reps, config$B, length(jobs), cores)); flush.console()
started <- proc.time()[["elapsed"]]
if (length(jobs)) {
  if (cores > 1L) {
    cl <- parallel::makePSOCKcluster(min(cores, length(jobs)), outfile = file.path(config$outdir, "workers.log"))
    parallel::clusterCall(cl, load_runtime, root = root, src = src, dll = dll)
  } else cl <- NULL
  batches <- split(seq_along(jobs), ceiling(seq_along(jobs) / cores))
  done <- 0L
  tryCatch(for (batch in batches) {
    if (is.null(cl)) lapply(jobs[batch], panel_run_job, config = config, masks = masks) else
      parallel::parLapplyLB(cl, jobs[batch], function(job, config, masks) panel_run_job(job, config, masks), config, masks)
    done <- done + length(batch)
    cat(sprintf("Progress %d/%d chunks | elapsed %.1f minutes\n", done, length(jobs), (proc.time()[["elapsed"]] - started) / 60)); flush.console()
  }, finally = {if (!is.null(cl)) parallel::stopCluster(cl)})
}
files <- list.files(file.path(config$outdir, "chunks"), pattern = "\\.csv$", full.names = TRUE)
combined <- do.call(rbind, lapply(files, read.csv))
write.csv(combined, file.path(config$outdir, "all_replications.csv"), row.names = FALSE)
cat(sprintf("Validation run complete: %d recorded rows, %.1f minutes in this process\n", nrow(combined), (proc.time()[["elapsed"]] - started) / 60))
