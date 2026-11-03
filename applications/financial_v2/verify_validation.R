# Independent duplicated-sample rank checks on the saved end-to-end audits.
args <- commandArgs(TRUE)
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
folder <- if (length(args)) args[1] else file.path(root, "validation_v1_d5")
source(file.path(root, "../../R/learners.R"))
source(file.path(root, "../../R/tests.R"))
source(file.path(root, "panel_calibration.R"))
dyn.load(file.path(root, "native/panel_auc.dll"))
files <- list.files(file.path(folder, "audits"), pattern = "^.*_r00001\\.rds$", full.names = TRUE)
rows <- list()
for (path in files) {
  audit <- readRDS(path)
  for (len in names(audit$retained)) {
    retained <- audit$retained[[len]]
    columns <- unique(c(1L, 2L, 5L, 17L, ncol(retained$counts)))
    for (q in seq_along(audit$questions)) {
      question <- audit$questions[[q]]
      a <- retained$calibrated[[q]]$auc_matrix[, columns, drop = FALSE]
      direct <- matrix(NA_real_, nrow(a), ncol(a))
      for (b in seq_along(columns)) {
        index <- rep(seq_len(nrow(question$S)), retained$counts[question$time_id, columns[b]])
        for (c in seq_len(ncol(question$S)))
          direct[c, b] <- auc_vec(question$S[index, c], question$groups[index] == question$targets[c])
      }
      error <- max(abs(a - direct))
      stopifnot(error < 2e-13)
      rows[[length(rows) + 1L]] <- data.frame(scenario = audit$scenario, block_length = as.integer(len),
                                            question = q, draws_checked = length(columns), max_auc_error = error)
    }
  }
}
slow_path <- file.path(root, "validation_smoke_v2/audits/end_to_end.rds")
fast_path <- file.path(root, "validation_smoke_v3/audits/end_to_end.rds")
if (file.exists(slow_path) && file.exists(fast_path)) {
  slow <- readRDS(slow_path); fast <- readRDS(fast_path)
  for (q in seq_along(slow$questions)) for (key in c("S", "groups", "targets", "time_id", "scale"))
    stopifnot(identical(slow$questions[[q]][[key]], fast$questions[[q]][[key]]))
  cat("Original and accelerated kNN: four complete fitted score matrices bit-identical.\n")
}
result <- do.call(rbind, rows)
write.csv(result, file.path(folder, "independent_auc_audit.csv"), row.names = FALSE)
cat(sprintf("Passed explicit duplicated-sample rank audit: %d scenarios, %d question/block cells, max AUC error %.3g\n",
            length(files), nrow(result), max(result$max_auc_error)))

circular_files <- list.files(file.path(folder,"audits"),pattern="^circular_.*_L[0-9]+\\.rds$",full.names=TRUE)
if (length(circular_files)) {
  circular_rows <- list()
  for (path in circular_files) {
    case <- sub("^circular_(.*)_L[0-9]+\\.rds$","\\1",basename(path))
    audit <- readRDS(file.path(folder,"audits",paste0(case,"_r00001.rds")))
    retained <- readRDS(path)
    columns <- unique(c(1L,2L,5L,17L,ncol(retained$counts)))
    stopifnot(all(colSums(retained$counts)==nrow(retained$counts)))
    for (q in seq_along(audit$questions)) {
      question <- audit$questions[[q]]
      direct <- matrix(NA_real_,length(question$targets),length(columns))
      for (b in seq_along(columns)) {
        index <- rep(seq_len(nrow(question$S)),retained$counts[question$time_id,columns[b]])
        for (c in seq_len(ncol(question$S)))
          direct[c,b] <- auc_vec(question$S[index,c],question$groups[index]==question$targets[c])
      }
      error <- max(abs(direct-retained$calibrated[[q]]$auc_matrix[,columns,drop=FALSE]))
      stopifnot(error<2e-13)
      circular_rows[[length(circular_rows)+1L]] <- data.frame(scenario=case,source=basename(path),question=q,max_auc_error=error)
    }
    previous <- file.path(root,"validation_v1_d5/audits",paste0(case,"_r00001.rds"))
    if (file.exists(previous)) {
      original <- readRDS(previous)
      for (q in seq_along(audit$questions)) stopifnot(identical(original$questions[[q]]$S,audit$questions[[q]]$S))
    }
  }
  circular <- do.call(rbind,circular_rows)
  write.csv(circular,file.path(folder,"independent_circular_auc_audit.csv"),row.names=FALSE)
  cat(sprintf("Circular-block explicit rank audit: %d question/block cells, max error %.3g; pilot score matrices match original seeds.\n",nrow(circular),max(circular$max_auc_error)))
}
