# =============================================================================
# Rejection rates per cell and method, with Monte Carlo standard errors
# =============================================================================

P_COLS <- c("mca", "mca_l1", "mca_nn", "mca_2s", "mca_sum", "mca_comb",
            "shared_nn", "shared_lib", "macc", "pair", "pair_l1", "pair_nn",
            "disco", "mmd", "mst", "schott", "mca_ms_med", "mca_ms_mean",
            # third review round (1 Nov 2026): G2 accuracy, G1 competitors
            "macc_ovr", "maxmean", "clx_mean", "clx_cov", "schott_perm", "mmd_max",
            "kmd", "songchen", "songchen_A",
            # one-vs-rest max-type competitors (scenarios S3)
            "maxmean_ovr", "clx_cov_ovr")

read_exp <- function(exp, outdir = "results") {
  fs <- list.files(file.path(outdir, exp), pattern = "\\.csv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  # check.names = FALSE keeps library-set columns such as mca_set_l1+nn
  dl <- lapply(fs, read.csv, stringsAsFactors = FALSE, check.names = FALSE)
  keys <- unique(unlist(lapply(dl, names)))
  do.call(rbind, lapply(dl, function(d) { d[setdiff(keys, names(d))] <- NA; d[keys] }))
}

summarize_exp <- function(exp, outdir = "results", alpha = 0.05) {
  d <- read_exp(exp, outdir)
  if (is.null(d)) return(data.frame())
  d$err[is.na(d$err)] <- ""
  # fixed columns first, then single learners (mca_<name>) and library sets
  extra <- grep("^(mca_([a-z0-9]+|set_.+)|pair_[a-z0-9]+|macc_[a-z0-9]+)$", names(d), value = TRUE)
  pc <- union(intersect(P_COLS, names(d)), extra)
  tc <- grep("^t_", names(d), value = TRUE)
  dc <- union(grep("^(auc_.+_dev|acc_[a-z0-9]+_obs)$", names(d), value = TRUE),
              intersect(c("acc_obs", "ep_ovr_mean"), names(d)))
  dc <- union(dc, grep("^(l1|l1th|l1se|quad|scrm|scr2|knn|rf|nn)_(K|m|log10C|nnz|epochs|nonconv)$", names(d), value = TRUE))
  key <- c("cell", "model", "beta", "k", "n", "d", "tau", "nsplit", "cA", "cB", "sC", "cE",
           intersect("rot", names(d)))
  sp <- split(d, d$cell)
  out <- do.call(rbind, lapply(sp, function(x) {
    ok <- x[x$err == "", , drop = FALSE]
    base <- x[1L, key, drop = FALSE]
    base$reps <- nrow(ok); base$errors <- nrow(x) - nrow(ok)
    for (p in pc) base[[p]] <- if (all(is.na(ok[[p]]))) NA else mean(ok[[p]] <= alpha, na.rm = TRUE)
    for (t in dc) base[[t]] <- mean(ok[[t]], na.rm = TRUE)
    base$sec_per_rep <- mean(ok$t_total, na.rm = TRUE)
    base
  }))
  out <- out[order(out$model, out$k, out$d, out$beta, out$cell), ]
  # drop all-NA method columns
  keep <- vapply(out, function(v) !all(is.na(v)), logical(1))
  out <- out[, keep, drop = FALSE]
  write.csv(out, file.path(outdir, paste0(exp, "_summary.csv")), row.names = FALSE)
  # timing table (E7)
  if (length(tc)) {
    tt <- aggregate(d[tc], by = list(cell = d$cell, k = d$k), FUN = function(v) mean(v, na.rm = TRUE))
    write.csv(tt, file.path(outdir, paste0(exp, "_timing.csv")), row.names = FALSE)
  }
  rownames(out) <- NULL
  out
}
