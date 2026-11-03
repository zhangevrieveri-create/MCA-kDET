# =============================================================================
# Unit tests for the learner library (TASK_学习器库.md, section 4)
#
#   Rscript tests/test_learners.R --part=1               # direction and signal (serial, fast)
#   Rscript tests/test_learners.R --part=2 --cores=30    # size, 200 null reps
#   Rscript tests/test_learners.R --part=3               # reproducibility (serial)
#
# Output is printed; part 2 also writes tests/out/level_reps.csv.
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
setwd(file.path(script_dir, ".."))
SRC_DIR <- file.path(getwd(), "R")
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "experiments.R"))
  source(file.path(SRC_DIR, f))
suppressPackageStartupMessages({ requireNamespace("glmnet"); requireNamespace("ranger") })

part <- as.integer(get_arg("part", 1L))
ALL <- c("l1", "quad", "knn", "rf", "nn")

# ---------------------------------------------------------------- part 1
# Toy two-class problems, n = 200 per class for fitting and 200 per class for
# testing, unit weights. A learner "should have signal" if its test AUC > 0.6.
if (part == 1L) {
  set.seed(20261021L)
  n <- 200L
  toy <- list(
    mean = list(d = 2L, expect = ALL,
                gen = function(n, cls) { X <- matrix(rnorm(n * 2), n, 2); if (cls) X[, 1] <- X[, 1] + 1; X }),
    var  = list(d = 2L, expect = c("quad", "rf"),
                gen = function(n, cls) { X <- matrix(rnorm(n * 2), n, 2); if (cls) X[, 1] <- 2 * X[, 1]; X }),
    Cprime = list(d = 2L, expect = c("knn", "rf", "nn"),        # x2 = S x1, S = +-1 (model C', beta = 1)
                gen = function(n, cls) gen_population("C", if (cls) 1 else 0, n, 2L, if (cls) 2L else 1L, 2L, list(sC = 1L))),
    cov  = list(d = 5L, expect = "quad",                         # corr(x1, x2) = 0.8 vs 0
                gen = function(n, cls) {
                  X <- matrix(rnorm(n * 5), n, 5)
                  if (cls) X[, 2] <- 0.8 * X[, 1] + 0.6 * X[, 2]
                  X }))
  tab <- matrix(NA_real_, length(toy), length(ALL), dimnames = list(names(toy), ALL))
  info <- character(0)
  for (p in names(toy)) {
    tp <- toy[[p]]
    Xf <- rbind(tp$gen(n, FALSE), tp$gen(n, TRUE)); yf <- rep(0:1, each = n)
    Xt <- rbind(tp$gen(n, FALSE), tp$gen(n, TRUE)); yt <- rep(0:1, each = n)
    st <- standardize_fit(Xf)
    Zf <- apply_std(Xf, st); Zt <- apply_std(Xt, st)
    for (a in ALL) {
      o <- get(paste0("learner_", a))(Zf, yf, rep(1, 2 * n), yf + 1L, Zt)
      tab[p, a] <- auc_vec(o$score, yt == 1)
      if (length(o$info)) info <- c(info, sprintf("%-7s %-5s %s", p, a,
        paste(names(o$info), signif(unlist(o$info), 3), sep = "=", collapse = " ")))
    }
  }
  cat("\nTest AUC (* = signal expected; FAIL if expected and AUC <= 0.6)\n\n")
  show <- matrix(sprintf("%.3f", tab), nrow(tab), dimnames = dimnames(tab))
  for (p in names(toy)) for (a in toy[[p]]$expect)
    show[p, a] <- paste0(show[p, a], if (tab[p, a] > 0.6) "*" else "* FAIL")
  print(noquote(show))
  cat("\nSelected tuning values:\n"); cat(info, sep = "\n")
  ok <- all(vapply(names(toy), function(p) all(tab[p, toy[[p]]$expect] > 0.6), logical(1)))
  cat(if (ok) "\nPART 1 OK\n" else "\nPART 1: some expected signals missing\n")
}

# ---------------------------------------------------------------- part 2
# Size: k = 3, n = 60, d = 5, null of model A (AR(0.5) Gaussian), B = 99,
# 200 reps. All five learners are fitted in every rep; mca_<name> is the
# p-value of the library {name} alone (same permutations), mca the whole lib.
if (part == 2L) {
  cores <- as.integer(get_arg("cores", 30L))
  reps <- as.integer(get_arg("reps", 200L))
  cell <- as.list(cell_row("unittest", "A", 0, k = 3L, n = 60L, d = 5L, B = 99L,
                           groups = "ovr", reps = reps, lib = paste(ALL, collapse = ","),
                           libsets = "", tag = "_level"))
  cl <- make_cluster(cores, SRC_DIR)
  t0 <- proc.time()[["elapsed"]]
  rows <- parallel::parLapplyLB(cl, seq_len(reps), function(r, cell) safe_rep(cell, r), cell = cell)
  parallel::stopCluster(cl)
  el <- proc.time()[["elapsed"]] - t0
  df <- rows_to_df(rows)
  dir.create("tests/out", showWarnings = FALSE)
  write.csv(df, "tests/out/level_reps.csv", row.names = FALSE)
  ok <- df$err == ""
  cat(sprintf("\n%d reps in %.0f s, %d errors\n", reps, el, sum(!ok)))
  if (any(!ok)) cat("first error:", df$err[!ok][1L], "\n")
  cols <- c(paste0("mca_", ALL), "mca", "mca_2s", "mca_sum", "mca_comb")
  rr <- vapply(cols, function(p) mean(df[[p]][ok] <= 0.05), numeric(1))
  cat("\nRejection rate at alpha = 0.05 (target band [0.02, 0.08])\n")
  print(data.frame(p = cols, rate = round(rr, 3),
                   in_band = rr >= 0.02 & rr <= 0.08), row.names = FALSE)
  tc <- grep("^t_ovr_", names(df), value = TRUE)
  cat("\nMean seconds per rep per learner:\n"); print(round(colMeans(df[ok, tc, drop = FALSE]), 3))
}

# ---------------------------------------------------------------- part 3
# Reproducibility: the same cell and rep run twice give identical rows
# (timing columns excluded).
if (part == 3L) {
  cell <- as.list(cell_row("unittest", "D", 1, k = 3L, n = 60L, d = 5L, B = 19L,
                           groups = "ovr,dist,schott", reps = 1L,
                           lib = paste(ALL, collapse = ","),
                           libsets = "l1+nn;l1+quad+knn+nn;l1+quad+knn+rf+nn", tag = "_repro"))
  a <- one_rep(cell, 1L); b <- one_rep(cell, 1L)
  drop <- function(x) x[!grepl("^t_", names(x))]
  same <- identical(drop(a), drop(b))
  cat("\nColumns compared:", length(drop(a)), "\n")
  print(unlist(drop(a)[grepl("^(mca|auc_|knn_|quad_|l1_)", names(drop(a)))]))
  cat(if (same) "\nPART 3 OK: identical\n" else "\nPART 3 FAILED: rows differ\n")
}
