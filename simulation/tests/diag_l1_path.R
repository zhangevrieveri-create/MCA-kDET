# =============================================================================
# Diagnostic (no code changed): time the ACTUAL learners (learner_l1,
# learner_quad, which fit the whole C path with warm starts) on the slowest E1
# cells, and count glmnet's non-convergence warnings instead of suppressing them.
#   Rscript tests/diag_l1_path.R
# Rebuilds replication 1 of each cell exactly as one_rep does (same seed,
# data, split, standardization, weights) and calls the learner for every
# population i in turn, as fit_scores does. Prints only; writes nothing.
# =============================================================================
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- "R"; PROFILE <- "quick"
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "experiments.R"))
  source(file.path("R", f))
suppressPackageStartupMessages(requireNamespace("glmnet"))

# count (and do not suppress) glmnet warnings inside the learner
count_warn <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  list(val = val, warn = w)
}
# Since 22 Oct 2026 l1_fit_score counts non-convergence itself (info$nonconv);
# with the earlier code (archive/code_v2_20261022) the warnings are caught here.
if (!exists("glmnet_counted"))
  body(l1_fit_score) <- parse(text = paste(gsub("suppressWarnings", "identity", deparse(body(l1_fit_score))), collapse = "\n"))[[1]]
cat("L1_MAXIT =", if (exists("L1_MAXIT")) L1_MAXIT else "1e5 (old code)", "\n")

run_cell <- function(model, k, d, learner) {
  g <- make_grid("E1"); cell <- as.list(g[g$model == model & g$k == k & g$d == d & g$n == "300", ][1, ])
  cat("\n==", cell$cell, "| learner:", learner, "\n")
  set.seed(as.integer(cell$seed0) + 1L)
  data <- gen_data(cell$model, cell$beta, parse_n(cell$n, k), d,
                   list(cA = cell$cA, cB = cell$cB, sC = cell$sC, cE = cell$cE))
  sp <- split_data(data, cell$tau)
  f <- vapply(sp, function(s) nrow(s$fit), integer(1)); m <- vapply(sp, function(s) nrow(s$test), integer(1))
  Xfit <- do.call(rbind, lapply(sp, `[[`, "fit")); Xte <- do.call(rbind, lapply(sp, `[[`, "test"))
  gfit <- rep(seq_len(k), f); st <- standardize_fit(Xfit)
  Zf <- apply_std(Xfit, st); Zt <- apply_std(Xte, st)
  fn <- get(paste0("learner_", learner))
  out <- lapply(seq_len(k), function(i) {
    y <- as.integer(gfit == i); w <- ovr_weights(gfit, i, m / sum(m), f)
    t <- system.time(r <- count_warn(fn(Zf, y, w, gfit, Zt)))[["elapsed"]]
    nc <- if (!is.null(r$val$info$nonconv)) r$val$info$nonconv else sum(grepl("onvergence", r$warn))
    row <- data.frame(i = i, sec = round(t, 1), nonconv = nc, log10C = r$val$info$log10C,
                      nnz = r$val$info$nnz, score_constant = length(unique(r$val$score)) == 1L)
    cat(sprintf("  i = %2d  %6.1f s  nonconv = %d  log10C = %5.1f  nnz = %4d  constant score = %s\n",
                i, row$sec, nc, row$log10C, row$nnz, row$score_constant)); flush.console()
    row
  })
  tab <- do.call(rbind, out)
  cat(sprintf("  total %.1f s, %d non-converged fits\n", sum(tab$sec), sum(tab$nonconv)))
}
run_cell("A", 10L, 300L, "l1")
run_cell("D", 10L, 100L, "quad")
