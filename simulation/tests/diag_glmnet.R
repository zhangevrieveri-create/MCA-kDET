# =============================================================================
# Diagnostic (no code changed): where does glmnet spend its time in the slow
# E1 cells? Rebuilds the first one-vs-rest problem of replication 1 of a cell
# (same seed, data, split and standardization as one_rep) and times glmnet at
# every value of the C grid, on all fitting data and on one CV training fold.
#
#   Rscript tests/diag_glmnet.R A 10 300 l1      # E1 A null, k = 10, d = 300, l1
#   Rscript tests/diag_glmnet.R D 10 100 quad    # E1 D null, k = 10, d = 100, quad
# =============================================================================
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
a <- commandArgs(TRUE); model <- a[1]; k <- as.integer(a[2]); d <- as.integer(a[3]); which <- a[4]
SRC_DIR <- "R"; PROFILE <- "quick"
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "experiments.R"))
  source(file.path("R", f))
suppressPackageStartupMessages(requireNamespace("glmnet"))
g <- make_grid("E1"); cell <- as.list(g[g$model == model & g$k == k & g$d == d & g$n == "300", ][1, ])
cat("cell:", cell$cell, "\n")

# replicate one_rep up to the one-vs-rest problem of population 1
set.seed(as.integer(cell$seed0) + 1L)
nvec <- parse_n(cell$n, k)
data <- gen_data(cell$model, cell$beta, nvec, d, list(cA = cell$cA, cB = cell$cB, sC = cell$sC, cE = cell$cE))
sp <- split_data(data, cell$tau)
f <- vapply(sp, function(s) nrow(s$fit), integer(1)); m <- vapply(sp, function(s) nrow(s$test), integer(1))
Xfit <- do.call(rbind, lapply(sp, `[[`, "fit")); gfit <- rep(seq_len(k), f)
st <- standardize_fit(Xfit); Z <- apply_std(Xfit, st)
i <- 1L; y <- as.integer(gfit == i); w <- ovr_weights(gfit, i, m / sum(m), f); w <- w / mean(w)
if (which == "quad") {               # quad features as in learner_quad
  dd <- ncol(Z); npair <- dd * (dd - 1) / 2
  if (npair <= QUAD_MAX_PAIRS) pairs <- which(upper.tri(diag(dd)), arr.ind = TRUE) else {
    one <- y == 1L
    D <- abs(wcov(Z[one, , drop = FALSE], w[one]) - wcov(Z[!one, , drop = FALSE], w[!one]))
    ut <- which(upper.tri(D), arr.ind = TRUE); pairs <- ut[order(-D[ut])[seq_len(QUAD_MAX_PAIRS)], , drop = FALSE] }
  F <- quad_features(Z, pairs); Z <- apply_std(F, standardize_fit(F))
}
cat(sprintf("design: n = %d, p = %d, positives = %d\n", nrow(Z), ncol(Z), sum(y)))

time_fit <- function(X, yy, ww, label) {
  n <- length(yy); ww <- ww / mean(ww)
  out <- lapply(sort(L1_C_GRID), function(C) {
    t <- system.time(fit <- suppressWarnings(glmnet::glmnet(X, yy, family = "binomial", weights = ww,
           lambda = 1 / (C * n), standardize = FALSE, intercept = TRUE, maxit = 1e5)))[["elapsed"]]
    data.frame(set = label, log10C = log10(C), lambda = signif(1 / (C * n), 3), sec = round(t, 2),
               npasses = fit$npasses, dev_ratio = round(fit$dev.ratio, 4), nonzero = fit$df,
               jerr = fit$jerr)
  })
  do.call(rbind, out)
}
res <- time_fit(Z, y, w, "all fitting data (one lambda per call)")
fold <- strat_folds(gfit, 3L); trn <- fold != 1L
res <- rbind(res, time_fit(Z[trn, ], y[trn], w[trn], "CV training fold 1"))
print(res, row.names = FALSE)
# the path as l1_fit_score fits it (all six lambdas in one call, strongest penalty first)
t_path <- system.time(fp <- suppressWarnings(glmnet::glmnet(Z, y, family = "binomial", weights = w,
  lambda = 1 / (sort(L1_C_GRID) * length(y)), standardize = FALSE, intercept = TRUE, maxit = 1e5)))[["elapsed"]]
cat(sprintf("\nwhole path in one call (as in l1_fit_score): %.1f s; lambdas fitted: %d of 6; dev.ratio: %s\n",
            t_path, length(fp$lambda), paste(round(fp$dev.ratio, 4), collapse = ", ")))
