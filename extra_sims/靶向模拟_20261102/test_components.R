# Unit tests of the targeted-simulation code (deprf, adk_max, replication).
# Run from sim/:  Rscript 靶向模拟_20261102/test_components.R
# The seeds here (99001...) are not used by the formal grid (165e6 + ...).
source("靶向模拟_20261102/common.R")
ok <- function(cond, msg) { cat(if (isTRUE(cond)) "OK  " else "FAIL", msg, "\n"); flush.console(); invisible(cond) }
set.seed(99001)
mk_sc <- function(data, lib, tau = 0.5) fit_scores(split_data(data, tau), "ovr", lib = lib)

# 1. wcor equals cov.wt-based correlation
X <- matrix(rnorm(300), 100, 3); w <- runif(100)
ref <- cov.wt(X, wt = w / sum(w), cor = TRUE, method = "ML")$cor
ok(max(abs(wcor(X, X, w) - ref)) < 1e-10, "wcor equals cov.wt correlation")

# 2. screening: true pairs found in the deviating population (ideal AUC 0.75)
for (d in c(20L, 50L, 200L)) {
  data <- tgt_generate("SM3", 0.75, c(300L, 300L, 300L), d)
  sc <- mk_sc(data, "rf")
  f <- as.integer(table(sc$gfit)); omega <- sc$m / sum(sc$m)
  w3 <- ovr_weights(sc$gfit, 3L, omega, f)
  s <- dep_screen(sc$Zfit, as.integer(sc$gfit == 3L), w3 / mean(w3))
  key <- paste(s$pairs[, 1], s$pairs[, 2], sep = "-")
  cat(sprintf("     d = %3d: true pairs among the 10 selected: %d of 3 (%s)\n", d,
              sum(c("1-2", "3-4", "5-6") %in% key), paste(head(key, 5), collapse = " ")))
}
ok(TRUE, "screening run at d = 20, 50, 200 (counts above are informative, not a pass/fail criterion)")

# 3. deprf learns: held-out AUC of population 3 above 0.55 at d = 50, ideal AUC 0.75
data <- tgt_generate("SM3", 0.75, c(300L, 300L, 300L), 50L)
sc <- mk_sc(data, c("rf", "deprf"))
a_rf <- auc_vec(sc$Sovr$rf[, 3], sc$L == 3); a_dp <- auc_vec(sc$Sovr$deprf[, 3], sc$L == 3)
cat(sprintf("     d = 50: held-out AUC of population 3, rf %.3f, deprf %.3f\n", a_rf, a_dp))
ok(a_dp > 0.55, "deprf held-out AUC > 0.55")

# 4. adk_cols against an independent loop implementation of Scholz-Stephens (1987)
adk_loop <- function(x, lab) {
  N <- length(x); o <- order(x); lo <- lab[o]; k <- max(lab); A <- 0
  for (i in 1:k) { ni <- sum(lab == i); Mij <- 0
    for (j in 1:(N - 1)) { Mij <- Mij + (lo[j] == i); A <- A + (N * Mij - j * ni)^2 / (j * (N - j)) / ni } }
  A / N
}
x <- rnorm(60); lab <- rep(1:3, c(15, 20, 25)); Lm <- cbind(lab, sample(lab), sample(lab))
ok(max(abs(adk_cols(x, Lm, 3L) - apply(Lm, 2, function(l) adk_loop(x, l)))) < 1e-10, "adk_cols equals loop implementation")

# 5. null size of deprf alone and adk_max (200 reps, n = 100, d = 20, B = 99)
rej <- t(replicate(200, {
  data <- tgt_generate("SM3", 0.5, c(100L, 100L, 100L), 20L)
  sc <- mk_sc(data, "deprf")
  ct <- classifier_tests(sc, 99L)
  c(ct$p$mca_deprf, adk_max_test(data, 99L)$p) <= 0.05
}))
r <- colMeans(rej)
ok(all(r >= 0.02 & r <= 0.09), sprintf("null rejection deprf %.3f, adk_max %.3f in [0.02, 0.09]", r[1], r[2]))

# 6. one full replication twice with the same seed: identical; timing at d = 20 and d = 200
grid <- readRDS("靶向模拟_20261102/grid.rds")
for (d in c(20L, 200L)) {
  cell <- as.list(grid[grid$model == "SM3" & grid$target == 0.75 & grid$d == d, ])
  cell$B <- 19L
  t1 <- system.time(a <- tgt_one(cell, 1L))[["elapsed"]]
  if (d == 20L) {
    b <- tgt_one(cell, 1L)
    ok(identical(a$p, b$p) && identical(a$aucs, b$aucs) && identical(a$deprf_screen, b$deprf_screen),
       "same seed, identical replication")
    ok(!is.null(a$p$adk_max) && !is.null(a$p[["mca_set_l1+quad+knn+rf+deprf"]]), "p-values include adk_max and L0P")
  }
  cat(sprintf("     one replication, d = %d, B = 19: %.1f s\n", d, t1)); print(round(unlist(a$timing), 2))
}
