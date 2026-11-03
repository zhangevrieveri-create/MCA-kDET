# =============================================================================
# Unit tests of R/garch_ext.R (E8, GARCH null). Single process, writes nothing.
#   "D:/R-4.6.1/bin/Rscript.exe" tests/test_garch.R
# Expected: every line ends in OK; about 1-2 minutes.
# =============================================================================
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- "R"; PROFILE <- "quick"
suppressPackageStartupMessages({ requireNamespace("glmnet"); requireNamespace("ranger") })
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "garch_ext.R",
            "experiments.R"))
  source(file.path("R", f))
ok <- function(cond, msg) { cat(sprintf("%-62s %s\n", msg, if (isTRUE(cond)) "OK" else "FAILED")); invisible(cond) }

g <- make_grid("E8")
cat("E8 cells:\n"); print(g[, c("cell", "d", "reps", "seed0", "mpar")], row.names = FALSE)
ok(nrow(g) == 4L && !anyDuplicated(g$seed0), "1  grid: 4 cells, distinct seeds")

# 2 regime path: fixed per cell, every regime >= minreg weeks
mp <- parse_mpar(g$mpar[1L])
set.seed(g$seed0[1L]); r1 <- gen_regime_path(mp$T, 3L, mp$dur, mp$minreg)
set.seed(g$seed0[1L]); r2 <- gen_regime_path(mp$T, 3L, mp$dur, mp$minreg)
cat("   regime sizes:", tabulate(r1, 3L), "| episodes:", sum(diff(r1) != 0) + 1L, "\n")
ok(identical(r1, r2) && all(tabulate(r1, 3L) >= mp$minreg), "2  regime path reproducible, >= minreg weeks per regime")

# 3 GARCH generator: unit variance, volatility clustering, correlation.
#   Unit variance is checked with Gaussian innovations over 50 columns: with
#   t(5) innovations and (0.10, 0.85) the fourth moment is barely finite
#   (beta^2 + 2 alpha beta + 9 alpha^2 = 0.98), so single-column sample
#   variances scatter widely and converge slowly (2026-10-25: 1.18 in one column).
set.seed(1); Xg <- gen_ccc_garch(20000L, 50L, 0.10, 0.85, Inf, 0, 500L)
set.seed(1); X <- gen_ccc_garch(20000L, 3L, 0.10, 0.85, 5, 0.3, 500L)
cat(sprintf("   Gaussian: mean var %.3f | t(5): var %.3f %.3f %.3f | cor(1,2) %.3f | acf1 of x^2 %.3f\n",
            mean(apply(Xg, 2L, var)), var(X[, 1]), var(X[, 2]), var(X[, 3]), cor(X[, 1], X[, 2]), sq_acf1(X)))
ok(abs(mean(apply(Xg, 2L, var)) - 1) < 0.03 && sq_acf1(X) > 0.05 && abs(cor(X[, 1], X[, 2]) - 0.3) < 0.1,
   "3  CCC-GARCH: variance about 1, clustering, correlation")

# 4 QML filter recovers the persistence, removes the clustering
for (pp in list(c(0.10, 0.85), c(0.05, 0.94))) {
  set.seed(2); x <- gen_ccc_garch(5000L, 1L, pp[1], pp[2], 5, 0, 500L)[, 1]
  fl <- garch11_filter(x)
  cat(sprintf("   true persistence %.2f -> estimate %.3f | acf1 of x^2 %.3f -> z^2 %.3f\n",
              sum(pp), fl$persist, sq_acf1(matrix(x)), sq_acf1(matrix(fl$z))))
  ok(abs(fl$persist - sum(pp)) < 0.05 && abs(sq_acf1(matrix(fl$z))) < 0.05,
     sprintf("4  filter, persistence %.2f", sum(pp)))
}

# 5 block split: buffer respected, enough weeks per regime and role
set.seed(3); rb <- block_roles(r1, 3L, 0.5, 26L, 4L, 30L)
fitw <- which(rb == "fit"); testw <- which(rb == "test")
mind <- min(vapply(testw, function(t) min(abs(fitw - t)), numeric(1)))
cat(sprintf("   fit %d | test %d | dropped %d | min distance test-fit %d weeks\n",
            length(fitw), length(testw), sum(rb == "drop"), mind))
ok(mind > 4 && all(tabulate(r1[rb == "fit"], 3L) >= 30) && all(tabulate(r1[rb == "test"], 3L) >= 30),
   "5  block split with buffer")
ri <- iid_roles(r1, 3L, 0.5)
ok(all(abs(tabulate(r1[ri == "fit"], 3L) - floor(0.5 * tabulate(r1, 3L))) == 0), "6  random within-regime split")

# 7 one replication (B = 19): runs, p-values in (0, 1], reproducible
cell <- as.list(g[1L, ]); cell$B <- 19L
t0 <- proc.time()[["elapsed"]]; a <- one_rep_garch(cell, 1L); t1 <- proc.time()[["elapsed"]] - t0
b <- one_rep_garch(cell, 1L)
pc <- grep("^((flt_iid_|raw_blk_|raw_iid_|raw_)?(mca|disco|mmd|mst))$", names(a), value = TRUE)
cat("   p-values:", paste(sprintf("%s=%.2f", pc, unlist(a[pc])), collapse = " "), "\n")
cat(sprintf("   %.1f s per replication (B = 19); sq-acf raw %.3f, filtered %.3f\n",
            t1, a$sqacf_raw, a$sqacf_flt))
ok(length(pc) == 10L && all(unlist(a[pc]) > 0 & unlist(a[pc]) <= 1), "7  one replication: 10 p-values in (0, 1]")
ok(identical(a[pc], b[pc]), "8  same seed, same p-values")

# 9 the other models are untouched by the dispatch in one_rep()
cs <- cell_row("smoke", "A", 0, n = 40, d = 10, B = 19, reps = 1)
ok(is.list(one_rep(as.list(cs), 1L)), "9  one_rep() for model A still runs")
