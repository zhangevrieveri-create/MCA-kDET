# =============================================================================
# Checks of the code added for the third review round (1 Nov 2026):
# R/competitors2.R (G1), R/rotation.R (G1), the accuracy columns of
# classifier_tests() (G2) and the paired grids E9, E10. Writes nothing.
#
#   "D:/R-4.6.1/bin/Rscript.exe" tests/test_review3.R
#
# About 2-4 minutes. Every line starts with OK or FAIL; the last line says
# REVIEW3 TESTS OK or REVIEW3 TESTS FAILED. Optional packages: gTestsMulti
# (Song-Chen, needed for E9/E10) and KMD (only for check 4).
# =============================================================================

script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R",
            "competitors2.R", "rotation.R", "experiments.R"))
  source(file.path(SRC_DIR, f))
source(file.path(SRC_DIR, "summarize.R"))
suppressPackageStartupMessages(requireNamespace("glmnet"))
suppressPackageStartupMessages(requireNamespace("ranger", quietly = TRUE))

all_ok <- TRUE
chk <- function(cond, msg) {
  cond <- isTRUE(cond)
  cat(if (cond) "OK   " else "FAIL ", msg, "\n", sep = "")
  all_ok <<- all_ok && cond
}
par0 <- list(cA = 0.6, cB = 4, sC = 2L, cE = 0.2)
has_sc <- requireNamespace("gTestsMulti", quietly = TRUE)

# ---- 1 Schott statistic at the observed labels equals schott_test()
set.seed(1)
dat <- gen_data("A", 0, c(40L, 50L, 60L), 8L, par0)
X <- do.call(rbind, dat); L <- rep(1:3, c(40L, 50L, 60L))
cs <- cov_perm_stats(X, cbind(L), 3L)
chk(abs(cs[1L, "schott"] - schott_test(dat)[["stat"]]) < 1e-8,
    sprintf("1 Schott z at observed labels: %.10f vs schott_test %.10f",
            cs[1L, "schott"], schott_test(dat)[["stat"]]))

# ---- 2 CLX covariance statistic against a direct double loop
direct_clx <- function(dat) {
  k <- length(dat); d <- ncol(dat[[1L]]); best <- -Inf
  for (a in 1:(k - 1)) for (b in (a + 1):k) for (j in 1:d) for (l in j:d) {
    f <- function(Y) {
      y1 <- Y[, j] - mean(Y[, j]); y2 <- Y[, l] - mean(Y[, l])
      s <- mean(y1 * y2); c(s, mean((y1 * y2 - s)^2), nrow(Y))
    }
    u <- f(dat[[a]]); v <- f(dat[[b]])
    best <- max(best, (u[1] - v[1])^2 / (u[2] / u[3] + v[2] / v[3]))
  }
  best
}
chk(abs(cs[1L, "clx"] - direct_clx(dat)) < 1e-8 * max(1, abs(cs[1L, "clx"])),
    sprintf("2 CLX covariance statistic: %.8f vs direct %.8f", cs[1L, "clx"], direct_clx(dat)))

# ---- 3 max-type mean statistic against a direct computation
G <- group_indicator(cbind(L), 3L)
v <- apply(X, 2L, var)
mm <- maxmean_stat(X, G, c(40, 50, 60), 3L, 1L, v)
dm <- max(sapply(list(c(1, 2), c(1, 3), c(2, 3)), function(p)
  max((colMeans(dat[[p[1]]]) - colMeans(dat[[p[2]]]))^2 / v /
      (1 / nrow(dat[[p[1]]]) + 1 / nrow(dat[[p[2]]])))))
chk(abs(mm - dm) < 1e-10 * max(1, dm), sprintf("3 max-type mean statistic: %.8f vs direct %.8f", mm, dm))

# ---- 4 KMD against the KMD package (if installed)
D2 <- sqdist(X)
eta <- kmd_stat(D2, cbind(L), ceiling(nrow(X) / 10))
if (requireNamespace("KMD", quietly = TRUE)) {
  ek <- KMD::KMD(X, L, M = 3, Knn = ceiling(nrow(X) / 10), Kernel = "discrete")
  chk(abs(eta - ek) < 1e-8, sprintf("4 KMD: %.8f vs KMD::KMD %.8f", eta, ek))
} else cat("SKIP 4 KMD package not installed (own value", round(eta, 5), ")\n")

# ---- 5 rotation: orthogonal, fixed, and the random-number stream is unchanged
Q1 <- fixed_rotation(20L); Q2 <- fixed_rotation(20L)
chk(max(abs(crossprod(Q1) - diag(20))) < 1e-10 && identical(Q1, Q2),
    "5a rotation matrix orthogonal and reproducible")
set.seed(7); d0 <- gen_data("C", 1, c(30L, 30L, 30L), 20L, par0); u0 <- runif(3)
set.seed(7); d1 <- gen_data("C", 1, c(30L, 30L, 30L), 20L, par0); d1 <- rotate_data(d1); u1 <- runif(3)
chk(identical(u0, u1) && max(abs(d1[[3]] - d0[[3]] %*% Q1)) < 1e-12,
    "5b rotated data = data %*% Q, random-number stream unchanged")

# ---- 6 size of the comp2 tests (A null, k = 3, n = 40, d = 10, B = 99, 300 reps)
reps <- 300L
P <- t(vapply(seq_len(reps), function(r) {
  set.seed(5000L + r)
  dd <- gen_data("A", 0, c(40L, 40L, 40L), 10L, par0)
  p <- comp2_tests(dd, 99L, "comp2")$p
  unlist(p[c("maxmean", "clx_mean", "clx_cov", "schott_perm", "mmd_max", "kmd")])
}, numeric(6)))
rr <- colMeans(P <= 0.05)
chk(all(rr >= 0.01 & rr <= 0.095),
    paste("6 size in [0.01, 0.095] over 300 null reps:",
          paste(sprintf("%s %.3f", names(rr), rr), collapse = ", ")))

# ---- 7 power sanity checks (strong signals, 20 reps)
pw <- function(model, beta, cpar, cols) {
  pp <- t(vapply(1:20, function(r) {
    set.seed(9000L + r)
    dd <- gen_data(model, beta, c(60L, 60L, 60L), 10L, modifyList(par0, cpar))
    unlist(comp2_tests(dd, 99L, "comp2")$p[cols])
  }, numeric(length(cols))))
  colMeans(pp <= 0.05)
}
p_mean <- pw("A", 1, list(cA = 3), c("maxmean", "clx_mean"))
chk(all(p_mean >= 0.8), paste("7a mean shift detected:", paste(sprintf("%.2f", p_mean), collapse = " ")))
p_cov <- pw("B", 1, list(cB = 20), c("clx_cov", "schott_perm"))
chk(all(p_cov >= 0.8), paste("7b variance shift detected:", paste(sprintf("%.2f", p_cov), collapse = " ")))

# ---- 8 Song-Chen tests run and return p-values
if (has_sc) {
  set.seed(11)
  dd <- gen_data("A", 1, c(40L, 40L, 40L), 10L, modifyList(par0, list(cA = 3)))
  r8 <- comp2_tests(dd, 99L, "songchen")$p
  chk(is.numeric(r8$songchen) && is.numeric(r8$songchen_A) &&
        all(c(r8$songchen, r8$songchen_A) > 0 & c(r8$songchen, r8$songchen_A) <= 1),
      sprintf("8 Song-Chen p-values: S %.3f, S^A %.3f", r8$songchen, r8$songchen_A))
  E5 <- kmst_edges(sqrt(sqdist(do.call(rbind, dd))), SC_KMST)
  chk(nrow(E5) == SC_KMST * (120L - 1L) && !anyDuplicated(t(apply(E5, 1L, sort))),
      sprintf("8b %d-MST: %d distinct edges", SC_KMST, nrow(E5)))
} else cat("FAIL 8 gTestsMulti not installed: install.packages(\"gTestsMulti\") before E9/E10\n")
if (!has_sc) all_ok <- FALSE

# ---- 9 paired grids: E9 and E10 use the seeds of their E1/E2 source cells
src <- rbind(make_grid("E1"), make_grid("E2")); src$key <- sub("^E[0-9]+_", "", src$cell)
for (x in c("E9", "E10")) {
  g <- make_grid(x); m <- match(sub("^E[0-9]+_", "", g$cell), src$key)
  chk(!anyNA(m) && all(g$seed0 == src$seed0[m]) && all(g$reps == src$reps[m]) &&
        all(g$lib == src$lib[m]),
      sprintf("9 %s: %d cells paired with E1/E2 (same seed0, reps, lib)", x, nrow(g)))
}

# ---- 10 compatibility: an E2 replication of the quick run is reproduced,
#         every old column unchanged, the new accuracy columns added
same15 <- function(a, b) isTRUE(all(abs(a - b) <= 1e-13 * pmax(abs(a), abs(b)) | a == b))
e2 <- make_grid("E2")
# E2_D_b0.03 was run last (22 Oct 2026, after the glmnet change of that day,
# which left every p-value unchanged), so its chunk reflects the current learners
id <- e2$cell[e2$model == "D" & abs(e2$beta - 0.03) < 1e-9]
f10 <- file.path("results_v2_quick", "E2", sprintf("%s__r00001-00050.csv", id))
if (length(id) == 1L && file.exists(f10)) {
  old <- read.csv(f10, stringsAsFactors = FALSE, check.names = FALSE)
  o <- old[old$rep == 1L, ]
  new <- one_rep(as.list(e2[e2$cell == id, ]), 1L)
  common <- intersect(names(new), names(o))
  pcols <- intersect(union(P_COLS, grep("^(mca_|pair_)", common, value = TRUE)), common)
  p_same <- vapply(pcols, function(x)
    if (x == "schott") same15(as.numeric(new[[x]]), o[[x]])
    else same15(as.numeric(new[[x]]), o[[x]]), logical(1))   # csv keeps 15 digits
  chk(all(p_same), sprintf("10a %s rep 1: %d old p-value columns identical", id, length(pcols)))
  if (!all(p_same)) print(data.frame(col = pcols[!p_same], old = unlist(o[pcols[!p_same]]),
                                     new = unlist(new[pcols[!p_same]])), row.names = FALSE)
  chk(all(c("macc_ovr", "macc_l1", "macc_quad", "macc_knn", "macc_rf") %in% names(new)),
      sprintf("10b accuracy columns added: macc_ovr = %.3f", new$macc_ovr))
} else cat("SKIP 10 results_v2_quick/E2 chunk of", id, "not found\n")

cat(if (all_ok) "\nREVIEW3 TESTS OK\n" else "\nREVIEW3 TESTS FAILED\n")
