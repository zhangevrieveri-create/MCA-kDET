# =============================================================================
# Checks of the quick E1-E7 results against the theory of Sections 2-3
# (Theorem 1, Remark 12, Theorem 2 / Corollary 1, Remark 11, Remark 10).
# Reads results_v2_quick/ and results/libsize only; writes nothing else.
#   Rscript analysis_theory_checks.R > results_v2_quick/theory_checks.txt
# =============================================================================
source("R/summarize.R")
rd <- function(e, o = "results_v2_quick") { d <- read_exp(e, o); d$err[is.na(d$err)] <- ""; d[d$err == "", ] }
E1 <- rd("E1"); E2 <- rd("E2"); E3 <- rd("E3"); E4 <- rd("E4"); E5 <- rd("E5"); E6 <- rd("E6")
# E2_D_b0.01 (s = 1) is a null case, replaced by s = 3 on 22 Oct 2026 (R/calib_params_v3.R);
# its quick results stay on disk but are not an alternative
E2 <- E2[!grepl("^E2_D_b0\\.01_", E2$cell), ]
cp <- function(x, n) binom.test(x, n)$conf.int

cat("=== 1. Theorem 1: pooled size (all null replications of the quick runs)\n")
B <- 200; bound <- floor(0.05 * (B + 1)) / (B + 1)
cat(sprintf("bound floor(alpha(B+1))/(B+1) = %.4f\n", bound))
nulls <- list(E1 = E1, E4 = E4[E4$beta == 0, ], E5 = E5[E5$beta == 0, ])
vars <- c("mca","mca_l1","mca_quad","mca_knn","mca_rf","mca_2s","mca_sum","mca_comb","pair",
          "shared_nn","macc","mca_ms_med","mca_ms_mean")
for (v in vars) {
  x <- unlist(lapply(nulls, function(d) if (v %in% names(d)) d[[v]] else NULL)); x <- x[!is.na(x)]
  k <- sum(x <= 0.05); ci <- cp(k, length(x))
  cat(sprintf("  %-12s n = %5d  rate = %.4f  95%% CI [%.4f, %.4f]  %s\n", v, length(x), k / length(x),
              ci[1], ci[2], if (ci[1] > bound) "ABOVE BOUND" else "consistent"))
}
x <- E1$mca
cat(sprintf("  MCA in E1 at alpha = 0.01 / 0.05 / 0.10: %.4f / %.4f / %.4f (n = %d)\n",
            mean(x <= .01), mean(x <= .05), mean(x <= .10), length(x)))
cat(sprintf("  KS test of MCA p-values (E1) against U(0,1): p = %.3f (p-values lie on a grid of 1/201)\n",
            suppressWarnings(ks.test(x, "punif")$p.value)))

cat("\n=== 2. Remark 12: null distribution of the standardized statistic T (deviant population, E1)\n")
for (a in c("l1", "quad", "knn", "rf")) {
  t <- E1[[paste0("T_", a, "_dev")]]
  cat(sprintf("  T_%-5s mean %+.3f  sd %.3f  P(T = 0) %.3f  P(T > 1.645) %.3f\n", a, mean(t), sd(t),
              mean(t == 0), mean(t > 1.645)))
}

cat("\n=== 3. Theorem 2 / Corollary 1: power is driven by the learned signal\n")
kc <- c("mca", "auc_l1_dev", "auc_quad_dev", "auc_knn_dev", "auc_rf_dev")
alt <- rbind(E2[E2$beta > 0, kc], E3[E3$k == 3, kc], E6[E6$d %in% c(20, 100), kc])
auc <- pmax(alt$auc_l1_dev, alt$auc_quad_dev, alt$auc_knn_dev, alt$auc_rf_dev)
m <- 150; M <- 450; rho <- pmax(auc - 0.5, 0)^2 * m * (M - m) / M
kappa <- (M - m) * (M + 1) / (6 * (M - 1)^2); kq <- 3 * 4
br <- c(0.4, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.65, 0.70, 0.80, 1)
bin <- cut(auc, br, right = FALSE)
tab <- do.call(rbind, lapply(split(seq_along(auc), bin), function(ix) {
  r <- median(rho[ix]); bnd <- min(1, exp(-r / 2) + B / floor(0.05 * (B + 1)) * kq * exp(-3 * kappa * M * r / (M + 1)))
  data.frame(n = length(ix), power = round(mean(alt$mca[ix] <= 0.05), 3),
             median_rho = round(r, 1), thm2_lower_bound_on_power = round(1 - bnd, 3))
}))
cat("  best deviant test AUC over the 4 learners (k = 3 cells of E2, E3, E6):\n"); print(tab)
cat(sprintf("  null sd of the AUC: %.4f, so AUC 0.58 is about %.1f null sd above 1/2\n",
            sqrt((M + 1) / (12 * m * (M - m))), 0.08 / sqrt((M + 1) / (12 * m * (M - m)))))

cat("\n=== 4. Remark 11: given the learned signal, power does not depend on d (E6, model A)\n")
a6 <- E6[E6$model == "A", ]; auc6 <- pmax(a6$auc_l1_dev, a6$auc_quad_dev, a6$auc_knn_dev, a6$auc_rf_dev)
b6 <- cut(auc6, c(0.4, 0.55, 0.57, 0.59, 0.62, 1), right = FALSE)
print(round(tapply(a6$mca <= 0.05, list(auc_bin = b6, d = a6$d), mean), 2))
print(table(auc_bin = b6, d = a6$d))

cat("\n=== 5. Remark 10: price of the library in model A, beta = 0.8 (E2)\n")
a2 <- E2[E2$model == "A" & E2$beta == 0.8, ]
sdA <- sqrt((M + 1) / (12 * m * (M - m)))
for (z in c(qnorm(1 - 0.05 / 3), qnorm(1 - 0.05 / 12)))
  cat(sprintf("  P(AUC_l1_dev > 1/2 + %.2f null sd) = %.2f\n", z, mean(a2$auc_l1_dev > 0.5 + z * sdA)))
cat(sprintf("  observed: l1 alone %.2f, library %.2f\n", mean(a2$mca_l1 <= .05), mean(a2$mca <= .05)))
