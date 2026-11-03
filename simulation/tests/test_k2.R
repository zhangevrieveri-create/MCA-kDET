# =============================================================================
# 检查 k = 2 的实现（注 4）：每个学习器只拟合一个分类器，且 T_2 = T_1。
# 在 sim\ 目录下运行：
#   "D:/R-4.6.1/bin/Rscript.exe" tests/test_k2.R
# =============================================================================
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R")) source(file.path("R", f))
suppressPackageStartupMessages(requireNamespace("glmnet"))
set.seed(20261031)
par <- list(cA = 0.6, cB = 4, sC = 10, cE = 0.2)
data <- gen_data("A", 1, c(80L, 80L), 10L, par)
sp <- split_data(data, 0.5)

calls <- 0L
orig_l1 <- learner_l1
learner_l1 <- function(...) { calls <<- calls + 1L; orig_l1(...) }
sc <- fit_scores(sp, "ovr", lib = c("l1", "knn"))
cat("learner_l1 calls:", calls, "(expected 1)\n")
stopifnot(calls == 1L)

L <- perm_labels(sc$L, 19L)
for (a in sc$lib) {
  Tm <- ovr_T(sc$Sovr[[a]], 1:2, L)
  cat(a, ": max |T_1 - T_2| over 20 labellings =", max(abs(Tm[1, ] - Tm[2, ])), "\n")
  stopifnot(isTRUE(all.equal(Tm[1, ], Tm[2, ])))
}
ct <- classifier_tests(sc, 19L)
print(unlist(ct$p))
cat("k = 2 test passed\n")
