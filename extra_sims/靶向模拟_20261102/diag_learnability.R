# Learnability diagnostic (seeds 99101+, not in the formal grid): held-out AUC of the
# deviating population for RF / kNN on (a) the 6 true coordinates, (b) deprf-selected
# coordinates, (c) all coordinates. Informative only; no parameter is chosen from it.
source("靶向模拟_20261102/common.R")
res <- list()
for (t in c(0.75)) for (d in c(20L, 50L)) for (r in 1:10) {
  set.seed(99100 + r + 1000 * d)
  data <- tgt_generate("SM3", t, c(300L, 300L, 300L), d)
  sp <- split_data(data, 0.5)
  sc <- fit_scores(sp, "ovr", lib = c("rf", "knn", "deprf"))
  f <- as.integer(table(sc$gfit)); omega <- sc$m / sum(sc$m)
  y <- as.integer(sc$gfit == 3L); w <- ovr_weights(sc$gfit, 3L, omega, f); w <- w / mean(w)
  sel <- sort(unique(as.integer(dep_screen(sc$Zfit, y, w)$pairs)))
  a <- function(s) auc_vec(s, sc$L == 3L)
  res[[length(res) + 1L]] <- data.frame(d = d, rep = r,
    rf_all = a(sc$Sovr$rf[, 3]), knn_all = a(sc$Sovr$knn[, 3]), deprf = a(sc$Sovr$deprf[, 3]),
    rf_true6 = a(learner_rf(sc$Zfit[, 1:6], y, w, sc$gfit, sc$Zte[, 1:6])$score),
    knn_true6 = a(learner_knn(sc$Zfit[, 1:6], y, w, sc$gfit, sc$Zte[, 1:6])$score),
    knn_sel = a(learner_knn(sc$Zfit[, sel], y, w, sc$gfit, sc$Zte[, sel])$score),
    oracle = a({X <- sc$Zte; s <- 0; rho <- TGT_CAL[["0.75"]]
      for (p in 1:3) { x1 <- X[, 2*p-1]; x2 <- X[, 2*p]
        s <- s - rho^2*(x1^2+x2^2)/(2*(1-rho^2)) + log(cosh(rho*x1*x2/(1-rho^2))) }; s}))
}
res <- do.call(rbind, res)
print(aggregate(. ~ d, data = res[, -2], FUN = function(x) round(mean(x), 3)))
