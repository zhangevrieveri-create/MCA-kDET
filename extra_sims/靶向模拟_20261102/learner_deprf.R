# =============================================================================
# deprf: dependence-screened random forest (task brief TASK_第三轮_功效与场景_
# 20261101.md, section 2.5, specified on 1 Nov 2026). New file, 2 Nov 2026.
#
# On the fitting data only: for every pair j < l, the absolute difference
# between class 1 and class 0 (class 0 weighted by w, so that it represents
# the self-excluded reference) of each of
#   corr(x_j, x_l), corr(|x_j|, |x_l|), corr(|x_j|, x_l), corr(x_j, |x_l|);
# the pair score is the largest of the four. The DEP_NPAIR = 10 pairs with the
# largest scores are kept (ties: first in the order of upper.tri), and a
# probability forest with the settings of learner_rf is fitted on the
# coordinates of these pairs. It is an ordinary learner: it uses D_fit and one
# seed from the current random-number stream (Assumption 2).
# =============================================================================

DEP_NPAIR <- 10L

wcor <- function(A, B, wt) {               # weighted correlation matrix, columns of A x columns of B
  wt <- wt / sum(wt)
  Ac <- sweep(A, 2L, colSums(A * wt)); Bc <- sweep(B, 2L, colSums(B * wt))
  sa <- sqrt(colSums(Ac^2 * wt)); sb <- sqrt(colSums(Bc^2 * wt))
  sa[sa < 1e-12] <- Inf; sb[sb < 1e-12] <- Inf
  crossprod(Ac * wt, Bc) / outer(sa, sb)
}

# pairs (2-column matrix) selected on (X, y, w); also returns their scores
dep_screen <- function(X, y, w, npair = DEP_NPAIR) {
  d <- ncol(X); A <- abs(X)
  one <- y == 1L; zero <- !one
  delta <- function(f) abs(f(X[one, , drop = FALSE], A[one, , drop = FALSE], w[one]) -
                           f(X[zero, , drop = FALSE], A[zero, , drop = FALSE], w[zero]))
  D1 <- delta(function(X, A, w) wcor(X, X, w))
  D2 <- delta(function(X, A, w) wcor(A, A, w))
  D3 <- delta(function(X, A, w) wcor(A, X, w))           # [j, l] = corr(|x_j|, x_l)
  S <- pmax(D1, D2, D3, t(D3))
  pr <- which(upper.tri(S), arr.ind = TRUE)
  sc <- S[pr]
  top <- order(-sc)[seq_len(min(npair, nrow(pr)))]           # stable: first on ties
  list(pairs = pr[top, , drop = FALSE], score = sc[top])
}

learner_deprf <- function(X, y, w, grp, Xnew, num.trees = 500L, ...) {
  w <- w / mean(w)
  scr <- dep_screen(X, y, w)
  cols <- sort(unique(as.integer(scr$pairs)))
  o <- learner_rf(X[, cols, drop = FALSE], y, w, grp, Xnew[, cols, drop = FALSE], num.trees = num.trees)
  list(score = o$score, info = list(deprf_ncol = length(cols)))
}

LEARNERS2 <- union(LEARNERS2, "deprf")
