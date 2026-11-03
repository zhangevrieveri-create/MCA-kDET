# =============================================================================
# Third review round (G1/G2 follow-up): sparser learners and screening
# components (TASK_第三轮_功效与场景_20261101.md, section 2). Written 1 Nov 2026.
#
# Learners (same interface as R/learners_lib.R, appear in lib = "..."):
#   l1th  L1 logistic with the theoretical penalty lambda = c (2 log d / N_fit)^(1/2)
#   l1se  L1 logistic, penalty by the 1-SE rule (strongest penalty within one
#         standard error of the best cross-validated AUC)
#   scrm  mean screening: top-m standardized mean differences, m in {1, 3, 10}
#   scr2  second-order screening on x_j^2 and screened x_j x_l, same rule
#
# Components (appear in lib = "..." but are NOT one-vs-rest classifiers; they
# are computed in classifier_tests() from the test labels of each permutation
# and enter the same maximum, with the same permutations):
#   xscr1, xscr3  cross-fitted screening: select the m best of {x_j, x_j^2,
#                 x_j x_l} on one half (labels of that half), AUC on the other,
#                 both directions, T = (T1 + T2)/sqrt(2). For every permuted
#                 labelling the selection in direction 2 is redone from the
#                 permuted test labels (task section 2.6, validity: given
#                 D_fit, T is a function of the test points and test labels)
#   g5t           within-population edge count of a 5-MST on the test points
#   g5f           same on fit + test points, fit labels fixed, edges counted if
#                 at least one end is a test point
# A learner uses only the fitting data and the random-number stream; a
# component is a function of D_fit and the test points / test labels only.
# =============================================================================

LEARNERS2 <- c("l1th", "l1se", "scrm", "scr2")
COMPONENTS <- c("xscr1", "xscr3", "g5t", "g5f")
SCREEN_MS <- c(1L, 3L, 10L)
XSCR_MAX_D_CROSS <- 100L        # cross products x_j x_l only for d <= 100

# ---------------------------------------------------------------- l1th
learner_l1th <- function(X, y, w, grp, Xnew, cth = 1, ...) {
  N <- nrow(X); d <- ncol(X)
  w <- w / mean(w)
  lam <- cth * sqrt(2 * log(max(d, 2L)) / N)
  g <- glmnet_counted(X, y, family = "binomial", weights = w, lambda = lam,
                      standardize = FALSE, intercept = TRUE, maxit = L1_MAXIT)
  P <- predict(g$fit, Xnew, type = "link")
  list(score = as.numeric(P[, ncol(P)]),
       info = list(nnz = sum(g$fit$beta[, ncol(P)] != 0), nonconv = g$nonconv))
}

# ---------------------------------------------------------------- l1se
# same grid and folds as l1_fit_score(); the penalty is chosen by the 1-SE rule
l1_fit_score_1se <- function(X, y, w, grp, Xnew, Cs = L1_C_GRID, nfold = 3L) {
  nonconv <- 0L
  w <- w / mean(w)
  fold <- integer(length(y))
  for (g in unique(grp)) {
    ix <- which(grp == g)
    fold[ix] <- sample(rep_len(seq_len(nfold), length(ix)))
  }
  Cs <- Cs[order(Cs)]
  cv <- matrix(NA_real_, nfold, length(Cs))
  for (f in seq_len(nfold)) {
    trn <- fold != f; val <- !trn
    if (length(unique(y[trn])) < 2L || length(unique(y[val])) < 2L) next
    wt <- w[trn] / mean(w[trn])
    lam <- 1 / (Cs * sum(trn))
    g <- glmnet_counted(X[trn, , drop = FALSE], y[trn], family = "binomial", weights = wt,
                        lambda = lam, standardize = FALSE, intercept = TRUE, maxit = L1_MAXIT)
    fit <- g$fit; nonconv <- nonconv + g$nonconv
    P <- predict(fit, X[val, , drop = FALSE], type = "link")
    for (c in seq_len(min(ncol(P), length(Cs)))) cv[f, c] <- auc_vec(P[, c], y[val] == 1)
  }
  mcv <- colMeans(cv, na.rm = TRUE)
  mcv[!is.finite(mcv)] <- -Inf
  cbest <- which(mcv == max(mcv))[1L]
  se <- stats::sd(cv[, cbest], na.rm = TRUE) / sqrt(sum(!is.na(cv[, cbest])))
  if (!is.finite(se)) se <- 0
  c1se <- which(mcv >= mcv[cbest] - se)[1L]       # smallest C = strongest penalty
  lam_all <- 1 / (Cs * length(y))
  g <- glmnet_counted(X, y, family = "binomial", weights = w, lambda = lam_all[seq_len(c1se)],
                      standardize = FALSE, intercept = TRUE, maxit = L1_MAXIT)
  fit <- g$fit; nonconv <- nonconv + g$nonconv
  P <- predict(fit, Xnew, type = "link")
  list(score = as.numeric(P[, ncol(P)]), C = Cs[ncol(P)],
       nnz = sum(fit$beta[, ncol(P)] != 0), nonconv = nonconv)
}
learner_l1se <- function(X, y, w, grp, Xnew, ...) {
  o <- l1_fit_score_1se(X, y, w, grp, Xnew)
  list(score = o$score, info = list(log10C = log10(o$C), nnz = o$nnz, nonconv = o$nonconv))
}

# ---------------------------------------------------------------- screening
# weighted standardized difference of the class means, class 1 minus class 0
wt_tstat <- function(F, y, w) {
  one <- y == 1
  st <- function(idx) {
    ww <- w[idx]; s <- sum(ww)
    Fi <- F[idx, , drop = FALSE]
    mu <- colSums(Fi * ww) / s
    v <- colSums(sweep(Fi, 2L, mu)^2 * ww) / s
    list(mu = mu, v = v, ne = s^2 / sum(ww^2))
  }
  a <- st(one); b <- st(!one)
  (a$mu - b$mu) / sqrt(a$v / a$ne + b$v / b$ne + 1e-12)
}

# top-m screening with m chosen by stratified 3-fold CV AUC (ties: smaller m)
screen_fit <- function(F, y, w, grp, Fnew, Ms = SCREEN_MS, nfold = 3L) {
  p <- ncol(F)
  Ms <- Ms[Ms <= p]
  if (!length(Ms)) Ms <- p
  fold <- strat_folds(grp, nfold)
  cv <- matrix(NA_real_, nfold, length(Ms))
  for (f in seq_len(nfold)) {
    trn <- fold != f; val <- !trn
    if (length(unique(y[trn])) < 2L || length(unique(y[val])) < 2L) next
    t <- wt_tstat(F[trn, , drop = FALSE], y[trn], w[trn])
    ord <- order(-abs(t))
    for (c in seq_along(Ms)) {
      sel <- ord[seq_len(Ms[c])]
      cv[f, c] <- auc_vec(as.numeric(F[val, sel, drop = FALSE] %*% t[sel]), y[val] == 1)
    }
  }
  mcv <- colMeans(cv, na.rm = TRUE)
  mcv[!is.finite(mcv)] <- -Inf
  m <- Ms[which.max(mcv)]                          # first maximum = smaller m
  t <- wt_tstat(F, y, w)
  sel <- order(-abs(t))[seq_len(m)]
  list(score = as.numeric(Fnew[, sel, drop = FALSE] %*% t[sel]), m = m)
}

learner_scrm <- function(X, y, w, grp, Xnew, ...) {
  o <- screen_fit(X, y, w, grp, Xnew)
  list(score = o$score, info = list(m = o$m))
}

learner_scr2 <- function(X, y, w, grp, Xnew, max_pairs = QUAD_MAX_PAIRS, ...) {
  d <- ncol(X)
  npair <- d * (d - 1L) / 2
  if (d < 2L) {
    pairs <- matrix(0L, 0L, 2L)
  } else if (npair <= max_pairs) {
    pairs <- which(upper.tri(diag(d)), arr.ind = TRUE)
  } else {                                         # same screening as learner_quad
    one <- y == 1L
    D <- abs(wcov(X[one, , drop = FALSE], w[one]) - wcov(X[!one, , drop = FALSE], w[!one]))
    ut <- which(upper.tri(D), arr.ind = TRUE)
    pairs <- ut[order(-D[ut])[seq_len(max_pairs)], , drop = FALSE]
  }
  pairs <- pairs[order(pairs[, 2L], pairs[, 1L]), , drop = FALSE]
  feat <- function(Z) cbind(Z^2, if (nrow(pairs)) Z[, pairs[, 1L], drop = FALSE] * Z[, pairs[, 2L], drop = FALSE])
  F <- feat(X); Fnew <- feat(Xnew)
  # features centered and scaled with the fitting data (weights w), so that the
  # t-weighted sum is a diagonal-LDA score on comparable scales
  ws <- w / sum(w)
  mu <- colSums(F * ws)
  sdv <- sqrt(colSums(sweep(F, 2L, mu)^2 * ws)); sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
  F <- sweep(sweep(F, 2L, mu), 2L, sdv, "/")
  Fnew <- sweep(sweep(Fnew, 2L, mu), 2L, sdv, "/")
  o <- screen_fit(F, y, w, grp, Fnew)
  list(score = o$score, info = list(m = o$m))
}

# ---------------------------------------------------------------- xscr (2.6)
# standardized AUC of the columns of S against ONE fixed labelling (pos)
std_auc_cols <- function(S, pos) {
  n <- nrow(S); n1 <- sum(pos); n0 <- n - n1
  R <- apply(S, 2L, rank, ties.method = "average")
  if (!is.matrix(R)) R <- matrix(R, ncol = ncol(S))
  auc <- (colSums(R[pos, , drop = FALSE]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  (auc - 0.5) / sqrt((n + 1) / (12 * n1 * n0))
}

# lighter than apply(M, 2, sd) / order(): this code runs inside big loops, and on
# this machine R occasionally dies in base functions (CLAUDE.md), the more the
# more temporary objects there are
colsd <- function(M) sqrt(pmax(colSums(M^2) - nrow(M) * colMeans(M)^2, 0) / (nrow(M) - 1))
top_m <- function(v, m) {                     # indices of the m largest, ties: first
  ix <- integer(m)
  for (j in seq_len(m)) { ix[j] <- which.max(v); v[ix[j]] <- -Inf }
  ix
}

xscr_features <- function(Z, pairs)
  cbind(Z, Z^2, if (nrow(pairs)) Z[, pairs[, 1L], drop = FALSE] * Z[, pairs[, 2L], drop = FALSE])

# k x (B + 1) matrix of cross-fitted screening statistics for top-m selection
xscr_T <- function(Zfit, gfit, Zte, Lall, k, m) {
  d <- ncol(Zfit)
  pairs <- if (d <= XSCR_MAX_D_CROSS && d >= 2L) which(upper.tri(diag(d)), arr.ind = TRUE)
           else matrix(0L, 0L, 2L)
  Ff <- xscr_features(Zfit, pairs); Ft <- xscr_features(Zte, pairs)
  nf <- nrow(Ff); M <- nrow(Ft); nb <- ncol(Lall)
  mf <- colMeans(Ff); sf <- colsd(Ff); sf[!is.finite(sf) | sf < 1e-12] <- 1
  mt <- colMeans(Ft); st <- colsd(Ft); st[!is.finite(st) | st < 1e-12] <- 1
  Ffz <- sweep(sweep(Ff, 2L, mf), 2L, sf, "/")        # fit half, own statistics
  Ftz <- sweep(sweep(Ft, 2L, mt), 2L, st, "/")        # test half, own statistics
  Ftf <- sweep(sweep(Ft, 2L, mf), 2L, sf, "/")        # test half, fit statistics
  Fft <- sweep(sweep(Ff, 2L, mt), 2L, st, "/")        # fit half, test statistics
  m <- min(m, ncol(Ff))
  tot_f <- colSums(Ffz); tot_t <- colSums(Ftz)
  S1 <- matrix(0, M, k)
  T2 <- matrix(0, k, nb)
  for (i in seq_len(k)) {
    # direction 1: selection on the fit half, true labels
    ii <- gfit == i; n1 <- sum(ii); n0 <- nf - n1
    s1 <- colSums(Ffz[ii, , drop = FALSE])
    t1 <- (s1 / n1 - (tot_f - s1) / n0) / sqrt(1 / n1 + 1 / n0)
    sel <- top_m(abs(t1), m)
    S1[, i] <- as.numeric(Ftf[, sel, drop = FALSE] %*% t1[sel])
    # direction 2: selection on the test half for EVERY labelling b
    I <- (Lall == i) * 1
    ni <- colSums(I)
    si <- crossprod(I, Ftz)                              # nb x p
    tt <- (si / ni - (rep(tot_t, each = nb) - si) / (M - ni)) / sqrt(1 / ni + 1 / (M - ni))
    S2 <- matrix(0, nf, nb)
    if (m == 1L) {
      j <- max.col(abs(tt), ties.method = "first")
      S2 <- sweep(Fft[, j, drop = FALSE], 2L, tt[cbind(seq_len(nb), j)], "*")
    } else {
      for (b in seq_len(nb)) {
        sel <- top_m(abs(tt[b, ]), m)
        S2[, b] <- Fft[, sel, drop = FALSE] %*% tt[b, sel]
      }
    }
    T2[i, ] <- std_auc_cols(S2, gfit == i)
  }
  T1 <- ovr_T(S1, seq_len(k), Lall)
  (T1 + T2) / sqrt(2)
}

# ---------------------------------------------------------------- g5t, g5f
g5_standardize <- function(W) {
  mu <- rowMeans(W)
  sdv <- apply(W, 1L, sd)
  T <- (W - mu) / ifelse(sdv > 0, sdv, 1)
  T[sdv == 0, ] <- 0
  T
}

g5_counts <- function(E, lab_mat, keep, k) {
  LA <- lab_mat[E[, 1L], , drop = FALSE]; LB <- lab_mat[E[, 2L], , drop = FALSE]
  same <- (LA == LB) * keep
  W <- matrix(0, k, ncol(lab_mat))
  for (i in seq_len(k)) W[i, ] <- colSums(same * (LA == i))
  W
}

g5t_T <- function(Zte, Lall, k) {
  E <- kmst_edges(sqrt(sqdist(Zte)), SC_KMST)
  g5_standardize(g5_counts(E, Lall, rep(1, nrow(E)), k))
}

g5f_T <- function(Zfit, gfit, Zte, Lall, k) {
  nf <- nrow(Zfit)
  Z <- rbind(Zfit, Zte)
  E <- kmst_edges(sqrt(sqdist(Z)), SC_KMST)
  lab <- rbind(matrix(gfit, nf, ncol(Lall)), Lall)
  keep <- as.numeric(E[, 1L] > nf | E[, 2L] > nf)
  g5_standardize(g5_counts(E, lab, keep, k))
}

component_T <- function(name, sc, Lall) {
  switch(name,
    xscr1 = xscr_T(sc$Zfit, sc$gfit, sc$Zte, Lall, sc$k, 1L),
    xscr3 = xscr_T(sc$Zfit, sc$gfit, sc$Zte, Lall, sc$k, 3L),
    g5t   = g5t_T(sc$Zte, Lall, sc$k),
    g5f   = g5f_T(sc$Zfit, sc$gfit, sc$Zte, Lall, sc$k),
    stop("unknown component: ", name))
}
