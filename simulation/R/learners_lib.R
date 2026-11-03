# =============================================================================
# Learner library for the one-vs-rest classifiers (Section 2.4, eq. (2.11))
#
# Common interface
#   learner_<name>(X, y, w, grp, Xnew, ...) -> list(score, info)
#     X, Xnew : fitting / new data, already standardized with the fitting data
#     y       : 0/1, 1 = population i
#     w       : observation weights of Section 2.4 (mean one)
#     grp     : population of each fitting observation (stratified CV)
#     score   : length nrow(Xnew), larger = more like population i
#     info    : named list of scalars (selected tuning values etc.)
# A learner may use only X, y, w, grp and the current R random-number stream,
# never the test labels; this is what Theorem 1 requires.
# =============================================================================

LIB_DEFAULT <- "l1,nn"
KNN_K_GRID <- c(5L, 10L, 20L, 40L, 80L)   # plus K_N = ceil(sqrt(N_fit)), see learner_knn
# Upper bound on the screened cross products of quad: 5000 in the pilots
# (libpilot, libpilot2, libweak, libsize), 1000 from the revised E1-E7 on
# (21 Oct 2026). Linear and squared terms are always kept; in the revised
# models the signal-bearing cross products number at most 190 (C', C2 at
# d = 20; D with s = ceil(100 beta) <= 20), so the cap does not bind on them.
QUAD_MAX_PAIRS <- 1000L

# folds stratified by population (same construction as in l1_fit_score)
strat_folds <- function(grp, nfold) {
  fold <- integer(length(grp))
  for (g in unique(grp)) {
    ix <- which(grp == g)
    fold[ix] <- sample(rep_len(seq_len(nfold), length(ix)))
  }
  fold
}

# ---------------------------------------------------------------- l1
learner_l1 <- function(X, y, w, grp, Xnew, ...) {
  o <- l1_fit_score(X, y, w, grp, Xnew)
  list(score = o$score, info = list(log10C = log10(o$C), nnz = o$nnz, nonconv = o$nonconv))
}

# ---------------------------------------------------------------- quad
# L1 logistic regression on [x_j, x_j^2, selected x_j x_l]. The log density
# ratio of two Gaussians is quadratic, so this is the right function class for
# variance and covariance changes (models B, D).
# The cross products are screened by |Delta Sigma_jl| (weighted class
# covariances) on ALL fitting data, outside the CV folds. Screening uses the
# fitting data only, so it affects tuning, not validity (Theorem 1). When
# d(d-1)/2 <= QUAD_MAX_PAIRS every pair is kept and no screening takes place.
wcov <- function(X, w) {
  w <- w / sum(w)
  mu <- colSums(X * w)
  Xc <- sweep(X, 2L, mu)
  crossprod(Xc * sqrt(w))
}

quad_features <- function(X, pairs) {
  cbind(X, X^2, if (nrow(pairs)) X[, pairs[, 1L], drop = FALSE] * X[, pairs[, 2L], drop = FALSE])
}

learner_quad <- function(X, y, w, grp, Xnew, max_pairs = QUAD_MAX_PAIRS, ...) {
  d <- ncol(X)
  npair <- d * (d - 1L) / 2
  if (d < 2L) {
    pairs <- matrix(0L, 0L, 2L)
  } else if (npair <= max_pairs) {
    pairs <- which(upper.tri(diag(d)), arr.ind = TRUE)
  } else {
    one <- y == 1L
    D <- abs(wcov(X[one, , drop = FALSE], w[one]) - wcov(X[!one, , drop = FALSE], w[!one]))
    ut <- which(upper.tri(D), arr.ind = TRUE)
    pairs <- ut[order(-D[ut])[seq_len(max_pairs)], , drop = FALSE]
  }
  pairs <- pairs[order(pairs[, 2L], pairs[, 1L]), , drop = FALSE]   # fixed column order
  F <- quad_features(X, pairs); Fnew <- quad_features(Xnew, pairs)
  st <- standardize_fit(F)
  o <- l1_fit_score(apply_std(F, st), y, w, grp, apply_std(Fnew, st))
  list(score = o$score, info = list(log10C = log10(o$C), nnz = o$nnz, nonconv = o$nonconv))
}

# ---------------------------------------------------------------- knn
# weighted kNN vote with rank kernel kappa_r = (K - r + 1) / K:
#   s(z) = sum_{r<=K} kappa_r w_(r) y_(r) / sum_{r<=K} kappa_r w_(r)
# Returns, for every K in Ks, the scores of the rows of Xnew.
knn_scores <- function(X, y, w, Xnew, Ks) {
  D2 <- outer(rowSums(Xnew^2), rowSums(X^2), "+") - 2 * tcrossprod(Xnew, X)
  Kmax <- max(Ks)
  # one sort per new point; columns of NB = the Kmax nearest fitting points
  NB <- t(apply(D2, 1L, function(v) order(v)[seq_len(Kmax)]))
  if (Kmax == 1L) NB <- matrix(NB, ncol = 1L)
  WN <- matrix(w[NB], nrow(NB)); YN <- matrix(y[NB], nrow(NB))
  vapply(Ks, function(K) {
    kap <- (K - seq_len(K) + 1) / K
    a <- WN[, seq_len(K), drop = FALSE] %*% kap
    b <- (WN[, seq_len(K), drop = FALSE] * YN[, seq_len(K), drop = FALSE]) %*% kap
    as.numeric(b / a)
  }, numeric(nrow(Xnew)))
}

# The grid contains K_N = ceil(sqrt(N_fit)), which grows with the sample size
# (K_N -> Inf, K_N / N -> 0), as the universal-consistency remark requires.
learner_knn <- function(X, y, w, grp, Xnew, Ks = KNN_K_GRID, nfold = 3L, ...) {
  Ks <- sort(unique(c(Ks, as.integer(ceiling(sqrt(nrow(X)))))))
  w <- w / mean(w)
  fold <- strat_folds(grp, nfold)
  ntr_min <- min(vapply(seq_len(nfold), function(f) sum(fold != f), integer(1)))
  Ks <- Ks[Ks <= ntr_min - 1L]
  if (!length(Ks)) Ks <- max(1L, ntr_min - 1L)
  cv <- matrix(NA_real_, nfold, length(Ks))
  for (f in seq_len(nfold)) {
    trn <- fold != f; val <- !trn
    if (length(unique(y[trn])) < 2L || length(unique(y[val])) < 2L) next
    S <- knn_scores(X[trn, , drop = FALSE], y[trn], w[trn], X[val, , drop = FALSE], Ks)
    if (!is.matrix(S)) S <- matrix(S, ncol = length(Ks))
    for (c in seq_along(Ks)) cv[f, c] <- auc_vec(S[, c], y[val] == 1)
  }
  mcv <- colMeans(cv, na.rm = TRUE)
  mcv[!is.finite(mcv)] <- -Inf
  K <- Ks[max(which(mcv == max(mcv)))]     # largest K among ties
  s <- knn_scores(X, y, w, Xnew, K)
  list(score = as.numeric(s), info = list(K = K))
}

# ---------------------------------------------------------------- rf
# probability forest (ranger), no tuning; the seed is drawn from the current
# R stream so the fit is reproducible from the replication seed.
learner_rf <- function(X, y, w, grp, Xnew, num.trees = 500L, ...) {
  seed <- sample.int(.Machine$integer.max, 1L)
  colnames(X) <- colnames(Xnew) <- paste0("x", seq_len(ncol(X)))
  fit <- ranger::ranger(x = X, y = factor(y, levels = c(0L, 1L)), probability = TRUE,
                        num.trees = num.trees, mtry = max(1L, floor(sqrt(ncol(X)))),
                        min.node.size = 5L, case.weights = w, num.threads = 1L,
                        seed = seed, verbose = FALSE)
  P <- predict(fit, data = Xnew, num.threads = 1L)$predictions
  list(score = as.numeric(P[, "1"]), info = list())
}

# ---------------------------------------------------------------- nn
learner_nn <- function(X, y, w, grp, Xnew, nn_cfg = NN_DEFAULT, ...) {
  fit <- nn_fit(X, y, w, grp, out = 1L, cfg = nn_cfg)
  list(score = nn_score(fit, Xnew), info = list(epochs = fit$epochs))
}

# ---------------------------------------------------------------- helpers
LEARNERS <- c("l1", "quad", "knn", "rf", "nn")

parse_lib <- function(lib) {
  if (is.null(lib) || length(lib) == 0L || is.na(lib) || !nzchar(lib)) lib <- LIB_DEFAULT
  v <- trimws(strsplit(lib, ",", fixed = TRUE)[[1L]])
  known <- c(LEARNERS, if (exists("LEARNERS2")) LEARNERS2, if (exists("COMPONENTS")) COMPONENTS)
  bad <- setdiff(v, known)
  if (length(bad)) stop("unknown learner(s): ", paste(bad, collapse = ", "))
  if (anyDuplicated(v)) stop("duplicated learner in lib: ", lib)
  v
}

# "l1+nn;l1+quad+knn+nn" -> list(c("l1","nn"), c("l1","quad","knn","nn"));
# default: one set equal to lib
parse_libsets <- function(libsets, lib) {
  if (is.null(libsets) || length(libsets) == 0L || is.na(libsets) || !nzchar(libsets))
    return(list(lib))
  sets <- lapply(strsplit(libsets, ";", fixed = TRUE)[[1L]],
                 function(s) trimws(strsplit(s, "+", fixed = TRUE)[[1L]]))
  for (s in sets) if (!all(s %in% lib))
    stop("libset '", paste(s, collapse = "+"), "' is not contained in lib = ", paste(lib, collapse = ","))
  sets
}
