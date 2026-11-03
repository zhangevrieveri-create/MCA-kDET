# =============================================================================
# Learners of Section 2.4
#   * feedforward ReLU network (2.9), q = (d, 64, 64, out), Glorot-uniform
#     init, Adam (lr 1e-3, batch 64, <= 200 epochs), first-layer spectral
#     projection W0 <- W0 / max(1, ||W0||_2) after every epoch (Gouk et al.),
#     weighted cross-entropy, early stopping on weighted validation CE
#     (20% of each population's fitting observations, patience 10,
#     min decrease 1e-6), weights of the best epoch restored.
#   * L1-penalized logistic regression, inverse penalty C in
#     10^{-4, -2.4, -0.8, 0.8, 2.4, 4}, 3-fold CV stratified by population,
#     largest mean AUC, ties -> strongest penalty, refit on all fitting data.
#     Implemented with glmnet: the LIBLINEAR objective ||w||_1 + C sum_i w_i l_i
#     equals (C n) * [ (1/n) sum_i w_i l_i + lambda ||w||_1 ] with
#     lambda = 1 / (C n) when the weights have mean one (glmnet's scaling);
#     glmnet leaves the intercept unpenalized.
# Written in base R so that the default network is reproduced exactly.
# =============================================================================

NN_DEFAULT <- list(width = c(64L, 64L), lr = 1e-3, batch = 64L, epochs = 200L,
                   patience = 10L, min_delta = 1e-6, val_frac = 0.2,
                   beta1 = 0.9, beta2 = 0.999, eps = 1e-8, spectral = TRUE)

L1_C_GRID <- 10^c(-4, -2.4, -0.8, 0.8, 2.4, 4)

# ---------------------------------------------------------------- utilities
standardize_fit <- function(Xfit) {
  mu <- colMeans(Xfit)
  sdv <- sqrt(colSums(sweep(Xfit, 2L, mu)^2) / max(1, nrow(Xfit) - 1))
  sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
  list(mu = mu, sd = sdv)
}
apply_std <- function(X, st) sweep(sweep(X, 2L, st$mu, "-"), 2L, st$sd, "/")

softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))

# AUC with ties counted as one half (Mann-Whitney)
auc_vec <- function(score, pos) {
  n1 <- sum(pos); n0 <- length(pos) - n1
  r <- rank(score, ties.method = "average")
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# validation subset: val_frac of the fitting observations of each population
val_split <- function(grp, frac) {
  v <- logical(length(grp))
  for (g in unique(grp)) {
    ix <- which(grp == g)
    nv <- max(1L, round(frac * length(ix)))
    if (length(ix) >= 2L) v[ix[sample.int(length(ix), nv)]] <- TRUE
  }
  v
}

# ---------------------------------------------------------------- network
# X: n x d (already standardized); y: 0/1 vector (out = 1) or integer
# labels 1..K (out = K, softmax); w: observation weights; grp: population of
# each observation (for the stratified validation split).
nn_fit <- function(X, y, w, grp, out = 1L, cfg = NN_DEFAULT) {
  n <- nrow(X); d <- ncol(X)
  dims <- c(d, cfg$width, out)
  L <- length(dims) - 1L
  W <- vector("list", L); b <- vector("list", L)
  for (l in seq_len(L)) {
    a <- sqrt(6 / (dims[l] + dims[l + 1L]))
    W[[l]] <- matrix(runif(dims[l] * dims[l + 1L], -a, a), dims[l], dims[l + 1L])
    b[[l]] <- numeric(dims[l + 1L])
  }
  mW <- lapply(W, function(x) x * 0); vW <- mW
  mb <- lapply(b, function(x) x * 0); vb <- mb
  t_adam <- 0L

  w <- w / mean(w)
  isval <- val_split(grp, cfg$val_frac)
  tr <- which(!isval); va <- which(isval)
  if (out > 1L) {
    Y <- matrix(0, n, out); Y[cbind(seq_len(n), y)] <- 1
  } else {
    Y <- matrix(as.numeric(y), n, 1L)
  }

  forward <- function(Xb, W, b) {
    A <- vector("list", L)        # post-activation of hidden layers
    H <- Xb
    for (l in seq_len(L - 1L)) {
      Z <- H %*% W[[l]]; Z <- Z + rep(b[[l]], each = nrow(Z))
      H <- Z * (Z > 0)            # ReLU
      A[[l]] <- H
    }
    Z <- H %*% W[[L]]; Z <- Z + rep(b[[L]], each = nrow(Z))
    list(A = A, Z = Z)
  }
  loss_of <- function(Z, Yb, wb) {
    if (out == 1L) {
      li <- softplus(Z[, 1L]) - Yb[, 1L] * Z[, 1L]
    } else {
      mx <- apply(Z, 1L, max)
      lse <- mx + log(rowSums(exp(Z - mx)))
      li <- lse - rowSums(Z * Yb)
    }
    sum(wb * li) / sum(wb)
  }
  # gradient of the dZ wrt output: w * (p - y) / batch size
  dout <- function(Z, Yb, wb) {
    if (out == 1L) {
      P <- 1 / (1 + exp(-Z))
    } else {
      P <- exp(Z - apply(Z, 1L, max)); P <- P / rowSums(P)
    }
    (P - Yb) * (wb / length(wb))
  }

  best <- Inf; bestW <- W; bestb <- b; wait <- 0L; ep_done <- 0L
  Xva <- X[va, , drop = FALSE]; Yva <- Y[va, , drop = FALSE]; wva <- w[va]
  for (ep in seq_len(cfg$epochs)) {
    ord <- tr[sample.int(length(tr))]
    nb <- ceiling(length(ord) / cfg$batch)
    for (bi in seq_len(nb)) {
      ix <- ord[((bi - 1L) * cfg$batch + 1L):min(bi * cfg$batch, length(ord))]
      Xb <- X[ix, , drop = FALSE]
      fw <- forward(Xb, W, b)
      G <- dout(fw$Z, Y[ix, , drop = FALSE], w[ix])
      gW <- vector("list", L); gb <- vector("list", L)
      for (l in L:1L) {
        Hin <- if (l == 1L) Xb else fw$A[[l - 1L]]
        gW[[l]] <- crossprod(Hin, G)
        gb[[l]] <- colSums(G)
        if (l > 1L) G <- tcrossprod(G, W[[l]]) * (fw$A[[l - 1L]] > 0)
      }
      t_adam <- t_adam + 1L
      bc1 <- 1 - cfg$beta1^t_adam; bc2 <- 1 - cfg$beta2^t_adam
      for (l in seq_len(L)) {
        mW[[l]] <- cfg$beta1 * mW[[l]] + (1 - cfg$beta1) * gW[[l]]
        vW[[l]] <- cfg$beta2 * vW[[l]] + (1 - cfg$beta2) * gW[[l]]^2
        W[[l]] <- W[[l]] - cfg$lr * (mW[[l]] / bc1) / (sqrt(vW[[l]] / bc2) + cfg$eps)
        mb[[l]] <- cfg$beta1 * mb[[l]] + (1 - cfg$beta1) * gb[[l]]
        vb[[l]] <- cfg$beta2 * vb[[l]] + (1 - cfg$beta2) * gb[[l]]^2
        b[[l]] <- b[[l]] - cfg$lr * (mb[[l]] / bc1) / (sqrt(vb[[l]] / bc2) + cfg$eps)
      }
    }
    if (cfg$spectral) {
      s1 <- svd(W[[1L]], nu = 0L, nv = 0L)$d[1L]
      if (s1 > 1) W[[1L]] <- W[[1L]] / s1
    }
    ep_done <- ep
    vl <- loss_of(forward(Xva, W, b)$Z, Yva, wva)
    if (vl < best - cfg$min_delta) {
      best <- vl; bestW <- W; bestb <- b; wait <- 0L
    } else {
      wait <- wait + 1L
      if (wait >= cfg$patience) break
    }
  }
  structure(list(W = bestW, b = bestb, out = out, epochs = ep_done,
                 val_loss = best, forward = forward), class = "mca_nn")
}

# binary: returns the logit (monotone in the sigmoid, avoids saturation ties)
# multiclass: returns the matrix of log-softmax outputs (monotone in h_i)
nn_score <- function(fit, X) {
  Z <- fit$forward(X, fit$W, fit$b)$Z
  if (fit$out == 1L) return(Z[, 1L])
  mx <- apply(Z, 1L, max)
  Z - (mx + log(rowSums(exp(Z - mx))))
}

# ---------------------------------------------------------------- L1 logistic
# glmnet fits the C path from the strongest penalty down, with warm starts.
# When the fit at some C does not converge (nearly separable data, weak
# penalty), glmnet returns the path up to the last converged C and warns; the
# CV then ignores the missing values of C and the final fit falls back to the
# last converged C. maxit was 1e5 up to 22 Oct 2026; each non-converged fit
# then ran 1e5 passes (about 10 s at k = 10, d = 300), which made E1 slow.
# It is now 1e4; non-convergence is counted (info: nonconv), other warnings
# are muffled as before.
L1_MAXIT <- 1e4
glmnet_counted <- function(...) {
  nc <- 0L
  fit <- withCallingHandlers(glmnet::glmnet(...), warning = function(w) {
    if (grepl("onvergence", conditionMessage(w))) nc <<- nc + 1L
    invokeRestart("muffleWarning")
  })
  list(fit = fit, nonconv = nc)
}

# returns the linear predictor on Xnew
l1_fit_score <- function(X, y, w, grp, Xnew, Cs = L1_C_GRID, nfold = 3L) {
  nonconv <- 0L
  w <- w / mean(w)
  # folds stratified by population
  fold <- integer(length(y))
  for (g in unique(grp)) {
    ix <- which(grp == g)
    fold[ix] <- sample(rep_len(seq_len(nfold), length(ix)))
  }
  ordC <- order(Cs)                       # increasing C = decreasing lambda
  Cs <- Cs[ordC]
  cv <- matrix(NA_real_, nfold, length(Cs))
  for (f in seq_len(nfold)) {
    trn <- fold != f; val <- !trn
    if (length(unique(y[trn])) < 2L || length(unique(y[val])) < 2L) next
    wt <- w[trn] / mean(w[trn])
    lam <- 1 / (Cs * sum(trn))
    g <- glmnet_counted(X[trn, , drop = FALSE], y[trn],
              family = "binomial", weights = wt, lambda = lam,
              standardize = FALSE, intercept = TRUE, maxit = L1_MAXIT)
    fit <- g$fit; nonconv <- nonconv + g$nonconv
    P <- predict(fit, X[val, , drop = FALSE], type = "link")
    # glmnet may stop the path early; missing columns get NA
    for (c in seq_len(min(ncol(P), length(Cs)))) cv[f, c] <- auc_vec(P[, c], y[val] == 1)
  }
  mcv <- colMeans(cv, na.rm = TRUE)
  mcv[!is.finite(mcv)] <- -Inf
  cbest <- which(mcv == max(mcv))[1L]     # smallest C (strongest penalty) among ties
  lam_all <- 1 / (Cs * length(y))
  g <- glmnet_counted(X, y, family = "binomial", weights = w,
            lambda = lam_all[seq_len(cbest)], standardize = FALSE,
            intercept = TRUE, maxit = L1_MAXIT)
  fit <- g$fit; nonconv <- nonconv + g$nonconv
  P <- predict(fit, Xnew, type = "link")
  # C actually used: the selected one, or the last converged one before it
  list(score = as.numeric(P[, ncol(P)]), C = Cs[ncol(P)],
       nnz = sum(fit$beta[, ncol(P)] != 0),   # nonzero slopes (diagnostic only)
       nonconv = nonconv)                     # non-converged glmnet fits (diagnostic)
}
