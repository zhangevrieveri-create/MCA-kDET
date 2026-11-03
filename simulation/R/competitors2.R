# =============================================================================
# Further competitors (third review round, item G1; written 1 Nov 2026)
#
# method groups (column "groups" of a cell):
#   comp2     maxmean     max-type mean test, coordinatewise standardized
#                         differences of means, largest over coordinates and
#                         over pairs of populations (Cai, Liu and Xia 2014,
#                         statistic without the precision transform)
#             clx_mean    the same after the transform x -> Omega x with
#                         Omega = inverse of the label-free sample covariance
#                         and denominators omega_jj (Cai, Liu and Xia 2014);
#                         NA when N <= 2d
#             clx_cov     max-type covariance test (Cai, Liu and Xia 2013),
#                         largest over entries j <= l and over pairs of
#                         populations
#             schott_perm Schott's (2007) statistic, calibrated by permutation
#                         instead of its normal limit
#             mmd_max     largest pairwise unbiased MMD^2, same Gaussian kernel
#                         and bandwidth as mmd (instead of the average)
#             kmd         KMD of Huang and Sen (2024), discrete kernel, directed
#                         K-NN graph with K = ceiling(N/10), the default of
#                         KMD::KMD_test; Euclidean distances on the raw data
#   songchen  songchen    S of Song and Chen (2022), gTestsMulti::gtestsmulti,
#             songchen_A  and S^A, on the k-MST with k = SC_KMST, with the
#                         package's own permutation p-values (B permutations)
#
# All comp2 tests use ONE set of B uniform permutations of the pooled labels
# and the p-value (2.8). Every statistic is recomputed from the data and the
# permuted labels (variances, covariances and the KMD graph are label-free or
# recomputed per permutation), so each test has exact level under H0.
#
# The groups are evaluated inside competitor_tests() after "dist" and
# "schott", so adding them to a cell leaves every other column of the cell
# unchanged; "songchen" draws its permutations last, inside the package.
# Needs gTestsMulti for "songchen" (install it separately, see CLAUDE.md).
# =============================================================================

SC_KMST <- 5L   # k of the k-MST used by the Song-Chen tests (1 = plain MST)

# ---------------------------------------------------------------- means
# Z: N x d data (possibly transformed); G: indicator matrix of group_indicator();
# v: length-d variances under H0. Statistic for labelling b:
#   max_{i<j} max_l (zbar_il - zbar_jl)^2 / {v_l (1/n_i + 1/n_j)}
maxmean_stat <- function(Z, G, n, k, nb, v) {
  Sm <- crossprod(G, Z)                      # (k nb) x d, row (b - 1) k + i
  pr <- t(combn(k, 2L))
  vapply(seq_len(nb), function(b) {
    Mb <- Sm[(b - 1L) * k + seq_len(k), , drop = FALSE] / n
    max(vapply(seq_len(nrow(pr)), function(r) {
      i <- pr[r, 1L]; j <- pr[r, 2L]
      max((Mb[i, ] - Mb[j, ])^2 / v) / (1 / n[i] + 1 / n[j])
    }, numeric(1)))
  }, numeric(1))
}

# ---------------------------------------------------------------- covariances
# Schott's z statistic from k covariance matrices (divisor n_i - 1) and
# nn = n_i - 1; the same arithmetic as schott_test() in tests.R
schott_z <- function(S, nn) {
  g <- length(S)
  ahat <- function(Si, ni) ni^2 / ((ni - 1) * (ni + 2)) * (sum(Si^2) - sum(diag(Si))^2 / ni)
  a <- mapply(ahat, S, nn)
  t <- 0
  for (i in 1:(g - 1)) for (j in (i + 1):g) t <- t + a[i] + a[j] - 2 * sum(S[[i]] * S[[j]])
  npool <- sum(nn)
  Sp <- Reduce(`+`, Map(`*`, S, nn)) / npool
  a0 <- ahat(Sp, npool)
  pr <- t(combn(g, 2L))
  v <- 4 * a0^2 * (sum((1 / nn[pr[, 1L]] + 1 / nn[pr[, 2L]])^2) + (g - 1) * (g - 2) * sum(1 / nn^2))
  t / sqrt(v)
}

# For every labelling (column of Lall): the CLX (2013) statistic
#   max_{i<j} max_{l <= m} (s_i,lm - s_j,lm)^2 / (theta_i,lm / n_i + theta_j,lm / n_j),
#   s_i = covariance with divisor n_i, theta_i,lm = mean{(y_l y_m - s_i,lm)^2}
#   over the centred observations y of population i,
# and Schott's z statistic. Returns an nb x 2 matrix.
cov_perm_stats <- function(X, Lall, k) {
  nb <- ncol(Lall); d <- ncol(X)
  ut <- upper.tri(matrix(0, d, d), diag = TRUE)
  pr <- t(combn(k, 2L))
  res <- matrix(NA_real_, nb, 2L, dimnames = list(NULL, c("clx", "schott")))
  for (b in seq_len(nb)) {
    S <- vector("list", k); Th <- vector("list", k); ni <- integer(k)
    for (i in seq_len(k)) {
      Y <- X[Lall[, b] == i, , drop = FALSE]
      ni[i] <- nrow(Y)
      Y <- Y - rep(colMeans(Y), each = ni[i])
      S[[i]] <- crossprod(Y) / ni[i]
      Th[[i]] <- crossprod(Y * Y) / ni[i] - S[[i]]^2
    }
    res[b, 1L] <- max(vapply(seq_len(nrow(pr)), function(r) {
      i <- pr[r, 1L]; j <- pr[r, 2L]
      v <- (S[[i]] - S[[j]])^2 / (Th[[i]] / ni[i] + Th[[j]] / ni[j])
      max(v[ut], na.rm = TRUE)
    }, numeric(1)))
    res[b, 2L] <- schott_z(Map(function(Si, m) Si * m / (m - 1), S, ni), ni - 1)
  }
  res
}


# ---------------------------------------------------------------- one-vs-rest versions
# (third review round follow-up: in the scenarios S3 ONE population among
# k = 5, 10, 20 deviates, and the pairwise maximum over k(k-1)/2 pairs wastes
# power, so the max-type competitors also get a population-vs-rest version)
# maxmean_ovr: max_i max_l (zbar_il - zbar_rest,l)^2 / {v_l (1/n_i + 1/(N - n_i))}
maxmean_ovr_stat <- function(Z, G, n, k, nb, v) {
  Sm <- crossprod(G, Z); tot <- colSums(Z); N <- nrow(Z)
  vapply(seq_len(nb), function(b) {
    Sb <- Sm[(b - 1L) * k + seq_len(k), , drop = FALSE]
    Mi <- Sb / n
    Mr <- (rep(tot, each = k) - Sb) / (N - n)
    max(apply((Mi - Mr)^2 / rep(v, each = k), 1L, max) / (1 / n + 1 / (N - n)))
  }, numeric(1))
}

# Covariance (divisor n, centred at the mean of the set) S and
# Theta_lm = mean{(a_l a_m - S_lm)^2}, a = x - mean, from RAW moments that are
# additive over observations: m1 = sum x, m2 = sum x x', m21[l, m] = sum x_l^2 x_m,
# m22[l, m] = sum x_l^2 x_m^2. sum a_l^2 a_m^2 is expanded in these moments.
cov_from_moments <- function(n, m1, m2, m21, m22) {
  mu <- m1 / n
  q <- diag(m2)                                             # sum x_l^2
  A <- m22 - 2 * sweep(m21, 2L, mu, "*") - 2 * mu * t(m21) +
       outer(q, mu^2) + outer(mu^2, q) + 4 * outer(mu, mu) * m2 -
       3 * n * outer(mu^2, mu^2)
  S <- m2 / n - tcrossprod(mu)
  list(S = S, Theta = A / n - S^2)
}
raw_moments <- function(Y) {
  Y2 <- Y * Y
  list(n = nrow(Y), m1 = colSums(Y), m2 = crossprod(Y), m21 = crossprod(Y2, Y), m22 = crossprod(Y2))
}

# clx_cov_ovr: max_i max_{l <= m} (S_i,lm - S_rest,lm)^2 /
#              (Theta_i,lm / n_i + Theta_rest,lm / n_rest), every labelling
clx_cov_ovr_stat <- function(X, Lall, k) {
  d <- ncol(X); ut <- upper.tri(matrix(0, d, d), diag = TRUE)
  tot <- raw_moments(X); N <- nrow(X)
  vapply(seq_len(ncol(Lall)), function(b) {
    best <- -Inf
    for (i in seq_len(k)) {
      rm <- raw_moments(X[Lall[, b] == i, , drop = FALSE])
      ci <- cov_from_moments(rm$n, rm$m1, rm$m2, rm$m21, rm$m22)
      nr <- N - rm$n
      cr <- cov_from_moments(nr, tot$m1 - rm$m1, tot$m2 - rm$m2, tot$m21 - rm$m21, tot$m22 - rm$m22)
      v <- (ci$S - cr$S)^2 / (ci$Theta / rm$n + cr$Theta / nr)
      best <- max(best, max(v[ut], na.rm = TRUE))
    }
    best
  }, numeric(1))
}

# ---------------------------------------------------------------- kernels, graphs
# largest pairwise unbiased MMD^2 (mmd_stat in tests.R takes the average)
mmd_max_stat <- function(K, G, n, k, nb) {
  KG <- K %*% G
  pairs <- t(combn(k, 2L))
  vapply(seq_len(nb), function(b) {
    cols <- (b - 1L) * k + seq_len(k)
    S <- crossprod(G[, cols], KG[, cols])
    w <- (diag(S) - n) / (n * (n - 1))
    max(w[pairs[, 1L]] + w[pairs[, 2L]] - 2 * S[pairs] / (n[pairs[, 1L]] * n[pairs[, 2L]]))
  }, numeric(1))
}

# KMD with the discrete kernel on the directed K-NN graph (Huang and Sen 2024):
#   eta = {(N K)^{-1} sum_i sum_{j in N_K(i)} 1(Y_i = Y_j) - p0} / (1 - p0),
#   p0 = sum_c n_c (n_c - 1) / {N (N - 1)}.
# The label counts do not change under permutation, so p0 is a constant.
# D2: squared distances. Returns eta for every column of Lall.
kmd_stat <- function(D2, Lall, K) {
  N <- nrow(D2)
  Dd <- D2; diag(Dd) <- Inf
  NN <- do.call(rbind, lapply(seq_len(N), function(i) order(Dd[i, ])[seq_len(K)]))   # N x K
  ii <- rep(seq_len(N), K)                                    # column-major match of NN
  jj <- as.vector(NN)
  same <- vapply(seq_len(ncol(Lall)), function(b) {
    Lb <- Lall[, b]
    sum(Lb[ii] == Lb[jj])
  }, numeric(1))
  nc <- tabulate(Lall[, 1L])
  p0 <- sum(nc * (nc - 1)) / (N * (N - 1))
  (same / (N * K) - p0) / (1 - p0)
}

# union of kk successive minimum spanning trees (each one on the complete graph
# without the edges of the earlier ones); D: distances
kmst_edges <- function(D, kk) {
  E <- NULL; Dw <- D
  for (t in seq_len(kk)) {
    e <- mst_edges(Dw)
    E <- rbind(E, e)
    Dw[e] <- Inf; Dw[e[, 2:1, drop = FALSE]] <- Inf
  }
  E
}

# ---------------------------------------------------------------- driver
comp2_tests <- function(data, B, groups) {
  k <- length(data)
  n <- vapply(data, nrow, integer(1))
  X <- do.call(rbind, data)
  N <- nrow(X); d <- ncol(X)
  L <- rep(seq_len(k), n)
  out <- list(); tm <- list()
  D2 <- sqdist(X)
  Lall <- NULL
  if ("comp2" %in% groups) {
    Lall <- perm_labels(L, B); nb <- ncol(Lall)
    G <- group_indicator(Lall, k)
    tm$maxmean <- system.time({
      out$maxmean <- pval_perm(maxmean_stat(X, G, n, k, nb, apply(X, 2L, var)))
    })[["elapsed"]]
    tm$clx_mean <- system.time({
      out$clx_mean <- if (N > 2L * d) {
        Om <- solve(cov(X))
        pval_perm(maxmean_stat(X %*% Om, G, n, k, nb, diag(Om)))
      } else NA_real_
    })[["elapsed"]]
    tm$clx_cov <- system.time(cs <- cov_perm_stats(X, Lall, k))[["elapsed"]]
    out$clx_cov <- pval_perm(cs[, "clx"])
    out$schott_perm <- pval_perm(cs[, "schott"])
    tm$mmd_max <- system.time({
      h2 <- median(D2[upper.tri(D2)])
      out$mmd_max <- pval_perm(mmd_max_stat(exp(-D2 / (2 * h2)), G, n, k, nb))
    })[["elapsed"]]
    tm$kmd <- system.time({
      eta <- kmd_stat(D2, Lall, ceiling(N / 10))
      out$kmd <- pval_perm(eta)
    })[["elapsed"]]
    out$kmd_eta <- eta[1L]
  }
  if ("comp2ovr" %in% groups) {
    if (is.null(Lall)) Lall <- perm_labels(L, B)
    nb <- ncol(Lall)
    G <- group_indicator(Lall, k)
    tm$maxmean_ovr <- system.time(
      out$maxmean_ovr <- pval_perm(maxmean_ovr_stat(X, G, n, k, nb, apply(X, 2L, var))))[["elapsed"]]
    tm$clx_cov_ovr <- system.time(
      out$clx_cov_ovr <- pval_perm(clx_cov_ovr_stat(X, Lall, k)))[["elapsed"]]
  }
  if ("songchen" %in% groups) {
    tm$songchen <- system.time({
      E <- kmst_edges(sqrt(D2), SC_KMST)
      r <- gTestsMulti::gtestsmulti(E, data, perm = B)
      pv <- unlist(r$pval)
      # the package returns #{perm >= obs}/B, which can be exactly 0 and is not
      # a valid p-value; use (#{...} + 1)/(B + 1) as everywhere else
      out$songchen <- unname((round(pv["S_perm"] * B) + 1) / (B + 1))
      out$songchen_A <- unname((round(pv["S_A_perm"] * B) + 1) / (B + 1))
    })[["elapsed"]]
  }
  list(p = out, tm = tm)
}
