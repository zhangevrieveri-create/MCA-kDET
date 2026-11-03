# =============================================================================
# MCA-kDET (Algorithm 1) and its variants, plus the competitors of Section 4
# All permutation tests use the p-value (2.8): (1 + #{T_b >= T_0}) / (B + 1).
# =============================================================================

pval_perm <- function(v) {            # v[1] observed, v[-1] permuted
  (1 + sum(v[-1L] >= v[1L] - 1e-10)) / length(v)
}

# combined max/sum p-value (2.12)
pval_comb <- function(vmax, vsum) {
  nb <- length(vmax)
  rmax <- vapply(vmax, function(x) sum(vmax >= x - 1e-10), numeric(1)) / nb
  rsum <- vapply(vsum, function(x) sum(vsum >= x - 1e-10), numeric(1)) / nb
  C <- pmin(rmax, rsum)
  (1 + sum(C[-1L] <= C[1L] + 1e-10)) / nb
}

# Step 1: within-population split
split_data <- function(data, tau) {
  lapply(data, function(X) {
    n <- nrow(X); f <- floor(tau * n)
    if (f < 1L || f > n - 1L) stop("infeasible split")
    idx <- sample.int(n)
    list(fit = X[idx[seq_len(f)], , drop = FALSE],
         test = X[idx[(f + 1L):n], , drop = FALSE])
  })
}

# labels: observed labels in column 1, B uniform permutations in columns 2..B+1
perm_labels <- function(L, B) {
  M <- length(L)
  cbind(L, vapply(seq_len(B), function(b) L[sample.int(M)], integer(M)))
}

# Standardized one-vs-rest AUCs (2.5)-(2.6) for every labelling.
#   S   : M x C matrix of scores on the pooled test set
#   pop : population that column c is oriented toward
#   Lall: M x (B+1) labels
# returns C x (B+1) matrix of T_{i,l}
ovr_T <- function(S, pop, Lall) {
  M <- nrow(S)
  R <- apply(S, 2L, rank, ties.method = "average")
  if (!is.matrix(R)) R <- matrix(R, ncol = 1L)
  Tm <- matrix(NA_real_, ncol(S), ncol(Lall))
  for (i in unique(pop)) {
    cols <- which(pop == i)
    I <- (Lall == i) * 1
    mi <- sum(Lall[, 1L] == i)
    RS <- crossprod(R[, cols, drop = FALSE], I)
    auc <- (RS - mi * (mi + 1) / 2) / (mi * (M - mi))
    sig <- sqrt((M + 1) / (12 * mi * (M - mi)))
    Tm[cols, ] <- (auc - 0.5) / sig
  }
  Tm
}

# pairwise standardized AUCs: list of score vectors with (i, j) orientation
pair_T <- function(Slist, pairs, Lall) {
  nb <- ncol(Lall)
  Tm <- matrix(NA_real_, length(Slist), nb)
  for (c in seq_along(Slist)) {
    i <- pairs[c, 1L]; j <- pairs[c, 2L]; s <- Slist[[c]]
    for (b in seq_len(nb)) {
      lb <- Lall[, b]
      sel <- which(lb == i | lb == j)
      r <- rank(s[sel], ties.method = "average")
      pos <- lb[sel] == i
      n1 <- sum(pos); n0 <- length(sel) - n1
      auc <- (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
      Tm[c, b] <- (auc - 0.5) / sqrt((n1 + n0 + 1) / (12 * n1 * n0))
    }
  }
  Tm
}

# ---------------------------------------------------------------- weights
# one-vs-rest weights of Section 2.4: class i total 1/2, population j != i
# weight proportional to omega_j / f_j, classes equal total weight
ovr_weights <- function(grp, i, omega, f) {
  w <- numeric(length(grp))
  w[grp == i] <- 0.5 / f[i]
  for (j in setdiff(seq_along(f), i))
    w[grp == j] <- 0.5 * (omega[j] / (1 - omega[i])) / f[j]
  w / mean(w)
}
shared_weights <- function(grp, omega, f) {
  w <- (omega / f)[grp]
  w / mean(w)
}

# ---------------------------------------------------------------- classifier scores
# Fits every classifier needed by the requested method groups on one split.
# Returns score matrices on the pooled test set.
#   lib: one-vs-rest learners, fitted for each i in this order (for each i,
#        the learners in the order of lib), so lib = c("l1", "nn") consumes
#        the random-number stream exactly as the original two-learner code.
fit_scores <- function(sp, groups, nn_cfg = NN_DEFAULT, lib = c("l1", "nn")) {
  k <- length(sp)
  f <- vapply(sp, function(s) nrow(s$fit), integer(1))
  m <- vapply(sp, function(s) nrow(s$test), integer(1))
  M <- sum(m); omega <- m / M
  Xfit <- do.call(rbind, lapply(sp, `[[`, "fit"))
  gfit <- rep(seq_len(k), f)
  Xte <- do.call(rbind, lapply(sp, `[[`, "test"))
  L <- rep(seq_len(k), m)
  st <- standardize_fit(Xfit)
  Zfit <- apply_std(Xfit, st); Zte <- apply_std(Xte, st)
  # components of R/learners_lib2.R (xscr1, xscr3, g5t, g5f) are not
  # one-vs-rest classifiers: they are computed in classifier_tests(). The
  # learners keep their order, hence their random-number consumption.
  comp <- if (exists("COMPONENTS")) intersect(lib, COMPONENTS) else character(0)
  lib <- setdiff(lib, comp)
  res <- list(L = L, m = m, k = k, tm = list(), comp = comp, Zfit = Zfit, gfit = gfit, Zte = Zte)

  if ("ovr" %in% groups) {
    S <- setNames(lapply(lib, function(a) matrix(NA_real_, M, k)), lib)
    info <- setNames(lapply(lib, function(a) list()), lib)
    tl <- setNames(numeric(length(lib)), lib)
    for (i in seq_len(k)) {
      # k = 2 (Remark 4 of the paper): the comparison of population 2 is that of
      # population 1 with the class labels reversed, so its score is the negative
      # score of population 1 and T_2 = T_1, ties included. Only one classifier
      # per learner is fitted. For k >= 3 nothing changes, including the use of
      # the random-number stream.
      if (k == 2L && i == 2L) {
        for (a in lib) S[[a]][, 2L] <- -S[[a]][, 1L]
        next
      }
      y <- as.integer(gfit == i)
      w <- ovr_weights(gfit, i, omega, f)
      for (a in lib) {
        fn <- get(paste0("learner_", a), mode = "function")
        tl[[a]] <- tl[[a]] + system.time({
          o <- fn(Zfit, y, w, gfit, Zte, nn_cfg = nn_cfg)
        })[["elapsed"]]
        S[[a]][, i] <- o$score
        for (nm in names(o$info)) info[[a]][[nm]] <- c(info[[a]][[nm]], o$info[[nm]])
      }
    }
    res$Sovr <- S; res$info_ovr <- info; res$lib <- lib
    for (a in lib) res$tm[[paste0("ovr_", a)]] <- tl[[a]]
  }
  if ("shared" %in% groups) {
    w <- shared_weights(gfit, omega, f)
    res$tm$shared_nn <- system.time({
      nn <- nn_fit(Zfit, gfit, w, gfit, out = k, cfg = nn_cfg)
      LS <- nn_score(nn, Zte)
    })[["elapsed"]]
    res$SSH <- LS
    res$yhat <- max.col(LS, ties.method = "first")
    res$ep_shared <- nn$epochs
  }
  # pairwise variant: the same library as the one-vs-rest classifiers; for
  # each pair the learners are fitted in the order of lib, so lib = c("l1",
  # "nn") consumes the random-number stream exactly as the original code
  if ("pair" %in% groups) {
    pairs <- t(combn(k, 2L))
    P <- setNames(lapply(lib, function(a) vector("list", nrow(pairs))), lib)
    tl <- setNames(numeric(length(lib)), lib)
    for (r in seq_len(nrow(pairs))) {
      i <- pairs[r, 1L]; j <- pairs[r, 2L]
      Xp <- rbind(sp[[i]]$fit, sp[[j]]$fit)
      gp <- rep(c(i, j), c(f[i], f[j]))
      y <- as.integer(gp == i)
      w <- ifelse(gp == i, 0.5 / f[i], 0.5 / f[j]); w <- w / mean(w)
      stp <- standardize_fit(Xp)
      Zp <- apply_std(Xp, stp); Zt <- apply_std(Xte, stp)
      for (a in lib) {
        fn <- get(paste0("learner_", a), mode = "function")
        tl[[a]] <- tl[[a]] + system.time({
          P[[a]][[r]] <- fn(Zp, y, w, gp, Zt, nn_cfg = nn_cfg)$score
        })[["elapsed"]]
      }
    }
    res$pairs <- pairs; res$Ppair <- P
    for (a in lib) res$tm[[paste0("pair_", a)]] <- tl[[a]]
  }
  res
}

# ---------------------------------------------------------------- classifier tests
# Computes every classifier-based p-value from one split with ONE set of
# permutations (Step 6), shared by all statistics.
#   libsets: list of character vectors (subsets of sc$lib); each gives the
#            p-value mca_set_<a+b+...> of the library restricted to that set.
classifier_tests <- function(sc, B, libsets = NULL) {
  k <- sc$k; L <- sc$L
  t0 <- proc.time()[["elapsed"]]
  Lall <- perm_labels(L, B)
  out <- list(); diag <- list()
  T1 <- NULL
  if (!is.null(sc$Sovr)) {
    lib <- sc$lib
    Tl <- lapply(sc$Sovr, function(S) ovr_T(S, seq_len(k), Lall))   # same Lall for all
    # components of R/learners_lib2.R, same permutations, same maximum
    for (cn in sc$comp) {
      t0c <- proc.time()[["elapsed"]]
      Tl[[cn]] <- component_T(cn, sc, Lall)
      diag[[paste0("t_", cn)]] <- proc.time()[["elapsed"]] - t0c
    }
    lib_all <- names(Tl)
    Tlib <- do.call(rbind, Tl)
    vmax <- apply(Tlib, 2L, max)
    vsum <- colSums(Tlib)
    out$mca      <- pval_perm(vmax)                       # MCA-kDET (whole lib)
    for (a in lib_all) out[[paste0("mca_", a)]] <- pval_perm(apply(Tl[[a]], 2L, max))
    out$mca_2s   <- pval_perm(apply(abs(Tlib), 2L, max))  # two-sided
    out$mca_sum  <- pval_perm(vsum)
    out$mca_comb <- pval_comb(vmax, vsum)
    for (s in libsets)
      out[[paste0("mca_set_", paste(s, collapse = "+"))]] <-
        pval_perm(apply(do.call(rbind, Tl[s]), 2L, max))
    # G2 (third review round, 1 Nov 2026): multiclass accuracy from the SAME
    # one-vs-rest fits and the SAME permutations as MCA-kDET. A test point is
    # assigned to the population whose classifier of learner a gives the
    # largest score (the k scores of one learner are on one scale: logits,
    # vote shares or probabilities); macc_ovr takes the largest accuracy over
    # the learners. Ties go to the first population, so no random numbers
    # are used and every other column is unchanged.
    acc_l <- lapply(sc$Sovr, function(S) colMeans(Lall == max.col(S, ties.method = "first")))
    for (a in lib) out[[paste0("macc_", a)]] <- pval_perm(acc_l[[a]])
    if (length(acc_l)) out$macc_ovr <- pval_perm(do.call(pmax, unname(acc_l)))
    for (a in lib) diag[[paste0("acc_", a, "_obs")]] <- acc_l[[a]][1L]
    for (a in lib) diag[[paste0("T_", a, "_dev")]] <- Tl[[a]][k, 1L]
    for (a in lib) diag[[paste0("auc_", a, "_dev")]] <- auc_vec(sc$Sovr[[a]][, k], L == k)
    diag$argmax <- which.max(Tlib[, 1L])
    # scalar tuning information, averaged over the k binary problems
    for (a in lib) for (nm in names(sc$info_ovr[[a]]))
      diag[[paste0(a, "_", nm)]] <- mean(sc$info_ovr[[a]][[nm]])
    if ("nn" %in% lib) diag$ep_ovr_mean <- mean(sc$info_ovr$nn$epochs)
    if ("l1" %in% lib) T1 <- Tl$l1
  }
  if (!is.null(sc$SSH)) {
    TS <- ovr_T(sc$SSH, seq_len(k), Lall)
    out$shared_nn <- pval_perm(apply(TS, 2L, max))
    if (!is.null(T1)) out$shared_lib <- pval_perm(apply(rbind(T1, TS), 2L, max))
    acc <- colMeans(Lall == sc$yhat)
    out$macc <- pval_perm(acc)                            # multiclass accuracy
    diag$acc_obs <- acc[1L]; diag$T_sh_dev <- TS[k, 1L]; diag$ep_shared <- sc$ep_shared
  }
  if (!is.null(sc$Ppair)) {
    Pl <- lapply(sc$Ppair, function(S) pair_T(S, sc$pairs, Lall))   # same Lall
    out$pair <- pval_perm(apply(do.call(rbind, Pl), 2L, max))         # whole lib
    for (a in names(Pl)) out[[paste0("pair_", a)]] <- pval_perm(apply(Pl[[a]], 2L, max))
  }
  diag$t_perm <- proc.time()[["elapsed"]] - t0
  list(p = out, diag = diag)
}

# ---------------------------------------------------------------- competitors
sqdist <- function(X) {
  r <- rowSums(X^2)
  D2 <- outer(r, r, "+") - 2 * tcrossprod(X)
  D2[D2 < 0] <- 0
  diag(D2) <- 0
  D2
}

group_indicator <- function(Lall, k) {
  # N x (k * nb) matrix, columns ordered (b=1: groups 1..k), (b=2: ...), ...
  N <- nrow(Lall); nb <- ncol(Lall)
  G <- matrix(0, N, k * nb)
  for (b in seq_len(nb)) G[cbind(seq_len(N), (b - 1L) * k + Lall[, b])] <- 1
  G
}

# DISCO between-sample statistic (Rizzo & Szekely 2010):
#   S = T/(2N) - (1/2) sum_i S_ii / n_i, S_ii = sum of distances within sample i,
#   T = total sum (permutation invariant); the DISCO F ratio is monotone in S.
disco_stat <- function(D, G, n, k, nb) {
  DG <- D %*% G
  Sii <- matrix(colSums(G * DG), k, nb)
  -colSums(Sii / n)
}

# average pairwise unbiased MMD^2, Gaussian kernel exp(-|x-y|^2 / (2 h^2)),
# h = median pairwise distance
mmd_stat <- function(K, G, n, k, nb) {
  KG <- K %*% G
  pairs <- t(combn(k, 2L))
  vapply(seq_len(nb), function(b) {
    cols <- (b - 1L) * k + seq_len(k)
    S <- crossprod(G[, cols], KG[, cols])
    w <- (diag(S) - n) / (n * (n - 1))
    mean(w[pairs[, 1L]] + w[pairs[, 2L]] - 2 * S[pairs] / (n[pairs[, 1L]] * n[pairs[, 2L]]))
  }, numeric(1))
}

# minimum spanning tree (Prim, O(N^2)) on a distance matrix
mst_edges <- function(D) {
  N <- nrow(D)
  intree <- logical(N); intree[1L] <- TRUE
  key <- D[1L, ]; parent <- rep(1L, N); key[1L] <- Inf
  E <- matrix(0L, N - 1L, 2L)
  for (t in seq_len(N - 1L)) {
    key[intree] <- Inf
    v <- which.min(key)
    E[t, ] <- c(parent[v], v)
    intree[v] <- TRUE
    dv <- D[v, ]
    upd <- !intree & dv < key
    key[upd] <- dv[upd]; parent[upd] <- v
  }
  E
}

# Schott (2007b) k-sample covariance test, standard normal limit.
#   t = sum_{i<j} [ a_i + a_j - 2 tr(S_i S_j) ],
#   a_i = n_i^2/((n_i-1)(n_i+2)) [tr S_i^2 - (tr S_i)^2 / n_i]  (unbiased tr Sigma^2),
#   Var_0(t) = 4 a^2 [ sum_{i<j} (1/n_i + 1/n_j)^2 + (g-1)(g-2) sum_i 1/n_i^2 ],
#   n_i = N_i - 1, a estimated from the pooled covariance matrix.
schott_test <- function(data) {
  g <- length(data)
  S <- lapply(data, cov)
  nn <- vapply(data, nrow, integer(1)) - 1
  ahat <- function(Si, ni) ni^2 / ((ni - 1) * (ni + 2)) * (sum(Si^2) - sum(diag(Si))^2 / ni)
  a <- mapply(ahat, S, nn)
  t <- 0
  for (i in 1:(g - 1)) for (j in (i + 1):g) t <- t + a[i] + a[j] - 2 * sum(S[[i]] * S[[j]])
  npool <- sum(nn)
  Sp <- Reduce(`+`, Map(`*`, S, nn)) / npool
  a0 <- ahat(Sp, npool)
  pr <- t(combn(g, 2L))
  v <- 4 * a0^2 * (sum((1 / nn[pr[, 1L]] + 1 / nn[pr[, 2L]])^2) + (g - 1) * (g - 2) * sum(1 / nn^2))
  z <- t / sqrt(v)
  c(stat = z, p = pnorm(z, lower.tail = FALSE))
}

competitor_tests <- function(data, B, groups) {
  k <- length(data)
  n <- vapply(data, nrow, integer(1))
  X <- do.call(rbind, data)
  L <- rep(seq_len(k), n)
  out <- list(); tm <- list()
  if (any(c("dist") %in% groups)) {
    t0 <- proc.time()[["elapsed"]]
    Lall <- perm_labels(L, B); nb <- ncol(Lall)
    G <- group_indicator(Lall, k)
    D2 <- sqdist(X); D <- sqrt(D2)
    tm$dist_prep <- proc.time()[["elapsed"]] - t0
    tm$disco <- system.time(out$disco <- pval_perm(disco_stat(D, G, n, k, nb)))[["elapsed"]]
    tm$mmd <- system.time({
      h2 <- median(D2[upper.tri(D2)])
      K <- exp(-D2 / (2 * h2))
      out$mmd <- pval_perm(mmd_stat(K, G, n, k, nb))
    })[["elapsed"]]
    tm$mst <- system.time({
      E <- mst_edges(D)
      out$mst <- pval_perm(colSums(Lall[E[, 1L], , drop = FALSE] == Lall[E[, 2L], , drop = FALSE]))
    })[["elapsed"]]
  }
  if ("schott" %in% groups) {
    tm$schott <- system.time(out$schott <- schott_test(data)[["p"]])[["elapsed"]]
  }
  # further competitors of R/competitors2.R (G1, 1 Nov 2026), evaluated after
  # the groups above so that these keep their random numbers
  if (any(c("comp2", "comp2ovr", "songchen") %in% groups)) {
    r2 <- comp2_tests(data, B, groups)
    out <- c(out, r2$p); tm <- c(tm, r2$tm)
  }
  list(p = out, tm = tm)
}
