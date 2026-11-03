# =============================================================================
# Targeted simulation (2026-11-02): dependence-screening component and a
# coordinatewise k-sample Anderson--Darling competitor. New file; nothing in
# R/ is modified. Loaded after R/learners_lib2.R by common.R.
#
# xscrd1, xscrd3: cross-fitted screening exactly as xscr1/xscr3 (task 2.6),
#   with the dictionary extended to all pairwise mixed moments of bidegree at
#   most (2, 2):
#       z_j, z_j^2                                   (j = 1..d)
#       z_j z_l, z_j^2 z_l, z_j z_l^2, z_j^2 z_l^2   (all j < l, no cap on d)
#   Direction 1 selects on the fit half with the true labels and evaluates the
#   AUC on the test half; direction 2 selects on the test half separately for
#   EVERY labelling b and evaluates on the fit half with the fit labels.
#   T = (T1 + T2)/sqrt(2). Given D_fit, T is a function of the test points and
#   the test labels only, so the permutation p-value is exact (Theorem 1
#   argument). No random numbers are used. Pairs are processed in chunks so
#   that d = 200 (about 8e4 features) fits in memory.
#
# adk_max: for each coordinate the k-sample Anderson--Darling statistic of
#   Scholz and Stephens (1987, continuous version) on the FULL sample; the test
#   statistic is the maximum over coordinates, calibrated by B permutations of
#   all N labels (own permutations, drawn after every other test).
# =============================================================================

XSCRD_CHUNK <- 2000L          # pairs per chunk (4 features each)

# columns of the extended dictionary, identified by (type, j, l):
#   type 1: z_j; 2: z_j^2; 3: z_j z_l; 4: z_j^2 z_l; 5: z_j z_l^2; 6: z_j^2 z_l^2
xscrd_cols <- function(Z, type, j, l) {
  out <- matrix(0, nrow(Z), length(type))
  for (t in unique(type)) {
    w <- which(type == t); a <- Z[, j[w], drop = FALSE]
    b <- if (t >= 3L) Z[, l[w], drop = FALSE] else NULL
    out[, w] <- switch(t, a, a^2, a * b, a^2 * b, a * b^2, a^2 * b^2)
  }
  out
}

# dictionary index -> (type, j, l); singles first, then chunks of pairs with
# the four pair types of one chunk in order 3, 4, 5, 6
xscrd_index <- function(d, types = 1:6) {
  idx <- list()
  if (1L %in% types) idx[[length(idx) + 1L]] <- cbind(1L, seq_len(d), 0L)
  if (2L %in% types) idx[[length(idx) + 1L]] <- cbind(2L, seq_len(d), 0L)
  if (d >= 2L) {
    pr <- which(upper.tri(diag(d)), arr.ind = TRUE)      # same order as xscr
    pt <- intersect(3:6, types)
    starts <- seq(1L, nrow(pr), by = XSCRD_CHUNK)
    for (s in starts) {
      rr <- s:min(nrow(pr), s + XSCRD_CHUNK - 1L)
      for (t in pt) idx[[length(idx) + 1L]] <- cbind(t, pr[rr, 1L], pr[rr, 2L])
    }
  }
  out <- do.call(rbind, idx); colnames(out) <- c("type", "j", "l")
  storage.mode(out) <- "integer"
  out
}

# blocks of the index processed together (singles, then one block per chunk)
xscrd_blocks <- function(index, size = 4L * XSCRD_CHUNK) {
  split(seq_len(nrow(index)), ceiling(seq_len(nrow(index)) / size))
}

zstd <- function(F, mu, s) sweep(sweep(F, 2L, mu), 2L, s, "/")
safe_sd <- function(F) { s <- colsd(F); s[!is.finite(s) | s < 1e-12] <- 1; s }

# k x (B + 1) matrix; types = 1:6 is xscrd, types = 1:3 reproduces xscr for d <= 100
xscrd_T <- function(Zfit, gfit, Zte, Lall, k, m, types = 1:6, keep_diag = FALSE) {
  d <- ncol(Zfit); nf <- nrow(Zfit); M <- nrow(Zte); nb <- ncol(Lall)
  index <- xscrd_index(d, types)
  m <- min(m, nrow(index))
  # running top-m: direction 1 per population (k x m), direction 2 per (i, b)
  v1 <- matrix(-Inf, k, m); id1 <- matrix(0L, k, m); t1v <- matrix(0, k, m)
  # rows of direction 2 follow the columns of I: r = (i - 1) nb + b
  v2 <- array(-Inf, c(nb, k, m)); id2 <- array(0L, c(nb, k, m)); t2v <- array(0, c(nb, k, m))
  I <- do.call(cbind, lapply(seq_len(k), function(i) (Lall == i) * 1))   # M x (k nb)
  ni <- colSums(I)
  merge_top <- function(val, ids, tv, cand_abs, cand_id, cand_t) {
    # val, ids, tv: R x m current; cand_*: R x c candidates (c >= 1); keeps first on ties
    allv <- cbind(val, cand_abs); alli <- cbind(ids, cand_id); allt <- cbind(tv, cand_t)
    for (q in seq_len(m)) {
      j <- max.col(allv, ties.method = "first")
      rr <- cbind(seq_len(nrow(allv)), j)
      val[, q] <- allv[rr]; ids[, q] <- alli[rr]; tv[, q] <- allt[rr]
      allv[rr] <- -Inf
    }
    list(val = val, ids = ids, tv = tv)
  }
  for (blk in xscrd_blocks(index)) {
    ix <- index[blk, , drop = FALSE]
    Ff <- xscrd_cols(Zfit, ix[, 1], ix[, 2], ix[, 3]); Ft <- xscrd_cols(Zte, ix[, 1], ix[, 2], ix[, 3])
    Ffz <- zstd(Ff, colMeans(Ff), safe_sd(Ff)); Ftz <- zstd(Ft, colMeans(Ft), safe_sd(Ft))
    rm(Ff, Ft)
    p <- length(blk); ids <- matrix(blk, 1L, p)
    tot_f <- colSums(Ffz); tot_t <- colSums(Ftz)
    # direction 1
    for (i in seq_len(k)) {
      ii <- gfit == i; n1 <- sum(ii); n0 <- nf - n1
      s1 <- colSums(Ffz[ii, , drop = FALSE])
      t1 <- (s1 / n1 - (tot_f - s1) / n0) / sqrt(1 / n1 + 1 / n0)
      cand <- order(-abs(t1))[seq_len(min(m, p))]          # order() is stable: first on ties
      r <- merge_top(v1[i, , drop = FALSE], id1[i, , drop = FALSE], t1v[i, , drop = FALSE],
                     matrix(abs(t1[cand]), 1L), matrix(blk[cand], 1L), matrix(t1[cand], 1L))
      v1[i, ] <- r$val; id1[i, ] <- r$ids; t1v[i, ] <- r$tv
    }
    # direction 2, every labelling
    si <- crossprod(I, Ftz)                                 # (k nb) x p
    tt <- (si / ni - (rep(tot_t, each = k * nb) - si) / (M - ni)) / sqrt(1 / ni + 1 / (M - ni))
    rm(si)
    at <- abs(tt); cid <- matrix(0L, k * nb, m); cv <- matrix(-Inf, k * nb, m); ct <- matrix(0, k * nb, m)
    for (q in seq_len(min(m, p))) {
      j <- max.col(at, ties.method = "first"); rr <- cbind(seq_len(k * nb), j)
      cv[, q] <- at[rr]; cid[, q] <- blk[j]; ct[, q] <- tt[rr]; at[rr] <- -Inf
    }
    rm(at, tt)
    cur_v <- matrix(v2, k * nb, m); cur_i <- matrix(id2, k * nb, m); cur_t <- matrix(t2v, k * nb, m)
    r <- merge_top(cur_v, cur_i, cur_t, cv, cid, ct)
    v2[] <- r$val; id2[] <- r$ids; t2v[] <- r$tv
    rm(Ffz, Ftz)
  }
  # standardization statistics of the selected columns (label-free, per half)
  used <- sort(unique(c(id1, id2)))
  ix <- index[used, , drop = FALSE]
  Ff <- xscrd_cols(Zfit, ix[, 1], ix[, 2], ix[, 3]); Ft <- xscrd_cols(Zte, ix[, 1], ix[, 2], ix[, 3])
  mf <- colMeans(Ff); sf <- safe_sd(Ff); mt <- colMeans(Ft); st <- safe_sd(Ft)
  Ftf <- zstd(Ft, mf, sf)            # test half, fit statistics (direction 1 scores)
  Fft <- zstd(Ff, mt, st)            # fit half, test statistics (direction 2 scores)
  pos <- function(id) match(id, used)
  S1 <- matrix(0, M, k); T2 <- matrix(0, k, nb)
  for (i in seq_len(k)) {
    S1[, i] <- as.numeric(Ftf[, pos(id1[i, ]), drop = FALSE] %*% t1v[i, ])
    S2 <- matrix(0, nf, nb)
    for (q in seq_len(m)) S2 <- S2 + sweep(Fft[, pos(id2[, i, q]), drop = FALSE], 2L, t2v[, i, q], "*")
    T2[i, ] <- std_auc_cols(S2, gfit == i)
  }
  T1 <- ovr_T(S1, seq_len(k), Lall)
  ans <- (T1 + T2) / sqrt(2)
  if (keep_diag) {
    lab <- function(id) { r <- index[id, ]
      c("z%d", "z%d^2", "z%d*z%d", "z%d^2*z%d", "z%d*z%d^2", "z%d^2*z%d^2")[r[1]] |>
        (\(f) if (r[1] <= 2L) sprintf(f, r[2]) else sprintf(f, r[2], r[3]))() }
    attr(ans, "selected") <- data.frame(
      population = rep(seq_len(k), 2L), direction = rep(1:2, each = k),
      features = c(vapply(seq_len(k), function(i) paste(vapply(id1[i, ], lab, ""), collapse = ";"), ""),
                   vapply(seq_len(k), function(i) paste(vapply(id2[1L, i, ], lab, ""), collapse = ";"), "")),
      T1 = c(T1[, 1L], rep(NA, k)), T2 = c(rep(NA, k), T2[, 1L]))
  }
  ans
}

# register the new components without editing R/learners_lib2.R
if (!exists(".component_T_frozen")) .component_T_frozen <- component_T
COMPONENTS <- union(COMPONENTS, c("xscrd1", "xscrd3"))
component_T <- function(name, sc, Lall) {
  switch(name,
    xscrd1 = xscrd_T(sc$Zfit, sc$gfit, sc$Zte, Lall, sc$k, 1L),
    xscrd3 = xscrd_T(sc$Zfit, sc$gfit, sc$Zte, Lall, sc$k, 3L),
    .component_T_frozen(name, sc, Lall))
}

# ---------------------------------------------------------------- adk_max
# A^2_akN of Scholz & Stephens (1987), continuous data, for every column of
# the label matrix; x is one coordinate of the pooled sample
adk_cols <- function(x, Lmat, k) {
  N <- length(x); o <- order(x); Lo <- Lmat[o, , drop = FALSE]
  j <- seq_len(N - 1L); w <- 1 / (j * (N - j))
  tot <- 0
  for (i in seq_len(k)) {
    Mi <- apply((Lo == i) * 1, 2L, cumsum)[j, , drop = FALSE]
    n_i <- sum(Lmat[, 1L] == i)
    D <- N * Mi - j * n_i
    tot <- tot + colSums(D^2 * w) / n_i
  }
  tot / N
}

adk_max_test <- function(data, B) {
  k <- length(data); X <- do.call(rbind, data)
  L <- rep(seq_len(k), vapply(data, nrow, integer(1)))
  Lmat <- cbind(L, vapply(seq_len(B), function(b) sample(L), integer(length(L))))
  stat <- rep(-Inf, B + 1L)
  for (c in seq_len(ncol(X))) stat <- pmax(stat, adk_cols(X[, c], Lmat, k))
  list(p = (1 + sum(stat[-1L] >= stat[1L])) / (B + 1), stat = stat[1L])
}
