# =============================================================================
# Ideal AUC ("oracle AUC") of the scenarios (task section 3, written 1 Nov 2026)
#
# The ideal AUC of the deviant population is the AUC of its true log density ratio
#   Lambda(x) = log f_dev(x) - log f_ref(x),
# f_ref = the base population N(0, I) (S1, S3) or, in S2, each of the three
# departures against N(0, I) (the other populations are not used), computed by Monte
# Carlo (2e5 points per sample, fixed seed, common random numbers so that the
# AUC is a smooth monotone function of the parameter). It does not involve any
# test. oracle_solve() inverts it by bisection; tests/make_calib_s.R writes
# the result to R/calib_params_s.R, which is then frozen.
# =============================================================================

ORACLE_NMC <- 2e5L
ORACLE_SEED <- 20261101L

mc_auc <- function(sd, sr) {
  nd <- as.numeric(length(sd)); nr <- as.numeric(length(sr))
  r <- rank(c(sd, sr))
  (sum(r[seq_len(nd)]) - nd * (nd + 1) / 2) / (nd * nr)
}
logsumexp2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }

# ---------------------------------------------------------------- common random numbers
oracle_noise <- function(ncol, nmc = ORACLE_NMC, seed = ORACLE_SEED) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv)
          else rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  list(Zd = matrix(rnorm(nmc * ncol), nmc, ncol), Zr = matrix(rnorm(nmc * ncol), nmc, ncol),
       Sd = matrix(sample(c(-1, 1), nmc * ncol, replace = TRUE), nmc, ncol),
       Ud = matrix(runif(nmc * ncol), nmc, ncol), nmc = nmc)
}

# ---------------------------------------------------------------- S1a, S1b (k populations, rest = base)
oracle_S1a <- function(a, nz) {
  s <- sqrt(1 - a^2)
  x <- a * nz$Sd[, 1:S1_NDEV] + s * nz$Zd[, 1:S1_NDEV]
  l <- function(x) log(0.5) + logsumexp2(dnorm(x, a, s, log = TRUE), dnorm(x, -a, s, log = TRUE)) -
                   dnorm(x, log = TRUE)
  mc_auc(rowSums(l(x)), rowSums(l(nz$Zr[, 1:S1_NDEV])))
}
oracle_S1b <- function(u, nz) {                      # u = 1 / nu
  nu <- 1 / u; cc <- sqrt(nu / (nu - 2))
  x <- qt(nz$Ud[, 1:S1_NDEV], df = nu) / cc
  l <- function(x) log(cc) + dt(cc * x, df = nu, log = TRUE) - dnorm(x, log = TRUE)
  mc_auc(rowSums(l(x)), rowSums(l(nz$Zr[, 1:S1_NDEV])))
}

# ---------------------------------------------------------------- S2 (k = 4, equal n)
# log f_p / f_0 of populations 1..4 on the first 20 coordinates; X is n x 20
s2_logr <- function(X, th) {
  delta <- th[1L]; gamma <- th[2L]; rho <- th[3L]
  l2 <- delta * rowSums(X[, S2_BLOCKS$mean, drop = FALSE]) - length(S2_BLOCKS$mean) * delta^2 / 2
  Xv <- X[, S2_BLOCKS$var, drop = FALSE]
  l3 <- rowSums(-0.5 * log1p(gamma) + 0.5 * Xv^2 * (gamma / (1 + gamma)))
  Xp <- X[, S2_BLOCKS$pairs, drop = FALSE]
  l4 <- 0
  for (p in seq_len(ncol(Xp) %/% 2L)) {
    x1 <- Xp[, 2L * p - 1L]; x2 <- Xp[, 2L * p]
    q <- rho^2 * (x1^2 + x2^2) / (2 * (1 - rho^2))
    c <- rho * x1 * x2 / (1 - rho^2)
    l4 <- l4 - 0.5 * log1p(-rho^2) - q + log(0.5) + logsumexp2(c, -c)
  }
  cbind(0, l2, l3, l4)
}
# samples of population p (2..4) from N(0, I_20) noise Z (and signs S)
s2_sample <- function(p, Z, S, th) {
  X <- Z
  if (p == 2L) X[, S2_BLOCKS$mean] <- X[, S2_BLOCKS$mean] + th[1L]
  if (p == 3L) X[, S2_BLOCKS$var] <- X[, S2_BLOCKS$var] * sqrt(1 + th[2L])
  if (p == 4L) for (q in seq_len(length(S2_BLOCKS$pairs) %/% 2L)) {
    a <- S2_BLOCKS$pairs[2L * q - 1L]; b <- S2_BLOCKS$pairs[2L * q]
    X[, b] <- S[, q] * th[3L] * Z[, a] + sqrt(1 - th[3L]^2) * Z[, b]
  }
  X
}
# ideal AUCs of populations 2, 3, 4, EACH AGAINST N(0, I) (decision of 1 Nov 2026,
# after the first version had used the mixture of the other three populations,
# which made the three parameters depend on each other): the score is
# log f_p / f_0 and the reference sample is the base population. The strength of
# the three kinds of departure is then directly comparable and the parameters
# are solved one at a time.
oracle_S2 <- function(th, nz, which = 2:4) {
  out <- rep(NA_real_, 3L)
  for (p in which) {
    Xd <- s2_sample(p, nz$Zd, nz$Sd[, 1:5], th)
    out[p - 1L] <- mc_auc(s2_logr(Xd, th)[, p], s2_logr(nz$Zr, th)[, p])
  }
  out
}

# ---------------------------------------------------------------- S3b (pairs, rest = base)
oracle_S3b <- function(rho, nz, npair = 10L) {
  x <- nz$Zd[, 1:(2 * npair)]
  for (p in seq_len(npair)) {
    a <- 2L * p - 1L; b <- 2L * p
    x[, b] <- nz$Sd[, p] * rho * x[, a] + sqrt(1 - rho^2) * x[, b]
  }
  lr <- function(X) {
    s <- 0
    for (p in seq_len(npair)) {
      x1 <- X[, 2L * p - 1L]; x2 <- X[, 2L * p]
      q <- rho^2 * (x1^2 + x2^2) / (2 * (1 - rho^2)); c <- rho * x1 * x2 / (1 - rho^2)
      s <- s - 0.5 * log1p(-rho^2) - q + log(0.5) + logsumexp2(c, -c)
    }
    s
  }
  mc_auc(lr(x), lr(nz$Zr[, 1:(2 * npair)]))
}

# ---------------------------------------------------------------- models A, B, D, E (report only)
ar1_cov <- function(d, rho = 0.5) rho^abs(outer(seq_len(d), seq_len(d), "-"))

# par: list with cA, cB, cE. Returns the ideal AUC (mean over populations for E).
oracle_auc_old <- function(model, beta, par, d = 100L, k = 3L, nmc = ORACLE_NMC, seed = ORACLE_SEED) {
  if (beta == 0) return(0.5)
  nz <- oracle_noise(d, nmc, seed)
  j <- seq_len(d); S0 <- ar1_cov(d); Ch <- chol(S0)            # t(Ch) %*% Ch = S0
  switch(model,
    A = {
      mu <- par$cA * beta * j^(-3)
      Si <- solve(S0); w <- as.numeric(Si %*% mu)
      mc_auc(as.numeric((nz$Zd %*% Ch + rep(mu, each = nmc)) %*% w), as.numeric((nz$Zr %*% Ch) %*% w))
    },
    B = {
      S1 <- S0 + diag(par$cB * beta * j^(-3), nrow = d); Q <- solve(S0) - solve(S1)
      sc <- function(X) 0.5 * rowSums((X %*% Q) * X)
      Xd <- nz$Zd %*% chol(S1)
      mc_auc(sc(Xd), sc(nz$Zr %*% Ch))
    },
    D = {
      s <- min(d, ceiling(100 * beta))
      SD <- diag(d); SD[seq_len(s), seq_len(s)] <- ar1_cov(s, 0.9)
      SDi <- solve(SD)
      lg <- function(X, Si) {                                  # log of the +-1 mixture, up to constants
        o <- rep(1, d)
        a <- -0.5 * rowSums(((X - rep(o, each = nrow(X))) %*% Si) * (X - rep(o, each = nrow(X))))
        b <- -0.5 * rowSums(((X + rep(o, each = nrow(X))) %*% Si) * (X + rep(o, each = nrow(X))))
        logsumexp2(a, b)
      }
      sgn <- nz$Sd[, 1]
      Xd <- nz$Zd; Xd[, seq_len(s)] <- Xd[, seq_len(s), drop = FALSE] %*% chol(ar1_cov(s, 0.9))
      Xd <- Xd + sgn
      Xr <- nz$Zr + sgn
      ld <- -0.5 * as.numeric(determinant(SD)$modulus)
      sc <- function(X) lg(X, SDi) + ld - lg(X, diag(d))
      mc_auc(sc(Xd), sc(Xr))
    },
    E = {
      Si <- solve(S0)
      mus <- lapply(seq_len(k), function(g) { m <- numeric(d); blk <- (5L * (g - 1L) + 1L):(5L * g)
                                              m[blk[blk <= d]] <- par$cE * beta; m })
      lN <- function(X, m) { Z <- X - rep(m, each = nrow(X)); -0.5 * rowSums((Z %*% Si) * Z) }
      mean(vapply(seq_len(k), function(g) {
        others <- setdiff(seq_len(k), g)
        lam <- function(X) lN(X, mus[[g]]) - (log(rowSums(exp(sapply(others, function(h) lN(X, mus[[h]]))))) - log(length(others)))
        Xd <- nz$Zd %*% Ch + rep(mus[[g]], each = nmc)
        Xr <- nz$Zr %*% Ch + rep(mus[[others[1L]]], each = nmc)
        # reference: equal parts of the other populations
        part <- nmc %/% length(others)
        for (q in seq_along(others)) { ix <- ((q - 1L) * part + 1L):(q * part)
          Xr[ix, ] <- nz$Zr[ix, ] %*% Ch + rep(mus[[others[q]]], each = length(ix)) }
        mc_auc(lam(Xd), lam(Xr))
      }, numeric(1)))
    },
    stop("no oracle for model ", model))
}

# ---------------------------------------------------------------- inversion
# increasing function f on [lo, hi]; returns the root of f(x) = target, or NA
oracle_solve <- function(f, target, lo, hi, tol = 1e-5, maxit = 40L) {
  started <- proc.time()[["elapsed"]]; neval <- 0L
  evaluate <- function(x) {
    neval <<- neval + 1L
    verbose <- isTRUE(getOption("power3.progress", FALSE))
    if (verbose) { cat(sprintf("  oracle target %.2f | eval %d | parameter %.7g ...\n", target, neval, x)); flush.console() }
    value <- f(x)
    if (verbose) { cat(sprintf("    AUC %.6f | elapsed %.1fs\n", value, proc.time()[["elapsed"]] - started)); flush.console() }
    value
  }
  flo <- evaluate(lo); fhi <- evaluate(hi)
  if (target <= flo) return(lo)
  if (target > fhi) { warning(sprintf("target %.3f above the maximum %.3f in [%g, %g]", target, fhi, lo, hi)); return(NA_real_) }
  for (it in seq_len(maxit)) {
    mid <- (lo + hi) / 2
    if (evaluate(mid) < target) lo <- mid else hi <- mid
    if (hi - lo < tol) break
  }
  (lo + hi) / 2
}

# natural parameters for the ideal AUC `target`; S2 returns c(delta, gamma, rho)
oracle_param <- function(model, target, par = list(cB = 3)) {
  # Generate common random numbers once per root, not once per bisection step.
  nc <- if (model %in% c("S1a", "S1b")) S1_NDEV else 20L
  if (model != "S3a") nz <- oracle_noise(nc)
  switch(model,
    S1a = oracle_solve(function(a) oracle_S1a(a, nz), target, 1e-6, 0.999),
    S1b = oracle_solve(function(u) oracle_S1b(u, nz), target, 1e-6, 0.49),
    S3b = oracle_solve(function(r) oracle_S3b(r, nz), target, 1e-6, 0.999),
    S3a = {
      # In whitened eigen-coordinates the reference variances are 1 and the
      # alternative variances are 1 + b*ev. These are exactly the model-B
      # density-ratio distributions, without repeated 200000 x 100 products.
      d <- 100L
      Cinv <- solve(chol(ar1_cov(d)))
      ev <- eigen(t(Cinv) %*% diag(par$cB * seq_len(d)^(-3)) %*% Cinv,
                  symmetric = TRUE, only.values = TRUE)$values
      noise <- oracle_noise(d)
      Zd2 <- noise$Zd^2; Zr2 <- noise$Zr^2
      rm(noise)
      f <- function(b) mc_auc(as.numeric(Zd2 %*% (b * ev)),
                              as.numeric(Zr2 %*% (b * ev / (1 + b * ev))))
      oracle_solve(f, target, 1e-6, 1)
    },
    S2 = {
      nz <- oracle_noise(20L)
      th <- c(0.1, 0.1, 0.1)                       # the other entries do not matter
      ub <- c(5, 20, 0.999)
      for (q in 1:3) {
        f <- function(v) { t2 <- th; t2[q] <- v; oracle_S2(t2, nz, which = q + 1L)[q] }
        th[q] <- oracle_solve(f, target, 1e-6, ub[q])
      }
      th
    },
    stop("unknown model: ", model))
}
