# =============================================================================
# MCA-kDET simulation study -- data-generating models (Section 4 design)
#
# Populations 1, ..., k-1 follow P_base, population k follows P_dev
# (models A, B, C', D); beta = 0 is the null hypothesis.
# Model E is the "all populations differ slightly" model used in E4.
# =============================================================================

# n x d draws of N(0, R), R_ij = rho^|i-j|, by the AR(1) recursion
rar1 <- function(n, d, rho) {
  Z <- matrix(rnorm(n * d), n, d)
  if (d == 1L || rho == 0) return(Z)
  cc <- sqrt(1 - rho^2)
  for (j in 2:d) Z[, j] <- rho * Z[, j - 1L] + cc * Z[, j]
  Z
}

# One population.
#   g      : population index (1..k)
#   k      : number of populations
#   par    : list with cA, cB, sC (number of pairs in C'), cE
gen_population <- function(model, beta, n, d, g, k, par) {
  dev <- (g == k)
  j <- seq_len(d)
  switch(model,
    # A: sparse mean shift, Sigma_ij = 0.5^|i-j|, mu_j = cA * beta * j^-3
    A = {
      X <- rar1(n, d, 0.5)
      if (dev && beta > 0) X <- sweep(X, 2L, par$cA * beta * j^(-3), "+")
      X
    },
    # B: sparse variance shift, Sigma'_jj = 1 + cB * beta * j^-3 (off-diagonals unchanged)
    B = {
      X <- rar1(n, d, 0.5)
      if (dev && beta > 0) {
        extra <- sqrt(par$cB * beta * j^(-3))
        X <- X + sweep(matrix(rnorm(n * d), n, d), 2L, extra, "*")
      }
      X
    },
    # C': equal marginals N(0,1) and equal covariance I_d; the first 2*sC
    #     coordinates of the deviant population form sC pairs that, with
    #     probability beta, equal (U, S*U) with S = +-1 independent.
    C = {
      X <- matrix(rnorm(n * d), n, d)
      if (dev && beta > 0) {
        s <- min(par$sC, floor(d / 2))
        for (p in seq_len(s)) {
          a <- 2L * p - 1L; b <- 2L * p
          hit <- runif(n) < beta
          S <- sample(c(-1, 1), n, replace = TRUE)
          X[hit, b] <- S[hit] * X[hit, a]
        }
      }
      X
    },
    # C2: like C', but the second coordinate of each pair is, with probability
    #     beta, V = Phi^{-1}(2 Phi(|U|) - 1): still N(0,1) and uncorrelated
    #     with U (covariance I_d), but E(U^2 V) > 0, a third-order signal.
    C2 = {
      X <- matrix(rnorm(n * d), n, d)
      if (dev && beta > 0) {
        s <- min(par$sC, floor(d / 2))
        for (p in seq_len(s)) {
          a <- 2L * p - 1L; b <- 2L * p
          hit <- runif(n) < beta
          u <- X[hit, a]
          X[hit, b] <- qnorm(pmin(2 * pnorm(abs(u)) - 1, 1 - 1e-12))
        }
      }
      X
    },
    # D: Gaussian mixture with a shared +-1 component; the deviant population
    #    replaces the first s = ceil(100 beta) coordinates of the noise by an
    #    AR(0.9) block (identical marginals, different copula and covariance).
    D = {
      sgn <- sample(c(-1, 1), n, replace = TRUE)
      E <- matrix(rnorm(n * d), n, d)
      if (dev && beta > 0) {
        s <- min(d, ceiling(100 * beta))
        E[, seq_len(s)] <- rar1(n, s, 0.9)
      }
      E + sgn
    },
    # E: dense alternative -- every population shifts its own block of 5
    #    coordinates by cE * beta (AR(0.5) covariance); all k populations
    #    depart from their references by the same amount.
    E = {
      X <- rar1(n, d, 0.5)
      if (beta > 0) {
        blk <- (5L * (g - 1L) + 1L):(5L * g)
        blk <- blk[blk <= d]
        X[, blk] <- X[, blk] + par$cE * beta
      }
      X
    },
    stop("unknown model: ", model))
}

# n: vector of k sample sizes
gen_data <- function(model, beta, n, d, par) {
  k <- length(n)
  lapply(seq_len(k), function(g) gen_population(model, beta, n[g], d, g, k, par))
}
