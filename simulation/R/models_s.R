# =============================================================================
# New scenarios S1-S3 (third review round, TASK_第三轮_功效与场景_20261101.md,
# section 3). Written 1 Nov 2026. Models A-E in R/models.R are untouched.
#
# For every S model the cell column `beta` is the IDEAL AUC minus 0.5: the AUC
# of the true log density ratio of the deviant population against the
# other populations (S1/S3); S2 calibrates each departure against N(0,I),
# independently of the other departures. beta = 0 is the null. The natural
# parameters of the models (a, 1/nu, delta, gamma, rho, ...) belong to the
# table S_CALIB in R/calib_params_s.R, produced by tests/make_calib_s.R from
# R/oracle_s.R and frozen; it does not depend on any test.
#
#   S1a  k populations, deviant = population k: the first S1_NDEV = 3
#        coordinates are bimodal, X_j = a S + sqrt(1 - a^2) Z, S = +-1, Z ~ N(0,1);
#        mean 0, variance 1, independent coordinates. Base: N(0, I_d).
#   S1b  same coordinates standardized t(nu), signal parameter 1/nu.
#   S2   k = 4, base N(0, I_d). Pop 1 base; pop 2 mean shift delta on
#        coordinates 1-5; pop 3 variance 1 + gamma on coordinates 6-10; pop 4
#        sign-mixed correlation (X1, X2) ~ 1/2 N(0, R_rho) + 1/2 N(0, R_-rho) on
#        the 5 pairs (11,12), ..., (19,20).
#   S3a  k in {5, 10, 20}: deviant = (M2)-type sparse variance change (model B).
#   S3b  k in {5, 10, 20}, d = 20: deviant = sign-mixed correlation on all 10 pairs.
# =============================================================================

S_MODELS <- c("S1a", "S1b", "S2", "S3a", "S3b")
S1_NDEV <- 3L
S2_BLOCKS <- list(mean = 1:5, var = 6:10, pairs = 11:20)
S_LEVELS <- c(0.04, 0.07, 0.10, 0.14, 0.18)      # ideal AUC 0.54 ... 0.68

if (exists("SRC_DIR") && file.exists(file.path(SRC_DIR, "calib_params_s.R")))
  source(file.path(SRC_DIR, "calib_params_s.R"))
if (!exists("S_CALIB")) S_CALIB <- NULL

# natural parameters of model at level beta; S2 returns c(delta, gamma, rho)
s_param <- function(model, beta) {
  key <- sprintf("%.2f", beta)
  v <- if (is.null(S_CALIB)) NULL else S_CALIB[[model]][[key]]
  if (is.null(v) || anyNA(v))
    stop("scenario ", model, " is not calibrated at beta = ", key,
         ": run tests/make_calib_s.R and fill in R/calib_params_s.R")
  v
}

# sign-mixed correlation on the pairs (cols[1], cols[2]), (cols[3], cols[4]), ...
# X is n x d with N(0, I) entries; the second coordinate of every pair is
# replaced by S rho X1 + sqrt(1 - rho^2) X2 with S = +-1 (independent per row and pair)
sign_mix_pairs <- function(X, cols, rho) {
  n <- nrow(X)
  for (p in seq_len(length(cols) %/% 2L)) {
    a <- cols[2L * p - 1L]; b <- cols[2L * p]
    S <- sample(c(-1, 1), n, replace = TRUE)
    X[, b] <- S * rho * X[, a] + sqrt(1 - rho^2) * X[, b]
  }
  X
}

gen_population_s <- function(model, beta, n, d, g, k, par) {
  dev <- (g == k)
  switch(model,
    S1a = {
      X <- matrix(rnorm(n * d), n, d)
      if (dev && beta > 0) {
        a <- s_param(model, beta)[1L]
        s <- min(S1_NDEV, d)
        S <- matrix(sample(c(-1, 1), n * s, replace = TRUE), n, s)
        X[, seq_len(s)] <- a * S + sqrt(1 - a^2) * X[, seq_len(s), drop = FALSE]
      }
      X
    },
    S1b = {
      X <- matrix(rnorm(n * d), n, d)
      if (dev && beta > 0) {
        nu <- 1 / s_param(model, beta)[1L]
        s <- min(S1_NDEV, d)
        X[, seq_len(s)] <- matrix(rt(n * s, df = nu), n, s) / sqrt(nu / (nu - 2))
      }
      X
    },
    S2 = {
      if (k != 4L) stop("S2 has k = 4 populations")
      X <- matrix(rnorm(n * d), n, d)
      if (beta > 0 && g > 1L) {
        th <- s_param(model, beta)                    # delta, gamma, rho
        if (g == 2L) X[, S2_BLOCKS$mean] <- X[, S2_BLOCKS$mean] + th[1L]
        if (g == 3L) X[, S2_BLOCKS$var] <- X[, S2_BLOCKS$var] * sqrt(1 + th[2L])
        if (g == 4L) X <- sign_mix_pairs(X, S2_BLOCKS$pairs, th[3L])
      }
      X
    },
    S3a = gen_population("B", if (dev && beta > 0) s_param(model, beta)[1L] else 0,
                         n, d, g, k, par),
    S3b = {
      X <- matrix(rnorm(n * d), n, d)
      if (dev && beta > 0) X <- sign_mix_pairs(X, seq_len(2L * (d %/% 2L)), s_param(model, beta)[1L])
      X
    },
    stop("unknown scenario: ", model))
}

gen_data_s <- function(model, beta, n, d, par) {
  k <- length(n)
  lapply(seq_len(k), function(g) gen_population_s(model, beta, n[g], d, g, k, par))
}

# one entry point for one_rep()
gen_any <- function(model, beta, n, d, par)
  if (model %in% S_MODELS) gen_data_s(model, beta, n, d, par) else gen_data(model, beta, n, d, par)
