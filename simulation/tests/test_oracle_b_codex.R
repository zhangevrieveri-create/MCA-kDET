# Focused deterministic checks for the B-oracle sign and spectral S3a solver.
# Run from sim/: Rscript --vanilla tests/test_oracle_b_codex.R
source("R/models_s.R")
source("R/oracle_s.R")
set.seed(90210)
d <- 10L; beta <- 0.4; cB <- 3
S0 <- ar1_cov(d); Cinv <- solve(chol(S0))
ee <- eigen(t(Cinv) %*% diag(cB * seq_len(d)^(-3)) %*% Cinv, symmetric = TRUE)
Z <- matrix(rnorm(2000), 200, d)
S1 <- S0 + diag(cB * beta * seq_len(d)^(-3))
Q <- solve(S0) - solve(S1)
# X = Z V' C has covariance S0. Under the alternative replace Z by
# Z diag(sqrt(1+beta*ev)); both scores omit the same log-determinant constant.
Xr <- Z %*% t(ee$vectors) %*% chol(S0)
Xd <- sweep(Z, 2, sqrt(1 + beta * ee$values), "*") %*% t(ee$vectors) %*% chol(S0)
direct <- function(X) rowSums((X %*% Q) * X)
stopifnot(max(abs(direct(Xr) - as.numeric(Z^2 %*% (beta * ee$values / (1 + beta * ee$values))))) < 1e-10,
          max(abs(direct(Xd) - as.numeric(Z^2 %*% (beta * ee$values)))) < 1e-10)
cat("OK spectral B scores equal direct Gaussian density-ratio scores\n")
auc <- oracle_auc_old("B", beta, list(cB = cB), d = 1L)
exact <- 2 / pi * atan(sqrt(1 + cB * beta))
stopifnot(abs(auc - exact) < 0.002)
cat(sprintf("OK B oracle AUC %.6f; one-dimensional closed form %.6f\n", auc, exact))
cat("ORACLE B CHECKS OK\n")
