# Ideal-AUC calibration of the T1/T2 sign-mixed-correlation model (3 pairs).
# Uses oracle_S3b() of R/oracle_s.R (frozen) with npair = 3 and a NEW seed,
# not used by any test. Run from sim/:  Rscript 靶向模拟_20261102/calibrate.R
SRC_DIR <- file.path(getwd(), "R")
source(file.path(SRC_DIR, "oracle_s.R"))
OUT <- "靶向模拟_20261102"
NPAIR <- 3L
SEED_CAL <- 20261103L; SEED_CHECK <- 20261104L
targets <- c(0.65, 0.70, 0.75)

nz <- oracle_noise(2L * NPAIR, seed = SEED_CAL)
rho <- vapply(targets, function(t)
  oracle_solve(function(r) oracle_S3b(r, nz, npair = NPAIR), t, 1e-6, 0.999), numeric(1))
nz2 <- oracle_noise(2L * NPAIR, seed = SEED_CHECK)
check <- vapply(rho, function(r) oracle_S3b(r, nz2, npair = NPAIR), numeric(1))

# moment check of the generator: marginals N(0,1), covariance I, |X1|,|X2| correlated
set.seed(SEED_CHECK + 1L)
mom <- t(vapply(rho, function(r) {
  X <- sign_mix_check <- matrix(rnorm(4e5), 2e5, 2)
  S <- sample(c(-1, 1), 2e5, replace = TRUE)
  X[, 2] <- S * r * X[, 1] + sqrt(1 - r^2) * X[, 2]
  c(mean2 = mean(X[, 2]), var2 = var(X[, 2]), cor12 = cor(X[, 1], X[, 2]),
    ks2 = suppressWarnings(ks.test(X[, 2], "pnorm")$p.value),
    cor_abs = cor(abs(X[, 1]), abs(X[, 2])))
}, numeric(5)))

tab <- data.frame(target = targets, rho = rho, auc_check_new_seed = check, mom)
print(tab, digits = 6)
write.csv(tab, file.path(OUT, "calibration_check.csv"), row.names = FALSE)
saveRDS(setNames(as.list(rho), sprintf("%.2f", targets)), file.path(OUT, "calibration.rds"))
