# =============================================================================
# Checks of the code of TASK_第三轮_功效与场景_20261101.md (written 1 Nov 2026):
# R/learners_lib2.R (l1th, l1se, scrm, scr2, xscr1, xscr3, g5t, g5f),
# R/competitors2.R (one-vs-rest max-type tests), R/models_s.R (S1-S3),
# R/oracle_s.R, the grids libpilot3, E11-E13. Writes nothing.
#
#   "C:/Program Files/R/R-4.3.3/bin/Rscript.exe" tests/test_power3.R
#
# Every line starts with OK / FAIL / SKIP; the last line says POWER3 TESTS OK
# or POWER3 TESTS FAILED. Needs ranger, gTestsMulti is not used here.
# =============================================================================
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "competitors2.R",
            "rotation.R", "learners_lib2.R", "models_s.R", "oracle_s.R", "experiments.R"))
  source(file.path(SRC_DIR, f))
source(file.path(SRC_DIR, "summarize.R"))
suppressPackageStartupMessages(requireNamespace("glmnet"))
suppressPackageStartupMessages(requireNamespace("ranger", quietly = TRUE))
if (is.na(LIB_MAIN_V3)) LIB_MAIN_V3 <- "l1,quad,knn,rf"   # only ids and seeds are checked here
T3_PROGRESS_START <- proc.time()[["elapsed"]]

# R on this machine sometimes fails inside base functions ("locked binding", "could not
# find function"; seen with project-free base R code too, CLAUDE.md). As safe_rep() in
# sim.R, a failed call is repeated; the same seed gives the same result.
retry_f <- function(f, tries = 6L) {
  for (a in seq_len(tries)) {
    r <- tryCatch(f(), error = function(e) e)
    if (!inherits(r, "error")) return(r)
    cat("    (retry", a, ":", conditionMessage(r), ")
")
  }
  stop("failed ", tries, " times: ", conditionMessage(r))
}
# cached loop: the results of the first r replications survive a crash of R
# (delete tests/out/t3cache_*.rds after changing the code under test)
cvapply <- function(name, n, f, len) {
  cf <- file.path("tests", "out", paste0("t3cache_", name, ".rds"))
  res <- if (file.exists(cf)) tryCatch(readRDS(cf), error = function(e) list()) else list()
  length(res) <- n
  n_cached <- sum(!vapply(res, is.null, logical(1)))
  cat(sprintf("[%s] CACHE %s: %d/%d already saved\n", format(Sys.time(), "%H:%M:%S"), name, n_cached, n)); flush.console()
  for (r in seq_len(n)) {
    if (is.null(res[[r]])) {
      cat(sprintf("[%s] START %s %d/%d\n", format(Sys.time(), "%H:%M:%S"), name, r, n)); flush.console()
      res[[r]] <- f(r)
      if (r %% 10L == 0L || r == n) saveRDS(res, cf)
      pbar(r, n, name)
    }
  }
  pbar(n, n, paste(name, "including cache"))
  t(vapply(res, identity, numeric(len)))
}
# progress bar for the long loops (printed on one line per update)
pbar <- function(i, n, what = "") {
  if (TRUE) {
    nb <- round(30 * i / n)
    cat(sprintf("  [%s%s] %3.0f%% %s %d/%d
", strrep("#", nb), strrep("-", 30 - nb), 100 * i / n, what, i, n))
    cat(sprintf("    elapsed %.1fs\n", proc.time()[["elapsed"]] - T3_PROGRESS_START)); flush.console()
  }
}
all_ok <- TRUE
chk <- function(cond, msg) {
  cond <- isTRUE(cond)
  cat(if (cond) "OK   " else "FAIL ", msg, "\n", sep = "")
  all_ok <<- all_ok && cond
}
par0 <- list(cA = 0.6, cB = 4, sC = 2L, cE = 0.2)
mkcell <- function(model, beta, k, n, d, lib, groups = "ovr", B = 99L, ...) {
  g <- cell_row("t3", model, beta, k = k, n = n, d = d, groups = groups, B = B, reps = 1L,
                lib = lib, seed_scheme = "v2", ...)
  as.list(g)
}

# ---- 1 toy problems: direction of the new learners (strong signals, fresh test data)
set.seed(101)
toy <- function(n, d, kind) {
  y <- rep(0:1, each = n); X <- matrix(rnorm(2 * n * d), 2 * n, d)
  one <- y == 1
  if (kind == "mean") X[one, 1] <- X[one, 1] + 1.5
  if (kind == "var") X[one, 1] <- X[one, 1] * 2.5
  list(X = X, y = y)
}
auc_test <- function(name, kind, n = 150, d = 10) {
  tr <- toy(n, d, kind); te <- toy(500, d, kind)
  w <- rep(1, nrow(tr$X)); grp <- tr$y + 1L
  o <- get(paste0("learner_", name))(tr$X, tr$y, w, grp, te$X)
  auc_vec(o$score, te$y == 1)
}
for (nm in c("l1th", "l1se", "scrm")) {
  a <- auc_test(nm, "mean")
  chk(a > 0.75, sprintf("1 %s detects a mean shift: test AUC %.3f", nm, a))
}
for (nm in c("scr2")) {
  a <- auc_test(nm, "var")
  chk(a > 0.65, sprintf("1 %s detects a variance change: test AUC %.3f", nm, a))
}
# screening picks the right coordinate
tr <- toy(150, 10, "mean")
o <- learner_scrm(tr$X, tr$y, rep(1, nrow(tr$X)), tr$y + 1L, tr$X)
chk(o$info$m %in% SCREEN_MS, sprintf("1 scrm selects m = %d", o$info$m))

# ---- 2 xscr_T against a naive implementation (all selections redone per labelling)
naive_xscr <- function(Zfit, gfit, Zte, Lall, k, m) {
  d <- ncol(Zfit)
  pairs <- if (d <= XSCR_MAX_D_CROSS) which(upper.tri(diag(d)), arr.ind = TRUE) else matrix(0L, 0, 2)
  fe <- function(Z) cbind(Z, Z^2, Z[, pairs[, 1]] * Z[, pairs[, 2]])
  Ff <- fe(Zfit); Ft <- fe(Zte)
  tstat <- function(F, lab, i) {
    Fz <- scale(F); pos <- lab == i
    (colMeans(Fz[pos, , drop = FALSE]) - colMeans(Fz[!pos, , drop = FALSE])) /
      sqrt(1 / sum(pos) + 1 / sum(!pos))
  }
  stdz <- function(F, mu, sdv) sweep(sweep(F, 2, mu), 2, sdv, "/")
  stdauc <- function(s, pos) {
    n <- length(s); n1 <- sum(pos); n0 <- n - n1; r <- rank(s)
    ((sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0) - 0.5) / sqrt((n + 1) / (12 * n1 * n0))
  }
  mf <- colMeans(Ff); sf <- apply(Ff, 2, sd); mt <- colMeans(Ft); st <- apply(Ft, 2, sd)
  out <- matrix(NA_real_, k, ncol(Lall))
  for (b in seq_len(ncol(Lall))) for (i in seq_len(k)) {
    t1 <- tstat(Ff, gfit, i); s1 <- order(-abs(t1))[seq_len(m)]
    sc1 <- stdz(Ft, mf, sf)[, s1, drop = FALSE] %*% t1[s1]
    T1 <- stdauc(as.numeric(sc1), Lall[, b] == i)
    t2 <- tstat(Ft, Lall[, b], i); s2 <- order(-abs(t2))[seq_len(m)]
    sc2 <- stdz(Ff, mt, st)[, s2, drop = FALSE] %*% t2[s2]
    T2 <- stdauc(as.numeric(sc2), gfit == i)
    out[i, b] <- (T1 + T2) / sqrt(2)
  }
  out
}
set.seed(7)
k <- 3; d <- 4
Zfit <- matrix(rnorm(45 * d), 45, d); gfit <- rep(1:3, each = 15)
Zte <- matrix(rnorm(48 * d), 48, d); Lte <- rep(1:3, each = 16)
Zte[Lte == 3, 1] <- Zte[Lte == 3, 1] * 2
Lall <- perm_labels(Lte, 6L)
for (m in c(1L, 3L)) {
  a <- xscr_T(Zfit, gfit, Zte, Lall, k, m); b <- naive_xscr(Zfit, gfit, Zte, Lall, k, m)
  chk(max(abs(a - b)) < 1e-8, sprintf("2 xscr%d equals the naive implementation (max diff %.1e)", m, max(abs(a - b))))
  rng <- .Random.seed
  chk(identical(a, xscr_T(Zfit, gfit, Zte, Lall, k, m)) && identical(rng, .Random.seed),
      sprintf("2 xscr%d repeated statistics are bit-identical and use no RNG", m))
  cols <- rev(seq_len(ncol(Lall)))
  chk(max(abs(xscr_T(Zfit, gfit, Zte, Lall[, cols], k, m) - a[, cols])) < 1e-8,
      sprintf("2 xscr%d is equivariant to permutation-column reordering", m))
}

# ---- 3 g5t / g5f counts against a direct count
set.seed(8)
Z1 <- matrix(rnorm(30 * 3), 30, 3); Z2 <- matrix(rnorm(24 * 3), 24, 3)
g1 <- rep(1:3, each = 10); L2 <- rep(1:3, each = 8); La <- perm_labels(L2, 3L)
Et <- kmst_edges(sqrt(sqdist(Z2)), SC_KMST)
Wt <- g5_counts(Et, La, rep(1, nrow(Et)), 3L)
direct <- sapply(1:ncol(La), function(b) sapply(1:3, function(i) sum(La[Et[, 1], b] == i & La[Et[, 2], b] == i)))
chk(all(Wt == direct), "3a g5t within-population edge counts equal a direct count")
Ef <- kmst_edges(sqrt(sqdist(rbind(Z1, Z2))), SC_KMST); nf <- nrow(Z1)
labs <- function(b) c(g1, La[, b])
direct <- sapply(1:ncol(La), function(b) sapply(1:3, function(i) {
  l <- labs(b); sum(l[Ef[, 1]] == i & l[Ef[, 2]] == i & (Ef[, 1] > nf | Ef[, 2] > nf)) }))
Wf <- g5_counts(Ef, rbind(matrix(g1, nf, ncol(La)), La), as.numeric(Ef[, 1] > nf | Ef[, 2] > nf), 3L)
chk(all(Wf == direct), "3b g5f counts (at least one test end, fit labels fixed) equal a direct count")

# ---- 4 one-vs-rest max-type statistics against direct computation
set.seed(9)
X <- matrix(rnorm(60 * 5), 60, 5); L <- rep(1:3, each = 20); La <- perm_labels(L, 4L)
rm1 <- raw_moments(X[1:20, ]); cf <- cov_from_moments(rm1$n, rm1$m1, rm1$m2, rm1$m21, rm1$m22)
Y <- sweep(X[1:20, ], 2, colMeans(X[1:20, ]))
S_dir <- crossprod(Y) / 20
Th_dir <- matrix(0, 5, 5); for (l in 1:5) for (m in 1:5) Th_dir[l, m] <- mean((Y[, l] * Y[, m] - S_dir[l, m])^2)
chk(max(abs(cf$S - S_dir)) < 1e-10 && max(abs(cf$Theta - Th_dir)) < 1e-10,
    "4a covariance and theta from raw moments equal the direct computation")
ut <- upper.tri(diag(5), diag = TRUE)
dir_cov <- function(b) {
  best <- -Inf
  for (i in 1:3) {
    pos <- La[, b] == i
    f <- function(Z) { Y <- sweep(Z, 2, colMeans(Z)); n <- nrow(Z); S <- crossprod(Y) / n
      Th <- matrix(0, 5, 5); for (l in 1:5) for (m in 1:5) Th[l, m] <- mean((Y[, l] * Y[, m] - S[l, m])^2)
      list(S = S, Th = Th, n = n) }
    a <- f(X[pos, ]); r <- f(X[!pos, ])
    v <- (a$S - r$S)^2 / (a$Th / a$n + r$Th / r$n)
    best <- max(best, max(v[ut]))
  }
  best
}
fast <- clx_cov_ovr_stat(X, La, 3L)
chk(max(abs(fast - sapply(seq_len(ncol(La)), dir_cov))) < 1e-8, "4b clx_cov_ovr equals the direct computation")
G <- group_indicator(La, 3L); vv <- apply(X, 2, var)
dir_mean <- function(b) max(sapply(1:3, function(i) { pos <- La[, b] == i
  max((colMeans(X[pos, ]) - colMeans(X[!pos, ]))^2 / vv) / (1 / sum(pos) + 1 / sum(!pos)) }))
chk(max(abs(maxmean_ovr_stat(X, G, rep(20, 3), 3L, ncol(La), vv) - sapply(seq_len(ncol(La)), dir_mean))) < 1e-10,
    "4c maxmean_ovr equals the direct computation")

# ---- 5 size of every new learner and component as a library of its own (200 null reps)
lib5 <- "l1th,l1se,scrm,scr2,xscr1,xscr3,g5t,g5f"
cell5 <- mkcell("A", 0, 3L, 60L, 10L, lib5)
pp <- cvapply("null_reps", 200L, function(r) {
  one <- retry_f(function() one_rep(cell5, r))
  c(unlist(one[paste0("mca_", strsplit(lib5, ",")[[1L]])]), l1th_nnz = one$l1th_nnz)
}, 9L)
nnz5 <- pp[, 9L]; pp <- pp[, 1:8]
colnames(pp) <- strsplit(lib5, ",")[[1L]]; rr <- colMeans(pp <= 0.05)
# l1th with the theoretical penalty empties the model under H0 (nnz = 0): its scores are then
# constant, T = 0 for every labelling and p = 1, so it is conservative BY DESIGN; its lower
# bound is 0, the upper bound is the same as for the others
lo5 <- ifelse(names(rr) == "l1th", 0, 0.02)
cat(sprintf("  l1th: mean number of nonzero slopes under H0 %.2f, share of replications with an empty model %.2f
",
            mean(nnz5), mean(nnz5 == 0)))
chk(all(rr >= lo5 & rr <= 0.09),
    paste("5 size in [0.02, 0.09] over 200 null reps:", paste(sprintf("%s %.3f", names(rr), rr), collapse = ", ")))
set.seed(1); cell5b <- mkcell("A", 0, 3L, 60L, 10L, paste("l1,quad,knn,rf", lib5, sep = ","))
pb <- retry_f(function() one_rep(cell5b, 3L)); pb2 <- retry_f(function() one_rep(cell5b, 3L))
pc6 <- grep("^(mca|macc|pair|shared)", names(pb), value = TRUE)
pb <- unlist(pb[pc6]); pb2 <- unlist(pb2[pc6])
chk(identical(pb, pb2), "6 the same cell and replication gives bit-identical p-values")

# ---- 7 old library: an E2 replication of the quick run is reproduced
same15 <- function(a, b) isTRUE(all(abs(a - b) <= 1e-13 * pmax(abs(a), abs(b)) | a == b))
e2 <- make_grid("E2")
id <- e2$cell[e2$model == "D" & abs(e2$beta - 0.03) < 1e-9]
f7 <- file.path("results_v2_quick", "E2", sprintf("%s__r00001-00050.csv", id))
if (length(id) == 1L && file.exists(f7)) {
  old <- read.csv(f7, stringsAsFactors = FALSE, check.names = FALSE); o <- old[old$rep == 1L, ]
  new <- retry_f(function() one_rep(as.list(e2[e2$cell == id, ]), 1L))
  common <- intersect(names(new), names(o))
  pcols <- intersect(union(P_COLS, grep("^(mca_|pair_)", common, value = TRUE)), common)
  ok7 <- vapply(pcols, function(x) same15(as.numeric(new[[x]]), o[[x]]), logical(1))
  chk(all(ok7), sprintf("7 %s rep 1: %d old p-value columns unchanged with the old library", id, length(pcols)))
  if (!all(ok7)) print(data.frame(col = pcols[!ok7], old = unlist(o[pcols[!ok7]]), new = unlist(new[pcols[!ok7]])))
} else chk(FALSE, "7 required results_v2_quick/E2 chunk not found")

# ---- 8 sign-mixed correlation is detected by the graph components (d = 2, strong)
set.seed(88)
pg <- cvapply("sign_mix", 20L, function(r) {
  retry_f(function() {
  set.seed(88L + r)
  data <- list(matrix(rnorm(200), 100, 2), matrix(rnorm(200), 100, 2), matrix(rnorm(200), 100, 2))
  data[[3]] <- sign_mix_pairs(data[[3]], 1:2, 0.95)
  sp <- split_data(data, 0.5)
  sc <- fit_scores(sp, "ovr", lib = c("g5t", "g5f", "knn"))
  unlist(classifier_tests(sc, 99L)$p[c("mca_g5t", "mca_g5f")])
  })
}, 2L)
rg <- colMeans(pg <= 0.05)
chk(all(rg >= 0.5), sprintf("8 g5t / g5f detect a sign-mixed correlation: %.2f %.2f", rg[1], rg[2]))

# ---- 9 grids: ids unique, seeds not shared with any other grid, libpilot3 is fresh
ranges <- function(g) data.frame(cell = g$cell, a = g$seed0 + 1, b = g$seed0 + g$reps)
nx <- do.call(rbind, lapply(V3_EXPS, function(x) ranges(make_grid(x))))
ox <- do.call(rbind, lapply(c(paste0("E", 1:8), "calib2", "libpilot", "libpilot2", "libsize", "libweak",
                              "smoke", "quick", "E0", "Cdiag", "nndiag", "nndiag2"),
                            function(x) ranges(make_grid(x))))
hit <- 0L
for (i in seq_len(nrow(nx))) {
  oo <- rbind(nx[-i, ], ox); oo <- oo[oo$a <= nx$b[i] & oo$b >= nx$a[i], ]
  if (nrow(oo)) hit <- hit + 1L
}
chk(hit == 0L, sprintf("9 %d cells of libpilot3, E11-E13: no seed range overlaps another grid", nrow(nx)))

# ---- 10 oracle AUC and the new generators (need R/calib_params_s.R)
# This closed-form check must run even before calibration exists.
ab <- oracle_auc_old("B", 0.4, list(cB = 3), d = 1L)
chk(abs(ab - 2 / pi * atan(sqrt(1 + 3 * 0.4))) < 0.002,
    sprintf("10pre model B oracle AUC %.4f agrees with 1D closed form", ab))
if (!is.null(S_CALIB)) {
  # (M1): the oracle equals Phi(Delta / sqrt(2))
  mu <- 0.6 * 0.8 * (1:100)^(-3); Delta <- sqrt(drop(t(mu) %*% solve(ar1_cov(100)) %*% mu))
  chk(abs(oracle_auc_old("A", 0.8, list(cA = 0.6), d = 100L) - pnorm(Delta / sqrt(2))) < 0.002,
      sprintf("10a oracle AUC of (M1) %.4f vs Phi(Delta/sqrt 2) = %.4f",
              oracle_auc_old("A", 0.8, list(cA = 0.6), d = 100L), pnorm(Delta / sqrt(2))))
  # the calibrated parameters give the target ideal AUC
  nz1 <- oracle_noise(S1_NDEV)
  a1 <- oracle_S1a(s_param("S1a", 0.10), nz1); a2 <- oracle_S1b(s_param("S1b", 0.10), nz1)
  nz2 <- oracle_noise(20L); a3 <- oracle_S2(s_param("S2", 0.10), nz2)
  chk(all(abs(c(a1, a2, a3) - 0.6) < 0.003),
      sprintf("10b calibrated ideal AUC at 0.60: S1a %.4f S1b %.4f S2 %.4f %.4f %.4f", a1, a2, a3[1], a3[2], a3[3]))
  # mean 0, identity covariance of the deviant populations (large sample)
  set.seed(1010)
  chkcov <- function(X, nm) {
    ok <- max(abs(colMeans(X))) < 0.02 && max(abs(cov(X) - diag(ncol(X)))) < 0.03
    chk(ok, sprintf("10c %s: mean 0 and identity covariance (max |mean| %.3f, max |cov - I| %.3f)",
                    nm, max(abs(colMeans(X))), max(abs(cov(X) - diag(ncol(X))))))
  }
  chkcov(gen_population_s("S1a", 0.10, 4e5L, 6L, 3L, 3L, par0)[, 1:3], "S1a deviant, coordinates 1-3")
  Xt <- gen_population_s("S1b", 0.10, 4e5L, 6L, 3L, 3L, par0)[, 1:3]
  nu <- 1 / s_param("S1b", 0.10); cc <- sqrt(nu / (nu - 2))
  chk(abs(mean(abs(Xt) < 1) - (2 * pt(cc, nu) - 1)) < 0.005 && max(abs(colMeans(Xt))) < 0.02,
      sprintf("10c S1b deviant: standardized t(%.2f) scale (P(|X|<1) %.4f vs %.4f)", nu,
              mean(abs(Xt) < 1), 2 * pt(cc, nu) - 1))
  chkcov(gen_population_s("S2", 0.10, 4e5L, 22L, 4L, 4L, par0)[, 11:20], "S2 population 4, pair block")
  chkcov(gen_population_s("S3b", 0.10, 4e5L, 20L, 10L, 10L, par0), "S3b deviant")
  # S2: pop 2 shifted mean, pop 3 changed variance, others N(0, I)
  th <- s_param("S2", 0.10)
  X2 <- gen_population_s("S2", 0.10, 4e5L, 22L, 2L, 4L, par0); X3 <- gen_population_s("S2", 0.10, 4e5L, 22L, 3L, 4L, par0)
  chk(abs(mean(X2[, 1:5]) - th[1]) < 0.01 && abs(mean(X3[, 6:10]^2) - (1 + th[2])) < 0.02 &&
        abs(mean(X2[, 6:22])) < 0.01,
      "10d S2 populations 2 and 3: mean shift on coordinates 1-5, variance 1 + gamma on 6-10")
} else chk(FALSE, "10 calibration missing: full validation is incomplete")

cat(if (all_ok) "\nPOWER3 TESTS OK\n" else "\nPOWER3 TESTS FAILED\n")
# marker for tests/retry_r.bat (R sometimes exits silently on this machine)
writeLines(if (all_ok) "OK" else "FAILED", file.path("tests", "out", "test_power3.done"))
if (!all_ok) quit(status = 1L)
