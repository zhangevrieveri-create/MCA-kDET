# Targeted simulation T1/T2 (2026-11-02). Frozen files in R/ are read only.
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for (p in c("glmnet", "ranger", "gTestsMulti")) stopifnot(requireNamespace(p, quietly = TRUE))
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "garch_ext.R",
            "competitors2.R", "rotation.R", "learners_lib2.R", "models_s.R", "experiments.R"))
  source(file.path(SRC_DIR, f))
TGT_DIR <- "靶向模拟_20261102"
source(file.path(TGT_DIR, "components_dep.R"))     # adk_max (xscrd is registered but not used)
source(file.path(TGT_DIR, "learner_deprf.R"))
TGT_CAL <- readRDS(file.path(TGT_DIR, "calibration.rds"))     # rho by ideal AUC, 3 pairs
TGT_PAIRS <- 1:6                                                  # coordinates of the 3 pairs

tgt_arg <- function(key, default = NULL) {
  a <- grep(paste0("^--", key, "="), commandArgs(TRUE), value = TRUE)
  if (!length(a)) return(default)
  stopifnot(length(a) == 1L); sub(paste0("^--", key, "="), "", a)
}
tgt_manifest_check <- function(out) {
  z <- read.csv(file.path(out, "source_md5.csv"), stringsAsFactors = FALSE)
  actual <- as.character(tools::md5sum(z$file))
  if (any(is.na(actual)) || !identical(actual, z$md5))
    stop("Source fingerprint mismatch: ", paste(z$file[is.na(actual) | actual != z$md5], collapse = ", "))
  invisible(TRUE)
}

# SM3: N(0, I_d) everywhere; population k has sign-mixed correlation rho on the
#      pairs (1,2), (3,4), (5,6) (frozen sign_mix_pairs of R/models_s.R).
#      target = ideal AUC; target 0.5 is the null.
# C3:  paper model (M3) (frozen gen_population "C") with sC = 3 pairs, beta fixed.
tgt_generate <- function(model, target, nvec, d) {
  k <- length(nvec)
  if (model == "SM3") {
    return(lapply(seq_len(k), function(g) {
      X <- matrix(rnorm(nvec[g] * d), nvec[g], d)
      if (g == k && target > 0.5) {
        rho <- TGT_CAL[[sprintf("%.2f", target)]]
        stopifnot(is.finite(rho), rho > 0, rho < 1)
        X <- sign_mix_pairs(X, TGT_PAIRS, rho)
      }
      X
    }))
  }
  if (model == "C3") return(gen_data("C", target, nvec, d, list(sC = 3L)))
  stop("unknown model ", model)
}

tgt_one <- function(cell, r) {
  seed <- as.integer(cell$seed0 + r); set.seed(seed); t0 <- proc.time()[["elapsed"]]
  k <- as.integer(cell$k); B <- as.integer(cell$B); lib <- parse_lib(cell$lib)
  nvec <- parse_n(cell$n, k)
  data <- tgt_generate(cell$model, cell$target, nvec, as.integer(cell$d))
  stopifnot(length(data) == k, identical(vapply(data, nrow, integer(1)), nvec),
            all(vapply(data, ncol, integer(1)) == cell$d))
  sc <- fit_scores(split_data(data, cell$tau), "ovr", lib = lib)
  ct <- classifier_tests(sc, B, parse_libsets(cell$libsets, lib))
  comp <- competitor_tests(data, B, strsplit(cell$groups, ",", fixed = TRUE)[[1L]])
  t_ad <- proc.time()[["elapsed"]]
  ad <- adk_max_test(data, B)                       # last: own permutations
  t_ad <- proc.time()[["elapsed"]] - t_ad
  p <- c(ct$p, comp$p[setdiff(names(comp$p), "kmd_eta")], list(adk_max = ad$p))
  pv <- unlist(p); stopifnot(!anyDuplicated(names(p)), all(is.finite(pv)), all(pv >= 0 & pv <= 1))
  aucs <- do.call(rbind, lapply(names(sc$Sovr), function(a) data.frame(learner = a, population = seq_len(k),
    auc = vapply(seq_len(k), function(i) auc_vec(sc$Sovr[[a]][, i], sc$L == i), numeric(1)))))
  # deprf screening record (deterministic, same inputs as in fit_scores)
  f <- as.integer(table(factor(sc$gfit, levels = seq_len(k)))); omega <- sc$m / sum(sc$m)
  truth <- apply(matrix(TGT_PAIRS, ncol = 2L, byrow = TRUE), 1L, paste, collapse = "-")
  sel <- do.call(rbind, lapply(seq_len(k), function(i) {
    s <- dep_screen(sc$Zfit, as.integer(sc$gfit == i), ovr_weights(sc$gfit, i, omega, f) / mean(ovr_weights(sc$gfit, i, omega, f)))
    key <- paste(s$pairs[, 1L], s$pairs[, 2L], sep = "-")
    data.frame(population = i, pairs = paste(key, collapse = ";"), hits = sum(truth %in% key),
               top_score = s$score[1L])
  }))
  list(cell = cell$cell, model = cell$model, target = cell$target, rep = as.integer(r), seed = seed,
       B = B, k = k, n = cell$n, d = cell$d, tau = cell$tau, lib = cell$lib, libsets = cell$libsets,
       p = p, adk_stat = ad$stat, diag = ct$diag, aucs = aucs, deprf_screen = sel, err = "",
       timing = c(sc$tm, comp$tm, list(adk = t_ad)), seconds = proc.time()[["elapsed"]] - t0)
}
