# Application helpers; frozen sim/R learners and statistic engine are not modified.
install_panel_knn <- function() {
  original <- knn_scores
  cache <- new.env(parent = emptyenv())
  fast <- function(X, y, w, Xnew, Ks) {
    Kmax <- max(Ks)
    hit <- exists("X", envir = cache, inherits = FALSE) &&
      identical(X, cache$X) && identical(Xnew, cache$Xnew) && identical(Kmax, cache$Kmax)
    if (hit) NB <- cache$NB else {
      D2 <- outer(rowSums(Xnew^2), rowSums(X^2), "+") - 2 * tcrossprod(Xnew, X)
      NB <- matrix(.C("mca_panel_neighbors", as.double(D2), as.integer(nrow(D2)),
                       as.integer(ncol(D2)), as.integer(Kmax),
                       answer = integer(nrow(D2) * Kmax))$answer, nrow(D2))
      if (nrow(Xnew) > nrow(X)) {
        cache$X <- X; cache$Xnew <- Xnew; cache$Kmax <- Kmax; cache$NB <- NB
      }
    }
    WN <- matrix(w[NB], nrow(NB)); YN <- matrix(y[NB], nrow(NB))
    vapply(Ks, function(K) {
      kap <- (K - seq_len(K) + 1) / K
      a <- WN[, seq_len(K), drop = FALSE] %*% kap
      b <- (WN[, seq_len(K), drop = FALSE] * YN[, seq_len(K), drop = FALSE]) %*% kap
      as.numeric(b / a)
    }, numeric(nrow(Xnew)))
  }
  set.seed(918)
  for (nnew in c(17L, 100L)) for (tied in c(FALSE, TRUE)) {
    X <- matrix(rnorm(60L * 5L), 60L); new <- matrix(rnorm(nnew * 5L), nnew)
    if (tied) {X <- round(X); new <- round(new)}
    y <- rep(0:1, 30L); w <- seq(.4, 1.6, length.out = 60L)
    expected <- original(X, y, w, new, c(1L, 5L, 13L))
    actual <- fast(X, y, w, new, c(1L, 5L, 13L))
    stopifnot(identical(actual, expected))
    stopifnot(identical(fast(X, 1L - y, rev(w), new, c(1L, 5L, 13L)),
                        original(X, 1L - y, rev(w), new, c(1L, 5L, 13L))))
  }
  assign("knn_scores", fast, envir = .GlobalEnv)
  invisible(TRUE)
}

panel_auc <- function(S, groups, targets, time_id, counts) {
  stopifnot(nrow(S) == length(groups), length(time_id) == nrow(S),
            ncol(S) == length(targets), all(is.finite(S)), min(time_id) >= 1L,
            max(time_id) <= nrow(counts), all(counts >= 0))
  result <- .C("mca_panel_auc", as.double(S), as.integer(nrow(S)), as.integer(ncol(S)),
               as.integer(groups), as.integer(targets), as.integer(time_id - 1L),
               as.integer(nrow(counts)), as.integer(ncol(counts)), as.integer(counts),
               answer = double(ncol(S) * ncol(counts)))$answer
  matrix(result, nrow = ncol(S))
}

panel_counts <- function(Tn, length, B, seed) {
  set.seed(seed)
  out <- matrix(0L, Tn, B + 1L); out[, 1L] <- 1L
  for (b in seq_len(B)) {
    starts <- sample.int(Tn - length + 1L, ceiling(Tn / length), replace = TRUE)
    indices <- as.vector(t(outer(starts, seq.int(0L, length - 1L), "+")))[seq_len(Tn)]
    out[, b + 1L] <- tabulate(indices, Tn)
  }
  out
}

panel_question <- function(fitting, testing, fit_valid, test_valid, members, groups, seed) {
  # Array axes: calendar block, stock, feature. Rows group -> stock -> calendar.
  k <- max(groups); d <- dim(fitting)[3L]
  data_for <- function(A, valid, stocks) do.call(rbind, lapply(stocks, function(j)
    matrix(A[valid, j, , drop = FALSE], ncol = d)))
  sp <- lapply(seq_len(k), function(g) {
    ids <- members[groups == g]
    list(fit = data_for(fitting, fit_valid, ids), test = data_for(testing, test_valid, ids))
  })
  set.seed(seed)
  warnings_seen <- character()
  sc <- withCallingHandlers(fit_scores(sp, "ovr", lib = c("l1", "quad", "knn", "rf")),
    warning = function(w) {warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning")})
  S <- do.call(cbind, sc$Sovr)
  targets <- rep(seq_len(k), length(sc$lib))
  time_id <- rep(which(test_valid), length(members))
  mi <- sc$m[targets]; M <- sum(sc$m)
  scale <- sqrt((M + 1) / (12 * mi * (M - mi)))
  list(S = S, groups = sc$L, targets = targets, time_id = time_id,
       scale = scale, fit_n = sum(vapply(sp, function(s) nrow(s$fit), integer(1))),
       test_n = M, warnings = unique(warnings_seen), tuning = sc$info_ovr, timing = sc$tm,
       nonconv = sum(unlist(lapply(sc$info_ovr, function(o) o$nonconv))))
}

panel_calibrate <- function(question, counts) {
  A <- panel_auc(question$S, question$groups, question$targets, question$time_id, counts)
  if (!all(is.finite(A))) stop("A bootstrap draw has no complete panel")
  observed <- (A[, 1L] - .5) / question$scale
  centered <- sweep(sweep(A[, -1L, drop = FALSE], 1L, A[, 1L], "-"), 1L, question$scale, "/")
  boot <- apply(centered, 2L, max)
  statistic <- max(observed)
  p <- (1 + sum(boot >= statistic - 1e-10)) / ncol(counts)
  # Frozen engine check at every replication, using original observed scores.
  engine <- ovr_T(question$S, question$targets, matrix(question$groups, ncol = 1L))[, 1L]
  stopifnot(max(abs(engine - observed)) < 1e-10)
  list(p = p, statistic = statistic, observed_auc = A[, 1L], bootstrap_max = boot,
       critical95 = as.numeric(quantile(boot, .95, type = 1)), auc_matrix = A)
}

panel_generate <- function(scenario, rep, fit_blocks, test_blocks, d, seed) {
  # Every coordinate has the same marginal time-process under every listed null.
  set.seed(seed); J <- 15L; burn <- 2000L; gap <- ceiling(20L / d)
  blocks <- fit_blocks + gap + test_blocks; days <- blocks * d
  N <- days + burn
  loading <- seq(.15, .65, length.out = J)
  sector <- rep(seq_len(3L), each = 5L); industry_loading <- .4
  eps <- matrix(rnorm(N * J), N, J) * rep(sqrt(1 - loading^2 - industry_loading^2), each = N)
  common <- rnorm(N); industry <- matrix(rnorm(N * 3L), N, 3L)
  eps <- eps + outer(common, loading) + industry_loading * industry[, sector]
  if (scenario == "null_t5_ar") eps <- eps * sqrt(3 / rchisq(N, 5))
  phi <- switch(scenario, null_ar50 = .5, null_ar90 = .9, null_t5_ar = .5, 0)
  phis <- rep(phi, J)
  if (scenario == "alt_sector_dependence") phis <- c(0, .4, .7)[sector]
  X <- matrix(0, N, J)
  if (scenario %in% c("null_garch95", "null_garch99")) {
    a <- .1; b <- if (scenario == "null_garch95") .85 else .89
    h <- previous <- rep(1, J)
    for (t in seq_len(N)) {
      h <- 1 - a - b + a * previous + b * h
      X[t, ] <- sqrt(h) * eps[t, ]; previous <- X[t, ]^2
    }
  } else {
    previous <- rep(0, J)
    for (t in seq_len(N)) {
      X[t, ] <- phis * previous + sqrt(1 - phis^2) * eps[t, ]
      previous <- X[t, ]
    }
  }
  X <- X[seq.int(burn + 1L, N), , drop = FALSE]
  if (scenario == "alt_sector_mean") X <- X + rep(c(-.2, 0, .2)[sector], each = days)
  if (scenario == "alt_sector_scale") X <- X * rep(c(.7, 1, 1.4)[sector], each = days)
  if (scenario == "alt_one_bank") X[, 1L] <- X[, 1L] + .35
  A <- aperm(array(X, c(d, blocks, J)), c(2L, 3L, 1L))
  list(fit = A[seq_len(fit_blocks), , , drop = FALSE],
       test = A[seq.int(fit_blocks + gap + 1L, blocks), , , drop = FALSE])
}

panel_rep <- function(scenario, r, d, B, lengths, masks, seed0, keep = FALSE) {
  id <- match(scenario, c("null_iid", "null_ar50", "null_ar90", "null_garch95", "null_garch99",
                          "null_t5_ar", "alt_sector_mean", "alt_sector_scale", "alt_sector_dependence", "alt_one_bank"))
  seed <- as.integer(seed0 + id * 1000000L + r * 100L)
  started <- proc.time()[["elapsed"]]
  generated <- panel_generate(scenario, r, length(masks$fit), length(masks$test), d, seed)
  qnames <- c("within_bank", "within_pharma", "within_machinery", "between_sectors")
  questions <- lapply(seq_along(qnames), function(q) {
    members <- if (q <= 3L) seq.int((q - 1L) * 5L + 1L, q * 5L) else seq_len(15L)
    groups <- if (q <= 3L) seq_len(5L) else rep(seq_len(3L), each = 5L)
    panel_question(generated$fit, generated$test, masks$fit, masks$test, members, groups, seed + q)
  })
  rows <- list(); retained <- list(); z <- 1L
  for (block_len in lengths) {
    counts <- panel_counts(length(masks$test), block_len, B, seed + 50L + block_len)
    calibrated <- lapply(questions, panel_calibrate, counts = counts)
    ps <- vapply(calibrated, `[[`, numeric(1), "p"); adjusted <- p.adjust(ps, "holm")
    for (q in seq_along(qnames)) {
      o <- calibrated[[q]]
      rows[[z]] <- data.frame(scenario = scenario, rep = r, seed = seed, d = d, block_length = block_len,
                              B = B, question = qnames[q], p = o$p, p_holm = adjusted[q],
                              statistic = o$statistic, critical95 = o$critical95,
                              fitting_complete_blocks = sum(masks$fit), test_complete_blocks = sum(masks$test),
                              fitting_rows = questions[[q]]$fit_n, testing_rows = questions[[q]]$test_n,
                              learner_warning_count = length(questions[[q]]$warnings),
                              fit_nonconvergence_flags = questions[[q]]$nonconv,
                              elapsed_seconds = proc.time()[["elapsed"]] - started)
      z <- z + 1L
    }
    if (keep) retained[[as.character(block_len)]] <- list(counts = counts, calibrated = calibrated)
  }
  list(rows = do.call(rbind, rows), audit = if (keep) list(questions = questions, retained = retained,
                                                         scenario = scenario, rep = r, seed = seed) else NULL)
}

panel_run_job <- function(job, config, masks) {
  destination <- file.path(config$outdir, "chunks", sprintf("%s__r%05d-%05d.csv", job$scenario, job$start, job$end))
  if (file.exists(destination)) return(list(file = destination, skipped = TRUE))
  rows <- list()
  for (r in seq.int(job$start, job$end)) {
    attempt <- NULL
    for (a in seq_len(3L)) {
      attempt <- tryCatch(panel_rep(job$scenario, r, config$d, config$B, config$lengths,
                                     masks, config$seed0, keep = r <= config$audit_reps), error = function(e) e)
      if (!inherits(attempt, "error")) break
    }
    if (inherits(attempt, "error")) {
      saveRDS(list(job = job, rep = r, message = conditionMessage(attempt)),
              file.path(config$outdir, "audits", sprintf("error_%s_r%05d.rds", job$scenario, r)))
      stop(job$scenario, " rep ", r, ": ", conditionMessage(attempt))
    }
    rows[[length(rows) + 1L]] <- attempt$rows
    if (!is.null(attempt$audit)) saveRDS(attempt$audit, file.path(config$outdir, "audits",
                                                sprintf("%s_r%05d.rds", job$scenario, r)))
  }
  tab <- do.call(rbind, rows)
  part <- paste0(destination, ".part")
  write.csv(tab, part, row.names = FALSE)
  if (!file.rename(part, destination)) stop("Could not publish completed chunk")
  cat(sprintf("Completed %s reps %d-%d\n", job$scenario, job$start, job$end)); flush.console()
  list(file = destination, skipped = FALSE)
}
