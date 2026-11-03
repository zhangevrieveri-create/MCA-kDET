# =============================================================================
# E8: size under serial dependence (GARCH null), Section 3.1 "Dependent
# observations" of the draft and the financial application of Section 5.
#
# One replication generates a T x d series of weekly returns from a
# constant-conditional-correlation GARCH(1,1) with standardized t innovations.
# The regimes come from an external indicator (a k-state Markov chain) that is
# independent of the returns and fixed in advance: it is drawn once per cell
# from seed0 and is the same in every replication. Hence the k regimes have
# the same marginal distribution (H_0 holds) but the observations are serially
# dependent.
#
# Four versions of MCA-kDET are applied to the same replication:
#   flt_blk  coordinates standardized by univariate GARCH(1,1) filters
#            (Gaussian QML, fitted on the whole series, label-blind) and a
#            split into contiguous blocks with a buffer: the procedure of
#            Section 5. Its p-values carry the usual names (mca, mca_l1, ...).
#   flt_iid  filtered, random within-regime split (Step 1 of Algorithm 1)
#   raw_blk  unfiltered returns, block split with buffer
#   raw_iid  unfiltered returns, random within-regime split
# DISCO, MMD and MST are applied to the filtered data (disco, mmd, mst) and
# to the raw data (raw_disco, raw_mmd, raw_mst).
#
# Model parameters come from the cell column mpar ("key=value;..."), set in
# make_grid("E8"); see GARCH_DEFAULT for the keys.
# R/models.R and the permutation engine in R/tests.R are not touched.
# =============================================================================

GARCH_DEFAULT <- list(
  alpha = 0.10, beta = 0.85,   # GARCH(1,1); omega = 1 - alpha - beta (unit variance)
  nu = 5,                      # degrees of freedom of the t innovations (Inf = Gaussian)
  rho = 0.3,                   # constant conditional correlation (equicorrelation)
  T = 900,                     # weeks
  burn = 500,                  # burn-in weeks, discarded
  dur = 26,                    # mean duration (weeks) of a regime episode
  minreg = 150,                # every regime has at least this many weeks
  block = 26,                  # block length (weeks) of the block split
  buf = 4,                     # buffer (weeks): test weeks this close to a fit week are dropped
  minobs = 30                  # every regime has >= minobs fit and test weeks
)

parse_mpar <- function(s, default = GARCH_DEFAULT) {
  out <- default
  if (is.null(s) || length(s) == 0L || is.na(s) || !nzchar(s)) return(out)
  for (kv in strsplit(strsplit(s, ";", fixed = TRUE)[[1L]], "=", fixed = TRUE)) {
    if (!kv[1L] %in% names(default)) stop("unknown mpar key: ", kv[1L])
    out[[kv[1L]]] <- as.numeric(kv[2L])
  }
  out
}

# k-state Markov chain: stays with probability 1 - 1/dur, otherwise moves to
# one of the other states uniformly. Redrawn until every state has >= minreg
# weeks. Uses the current random-number stream.
gen_regime_path <- function(Tn, k, dur, minreg, tries = 1000L) {
  for (a in seq_len(tries)) {
    s <- integer(Tn)
    s[1L] <- sample.int(k, 1L)
    sw <- runif(Tn) < 1 / dur
    for (t in 2:Tn) {
      if (sw[t]) {
        o <- setdiff(seq_len(k), s[t - 1L])
        s[t] <- o[sample.int(length(o), 1L)]
      } else s[t] <- s[t - 1L]
    }
    if (all(tabulate(s, k) >= minreg)) return(s)
  }
  stop("gen_regime_path: no path with >= ", minreg, " weeks per regime")
}

# CCC-GARCH(1,1): x_tj = sigma_tj z_tj, sigma^2_tj = omega + alpha x^2_{t-1,j}
# + beta sigma^2_{t-1,j}, z_t standardized multivariate t (correlation R).
gen_ccc_garch <- function(Tn, d, alpha, beta, nu, rho, burn) {
  N <- Tn + burn
  R <- matrix(rho, d, d); diag(R) <- 1
  E <- matrix(rnorm(N * d), N, d) %*% chol(R)
  if (is.finite(nu)) E <- E * sqrt((nu - 2) / rchisq(N, nu))   # unit variance
  om <- 1 - alpha - beta
  X <- matrix(0, N, d)
  s2 <- rep(1, d); x2 <- rep(1, d)
  for (t in seq_len(N)) {
    s2 <- om + alpha * x2 + beta * s2
    X[t, ] <- sqrt(s2) * E[t, ]
    x2 <- X[t, ]^2
  }
  X[(burn + 1L):N, , drop = FALSE]
}

# Univariate GARCH(1,1) by Gaussian quasi-maximum likelihood; returns the
# standardized residuals of the demeaned series. Parametrization:
# omega = exp(th1), persistence p = plogis(th2), alpha = p plogis(th3),
# beta = p - alpha. sigma^2 starts at the sample variance.
garch11_filter <- function(x) {
  x <- x - mean(x)
  Tn <- length(x); x2 <- x^2; v0 <- mean(x2)
  s2_of <- function(th) {
    om <- exp(th[1L]); p <- plogis(th[2L]); a <- p * plogis(th[3L]); b <- p - a
    u <- om + a * c(v0, x2[-Tn])                       # omega + alpha x^2_{t-1}
    as.numeric(stats::filter(u, b, method = "recursive", init = v0))
  }
  nll <- function(th) {
    s2 <- s2_of(th)
    if (any(!is.finite(s2)) || any(s2 <= 0)) return(1e100)
    0.5 * sum(log(s2) + x2 / s2)
  }
  th0 <- c(log(0.05 * v0), qlogis(0.95), qlogis(0.1 / 0.95))
  o <- optim(th0, nll, method = "Nelder-Mead", control = list(maxit = 2000L, reltol = 1e-10))
  s2 <- s2_of(o$par)
  list(z = x / sqrt(s2), persist = plogis(o$par[2L]), conv = o$convergence)
}

# Random within-regime split (Step 1 of Algorithm 1), as a role vector.
iid_roles <- function(reg, k, tau) {
  role <- character(length(reg))
  for (i in seq_len(k)) {
    ix <- which(reg == i)
    ix <- ix[sample.int(length(ix))]
    nf <- floor(tau * length(ix))
    role[ix[seq_len(nf)]] <- "fit"
    role[ix[(nf + 1L):length(ix)]] <- "test"
  }
  role
}

# Block split: the weeks are cut into contiguous blocks of length `block`; a
# random set of round(tau * #blocks) blocks goes to D_fit, the rest to D_test;
# a test week within `buf` weeks of a fit week is dropped. Redrawn (independently
# of the returns) until every regime has >= minobs fit and test weeks.
block_roles <- function(reg, k, tau, block, buf, minobs, tries = 1000L) {
  Tn <- length(reg)
  bid <- ceiling(seq_len(Tn) / block); nb <- max(bid)
  for (a in seq_len(tries)) {
    isfit <- logical(nb); isfit[sample.int(nb, round(tau * nb))] <- TRUE
    role <- ifelse(isfit[bid], "fit", "test")
    if (buf > 0) {
      cs <- c(0L, cumsum(role == "fit"))
      lo <- pmax(1L, seq_len(Tn) - buf); hi <- pmin(Tn, seq_len(Tn) + buf)
      near <- (cs[hi + 1L] - cs[lo]) > 0
      role[role == "test" & near] <- "drop"
    }
    if (all(tabulate(reg[role == "fit"], k) >= minobs) &&
        all(tabulate(reg[role == "test"], k) >= minobs)) return(role)
  }
  stop("block_roles: no admissible block split")
}

# the list(fit, test) per regime that fit_scores() expects
roles_to_split <- function(M, reg, role, k)
  lapply(seq_len(k), function(i)
    list(fit = M[reg == i & role == "fit", , drop = FALSE],
         test = M[reg == i & role == "test", , drop = FALSE]))

sq_acf1 <- function(M) mean(apply(M^2, 2L, function(v) cor(v[-1L], v[-length(v)])))

one_rep_garch <- function(cell, r) {
  mp <- parse_mpar(cell$mpar)
  k <- as.integer(cell$k); d <- as.integer(cell$d); Tn <- as.integer(mp$T)
  B <- as.integer(cell$B); tau <- as.numeric(cell$tau)
  lib <- parse_lib(cell$lib)
  # regime path: fixed in advance, the same in every replication of the cell
  set.seed(as.integer(cell$seed0))
  reg <- gen_regime_path(Tn, k, mp$dur, mp$minreg)
  seed <- as.integer(cell$seed0) + as.integer(r)
  set.seed(seed)
  t_all <- proc.time()[["elapsed"]]

  X <- gen_ccc_garch(Tn, d, mp$alpha, mp$beta, mp$nu, mp$rho, as.integer(mp$burn))
  t_f <- system.time(fl <- lapply(seq_len(d), function(j) garch11_filter(X[, j])))[["elapsed"]]
  Z <- do.call(cbind, lapply(fl, `[[`, "z"))

  rb <- block_roles(reg, k, tau, as.integer(mp$block), as.integer(mp$buf), as.integer(mp$minobs))
  ri <- iid_roles(reg, k, tau)

  row <- list(exp = cell$exp, cell = cell$cell, rep = r, seed = seed, err = "")
  variants <- list(flt_blk = list(Z, rb), flt_iid = list(Z, ri),
                   raw_blk = list(X, rb), raw_iid = list(X, ri))
  for (v in names(variants)) {
    sp <- roles_to_split(variants[[v]][[1L]], reg, variants[[v]][[2L]], k)
    sc <- fit_scores(sp, "ovr", lib = lib)
    ct <- classifier_tests(sc, B, libsets = NULL)
    pre <- if (v == "flt_blk") "" else paste0(v, "_")
    names(ct$p) <- paste0(pre, names(ct$p))
    dg <- ct$diag[grep("^(auc_.+_dev|t_perm)$", names(ct$diag))]
    names(dg) <- paste0(pre, names(dg))
    row <- c(row, ct$p, dg)
    if (v == "flt_blk") {
      row <- c(row, ct$diag[setdiff(names(ct$diag), c(grep("^auc_.+_dev$", names(ct$diag), value = TRUE), "t_perm"))],
               setNames(sc$tm, paste0("t_", names(sc$tm))))
    }
  }
  cf <- competitor_tests(lapply(seq_len(k), function(i) Z[reg == i, , drop = FALSE]), B, "dist")
  cr <- competitor_tests(lapply(seq_len(k), function(i) X[reg == i, , drop = FALSE]), B, "dist")
  row <- c(row, cf$p, setNames(cr$p, paste0("raw_", names(cr$p))))

  # diagnostics
  row$n_fit_blk <- sum(rb == "fit"); row$n_test_blk <- sum(rb == "test")
  row$n_drop_blk <- sum(rb == "drop")
  row$sqacf_raw <- sq_acf1(X); row$sqacf_flt <- sq_acf1(Z)
  row$persist_hat <- mean(vapply(fl, `[[`, numeric(1), "persist"))
  row$filter_nonconv <- sum(vapply(fl, `[[`, numeric(1), "conv") != 0)
  row$t_filter <- t_f
  row$t_total <- proc.time()[["elapsed"]] - t_all
  row
}
