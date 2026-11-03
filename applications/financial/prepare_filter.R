# Pre-2015 GARCH(1,1) QML fit; forward recursion without using future returns.
argv <- commandArgs(FALSE)
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)[1]), winslash = "/"))
input <- file.path(root, "prepared")
X <- read.csv(file.path(input, "weekly_returns.csv"), check.names = FALSE)
weeks <- read.csv(file.path(input, "weeks.csv"), check.names = FALSE)
stopifnot(identical(X$week_end, weeks$week_end))
M <- as.matrix(X[, -1L, drop = FALSE]); baseline <- weeks$period == "baseline"
recursion <- function(x, omega, alpha, beta, v0) {
  s <- numeric(length(x)); previous_h <- v0; previous_x2 <- v0
  for (t in seq_along(x)) {
    s[t] <- omega + alpha * previous_x2 + beta * previous_h
    previous_h <- s[t]
    previous_x2 <- if (is.finite(x[t])) x[t]^2 else s[t]
  }
  s
}
fit_one <- function(x) {
  mu <- mean(x[baseline], na.rm = TRUE); centered <- x - mu
  xfit <- centered[baseline]; v0 <- mean(xfit^2, na.rm = TRUE)
  decode <- function(th) {
    persistence <- plogis(th[2L]); alpha <- persistence * plogis(th[3L])
    c(omega = exp(th[1L]), alpha = alpha, beta = persistence - alpha)
  }
  objective <- function(th) {
    p <- decode(th)
    if (any(!is.finite(p))) return(1e100)
    h <- recursion(xfit, p[1L], p[2L], p[3L], v0)
    if (any(!is.finite(h)) || any(h <= 0)) return(1e100)
    good <- is.finite(xfit)
    .5 * sum(log(h[good]) + xfit[good]^2 / h[good])
  }
  starts <- c(.85, .95, .985)
  fits <- lapply(starts, function(p) optim(c(log((1 - p) * v0), qlogis(p), qlogis(.1 / p)),
                                           objective, method = "Nelder-Mead",
                                           control = list(maxit = 4000L, reltol = 1e-10)))
  converged <- which(vapply(fits, function(o) o$convergence == 0L && is.finite(o$value), logical(1)))
  if (!length(converged)) stop("No converged presample GARCH fit")
  best <- fits[[converged[which.min(vapply(fits[converged], `[[`, numeric(1), "value"))]]]
  p <- decode(best$par)
  h <- recursion(centered, p[1L], p[2L], p[3L], v0)
  list(z = centered / sqrt(h), h = h,
       parameters = c(mean = mu, p, persistence = p[2L] + p[3L], initial_variance = v0,
                      convergence = best$convergence, negative_log_likelihood = best$value,
                      baseline_observations = sum(is.finite(xfit))))
}
fit <- vector("list", ncol(M))
for (j in seq_len(ncol(M))) {
  fit[[j]] <- fit_one(M[, j])
  cat("Fitted GARCH ", j, "/", ncol(M), ": ", colnames(M)[j], "\n", sep = ""); flush.console()
}
Z <- do.call(cbind, lapply(fit, `[[`, "z")); H <- do.call(cbind, lapply(fit, `[[`, "h"))
colnames(Z) <- colnames(H) <- colnames(M)
parameters <- data.frame(stock = colnames(M), do.call(rbind, lapply(fit, `[[`, "parameters")), row.names = NULL)
names(parameters) <- sub("\\.alpha$", "", names(parameters))
write.csv(parameters, file.path(input, "garch_parameters.csv"), row.names = FALSE)
write_exact <- function(mat, filename) {
  tab <- data.frame(week_end = weeks$week_end, check.names = FALSE)
  for (j in seq_len(ncol(mat))) tab[[colnames(mat)[j]]] <- ifelse(is.finite(mat[, j]), sprintf("%.17g", mat[, j]), NA_character_)
  write.csv(tab, file.path(input, filename), row.names = FALSE, na = "")
}
write_exact(Z, "filtered_returns.csv"); write_exact(H, "conditional_variance.csv")
acf_pair <- function(v, lag) {
  u <- v[seq_len(length(v) - lag)]; w <- v[seq.int(lag + 1L, length(v))]
  good <- is.finite(u) & is.finite(w)
  if (sum(good) < 5L) return(NA_real_)
  cor(u[good], w[good])
}
diagnostics <- list(); z <- 1L
for (j in seq_len(ncol(M))) for (representation in c("raw", "filtered")) for (power in 1:2) for (lag in 1:4) {
  v <- if (representation == "raw") M[weeks$period == "analysis", j] else Z[weeks$period == "analysis", j]
  diagnostics[[z]] <- data.frame(stock = colnames(M)[j], representation = representation,
                                  power = power, lag = lag, autocorrelation = acf_pair(v^power, lag))
  z <- z + 1L
}
write.csv(do.call(rbind, diagnostics), file.path(input, "autocorrelation.csv"), row.names = FALSE)
saveRDS(list(fit = fit, parameters = parameters, weeks = weeks, X = M, Z = Z, H = H,
             baseline = baseline, session = sessionInfo()), file.path(input, "filters.rds"))
cat("Pre-2015-only GARCH filters completed.\n")
