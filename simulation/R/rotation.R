# =============================================================================
# Rotated models (third review round, item G1; written 1 Nov 2026)
#
# A cell with rot = r > 0 multiplies every observation of every population by
# one fixed d x d orthogonal matrix Q (x -> Q'x, i.e. X %*% Q), drawn from the
# Haar distribution with the seed 910000 + 1000 r + d. Q is the same in every
# replication and every cell with the same (r, d). Q is drawn with its own
# seed and the random-number state is restored afterwards, so the data of
# replication j of a rotated cell are exactly those of replication j of the
# unrotated cell with the same seed, multiplied by Q.
# DISCO, MMD, MST, KMD, the Song-Chen tests and Schott's statistic are
# invariant under x -> Q'x; the max-type tests, the coordinatewise
# standardization and the l1, quad and rf learners are not. The rotation
# removes the alignment of the departure with the coordinates.
# =============================================================================

fixed_rotation <- function(d, key = 1L) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv)
          else rm(".Random.seed", envir = .GlobalEnv))
  set.seed(910000L + 1000L * as.integer(key) + as.integer(d))
  qz <- qr(matrix(rnorm(d * d), d, d))
  qr.Q(qz) %*% diag(sign(diag(qr.R(qz))), nrow = d)
}

rotate_data <- function(data, key = 1L) {
  Q <- fixed_rotation(ncol(data[[1L]]), key)
  lapply(data, function(X) X %*% Q)
}
