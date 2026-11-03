# =============================================================================
# Seed-range check: no two cells may share a seed (seed0 + 1 .. seed0 + reps).
# Checks every revised grid (E1-E8 in the FULL version, calib2)
# against each other and against all earlier grids. Writes nothing.
#
#   Rscript tests/check_seeds.R                       # uses R/calib_params_v2.R
#   Rscript tests/check_seeds.R 0.1 0.5 0.8           # trial values: max of the
#                                                     # D grid, max of the C2 grid,
#                                                     # beta_fix_C (grids = 0.2-1.0 x max)
# Re-run after filling in calib_params_v2.R: the ids, hence the seeds, of the
# D, C2 and C' cells depend on the calibrated betas.
# =============================================================================
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- "R"; PROFILE <- "full"
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "competitors2.R",
            "rotation.R", "learners_lib2.R", "models_s.R", "experiments.R"))
  source(file.path("R", f))
# only the cell ids, hence the seeds, matter here; the library of E11-E13 is chosen later
if (is.na(LIB_MAIN_V3)) LIB_MAIN_V3 <- "l1,quad,knn,rf"
a <- as.numeric(commandArgs(TRUE))
if (length(a) == 3L) {
  f <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  PAR2$beta_grid_D <- a[1] * f; PAR2$beta_grid_C2 <- a[2] * f; PAR2$beta_fix_C <- a[3]
  cat("trial calibration values:", a, "\n")
} else {
  cat("calibration values from R/calib_params_v2.R: beta_fix_C =", PAR2$beta_fix_C,
      "| C2:", PAR2$beta_grid_C2, "| D:", PAR2$beta_grid_D, "\n")
}

ranges <- function(g) data.frame(cell = g$cell, a = g$seed0 + 1, b = g$seed0 + g$reps)
newx <- c(paste0("E", 1:8), "calib2", V3_EXPS)
oldx <- c("libpilot", "libpilot2", "libsize", "libweak", "smoke", "quick", "E0",
          "Cdiag", "nndiag", "nndiag2")
newg <- do.call(rbind, lapply(newx, function(x) ranges(make_grid(x))))
oldg <- do.call(rbind, lapply(oldx, function(x) ranges(make_grid(x))))
allg <- rbind(newg, oldg)
hit <- 0L
for (i in seq_len(nrow(newg))) {
  o <- allg[allg$cell != newg$cell[i] & allg$a <= newg$b[i] & allg$b >= newg$a[i], ]
  if (nrow(o)) { hit <- hit + 1L; cat("overlap:", newg$cell[i], "<->", unique(o$cell), "\n") }
}
cat(sprintf("revised cells: %d (%d reps in the full version); seed range %.0f - %.0f\n",
            nrow(newg), sum(newg$b - newg$a + 1), min(newg$a), max(newg$b)))
cat(sprintf("earlier cells: %d; seed range %.0f - %.0f\n", nrow(oldg), min(oldg$a), max(oldg$b)))
cat(if (hit == 0L) "SEEDS OK: no overlaps\n" else sprintf("SEEDS FAILED: %d cells overlap\n", hit))

# E9 and E10 are paired with E1/E2 on purpose (same seeds, same data); they
# are not in newx above. Check that every paired cell has exactly the seed0
# and reps of its source cell.
src <- rbind(make_grid("E1"), make_grid("E2"))
src$key <- sub("^E[0-9]+_", "", src$cell)
bad <- 0L
for (x in c("E9", "E10")) {
  g <- make_grid(x)
  m <- match(sub("^E[0-9]+_", "", g$cell), src$key)
  bad <- bad + sum(is.na(m)) + sum(g$seed0 != src$seed0[m] | g$reps != src$reps[m], na.rm = TRUE)
  cat(sprintf("%s: %d cells paired with E1/E2\n", x, nrow(g)))
}
cat(if (bad == 0L) "PAIRED OK: E9/E10 use the seeds of their source cells\n"
    else sprintf("PAIRED FAILED: %d cells\n", bad))
