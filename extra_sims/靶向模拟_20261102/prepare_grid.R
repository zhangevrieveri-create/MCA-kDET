# Fixed grid of the targeted simulation (written before any formal run).
# Run from sim/:  Rscript 靶向模拟_20261102/prepare_grid.R
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"; source("R/experiments.R")
folder <- "靶向模拟_20261102"
LIB <- paste0(LIBPILOT3_LIB, ",deprf")
SETS <- paste("l1+quad+knn+rf",                      # L0, paper default
              "l1+quad+knn+rf+deprf",                # L0P, primary new library
              "xscr1+xscr3+knn+rf+g5f",              # L6
              "xscr1+xscr3+knn+rf+g5f+deprf", sep = ";")
GROUPS <- "ovr,dist,schott,comp2,comp2ovr,songchen"   # + adk_max, computed in common.R

specs <- list()
add <- function(...) specs[[length(specs) + 1L]] <<- list(...)
# T1: dimension, k = 3, n = 300
for (t in c(0.65, 0.75)) for (d in c(20L, 50L, 100L, 200L)) add(part = "T1", model = "SM3", target = t, k = 3L, n = 300L, d = d)
for (d in c(20L, 100L)) add(part = "T1", model = "C3", target = 0.8, k = 3L, n = 300L, d = d)   # target = beta of (M3)
for (d in c(20L, 50L, 100L, 200L)) add(part = "T1", model = "SM3", target = 0.5, k = 3L, n = 300L, d = d)
# T2: number of populations, d = 20, n = 100 per population
for (k in c(3L, 5L, 10L, 20L)) add(part = "T2", model = "SM3", target = 0.70, k = k, n = 100L, d = 20L)
for (k in c(10L, 20L)) add(part = "T2", model = "SM3", target = 0.5, k = k, n = 100L, d = 20L)

grid <- do.call(rbind, lapply(seq_along(specs), function(i) {
  s <- specs[[i]]
  null <- s$model == "SM3" && s$target == 0.5
  beta <- if (s$model == "SM3") s$target - 0.5 else s$target
  g <- cell_row("tgt", s$model, beta, k = s$k, n = s$n, d = s$d, B = 200L,
                reps = if (null) 500L else 100L, groups = GROUPS, lib = LIB, libsets = SETS,
                seed_scheme = "v2", tag = paste0("_", s$part))
  g$target <- s$target; g$part <- s$part
  g$seed0 <- 165000000L + i * 10000L
  g
}))
stopifnot(nrow(grid) == 20L, !anyDuplicated(grid$cell), !anyDuplicated(grid$seed0),
          min(grid$seed0) > 160300000L, max(grid$seed0) + 3000L < 200000000L,
          all(grid$reps <= 2000L))
path <- file.path(folder, "grid.rds")
if (file.exists(path)) stopifnot(identical(grid, readRDS(path))) else saveRDS(grid, path)
csv <- file.path(folder, "grid.csv")
if (!file.exists(csv)) write.csv(grid, csv, row.names = FALSE)
cat("FIXED GRID:", nrow(grid), "cells;", sum(grid$reps), "replications; B = 200\n")
print(grid[, c("cell", "model", "target", "k", "n", "d", "reps", "seed0")], row.names = FALSE)
