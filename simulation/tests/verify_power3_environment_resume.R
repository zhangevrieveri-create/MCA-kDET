# Verify identical p-values before copying ONLY missing successful chunks.
# Run from sim/. Raw chunks already in the destination are never overwritten.
source("R/summarize.R")
old <- "results_trial_power3_progress_fixed_20261101"
fresh <- "results_trial_power3_cpu16_jit0_20261101"
envdir <- "tests/out/environment_20261101"
a <- readRDS(file.path(envdir, "r43_cpu16_jit0_a.rds"))
b <- readRDS(file.path(envdir, "r43_cpu16_jit0_b.rds"))
stopifnot(identical(a$p, b$p), identical(a$seed, b$seed), identical(a$packages, b$packages))
cat("OK independent processes give bit-identical p-values\n")
fs <- list.files(file.path(fresh, "libpilot3"), pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(fs) == 6L)
dest <- file.path(old, "libpilot3", basename(fs))
overlaps <- 0L
for (i in seq_along(fs)) {
  n <- read.csv(fs[i], check.names = FALSE)
  stopifnot(nrow(n) == 1L, all(is.na(n$err) | n$err == ""), n$B == 19, n$model == "S2")
  if (file.exists(dest[i])) {
    o <- read.csv(dest[i], check.names = FALSE)
    pc <- union(intersect(P_COLS, names(n)), grep("^(mca|macc)_", names(n), value = TRUE))
    stopifnot(length(pc) > 0L, identical(n$seed, o$seed), identical(n[, pc], o[, pc]))
    overlaps <- overlaps + 1L
  }
}
stopifnot(overlaps == 2L)
cat("OK both existing S2 replications have identical seeds and p-values\n")
missing <- !file.exists(dest)
stopifnot(sum(missing) == 4L, all(file.copy(fs[missing], dest[missing], overwrite = FALSE)))
d <- read_exp("libpilot3", old)
g <- read.csv(file.path(old, "driver_grid.csv"))
stopifnot(nrow(d) == 81L, all(is.na(d$err) | d$err == ""), all(d$B == 19), length(unique(d$cell)) == 27L)
for (id in g$cell) {
  z <- d[d$cell == id, ]
  stopifnot(nrow(z) == 3L, identical(sort(z$rep), 1:3),
            all(z$seed == g$seed0[match(id, g$cell)] + z$rep),
            all(z$lib == g$lib[match(id, g$cell)]))
}
s <- summarize_exp("libpilot3", old)
writeLines(character(), file.path(old, "pending_cells.txt"))
writeLines(c("TIMING TRIAL COMPLETE: 27 cells, 81 successful replications, B=19.",
             "Four missing raw chunks copied; original 77 chunks retained unchanged.",
             "Recovery: R 4.3.3, CPU affinity index 16, BLAS/OMP threads=1, R_ENABLE_JIT=0.",
             "Four same-seed evaluations in two fresh processes had identical p-values.",
             "Two overlap replications matched existing seeds and p-values exactly.",
             "Formal pilot has not been started; runtime stability is not fully established."),
           file.path(old, "environment_resume_status.txt"))
cat("OK all 81 rows: seeds, cell libraries, repetitions and error status verified\n")
cat(sprintf("COMPLETE: %d cells, %d successful repetitions, zero recorded errors\n", nrow(s), sum(s$reps)))
