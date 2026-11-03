# Same formal B=200 S2 replication, repeated twice, for concurrent CPU probes.
args <- commandArgs(TRUE); stopifnot(length(args) == 1L)
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for (p in c("glmnet", "ranger", "gTestsMulti")) stopifnot(requireNamespace(p, quietly = TRUE))
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "garch_ext.R",
            "competitors2.R", "rotation.R", "learners_lib2.R", "models_s.R", "experiments.R", "summarize.R"))
  source(file.path(SRC_DIR, f))
g <- make_grid("libpilot3")
cell <- as.list(g[g$model == "S2" & abs(g$beta - 0.07) < 1e-9, ])
stopifnot(cell$B == 200L, length(cell$cell) == 1L)
rows <- vector("list", 2L)
for (i in 1:2) {
  cat(sprintf("[%s] formal B=200 S2 probe %d/2\n", format(Sys.time(), "%H:%M:%S"), i)); flush.console()
  rows[[i]] <- one_rep(cell, 1L)
  cat("completed in", rows[[i]]$t_total, "seconds\n"); flush.console()
}
pc <- union(intersect(P_COLS, names(rows[[1]])), grep("^(mca|macc)_", names(rows[[1]]), value = TRUE))
stopifnot(identical(rows[[1]][pc], rows[[2]][pc]))
saveRDS(list(p=rows[[1]][pc], seed=rows[[1]]$seed, seconds=vapply(rows, function(x) x$t_total, 0)), args[1L])
cat("PARALLEL CPU PROBE OK: identical repeated p-values\n")
