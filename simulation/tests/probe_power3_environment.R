# Environment-only validation; the scientific algorithms are sourced unchanged.
# Run from sim/: Rscript --vanilla tests/probe_power3_environment.R <output.rds>
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
cat("R:", R.version.string, "\nLibraries:", paste(.libPaths(), collapse = "; "), "\n")
for (pkg in c("glmnet", "ranger", "gTestsMulti")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Cannot load package: ", pkg)
  cat(pkg, as.character(utils::packageVersion(pkg)), "\n")
}
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "garch_ext.R",
            "competitors2.R", "rotation.R", "learners_lib2.R", "models_s.R", "experiments.R", "summarize.R"))
  source(file.path(SRC_DIR, f))
grid <- read.csv("results_trial_power3_progress_fixed_20261101/driver_grid.csv", stringsAsFactors = FALSE)
cell <- as.list(grid[grid$model == "S2" & abs(grid$beta - 0.07) < 1e-9, ])
stopifnot(length(cell$cell) == 1L)
results <- vector("list", 2L)
for (i in 1:2) {
  cat(sprintf("[%s] [%-20s] %d/2: SAME S2 cell, SAME replication 1\n", format(Sys.time(), "%H:%M:%S"),
              strrep("#", 10L * (i - 1L)), i - 1L)); flush.console()
  results[[i]] <- one_rep(cell, 1L)
  cat("completed evaluation", i, "in", results[[i]]$t_total, "seconds\n"); flush.console()
}
pc <- union(intersect(P_COLS, names(results[[1]])), grep("^(mca|macc)_", names(results[[1]]), value = TRUE))
stopifnot(length(pc) > 0L, identical(results[[1]][pc], results[[2]][pc]))
payload <- list(R = R.version.string, libraries = .libPaths(), cell = cell$cell, seed = results[[1]]$seed,
                p = results[[1]][pc], rows = results,
                packages = vapply(c("glmnet", "ranger", "gTestsMulti"), function(p) as.character(utils::packageVersion(p)), ""))
saveRDS(payload, args[1L])
cat("[####################] 2/2 ENVIRONMENT CHECK OK; identical p-values\n")
