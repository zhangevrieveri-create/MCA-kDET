# =============================================================================
# Step 4 of the glmnet-fix check: compare the replications rerun with the new
# code (maxit 1e4) with the same replications of the quick E1 run (maxit 1e5).
#   Rscript tests/compare_glmnet_fix.R
# Reads results_check_glmnet/E1 and results_v2_quick/E1; writes nothing.
# =============================================================================
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
source("R/summarize.R")
new <- read_exp("E1", "results_check_glmnet"); old <- read_exp("E1", "results_v2_quick")
old <- old[old$cell %in% new$cell, ]
key <- paste(new$cell, new$rep); okey <- paste(old$cell, old$rep)
old <- old[match(key, okey), ]
stopifnot(all(paste(old$cell, old$rep) == key))
pc <- intersect(grep("^(mca|pair|shared|macc|disco|mmd|mst|schott)", names(new), value = TRUE), names(old))
for (cl in unique(new$cell)) {
  ix <- new$cell == cl
  same <- sapply(pc, function(p) sum(new[[p]][ix] == old[[p]][ix], na.rm = TRUE))
  cat(sprintf("\n== %s: %d reps\n", cl, sum(ix)))
  cat("  p-value columns identical in all reps:", sum(same == sum(ix)), "of", length(pc), "\n")
  diffcols <- names(same)[same < sum(ix)]
  if (length(diffcols)) for (p in diffcols)
    cat(sprintf("  %-10s differs in %d reps; rejections old %d -> new %d\n", p, sum(ix) - same[[p]],
                sum(old[[p]][ix] <= .05), sum(new[[p]][ix] <= .05)))
  for (a in c("l1", "quad")) {
    nc <- new[[paste0(a, "_nonconv")]][ix]
    if (!is.null(nc)) cat(sprintf("  %s: mean non-converged fits per binary problem %.2f\n", a, mean(nc)))
    t0 <- old[[paste0("t_ovr_", a)]][ix]; t1 <- new[[paste0("t_ovr_", a)]][ix]
    cat(sprintf("  t_ovr_%s: old %.1f s -> new %.1f s per rep\n", a, mean(t0), mean(t1)))
  }
  cat(sprintf("  t_total: old %.1f s -> new %.1f s per rep\n", mean(old$t_total[ix]), mean(new$t_total[ix])))
}
