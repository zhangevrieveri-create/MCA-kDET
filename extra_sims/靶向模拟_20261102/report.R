# Report of the targeted simulation: reads the saved replications only (no recomputation).
# Run from sim/:  Rscript 靶向模拟_20261102/report.R --out=<output folder>
args <- commandArgs(TRUE)
out <- sub("^--out=", "", grep("^--out=", args, value = TRUE))
stopifnot(length(out) == 1L, dir.exists(out))
grid <- readRDS(file.path(out, "grid.rds"))
files <- list.files(file.path(out, "rows"), pattern = "[.]rds$", full.names = TRUE)
rows <- lapply(files, readRDS)
tabdir <- file.path(out, "report_tables"); dir.create(tabdir, showWarnings = FALSE)

# integrity
key <- vapply(rows, function(r) paste(r$cell, r$rep), "")
seeds <- vapply(rows, function(r) r$seed, 1)
integ <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; rr <- Filter(function(r) r$cell == g$cell, rows)
  data.frame(cell = g$cell, expected = g$reps, saved = length(rr),
             unique_reps = length(unique(vapply(rr, `[[`, 1L, "rep"))),
             seeds_ok = all(vapply(rr, function(r) r$seed == g$seed0 + r$rep && r$B == g$B && r$err == "", TRUE)))
}))
write.csv(integ, file.path(tabdir, "integrity.csv"), row.names = FALSE)
stopifnot(!anyDuplicated(key), !anyDuplicated(seeds), all(integ$saved == integ$expected),
          all(integ$unique_reps == integ$expected), all(integ$seeds_ok))

# p-values, long
pv <- do.call(rbind, lapply(rows, function(r) data.frame(cell = r$cell, rep = r$rep, seed = r$seed,
  method = names(r$p), p = unlist(r$p), stringsAsFactors = FALSE)))
write.csv(pv, file.path(tabdir, "pvalues.csv"), row.names = FALSE)
rates <- do.call(rbind, lapply(split(pv, list(pv$cell, pv$method), drop = TRUE), function(z) {
  x <- sum(z$p <= 0.05); n <- nrow(z); ci <- binom.test(x, n)$conf.int
  data.frame(cell = z$cell[1], method = z$method[1], rejections = x, reps = n, rate = x / n,
             lower95 = ci[1], upper95 = ci[2])
}))
rates <- merge(grid[, c("cell", "part", "model", "target", "k", "n", "d")], rates, by = "cell")
write.csv(rates, file.path(tabdir, "all_method_rates_and_intervals.csv"), row.names = FALSE)

# learner AUCs of the deviating population (population k) and of all populations
au <- do.call(rbind, lapply(rows, function(r) cbind(cell = r$cell, rep = r$rep, k = r$k, r$aucs)))
aus <- aggregate(auc ~ cell + learner + population, data = au, FUN = mean)
write.csv(aus, file.path(tabdir, "learner_auc_summary.csv"), row.names = FALSE)
scr <- do.call(rbind, lapply(rows, function(r) cbind(cell = r$cell, rep = r$rep, k = r$k, r$deprf_screen)))
write.csv(scr, file.path(tabdir, "deprf_screen.csv"), row.names = FALSE)
tim <- do.call(rbind, lapply(rows, function(r) data.frame(cell = r$cell, rep = r$rep, seconds = r$seconds)))
write.csv(aggregate(seconds ~ cell, data = tim, FUN = mean), file.path(tabdir, "timing_summary.csv"), row.names = FALSE)

# main wide table
main <- c(L0 = "mca_set_l1+quad+knn+rf", L0P = "mca_set_l1+quad+knn+rf+deprf", deprf = "mca_deprf",
          rf = "mca_rf", knn = "mca_knn", L6 = "mca_set_xscr1+xscr3+knn+rf+g5f",
          L6P = "mca_set_xscr1+xscr3+knn+rf+g5f+deprf",
          SC_S = "songchen", SC_SA = "songchen_A", MST = "mst", KMD = "kmd", MMDmax = "mmd_max",
          DISCO = "disco", MMD = "mmd", CLXcov = "clx_cov", maxmean = "maxmean", ADmax = "adk_max")
wide <- grid[, c("cell", "part", "model", "target", "k", "n", "d", "reps")]
for (nm in names(main)) {
  z <- rates[rates$method == main[[nm]], c("cell", "rate")]; names(z)[2] <- nm
  wide <- merge(wide, z, by = "cell", all.x = TRUE)
}
dev <- aus[aus$population == vapply(aus$cell, function(c) grid$k[grid$cell == c], 1), ]
for (a in c("rf", "deprf", "knn")) {
  z <- dev[dev$learner == a, c("cell", "auc")]; names(z)[2] <- paste0("auc_", a)
  wide <- merge(wide, z, by = "cell", all.x = TRUE)
}
hk <- aggregate(hits ~ cell, data = scr[scr$population == scr$k, ], FUN = mean); names(hk)[2] <- "true_pairs_found"
wide <- merge(wide, hk, by = "cell", all.x = TRUE)
wide <- wide[order(wide$part, wide$model, wide$target == 0.5, wide$target, wide$k, wide$d), ]
write.csv(wide, file.path(tabdir, "main_table.csv"), row.names = FALSE)

# null diagnostics
nul <- rates[rates$model == "SM3" & rates$target == 0.5, ]
nul$ci_excludes_005 <- nul$lower95 > 0.05 | nul$upper95 < 0.05
write.csv(nul, file.path(tabdir, "null_rates_and_flags.csv"), row.names = FALSE)

fmt <- function(df) {
  df[] <- lapply(df, function(x) if (is.numeric(x)) format(round(x, 3), nsmall = 0) else x)
  c(paste("|", paste(names(df), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
lines <- c("# 靶向模拟 T1/T2 自动汇总", "", paste("生成时间：", format(Sys.time())), "",
  sprintf("%d 格，%d 次重复，完整性核对通过（次数、种子、B、无错误）。", nrow(grid), nrow(integ) * 0 + sum(integ$saved)), "",
  "## 主表（拒绝率，α = 0.05；auc_* 为偏离总体的平均留出 AUC；true_pairs_found 为 deprf 筛选命中的真依赖对数，满分 3）", "",
  fmt(wide[, setdiff(names(wide), "cell")]), "",
  "## 原假设：95% 区间不含 0.05 的方法", "",
  if (any(nul$ci_excludes_005)) fmt(nul[nul$ci_excludes_005, c("k", "d", "method", "rejections", "reps", "rate", "lower95", "upper95")]) else "无。", "")
writeLines(lines, file.path(out, "自动汇总.md"), useBytes = FALSE)
cat("REPORT DONE\n")
