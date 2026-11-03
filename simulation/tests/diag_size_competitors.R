# =============================================================================
# 诊断：在一个原假设格子里，用全新种子只跑对比方法（DISCO、MMD、MST），
# 检查它们的拒绝率是否落在精确置换检验应有的范围内。
# 只读调用冻结代码（R/models.R 的 gen_data，R/tests.R 的 competitor_tests），不修改任何文件。
#
# 用法（在 sim\ 目录下，先小样本试跑）：
#   "D:/R-4.6.1/bin/Rscript.exe" tests/diag_size_competitors.R --model=A --k=5 --d=100 --reps=20 --cores=4
#   "D:/R-4.6.1/bin/Rscript.exe" tests/diag_size_competitors.R --model=A --k=5 --d=100 --reps=2000 --cores=30
# 结果：results_diag\size_<model>_k<k>_d<d>_seed<seed>.csv（逐次 p 值）和同名 _summary.txt
# =============================================================================
args <- commandArgs(TRUE)
opt <- list(model = "A", k = "5", d = "100", n = "300", reps = "20", B = "200",
            cores = "4", seed = "987650000")
for (a in args) { kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]; opt[[kv[1]]] <- kv[2] }
model <- opt$model; k <- as.integer(opt$k); d <- as.integer(opt$d); n <- as.integer(opt$n)
reps <- as.integer(opt$reps); B <- as.integer(opt$B); cores <- as.integer(opt$cores)
seed0 <- as.numeric(opt$seed)
par <- list(cA = 0.6, cB = 4, sC = 10, cE = 0.2)   # beta = 0: these constants are not used

one <- function(r) {
  set.seed(seed0 + r)
  data <- gen_data(model, 0, rep(n, k), d, par)
  p <- competitor_tests(data, B, "dist")$p
  c(rep = r, seed = seed0 + r, disco = p$disco, mmd = p$mmd, mst = p$mst)
}

library(parallel)
cl <- makeCluster(cores)
clusterExport(cl, c("model", "k", "d", "n", "B", "seed0", "par"))
invisible(clusterEvalQ(cl, { source("R/models.R"); source("R/tests.R"); NULL }))
t0 <- Sys.time()
res <- parLapplyLB(cl, seq_len(reps), one)
stopCluster(cl)
D <- as.data.frame(do.call(rbind, res))

dir.create("results_diag", showWarnings = FALSE)
tag <- sprintf("size_%s_k%d_d%d_seed%.0f", model, k, d, seed0)
write.csv(D, file.path("results_diag", paste0(tag, ".csv")), row.names = FALSE)

a <- floor(0.05 * (B + 1)) / (B + 1)       # exact size of the permutation test
lines <- c(sprintf("model %s, k = %d, d = %d, n = %d, reps = %d, B = %d, seeds %.0f+1..%d, %.1f min",
                   model, k, d, n, reps, B, seed0, reps, as.numeric(difftime(Sys.time(), t0, units = "mins"))),
           sprintf("exact size of each test: %.4f; 95%% Monte Carlo range for %d reps: [%.3f, %.3f]",
                   a, reps, qbinom(0.025, reps, a) / reps, qbinom(0.975, reps, a) / reps))
for (m in c("disco", "mmd", "mst")) {
  x <- sum(D[[m]] <= 0.05)
  lines <- c(lines, sprintf("%-6s rejections %4d / %d = %.3f   P(Bin >= x) = %.2g   P(Bin <= x) = %.2g   mean p = %.3f",
                            m, x, reps, x / reps, pbinom(x - 1, reps, a, lower.tail = FALSE),
                            pbinom(x, reps, a), mean(D[[m]])))
}
lines <- c(lines, sprintf("unique seeds: %d of %d", length(unique(D$seed)), reps))
writeLines(lines, file.path("results_diag", paste0(tag, "_summary.txt")))
cat(lines, sep = "\n")
