# =============================================================================
# 只读审计：E1（水平）结果的重复是否独立、各对比方法的拒绝率、p 值是否在 k/(B+1) 网格上。
# 不重算、不修改任何结果文件。
#
# 用法（在 sim\ 目录下）：
#   "D:/R-4.6.1/bin/Rscript.exe" tests/audit_E1_size.R results_v2_quick/E1 [另一份结果目录 ...]
# 如果表 1 的 1000 次结果在别的机器或别的目录，把那个目录也作为参数传进来。
# 报告写到 audit_E1_size_report.txt。
# =============================================================================
dirs <- commandArgs(TRUE)
if (!length(dirs)) dirs <- c("results_v2_full/E1", "results_v2_quick/E1")
pcols <- c("mca", "macc", "disco", "mmd", "mst", "schott")
out <- character()
say <- function(...) { s <- sprintf(...); cat(s, "\n"); out <<- c(out, s) }
seed0_by_dir <- list()

for (dd in dirs) {
  if (!dir.exists(dd)) { say("== %s: directory not found", dd); next }
  fs <- list.files(dd, pattern = "\\.csv$", full.names = TRUE)
  say("== %s: %d chunk files", dd, length(fs))
  if (!length(fs)) next
  dl <- lapply(fs, function(f) read.csv(f, stringsAsFactors = FALSE))
  allc <- Reduce(union, lapply(dl, names))
  dl <- lapply(dl, function(x) { for (m in setdiff(allc, names(x))) x[[m]] <- NA; x[allc] })
  D <- do.call(rbind, dl)
  cells <- sort(unique(D$cell))
  s0 <- setNames(numeric(length(cells)), cells)
  for (cl in cells) {
    x <- D[D$cell == cl, ]
    Bc <- suppressWarnings(as.integer(x$B[1]))
    nrow0 <- nrow(x); nseed <- length(unique(x$seed)); nrep <- length(unique(x$rep))
    s0[cl] <- min(x$seed - x$rep)
    flag <- if (nrow0 != nseed || nseed != nrep || length(unique(x$seed - x$rep)) != 1) "  <-- 重复或种子不一致" else ""
    say("%s: rows %d, unique seeds %d, unique reps %d, B = %s, seed0 %.0f%s",
        cl, nrow0, nseed, nrep, Bc, s0[cl], flag)
    xs <- x[!duplicated(x$seed), ]
    a <- if (is.na(Bc)) 0.05 else floor(0.05 * (Bc + 1)) / (Bc + 1)
    for (p in intersect(pcols, names(xs))) {
      v <- xs[[p]]; v <- v[!is.na(v)]
      if (!length(v)) next
      xr <- sum(v <= 0.05); nn <- length(v)
      grid <- if (is.na(Bc) || p == "schott") NA else mean(abs(v * (Bc + 1) - round(v * (Bc + 1))) < 1e-6)
      say("    %-6s rej %.3f (%d/%d)  P(Bin>=x) %.2g  P(Bin<=x) %.2g  share on k/(B+1) grid %s",
          p, xr / nn, xr, nn, pbinom(xr - 1, nn, a, lower.tail = FALSE), pbinom(xr, nn, a),
          if (is.na(grid)) "-" else sprintf("%.2f", grid))
    }
    # the same seed appearing twice with different p-values would mean the code changed between runs
    dseed <- x$seed[duplicated(x$seed)]
    if (length(dseed)) {
      y <- x[x$seed %in% dseed, ]
      ndiff <- sum(tapply(y$disco, y$seed, function(z) length(unique(round(z, 8))) > 1))
      say("    %d duplicated seeds; %d of them with different DISCO p-values", length(dseed), ndiff)
    }
  }
  seed0_by_dir[[dd]] <- s0
}

# overlap of seeds between directories (e.g. a "second run" that reuses most seeds of the first)
if (length(seed0_by_dir) >= 2) {
  nm <- names(seed0_by_dir)
  for (i in 1:(length(nm) - 1)) for (j in (i + 1):length(nm)) {
    a <- seed0_by_dir[[nm[i]]]; b <- seed0_by_dir[[nm[j]]]
    common <- intersect(names(a), names(b))
    same <- sum(a[common] == b[common])
    say("seed0 shared between %s and %s: %d of %d common cells", nm[i], nm[j], same, length(common))
  }
}
writeLines(out, "audit_E1_size_report.txt")
