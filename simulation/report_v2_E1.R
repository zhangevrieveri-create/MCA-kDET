# Quick summary of E1 (revised design, size): Rscript report_v2_E1.R [--out=results_v2_quick]
# Rejection rates at alpha = 0.05 (and 0.01, 0.10 for the main statistic);
# "!" marks a rate outside the 95% Monte Carlo range of its cell.
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_quick")[1L])
source("R/summarize.R")
d <- read_exp("E1", outdir); d$err[is.na(d$err)] <- ""
cat("rows:", nrow(d), "| errors:", sum(d$err != ""), "| retries:", sum(d$retries, na.rm = TRUE), "\n")
ok <- d[d$err == "", ]
cols <- c("mca","mca_l1","mca_quad","mca_knn","mca_rf","mca_2s","mca_sum","mca_comb","pair",
          "shared_nn","shared_lib","macc","disco","mmd","mst","schott")
lab <- c("MCA","l1","quad","knn","rf","2s","sum","comb","pair","shared","sh_lib","macc",
         "DISCO","MMD","MST","Schott")
g <- read.csv(file.path(outdir, "E1_grid.csv"), stringsAsFactors = FALSE)
cells <- g$cell[g$cell %in% ok$cell]
flags <- 0L; tests <- 0L
rows <- lapply(cells, function(cl) {
  x <- ok[ok$cell == cl, ]; n <- nrow(x)
  lo <- qbinom(0.025, n, 0.05) / n; hi <- qbinom(0.975, n, 0.05) / n
  v <- sapply(cols, function(p) if (all(is.na(x[[p]]))) NA else mean(x[[p]] <= 0.05, na.rm = TRUE))
  out <- ifelse(is.na(v), "-", sprintf("%.3f%s", v, ifelse(v < lo | v > hi, "!", "")))
  flags <<- flags + sum(v < lo | v > hi, na.rm = TRUE); tests <<- tests + sum(!is.na(v))
  c(model = sub("^C$", "C(N(0,I))", x$model[1]), k = x$k[1], n = x$n[1], d = x$d[1], reps = n,
    range = sprintf("[%.3f,%.3f]", lo, hi), out,
    mca_a01 = sprintf("%.3f", mean(x$mca <= 0.01)), mca_a10 = sprintf("%.3f", mean(x$mca <= 0.10)))
})
t <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(t) <- c("model","k","n","d","reps","95% MC range",lab,"MCA a=.01","MCA a=.10")
print(t, row.names = FALSE)
cat(sprintf("\noutside the 95%% MC range: %d of %d (cell, method) pairs; about %.1f expected by chance\n",
            flags, tests, 0.05 * tests))
cat("\nmean seconds per rep by k:", round(tapply(ok$t_total, ok$k, mean), 1), "\n")
