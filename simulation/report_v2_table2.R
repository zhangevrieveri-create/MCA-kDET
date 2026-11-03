# Table 2 of Section 4.3: power of the variants at selected cells of E2 (k = 3).
#   Rscript report_v2_table2.R [--out=results_v2_full]
# Writes <out>/E2_table2_variants.txt as well.
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_full")[1L])
source("R/summarize.R")
d <- read_exp("E2", outdir); d$err[is.na(d$err)] <- ""
ok <- d[d$err == "", ]
cells <- data.frame(label = c("(M1), beta = 0.8", "(M2), beta = 0.6", "(M3), beta = 0.8",
                              "(M4), beta = 0.6", "(M5), s = 3"),
                    model = c("A", "B", "C", "C2", "D"),   # model C of E2 is C' (sC = 10)
                    beta = c(0.8, 0.6, 0.8, 0.6, 0.03))
cols <- c(MCA = "mca", Sum = "mca_sum", Comb = "mca_comb", `Two-sided` = "mca_2s",
          Pair = "pair", Shared = "shared_nn")
tab <- do.call(rbind, lapply(seq_len(nrow(cells)), function(j) {
  x <- ok[ok$model == cells$model[j] & abs(ok$beta - cells$beta[j]) < 1e-8, ]
  data.frame(cell = cells$label[j], reps = nrow(x),
             t(sprintf("%.2f", vapply(cols, function(p) mean(x[[p]] <= 0.05, na.rm = TRUE), numeric(1)))),
             check.names = FALSE)
}))
names(tab)[3:(2 + length(cols))] <- names(cols)
sink(file.path(outdir, "E2_table2_variants.txt"), split = TRUE)
cat("Table 2: power at alpha = 0.05 of the variants, E2, k = 3 |", outdir, "| errors:", sum(d$err != ""), "\n\n")
print(tab, row.names = FALSE)
sink()
