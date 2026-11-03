# Summary of E8 (size under a GARCH null): Rscript report_v2_E8.R [--out=results_v2_full]
# Writes the table to <out>/E8_summary_table.txt as well.
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_full")[1L])
source("R/summarize.R")
d <- read_exp("E8", outdir); d$err[is.na(d$err)] <- ""
sink(file.path(outdir, "E8_summary_table.txt"), split = TRUE)
cat("E8, GARCH null |", outdir, "| rows:", nrow(d), "| errors:", sum(d$err != ""),
    "| retries:", sum(d$retries, na.rm = TRUE), "\n\n")
ok <- d[d$err == "", ]
rr <- function(x, p, a = 0.05) mean(x[[p]] <= a, na.rm = TRUE)

cat("Rejection rates at alpha = 0.05 (flt = GARCH-filtered, raw = unfiltered;\n",
    "blk = block split with buffer, iid = random within-regime split)\n", sep = "")
cols <- c(flt_blk = "mca", flt_iid = "flt_iid_mca", raw_blk = "raw_blk_mca", raw_iid = "raw_iid_mca",
          flt_blk_sum = "mca_sum", DISCO_flt = "disco", MMD_flt = "mmd", MST_flt = "mst",
          DISCO_raw = "raw_disco", MMD_raw = "raw_mmd", MST_raw = "raw_mst")
tab <- do.call(rbind, lapply(split(ok, ok$cell), function(x) {
  R <- nrow(x)
  data.frame(cell = sub("^E8_GARCH_b0.00_k3_n300_", "", x$cell[1L]), reps = R,
             t(sprintf("%.3f", vapply(cols, function(p) rr(x, p), numeric(1)))),
             range95 = sprintf("[%.3f, %.3f]", qbinom(.025, R, .05) / R, qbinom(.975, R, .05) / R),
             check.names = FALSE)
}))
names(tab)[3:(2 + length(cols))] <- names(cols)
print(tab, row.names = FALSE)

cat("\nMCA-kDET of the Section 5 procedure (flt_blk): alpha = 0.01 / 0.05 / 0.10, KS test of uniform p-values:\n")
ks <- do.call(rbind, lapply(split(ok, ok$cell), function(x)
  data.frame(cell = sub("^E8_GARCH_b0.00_k3_n300_", "", x$cell[1L]),
             a01 = sprintf("%.3f", rr(x, "mca", .01)), a05 = sprintf("%.3f", rr(x, "mca", .05)),
             a10 = sprintf("%.3f", rr(x, "mca", .10)),
             KS_p = sprintf("%.3f", suppressWarnings(ks.test(x$mca, "punif")$p.value)))))
print(ks, row.names = FALSE)

cat("\nSingle learners and variants at alpha = 0.05, all four versions:\n")
vers <- c(flt_blk = "", flt_iid = "flt_iid_", raw_blk = "raw_blk_", raw_iid = "raw_iid_")
sl <- do.call(rbind, lapply(split(ok, ok$cell), function(x) do.call(rbind, lapply(names(vers), function(v)
  data.frame(cell = sub("^E8_GARCH_b0.00_k3_n300_", "", x$cell[1L]), version = v,
             t(sprintf("%.3f", vapply(paste0(vers[[v]], c("mca", "mca_l1", "mca_quad", "mca_knn", "mca_rf",
                                                          "mca_2s", "mca_sum", "mca_comb")),
                                      function(p) rr(x, p), numeric(1)))))))))
names(sl)[-(1:2)] <- c("MCA", "l1", "quad", "knn", "rf", "2s", "sum", "comb")
print(sl, row.names = FALSE)

cat("\nDiagnostics (means): lag-1 autocorrelation of squares raw -> filtered, estimated persistence,\n",
    "weeks fit / test / dropped in the block split, filter non-convergence, seconds per replication\n", sep = "")
dg <- aggregate(ok[c("sqacf_raw", "sqacf_flt", "persist_hat", "n_fit_blk", "n_test_blk", "n_drop_blk",
                     "filter_nonconv", "t_total")], list(cell = ok$cell), mean)
dg$cell <- sub("^E8_GARCH_b0.00_k3_n300_", "", dg$cell)
print(cbind(dg[1], round(dg[-1], 3)), row.names = FALSE)
sink()
