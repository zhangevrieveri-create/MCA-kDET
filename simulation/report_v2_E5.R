# Quick summary of E5 (revised design, split fraction and multiple splits): Rscript report_v2_E5.R [--out=results_v2_quick]
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_quick")[1L])
source("R/summarize.R")
d <- read_exp("E5", outdir); d$err[is.na(d$err)] <- ""
cat("rows:", nrow(d), "| errors:", sum(d$err != ""), "| retries:", sum(d$retries, na.rm = TRUE), "\n")
ok <- d[d$err == "", ]
cols <- c("mca","mca_l1","mca_quad","mca_knn","mca_rf","mca_sum","mca_comb","mca_ms_med","mca_ms_mean")
lab <- c("MCA(1 split)","l1","quad","knn","rf","sum","comb","2x median(5)","2x mean(5)")
key <- unique(ok[order(ok$model, ok$beta, ok$tau), c("model", "beta", "tau")])
rows <- lapply(seq_len(nrow(key)), function(i) {
  x <- ok[ok$model == key$model[i] & ok$beta == key$beta[i] & ok$tau == key$tau[i], ]
  c(model = sub("^C$", "C'", key$model[i]), beta = key$beta[i], tau = key$tau[i], reps = nrow(x),
    sprintf("%.3f", sapply(cols, function(p) mean(x[[p]] <= 0.05, na.rm = TRUE)))) })
t <- as.data.frame(do.call(rbind, rows)); names(t) <- c("model","beta","tau","reps",lab); print(t, row.names = FALSE)
cat(sprintf("\n95%% MC range for the null cells (200 reps): [%.3f, %.3f]\n", qbinom(.025, 200, .05) / 200, qbinom(.975, 200, .05) / 200))
auc <- aggregate(ok[c("auc_l1_dev","auc_quad_dev","auc_knn_dev","auc_rf_dev")], list(model = ok$model, beta = ok$beta, tau = ok$tau), mean)
cat("\nmean deviant AUC (first split):\n"); print(cbind(auc[1:3], round(auc[-(1:3)], 3)), row.names = FALSE)
cat("\nmean seconds per rep (5 splits):", round(tapply(ok$t_total, ok$model, mean), 1), "\n")
