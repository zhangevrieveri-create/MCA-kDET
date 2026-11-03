# Quick summary of E4 (revised design, model E): Rscript report_v2_E4.R [--out=results_v2_quick]
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_quick")[1L])
source("R/summarize.R")
d <- read_exp("E4", outdir); d$err[is.na(d$err)] <- ""
cat("rows:", nrow(d), "| errors:", sum(d$err != ""), "| retries:", sum(d$retries, na.rm = TRUE), "\n")
ok <- d[d$err == "", ]
cols <- c("mca_l1","mca_quad","mca_knn","mca_rf","mca","mca_2s","mca_sum","mca_comb","shared_nn","macc","disco","mmd","mst")
lab <- c("l1","quad","knn","rf","MCA","2s","sum","comb","shared","macc","DISCO","MMD","MST")
rows <- lapply(sort(unique(ok$beta)), function(b) { x <- ok[ok$beta == b, ]
  c(beta = sprintf("%.1f", b), reps = nrow(x), sprintf("%.3f", sapply(cols, function(p) mean(x[[p]] <= 0.05, na.rm = TRUE)))) })
t <- as.data.frame(do.call(rbind, rows)); names(t) <- c("beta","reps",lab); print(t, row.names = FALSE)
n0 <- sum(ok$beta == 0); cat(sprintf("\n95%% MC range at beta = 0 (%d reps): [%.3f, %.3f]\n", n0, qbinom(.025, n0, .05) / n0, qbinom(.975, n0, .05) / n0))
auc <- aggregate(ok[c("auc_l1_dev","auc_quad_dev","auc_knn_dev","auc_rf_dev")], list(beta = ok$beta), mean)
cat("\nmean deviant AUC:\n"); print(cbind(auc[1], round(auc[-1], 3)), row.names = FALSE)
cat("\nmean seconds per rep:", round(mean(ok$t_total), 2), "\n")
