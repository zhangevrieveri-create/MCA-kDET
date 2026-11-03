# Quick summary of E6 (revised design, dimension): Rscript report_v2_E6.R [--out=results_v2_quick]
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_quick")[1L])
source("R/summarize.R")
d <- read_exp("E6", outdir); d$err[is.na(d$err)] <- ""
cat("rows:", nrow(d), "| errors:", sum(d$err != ""), "| retries:", sum(d$retries, na.rm = TRUE), "\n")
ok <- d[d$err == "", ]
cols <- c("mca_l1","mca_quad","mca_knn","mca_rf","mca","mca_sum","mca_comb","shared_nn","macc","disco","mmd","mst","schott")
lab <- c("l1","quad","knn","rf","MCA","sum","comb","shared","macc","DISCO","MMD","MST","Schott")
key <- unique(ok[order(ok$model, ok$d), c("model", "d")])
rows <- lapply(seq_len(nrow(key)), function(i) { x <- ok[ok$model == key$model[i] & ok$d == key$d[i], ]
  c(model = sub("^C$", "C'", key$model[i]), d = key$d[i], beta = x$beta[1], sC = x$sC[1], reps = nrow(x),
    sprintf("%.2f", sapply(cols, function(p) mean(x[[p]] <= 0.05, na.rm = TRUE)))) })
t <- as.data.frame(do.call(rbind, rows)); names(t) <- c("model","d","beta","sC","reps",lab); print(t, row.names = FALSE)
auc <- aggregate(ok[c("auc_l1_dev","auc_quad_dev","auc_knn_dev","auc_rf_dev")], list(model = ok$model, d = ok$d), mean)
cat("\nmean deviant AUC:\n"); print(cbind(auc[1:2], round(auc[-(1:2)], 3)), row.names = FALSE)
tt <- aggregate(ok[c("t_ovr_l1","t_ovr_quad","t_ovr_knn","t_ovr_rf","t_total")], list(model = ok$model, d = ok$d), mean)
cat("\nmean seconds per rep:\n"); print(cbind(tt[1:2], round(tt[-(1:2)], 2)), row.names = FALSE)
