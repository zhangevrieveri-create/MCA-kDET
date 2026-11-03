# Quick summary of E3 (revised design, effect of k): Rscript report_v2_E3.R [--out=results_v2_quick]
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_quick")[1L])
source("R/summarize.R")
d <- read_exp("E3", outdir); d$err[is.na(d$err)] <- ""
cat("rows:", nrow(d), "| errors:", sum(d$err != ""), "| retries:", sum(d$retries, na.rm = TRUE), "\n")
ok <- d[d$err == "", ]
cols <- c("mca","mca_l1","mca_quad","mca_knn","mca_rf","mca_sum","mca_comb","pair","pair_l1","pair_quad","pair_knn","pair_rf","shared_nn","shared_lib","macc")
lab <- c("MCA","l1","quad","knn","rf","sum","comb","PAIR","p_l1","p_quad","p_knn","p_rf","shared","sh_lib","macc")
key <- unique(ok[order(ok$model, ok$k), c("model", "k")])
rows <- lapply(seq_len(nrow(key)), function(i) { x <- ok[ok$model == key$model[i] & ok$k == key$k[i], ]
  c(model = sub("^C$", "C'", key$model[i]), k = key$k[i], reps = nrow(x),
    sprintf("%.2f", sapply(cols, function(p) mean(x[[p]] <= 0.05, na.rm = TRUE)))) })
t <- as.data.frame(do.call(rbind, rows)); names(t) <- c("model","k","reps",lab); print(t, row.names = FALSE)
auc <- aggregate(ok[c("auc_l1_dev","auc_quad_dev","auc_knn_dev","auc_rf_dev")], list(model = ok$model, k = ok$k), mean)
cat("\nmean deviant AUC (one-vs-rest):\n"); print(cbind(auc[1:2], round(auc[-(1:2)], 3)), row.names = FALSE)
tc <- intersect(c("t_ovr_l1","t_ovr_quad","t_ovr_knn","t_ovr_rf","t_pair_l1","t_pair_quad","t_pair_knn","t_pair_rf","t_shared_nn","t_perm","t_total"), names(ok))
tt <- aggregate(ok[tc], list(model = ok$model, k = ok$k), mean)
cat("\nmean seconds per rep:\n"); print(cbind(tt[1:2], round(tt[-(1:2)], 1)), row.names = FALSE)
