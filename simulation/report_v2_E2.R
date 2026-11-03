# Quick summary of E2 (revised design): Rscript report_v2_E2.R  (reads results_v2_quick/E2)
source("R/summarize.R")
d <- read_exp("E2", "results_v2_quick"); d$err[is.na(d$err)] <- ""
cat("rows:", nrow(d), "| errors:", sum(d$err != ""), "| retries:", sum(d$retries, na.rm = TRUE), "\n")
ok <- d[d$err == "", ]
cols <- c("mca_l1","mca_quad","mca_knn","mca_rf","mca","mca_sum","mca_comb","pair","shared_nn","macc","disco","mmd","mst","schott")
lab <- c("l1","quad","knn","rf","MCA","sum","comb","pair","shared","macc","DISCO","MMD","MST","Schott")
g <- read.csv("results_v2_quick/E2_grid.csv"); cells <- g$cell[g$cell %in% ok$cell]
rows <- lapply(cells, function(cl) { x <- ok[ok$cell == cl, ]
  c(model = sub("^C$", "C'", x$model[1]), beta = sprintf("%.2f", x$beta[1]), reps = nrow(x),
    sprintf("%.2f", sapply(cols, function(p) mean(x[[p]] <= 0.05, na.rm = TRUE)))) })
t <- as.data.frame(do.call(rbind, rows)); names(t) <- c("model","beta","reps",lab)
print(t, row.names = FALSE)
auc <- aggregate(ok[c("auc_l1_dev","auc_quad_dev","auc_knn_dev","auc_rf_dev")], list(model = ok$model, beta = ok$beta), mean)
cat("\nmean deviant AUC:\n"); print(cbind(auc[1:2], round(auc[-(1:2)], 3)), row.names = FALSE)
tc <- c("t_ovr_l1","t_ovr_quad","t_ovr_knn","t_ovr_rf","t_pair_l1","t_pair_quad","t_pair_knn","t_pair_rf","t_shared_nn","t_total")
tt <- aggregate(ok[intersect(tc, names(ok))], list(model = ok$model), mean)
cat("\nmean seconds per rep:\n"); print(cbind(tt[1], round(tt[-1], 2)), row.names = FALSE)
