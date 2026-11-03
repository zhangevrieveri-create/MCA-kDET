# Quick summary of E7 (revised design, timing): Rscript report_v2_E7.R [--out=results_v2_quick]
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_quick")[1L])
source("R/summarize.R")
d <- read_exp("E7", outdir); d$err[is.na(d$err)] <- ""
cat("rows:", nrow(d), "| errors:", sum(d$err != ""), "| retries:", sum(d$retries, na.rm = TRUE), "\n")
ok <- d[d$err == "", ]
grp <- list("MCA one-vs-rest (4 learners)" = c("t_ovr_l1","t_ovr_quad","t_ovr_knn","t_ovr_rf"),
            "  of which l1" = "t_ovr_l1", "  of which quad" = "t_ovr_quad", "  of which knn" = "t_ovr_knn",
            "  of which rf" = "t_ovr_rf", "MCA permutations (all classifier stats)" = "t_perm",
            "pairwise (4 learners)" = c("t_pair_l1","t_pair_quad","t_pair_knn","t_pair_rf"),
            "shared network" = "t_shared_nn", "distance matrix" = "t_dist_prep",
            "DISCO" = "t_disco", "MMD" = "t_mmd", "MST" = "t_mst", "Schott" = "t_schott", "total" = "t_total")
tab <- sapply(sort(unique(ok$k)), function(k) { x <- ok[ok$k == k, ]
  sapply(grp, function(cc) mean(rowSums(x[intersect(cc, names(x))])))})
colnames(tab) <- paste0("k=", sort(unique(ok$k)))
cat("\nmean seconds per replication (5 reps per k, model A, d = 100, n = 300 per population, B = 200):\n")
print(round(tab, 2))
