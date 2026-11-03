# =============================================================================
# Tables for sim/试点报告_学习器库.md from results/<exp>/
#   Rscript report_libpilot.R [--exp=libpilot] [--out=results]   -> <out>/<exp>_tables.md
# The library sets are taken from the mca_set_* columns of the results.
# Reads the per-replication csv files only; writes nothing else.
# =============================================================================

script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(script_dir)
source(file.path("R", "summarize.R"))

args <- commandArgs(trailingOnly = TRUE)
exp <- sub("^--exp=", "", c(grep("^--exp=", args, value = TRUE), "--exp=libpilot")[1L])
alpha <- 0.05
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results")[1L])
d <- read_exp(exp, outdir)
d$err[is.na(d$err)] <- ""
ok <- d[d$err == "", ]

LRN  <- c("l1", "quad", "knn", "rf", "nn")
LRN  <- LRN[paste0("mca_", LRN) %in% names(d)]     # learners present in the results
SETS <- grep("^mca_set_", names(d), value = TRUE)
SETSIZE <- lengths(strsplit(SETS, "+", fixed = TRUE))
SETS <- SETS[order(SETSIZE)]; SETSIZE <- sort(SETSIZE)         # small libraries first
SETLAB <- paste0("{", gsub("+", ",", sub("^mca_set_", "", SETS), fixed = TRUE), "}")
COMP <- c("disco", "mmd", "mst", "schott")

lab <- function(x) {
  s <- sprintf("%s, d = %d, β = %.1f", x$model[1L], x$d[1L], x$beta[1L])
  sub("^C,", "C′,", s)
}
cells <- unique(d$cell)
# keep the order of the grid file
g <- read.csv(file.path(outdir, paste0(exp, "_grid.csv")), stringsAsFactors = FALSE)
cells <- g$cell[g$cell %in% cells]

out <- character(0)
md <- function(...) out <<- c(out, paste0(...))
fmt <- function(v, dig = 2) formatC(v, format = "f", digits = dig)
table_md <- function(df) {
  md("| ", paste(names(df), collapse = " | "), " |")
  md("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  for (r in seq_len(nrow(df))) md("| ", paste(unlist(df[r, ]), collapse = " | "), " |")
  md("")
}
rej <- function(x, p) mean(x[[p]] <= alpha, na.rm = TRUE)

md("# ", exp, " tables (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")\n")
md("Replications: ", nrow(d), "; errors: ", sum(d$err != ""), "; same-seed retries: ",
   sum(d$retries, na.rm = TRUE), "\n")

# ---- power / size table: every cell
md("## Rejection rates at α = 0.05 (all ", length(cells), " cells)\n")
rows <- lapply(cells, function(cl) {
  x <- ok[ok$cell == cl, ]
  v <- c(sapply(paste0("mca_", LRN), function(p) rej(x, p)),
         sapply(SETS, function(p) rej(x, p)), sapply(COMP, function(p) rej(x, p)))
  best <- max(v[seq_along(LRN)])
  c(cell = lab(x), reps = nrow(x), fmt(v), best_single = fmt(best),
    fmt(best - v[length(LRN) + seq_along(SETS)]))
})
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(tab) <- c("cell", "reps", LRN, SETLAB, "DISCO", "MMD", "MST", "Schott",
                "best single", paste("gap", SETLAB))
table_md(tab)

# ---- size table with 95% Monte Carlo range
nulls <- cells[sapply(cells, function(cl) ok$beta[ok$cell == cl][1L] == 0)]
if (length(nulls)) {
md("## Size (null cells) with the 95% Monte Carlo range\n")
rows <- lapply(nulls, function(cl) {
  x <- ok[ok$cell == cl, ]; n <- nrow(x)
  lo <- qbinom(0.025, n, alpha) / n; hi <- qbinom(0.975, n, alpha) / n
  v <- c(sapply(paste0("mca_", LRN), function(p) rej(x, p)), sapply(SETS, function(p) rej(x, p)),
         sapply(c("mca_2s", "mca_sum", "mca_comb"), function(p) rej(x, p)), sapply(COMP, function(p) rej(x, p)))
  flag <- ifelse(v < lo | v > hi, "!", "")
  c(cell = lab(x), reps = n, range = sprintf("[%.3f, %.3f]", lo, hi), paste0(fmt(v, 3), flag))
})
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(tab) <- c("cell", "reps", "95% MC range", LRN, SETLAB, "two-sided", "sum", "comb",
                "DISCO", "MMD", "MST", "Schott")
table_md(tab)
md("(! = outside the 95% Monte Carlo range)\n")
}

# ---- learned signal: mean test AUC of the deviant population
md("## Mean test AUC of the deviant population, per learner\n")
rows <- lapply(cells, function(cl) {
  x <- ok[ok$cell == cl, ]
  c(cell = lab(x), fmt(sapply(paste0("auc_", LRN, "_dev"), function(p) mean(x[[p]])), 3))
})
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE); names(tab) <- c("cell", LRN)
table_md(tab)

# ---- tuning values (each value is the mean over the k = 3 binary problems of a rep)
md("## Selected tuning values (per rep: mean over the k = 3 binary problems; median [min, max] over reps)\n")
q3 <- function(v, dig = 1) {
  if (is.null(v) || all(is.na(v))) return("—")                # learner not in this experiment
  sprintf("%s [%s, %s]", fmt(median(v, na.rm = TRUE), dig), fmt(min(v, na.rm = TRUE), dig),
          fmt(max(v, na.rm = TRUE), dig))
}
rows <- lapply(cells, function(cl) {
  x <- ok[ok$cell == cl, ]
  c(cell = lab(x), knn_K = q3(x$knn_K, 0), l1_log10C = q3(x$l1_log10C), quad_log10C = q3(x$quad_log10C),
    quad_nnz = q3(x$quad_nnz, 0), l1_nnz = q3(x$l1_nnz, 0), nn_epochs = q3(x$nn_epochs, 0))
})
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
table_md(tab)

# ---- timing
md("## Mean seconds per replication\n")
tc <- c(paste0("t_ovr_", LRN), "t_perm", "t_total")
rows <- lapply(cells, function(cl) {
  x <- ok[ok$cell == cl, ]
  c(cell = lab(x), fmt(colMeans(x[tc], na.rm = TRUE), 2))
})
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE); names(tab) <- c("cell", LRN, "perm", "total")
table_md(tab)
md("All cells: ", paste(paste(c(LRN, "total"), fmt(colMeans(ok[tc[-6]], na.rm = TRUE), 2)),
                         collapse = ", "), "\n")

# ---- library rule
md("## Library rule: max over alternative cells of (best single learner − library)\n")
alts <- setdiff(cells, nulls)
gap <- sapply(SETS, function(s) max(sapply(alts, function(cl) {
  x <- ok[ok$cell == cl, ]
  max(sapply(paste0("mca_", LRN), function(p) rej(x, p))) - rej(x, s)
})))
worst <- sapply(SETS, function(s) {
  gg <- sapply(alts, function(cl) {
    x <- ok[ok$cell == cl, ]
    max(sapply(paste0("mca_", LRN), function(p) rej(x, p))) - rej(x, s)
  })
  lab(ok[ok$cell == alts[which.max(gg)], ])
})
mean_pow <- sapply(SETS, function(s) mean(sapply(alts, function(cl) rej(ok[ok$cell == cl, ], s))))
table_md(data.frame(library = SETLAB, size = SETSIZE, max_gap = fmt(gap), worst_cell = unname(worst),
                    mean_power = fmt(mean_pow), meets_0.05 = unname(gap <= 0.05)))

writeLines(out, file.path(outdir, paste0(exp, "_tables.md")), useBytes = TRUE)
cat(out, sep = "\n")
