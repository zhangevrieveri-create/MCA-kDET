# =============================================================================
# Tables for the third review round (G1, G2), 1 Nov 2026
#
#   "D:/R-4.6.1/bin/Rscript.exe" report_review3.R --out=results_v2_quick
#
# Reads E2, E9 and E10 from the output directory and writes
#   <out>/review3_power.csv     E2 cells: MCA-kDET and the old competitors
#                               (from E2) next to the new competitors (E9),
#                               paired replication by replication
#   <out>/review3_size.csv      E9 null cells (k = 3, equal n)
#   <out>/review3_rotated.csv   E10: rotated (M3)-(M5), all methods, with the
#                               unrotated MCA-kDET of E2 for comparison
#   <out>/review3_accuracy.csv  G2: AUC vs accuracy from the same fits
#                               (E10, and E2 when it was run with the new code)
#   <out>/review3_tables.txt    the four tables as text
# Paper names: A = (M1), B = (M2), C = (M3), C2 = (M4), D = (M5).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1L]) else default
}
out <- get_arg("out", "results_v2_quick")
source("R/summarize.R")

MNAME <- c(A = "(M1)", B = "(M2)", C = "(M3)", C2 = "(M4)", D = "(M5)")
summ <- function(x) {
  s <- summarize_exp(x, out)
  if (!nrow(s)) stop("no results for ", x, " in ", out)
  s$key <- sub("^E[0-9]+_", "", s$cell)
  s$M <- MNAME[s$model]
  s
}
pick <- function(s, cols) {
  cols <- intersect(cols, names(s))
  s[, cols, drop = FALSE]
}
fmt <- function(df) {
  for (v in names(df)) if (is.numeric(df[[v]]) && !v %in% c("reps", "d", "k"))
    df[[v]] <- ifelse(is.na(df[[v]]), "", sprintf("%.2f", df[[v]]))
  df
}
OLD <- c("mca", "mca_l1", "mca_quad", "mca_knn", "mca_rf", "pair", "macc", "macc_ovr",
         "disco", "mmd", "mst", "schott")
NEW <- c("maxmean", "clx_mean", "clx_cov", "schott_perm", "mmd_max", "kmd",
         "songchen", "songchen_A")

e2 <- summ("E2"); e9 <- summ("E9")
lines <- character(0)
add <- function(title, df) {
  lines <<- c(lines, "", title, strrep("-", nchar(title)),
              utils::capture.output(print(fmt(df), row.names = FALSE)))
}

# ---- power: E2 cells, paired
p9 <- e9[e9$beta > 0, ]
pw <- merge(pick(e2, c("key", "M", "beta", "d", "reps", OLD)),
            pick(p9, c("key", "reps", NEW)), by = "key", suffixes = c("", "_E9"))
pw <- pw[order(pw$M, pw$beta), setdiff(names(pw), "key")]
write.csv(pw, file.path(out, "review3_power.csv"), row.names = FALSE)
add("Power at alpha = 0.05: MCA-kDET and old competitors (E2) | new competitors (E9), same data", pw)

# ---- size: E9 null cells
s9 <- e9[e9$beta == 0, ]
s9 <- s9[order(s9$model, s9$d), ]
s9$M <- c(A = "(M1)/(M2) null", C = "(M3)/(M4) null", D = "(M5) null")[s9$model]
s9 <- pick(s9, c("M", "d", "reps", NEW))
write.csv(s9, file.path(out, "review3_size.csv"), row.names = FALSE)
r <- if (nrow(s9)) s9$reps[1L] else NA
add(sprintf("Size of the new competitors (E9, k = 3, n = 300); 95%% MC range for %s reps: [%.3f, %.3f]",
            r, 0.05 - 1.96 * sqrt(0.0475 / r), 0.05 + 1.96 * sqrt(0.0475 / r)), s9)

# ---- rotated models
if (dir.exists(file.path(out, "E10"))) {
  e10 <- summ("E10")
  rt <- merge(pick(e10, c("key", "M", "beta", "d", "reps", OLD, NEW)),
              setNames(pick(e2, c("key", "mca")), c("key", "mca_unrotated")), by = "key")
  rt <- rt[order(rt$M, rt$beta), setdiff(names(rt), "key")]
  write.csv(rt, file.path(out, "review3_rotated.csv"), row.names = FALSE)
  add("Rotated (M3)-(M5) (E10): all methods; mca_unrotated = MCA-kDET on the same data before rotation (E2)", rt)

  # ---- G2: AUC vs accuracy from the same fits
  AC <- c("mca", "macc_ovr", "mca_l1", "macc_l1", "mca_quad", "macc_quad",
          "mca_knn", "macc_knn", "mca_rf", "macc_rf", "macc")
  ac <- rbind(if ("macc_ovr" %in% names(e2)) cbind(rot = "no", pick(e2, c("M", "beta", "reps", AC))),
              cbind(rot = "yes", pick(e10, c("M", "beta", "reps", AC))))
  ac <- ac[order(ac$rot, ac$M, ac$beta), ]
  write.csv(ac, file.path(out, "review3_accuracy.csv"), row.names = FALSE)
  add("G2: standardized AUC (mca_*) vs accuracy (macc_*) from the same one-vs-rest fits and permutations; macc = shared network", ac)
}

writeLines(lines, file.path(out, "review3_tables.txt"))
cat(lines, sep = "\n")
