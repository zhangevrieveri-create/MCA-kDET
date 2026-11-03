# =============================================================================
# Figures 1 and 2 of Section 4 (base R, no extra packages).
#
#   "D:\R-4.6.1\bin\Rscript.exe" plot_v2.R --out=results_v2_full
#   "D:\R-4.6.1\bin\Rscript.exe" plot_v2.R --out=results_v2_quick   (preview)
#
# Reads <out>/E2, <out>/E3 and <out>/E6 (per-replication csv chunks) and writes
#   <out>/figures/fig1_power.pdf   power against beta in (M1)-(M5), k = 3
#   <out>/figures/fig2_auc.pdf     left: power against the largest learned AUC
#                                  right: power in (M1), beta = 0.8, against d
#   <out>/figures/fig_numbers.txt  the numbers quoted in Section 4.3
# Nothing else is written; results are not modified.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_full")[1L])
source("R/summarize.R")
figdir <- file.path(outdir, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

rd <- function(e) {
  d <- read_exp(e, outdir)
  if (is.null(d)) stop("no results for ", e, " in ", outdir)
  d$err[is.na(d$err)] <- ""
  d[d$err == "", ]
}
E2 <- rd("E2"); E3 <- rd("E3"); E6 <- rd("E6")
alpha <- 0.05
pw <- function(p) mean(p <= alpha, na.rm = TRUE)

# ---------------------------------------------------------------- style
COL <- c(mca = "black", mca_l1 = "#1f5fa8", mca_quad = "#d0702a", mca_knn = "#2e8b57",
         mca_rf = "#8a3fa0", disco = "grey35", mmd = "grey55", mst = "grey20",
         schott = "grey45", macc = "grey70")
LTY <- c(mca = 1, mca_l1 = 2, mca_quad = 2, mca_knn = 2, mca_rf = 2,
         disco = 1, mmd = 3, mst = 4, schott = 5, macc = 6)
LWD <- c(mca = 2.4, mca_l1 = 1.2, mca_quad = 1.2, mca_knn = 1.2, mca_rf = 1.2,
         disco = 1.2, mmd = 1.2, mst = 1.2, schott = 1.2, macc = 1.2)
PCH <- c(mca = 19, mca_l1 = NA, mca_quad = NA, mca_knn = NA, mca_rf = NA,
         disco = 1, mmd = 2, mst = 0, schott = 5, macc = 4)
LAB <- c(mca = "MCA-kDET", mca_l1 = "l1", mca_quad = "quad", mca_knn = "knn", mca_rf = "rf",
         disco = "DISCO", mmd = "MMD", mst = "MST", schott = "Schott", macc = "Accuracy")
METH <- names(LAB)

# ---------------------------------------------------------------- Figure 1
panels <- list(A = "(M1)", B = "(M2)", C = "(M3)", C2 = "(M4)", D = "(M5)")
num <- character()
pdf(file.path(figdir, "fig1_power.pdf"), width = 11, height = 3.4)
layout(matrix(c(1:5, rep(6, 5)), 2, byrow = TRUE), heights = c(1, 0.22))
par(mar = c(3.6, 3.2, 2, 0.6), mgp = c(2.1, 0.6, 0), las = 1, cex.axis = 0.85)
for (m in names(panels)) {
  x <- E2[E2$model == m & E2$beta > 0, ]
  b <- sort(unique(x$beta))
  xv <- if (m == "D") ceiling(100 * b) else b
  P <- sapply(METH, function(p) sapply(b, function(v) pw(x[x$beta == v, p])))
  if (length(b) == 1L) P <- matrix(P, nrow = 1L, dimnames = list(NULL, METH))
  plot(range(xv), c(0, 1), type = "n", xlab = if (m == "D") expression(italic(s)) else expression(beta),
       ylab = if (m == "A") "Power" else "", main = panels[[m]], xaxt = "n")
  axis(1, at = xv, labels = format(xv))
  abline(h = alpha, lty = 3, col = "grey60")
  for (p in rev(METH)) {
    lines(xv, P[, p], col = COL[p], lty = LTY[p], lwd = LWD[p])
    if (!is.na(PCH[p])) points(xv, P[, p], col = COL[p], pch = PCH[p], cex = 0.8)
  }
  num <- c(num, sprintf("%s: %s", panels[[m]],
           paste(sprintf("%s=%s", format(xv), apply(round(P, 2), 1, function(r) paste(r, collapse = "/"))), collapse = "; ")))
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = LAB[METH], col = COL[METH], lty = LTY[METH], lwd = LWD[METH],
       pch = PCH[METH], ncol = 5, bty = "n", cex = 0.9, seg.len = 2.5)
dev.off()

# ---------------------------------------------------------------- Figure 2
aucc <- c("auc_l1_dev", "auc_quad_dev", "auc_knn_dev", "auc_rf_dev")
pool <- rbind(E2[E2$beta > 0, c("mca", aucc)], E3[E3$k == 3, c("mca", aucc)], E6[, c("mca", aucc)])
best <- apply(pool[aucc], 1, max, na.rm = TRUE)
rej <- pool$mca <= alpha
m <- 150; M <- 450                       # k = 3, n = 300, tau = 1/2
sdA <- sqrt((M + 1) / (12 * m * (M - m)))
thr <- 0.5 + qnorm(1 - alpha / 12) * sdA # Bonferroni value for 3 populations x 4 learners
br <- c(-Inf, seq(0.50, 0.70, by = 0.01), Inf)
bin <- cut(best, br, right = FALSE)
nb <- tapply(rej, bin, length); pb <- tapply(rej, bin, mean)
mid <- c(0.495, seq(0.505, 0.695, by = 0.01), 0.71)
ok <- !is.na(pb) & nb >= 10
ci <- t(sapply(which(ok), function(j) binom.test(sum(rej[which(bin == levels(bin)[j])]), nb[j])$conf.int))

pdf(file.path(figdir, "fig2_auc.pdf"), width = 9, height = 3.6)
par(mfrow = c(1, 2), mar = c(3.6, 3.4, 1.2, 0.8), mgp = c(2.2, 0.6, 0), las = 1, cex.axis = 0.85)
plot(mid[ok], pb[ok], ylim = c(0, 1), xlim = c(0.49, 0.72), pch = 19, cex = 0.8,
     xlab = "Largest test AUC of the four learners", ylab = "Power of MCA-kDET")
arrows(mid[ok], ci[, 1], mid[ok], ci[, 2], angle = 90, code = 3, length = 0.02, col = "grey40")
lines(mid[ok], pb[ok])
abline(v = thr, lty = 2); abline(h = alpha, lty = 3, col = "grey60")
mtext(sprintf("%.3f", thr), side = 1, at = thr, line = -1.2, cex = 0.7, adj = -0.1)

a6 <- E6[E6$model == "A", ]
dd <- sort(unique(a6$d))
mm <- c("mca", "mca_l1", "disco", "mmd", "mst")
P6 <- sapply(mm, function(p) sapply(dd, function(v) pw(a6[a6$d == v, p])))
plot(range(dd), c(0, 1), type = "n", log = "x", xaxt = "n", xlab = expression(italic(d)), ylab = "Power")
axis(1, at = dd, labels = dd)
abline(h = alpha, lty = 3, col = "grey60")
for (p in mm) { lines(dd, P6[, p], col = COL[p], lty = LTY[p], lwd = LWD[p])
  points(dd, P6[, p], col = COL[p], pch = ifelse(is.na(PCH[p]), 17, PCH[p]), cex = 0.8) }
legend("topright", legend = LAB[mm], col = COL[mm], lty = LTY[mm], lwd = LWD[mm],
       pch = ifelse(is.na(PCH[mm]), 17, PCH[mm]), bty = "n", cex = 0.8)
dev.off()

# ---------------------------------------------------------------- numbers for the text
cut3 <- cut(best, c(-Inf, 0.56, 0.58, Inf), right = FALSE)
auc6 <- apply(a6[aucc], 1, max, na.rm = TRUE)
b6 <- cut(auc6, c(0.57, 0.59), right = FALSE)
txt <- c(sprintf("results: %s   replications pooled for Figure 2 (left): %d", outdir, length(rej)),
         sprintf("null sd of the AUC = %.4f; threshold 1/2 + %.3f sd = %.3f", sdA, qnorm(1 - alpha / 12), thr),
         "power by largest AUC (<0.56 / [0.56,0.58) / >=0.58):",
         paste(sprintf("%.3f", tapply(rej, cut3, mean)), collapse = " / "),
         "power in (M1), beta = 0.8, by d (MCA, l1, DISCO, MMD, MST):",
         apply(cbind(dd, round(P6, 2)), 1, paste, collapse = "  "),
         "(M1): power among replications with AUC in [0.57, 0.59), by d:",
         paste(sprintf("d=%d: %.2f (n=%d)", dd, sapply(dd, function(v) mean(a6$mca[a6$d == v & !is.na(b6)] <= alpha)),
                       sapply(dd, function(v) sum(a6$d == v & !is.na(b6)))), collapse = "; "),
         "Figure 1 (methods in the order of the legend):", paste(METH, collapse = "/"), num)
writeLines(txt, file.path(figdir, "fig_numbers.txt"))
cat(txt, sep = "\n")
cat("\nwritten:", file.path(figdir, c("fig1_power.pdf", "fig2_auc.pdf", "fig_numbers.txt")), sep = "\n  ")
