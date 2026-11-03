# Supplementary Table S9 (E8, size under the GARCH null) as LaTeX.
#   Rscript report_v2_E8_supp.R [--out=results_v2_full]
# Writes <out>/supp_garch.tex; copy it to 写作/latex/supp_garch.tex.
# Panel (a): MCA-kDET, its single learners and variants for the four versions
# (filtered/raw x block/random split); panel (b): DISCO, MMD, MST on filtered
# and raw returns, and the diagnostics.
args <- commandArgs(trailingOnly = TRUE)
outdir <- sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=results_v2_full")[1L])
source("R/summarize.R")
d <- read_exp("E8", outdir); d$err[is.na(d$err)] <- ""
ok <- d[d$err == "", ]
quick <- grepl("quick", outdir)

cells <- sort(unique(ok$cell))
R <- min(table(ok$cell))
lo <- qbinom(.025, R, .05) / R; hi <- qbinom(.975, R, .05) / R
rr <- function(x, p, a = 0.05) mean(x[[p]] <= a, na.rm = TRUE)
f05 <- function(v) paste0(sprintf("%.3f", v), ifelse(v < lo | v > hi, "$^*$", ""))
f3 <- function(v) sub("^-", "$-$", sprintf("%.3f", v))
setting <- function(x) {
  d0 <- x$d[1L]
  ab <- if (grepl("p99$", x$cell[1L])) "0.99" else "0.95"
  c(as.character(d0), ab)
}

vers <- list(c("Filtered", "Blocks", ""), c("Filtered", "Random", "flt_iid_"),
             c("Raw", "Blocks", "raw_blk_"), c("Raw", "Random", "raw_iid_"))
lines_a <- character(); lines_b <- character()
for (j in seq_along(cells)) {
  x <- ok[ok$cell == cells[j], ]
  st <- setting(x)
  for (v in seq_along(vers)) {
    pre <- vers[[v]][3L]
    p <- function(s) paste0(pre, s)
    row <- c(if (v == 1L) st else c("", ""), vers[[v]][1:2],
             f05(rr(x, p("mca"))), f3(rr(x, p("mca"), .01)), f3(rr(x, p("mca"), .10)),
             vapply(c("mca_l1", "mca_quad", "mca_knn", "mca_rf", "mca_2s", "mca_sum", "mca_comb"),
                    function(s) f05(rr(x, p(s))), ""))
    lines_a <- c(lines_a, paste0(paste(row, collapse = " & "), " \\\\"))
  }
  if (j < length(cells)) lines_a <- c(lines_a, "\\addlinespace")
  rowb <- c(st,
            vapply(c("disco", "mmd", "mst", "raw_disco", "raw_mmd", "raw_mst"), function(s) f05(rr(x, s)), ""),
            f3(mean(x$sqacf_raw)), f3(mean(x$sqacf_flt)), f3(mean(x$persist_hat)),
            sprintf("%.0f", mean(x$n_fit_blk)), sprintf("%.0f", mean(x$n_test_blk)),
            sprintf("%.0f", mean(x$n_drop_blk)), sprintf("%d", as.integer(sum(x$filter_nonconv))),
            sprintf("%.1f", mean(x$t_total)))
  lines_b <- c(lines_b, paste0(paste(rowb, collapse = " & "), " \\\\"))
}

qn <- if (quick) sprintf(" \\Qnote{%d replications per setting; asterisks mark values outside $[%.3f, %.3f]$}", R, lo, hi) else
  sprintf(" Each entry is based on %d replications; asterisks mark values outside the 95\\%% Monte Carlo range $[%.4f, %.4f]$.", R, lo, hi)
tex <- c(
  "\\begin{table}[p]",
  "\\renewcommand{\\baselinestretch}{1}\\footnotesize",
  paste0("\\caption{Rejection rates under the GARCH null of Section 4.2: $T = 900$ weeks, $k = 3$ regimes, $B = 200$. ",
         "Data: coordinates standardized by GARCH(1,1) fits, or raw returns; split: blocks of 26 weeks with a buffer of 4 weeks, or random within regimes (Algorithm 1). ",
         "Panel (a): MCA-kDET at $\\alpha = 0.05$, 0.01 and 0.10, and at $\\alpha = 0.05$ its single learners and the two-sided (2s), sum-type and combined variants. ",
         "Panel (b): DISCO, MMD and MST at $\\alpha = 0.05$ on the filtered and the raw returns; lag-1 autocorrelation of the squared returns before and after filtering (mean over coordinates); ",
         "estimated persistence $\\hat a + \\hat b$ (mean over coordinates); mean numbers of fitting, test and dropped weeks in the block split; number of filter fits that did not converge; ",
         "seconds per replication for all tests.", qn, "\\label{tab:S-garch}}"),
  "\\begin{center}",
  "\\setlength{\\tabcolsep}{3pt}",
  "(a) MCA-kDET\\\\[2pt]",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{cccccccccccccc}",
  "\\toprule",
  "$d$ & $a + b$ & Data & Split & MCA & $\\alpha{=}.01$ & $\\alpha{=}.10$ & l1 & quad & knn & rf & 2s & Sum & Comb \\\\",
  "\\midrule",
  lines_a,
  "\\bottomrule",
  "\\end{tabular}}",
  "\\\\[10pt]",
  "(b) Competitors and diagnostics\\\\[2pt]",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{cccccccccccccccc}",
  "\\toprule",
  " & & \\multicolumn{3}{c}{Filtered} & \\multicolumn{3}{c}{Raw} & \\multicolumn{2}{c}{Autocorr.\\ of squares} & & \\multicolumn{3}{c}{Weeks (block split)} & & \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-10}\\cmidrule(lr){12-14}",
  "$d$ & $a + b$ & DISCO & MMD & MST & DISCO & MMD & MST & Raw & Filtered & $\\hat a + \\hat b$ & Fit & Test & Dropped & Non-conv. & s/rep \\\\",
  "\\midrule",
  lines_b,
  "\\bottomrule",
  "\\end{tabular}}",
  "\\end{center}",
  "\\end{table}")
writeLines(tex, file.path(outdir, "supp_garch.tex"), useBytes = TRUE)
cat("written:", file.path(outdir, "supp_garch.tex"), "| cells:", length(cells), "| reps per cell:", R, "\n")
