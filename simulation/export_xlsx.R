# =============================================================================
# Excel versions of the experiment summaries (E1-E8), for reading.
#
#   "D:\R-4.6.1\bin\Rscript.exe" export_xlsx.R --out=results_v2_quick
#   "D:\R-4.6.1\bin\Rscript.exe" export_xlsx.R --out=results_v2_full [--exp=E2,E8]
#
# Runs the existing report scripts (report_v2_E<n>.R; for E2 also
# report_v2_table2.R) unchanged, captures every table they print and writes
#   <out>/xlsx/E<n>_summary.xlsx   one sheet per printed table, plus "notes"
#                                  with the text lines (errors, MC ranges, ...)
#   <out>/xlsx/E1-E8_summary.xlsx  all of them in one workbook
# Numbers are stored as numbers; values the reports flag as outside the Monte
# Carlo range ("!" or "*") are shown in red. Results are not modified.
# Needs the package openxlsx (installed 25 Oct 2026).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
  a <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", key, "="), "", a[1L]) else default
}
outdir <- get_arg("out", "results_v2_quick")
exps <- strsplit(get_arg("exp", paste0("E", 1:8, collapse = ",")), ",")[[1L]]
exps <- exps[dir.exists(file.path(outdir, exps))]
suppressPackageStartupMessages(library(openxlsx))
xdir <- file.path(outdir, "xlsx"); dir.create(xdir, showWarnings = FALSE)

# ---- run one report script, capturing printed tables and cat() text
run_report <- function(script) {
  tabs <- list(); heads <- character(); notes <- character(); last <- ""
  env <- new.env(parent = globalenv())
  env$commandArgs <- function(trailingOnly = FALSE)
    if (trailingOnly) paste0("--out=", outdir) else c(base::commandArgs(FALSE), paste0("--out=", outdir))
  env$cat <- function(..., file = "", sep = " ", fill = FALSE, labels = NULL, append = FALSE) {
    if (identical(file, "")) {
      txt <- paste(vapply(list(...), function(z) paste(format(z), collapse = " "), ""), collapse = sep)
      ln <- trimws(strsplit(txt, "\n", fixed = TRUE)[[1L]]); ln <- ln[nzchar(ln)]
      if (length(ln)) { notes <<- c(notes, ln); last <<- paste(ln, collapse = " ") }
    }
    base::cat(..., file = file, sep = sep, fill = fill, labels = labels, append = append)
  }
  env$print <- function(x, ...) {
    if (is.data.frame(x) || is.matrix(x) || is.table(x)) {
      df <- if (is.data.frame(x)) x else {
        m <- as.matrix(unclass(x)); r <- data.frame(row = rownames(m), check.names = FALSE)
        cbind(r, as.data.frame(m, check.names = FALSE))
      }
      tabs[[length(tabs) + 1L]] <<- df
      heads <<- c(heads, sub(":$", "", last))   # text of the last cat() call before the table
      last <<- ""
    }
    base::print(x, ...)
  }
  txt <- readLines(script, warn = FALSE, encoding = "UTF-8")
  txt <- gsub('"results_v2_quick"', sprintf('"%s"', outdir), txt, fixed = TRUE)   # E2 has the path fixed
  for (e in parse(text = txt)) eval(e, env)
  list(tabs = tabs, heads = heads, notes = notes)
}

# ---- numeric conversion: a column is numeric if every entry is a number,
# optionally followed by a flag ("!", "*", "$^*$"), or a missing-value mark
num_re <- "^-?[0-9]*\\.?[0-9]+(e[-+]?[0-9]+)?$"
flag_re <- "(!|\\*|\\$\\^\\*\\$)$"
convert <- function(df) {
  flags <- matrix(FALSE, nrow(df), ncol(df))
  for (j in seq_along(df)) {
    v <- trimws(as.character(df[[j]]))
    core <- sub(flag_re, "", v)
    miss <- v %in% c("-", "--", "\u2013", "NA", "NaN", "")
    if (all(grepl(num_re, core) | miss) && any(!miss)) {
      flags[, j] <- grepl(flag_re, v) & !miss
      df[[j]] <- suppressWarnings(as.numeric(ifelse(miss, NA, core)))
    }
  }
  # E8 prints the Monte Carlo range per row ("range95" = "[lo, hi]") instead of
  # flagging: mark the rejection rates outside it (every numeric column but reps)
  if ("range95" %in% names(df)) {
    lim <- do.call(rbind, lapply(regmatches(df$range95, gregexpr("-?[0-9.]+", df$range95)), as.numeric))
    for (j in which(vapply(df, is.numeric, TRUE) & names(df) != "reps"))
      flags[, j] <- flags[, j] | (!is.na(df[[j]]) & (df[[j]] < lim[, 1L] | df[[j]] > lim[, 2L]))
  }
  list(df = df, flags = flags)
}

short_head <- function(h) {
  if (grepl("^MCA-kDET of the Section 5", h)) return("MCA flt_blk levels and KS")
  if (grepl("^Single learners and variants", h)) return("single learners, 4 versions")
  if (grepl("^Diagnostics", h)) return("diagnostics")
  trimws(sub("\\s*\\(.*$", "", h))                  # drop a trailing parenthetical
}

sheet_name <- function(s, used) {
  for (ch in c("[", "]", ":", "*", "?", "/", "\\")) s <- gsub(ch, " ", s, fixed = TRUE)   # not allowed in sheet names
  s <- trimws(gsub("\\s+", " ", s))
  s <- substr(s, 1L, 31L); b <- s; i <- 2L
  while (tolower(s) %in% tolower(used)) { s <- paste0(substr(b, 1L, 28L), " ", i); i <- i + 1L }
  s
}

st_head <- createStyle(textDecoration = "bold", fgFill = "#DDE7F3", border = "bottom")
st_flag <- createStyle(fontColour = "#B00020", textDecoration = "bold", fgFill = "#FDE7E9")
st_wrap <- createStyle(wrapText = FALSE)

add_tables <- function(wb, res, prefix, used) {
  for (i in seq_along(res$tabs)) {
    nm <- if (i == 1L) "main" else if (nzchar(res$heads[i])) short_head(res$heads[i]) else paste("table", i)
    nm <- sheet_name(paste0(prefix, nm), used); used <- c(used, nm)
    cv <- convert(res$tabs[[i]])
    addWorksheet(wb, nm)
    writeData(wb, nm, cv$df, headerStyle = st_head)
    freezePane(wb, nm, firstRow = TRUE)
    setColWidths(wb, nm, cols = seq_along(cv$df), widths = "auto")
    fl <- which(cv$flags, arr.ind = TRUE)
    if (nrow(fl)) addStyle(wb, nm, st_flag, rows = fl[, 1L] + 1L, cols = fl[, 2L], gridExpand = FALSE, stack = TRUE)
    if (!nzchar(prefix)) base::cat(sprintf("   sheet %-30s %3d x %2d, flagged cells: %d\n", nm, nrow(cv$df), ncol(cv$df), nrow(fl)))
  }
  if (length(res$notes)) {
    nm <- sheet_name(paste0(prefix, "notes"), used); used <- c(used, nm)
    addWorksheet(wb, nm)
    writeData(wb, nm, data.frame(notes = res$notes), headerStyle = st_head)
    setColWidths(wb, nm, cols = 1L, widths = 120)
  }
  used
}

wb_all <- createWorkbook(); used_all <- character()
for (e in exps) {
  scripts <- c(sprintf("report_v2_%s.R", e), if (e == "E2") "report_v2_table2.R")
  scripts <- scripts[file.exists(scripts)]
  res <- list(tabs = list(), heads = character(), notes = character())
  for (s in scripts) {
    cat(sprintf("\n===== %s: %s =====\n", e, s))
    r <- run_report(s)
    if (s == "report_v2_table2.R" && length(r$heads)) r$heads[] <- "Table 2 variants"
    res$tabs <- c(res$tabs, r$tabs); res$heads <- c(res$heads, r$heads); res$notes <- c(res$notes, r$notes)
  }
  wb <- createWorkbook(); add_tables(wb, res, "", character())
  saveWorkbook(wb, file.path(xdir, sprintf("%s_summary.xlsx", e)), overwrite = TRUE)
  used_all <- add_tables(wb_all, res, paste0(e, " "), used_all)
  cat(sprintf("-> %s: %d tables\n", file.path(xdir, sprintf("%s_summary.xlsx", e)), length(res$tabs)))
}
saveWorkbook(wb_all, file.path(xdir, "E1-E8_summary.xlsx"), overwrite = TRUE)
cat("\nXLSX DONE:", file.path(xdir, "E1-E8_summary.xlsx"), "| experiments:", paste(exps, collapse = ","), "\n")
