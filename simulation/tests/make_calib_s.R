# =============================================================================
# Calibration of the scenarios S1-S3 by the ideal AUC (task section 3).
# Writes R/calib_params_s.R (does NOT overwrite an existing one unless --force).
# No test is involved; the grid is fixed by the task: ideal AUC 0.54, 0.57,
# 0.60, 0.64, 0.68 (S1a, S1b, S2) and 0.60, 0.64 (S3a, S3b). beta = AUC - 0.5.
#   Rscript tests/make_calib_s.R [--force]
# =============================================================================
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", fa[1L])))
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- "R"; PROFILE <- "quick"
options(power3.progress = TRUE)
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "competitors2.R",
            "rotation.R", "models_s.R", "oracle_s.R", "learners_lib2.R", "experiments.R"))
  source(file.path("R", f))
out <- file.path("R", "calib_params_s.R")
if (file.exists(out) && !"--force" %in% commandArgs(TRUE)) stop(out, " exists; use --force to overwrite")
S_CALIB <- NULL                       # not needed to solve
lev <- S_LEVELS
# every finished (model, level) is cached: R sometimes exits silently on this machine
# (CLAUDE.md); simply run the same command again to resume
cache_f <- file.path("tests", "out", "calib_s_cache_codex_v2.rds")
cache <- if (file.exists(cache_f)) readRDS(cache_f) else list()
fmt <- function(v) paste0("c(", paste(sprintf("%.8g", v), collapse = ", "), ")")
lines <- c("# Ideal-AUC calibration of the scenarios S1-S3 (written by tests/make_calib_s.R, FROZEN).",
           sprintf("# written %s; cB used for S3a = %g; beta = ideal AUC - 0.5", format(Sys.time()), PAR$cB),
           "S_CALIB <- list(")
n_total <- 3L * length(lev) + 4L; n_done <- 0L
for (m in c("S1a", "S1b", "S2", "S3a", "S3b")) {
  lv <- if (m %in% c("S3a", "S3b")) c(0.10, 0.14) else lev
  ent <- character(0)
  for (b in lv) {
    t0 <- proc.time()[["elapsed"]]
    key <- sprintf("%s_%.2f", m, b)
    if (is.null(cache[[key]])) {
      cat(sprintf("START %s target %.2f\n", m, 0.5 + b)); flush.console()
      v <- oracle_param(m, 0.5 + b, par = list(cB = PAR$cB))
      if (any(!is.finite(v))) stop("non-finite calibration: ", key)
      cache[[key]] <- v
      saveRDS(cache, cache_f)
    }
    v <- cache[[key]]
    cat(sprintf("%s level %.2f (ideal AUC %.2f): %s  [%.0fs]\n", m, b, 0.5 + b,
                paste(sprintf("%.6g", v), collapse = " "), proc.time()[["elapsed"]] - t0))
    n_done <- n_done + 1L
    cat(sprintf("  [%s%s] %d/%d
", strrep("#", round(30 * n_done / n_total)),
                strrep("-", 30 - round(30 * n_done / n_total)), n_done, n_total))
    ent <- c(ent, sprintf("    \"%s\" = %s", sprintf("%.2f", b), if (length(v) == 1L) sprintf("%.8g", v) else fmt(v)))
  }
  lines <- c(lines, sprintf("  %s = list(\n%s)%s", m, paste(ent, collapse = ",\n"),
                            if (m == "S3b") "" else ","))
}
lines <- c(lines, ")")
writeLines(lines, out)
cat("written", out, "\n")
