# =============================================================================
# Backward-compatibility check (TASK_学习器库.md, section 3.6)
#
# Re-runs two smoke cells, reps 1-2, with lib = "l1,nn" through the new code
# and compares with results/smoke/ written by the original code:
#   * every p-value column must be bit-identical;
#   * the Schott p-value and every other numeric column (diagnostics) must
#     agree to the 15 significant digits that write.csv keeps (relative
#     tolerance 1e-14); timing columns are skipped.
# Serial (no cluster); writes nothing to results/.
#
#   Rscript tests/test_compat.R
# =============================================================================

script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- file.path(getwd(), "R")
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "experiments.R"))
  source(file.path(SRC_DIR, f))
source(file.path(SRC_DIR, "summarize.R"))   # for P_COLS

CELLS <- c("smoke_B_b1.00_k3_n40_d10_tau0.5_S1",
           "smoke_D_b1.00_k3_n40_d10_tau0.5_S1")
REPS <- 1:2

grid <- make_grid("smoke")
all_ok <- TRUE
# a value written by write.csv keeps 15 significant digits; reading it back
# gives the double nearest to that decimal string, so compare with a relative
# tolerance of 1e-14 (the rounding error is at most 5e-15)
same15 <- function(a, b) isTRUE(all(abs(a - b) <= 1e-14 * pmax(abs(a), abs(b)) | a == b))

for (id in CELLS) {
  cell <- as.list(grid[grid$cell == id, ])
  stopifnot(length(cell$cell) == 1L, cell$lib == "l1,nn", cell$libsets == "")
  old <- read.csv(file.path("results", "smoke", sprintf("%s__r00001-00002.csv", id)),
                  stringsAsFactors = FALSE, check.names = FALSE)
  for (r in REPS) {
    new <- one_rep(cell, r)
    o <- old[old$rep == r, ]
    # an all-empty err column is read back as NA
    if (nrow(o) != 1L || !(is.na(o$err) || o$err == "")) {
      cat(id, "rep", r, ": old row missing or had an error\n"); all_ok <- FALSE; next
    }
    common <- intersect(names(new), names(o))
    pcols <- intersect(P_COLS, common)
    dcols <- setdiff(common[vapply(common, function(x) is.numeric(o[[x]]), logical(1))],
                     c(pcols, grep("^t_", common, value = TRUE), "rep", "seed"))
    # permutation p-values (multiples of 1/(B+1)) must be bit-identical; the
    # Schott p-value is a continuous pnorm value, which the csv keeps to 15 digits
    p_same <- vapply(pcols, function(x)
      if (x == "schott") same15(as.numeric(new[[x]]), o[[x]])
      else identical(as.numeric(new[[x]]), as.numeric(o[[x]])), logical(1))
    d_same <- vapply(dcols, function(x) same15(as.numeric(new[[x]]), o[[x]]), logical(1))
    cat(sprintf("%s rep %d: %d p-value columns identical: %s | %d diagnostics equal (15 digits): %s\n",
                id, r, length(pcols), all(p_same), length(dcols), all(d_same)))
    if (!all(p_same)) {
      bad <- pcols[!p_same]
      print(data.frame(col = bad, old = unlist(o[bad]), new = unlist(new[bad])), row.names = FALSE)
    }
    if (!all(d_same)) {
      bad <- dcols[!d_same]
      print(data.frame(col = bad, old = unlist(o[bad]), new = unlist(new[bad])), row.names = FALSE)
    }
    missing <- setdiff(setdiff(names(o), grep("^t_", names(o), value = TRUE)), c(names(new), names(cell)))
    if (length(missing)) cat("   columns in old but not new:", paste(missing, collapse = ", "), "\n")
    all_ok <- all_ok && all(p_same) && all(d_same)
  }
}
cat(if (all_ok) "\nCOMPAT OK: new code with lib = l1,nn reproduces the old smoke results\n"
    else "\nCOMPAT FAILED\n")
