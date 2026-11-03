# Disjoint repetition ranges for parallel workers. Uses the frozen one_rep/safe_rep.
args <- commandArgs(TRUE)
parg <- function(name) {
  x <- grep(paste0("^--", name, "="), args, value=TRUE)
  stopifnot(length(x) == 1L)
  sub(paste0("^--", name, "="), "", x)
}
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for (p in c("glmnet", "ranger", "gTestsMulti")) stopifnot(requireNamespace(p, quietly=TRUE))
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "garch_ext.R",
            "competitors2.R", "rotation.R", "learners_lib2.R", "models_s.R", "experiments.R"))
  source(file.path(SRC_DIR, f))
out <- parg("out"); id <- parg("cell")
first <- as.integer(parg("start")); last <- as.integer(parg("end"))
g <- readRDS(file.path(out, "formal_grid.rds"))
stopifnot(identical(g, make_grid("libpilot3")))
cell <- as.list(g[g$cell == id, ])
stopifnot(length(cell$cell) == 1L, cell$B == 200L, first >= 1L, last <= cell$reps)
for (r in first:last) {
  target <- file.path(out, "libpilot3", sprintf("%s__r%05d-%05d.csv", id, r, r))
  if (file.exists(target)) {
    z <- read.csv(target, check.names=FALSE)
    stopifnot(nrow(z)==1L, all(is.na(z$err) | z$err==""), z$seed==cell$seed0+r, z$B==200L)
    next
  }
  cat(sprintf("[%s] START %s rep %d\n", format(Sys.time(), "%H:%M:%S"), id, r)); flush.console()
  result <- safe_rep(cell, r)
  if (!identical(result$err, "")) {
    saveRDS(result, paste0(target, ".failed.rds"))
    writeLines(c(id, as.character(r), result$err), file.path(out, "STOP_error.txt"))
    stop("Repeated calculation error; see failed.rds")
  }
  cpart <- cell[setdiff(names(cell), c("exp", "cell"))]
  df <- cbind(data.frame(cpart, stringsAsFactors=FALSE, check.names=FALSE), rows_to_df(list(result)))
  stopifnot(nrow(df)==1L, df$seed==cell$seed0+r, df$B==200L)
  write.csv(df, paste0(target, ".tmp"), row.names=FALSE)
  stopifnot(!file.exists(target), file.rename(paste0(target, ".tmp"), target))
  cat(sprintf("[%s] SAVED %s rep %d (%.1fs)\n", format(Sys.time(), "%H:%M:%S"), id, r, result$t_total)); flush.console()
}
cat("PARALLEL CHUNK OK\n")
