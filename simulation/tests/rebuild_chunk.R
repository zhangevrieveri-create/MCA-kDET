# =============================================================================
# Rebuild a chunk file from its .debug.rds (written by run_cell when the main
# process failed to assemble the rows), with the same code as run_cell and
# without recomputing anything. Checks reps, seeds, errors, columns and cell
# settings against a sibling chunk of the same cell before writing.
#
#   Rscript tests/rebuild_chunk.R <path to the .csv that is missing>
# =============================================================================
script_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1L]))) else getwd()
})
setwd(file.path(script_dir, ".."))
SRC_DIR <- "R"; source("R/sim.R")
f <- commandArgs(TRUE)[1L]
if (file.exists(f)) stop("chunk already exists: ", f)
dbg <- readRDS(paste0(f, ".debug.rds"))
cell <- dbg$cell; rows <- dbg$rows
rr <- as.integer(regmatches(f, regexpr("[0-9]{5}-[0-9]{5}", f)) |> strsplit("-") |> unlist())
cpart <- cell[setdiff(names(cell), c("exp", "cell"))]
cpart <- lapply(cpart, function(v) rep(unname(unlist(v))[1L], length(rows)))
df <- cbind(data.frame(cpart, stringsAsFactors = FALSE, check.names = FALSE), rows_to_df(rows))
df$err[is.na(df$err)] <- ""
sibs <- setdiff(list.files(dirname(f), pattern = paste0("^", cell$cell, "__r.*\\.csv$"), full.names = TRUE), f)
if (!length(sibs)) stop("no sibling chunk to compare with")
sib <- read.csv(sibs[1L], check.names = FALSE, stringsAsFactors = FALSE)
norm <- function(v) { v <- as.character(v); v[is.na(v)] <- ""; v }
ok <- c(
  reps    = identical(sort(as.integer(df$rep)), rr[1]:rr[2]),
  seeds   = all(df$seed == cell$seed0 + df$rep),
  noerr   = all(df$err == ""),
  columns = setequal(names(df), names(sib)),
  setting = all(vapply(names(cpart), function(k) norm(df[[k]][1]) == norm(sib[[k]][1]), logical(1))))
print(ok)
cat("retries:", sum(df$retries), "| compared with", basename(sibs[1L]), "\n")
if (!all(ok)) stop("checks failed; nothing written")
write.csv(df[, names(sib)], paste0(f, ".tmp"), row.names = FALSE)
file.rename(paste0(f, ".tmp"), f)
cat("REBUILT:", basename(f), "\n")
