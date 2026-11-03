out <- "motivated_v4_20261102/smoke_B19"
grid <- readRDS(file.path(out,"grid.rds"))
files <- list.files(file.path(out,"rows"),pattern="[.]rds$",full.names=TRUE)
rows <- lapply(files,readRDS)
stopifnot(length(rows)==60L)
for(z in rows) {
  cell <- grid[grid$cell==z$cell,,drop=FALSE]
  stopifnot(nrow(cell)==1L,z$rep%in%1:3,z$B==19L,z$seed==cell$seed0+z$rep)
  pv <- unlist(z$p[setdiff(names(z$p),"kmd_eta")])
  stopifnot(all(is.finite(pv)),all(pv>=0&pv<=1))
}
keys <- vapply(rows,function(z)paste(z$cell,z$rep),character(1))
stopifnot(!anyDuplicated(keys),all(table(vapply(rows,`[[`,character(1),"cell"))==3L))
writeLines(c("60 unique (cell,rep), 20 cells x 3 repeats; every B=19 and every seed=pre-fixed seed0+rep.",
"All actual p-values finite and in [0,1]; signed KMD effect statistic remains in raw records and is excluded from p-value checks.",
"CSV-level verification independently recomputed source MD5 and checked p-value export keys/ranges.",
"Smoke diagnostics are not power estimates."),file.path(out,"integrity_check.txt"))
cat("RAW RECORD SEEDS, COUNTS, B AND P-VALUES VERIFIED\n")
