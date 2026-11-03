probe <- "tests/out/parallel15_20261101"
cores <- seq(0L, 28L, 2L)
data <- lapply(cores, function(i) readRDS(file.path(probe, sprintf("cpu%d.rds", i))))
for (x in data) stopifnot(identical(x$p, data[[1L]]$p), identical(x$seed, data[[1L]]$seed))
cat("OK 15 cores x 2 evaluations: identical p-values within and across processes\n")
times <- unlist(lapply(data, `[[`, "seconds"))
cat(sprintf("Concurrent S2 B=200 seconds per evaluation: min %.1f, median %.1f, max %.1f\n",
            min(times), median(times), max(times)))
f <- list.files("tests/out/parallel15_writercheck/libpilot3", pattern="\\.csv$", full.names=TRUE)
stopifnot(length(f)==1L)
z <- read.csv(f, check.names=FALSE)
stopifnot(nrow(z)==1L, z$seed==data[[1L]]$seed, z$B==200L, all(is.na(z$err)|z$err==""))
for(p in names(data[[1L]]$p)) stopifnot(abs(z[[p]]-data[[1L]]$p[[p]]) < 1e-13)
cat("OK disjoint-range writer preserves reference seed and every p-value\n")
