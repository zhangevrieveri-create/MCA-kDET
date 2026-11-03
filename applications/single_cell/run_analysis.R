# Application wrapper: frozen learners/statistics, exhaustive animal-block calibration.
argv <- commandArgs(TRUE)
arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), argv, value = TRUE)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit[1]) else default
}
filearg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
root <- dirname(normalizePath(sub("^--file=", "", filearg), winslash = "/"))
src <- normalizePath(file.path(root, "..", "..", "R"), winslash = "/")
input <- arg("input", file.path(root, "prepared", "analysis"))
output <- arg("output", file.path(root, "results", "primary"))
d <- as.integer(arg("dimension", "100"))
seed <- as.integer(arg("seed", "20261027"))
lib <- strsplit(arg("lib", "l1,quad,knn,rf"), ",", fixed = TRUE)[[1]]
do_pair <- identical(arg("pair", "yes"), "yes")
for (pkg in c("glmnet", if ("rf" %in% lib) "ranger")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package not available: ", pkg)
}
for (name in c("learners.R", "learners_lib.R", "tests.R")) source(file.path(src, name))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
inputs <- c(file.path(input, c("samples.csv", "expression.csv.gz", "genes.csv",
                              "animal_label_assignments.csv", "design.json")),
            file.path(src, c("learners.R", "learners_lib.R", "tests.R")),
            file.path(root, "run_analysis.R"))
fingerprint <- list(md5 = tools::md5sum(inputs), dimension = d, seed = seed,
                    lib = lib, pair = do_pair)
result_file <- file.path(output, "analysis.rds")
if (file.exists(result_file)) {
  old <- readRDS(result_file)
  if (!identical(fingerprint, old$fingerprint)) stop("Existing result has a different fingerprint; use a new output directory")
  cat("Verified checkpoint: ", output, "\n", sep = "")
  quit(status = 0)
}
started <- proc.time()[["elapsed"]]
samples <- read.csv(file.path(input, "samples.csv"), check.names = FALSE)
con <- gzfile(file.path(input, "expression.csv.gz"), "rt")
dat <- read.csv(con, check.names = FALSE); close(con)
stopifnot(identical(samples$sampleID, dat$sampleID), d <= ncol(dat) - 1L)
X <- as.matrix(dat[, seq.int(2L, d + 1L), drop = FALSE])
stopifnot(all(is.finite(X)))
fit_ids <- unique(samples$verified_animal_block[samples$split == "fit"])
test_ids <- unique(samples$verified_animal_block[samples$split == "test"])
stopifnot(!length(intersect(fit_ids, test_ids)))
sp <- lapply(1:3, function(i) list(
  fit = X[samples$split == "fit" & samples$group_number == i, , drop = FALSE],
  test = X[samples$split == "test" & samples$group_number == i, , drop = FALSE]))
test_rows <- unlist(lapply(1:3, function(i) which(samples$split == "test" & samples$group_number == i)))
test <- samples[test_rows, , drop = FALSE]
assignments <- as.matrix(read.csv(file.path(input, "animal_label_assignments.csv"), check.names = FALSE))
units <- match(test$verified_animal_block, colnames(assignments))
stopifnot(!anyNA(units), nrow(unique(assignments)) == nrow(assignments))
Lall <- t(assignments[, units, drop = FALSE])
observed <- which(colSums(Lall == test$group_number) == nrow(test))
stopifnot(length(observed) == 1L)
Lall <- Lall[, c(observed, setdiff(seq_len(ncol(Lall)), observed)), drop = FALSE]
nb <- ncol(Lall); n <- vapply(sp, function(s) nrow(s$test), integer(1))
stopifnot(all(vapply(1:3, function(i) all(colSums(Lall == i) == n[i]), logical(1))))
stopifnot(all(Lall[, 1L] == test$group_number))
# pval_perm is exact here: observed once at position 1, remaining assignments once each.
stopifnot(identical(pval_perm(rep(1, nb)), 1))
set.seed(seed)
cat("Fit: d=", d, ", fit/test cells=", sum(vapply(sp, function(s) nrow(s$fit), integer(1))),
    "/", nrow(test), ", complete assignments=", nb, "\n", sep = ""); flush.console()
warnings_seen <- character()
sc <- withCallingHandlers(
  fit_scores(sp, groups = c("ovr", if (do_pair) "pair"), lib = lib),
  warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") })
stopifnot(identical(sc$L, as.integer(test$group_number)))
Tl <- lapply(sc$Sovr, function(S) ovr_T(S, 1:3, Lall))
Tlib <- do.call(rbind, Tl)
statistics <- list(MCA = apply(Tlib, 2L, max))
for (a in lib) statistics[[paste0("MCA_", a)]] <- apply(Tl[[a]], 2L, max)
statistics$Two_sided <- apply(abs(Tlib), 2L, max)
statistics$Sum <- colSums(Tlib)
p <- lapply(statistics, pval_perm)
p$Comb <- pval_comb(statistics$MCA, statistics$Sum)
if (do_pair) {
  Pl <- lapply(sc$Ppair, function(S) pair_T(S, sc$pairs, Lall))
  statistics$Pair <- apply(do.call(rbind, Pl), 2L, max)
  p$Pair <- pval_perm(statistics$Pair)
}
Zfit <- do.call(rbind, lapply(sp, `[[`, "fit"))
Ztest <- apply_std(do.call(rbind, lapply(sp, `[[`, "test")), standardize_fit(Zfit))
D2 <- sqdist(Ztest); D <- sqrt(D2)
G <- group_indicator(Lall, 3L)
statistics$DISCO <- disco_stat(D, G, n, 3L, nb)
h2 <- median(D2[upper.tri(D2)])
stopifnot(is.finite(h2), h2 > 0)
statistics$MMD <- mmd_stat(exp(-D2 / (2 * h2)), G, n, 3L, nb)
E <- mst_edges(D)
statistics$MST <- colSums(Lall[E[, 1L], , drop = FALSE] == Lall[E[, 2L], , drop = FALSE])
for (a in c("DISCO", "MMD", "MST")) p[[a]] <- pval_perm(statistics[[a]])
# Independent upper-tail count cross-check, without the +1 shortcut.
for (a in names(statistics)) stopifnot(abs(p[[a]] - mean(statistics[[a]] >= statistics[[a]][1] - 1e-10)) < 1e-14)
result <- data.frame(method = names(p), p_value = unlist(p), row.names = NULL)
result$statistic <- vapply(result$method, function(a) if (a %in% names(statistics)) statistics[[a]][1] else NA_real_, numeric(1))
result$dimension <- d; result$fitting_seed <- seed
result$test_animals <- length(test_ids); result$test_cells <- nrow(test)
result$assignments <- nb; result$calibration <- "complete_animal_block"
write.csv(result, file.path(output, "pvalues.csv"), row.names = FALSE)
write.csv(data.frame(assignment = seq_len(nb), observed = seq_len(nb) == 1L, statistics),
          file.path(output, "permutation_statistics.csv"), row.names = FALSE)
auc <- do.call(rbind, lapply(lib, function(a) do.call(rbind, lapply(1:3, function(i)
  data.frame(learner = a, target = i, auc = auc_vec(sc$Sovr[[a]][, i], sc$L == i),
             standardized_T = Tl[[a]][i, 1L])))))
write.csv(auc, file.path(output, "auc.csv"), row.names = FALSE)
score_table <- test[, c("sampleID", "verified_animal_block", "group", "group_number", "geo_tumor", "geo_Sex")]
for (a in lib) for (i in 1:3) score_table[[paste0(a, "_target", i)]] <- sc$Sovr[[a]][, i]
write.csv(score_table, file.path(output, "test_scores.csv"), row.names = FALSE)
timing <- data.frame(component = names(sc$tm), seconds = unlist(sc$tm), row.names = NULL)
timing <- rbind(timing, data.frame(component = "total_including_calibration", seconds = proc.time()[["elapsed"]] - started))
write.csv(timing, file.path(output, "timing.csv"), row.names = FALSE)
writeLines(unique(warnings_seen), file.path(output, "warnings.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
saveRDS(list(fingerprint = fingerprint, result = result, sc = sc, Lall = Lall,
             test_metadata = test, statistics = statistics, auc = auc,
             warnings = warnings_seen, h2 = h2, timing = timing), result_file)
print(result[, c("method", "p_value", "statistic")], row.names = FALSE)
cat("Completed in ", round(proc.time()[["elapsed"]] - started, 1), " seconds\n", sep = "")
