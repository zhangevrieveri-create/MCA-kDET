# Financial application; approximate week-label calibration, frozen learning/statistic code.
args <- commandArgs(TRUE)
arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit[1L]) else default
}
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), winslash = "/"))
src <- normalizePath(file.path(root, "..", "..", "R"), winslash = "/")
plan <- jsonlite::fromJSON(file.path(root, "study_plan.json"))
name <- arg("config", "primary"); smoke <- name == "smoke"
config <- as.list(plan$specifications[plan$specifications$name == if (smoke) "primary" else name, ])
stopifnot(length(config$name) == 1L)
if (smoke) { config$dimension <- 3L; config$permutations <- 19L }
lib <- if (smoke) "l1" else c("l1", "quad", "knn", "rf")
for (pkg in c("glmnet", if (!smoke) "ranger")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Package unavailable: ", pkg)
for (file in c("learners.R", "learners_lib.R", "tests.R")) source(file.path(src, file))
output <- file.path(root, "results", name); dir.create(output, recursive = TRUE, showWarnings = FALSE)
input <- file.path(root, "prepared")
fingerprint <- list(config = config, lib = lib, md5 = tools::md5sum(c(
  file.path(input, c("weeks.csv", "weekly_returns.csv", "filtered_returns.csv", "garch_parameters.csv")),
  file.path(root, c("study_plan.json", "run_analysis.R")), file.path(src, c("learners.R", "learners_lib.R", "tests.R")))))
result_file <- file.path(output, "analysis.rds")
if (file.exists(result_file)) {
  if (!identical(readRDS(result_file)$fingerprint, fingerprint)) stop("Checkpoint fingerprint differs; use a versioned output")
  cat("Verified checkpoint: ", name, "\n", sep = ""); quit(status = 0)
}
started <- proc.time()[["elapsed"]]
weeks <- read.csv(file.path(input, "weeks.csv"), check.names = FALSE)
filename <- if (config$representation == "filtered") "filtered_returns.csv" else "weekly_returns.csv"
dat <- read.csv(file.path(input, filename), check.names = FALSE, na.strings = c("", "NA"))
stopifnot(identical(weeks$week_end, dat$week_end))
selected <- weeks$period == "analysis"
weeks <- weeks[selected, , drop = FALSE]
X <- as.matrix(dat[selected, seq.int(2L, config$dimension + 1L), drop = FALSE])
reg <- if (config$regime == "fixed_15_25") weeks$regime_fixed else weeks$regime
usable <- rowSums(is.finite(X)) == ncol(X)
Tn <- nrow(X); set.seed(config$seed)
if (config$split == "random") {
  role <- rep("missing", Tn)
  for (i in 1:3) {
    ix <- which(usable & reg == i); ix <- ix[sample.int(length(ix))]
    nf <- floor(length(ix) / 2)
    role[ix[seq_len(nf)]] <- "fit"; role[ix[-seq_len(nf)]] <- "test"
  }
} else {
  bid <- ceiling(seq_len(Tn) / config$block); nb <- max(bid)
  for (attempt in 1:1000) {
    fit_blocks <- sample.int(nb, round(nb / 2))
    role <- ifelse(bid %in% fit_blocks, "fit", "test")
    original_role <- role
    cs <- c(0L, cumsum(role == "fit"))
    lo <- pmax(1L, seq_len(Tn) - config$buffer); hi <- pmin(Tn, seq_len(Tn) + config$buffer)
    role[role == "test" & (cs[hi + 1L] - cs[lo]) > 0] <- "buffer"
    role[!usable] <- "missing"
    if (all(tabulate(reg[role == "fit"], 3L) >= 30L) && all(tabulate(reg[role == "test"], 3L) >= 30L)) break
    if (attempt == 1000L) stop("No admissible block split")
  }
  # Calendar, not compressed-valid-row spacing. All offsets reuse the main split.
  if (config$thin > 1L)
    role[role == "test" & (seq_len(Tn) - 1L) %% config$thin != config$offset] <- "thinned"
}
stopifnot(all(tabulate(reg[role == "test"], 3L) >= 5L))
weeks$group_number <- reg; weeks$role <- role
write.csv(weeks, file.path(output, "week_allocation.csv"), row.names = FALSE)
sp <- lapply(1:3, function(i) list(fit = X[reg == i & role == "fit", , drop = FALSE],
                                   test = X[reg == i & role == "test", , drop = FALSE]))
test_indices <- unlist(lapply(1:3, function(i) which(reg == i & role == "test")))
test <- weeks[test_indices, , drop = FALSE]
f <- vapply(sp, function(s) nrow(s$fit), integer(1)); n <- vapply(sp, function(s) nrow(s$test), integer(1))
cat(name, ": fit=", paste(f, collapse = "/"), "; test=", paste(n, collapse = "/"), "\n", sep = ""); flush.console()
set.seed(config$seed)
warnings_seen <- character()
sc <- withCallingHandlers(fit_scores(sp, groups = "ovr", lib = lib),
  warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") })
set.seed(config$seed + 100000L)
Lall <- perm_labels(sc$L, config$permutations)
Tl <- lapply(sc$Sovr, function(S) ovr_T(S, 1:3, Lall)); Tlib <- do.call(rbind, Tl)
statistics <- list(MCA = apply(Tlib, 2L, max))
for (a in lib) statistics[[paste0("MCA_", a)]] <- apply(Tl[[a]], 2L, max)
statistics$Two_sided <- apply(abs(Tlib), 2L, max)
statistics$Sum <- colSums(Tlib)
p <- lapply(statistics, pval_perm); p$Comb <- pval_comb(statistics$MCA, statistics$Sum)
Xfit <- do.call(rbind, lapply(sp, `[[`, "fit"))
Ztest <- apply_std(do.call(rbind, lapply(sp, `[[`, "test")), standardize_fit(Xfit))
D2 <- sqdist(Ztest); D <- sqrt(D2); nb <- ncol(Lall)
G <- group_indicator(Lall, 3L)
statistics$DISCO <- disco_stat(D, G, n, 3L, nb)
h2 <- median(D2[upper.tri(D2)]); K <- exp(-D2 / (2 * h2))
statistics$MMD <- mmd_stat(K, G, n, 3L, nb)
edges <- mst_edges(D)
statistics$MST <- colSums(Lall[edges[, 1L], , drop = FALSE] == Lall[edges[, 2L], , drop = FALSE])
for (method in c("DISCO", "MMD", "MST")) p[[method]] <- pval_perm(statistics[[method]])
for (method in names(statistics)) stopifnot(abs(p[[method]] - mean(statistics[[method]] >= statistics[[method]][1L] - 1e-10)) < 1e-14)
values <- data.frame(method = names(p), p_value = unlist(p), row.names = NULL)
values$dimension <- config$dimension; values$fit_weeks <- sum(f); values$test_weeks <- sum(n)
values$permutations <- config$permutations; values$calibration <- "approximate_week_label"
values$MC_SE <- sqrt(values$p_value * (1 - values$p_value) / nb)
write.csv(values, file.path(output, "pvalues.csv"), row.names = FALSE)
write.csv(data.frame(assignment = seq_len(nb), statistics), file.path(output, "permutation_statistics.csv"), row.names = FALSE)
con <- gzfile(file.path(output, "permutation_labels.csv.gz"), "wt")
write.csv(t(Lall), con, row.names = FALSE); close(con)
score_table <- test[, c("week_end", "group_number", "lagged_vix")]
for (a in lib) for (i in 1:3) score_table[[paste0(a, "_target", i)]] <- sprintf("%.17g", sc$Sovr[[a]][, i])
write.csv(score_table, file.path(output, "test_scores_exact.csv"), row.names = FALSE)
auc <- do.call(rbind, lapply(lib, function(a) do.call(rbind, lapply(1:3, function(i)
  data.frame(learner = a, target = i, auc = auc_vec(sc$Sovr[[a]][, i], sc$L == i), standardized_T = Tl[[a]][i, 1L])))))
write.csv(auc, file.path(output, "auc.csv"), row.names = FALSE)
counts <- data.frame(regime = 1:3, fitting = f, testing = n)
write.csv(counts, file.path(output, "sample_counts.csv"), row.names = FALSE)
writeLines(unique(warnings_seen), file.path(output, "warnings.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
saveRDS(list(fingerprint = fingerprint, config = config, result = values, sc = sc,
             Lall = Lall, statistics = statistics, auc = auc, Xfit = Xfit, Ztest = Ztest,
             edges = edges, weeks = weeks, test = test, warnings = warnings_seen,
             seconds = proc.time()[["elapsed"]] - started), result_file)
print(values[values$method %in% c("MCA", "DISCO", "MMD", "MST"), ], row.names = FALSE)
cat("Completed ", name, " in ", round(proc.time()[["elapsed"]] - started, 1), " seconds\n", sep = "")
