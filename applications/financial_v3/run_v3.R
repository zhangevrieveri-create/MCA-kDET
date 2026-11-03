# =============================================================================
# Chapter 5, financial application v3: runs the jobs listed in a jobs file.
# Calls the frozen learners and statistics of sim/R (nothing there is changed).
#
#   "D:\R-4.6.1\bin\Rscript.exe" run_v3.R --jobs=jobs_calibration.csv --cores=30
#   "D:\R-4.6.1\bin\Rscript.exe" run_v3.R --jobs=jobs_main.csv --cores=30
#   "D:\R-4.6.1\bin\Rscript.exe" run_v3.R --jobs=jobs_main.csv --cores=1 --smoke   (2 jobs, B = 19, l1 only)
#
# A job whose output folder already contains pvalues.csv is skipped, so the
# same command resumes an interrupted run. Afterwards all p-values are
# collected in results/pvalues_<jobs file>.csv.
#
# Amendment 1 (protocol_amendment1.json), wrapper layer only:
#   - a job with strata = L > 0 permutes the test-week labels only within
#     calendar windows of L consecutive weeks of the analysis period
#     (strata = 0 or no strata column: the frozen perm_labels, as before);
#   - pvalues.csv also records strata and stat_obs, the observed statistic of
#     each method (Schott: its z statistic; Comb: NA), used by the rotation test.
# =============================================================================
args <- commandArgs(TRUE)
arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", key, "="), "", hit[1L]) else default
}
smoke <- "--smoke" %in% args
root <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), winslash = "/"))
src <- normalizePath(file.path(root, "..", "..", "R"), winslash = "/")
jobs_file <- arg("jobs", "jobs_main.csv")
cores <- as.integer(arg("cores", "30"))
jobs <- read.csv(file.path(root, jobs_file), stringsAsFactors = FALSE)
if (smoke) { jobs <- jobs[seq_len(min(2L, nrow(jobs))), ]; jobs$B <- 19L; jobs$out <- paste0("smoke/", jobs$out) }
for (pkg in c("glmnet", "ranger")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Package unavailable: ", pkg)
if (!"strata" %in% names(jobs)) jobs$strata <- 0L

# within-window permutation (amendment 1): column 1 is L itself, as in perm_labels;
# h gives the calendar window of each pooled test week
perm_labels_strata <- function(L, B, h) {
  ix <- Filter(function(s) length(s) > 1L, split(seq_along(L), h))
  cbind(L, vapply(seq_len(B), function(b) {
    Lb <- L
    for (s in ix) Lb[s] <- L[s[sample.int(length(s))]]
    Lb }, integer(length(L))))
}

run_job <- function(job, root, src, smoke) {
  out <- file.path(root, job$out)
  if (file.exists(file.path(out, "pvalues.csv"))) return("skipped")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  started <- proc.time()[["elapsed"]]
  lib <- if (smoke) "l1" else c("l1", "quad", "knn", "rf")
  X <- as.matrix(read.csv(file.path(root, job$x_file), check.names = FALSE))
  wk <- read.csv(file.path(root, job$role_file), stringsAsFactors = FALSE)
  lab <- read.csv(file.path(root, job$label_file), check.names = FALSE)
  reg <- as.integer(lab[[job$label_col]]); role <- wk$role
  stopifnot(nrow(X) == length(reg), length(role) == length(reg), all(is.finite(X)))
  k <- 3L
  sp <- lapply(seq_len(k), function(i) list(fit = X[reg == i & role == "fit", , drop = FALSE],
                                            test = X[reg == i & role == "test", , drop = FALSE]))
  f <- vapply(sp, function(s) nrow(s$fit), integer(1)); n <- vapply(sp, function(s) nrow(s$test), integer(1))
  if (any(f < 5L) || any(n < 5L)) stop("too few weeks in a group: fit ", paste(f, collapse = "/"), " test ", paste(n, collapse = "/"))
  set.seed(job$seed)
  warn <- character()
  sc <- withCallingHandlers(fit_scores(sp, groups = "ovr", lib = lib),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  set.seed(job$seed + 500000L)
  if (job$strata > 0L) {
    # analysis-period week index of each pooled test week, in the order of sc$L (group 1, 2, 3)
    wte <- unlist(lapply(seq_len(k), function(i) which(reg == i & role == "test")))
    stopifnot(length(wte) == length(sc$L), all(reg[wte] == sc$L))
    Lall <- perm_labels_strata(sc$L, job$B, (wte - 1L) %/% job$strata)
  } else Lall <- perm_labels(sc$L, job$B)
  Tl <- lapply(sc$Sovr, function(S) ovr_T(S, seq_len(k), Lall)); Tlib <- do.call(rbind, Tl)
  st <- list(MCA = apply(Tlib, 2L, max))
  for (a in lib) st[[paste0("MCA_", a)]] <- apply(Tl[[a]], 2L, max)
  st$Two_sided <- apply(abs(Tlib), 2L, max)
  st$Sum <- colSums(Tlib)
  p <- lapply(st, pval_perm); p$Comb <- pval_comb(st$MCA, st$Sum)
  if (job$competitors == 1L) {
    Xfit <- do.call(rbind, lapply(sp, `[[`, "fit"))
    Xte <- do.call(rbind, lapply(sp, `[[`, "test"))
    Zte <- apply_std(Xte, standardize_fit(Xfit))
    D2 <- sqdist(Zte); D <- sqrt(D2); nb <- ncol(Lall)
    G <- group_indicator(Lall, k)
    st$DISCO <- disco_stat(D, G, n, k, nb)
    h2 <- median(D2[upper.tri(D2)]); K <- exp(-D2 / (2 * h2))
    st$MMD <- mmd_stat(K, G, n, k, nb)
    E <- mst_edges(D)
    st$MST <- colSums(Lall[E[, 1L], , drop = FALSE] == Lall[E[, 2L], , drop = FALSE])
    for (m in c("DISCO", "MMD", "MST")) p[[m]] <- pval_perm(st[[m]])
    sch <- schott_test(lapply(seq_len(k), function(i) sp[[i]]$test))
    p$Schott <- unname(sch[["p"]])
  }
  stat_obs <- vapply(names(p), function(m) if (m %in% names(st)) st[[m]][1L] else if (m == "Schott") unname(sch[["stat"]]) else NA_real_, numeric(1))
  res <- data.frame(job = job$job, kind = job$kind, dataset = job$dataset, representation = job$representation,
                    label_col = job$label_col, strata = job$strata, method = names(p), p_value = unlist(p),
                    stat_obs = unname(stat_obs), B = job$B,
                    fit = paste(f, collapse = "/"), test = paste(n, collapse = "/"), d = ncol(X),
                    seconds = round(proc.time()[["elapsed"]] - started, 1), row.names = NULL)
  if (job$kind == "main") {
    auc <- do.call(rbind, lapply(lib, function(a) data.frame(learner = a, regime = seq_len(k),
      auc = vapply(seq_len(k), function(i) auc_vec(sc$Sovr[[a]][, i], sc$L == i), numeric(1)),
      T = Tl[[a]][, 1L])))
    write.csv(auc, file.path(out, "ovr_T.csv"), row.names = FALSE)
    info <- do.call(rbind, lapply(lib, function(a) { v <- sc$info_ovr[[a]]
      if (!length(v)) return(NULL)
      data.frame(learner = a, item = names(v), value = vapply(v, function(x) paste(signif(x, 6), collapse = "/"), ""))
    }))
    if (is.null(info)) info <- data.frame(learner = character(), item = character(), value = character())
    write.csv(info, file.path(out, "tuning.csv"), row.names = FALSE)
    saveRDS(list(job = job, statistics = st, L = sc$L, Sovr = sc$Sovr, Lall = Lall, warnings = warn), file.path(out, "detail.rds"))
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  }
  writeLines(unique(warn), file.path(out, "warnings.txt"))
  write.csv(res, file.path(out, "pvalues.csv"), row.names = FALSE)   # written last: marks the job as done
  "done"
}

cat(format(Sys.time()), " jobs:", nrow(jobs), " cores:", cores, "\n"); flush.console()
eng <- function(src) for (fl in c("learners.R", "learners_lib.R", "tests.R")) source(file.path(src, fl))
# progress bar after every finished job; the remaining time is extrapolated
# from jobs run in this session (skipped jobs count as done but not for the rate)
t0 <- proc.time()[["elapsed"]]; n_new <- 0L; n_err <- 0L
progress <- function(k, r, s) {
  if (s != "skipped") n_new <<- n_new + 1L
  if (grepl("^ERROR", s)) { n_err <<- n_err + 1L; cat("  ", jobs$job[r], jobs$out[r], s, "\n") }
  el <- (proc.time()[["elapsed"]] - t0) / 60; n <- nrow(jobs)
  eta <- if (n_new > 0L) el / n_new * (n - k) else NA
  w <- 30L; nb <- round(w * k / n)
  cat(sprintf("  [%s%s] %3.0f%% | %d/%d jobs | errors %d | elapsed %.1f min | remaining ~%s min\n",
              strrep("#", nb), strrep("-", w - nb), 100 * k / n, k, n, n_err, el,
              if (is.na(eta)) "?" else sprintf("%.1f", eta)))
  flush.console()
}
run_safe <- function(j, root, src, smoke)
  tryCatch(run_job(j, root, src, smoke), error = function(e) paste("ERROR:", conditionMessage(e)))
if (cores > 1L && nrow(jobs) > 1L) {
  cl <- parallel::makePSOCKcluster(min(cores, nrow(jobs)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, eng, src)
  parallel::clusterExport(cl, c("run_job", "perm_labels_strata"))
  todo <- split(jobs, seq_len(nrow(jobs)))
  # load-balanced dispatch as in clusterApplyLB, but reporting each job as it finishes
  status <- vector("list", length(todo)); nxt <- 1L
  send <- function(node) { parallel:::sendCall(cl[[node]], run_safe, list(todo[[nxt]], root, src, smoke), tag = nxt); nxt <<- nxt + 1L }
  for (node in seq_along(cl)) if (nxt <= length(todo)) send(node)
  for (k in seq_along(todo)) {
    d <- parallel:::recvOneResult(cl)
    s <- if (inherits(d$value, "try-error")) paste("ERROR:", as.character(d$value)) else d$value
    status[[d$tag]] <- s
    if (nxt <= length(todo)) send(d$node)
    progress(k, d$tag, s)
  }
} else {
  eng(src)
  status <- lapply(seq_len(nrow(jobs)), function(r) {
    s <- run_safe(jobs[r, ], root, src, smoke)
    cat(jobs$job[r], jobs$out[r], s, "\n"); progress(r, r, s); s })
}
status <- unlist(status)
cat("done:", sum(status == "done"), " skipped:", sum(status == "skipped"), " errors:", sum(grepl("^ERROR", status)), "\n")
if (any(grepl("^ERROR", status))) print(data.frame(job = jobs$job, out = jobs$out, status = status)[grepl("^ERROR", status), ], row.names = FALSE)
fs <- file.path(root, jobs$out, "pvalues.csv"); fs <- fs[file.exists(fs)]
if (length(fs)) {
  pl <- lapply(fs, read.csv, stringsAsFactors = FALSE)
  cols <- unique(unlist(lapply(pl, names)))   # files written before amendment 1 lack strata and stat_obs
  allp <- do.call(rbind, lapply(pl, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
  dir.create(file.path(root, if (smoke) "smoke/results" else "results"), recursive = TRUE, showWarnings = FALSE)
  target <- file.path(root, if (smoke) "smoke/results" else "results", paste0("pvalues_", sub("\\.csv$", "", jobs_file), ".csv"))
  write.csv(allp, target, row.names = FALSE)
  cat("collected", nrow(allp), "rows in", target, "\n")
}
