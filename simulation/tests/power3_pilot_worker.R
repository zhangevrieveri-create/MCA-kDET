# Isolated worker for the complete, pre-specified libpilot3 pilot.
# Scientific functions and seed definitions are sourced without modification.
args <- commandArgs(TRUE)
arg <- function(name, default = NULL) {
  x <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(x)) sub(paste0("^--", name, "="), "", x[1L]) else default
}
SRC_DIR <- file.path(getwd(), "R"); PROFILE <- "quick"
for (p in c("glmnet", "ranger", "gTestsMulti"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Package failed to load: ", p)
for (f in c("models.R", "learners.R", "learners_lib.R", "tests.R", "sim.R", "garch_ext.R",
            "competitors2.R", "rotation.R", "learners_lib2.R", "models_s.R", "experiments.R", "summarize.R"))
  source(file.path(SRC_DIR, f))
out <- arg("out"); action <- arg("action", "init")
stopifnot(!is.null(out))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
grid <- make_grid("libpilot3")
stopifnot(nrow(grid) == 27L, sum(grid$reps) == 5100L, all(grid$B == 200L))
manifest <- file.path(out, "formal_grid.rds")
if (file.exists(manifest)) stopifnot(identical(readRDS(manifest), grid))
emit <- function(...) { cat(paste0(...), "\n", sep = ""); flush.console() }
if (action == "init") {
  if (!file.exists(manifest)) {
    saveRDS(grid, manifest)
    write.csv(grid, file.path(out, "formal_grid.csv"), row.names = FALSE)
    writeLines(c(R.version.string, paste("Libraries:", paste(.libPaths(), collapse = "; ")),
                 vapply(c("glmnet", "ranger", "gTestsMulti"), function(p)
                   paste(p, utils::packageVersion(p)), "")), file.path(out, "R_environment.txt"))
  }
  emit("INIT OK: 27 cells; 5100 repetitions; B=200; original seeds and libraries.")
} else if (action == "chunk") {
  stopifnot(file.exists(manifest))
  id <- arg("cell"); end <- as.integer(arg("end")); start <- as.integer(arg("start"))
  row <- grid[grid$cell == id, ]; stopifnot(nrow(row) == 1L, start >= 1L, end <= row$reps)
  done <- as.integer(arg("done", "0")); total <- sum(grid$reps)
  progress <- function(n, new) {
    if (new) {
      done <<- done + n
      nb <- floor(30 * done / total)
      emit(sprintf("[%s] [%s%s] %.2f%% | %d/%d saved | %s",
                   format(Sys.time(), "%H:%M:%S"), strrep("#", nb), strrep("-", 30 - nb),
                   100 * done / total, done, total, id))
    }
  }
  emit(sprintf("CHUNK %s repetitions %d-%d, B=200", id, start, end))
  run_cell(as.list(row), reps = end, cl = NULL, outdir = out, chunk = 1L, log = emit, progress = progress)
  for (r in start:end) {
    f <- file.path(out, "libpilot3", sprintf("%s__r%05d-%05d.csv", id, r, r))
    z <- read.csv(f, check.names = FALSE)
    if (nrow(z) != 1L || any(!is.na(z$err) & z$err != "") || z$B != 200L || z$seed != row$seed0 + r) {
      writeLines(paste("Invalid saved repetition:", f), file.path(out, "STOP_error.txt"))
      stop("Invalid saved repetition: ", f)
    }
  }
  emit("CHUNK OK")
} else if (action == "check") {
  s <- summarize_exp("libpilot3", out)
  id <- arg("cell"); z <- s[s$cell == id, ]
  stopifnot(nrow(z) == 1L, z$errors == 0L, z$reps == 500L)
  cols <- c("mca_xscr1", "mca_xscr3", "mca_g5f")
  rates <- unlist(z[cols], use.names = TRUE)
  emit(paste("NULL SIZE", id, paste(sprintf("%s=%.4f", names(rates), rates), collapse = ", ")))
  if (any(!is.finite(rates) | rates < 0.031 | rates > 0.069)) {
    writeLines(c(id, paste(names(rates), rates)), file.path(out, "STOP_level.txt"))
    stop("Pre-specified null-size gate failed; stop and report to the user.")
  }
  emit("NULL SIZE GATE OK")
} else if (action == "summary") {
  s <- summarize_exp("libpilot3", out)
  stopifnot(sum(s$errors) == 0L)
  emit(sprintf("SUMMARY: %d cells, %d valid repetitions, %d errors", nrow(s), sum(s$reps), sum(s$errors)))
} else stop("Unknown action: ", action)
