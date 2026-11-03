# =============================================================================
# Experiment grids (Chapter 4 design, sections 2 and 4)
#
# PAR holds the quantities fixed by the calibration round (E0). Edit them
# only through calib_params.R, written after E0, then freeze.
# =============================================================================

PAR <- list(
  n = 300L, d = 100L, B = 200L,
  cA = 0.5, cB = 3, sC = 5L, cE = 0.3,     # placeholders until E0 is finished
  beta_fix = c(A = 0.6, C = 0.6, D = 0.6), # fixed beta for E3, E5, E6
  reps_null = 1000L, reps_alt = 500L, reps_calib = 50L,
  seed = 20261021L
)
if (file.exists(file.path(SRC_DIR, "calib_params.R"))) source(file.path(SRC_DIR, "calib_params.R"))

ALL_GROUPS <- "ovr,shared,pair,dist,schott"
BETAS <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)

# polynomial string hash in [0, 199998] (199999 is prime); exact in doubles
seed_hash <- function(codes) Reduce(function(h, c) (h * 131 + c) %% 199999, codes, 0)

cell_row <- function(exp, model, beta, k = 3L, n = PAR$n, d = PAR$d, tau = 0.5,
                     nsplit = 1L, groups = ALL_GROUPS, B = PAR$B,
                     cA = PAR$cA, cB = PAR$cB, sC = PAR$sC, cE = PAR$cE,
                     reps = if (beta == 0) PAR$reps_null else PAR$reps_alt, tag = "",
                     nn = "", seed_key = NULL, lib = "l1,nn", libsets = "",
                     seed_scheme = "v1") {
  nlab <- gsub("[,;: ]+", "-", as.character(n))
  id <- sprintf("%s_%s_b%.2f_k%d_n%s_d%d_tau%.1f_S%d%s", exp, model, beta, k, nlab, d,
                tau, nsplit, tag)
  # the seed depends on the cell only, so every cell is reproducible on its own
  #   v1 (all grids up to the pilots): seeds spaced by 1000, weak hash; cells
  #      with more than ~1000 reps can share seeds with other cells
  #   v2 (revised E1-E7, calib2): polynomial hash mod 199999,
  #      seeds 2e8 + 9000 h, i.e. above every v1 seed and 9000 apart; the
  #      grids are checked for duplicated hashes (make_grid)
  sk <- utf8ToInt(if (is.null(seed_key)) id else seed_key)
  seed0 <- switch(seed_scheme,
    v1 = PAR$seed + (sum(sk * seq_along(sk)) %% 100000L) * 1000L,
    v2 = 200000000 + seed_hash(sk) * 9000,
    stop("unknown seed_scheme: ", seed_scheme))
  data.frame(exp = exp, cell = id, model = model, beta = beta, k = k,
             n = as.character(n), d = d, tau = tau, nsplit = nsplit, groups = groups,
             B = B, cA = cA, cB = cB, sC = sC, cE = cE, nn = nn, reps = reps, seed0 = seed0,
             lib = lib, libsets = libsets, stringsAsFactors = FALSE)
}

# Default library chosen after the pilots (libpilot, libpilot2, libweak;
# decided 21 Oct 2026, see 试点报告_学习器库.md, section 10). cell_row()
# keeps lib = "l1,nn" as its default on purpose: the cell id does not contain
# the library, so changing the default would let a resumed E1-E7 run mix
# chunks of the old and the new library in the same results directory.
LIB_MAIN <- "l1,quad,knn,rf"

# learner-library pilot (TASK_学习器库.md, section 5)
LIBPILOT_LIB <- "l1,quad,knn,rf,nn"
LIBPILOT_SETS <- "l1+nn;l1+quad+knn+nn;l1+quad+knn+rf+nn"
lp_cell <- function(model, beta, d, reps, sC = PAR$sC, tag = "", exp = "libpilot",
                    libsets = LIBPILOT_SETS)
  cell_row(exp, model, beta, k = 3L, n = PAR$n, d = d, tau = 0.5,
           groups = "ovr,dist,schott", B = 200L, sC = sC, reps = reps, tag = tag,
           lib = LIBPILOT_LIB, libsets = libsets)

# second pilot: smaller candidate libraries, the 12 alternative cells of
# libpilot with fresh seeds (the cell id, hence the seed, starts with libpilot2)
LIBPILOT2_SETS <- "l1+nn;l1+rf;l1+quad+rf;l1+quad+knn+rf;l1+quad+rf+nn;l1+quad+knn+rf+nn"
lp2_cell <- function(model, beta, d, sC = PAR$sC, tag = "")
  lp_cell(model, beta, d, 100L, sC = sC, tag = tag, exp = "libpilot2", libsets = LIBPILOT2_SETS)

# =============================================================================
# Revised E1-E7 (E1-E7修订设计, 21 Oct 2026)
#
# TWO VERSIONS that differ ONLY in the number of replications:
#   --profile=full   precise version of the design (for the paper)
#   --profile=quick  fewer replications (used for the current results)
# Cells, cell ids and therefore seeds (scheme v2, see cell_row) are the same
# in both versions, so the
# replications 1..R of a quick run are bit-identical to the first R
# replications of the full run. run.R writes them to results_v2_full/ and
# results_v2_quick/ respectively. The profile is required for E1-E7.
# =============================================================================
if (!exists("PROFILE")) PROFILE <- NA_character_
V2_EXPS <- c(paste0("E", 1:10), "calib2")  # grids with seed scheme v2 (E9, E10: 1 Nov 2026)
# full version raised on 22 Oct 2026 to 2000 reps for every null cell and
# 1000 for every alternative cell (E7 measures time only and keeps 20).
# Seeds are spaced by 9000 per cell (scheme v2), so up to 9000 reps per cell
# never overlap. Earlier full values: E1 1000/500, E2-E3 250, E4 1000/500,
# E5 500/500, E6 500.
REPS_V2 <- list(                                   #   full   quick
  E1_k3   = c(full = 2000L, quick = 200L),         # E1, k = 3 (and unequal n)
  E1_k510 = c(full = 2000L, quick = 200L),         # E1, k = 5 and 10
  E2      = c(full = 1000L, quick = 100L),
  E3      = c(full = 1000L, quick = 100L),
  E4_null = c(full = 2000L, quick = 200L),
  E4_alt  = c(full = 1000L, quick = 100L),
  E5_null = c(full = 2000L, quick = 200L),
  E5_alt  = c(full = 1000L, quick = 100L),
  E6      = c(full = 1000L, quick = 100L),
  E7      = c(full =   20L, quick =   5L),
  E8      = c(full = 2000L, quick = 200L))         # GARCH null (size only)
nrep <- function(key) {
  if (length(PROFILE) != 1L || is.na(PROFILE) || !PROFILE %in% c("full", "quick"))
    stop("E1-E7 need --profile=full or --profile=quick")
  REPS_V2[[key]][[PROFILE]]
}

# second calibration (after the change of library): R/calib_params_v2.R
#   beta_fix_C   : fixed beta of C' (E3, E5, E6)
#   beta_grid_C2 : the 5 betas of C2 in E2
#   beta_grid_D  : the 5 betas of D in E2
PAR2 <- list(beta_fix_C = NA_real_, beta_grid_C2 = rep(NA_real_, 5L),
             beta_grid_D = rep(NA_real_, 5L))
if (file.exists(file.path(SRC_DIR, "calib_params_v2.R"))) source(file.path(SRC_DIR, "calib_params_v2.R"))
# grid change of model D after calib2 (s = 1 is a null case), 22 Oct 2026
if (file.exists(file.path(SRC_DIR, "calib_params_v3.R"))) source(file.path(SRC_DIR, "calib_params_v3.R"))
need_calib <- function(x, what, len = 1L) {
  if (is.null(x) || length(x) != len || anyNA(x))
    stop(what, " is not calibrated yet: run the experiment calib2 and fill in R/calib_params_v2.R")
  x
}

# E8: GARCH(1,1) parameters (alpha, beta) by persistence label; the other
# keys of GARCH_DEFAULT (R/garch_ext.R) keep their defaults
E8_PERSIST <- c(p95 = "alpha=0.10;beta=0.85", p99 = "alpha=0.05;beta=0.94")

SC_V2 <- 10L                          # dependent pairs of C' and C2
D_C_V2 <- 20L                         # dimension of C' and C2 in E2, E3, E5
GRID_FRAC <- c(0.2, 0.4, 0.6, 0.8, 1.0)
v2_dims_null <- function(m) if (m == "C") c(D_C_V2, 100L) else c(100L, 300L)
# a cell of the revised design (seed scheme v2); C'/C2 alternatives carry the
# number of pairs in the id
v2 <- function(exp, model, beta, ..., sC = SC_V2)
  cell_row(exp, model, beta, ..., sC = sC, lib = LIB_MAIN, seed_scheme = "v2",
           tag = if (model %in% c("C", "C2") && beta > 0) sprintf("_sC%d", sC) else "")
# second calibration round, 50 reps
c2_cell <- function(model, beta, d)
  cell_row("calib2", model, beta, d = d, groups = "ovr,dist,schott", reps = 50L, lib = LIB_MAIN,
           sC = SC_V2, seed_scheme = "v2", tag = if (model %in% c("C", "C2")) sprintf("_sC%d", SC_V2) else "")

# E9, E10 (third review round, G1, 1 Nov 2026): cells PAIRED with E1/E2.
# A paired cell keeps seed0 and reps of its source cell, so replication r
# uses exactly the data of replication r of the source cell (one_rep draws
# the data first); only the exp prefix of the id and the method groups
# change. E9 adds the competitors of R/competitors2.R to the E2 cells and to
# the k = 3, equal-n null cells of E1; E10 rotates (M3)-(M5) (R/rotation.R)
# and runs MCA-kDET, the accuracy tests and all competitors on them.
REVIEW3_COMP <- "comp2,songchen"
paired_grid <- function(src, exp, groups) {
  src$cell <- sub("^E[0-9]+_", paste0(exp, "_"), src$cell)
  src$exp <- exp
  src$groups <- groups
  src
}


# =============================================================================
# Third review round, power and scenarios (TASK_第三轮_功效与场景_20261101.md):
# libpilot3 (library pilot, fresh seeds, results_v3_pilot/), and E11-E13 (new
# scenarios S1-S3, results_v3_<profile>/). Nothing here changes an earlier grid.
# =============================================================================
# the library used by E11-E13 and, after the freeze, by E1-E10 of version 3;
# NA until the user has chosen it from the pilot report (stop point 2)
LIB_MAIN_V3 <- NA_character_
lib_v3 <- function() {
  if (is.na(LIB_MAIN_V3)) stop("LIB_MAIN_V3 is not set: the library is chosen after the pilot libpilot3")
  LIB_MAIN_V3
}
REPS_V3 <- list(S_alt = c(full = 1000L, quick = 100L), S_null = c(full = 2000L, quick = 200L))
nrep3 <- function(key) {
  if (length(PROFILE) != 1L || is.na(PROFILE) || !PROFILE %in% c("full", "quick"))
    stop("E11-E13 need --profile=full or --profile=quick")
  REPS_V3[[key]][[PROFILE]]
}
V3_EXPS <- c("libpilot3", "E11", "E12", "E13")

# pilot: every single learner / component of the task, and the 8 candidate
# libraries L0-L7 as library sets (columns mca_set_<a+b+...>):
#   L0 l1+quad+knn+rf                      L4 l1th+scrm+scr2+quad+knn+rf+g5f
#   L1 l1th+scrm+scr2+knn+rf               L5 xscr1+knn+rf
#   L2 L1+g5t                              L6 xscr1+xscr3+knn+rf+g5f
#   L3 L1+g5f                              L7 xscr1+l1th+quad+knn+rf
LIBPILOT3_LIB <- "l1,l1th,l1se,quad,scrm,scr2,knn,rf,xscr1,xscr3,g5t,g5f"
LIBPILOT3_SETS <- paste("l1+quad+knn+rf", "l1th+scrm+scr2+knn+rf", "l1th+scrm+scr2+knn+rf+g5t",
                        "l1th+scrm+scr2+knn+rf+g5f", "l1th+scrm+scr2+quad+knn+rf+g5f",
                        "xscr1+knn+rf", "xscr1+xscr3+knn+rf+g5f", "xscr1+l1th+quad+knn+rf", sep = ";")
LP3_GROUPS <- "ovr,dist,schott,comp2,songchen"
lp3_cell <- function(model, beta, d, reps, k = 3L, n = PAR$n, rot = 0L, groups = LP3_GROUPS) {
  tag <- paste0(if (model %in% c("C", "C2") && beta > 0) sprintf("_sC%d", SC_V2) else "",
                if (rot > 0L) sprintf("_rot%d", rot) else "")
  g <- cell_row("libpilot3", model, beta, k = k, n = n, d = d, groups = groups, reps = reps,
                sC = SC_V2, lib = LIBPILOT3_LIB, libsets = LIBPILOT3_SETS,
                seed_scheme = "v2", tag = tag)
  g$rot <- as.integer(rot)
  g
}
# scenario cell of E11-E13: beta = ideal AUC - 0.5
s_cell <- function(exp, model, beta, d, k, n, reps, groups) {
  g <- cell_row(exp, model, beta, k = k, n = n, d = d, groups = groups, reps = reps,
                lib = lib_v3(), seed_scheme = "v2")
  g
}
S_GROUPS <- "ovr,dist,schott,comp2,comp2ovr,songchen"

make_grid <- function(exp) {
  g <- switch(exp,
    # ---- smoke test: tiny problems, every method group, every model
    smoke = do.call(rbind, c(
      lapply(c("A", "B", "C", "D", "E"), function(m)
        cell_row("smoke", m, 0, n = 40, d = 10, B = 19, reps = 2, cA = 1, cB = 3, sC = 2, cE = 0.5)),
      lapply(c("A", "B", "C", "D", "E"), function(m)
        cell_row("smoke", m, 1, n = 40, d = 10, B = 19, reps = 2, cA = 1, cB = 3, sC = 2, cE = 0.5)),
      list(cell_row("smoke", "A", 0, n = "20,40,80", d = 10, B = 19, reps = 2),
           cell_row("smoke", "C", 1, k = 5, n = 30, d = 10, B = 19, reps = 2, sC = 2),
           cell_row("smoke", "A", 0.5, n = 40, d = 10, nsplit = 3, groups = "ovr", B = 19, reps = 2)))),

    # ---- quick signal check: simple -> complex dependence, moderate size
    quick = do.call(rbind, lapply(c("A", "B", "D", "C"), function(m) rbind(
      cell_row("quick", m, 0, n = 200, d = 20, B = 99, reps = 10,
               groups = "ovr,shared,dist,schott", cA = 1, cB = 3, sC = 2),
      cell_row("quick", m, 1, n = 200, d = 20, B = 99, reps = 10,
               groups = "ovr,shared,dist,schott", cA = 1, cB = 3, sC = 2)))),

    # ---- network diagnostics: why do B and C' fail at n = 200?
    nndiag = do.call(rbind, unlist(unlist(lapply(c("B", "C"), function(m)
      lapply(c(200L, 1000L), function(nn)
        lapply(c(default = "", pat30 = "patience=30", lr3 = "lr=0.003",
                 nospec = "spectral=FALSE", pat30lr3 = "patience=30;lr=0.003",
                 noES = "patience=1000"), function(v)
          cell_row("nndiag", m, 1, n = nn, d = 20, B = 99, reps = 10, groups = "ovr",
                   cB = 3, sC = 2, nn = v, seed_key = paste0("nndiag", m, nn),
                   tag = paste0("_", if (nzchar(v)) gsub("[=;.]", "", v) else "default"))))),
      recursive = FALSE), recursive = FALSE)),

    # ---- C' difficulty: sanity check vs sklearn (d = 4) and number of pairs
    nndiag2 = rbind(
      cell_row("nndiag2", "C", 1, n = 200, d = 4, sC = 2, B = 99, reps = 10, groups = "ovr", tag = "_s2"),
      do.call(rbind, unlist(lapply(c(300L, 1000L), function(nn) c(
        lapply(c(5L, 10L), function(s)
          cell_row("nndiag2", "C", 1, n = nn, d = 20, sC = s, B = 99, reps = 10, groups = "ovr", tag = sprintf("_s%d", s))),
        lapply(c(10L, 25L, 50L), function(s)
          cell_row("nndiag2", "C", 1, n = nn, d = 100, sC = s, B = 99, reps = 10, groups = "ovr", tag = sprintf("_s%d", s))))),
        recursive = FALSE))),

    # ---- where do the equal-covariance models become learnable?
    Cdiag = do.call(rbind, unlist(unlist(unlist(lapply(c("C", "C2"), function(m)
      lapply(c(6L, 10L, 20L), function(dd)
        lapply(c(1L, dd %/% 2L), function(s)
          lapply(c(300L, 1000L), function(nn)
            cell_row("Cdiag", m, 1, n = nn, d = dd, sC = s, groups = "ovr,dist,schott",
                     reps = 25L, tag = sprintf("_sC%d", s)))))),
      recursive = FALSE), recursive = FALSE), recursive = FALSE)),

    # ---- E0 calibration (50 reps per cell, beta = 1, k = 3, d = 100)
    E0 = rbind(
      do.call(rbind, lapply(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8), function(c)
        cell_row("E0", "A", 1, cA = c, groups = "ovr,dist,schott", reps = PAR$reps_calib, tag = sprintf("_cA%.2f", c)))),
      do.call(rbind, lapply(c(1, 2, 3, 4, 6), function(c)
        cell_row("E0", "B", 1, cB = c, groups = "ovr,dist,schott", reps = PAR$reps_calib, tag = sprintf("_cB%.1f", c)))),
      do.call(rbind, lapply(c(5L, 25L, 50L), function(s)
        cell_row("E0", "C", 1, sC = s, groups = "ovr,dist,schott", reps = PAR$reps_calib, tag = sprintf("_sC%d", s)))),
      cell_row("E0", "C", 1, n = 1000L, sC = 50L, groups = "ovr,dist,schott", reps = PAR$reps_calib, tag = "_sC50"),
      do.call(rbind, lapply(c(1L, 2L, 5L, 10L), function(s)
        cell_row("E0", "C2", 1, sC = s, groups = "ovr,dist,schott", reps = PAR$reps_calib, tag = sprintf("_sC%d", s)))),
      do.call(rbind, lapply(c(300L, 500L), function(nn)
        cell_row("E0", "D", 1, n = nn, groups = "ovr,dist,schott", reps = PAR$reps_calib))),
      do.call(rbind, lapply(c(0.1, 0.2, 0.3, 0.5), function(c)
        cell_row("E0", "E", 1, cE = c, groups = "ovr,dist,schott", reps = PAR$reps_calib, tag = sprintf("_cE%.2f", c))))),

    # ---- E1-E7: revised design (E1-E7修订设计, 21 Oct 2026). The grids of
    #      the original design are kept in archive/code_v1_20261021/.
    # ---- E1 size (Theorem 1)
    E1 = rbind(
      do.call(rbind, unlist(lapply(c("A", "C", "D"), function(m)
        unlist(lapply(c(3L, 5L, 10L), function(k)
          lapply(v2_dims_null(m), function(d)
            v2("E1", m, 0, k = k, d = d,
               reps = nrep(if (k == 3L) "E1_k3" else "E1_k510"),
               groups = if (k >= 10L) "ovr,shared,dist,schott" else ALL_GROUPS))),
          recursive = FALSE)),
        recursive = FALSE)),
      do.call(rbind, unlist(lapply(c("A", "C", "D"), function(m)
        lapply(v2_dims_null(m), function(d)
          v2("E1", m, 0, n = sprintf("%d,%d,%d", PAR$n %/% 2L, PAR$n, 2L * PAR$n), d = d,
             reps = nrep("E1_k3")))),
        recursive = FALSE))),

    # ---- E2 power, main figure: A, B (d = 100), C', C2 (d = 20), D (d = 100)
    E2 = rbind(
      do.call(rbind, lapply(GRID_FRAC, function(f) v2("E2", "A", f, reps = nrep("E2")))),
      do.call(rbind, lapply(GRID_FRAC, function(f) v2("E2", "B", f, reps = nrep("E2")))),
      do.call(rbind, lapply(GRID_FRAC, function(f)
        v2("E2", "C", f, d = D_C_V2, reps = nrep("E2")))),
      do.call(rbind, lapply(need_calib(PAR2$beta_grid_C2, "beta_grid_C2", 5L), function(b)
        v2("E2", "C2", b, d = D_C_V2, reps = nrep("E2")))),
      do.call(rbind, lapply(need_calib(PAR2$beta_grid_D, "beta_grid_D", 5L), function(b)
        v2("E2", "D", b, reps = nrep("E2"))))),

    # ---- E3 effect of k (Propositions 1-2), reduced version
    E3 = rbind(
      do.call(rbind, lapply(c(3L, 5L, 10L), function(k)
        v2("E3", "A", PAR$beta_fix[["A"]], k = k, groups = "ovr,shared,pair", reps = nrep("E3")))),
      do.call(rbind, lapply(c(3L, 5L, 10L), function(k)
        v2("E3", "C", need_calib(PAR2$beta_fix_C, "beta_fix_C"), k = k, d = D_C_V2,
           groups = "ovr,shared,pair", reps = nrep("E3"))))),

    # ---- E4 dense model E: max vs sum vs combined, and DISCO/MMD (the library
    #      ablation for A-D is read off E2)
    E4 = do.call(rbind, lapply(BETAS, function(b)
      v2("E4", "E", b, groups = "ovr,shared,dist",
         reps = nrep(if (b == 0) "E4_null" else "E4_alt")))),

    # ---- E5 split fraction and multiple splits (single split = first split)
    E5 = rbind(
      do.call(rbind, unlist(lapply(c(0, PAR$beta_fix[["A"]]), function(b)
        lapply(c(0.3, 0.5, 0.7), function(t)
          v2("E5", "A", b, tau = t, nsplit = 5L, groups = "ovr",
             reps = nrep(if (b == 0) "E5_null" else "E5_alt")))), recursive = FALSE)),
      do.call(rbind, unlist(lapply(c(0, need_calib(PAR2$beta_fix_C, "beta_fix_C")), function(b)
        lapply(c(0.3, 0.5, 0.7), function(t)
          v2("E5", "C", b, d = D_C_V2, tau = t, nsplit = 5L, groups = "ovr",
             reps = nrep(if (b == 0) "E5_null" else "E5_alt")))), recursive = FALSE))),

    # ---- E6 dimension: A in d = 50..1000; C' with a fixed number of dependent
    #      pairs and more noise coordinates, d = 10..100
    E6 = rbind(
      do.call(rbind, lapply(c(50L, 100L, 300L, 1000L), function(d)
        v2("E6", "A", PAR$beta_fix[["A"]], d = d, groups = "ovr,shared,dist,schott",
           reps = nrep("E6")))),
      do.call(rbind, lapply(c(10L, 20L, 50L, 100L), function(d)
        v2("E6", "C", need_calib(PAR2$beta_fix_C, "beta_fix_C"), d = d,
           sC = if (d == 10L) 5L else SC_V2, groups = "ovr,shared,dist,schott",
           reps = nrep("E6"))))),

    # ---- second calibration round (after the change of library), 50 reps
    calib2 = rbind(
      do.call(rbind, lapply(c(0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2), function(b)
        c2_cell("D", b, 100L))),
      do.call(rbind, lapply(c(0.2, 0.3, 0.4, 0.5), function(b) c2_cell("C2", b, D_C_V2))),
      do.call(rbind, lapply(c(0.7, 0.8, 0.9), function(b) c2_cell("C", b, D_C_V2)))),

    # ---- learner-library pilot: 12 alternative cells (100 reps), 3 null cells (200 reps)
    libpilot = rbind(
      lp_cell("A", 0.6, 100L, 100L), lp_cell("A", 1.0, 100L, 100L),
      lp_cell("B", 0.6, 100L, 100L), lp_cell("B", 1.0, 100L, 100L),
      lp_cell("D", 0.4, 100L, 100L), lp_cell("D", 1.0, 100L, 100L),
      lp_cell("C", 0.6, 20L, 100L, sC = 10L, tag = "_sC10"),
      lp_cell("C", 1.0, 20L, 100L, sC = 10L, tag = "_sC10"),
      lp_cell("C2", 0.6, 20L, 100L, sC = 10L, tag = "_sC10"),
      lp_cell("C2", 1.0, 20L, 100L, sC = 10L, tag = "_sC10"),
      lp_cell("C2", 1.0, 100L, 100L, sC = 10L, tag = "_sC10"),
      lp_cell("E", 1.0, 100L, 100L),
      lp_cell("A", 0, 100L, 200L),            # AR(0.5) Gaussian null (also B, E)
      lp_cell("C", 0, 20L, 200L),             # N(0, I) null
      lp_cell("D", 0, 100L, 200L)),           # Gaussian-mixture null

    # ---- size check of the D null (libpilot gave 0.095 for the full library
    #      in 200 reps): 1000 fresh reps, same configuration, own seed
    libsize = lp_cell("D", 0, 100L, 1000L, exp = "libsize"),

    # ---- second pilot: candidate smaller libraries, alternatives only
    libpilot2 = rbind(
      lp2_cell("A", 0.6, 100L), lp2_cell("A", 1.0, 100L),
      lp2_cell("B", 0.6, 100L), lp2_cell("B", 1.0, 100L),
      lp2_cell("D", 0.4, 100L), lp2_cell("D", 1.0, 100L),
      lp2_cell("C", 0.6, 20L, sC = 10L, tag = "_sC10"),
      lp2_cell("C", 1.0, 20L, sC = 10L, tag = "_sC10"),
      lp2_cell("C2", 0.6, 20L, sC = 10L, tag = "_sC10"),
      lp2_cell("C2", 1.0, 20L, sC = 10L, tag = "_sC10"),
      lp2_cell("C2", 1.0, 100L, sC = 10L, tag = "_sC10"),
      lp2_cell("E", 1.0, 100L)),

    # ---- weak covariance alternatives: does the library need quad when rf is
    #      present? (libpilot2 libraries, fresh seeds)
    libweak = rbind(
      lp_cell("D", 0.2, 100L, 100L, exp = "libweak", libsets = LIBPILOT2_SETS),
      lp_cell("B", 0.4, 100L, 100L, exp = "libweak", libsets = LIBPILOT2_SETS)),

    # ---- E7 timing: model A, all methods, every k
    E7 = do.call(rbind, lapply(c(3L, 5L, 10L), function(k)
      v2("E7", "A", PAR$beta_fix[["A"]], k = k, reps = nrep("E7")))),

    # ---- E8 size under a GARCH null (R/garch_ext.R): k = 3 regimes from an
    #      external Markov indicator, T = 900 weekly returns (about 300 per
    #      regime), d = 10, 30; moderate (0.10, 0.85) and high (0.05, 0.94)
    #      volatility persistence; t(5) innovations, conditional correlation
    #      0.3. The column mpar carries the model parameters (only in E8).
    E8 = {
      g8 <- do.call(rbind, unlist(lapply(c(10L, 30L), function(d)
        lapply(names(E8_PERSIST), function(p)
          cell_row("E8", "GARCH", 0, k = 3L, n = 300L, d = d, groups = "ovr,dist",
                   reps = nrep("E8"), lib = LIB_MAIN, seed_scheme = "v2",
                   tag = paste0("_", p)))), recursive = FALSE))
      g8$mpar <- E8_PERSIST[sub("^.*_", "", g8$cell)]
      g8
    },
    E9 = {
      e1 <- make_grid("E1")
      paired_grid(rbind(make_grid("E2"),
                        e1[e1$k == 3L & e1$n == as.character(PAR$n), , drop = FALSE]),
                  "E9", REVIEW3_COMP)
    },
    E10 = {
      e2 <- make_grid("E2")
      g10 <- paired_grid(e2[e2$model %in% c("C", "C2", "D"), , drop = FALSE], "E10",
                         paste0("ovr,shared,dist,schott,", REVIEW3_COMP))
      g10$rot <- 1L
      g10
    },

    libpilot3 = {
      alt <- 100L; nul <- 500L
      do.call(rbind, c(
        list(lp3_cell("A", 0.6, 100L, alt), lp3_cell("A", 1.0, 100L, alt),
             lp3_cell("B", 0.2, 100L, alt), lp3_cell("B", 0.6, 100L, alt),
             lp3_cell("C", 0.6, 20L, alt), lp3_cell("C", 1.0, 20L, alt),
             lp3_cell("C2", 0.36, 20L, alt), lp3_cell("C2", 0.60, 20L, alt),
             lp3_cell("D", 0.02, 100L, alt), lp3_cell("D", 0.04, 100L, alt),
             lp3_cell("E", 0.8, 100L, alt),
             lp3_cell("C", 1.0, 20L, alt, rot = 1L), lp3_cell("D", 0.04, 100L, alt, rot = 1L),
             lp3_cell("S1a", 0.07, 100L, alt), lp3_cell("S1a", 0.14, 100L, alt),
             lp3_cell("S1b", 0.07, 100L, alt), lp3_cell("S1b", 0.14, 100L, alt),
             lp3_cell("S2", 0.07, 100L, alt, k = 4L), lp3_cell("S2", 0.14, 100L, alt, k = 4L),
             lp3_cell("S3a", 0.10, 100L, alt, k = 10L, n = 100L, groups = paste0(LP3_GROUPS, ",comp2ovr")),
             lp3_cell("S3b", 0.10, 20L, alt, k = 10L, n = 100L, groups = paste0(LP3_GROUPS, ",comp2ovr")),
             # nulls
             lp3_cell("A", 0, 100L, nul), lp3_cell("C", 0, 20L, nul), lp3_cell("D", 0, 100L, nul),
             lp3_cell("S1a", 0, 100L, nul),
             lp3_cell("S3a", 0, 100L, nul, k = 10L, n = 100L, groups = paste0(LP3_GROUPS, ",comp2ovr")),
             lp3_cell("S3b", 0, 20L, nul, k = 10L, n = 100L, groups = paste0(LP3_GROUPS, ",comp2ovr")))))
    },
    E11 = do.call(rbind, unlist(lapply(c("S1a", "S1b"), function(m)
      lapply(c(0, S_LEVELS), function(b)
        s_cell("E11", m, b, 100L, 3L, PAR$n, if (b == 0) nrep3("S_null") else nrep3("S_alt"), S_GROUPS))),
      recursive = FALSE)),
    E12 = do.call(rbind, lapply(c(0, S_LEVELS), function(b)
      s_cell("E12", "S2", b, 100L, 4L, PAR$n, if (b == 0) nrep3("S_null") else nrep3("S_alt"), S_GROUPS))),
    E13 = do.call(rbind, unlist(lapply(c(5L, 10L, 20L), function(k)
      lapply(list(c("S3a", 100L), c("S3b", 20L)), function(ms)
        lapply(c(0, 0.10, 0.14), function(b)
          s_cell("E13", ms[1L], b, as.integer(ms[2L]), k, 100L,
                 if (b == 0) nrep3("S_null") else nrep3("S_alt"), S_GROUPS)))), recursive = FALSE)),
    stop("unknown experiment: ", exp))
  rownames(g) <- NULL
  if (anyDuplicated(g$cell)) stop("duplicated cell ids in grid ", exp)
  # revised grids only (seed scheme v2); some v1 grids are known to contain
  # duplicated seeds (e.g. libpilot) and must stay resumable as they are
  if ((exp %in% V2_EXPS || exp %in% V3_EXPS) && anyDuplicated(g$seed0)) stop("duplicated seeds in grid ", exp)
  g
}
