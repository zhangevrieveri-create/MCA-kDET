# =============================================================================
# Second calibration -- FROZEN 21 Oct 2026 (after the change of library to
# {l1, quad, knn, rf}; E1-E7修订设计, section 3). calib_params.R (first
# calibration, E0) stays frozen and unchanged; this file only adds the
# quantities below. Evidence: results/calib2/ (14 cells x 50 reps, k = 3,
# n = 300, B = 200, quad cap 1000) and, for C2, libpilot/libpilot2.
#
# Rule (as in E0): the best method (single learners and competitors included)
# has power >= 0.95 at the largest beta of a grid; beta_fix gives MCA-kDET
# power of about 0.6-0.8 at k = 3.
#
# beta_fix_C = 0.8   C' (d = 20, sC = 10), for E3, E5, E6.
#   MCA power 0.44 / 0.64 / 0.94 at beta = 0.7 / 0.8 / 0.9 (calib2).
#
# beta_grid_C2 = 0.12, 0.24, 0.36, 0.48, 0.60   C2 (d = 20, sC = 10), E2.
#   calib2: best method <= 0.56 for beta in {0.2, ..., 0.5} (rf 0.56 at 0.5).
#   At beta = 0.6 (libpilot, libpilot2; 100 reps each): rf 0.97 / 0.97,
#   library 0.95 / 0.94, so 0.6 is the smallest value meeting the rule;
#   the grid is 0.2, ..., 1.0 times 0.6.
#
# beta_grid_D = 0.01, 0.02, 0.04, 0.08, 0.16   D (d = 100), E2.
#   calib2: MCA 0.40 at 0.02 and 1.00 from 0.04 on (quad 0.66 / 1.00);
#   MST 0.08 / 0.38 / 0.64 / 0.88 / 0.98 / 1.00 at 0.04 / 0.06 / 0.08 / 0.10
#   / 0.15 / 0.20; DISCO, MMD <= 0.10, Schott 0. The rule gives 0.04, but
#   the deviating block has s = ceiling(100 beta) coordinates, so 0.2-1.0
#   times 0.04 would repeat s = 4. The grid doubles s instead (s = 1, 2, 4,
#   8, 16; chosen by the user, option (b)); the best method has power 1 at
#   its largest value, as the rule requires.
# =============================================================================
PAR2$beta_fix_C   <- 0.8
PAR2$beta_grid_C2 <- c(0.12, 0.24, 0.36, 0.48, 0.60)
PAR2$beta_grid_D  <- c(0.01, 0.02, 0.04, 0.08, 0.16)
