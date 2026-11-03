# =============================================================================
# Grid change after the second calibration -- 22 Oct 2026, decided by the user
# before the full run. calib_params.R and calib_params_v2.R stay frozen and
# unchanged; this file is sourced after them and overrides one quantity.
#
# beta_grid_D: 0.01, 0.02, 0.04, 0.08, 0.16  ->  0.02, 0.03, 0.04, 0.08, 0.16
#   In model D the deviating block has s = ceiling(100 beta) coordinates with
#   AR(0.9) correlations. For s = 1 the block is a single N(0, 1) coordinate,
#   so P_beta = P_0 and beta = 0.01 is a null case (quick run: every method
#   0.00-0.07, MCA 0.03). It is replaced by s = 3, so that the grid is
#   s = 2, 3, 4, 8, 16, all alternatives, denser where the power rises
#   (quick run: MCA 0.48 at s = 2 and 0.92 at s = 4).
#   The cells s = 2, 4, 8, 16 keep their ids and seeds (scheme v2), so their
#   results are unchanged; only the cell E2_D_b0.03 is new. The results of
#   E2_D_b0.01 stay on disk in results_v2_quick/E2 and are no longer part of
#   the design (copy_quick_to_full.bat does not copy them).
# =============================================================================
PAR2$beta_grid_D <- c(0.02, 0.03, 0.04, 0.08, 0.16)
stopifnot(identical(as.integer(ceiling(100 * PAR2$beta_grid_D)), c(2L, 3L, 4L, 8L, 16L)))
