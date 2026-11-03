# =============================================================================
# Frozen after the calibration round E0 (21 Oct 2026; 50 reps per cell,
# k = 3, d = 100, beta = 1, B = 200). Rule: the smallest grid value at which
# the best method has power >= 0.95 at beta = 1.
#   A: cA in {0.2,0.3,0.4,0.5,0.6,0.8}: best = MCA-L1 0.98 at 0.6        -> cA = 0.6
#   B: cB in {1,2,3,4,6}: best = Schott 0.94 at 3, 1.00 at 4              -> cB = 4
#   D: n in {300,500}: MCA 1.00 at n = 300                                -> n = 300
#   E: cE in {0.1,0.2,0.3,0.5}: best = DISCO 0.98 at 0.2                  -> cE = 0.2
#   C', C2: no method reaches 0.95 at d = 100 (see Cdiag); pending decision
# beta_fix (E3, E5, E6): beta at which MCA has power about 0.6-0.8 for k = 3.
# =============================================================================
PAR$n  <- 300L
PAR$cA <- 0.6
PAR$cB <- 4
PAR$cE <- 0.2
PAR$beta_fix <- c(A = 0.8, C = 0.6, D = 0.6)
