@echo off
rem ===========================================================================
rem  Copies the finished chunks of the quick version into the full version
rem  (E1-E6 only). Same cells, same seeds, same chunk boundaries (50 reps),
rem  so these files are bit-identical to the first 200 / 100 replications of
rem  the full run, which then skips them. E7 is NOT copied (different chunks).
rem  Nothing is deleted or overwritten in results_v2_quick. The cell
rem  E2_D_b0.01 (s = 1, a null case) is not copied: it left the design on
rem  22 Oct 2026 (R\calib_params_v3.R).
rem  1 Nov 2026: E2 is no longer copied. classifier_tests() now also writes the
rem  accuracy columns macc_ovr, macc_<learner> (third review round, G2), which
rem  the quick chunks of E2 do not have; the full run of E2 therefore computes
rem  all its replications (about 30 extra minutes). The other columns of these
rem  replications are unchanged (tests\test_review3.R, check 10).
rem ===========================================================================
chcp 65001 >nul
cd /d "%~dp0"
for %%E in (E1 E3 E4 E5 E6) do (
  robocopy results_v2_quick\%%E results_v2_full\%%E *.csv /XF E2_D_b0.01_* /XC /XN /XO /NJH /NJS /NP
)
echo copied quick chunks of E1, E3-E6 into results_v2_full
