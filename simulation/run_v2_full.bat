@echo off
rem ===========================================================================
rem  FULL VERSION (精确版) of the revised E1-E7: the replications of the
rem  design, for the paper. Results: results_v2_full\
rem  The first replications are bit-identical to the quick version (same
rem  seeds). Needs R\calib_params_v2.R filled in. 2000 reps per null cell,
rem  1000 per alternative cell (22 Oct 2026); about 75 h on 30 cores.
rem  Re-running resumes; finished chunks are skipped.
rem ===========================================================================
chcp 65001 >nul
cd /d "%~dp0"
set RS="D:\R-4.6.1\bin\Rscript.exe"
set OUT=results_v2_full
if not exist %OUT% mkdir %OUT%
for %%E in (E8 E2 E4 E6 E3 E5 E1 E7) do (
  %RS% run.R --exp=%%E --profile=full --cores=30 >> %OUT%\%%E_console.txt 2>&1
)
echo V2 FULL DONE >> %OUT%\DONE.txt
