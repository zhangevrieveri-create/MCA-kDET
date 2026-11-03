@echo off
rem ===========================================================================
rem  QUICK VERSION (快速版) of the revised E1-E7: fewer replications.
rem  Results: results_v2_quick\   (see 运行说明_精确版与快速版.md)
rem  Needs R\calib_params_v2.R filled in (run run_calib2.bat first).
rem  Re-running resumes; finished chunks are skipped.
rem ===========================================================================
chcp 65001 >nul
cd /d "%~dp0"
set RS="D:\R-4.6.1\bin\Rscript.exe"
set OUT=results_v2_quick
if not exist %OUT% mkdir %OUT%
for %%E in (E8 E2 E4 E6 E3 E5 E1 E7) do (
  %RS% run.R --exp=%%E --profile=quick --cores=30 >> %OUT%\%%E_console.txt 2>&1
)
echo V2 QUICK DONE >> %OUT%\DONE.txt
