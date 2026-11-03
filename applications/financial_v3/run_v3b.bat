@echo off
rem Amendment 1 (protocol_amendment1.json): dependence-robust inference.
rem Run from sim\applications\financial_v3 after run_v3_all.bat has finished.
set R="D:\R-4.6.1\bin\Rscript.exe"
set PYTHONIOENCODING=utf-8
rem 19 of 32 cores, about 60 percent, so the machine stays usable
set CORES=19
cd /d %~dp0
echo === 1 build jobs ===
python make_jobs_v3b.py || goto :err
echo === 2 smoke test of the within-era permutation (2 jobs, B = 19, l1 only) ===
%R% run_v3.R --jobs=jobs_strata_main.csv --cores=1 --smoke || goto :err
echo === 3 rotation test, primary (2586 jobs, B = 19) ===
%R% run_v3.R --jobs=jobs_rotation.csv --cores=%CORES% || goto :err
echo === 4 within-era permutation: placebo check (1200 jobs, B = 199) ===
%R% run_v3.R --jobs=jobs_strata_calib.csv --cores=%CORES% || goto :err
echo === 5 within-era permutation: real labels (B = 4999) ===
%R% run_v3.R --jobs=jobs_strata_main.csv --cores=%CORES% || goto :err
echo === 6 report ===
python report_robust.py > results\report_robust_console.txt || goto :err
echo Finished. See results\报告_v3_依赖稳健.md
goto :eof
:err
echo Stopped with an error. Re-running this file resumes: finished jobs are skipped.
exit /b 1
