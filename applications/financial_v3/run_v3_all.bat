@echo off
rem Chapter 5 financial application v3. Run from this folder: sim\applications\financial_v3
rem Step 0 (once): put 30_Industry_Portfolios_Daily_CSV.zip and 10_Industry_Portfolios_Daily_CSV.zip in raw\
set R="D:\R-4.6.1\bin\Rscript.exe"
rem the report prints Unicode symbols; the default GBK console encoding cannot encode them
set PYTHONIOENCODING=utf-8
cd /d %~dp0
echo === 1 prepare data and jobs ===
python prepare_v3.py || goto :err
python make_jobs_v3.py || goto :err
echo === 2 smoke test (2 jobs, B = 19, l1 only) ===
%R% run_v3.R --jobs=jobs_main.csv --cores=1 --smoke || goto :err
echo === 3 calibration checks (800 jobs, B = 199) ===
%R% run_v3.R --jobs=jobs_calibration.csv --cores=30 || goto :err
echo === 4 main analyses (B = 4999) ===
%R% run_v3.R --jobs=jobs_main.csv --cores=30 || goto :err
echo === 5 report ===
python report_v3.py > results\report_console.txt || goto :err
echo Finished. See results\报告_v3.md
goto :eof
:err
echo Stopped with an error. Re-running this file resumes: finished jobs are skipped.
exit /b 1
