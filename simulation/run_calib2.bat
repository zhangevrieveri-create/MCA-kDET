@echo off
rem second calibration round (after the change of library), 14 cells x 50 reps.
rem Same in both versions. Results: results\calib2\  -> fill in R\calib_params_v2.R
chcp 65001 >nul
cd /d "%~dp0"
set RS="D:\R-4.6.1\bin\Rscript.exe"
%RS% run.R --exp=calib2 --cores=30 >> results\calib2_console.txt 2>&1
echo CALIB2 DONE >> results\calib2_console.txt
