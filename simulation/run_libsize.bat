@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist results mkdir results
set RS="D:\R-4.6.1\bin\Rscript.exe"
rem libsize: D null, d=100, 1000 reps (R/experiments.R). Re-running this script
rem resumes: finished chunks are skipped, never overwritten. Keep --chunk=60.
%RS% run.R --exp=libsize --cores=30 --chunk=60 >> results\libsize_console.txt 2>&1
echo LIBSIZE DONE >> results\libsize_console.txt
