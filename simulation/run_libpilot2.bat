@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist results mkdir results
set RS="D:\R-4.6.1\bin\Rscript.exe"
rem libpilot2: smaller candidate libraries, 12 alternative cells x 100 reps,
rem fresh seeds (R/experiments.R). Re-running resumes; finished chunks are skipped.
%RS% run.R --exp=libpilot2 --cores=30 >> results\libpilot2_console.txt 2>&1
%RS% report_libpilot.R --exp=libpilot2 >> results\libpilot2_console.txt 2>&1
echo LIBPILOT2 DONE >> results\libpilot2_console.txt
