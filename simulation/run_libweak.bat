@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist results mkdir results
set RS="D:\R-4.6.1\bin\Rscript.exe"
rem libweak: weak covariance alternatives (D beta=0.2, B beta=0.4), 100 reps each,
rem libpilot2 libraries, fresh seeds (R/experiments.R). Re-running resumes.
%RS% run.R --exp=libweak --cores=30 >> results\libweak_console.txt 2>&1
%RS% report_libpilot.R --exp=libweak >> results\libweak_console.txt 2>&1
echo LIBWEAK DONE >> results\libweak_console.txt
