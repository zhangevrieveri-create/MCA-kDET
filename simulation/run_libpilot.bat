@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist results mkdir results
set RS="D:\R-4.6.1\bin\Rscript.exe"
rem learner-library pilot (TASK_学习器库.md, section 5); ranger must be installed beforehand
%RS% run.R --exp=libpilot --cores=30 > results\libpilot_console.txt 2>&1
echo LIBPILOT DONE >> results\libpilot_console.txt
