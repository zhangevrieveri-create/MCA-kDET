@echo off
rem Trial of libpilot3 (3 reps, B = 19) until the summary file exists; R crashes now and
rem then on this machine (CLAUDE.md) and run.R resumes from the finished chunks.
setlocal
cd /d "%~dp0.."
set RS="C:\Program Files\R\R-4.3.3\bin\Rscript.exe"
set OUT=results_trial_power3b
set /a n=0
:again
if exist %OUT%\libpilot3_summary.csv exit /b 0
set /a n+=1
if %n% GTR 40 (echo giving up & exit /b 1)
echo [try %n%]
%RS% run.R --exp=libpilot3 --profile=quick --reps=3 --B=19 --cores=4 --out=%OUT%
goto again
