@echo off
rem ===========================================================================
rem  Third review round (G1, G2), 1 Nov 2026: experiments E9 and E10.
rem    run_review3.bat quick   100 reps per alternative, 200 per null cell
rem    run_review3.bat full    1000 / 2000 (the precise version)
rem  E9  : new competitors on the E2 cells and the k = 3 null cells of E1,
rem        with the SAME data as E1/E2 (same seeds)
rem  E10 : rotated (M3)-(M5), MCA-kDET, accuracy tests and all competitors
rem  Then report_review3.R writes results_v2_<profile>\review3_*.csv/.txt
rem  Before the first run: install gTestsMulti in a separate R session and
rem  run tests\test_review3.R (see the notes of 1 Nov 2026 in sim\).
rem  Re-running resumes; finished chunks are skipped.
rem ===========================================================================
chcp 65001 >nul
cd /d "%~dp0"
rem  R 4.3.3: 4.6.1 / 4.4.3 crash on this machine (libdeflate), 2026-11-01
set RS="C:\Program Files\R\R-4.3.3\bin\Rscript.exe"
set PROF=%1
if "%PROF%"=="" set PROF=quick
set OUT=results_v2_%PROF%
if not exist %OUT% mkdir %OUT%
for %%E in (E9 E10) do (
  %RS% run.R --exp=%%E --profile=%PROF% --cores=30 >> %OUT%\%%E_console.txt 2>&1
)
%RS% report_review3.R --out=%OUT% >> %OUT%\review3_console.txt 2>&1
echo REVIEW3 %PROF% DONE >> %OUT%\DONE_review3.txt
