@echo off
rem Runs an R script again until a "done" file exists (max 40 tries). R sometimes
rem exits silently, even with exit code 0, on this machine (CLAUDE.md); scripts that
rem cache their results simply resume.
rem   tests\retry_r.bat <donefile> <script.R> [args]
setlocal
set RS="C:\Program Files\R\R-4.3.3\bin\Rscript.exe"
set DONE=%1
shift
set /a n=0
:again
if exist %DONE% exit /b 0
set /a n+=1
if %n% GTR 40 (
  echo giving up after 40 tries
  exit /b 1
)
echo [try %n%] %1 %2 %3 %4 %5 %6
%RS% %1 %2 %3 %4 %5 %6 %7 %8 %9
goto again
