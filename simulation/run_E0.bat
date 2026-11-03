@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist results mkdir results
"C:\Program Files\R\R-4.6.1\bin\Rscript.exe" run.R --exp=E0 --cores=28 --chunk=25 > results\E0_console.txt 2>&1
echo E0 DONE >> results\E0_console.txt
