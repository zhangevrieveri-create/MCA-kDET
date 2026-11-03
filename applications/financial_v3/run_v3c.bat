@echo off
rem Amendment 2 (protocol_amendment2.json): truncated time-shift (TTS) test.
rem Run from sim\applications\financial_v3 after run_v3_all.bat and run_v3b.bat.
set R="D:\R-4.6.1\bin\Rscript.exe"
set PYTHONIOENCODING=utf-8
rem 19 of 32 cores, about 60 percent, so the machine stays usable
set CORES=19
cd /d %~dp0
echo === 1 build jobs (labels for every shift, group-size check) ===
python make_jobs_v3c.py || goto :err
echo === 2 smoke test (2 jobs, B = 19, l1 only, written to smoke\) ===
%R% run_v3.R --jobs=jobs_tts.csv --cores=1 --smoke || goto :err
echo === 3 TTS test: main regimes first, then the lag-2 regimes (2236 jobs, B = 19) ===
%R% run_v3.R --jobs=jobs_tts.csv --cores=%CORES% || goto :err
echo === 4 lag-2 regimes, within-era permutation (2 jobs, B = 4999) ===
%R% run_v3.R --jobs=jobs_tts_strata.csv --cores=2 || goto :err
echo === 5 report ===
python report_tts.py > results\report_tts_console.txt || goto :err
echo Finished. See results\报告_v3_TTS.md
goto :eof
:err
echo Stopped with an error. Re-running this file resumes: finished jobs are skipped.
exit /b 1
