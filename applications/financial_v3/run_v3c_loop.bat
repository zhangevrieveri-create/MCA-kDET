@echo off
rem Amendment 2 (TTS) with R 4.3.3: on this machine R 4.4.3 and 4.6.1 fail to load packages
rem (libdeflate errors) in most starts, R 4.3.3 loads them every time. Same steps as run_v3c.bat.
rem Re-runs until every TTS and lag-2 within-era job has its pvalues.csv; finished jobs are skipped.
cd /d %~dp0
set R="C:\Program Files\R\R-4.3.3\bin\Rscript.exe"
set PYTHONIOENCODING=utf-8
set CORES=19
python make_jobs_v3c.py || exit /b 1
for /L %%i in (1,1,10) do (
  echo ##### attempt %%i #####
  %R% run_v3.R --jobs=jobs_tts.csv --cores=%CORES%
  %R% run_v3.R --jobs=jobs_tts_strata.csv --cores=2
  python -c "import pandas as pd,os,sys; j=pd.concat([pd.read_csv(f) for f in ['jobs_tts.csv','jobs_tts_strata.csv']]); m=[o for o in j.out if not os.path.exists(os.path.join(o,'pvalues.csv'))]; print('missing jobs:',len(m)); sys.exit(1 if m else 0)" && goto :done
)
echo Still incomplete after 10 attempts.
exit /b 1
:done
python report_tts.py > results\report_tts_console.txt
python profile_tts.py > results\tts_profile_console.txt
echo All jobs complete. See the TTS report and results\tts_profile.csv
