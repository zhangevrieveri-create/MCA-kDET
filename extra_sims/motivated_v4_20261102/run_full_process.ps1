param([int]$CpuIndex,[string]$Action='cell',[int]$Cell,[int]$First=1,[int]$Last=1,[string]$Out,[switch]$RepeatCheck)
$ErrorActionPreference='Stop'
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
$env:LANG='en_US.UTF-8';$env:LC_ALL='English_United States.utf8';$env:R_ENABLE_JIT='0'
$env:R_LIBS_USER='C:/Users/zhang/AppData/Local/R/win-library/4.3'
$env:OMP_NUM_THREADS='1';$env:OPENBLAS_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
(Get-Process -Id $PID).ProcessorAffinity=[IntPtr]([long]1 -shl $CpuIndex)
$r='C:/Program Files/R/R-4.3.3/bin/Rscript.exe'
if($Action -eq 'report') {
  & $r --vanilla motivated_v4_20261102/report_full.R "--out=$Out"
} else {
  $repeat=if($RepeatCheck){'true'}else{'false'}
  & $r --vanilla motivated_v4_20261102/full_worker.R cell "--out=$Out" "--cell=$Cell" "--start=$First" "--end=$Last" "--repeat_check=$repeat"
}
exit $LASTEXITCODE
