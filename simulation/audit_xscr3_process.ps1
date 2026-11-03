param([int]$CpuIndex,[int]$First,[int]$Last)
$ErrorActionPreference='Stop'
Set-Location -LiteralPath $PSScriptRoot
$env:LANG='en_US.UTF-8';$env:LC_ALL='English_United States.utf8'
$env:R_LIBS_USER='C:/Users/zhang/AppData/Local/R/win-library/4.3'
$env:OMP_NUM_THREADS='1';$env:OPENBLAS_NUM_THREADS='1';$env:MKL_NUM_THREADS='1';$env:R_ENABLE_JIT='0'
(Get-Process -Id $PID).ProcessorAffinity=[IntPtr]([long]1 -shl $CpuIndex)
$out='tests/out/xscr3_audit_20261101/rows'
New-Item -ItemType Directory -Path $out -Force|Out-Null
& 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla tests/audit_xscr3_null.R $First $Last $out
exit $LASTEXITCODE
