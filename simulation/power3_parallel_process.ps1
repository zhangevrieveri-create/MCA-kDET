param([int]$CpuIndex,[string]$Action,[string]$Cell,[int]$First=1,[int]$Last=1,[string]$Out)
$ErrorActionPreference='Stop'
Set-Location -LiteralPath $PSScriptRoot
$env:LANG='en_US.UTF-8'; $env:LC_ALL='English_United States.utf8'
$env:R_LIBS_USER='C:/Users/zhang/AppData/Local/R/win-library/4.3'
$env:OMP_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:R_ENABLE_JIT='0'
(Get-Process -Id $PID).ProcessorAffinity=[IntPtr]([long]1 -shl $CpuIndex)
$rPath='C:/Program Files/R/R-4.3.3/bin/Rscript.exe'
if($Action -eq 'chunk') {
    & $rPath --vanilla tests/power3_parallel_chunk.R "--out=$Out" "--cell=$Cell" "--start=$First" "--end=$Last"
} else {
    & $rPath --vanilla tests/power3_pilot_worker.R "--out=$Out" "--action=$Action" "--cell=$Cell"
}
exit $LASTEXITCODE
