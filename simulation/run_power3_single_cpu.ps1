# Pin this launcher and its R children to one logical CPU; restore the shell afterwards.
param(
    [string]$Filter = 'S2',
    [string]$Out = 'results_trial_power3_cpu16_jit0_20261101',
    [int]$CpuIndex = 16
)
$ErrorActionPreference = 'Stop'
if ($CpuIndex -lt 0 -or $CpuIndex -ge [Environment]::ProcessorCount) { throw 'Invalid logical CPU index.' }
$rscriptPath = 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe'
$settings = @{
    LANG='en_US.UTF-8'; LC_ALL='English_United States.utf8'
    R_LIBS_USER='C:/Users/zhang/AppData/Local/R/win-library/4.3'
    OMP_NUM_THREADS='1'; OPENBLAS_NUM_THREADS='1'; MKL_NUM_THREADS='1'; R_ENABLE_JIT='0'
}
$previousSettings = @{}
$launcherProcess = Get-Process -Id $PID
$previousAffinity = $launcherProcess.ProcessorAffinity
$previousDirectory = Get-Location
try {
    foreach ($key in $settings.Keys) {
        $previousSettings[$key] = [Environment]::GetEnvironmentVariable($key, 'Process')
        [Environment]::SetEnvironmentVariable($key, $settings[$key], 'Process')
    }
    $launcherProcess.ProcessorAffinity = [IntPtr]([long]1 -shl $CpuIndex)
    Set-Location -LiteralPath $PSScriptRoot
    Write-Host "R 4.3.3 | one logical CPU (index $CpuIndex) | one BLAS thread | JIT off"
    & $rscriptPath --vanilla run.R --exp=libpilot3 --profile=quick --reps=3 --B=19 --cores=1 --chunk=1 "--filter=$Filter" "--out=$Out" --print-summary=no
    $runExitCode = $LASTEXITCODE
} finally {
    $launcherProcess.ProcessorAffinity = $previousAffinity
    foreach ($key in $previousSettings.Keys) { [Environment]::SetEnvironmentVariable($key, $previousSettings[$key], 'Process') }
    Set-Location -LiteralPath $previousDirectory
}
exit $runExitCode
