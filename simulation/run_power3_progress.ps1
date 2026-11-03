# Visible, finite validation + three-replication timing run. No formal pilot.
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
$env:LANG = 'en_US.UTF-8'
$env:LC_ALL = 'English_United States.utf8'
$env:R_LIBS_USER = 'C:/Users/zhang/AppData/Local/R/win-library/4.3'
$rscriptPath = 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe'
$runStamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$runDir = Join-Path $PSScriptRoot ('results_trial_power3_progress_' + $runStamp)
New-Item -ItemType Directory -Path $runDir -Force | Out-Null
$transcriptPath = Join-Path $runDir 'console.txt'
Start-Transcript -LiteralPath $transcriptPath | Out-Null
try {
    Write-Host "Output: $runDir"
    Write-Host '[1/4] Calibration: reuse existing file, or generate with iteration progress.'
    if (-not (Test-Path -LiteralPath 'R/calib_params_s.R')) {
        & $rscriptPath --vanilla tests/make_calib_s.R
        if ($LASTEXITCODE -ne 0) { throw "Calibration failed (exit $LASTEXITCODE)." }
    } else { Write-Host 'Existing calibration found; validation follows.' }
    Write-Host '[2/4] Full validation (cached replications are explicitly counted).'
    for ($validationAttempt = 1; $validationAttempt -le 3; $validationAttempt++) {
        Write-Host "Validation process attempt $validationAttempt/3"
        & $rscriptPath --vanilla tests/test_power3.R
        if ($LASTEXITCODE -eq 0) { break }
        Write-Host "Validation process exited with code $LASTEXITCODE. No timing run started."
    }
    if ($LASTEXITCODE -ne 0) { throw "Validation failed after 3 attempts; timing trial was not started." }
    Write-Host '[3/4] Seed-range checks.'
    & $rscriptPath --vanilla tests/check_seeds.R
    if ($LASTEXITCODE -ne 0) { throw "Seed checks failed (exit $LASTEXITCODE)." }
    Write-Host '[4/4] Timing trial ONLY: 27 cells x 3 reps, B=19, 1 worker; save and update progress after each replication.'
    & $rscriptPath --vanilla run.R --exp=libpilot3 --profile=quick --reps=3 --B=19 --cores=1 --chunk=1 "--out=$runDir"
    if ($LASTEXITCODE -ne 0) { throw "Timing trial failed (exit $LASTEXITCODE). Completed chunks remain saved." }
    Write-Host "Finished. Summary: $(Join-Path $runDir 'libpilot3_summary.csv')"
    Write-Host 'STOP POINT 1: formal pilot has not been started.'
} finally { Stop-Transcript | Out-Null }
