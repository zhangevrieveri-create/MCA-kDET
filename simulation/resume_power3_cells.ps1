# One fresh R process per cell; preserve all existing result chunks.
param([string[]]$SkipCells = @())
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
$env:LANG = 'en_US.UTF-8'
$env:LC_ALL = 'English_United States.utf8'
$env:R_LIBS_USER = 'C:/Users/zhang/AppData/Local/R/win-library/4.3'
$resultDir = Join-Path $PSScriptRoot 'results_trial_power3_progress_fixed_20261101'
$gridPath = Join-Path $resultDir 'libpilot3_grid.csv'
$driverGridPath = Join-Path $resultDir 'driver_grid.csv'
if (Test-Path -LiteralPath $driverGridPath) { $gridPath = $driverGridPath }
$cells = @(Import-Csv -LiteralPath $gridPath)
if ($cells.Count -ne 27) { throw 'Expected the full 27-cell grid.' }
if (-not (Test-Path -LiteralPath $driverGridPath)) { Copy-Item -LiteralPath $gridPath -Destination $driverGridPath }
Start-Transcript -LiteralPath (Join-Path $resultDir 'cell_driver_console.txt') -Append | Out-Null
$pendingCells = @()
try {
    for ($cellIndex = 0; $cellIndex -lt $cells.Count; $cellIndex++) {
        $cellId = $cells[$cellIndex].cell
        if ($cellId -in $SkipCells) {
            Write-Host "PENDING cell $($cellIndex+1)/27: $cellId, explicitly deferred after earlier crashes."
            $pendingCells += $cellId
            continue
        }
        $cellFiles = @(Get-ChildItem -LiteralPath (Join-Path $resultDir 'libpilot3') -Filter ($cellId + '__r*.csv'))
        if ($cellFiles.Count -eq 3) {
            Write-Host "CACHE cell $($cellIndex+1)/27: $cellId, 3/3 saved"
            continue
        }
        $cellFilter = '^' + [regex]::Escape($cellId) + '$'
        for ($attemptIndex = 1; $attemptIndex -le 3; $attemptIndex++) {
            $savedCount = @(Get-ChildItem -LiteralPath (Join-Path $resultDir 'libpilot3') -Filter '*.csv').Count
            Write-Host "OVERALL $savedCount/81 saved | cell $($cellIndex+1)/27 | process attempt $attemptIndex/3"
            & 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe' --vanilla run.R --exp=libpilot3 --profile=quick --reps=3 --B=19 --cores=1 --chunk=1 "--filter=$cellFilter" "--out=$resultDir" --print-summary=no
            if ($LASTEXITCODE -eq 0) { break }
            Write-Host "R exit $LASTEXITCODE; completed repetitions retained."
        }
        if ($LASTEXITCODE -ne 0) {
            Write-Host "PENDING cell $cellId after 3 failed process attempts; continuing with the next cell."
            $pendingCells += $cellId
        }
    }
    $savedRows = @(Get-ChildItem -LiteralPath (Join-Path $resultDir 'libpilot3') -Filter '*.csv' | ForEach-Object { Import-Csv -LiteralPath $_.FullName })
    $errorRows = @($savedRows | Where-Object { $_.err -and $_.err -ne '' })
    Write-Host "OVERALL $($savedRows.Count)/81 saved; errors=$($errorRows.Count)."
    $pendingCells | Set-Content -LiteralPath (Join-Path $resultDir 'pending_cells.txt') -Encoding utf8
    if ($savedRows.Count -ne 81 -or $errorRows.Count -ne 0) {
        Write-Host "TIMING TRIAL PARTIAL: $($pendingCells.Count) pending cells; formal pilot not started."
        exit 2
    }
    Write-Host 'TIMING TRIAL COMPLETE. STOP POINT 1; formal pilot not started.'
} finally { Stop-Transcript | Out-Null }
