# Complete pilot only. Formal E1-E13 runs require the user's library decision.
param([int]$CpuIndex = 16, [int]$BlockSize = 5, [string]$Out = 'results_v3_pilot_20261101_cpu16')
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
if ($BlockSize -lt 1) { throw 'BlockSize must be positive.' }
$rscriptPath = 'C:/Program Files/R/R-4.3.3/bin/Rscript.exe'
$env:LANG='en_US.UTF-8'; $env:LC_ALL='English_United States.utf8'
$env:R_LIBS_USER='C:/Users/zhang/AppData/Local/R/win-library/4.3'
$env:OMP_NUM_THREADS='1'; $env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:R_ENABLE_JIT='0'
(Get-Process -Id $PID).ProcessorAffinity=[IntPtr]([long]1 -shl $CpuIndex)
$outPath = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot $Out))
New-Item -ItemType Directory -Path $outPath -Force | Out-Null
$resultPath = Join-Path $outPath 'libpilot3'
New-Item -ItemType Directory -Path $resultPath -Force | Out-Null
$statePath = Join-Path $outPath 'progress.json'
$timer = [Diagnostics.Stopwatch]::StartNew()
$initialCount = @(Get-ChildItem -LiteralPath $resultPath -Filter '*.csv').Count
function Save-State($state, $cellId, $done, $detail) {
    $newCount = $done - $initialCount
    $etaHours = if ($newCount -gt 0) { $timer.Elapsed.TotalHours / $newCount * (5100 - $done) } else { $null }
    [pscustomobject]@{ State=$state; Updated=(Get-Date -Format 'yyyy-MM-dd HH:mm:ss'); PID=$PID;
        CpuIndex=$CpuIndex; Workers=1; JIT=0; B=200; Done=$done; Total=5100; Cell=$cellId;
        ElapsedHours=$timer.Elapsed.TotalHours; RemainingHours=$etaHours; Detail=$detail } |
        ConvertTo-Json | Set-Content -LiteralPath $statePath -Encoding utf8
}
Start-Transcript -LiteralPath (Join-Path $outPath 'console.txt') -Append | Out-Null
try {
    Save-State 'initializing' '' $initialCount 'Complete pilot authorized by user; null cells first.'
    $sourcePaths = @(Get-ChildItem -LiteralPath (Join-Path $PSScriptRoot 'R') -Filter '*.R' -File | Select-Object -ExpandProperty FullName)
    $sourcePaths += Join-Path $PSScriptRoot 'tests/power3_pilot_worker.R'
    $hashPath = Join-Path $outPath 'source_sha256.csv'
    if (Test-Path -LiteralPath $hashPath) { $sourceHashes = @(Import-Csv -LiteralPath $hashPath) }
    else {
        $sourceHashes = @($sourcePaths | ForEach-Object { Get-FileHash -LiteralPath $_ -Algorithm SHA256 } | Select-Object Path,Hash)
        $sourceHashes | Export-Csv -LiteralPath $hashPath -NoTypeInformation -Encoding utf8
    }
    function Assert-Source {
        foreach ($h in $sourceHashes) {
            if ((Get-FileHash -LiteralPath $h.Path -Algorithm SHA256).Hash -ne $h.Hash) { throw "Source changed during pilot: $($h.Path)" }
        }
    }
    Assert-Source
    & $rscriptPath --vanilla tests/power3_pilot_worker.R --action=init "--out=$outPath"
    if ($LASTEXITCODE -ne 0) { throw 'Pilot initialization failed.' }
    $grid = @(Import-Csv -LiteralPath (Join-Path $outPath 'formal_grid.csv') | Sort-Object @{Expression={if([double]$_.beta -eq 0){0}else{1}}})
    foreach ($cell in $grid) {
        for ($start = 1; $start -le [int]$cell.reps; $start += $BlockSize) {
            $end = [Math]::Min([int]$cell.reps, $start + $BlockSize - 1)
            $expected = @($start..$end | ForEach-Object { Join-Path $resultPath ('{0}__r{1:D5}-{1:D5}.csv' -f $cell.cell,$_) })
            if (@($expected | Where-Object {-not (Test-Path -LiteralPath $_)}).Count -eq 0) { continue }
            for ($attempt = 1; $attempt -le 3; $attempt++) {
                Assert-Source
                $done = @(Get-ChildItem -LiteralPath $resultPath -Filter '*.csv').Count
                Save-State 'running' $cell.cell $done "Repetitions $start-$end; process attempt $attempt/3."
                Write-Host "[$(Get-Date -Format 'HH:mm:ss')] OVERALL $done/5100 | $($cell.cell) | $start-$end | attempt $attempt/3"
                & $rscriptPath --vanilla tests/power3_pilot_worker.R --action=chunk "--out=$outPath" "--cell=$($cell.cell)" "--start=$start" "--end=$end" "--done=$done"
                if ($LASTEXITCODE -eq 0) { break }
                if (Test-Path -LiteralPath (Join-Path $outPath 'STOP_error.txt')) { throw 'A saved repetition contains an error. Stop for review.' }
                Write-Host "R exit $LASTEXITCODE; same-seed resume on the next attempt."
            }
            if ($LASTEXITCODE -ne 0) { throw "Three process attempts failed at $($cell.cell), repetitions $start-$end." }
        }
        if ([double]$cell.beta -eq 0) {
            & $rscriptPath --vanilla tests/power3_pilot_worker.R --action=check "--out=$outPath" "--cell=$($cell.cell)"
            if ($LASTEXITCODE -ne 0) { throw "Null-size check failed or crashed for $($cell.cell). See console and STOP_level.txt." }
        }
        & $rscriptPath --vanilla tests/power3_pilot_worker.R --action=summary "--out=$outPath"
        if ($LASTEXITCODE -ne 0) { throw 'Summary failed; raw results are retained.' }
    }
    $done = @(Get-ChildItem -LiteralPath $resultPath -Filter '*.csv').Count
    if ($done -ne 5100) { throw "Expected 5100 saved repetitions; got $done." }
    Save-State 'complete' '' $done 'Pilot complete. STOP POINT 2: user selects the library; E1-E13 not started.'
    Write-Host 'PILOT COMPLETE. STOP POINT 2: library decision belongs to the user.'
} catch {
    $done = @(Get-ChildItem -LiteralPath $resultPath -Filter '*.csv').Count
    Save-State 'stopped' '' $done $_.Exception.Message
    Write-Host "STOPPED: $($_.Exception.Message)"
    throw
} finally { Stop-Transcript | Out-Null }
