param([int]$Workers=15,[int]$BlockSize=5,[string]$Out='results_v3_pilot_20261101_cpu16')
$ErrorActionPreference='Stop'
Set-Location -LiteralPath $PSScriptRoot
$outPath=Join-Path $PSScriptRoot $Out
$resultPath=Join-Path $outPath 'libpilot3'
$logPath=Join-Path $outPath 'parallel_console.txt'
$workerLogPath=Join-Path $outPath 'worker_logs'
New-Item -ItemType Directory -Path $workerLogPath -Force | Out-Null
$grid=@(Import-Csv -LiteralPath (Join-Path $outPath 'formal_grid.csv'))
$nulls=@($grid|Where-Object {[double]$_.beta -eq 0})
$cpus=@(0..($Workers-1)|ForEach-Object{2*$_})
if($Workers -lt 1 -or $Workers -gt 15 -or $BlockSize -lt 1){throw 'Invalid worker/block configuration.'}
$hashes=@(Import-Csv -LiteralPath (Join-Path $outPath 'source_sha256.csv'))
$extraHashPath=Join-Path $outPath 'parallel_source_sha256.csv'
if(Test-Path -LiteralPath $extraHashPath){$extraHashes=@(Import-Csv -LiteralPath $extraHashPath)}
else{
    $extraHashes=@('tests/power3_parallel_chunk.R','power3_parallel_process.ps1','run_power3_parallel15.ps1'|ForEach-Object{
        Get-FileHash -LiteralPath (Join-Path $PSScriptRoot $_) -Algorithm SHA256
    }|Select-Object Path,Hash)
    $extraHashes|Export-Csv -LiteralPath $extraHashPath -NoTypeInformation -Encoding utf8
}
function Check-Sources {
    foreach($h in @($hashes)+@($extraHashes)){
        if((Get-FileHash -LiteralPath $h.Path -Algorithm SHA256).Hash -ne $h.Hash){throw "Source changed: $($h.Path)"}
    }
}
$queue=[Collections.Generic.List[object]]::new()
function Queue-Phase($phaseRows){
    $maxReps=($phaseRows|Measure-Object reps -Maximum).Maximum
    for($s=1;$s -le $maxReps;$s+=$BlockSize){
        foreach($c in $phaseRows){
            if($s -gt [int]$c.reps){continue}
            $e=[Math]::Min([int]$c.reps,$s+$BlockSize-1)
            $missing=@($s..$e|Where-Object{-not(Test-Path -LiteralPath (Join-Path $resultPath ('{0}__r{1:D5}-{1:D5}.csv' -f $c.cell,$_)))})
            if($missing.Count){$queue.Add([pscustomobject]@{Action='chunk';Cell=$c.cell;First=$s;Last=$e;Attempts=0})}
        }
    }
}
$active=[Collections.Generic.List[object]]::new()
$checked=@{}
$checking=@{}
$phase='null'
$timer=[Diagnostics.Stopwatch]::StartNew()
$initialCount=@(Get-ChildItem -LiteralPath $resultPath -Filter '*.csv').Count
$lastLog=[datetime]::MinValue
function Show-State($state,$detail){
    $done=@(Get-ChildItem -LiteralPath $resultPath -Filter '*.csv').Count
    $eta=if($done-$initialCount -ge 15){$timer.Elapsed.TotalHours/($done-$initialCount)*(5100-$done)}else{$null}
    [pscustomobject]@{State=$state;Updated=(Get-Date -Format 'yyyy-MM-dd HH:mm:ss');PID=$PID;Workers=$Workers;
        ActiveWorkers=$active.Count;CpuIndices=$cpus;B=200;JIT=0;Done=$done;Total=5100;Phase=$phase;
        ElapsedHours=$timer.Elapsed.TotalHours;RemainingHours=$eta;Detail=$detail}|
        ConvertTo-Json|Set-Content -LiteralPath (Join-Path $outPath 'progress.json') -Encoding utf8
    Write-Host "[$(Get-Date -Format 'HH:mm:ss')] $state | $done/5100 | active=$($active.Count)/$Workers | $phase | remaining~$eta hours"
}
function Launch-Job($job,$cpu){
    Check-Sources
    $job.Attempts++
    $psi=[Diagnostics.ProcessStartInfo]::new()
    $psi.FileName=Join-Path $PSHOME 'pwsh.exe';$psi.UseShellExecute=$false;$psi.CreateNoWindow=$true
    $psi.RedirectStandardOutput=$true;$psi.RedirectStandardError=$true;$psi.WorkingDirectory=$PSScriptRoot
    foreach($a in @('-NoProfile','-File',(Join-Path $PSScriptRoot 'power3_parallel_process.ps1'),'-CpuIndex',[string]$cpu,
        '-Action',$job.Action,'-Cell',$job.Cell,'-First',[string]$job.First,'-Last',[string]$job.Last,'-Out',$outPath)){$psi.ArgumentList.Add($a)}
    $proc=[Diagnostics.Process]::Start($psi)
    $active.Add([pscustomobject]@{Job=$job;Cpu=$cpu;Process=$proc;Stdout=$proc.StandardOutput.ReadToEndAsync();Stderr=$proc.StandardError.ReadToEndAsync()})
    Write-Host "START CPU $cpu | $($job.Action) | $($job.Cell) | $($job.First)-$($job.Last) | attempt $($job.Attempts)/3"
}
Start-Transcript -LiteralPath $logPath -Append | Out-Null
try{
    Check-Sources
    Queue-Phase $nulls
    Show-State 'running' '15 independent, pinned, single-thread processes; disjoint repetition ranges; preserved earlier rows.'
    while($true){
        if((Test-Path -LiteralPath (Join-Path $outPath 'STOP_level.txt')) -or (Test-Path -LiteralPath (Join-Path $outPath 'STOP_error.txt'))){throw 'Pre-specified size gate or calculation error requires stopping.'}
        foreach($item in @($active|Where-Object{$_.Process.HasExited})){
            $item.Process.WaitForExit()
            $job=$item.Job
            $jobLog=Join-Path $workerLogPath ("$($job.Cell)_$($job.Action)_$($job.First)-$($job.Last)_try$($job.Attempts).log")
            @($item.Stdout.GetAwaiter().GetResult(),$item.Stderr.GetAwaiter().GetResult())|Set-Content -LiteralPath $jobLog -Encoding utf8
            $active.Remove($item)|Out-Null
            if($item.Process.ExitCode -eq 0){
                Write-Host "DONE CPU $($item.Cpu) | $($job.Cell) | $($job.Action) $($job.First)-$($job.Last)"
                if($job.Action -eq 'check'){$checked[$job.Cell]=$true;$checking.Remove($job.Cell)}
            }else{
                if($job.Attempts -ge 3){throw "Three process failures: $($job.Cell), $($job.First)-$($job.Last). See $jobLog"}
                Write-Host "RETRY exit=$($item.Process.ExitCode) | $($job.Cell), saved rows retained."
                $queue.Insert(0,$job)
            }
        }
        if($phase -eq 'null'){
            foreach($c in $nulls){
                if($checked.ContainsKey($c.cell) -or $checking.ContainsKey($c.cell)){continue}
                $queued=@($queue|Where-Object{$_.Cell -eq $c.cell}).Count
                $running=@($active|Where-Object{$_.Job.Cell -eq $c.cell}).Count
                if($queued -eq 0 -and $running -eq 0){
                    $count=@(Get-ChildItem -LiteralPath $resultPath -Filter ($c.cell+'__r*.csv')).Count
                    if($count -ne 500){throw "Null cell has $count saved rows, expected 500: $($c.cell)"}
                    $queue.Insert(0,[pscustomobject]@{Action='check';Cell=$c.cell;First=1;Last=500;Attempts=0})
                    $checking[$c.cell]=$true
                }
            }
            if($checked.Count -eq $nulls.Count -and $active.Count -eq 0 -and $queue.Count -eq 0){
                $phase='power';Queue-Phase @($grid|Where-Object{[double]$_.beta -ne 0})
                Write-Host 'All null-size gates passed. Starting every alternative cell.'
            }
        }
        foreach($cpu in $cpus){
            if($queue.Count -eq 0){break}
            if(@($active|Where-Object{$_.Cpu -eq $cpu}).Count){continue}
            $job=$queue[0];$queue.RemoveAt(0);Launch-Job $job $cpu
        }
        if((Get-Date)-$lastLog -gt [timespan]::FromSeconds(15)){Show-State 'running' 'Disjoint 5-repetition blocks; same seeds and statistical definitions.';$lastLog=Get-Date}
        if($phase -eq 'power' -and $active.Count -eq 0 -and $queue.Count -eq 0){break}
        Start-Sleep -Seconds 2
    }
    & (Join-Path $PSScriptRoot 'power3_parallel_process.ps1') -CpuIndex 16 -Action summary -Cell '' -Out $outPath
    if($LASTEXITCODE -ne 0){throw 'Final summary failed.'}
    if(@(Get-ChildItem -LiteralPath $resultPath -Filter '*.csv').Count -ne 5100){throw 'Expected 5100 saved repetitions.'}
    Show-State 'complete' 'Complete pilot. STOP POINT 2: user chooses library; E1-E13 not started.'
}catch{
    Show-State 'stopped' $_.Exception.Message
    throw
}finally{
    foreach($item in $active){if(-not $item.Process.HasExited){$item.Process.Kill($true)}}
    Stop-Transcript|Out-Null
}
