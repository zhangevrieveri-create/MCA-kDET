param([int]$Workers=20,[int]$TargetCpu=70,[int]$BlockSize=5,[int]$TimeoutMin=60,[switch]$Smoke,[switch]$PrepareOnly,[string]$Out='')
$ErrorActionPreference='Stop'
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
if(-not $Out){$Out=if($Smoke){'靶向模拟_20261102/smoke_B19'}else{'靶向模拟_20261102/full_B200'}}
$outPath=[IO.Path]::GetFullPath((Join-Path (Get-Location) $Out))
if($Workers -lt 1 -or $Workers -gt 26 -or $BlockSize -lt 1){throw 'Invalid worker/block configuration'}
$env:LANG='en_US.UTF-8';$env:LC_ALL='English_United States.utf8';$env:R_ENABLE_JIT='0'
$env:R_LIBS_USER='C:/Users/zhang/AppData/Local/R/win-library/4.3'
$env:OMP_NUM_THREADS='1';$env:OPENBLAS_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
$r='C:/Program Files/R/R-4.3.3/bin/Rscript.exe'
if(-not(Test-Path -LiteralPath (Join-Path $outPath 'grid.rds'))) {
  $s=if($Smoke){'true'}else{'false'}
  & $r --vanilla 靶向模拟_20261102/worker.R prepare "--out=$outPath" "--smoke=$s"
  if($LASTEXITCODE -ne 0){throw 'Preparation failed'}
}
if($PrepareOnly){return}
$manifest=@(Import-Csv -LiteralPath (Join-Path $outPath 'source_md5.csv'))
function Check-Sources {
  foreach($h in $manifest) {
    if((Get-FileHash -LiteralPath $h.file -Algorithm MD5).Hash.ToLowerInvariant() -ne $h.md5){throw "Source changed: $($h.file)"}
  }
}
Check-Sources
$lock=$null;$active=[Collections.Generic.List[object]]::new();$transcript=$false
try {
  $lock=[IO.File]::Open((Join-Path $outPath 'run.lock'),[IO.FileMode]::OpenOrCreate,[IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
  # A lost master can leave pinned workers alive. Never overlap their ranges.
  $orphans=@(Get-CimInstance Win32_Process | Where-Object {($_.Name -in @('pwsh.exe','powershell.exe','Rscript.exe')) -and $_.ProcessId -ne $PID -and $_.CommandLine -and $_.CommandLine.Contains($outPath) -and ($_.CommandLine.Contains('worker.R') -or $_.CommandLine.Contains('run_process.ps1'))})
  if($orphans.Count){throw "Existing workers for this output remain active: $($orphans.ProcessId -join ',')"}
  Start-Transcript -LiteralPath (Join-Path $outPath 'parallel_console.txt') -Append|Out-Null;$transcript=$true
  $grid=@(Import-Csv -LiteralPath (Join-Path $outPath 'grid.csv'))
  $total=[int](($grid|ForEach-Object{[int]$_.reps}|Measure-Object -Sum).Sum)
  $rowPath=Join-Path $outPath 'rows';$logPath=Join-Path $outPath 'worker_logs'
  $queue=[Collections.Generic.List[object]]::new()
  $maxReps=[int](($grid|ForEach-Object{[int]$_.reps}|Measure-Object -Maximum).Maximum)
  for($start=1;$start -le $maxReps;$start+=$BlockSize) {
    for($ci=0;$ci -lt $grid.Count;$ci++) {
      $c=$grid[$ci];if($start -gt [int]$c.reps){continue}
      $end=[Math]::Min([int]$c.reps,$start+$BlockSize-1)
      $missing=@($start..$end|Where-Object{-not(Test-Path -LiteralPath (Join-Path $rowPath ('{0}__r{1:D5}.rds' -f $c.cell,$_)))})
      if($missing.Count){$queue.Add([pscustomobject]@{Index=$ci+1;Cell=$c.cell;First=$start;Last=$end;Attempts=0})}
    }
  }
  if(-not ('MotivCpuClock' -as [type])) {
    Add-Type -TypeDefinition 'using System.Runtime.InteropServices; public static class MotivCpuClock { [DllImport("kernel32.dll")] public static extern bool GetSystemTimes(out long idle,out long kernel,out long user); [DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint flags); }'
  }
  function Cpu-Snapshot {
    $idle=0L;$kernel=0L;$user=0L
    if(-not [MotivCpuClock]::GetSystemTimes([ref]$idle,[ref]$kernel,[ref]$user)){throw 'CPU sampling failed'}
    return @($idle,$kernel,$user)
  }
  $clock=[Diagnostics.Stopwatch]::StartNew();$initialCount=@(Get-ChildItem -LiteralPath $rowPath -Filter '*.rds').Count
  $cpuPrevious=Cpu-Snapshot;$lastCpu=Get-Date;$lastState=[datetime]::MinValue;$lastTune=Get-Date;$cpuAverage=0.0
  $history=[Collections.Generic.Queue[double]]::new();$runStamp=Get-Date -Format 'yyyyMMdd_HHmmss'
  $cpuList=@(0..([Environment]::ProcessorCount-1))
  $Workers=[Math]::Min($Workers,$cpuList.Count)
  function Show-State($state,$detail) {
    $done=@(Get-ChildItem -LiteralPath $rowPath -Filter '*.rds').Count
    $eta=if($done-$initialCount -ge 10){$clock.Elapsed.TotalHours/($done-$initialCount)*($total-$done)}else{$null}
    $obj=[pscustomobject]@{State=$state;Updated=(Get-Date -Format 'yyyy-MM-dd HH:mm:ss');PID=$PID;
      Done=$done;Total=$total;B=if($Smoke){19}else{200};ActiveWorkers=$active.Count;WorkerLimit=$Workers;
      CpuBusyPercent=$cpuAverage;TargetCpuPercent=$TargetCpu;ElapsedHours=$clock.Elapsed.TotalHours;
      RemainingHours=$eta;Detail=$detail;Out=$outPath}
    $obj|ConvertTo-Json|Set-Content -LiteralPath (Join-Path $outPath 'progress.json') -Encoding utf8
    Write-Progress -Activity '靶向模拟T1T2' -Status "$done/$total; $($active.Count) workers; CPU~$([math]::Round($cpuAverage))%" -PercentComplete ($done*100/$total)
    Write-Host "[$(Get-Date -Format HH:mm:ss)] $state | $done/$total | active=$($active.Count)/$Workers | CPU~$([math]::Round($cpuAverage,1))% | remaining~$eta hours | $detail"
  }
  function Launch-Job($job,$cpu) {
    Check-Sources;$job.Attempts++
    $psi=[Diagnostics.ProcessStartInfo]::new()
    $psi.FileName=Join-Path $PSHOME 'pwsh.exe';$psi.UseShellExecute=$false;$psi.CreateNoWindow=$true
    $psi.RedirectStandardOutput=$true;$psi.RedirectStandardError=$true;$psi.WorkingDirectory=(Get-Location).Path
    foreach($a in @('-NoProfile','-File',(Join-Path $PSScriptRoot 'run_process.ps1'),'-CpuIndex',[string]$cpu,
      '-Cell',[string]$job.Index,'-First',[string]$job.First,'-Last',[string]$job.Last,'-Out',$outPath)){$psi.ArgumentList.Add($a)}
    if($Smoke -and $job.First -eq 1){$psi.ArgumentList.Add('-RepeatCheck')}
    $proc=[Diagnostics.Process]::Start($psi)
    $active.Add([pscustomobject]@{Job=$job;Cpu=$cpu;Process=$proc;Stdout=$proc.StandardOutput.ReadToEndAsync();Stderr=$proc.StandardError.ReadToEndAsync()})
    Write-Host "START CPU $cpu | $($job.Cell) | $($job.First)-$($job.Last) | process attempt $($job.Attempts)/3"
  }
  [MotivCpuClock]::SetThreadExecutionState([uint32]2147483649)|Out-Null
  Show-State 'running' 'All fixed methods; independent pinned single-thread processes; original rows preserved.'
  while($true) {
    foreach($item in @($active)) {
      if(-not $item.Process.HasExited -and ((Get-Date)-$item.Process.StartTime).TotalMinutes -gt $TimeoutMin) {
        Write-Host "HARD TIMEOUT ${TimeoutMin}min | $($item.Job.Cell) | $($item.Job.First)-$($item.Job.Last)"
        $item.Process.Kill($true)
      }
    }
    foreach($item in @($active|Where-Object{$_.Process.HasExited})) {
      $item.Process.WaitForExit();$job=$item.Job
      $jobLog=Join-Path $logPath ("$($job.Cell)_$($job.First)-$($job.Last)_${runStamp}_try$($job.Attempts).log")
      @($item.Stdout.GetAwaiter().GetResult(),$item.Stderr.GetAwaiter().GetResult())|Set-Content -LiteralPath $jobLog -Encoding utf8
      $active.Remove($item)|Out-Null
      if($item.Process.ExitCode -eq 0){Write-Host "DONE | $($job.Cell) | $($job.First)-$($job.Last)"}
      else {
        if($job.Attempts -ge 3){throw "Three process failures: $($job.Cell) $($job.First)-$($job.Last). See $jobLog"}
        Write-Host "RETRY exit=$($item.Process.ExitCode); same seeds; saved rows retained | $($job.Cell)"
        $queue.Insert(0,$job)
      }
    }
    if((Get-Date)-$lastCpu -ge [timespan]::FromSeconds(15)) {
      $now=Cpu-Snapshot;$dt=($now[1]-$cpuPrevious[1])+($now[2]-$cpuPrevious[2])
      if($dt -gt 0){$history.Enqueue(100*(1-($now[0]-$cpuPrevious[0])/$dt));while($history.Count -gt 3){$history.Dequeue()|Out-Null};$cpuAverage=($history|Measure-Object -Average).Average}
      $cpuPrevious=$now;$lastCpu=Get-Date
      if($history.Count -ge 2 -and (Get-Date)-$lastTune -ge [timespan]::FromSeconds(30)) {
        $old=$Workers
        if($cpuAverage -lt $TargetCpu-7 -and $queue.Count -gt 0 -and $active.Count -ge $Workers-2){$Workers=[Math]::Min([Math]::Min(26,$cpuList.Count),$Workers+2)}
        elseif($cpuAverage -gt $TargetCpu+7){$Workers=[Math]::Max(4,$Workers-2)}
        if($old -ne $Workers){Write-Host "CPU ADJUST $old -> $Workers; observed $([math]::Round($cpuAverage,1))%"}
        $lastTune=Get-Date
      }
      [MotivCpuClock]::SetThreadExecutionState([uint32]2147483649)|Out-Null
    }
    foreach($cpu in $cpuList) {
      if(-not $queue.Count -or $active.Count -ge $Workers){break}
      if(@($active|Where-Object{$_.Cpu -eq $cpu}).Count){continue}
      $job=$queue[0];$queue.RemoveAt(0);Launch-Job $job $cpu
    }
    if((Get-Date)-$lastState -ge [timespan]::FromSeconds(15)){Show-State 'running' 'Fixed grid; no outcome-based modification.';$lastState=Get-Date}
    if(-not $queue.Count -and -not $active.Count){break}
    Start-Sleep -Seconds 2
  }
  Check-Sources
  if(@(Get-ChildItem -LiteralPath $rowPath -Filter '*.rds').Count -ne $total){throw 'Saved repetition count mismatch'}
  Show-State 'reporting' 'All raw replications saved; generating integrity checks and comparison tables.'
  & (Join-Path $PSScriptRoot 'run_process.ps1') -CpuIndex 16 -Action report -Out $outPath
  if($LASTEXITCODE -ne 0){throw 'Raw results complete; report failed. Repair report only; do not recompute simulations.'}
  Show-State 'complete' 'All calculations and report checks complete; no final library or paper changes.'
  Write-Progress -Activity '靶向模拟T1T2' -Completed
} catch {
  if($clock){Show-State 'stopped' $_.Exception.Message}
  throw
} finally {
  foreach($item in @($active)) {
    if(-not $item.Process.HasExited){$item.Process.Kill($true)}
    $item.Process.WaitForExit()
    if($logPath) {
      $job=$item.Job
      $abortLog=Join-Path $logPath ("$($job.Cell)_$($job.First)-$($job.Last)_${runStamp}_try$($job.Attempts)_master_stopped.log")
      @($item.Stdout.GetAwaiter().GetResult(),$item.Stderr.GetAwaiter().GetResult())|Set-Content -LiteralPath $abortLog -Encoding utf8
    }
  }
  if('MotivCpuClock' -as [type]){[MotivCpuClock]::SetThreadExecutionState([uint32]2147483648)|Out-Null}
  if($lock){$lock.Dispose()}
  if($transcript){Stop-Transcript|Out-Null}
}
