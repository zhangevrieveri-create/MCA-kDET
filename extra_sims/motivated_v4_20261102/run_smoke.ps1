$ErrorActionPreference='Stop'
$env:LANG='en_US.UTF-8';$env:LC_ALL='English_United States.utf8';$env:R_ENABLE_JIT='0'
$env:R_LIBS_USER='C:/Users/zhang/AppData/Local/R/win-library/4.3'
$env:OMP_NUM_THREADS='1';$env:OPENBLAS_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
$r='C:/Program Files/R/R-4.3.3/bin/Rscript.exe'
$script='motivated_v4_20261102/smoke.R'
if(-not(Test-Path 'motivated_v4_20261102/smoke_B19/grid.rds')) {
  & $r --vanilla $script prepare
  if($LASTEXITCODE -ne 0){throw 'Preparation failed'}
}
$jobs=@();$attempts=@{}
for($i=1;$i -le 20;$i++) {
  $runStamp=Get-Date -Format 'yyyyMMdd_HHmmss'
  $log="motivated_v4_20261102/smoke_B19/logs/cell_${i}_$runStamp"
  $p=Start-Process -FilePath $r -ArgumentList @('--vanilla',$script,'cell',$i) -PassThru -WindowStyle Hidden -RedirectStandardOutput "$log.out" -RedirectStandardError "$log.err"
  $p.ProcessorAffinity=[IntPtr]([long]1 -shl ($i+5));$jobs+=@{id=$i;process=$p};$attempts[$i]=1
}
while($jobs.Count -gt 0) {
  Start-Sleep -Seconds 3
  $running=@()
  foreach($j in $jobs) {
    $j.process.Refresh()
    if(-not $j.process.HasExited){$running+=$j;continue}
    if($j.process.ExitCode -ne 0) {
      if($attempts[$j.id] -ge 3){throw "Cell $($j.id) failed after 3 process attempts; see logs"}
      $attempts[$j.id]++
      $log="motivated_v4_20261102/smoke_B19/logs/cell_$($j.id)_${runStamp}_attempt$($attempts[$j.id])"
      $p=Start-Process -FilePath $r -ArgumentList @('--vanilla',$script,'cell',$j.id) -PassThru -WindowStyle Hidden -RedirectStandardOutput "$log.out" -RedirectStandardError "$log.err"
      $p.ProcessorAffinity=[IntPtr]([long]1 -shl ($j.id+5));$running+=@{id=$j.id;process=$p}
    }
  }
  $jobs=$running
  $done=@(Get-ChildItem -LiteralPath 'motivated_v4_20261102/smoke_B19/rows' -Filter '*.rds').Count
  Write-Progress -Activity '新模型机械试运行（不作功效结论）' -Status "$done / 60 saved; $($jobs.Count) workers" -PercentComplete ($done*100/60)
  Write-Host "$(Get-Date -Format HH:mm:ss) $done/60 saved; $($jobs.Count) workers"
}
& $r --vanilla $script report
if($LASTEXITCODE -ne 0){throw 'Report failed'}
Write-Progress -Activity '新模型机械试运行（不作功效结论）' -Completed
