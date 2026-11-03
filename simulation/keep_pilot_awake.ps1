# Temporary sleep inhibition for this pilot only. The display may turn off.
$ErrorActionPreference='Stop'
Set-Location -LiteralPath $PSScriptRoot
Add-Type -TypeDefinition 'using System.Runtime.InteropServices; public static class PilotSleepRequest { [DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint flags); }'
$out=Join-Path $PSScriptRoot 'results_v3_pilot_20261101_cpu16'
$status=Join-Path $out 'sleep_guard.json'
$deadline=(Get-Date).Date.AddDays(1).AddHours(8)
$stoppedSince=$null
try{
    while((Get-Date) -lt $deadline){
        if([PilotSleepRequest]::SetThreadExecutionState([uint32]2147483649) -eq 0){throw 'Windows refused the temporary sleep request.'}
        $p=Get-Content -LiteralPath (Join-Path $out 'progress.json') -Raw|ConvertFrom-Json
        [pscustomobject]@{PID=$PID;Updated=(Get-Date -Format 'yyyy-MM-dd HH:mm:ss');PilotState=$p.State;
            PreventSystemSleep=$true;DisplayMayTurnOff=$true;Deadline=$deadline.ToString('yyyy-MM-dd HH:mm:ss')}|
            ConvertTo-Json|Set-Content -LiteralPath $status -Encoding utf8
        if($p.State -eq 'complete'){break}
        if($p.State -eq 'stopped'){
            if($null -eq $stoppedSince){$stoppedSince=Get-Date}
            if(((Get-Date)-$stoppedSince).TotalMinutes -ge 30){break}
        }else{$stoppedSince=$null}
        Start-Sleep -Seconds 30
    }
}finally{
    [PilotSleepRequest]::SetThreadExecutionState([uint32]2147483648)|Out-Null
    [pscustomobject]@{PID=$PID;Updated=(Get-Date -Format 'yyyy-MM-dd HH:mm:ss');PreventSystemSleep=$false}|
        ConvertTo-Json|Set-Content -LiteralPath $status -Encoding utf8
}
