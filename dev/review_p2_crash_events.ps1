# Review of Phase 2 (dev/review_p2.md): list the Windows "Application Error"
# events of the last hours that name Rterm.exe, with the exception code and
# the faulting module of each, to say how the R processes of the dup group
# of dev/review_p2_cases.R died. Reads the Application log; changes nothing.
#   pwsh -File dev/review_p2_crash_events.ps1 [-Hours 6]
param([int]$Hours = 6)
$since = (Get-Date).AddHours(-$Hours)
$events = Get-WinEvent -FilterHashtable @{
    LogName      = 'Application'
    ProviderName = 'Application Error'
    StartTime    = $since
} -ErrorAction SilentlyContinue | Where-Object { $_.Message -match 'Rterm\.exe' }
if (-not $events) {
    Write-Output "no Application Error event names Rterm.exe in the last $Hours hours"
    exit 0
}
$rows = foreach ($e in $events) {
    $code = if ($e.Message -match 'Exception code:\s*(0x[0-9a-fA-F]+)') { $Matches[1] } else { '?' }
    $module = if ($e.Message -match 'Faulting module name:\s*([^,\r\n]+)') { $Matches[1] } else { '?' }
    [pscustomobject]@{ Time = $e.TimeCreated; Code = $code; Module = $module }
}
Write-Output ("{0} events naming Rterm.exe since {1}" -f $rows.Count, $since)
$rows | Group-Object Code, Module | ForEach-Object {
    Write-Output ("{0,4} x  exception {1}" -f $_.Count, $_.Name)
}
Write-Output "first and last:"
$rows | Sort-Object Time | Select-Object -First 1 -Last 1 | Format-Table -AutoSize | Out-String | Write-Output
