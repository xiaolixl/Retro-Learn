. "$PSScriptRoot\agent_config.ps1"

param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$QueryParts
)

$query = ($QueryParts -join " ").Trim()
if (-not $query) {
    $query = Read-Host "User request"
}

python "$PSScriptRoot\agent_cli.py" -q $query
