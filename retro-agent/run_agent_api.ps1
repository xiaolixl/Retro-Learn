. "$PSScriptRoot\agent_config.ps1"

python -m uvicorn agent_api:app --host 127.0.0.1 --port 8000
