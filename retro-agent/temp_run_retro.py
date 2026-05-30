#!/usr/bin/env python3
import subprocess
import sys
import os

# UTF-8 wrapper to properly pass Chinese characters
cmd = [
    "conda", "run", "-n", "retro_env", 
    "python", 
    r"C:/Users/13558/.openclaw/workspace/skills/organic-chemistry-retrosynthesis/scripts/process_retro.py",
    "奎宁",
    "-n", "3"
]

print("Executing:", " ".join(cmd))
print("=" * 70, file=sys.stderr)
result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', cwd=os.path.expanduser("~"))
print(result.stdout)
print(result.stderr, file=sys.stderr)
sys.exit(result.returncode)