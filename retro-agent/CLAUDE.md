# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SimpRetro — a natural-language retrosynthesis agent. Users ask chemistry questions in Chinese or English, an LLM parses the request, and a template-based engine proposes synthetic routes with structure images.

## Architecture (two layers)

**Agent layer** (root `.py` files) — LLM orchestration:
- `agent_cli.py` — CLI entry point
- `agent_api.py` — FastAPI HTTP entry point (endpoints: `GET /`, `GET /health`, `POST /agent/query`)
- `agent_runtime.py` — core agent: `RetrosynthesisAgent` class handles parsing → resolution → planning → explanation → rendering
- `agent_rendering.py` — RDKit molecule structure image generation (grid images for routes/steps)
- `chem_resolution.py` — SMILES canonicalization and name-to-structure resolution (delegates to `SimpRetro4Learn.name2smiles` when available)

**Engine layer** (`SimpRetro4Learn/`) — single-step retrosynthesis:
- `SimpRetro4Learn/main.py` — loads templates and in-stock molecules, scores routes via CDScore/ASScore/RDScore/SiteSelectivity, outputs top-5 ranked routes as JSON
- Template files: `reaction_template.json`, `template_condition.json`
- In-stock databases: `emol_under_{1-6}_carbons.txt`

**Agent flow**: `parse_user_request` (LLM extracts JSON: target SMILES, step count, preferred reactants) → `resolve_request` (canonicalize SMILES, validate) → `plan_retrosynthesis` (call engine) → `explain_results` (LLM generates natural-language reply) → `render_agent_artifacts` (RDKit images).

## Key technical details

- LLM calls use OpenAI Responses API by default; auto-detects DeepSeek endpoints and falls back to Chat Completions API with `response_format: json_object`
- Default model: `OPENAI_MODEL` env var or `gpt-5.4-mini`
- Single-step mode returns top-3 routes; multi-step returns the single best cumulative route
- Outputs land in `user_output/agent_runs/<timestamp>_<id>/` with parsed_request.json, resolved_request.json, planning_result.json, agent_reply.md, agent_result.json, and structure PNGs
- Multi-step retrosynthesis: the agent repeatedly calls single-step engine, each step feeding the previous step's reactant as the new target

## Known issue: import path mismatch

`agent_runtime.py` imports from `SimpRetro4OrganicChemistryB.retro_engine` and `SimpRetro4OrganicChemistryB.route_planner`, but the actual directory on disk is `SimpRetro4Learn` and these modules (`retro_engine.py`, `route_planner.py`) do not exist. The engine's real implementation is in `SimpRetro4Learn/main.py` (functions: `load_global_data`, `perform_retrosynthesis`). Expect to create bridge modules or fix imports before the agent can run.

## Commands

```powershell
# Install dependencies
pip install -r requirements.txt

# CLI (requires OPENAI_API_KEY / OPENAI_BASE_URL / OPENAI_MODEL env vars)
python agent_cli.py -q "请对目标分子 CC(=O)C=C(C)C 做单步逆合成"
python agent_cli.py -q "..." --json          # full JSON output
python agent_cli.py -q "..." --model gpt-5.4  # model override

# API server
python -m uvicorn agent_api:app --host 127.0.0.1 --port 8000

# Standalone engine (no LLM)
python SimpRetro4Learn/main.py -s "CC(=O)C=C(C)C" -db emol_under_1_carbons -o result.json

# Standalone SMILES-to-image
python smiles_to_image.py -s "CCO" -o ./images

# Convenience scripts (PowerShell)
.\run_agent.ps1 -QueryParts "请对目标分子 CC(=O)C=C(C)C 做单步逆合成"
.\run_agent_api.ps1
```

## Environment variables

- `OPENAI_API_KEY` — required
- `OPENAI_MODEL` — optional, defaults to `gpt-5.4-mini`
- `OPENAI_BASE_URL` — optional, for API-compatible proxies (DeepSeek, etc.)

## Engine scoring weights

Four floats `[w1, w2, w3, w4]` default `[0.1, 0.2, 0.5, 0.0]`:
1. **CDScore** — favors convergent synthesis (cleaving into equal-sized fragments)
2. **ASScore** — favors commercially available reactants, penalizes organometallics (Mg/Li/Zn)
3. **RDScore** — bonus for ring-opening disconnections
4. **Site Selectivity** — penalizes templates that produce multiple isomers
