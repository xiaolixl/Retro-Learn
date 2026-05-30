# SimpRetro Retrosynthesis Agent

This is the root-level guide for the entire `local_retro` project.
If your goal is to enable the full agent, start from this root directory instead of directly running the lower-level scripts inside `SimpRetro4Learn/`.

The full workflow is:

1. The user submits a natural-language request
2. The LLM extracts the target molecule, preferred reactants, and requested route depth
3. The agent calls the retrosynthesis engine
4. The agent generates a natural-language explanation
5. The agent saves structure images and result files
6. (Via Web UI) visually displays synthetic routes with molecular structure drawings alongside SMILES

## 1. Prerequisites

Before enabling the agent, make sure you have:

- Python 3.9
- A working Python environment where dependencies can be installed
- A valid OpenAI API key (or compatible endpoint, e.g. DeepSeek)
- The current working directory set to the project root: `local_retro/`

Recommended structure:

```text
local_retro/
├── README.md
├── README_EN.md
├── SKILL.md
├── requirements.txt
├── agent_cli.py
├── agent_api.py
├── agent_runtime.py
├── chem_resolution.py
├── agent_rendering.py
├── agent_config.ps1
├── run_agent.ps1
├── run_agent_api.ps1
├── smiles_to_image.py
├── static/
│   └── index.html          ← Web frontend
├── user_output/
│   └── agent_runs/         ← per-run output
└── SimpRetro4Learn/        ← Retrosynthesis engine (template matching + scoring)
    ├── main.py
    ├── retro_engine.py     ← Engine constants (DB name, weights)
    ├── route_planner.py    ← Route planner (single-step + multi-step)
    ├── name2smiles.py      ← Molecule name → SMILES resolver
    ├── reaction_template.json
    ├── template_condition.json
    ├── emol_under_*.txt
    └── template/
```

## 2. Create an Environment

Using `conda` is recommended:

```powershell
conda create -n retro_agent python=3.9
conda activate retro_agent
```

If you prefer `venv`:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

## 3. Install Dependencies

From the project root, run:

```powershell
pip install -r requirements.txt
```

The root `requirements.txt` includes:

- root-level agent dependencies
- the retrosynthesis engine dependencies from `SimpRetro4Learn/`

If `rdkit` or `rdchiral` fails to install, fix those lower-level dependencies first and then re-run the install command. The agent cannot run correctly until the retrosynthesis engine dependencies are installed.

## 4. Configure LLM Environment Variables

Since your environment is PowerShell, use:

```powershell
$env:OPENAI_API_KEY="your_openai_api_key"
$env:OPENAI_MODEL="gpt-5.4-mini"
```

If you are using an OpenAI-compatible gateway, you can also set:

```powershell
$env:OPENAI_BASE_URL="https://your-compatible-endpoint"
```

You can also edit `agent_config.ps1` directly with your settings, then launch via `run_agent.ps1` / `run_agent_api.ps1`.

Notes:

- `OPENAI_API_KEY`: required
- `OPENAI_MODEL`: optional
- `OPENAI_BASE_URL`: optional, only for custom gateways

## 5. Three Ways to Enable the Agent

You can enable the agent in three ways:

- Option A: command-line mode
- Option B: HTTP API mode
- Option C: HTTP API + Web UI (recommended)

### Option A: Command-Line Mode

This is the simplest way to start the agent.

Run from the project root:

```powershell
python agent_cli.py -q "Run single-step retrosynthesis for target molecule CC(=O)C=C(C)C"
```

For a multi-step request:

```powershell
python agent_cli.py -q "Run a 3-step retrosynthesis for target molecule CC(=O)C=C(C)C"
```

To prefer certain starting materials:

```powershell
python agent_cli.py -q "Run a 2-step retrosynthesis for target molecule CC(=O)C=C(C)C and prefer ethanol and ethyl acetate as starting materials"
```

The agent will automatically:

1. detect the user language
2. extract the target molecule
3. extract the requested step count
4. extract preferred reactants
5. call the retrosynthesis planner
6. generate a natural-language explanation
7. save structure images to `user_output/agent_runs/...`

### Option B: HTTP API Mode

If you want to connect the agent to a frontend, chat UI, or another service, use the API mode.

Run from the project root:

```powershell
uvicorn agent_api:app --host 0.0.0.0 --port 8000
```

Available endpoints:

- `GET /` — Web frontend
- `GET /health` — health check
- `POST /agent/query` — submit a natural-language retrosynthesis query

First verify the service:

```powershell
Invoke-RestMethod -Uri "http://127.0.0.1:8000/health"
```

Then send a natural-language query:

```powershell
$body = @{
  message = "Run a 2-step retrosynthesis for target molecule CC(=O)C=C(C)C and prefer CCO as a starting material"
} | ConvertTo-Json

Invoke-RestMethod `
  -Uri "http://127.0.0.1:8000/agent/query" `
  -Method Post `
  -ContentType "application/json" `
  -Body $body
```

### Option C: Web UI (Recommended)

After starting the API, open `http://127.0.0.1:8000` in a browser to access the visual interface.

Web UI features:

- Natural-language query input (Chinese / English)
- Step count selector (1 step returns 3 routes; multi-step returns the single best route)
- Quick-example chips for one-click testing
- Routes displayed in **forward synthesis direction** (reactants → product), with reaction conditions written above the arrow and reaction type tag below
- Each route card shows molecular structure images (SVG/PNG) alongside reactant SMILES
- Score, stock status, and preferred reactant match info at a glance
- Best route (rank 1) highlighted with a green border
- Automatic reaction type classification with color-coded tags: Reduction / Elimination / Oxidation / Halogenation / Cycloaddition / Coupling

## 6. Default Agent Behavior

These rules control how the agent behaves.

### Step-Count Rules

- If the user does not specify the number of steps, the default is `1`
- If the user asks for `1` step, the agent returns the best `3` routes
- If the user asks for `2` or more steps, the agent repeatedly calls the single-step engine
- In multi-step mode, the agent returns only the single highest-ranked overall route

### Preferred Reactant Rules

- If the user provides preferred starting materials, matching routes are ranked first
- This is a ranking preference, not a hard filter
- If a preferred reactant name cannot be resolved reliably, the agent will ask for SMILES

### Molecule Representation Rules

- The most reliable input format is SMILES
- Molecule names are resolved via cirpy + pubchempy when possible
- If name resolution fails, the agent will ask the user to provide a structure explicitly

## 7. Recommended Input Style

Recommended:

```text
Please run single-step retrosynthesis for target molecule CC(=O)C=C(C)C.
```

Or:

```text
Please run a 3-step retrosynthesis for target molecule CC(=O)C=C(C)C and prefer ethanol as a starting material.
```

Chinese also works:

```text
请对目标分子 CC(=O)C=C(C)C 做 2 步逆合成，优先考虑乙醇作为原料。
```

If you only say "retro-synthesize aspirin", the agent may resolve it, but direct SMILES is more reliable.

## 8. What the Agent Produces

Each run writes output to:

```text
user_output/agent_runs/<timestamp>_<id>/
```

Typical files include:

- `parsed_request.json` — structured request from LLM
- `resolved_request.json` — validated execution parameters
- `planning_result.json` — retrosynthesis planning result
- `agent_reply.md` — natural-language explanation
- `agent_result.json` — complete output
- `target.png` / `target.svg` — target molecule structure (PNG + SVG)
- `route_x/` or `step_x/` — per-route / per-step reactant structures (PNG + SVG)
- `final_leaf_reactants/` — final starting material structures (multi-step mode)

So the agent does not only print text. It also saves structured outputs and structure images in both PNG and SVG formats. The Web UI embeds these SVG structures directly into synthetic route diagrams.

## 9. How to Confirm the Agent Is Enabled

The agent is considered successfully enabled if any of the following is true:

- `python agent_cli.py -q "..."` returns a natural-language answer
- `uvicorn agent_api:app ...` starts successfully
- `GET /health` returns `{"status": "ok"}`
- `POST /agent/query` returns `reply_markdown`
- a new directory appears under `user_output/agent_runs/`
- the Web UI loads at `http://127.0.0.1:8000`

## 10. Troubleshooting

### 1. `OPENAI_API_KEY is not set`

You have not configured the API key in the current PowerShell session. Run:

```powershell
$env:OPENAI_API_KEY="your_openai_api_key"
```

### 2. Missing `rdchiral`, `rdkit`, or related dependencies

The lower-level retrosynthesis engine is not installed correctly yet. Fix dependency installation first.

### 3. The user gives a molecule name and resolution fails

Provide the target molecule as SMILES.
For structure-based chemistry tasks, SMILES is more reliable than natural-language names.

### 4. The API starts but returns no useful route

Common reasons:

- invalid target structure
- template coverage does not support the molecule well
- preferred-reactant constraints reduce the chance of high-scoring matches

### 5. Why does multi-step mode return only one route?

That is the intended behavior:
multi-step mode searches multiple candidates internally, but only the single highest-ranked route is returned.

## 11. The Real Entry Points

If your goal is to enable the full agent, the correct entry points are:

- `agent_cli.py` — command line
- `agent_api.py` — HTTP API + Web UI

Not the lower-level file inside the subdirectory:

- `SimpRetro4Learn/main.py`

That is more suitable for algorithm-level testing. The root-level files are the real LLM agent entry points.

## 12. Current Limitations

- The core template engine is still fundamentally single-step
- Multi-step behavior is produced by upper-layer iterative search
- Outputs are heuristic suggestions, not experimentally validated plans
- If the target structure cannot be resolved reliably, the user still needs to provide SMILES

## 13. Shortest Activation Path

If you only want the fastest possible path to enable the agent, do these five steps:

```powershell
conda create -n retro_agent python=3.9
conda activate retro_agent
pip install -r requirements.txt
$env:OPENAI_API_KEY="your_openai_api_key"
python agent_cli.py -q "Run single-step retrosynthesis for target molecule CC(=O)C=C(C)C"
```

To use the Web UI, add one more step:

```powershell
python -m uvicorn agent_api:app --host 127.0.0.1 --port 8000
# open http://127.0.0.1:8000 in a browser
```

At that point, the agent is enabled.
