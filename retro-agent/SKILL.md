---
name: retro-learn
description: "SimpRetro organic chemistry retrosynthesis analysis skill. Triggers: retrosynthesis, synthetic route, retrosynthetic analysis, target molecule, route planning, single-step, multi-step, 逆合成, 合成路线, 逆合成分析, 倒推合成, 单步逆合成, 多步逆合成, 目标分子. Runs the SimpRetro template-based retrosynthesis engine, supports single-step and multi-step (beam search) route planning, and generates HTML visualizations with RDKit molecular structure drawings."
---

# RetroLearn — SimpRetro Retrosynthesis Skill

## Overview

Perform template-based retrosynthetic route planning for organic chemistry target molecules. The AI agent (you) parses the user's natural language request, invokes the SimpRetro engine, and generates an HTML visualization with molecular structure drawings.

**What this skill does**:
- Parse the user's request to extract target molecule SMILES, step count, and preferred reactants
- Execute single-step retrosynthesis (returns top routes) or multi-step retrosynthesis (beam search)
- Generate an HTML page with RDKit-rendered SVG molecular structures
- Present the visual result and a concise summary to the user

**What this skill does NOT need**:
- An external LLM API key — you (the AI agent) ARE the LLM
- A separate agent_runtime or agent_cli layer — you handle natural language understanding natively

---

## Repository Structure

This is the `retro-learn` repository — an AI-agent skill layer on top of the SimpRetro engine.

```
retro-learn/                           ← This repository
├── README.md                          ← Repository overview & setup guide
├── LICENSE
├── .gitignore
├── .gitmodules                        ← References SimpRetro4Learn
│
├── SimpRetro4Learn/       ← Engine (git submodule)
│   ├── __init__.py                    ←   ↳ https://github.com/wzhstat/SimpRetro4Learn
│   ├── retro_engine.py               ←   Template matching + scoring
│   ├── route_planner.py              ←   Multi-step beam search
│   ├── main.py                       ←   Standalone engine CLI
│   ├── name2smiles.py                ←   Molecule name → SMILES
│   ├── reaction_template.json        ←   Reaction template library (SMARTS)
│   ├── template_condition.json       ←   Reaction conditions per template
│   ├── requirements.txt              ←   Engine dependencies
│   └── emol_*.txt                     ←   Chemical stock databases
│
├── retro_agent/                       ← Helper scripts (invoked by the AI agent)
│   ├── run_retro.py                   ←   Unified CLI (single-step + multi-step)
│   ├── visualize.py                   ←   JSON → HTML with SVG molecules + reaction type labels
│   └── install.bat                    ←   One-click dependency installer (Windows)
│
└── retro-learn/                       ← WorkBuddy skill package
    └── SKILL.md                       ←   This file (skill definition)
```

### Two-repo architecture

| Repository | Purpose | URL |
|------------|---------|-----|
| **SimpRetro4Learn** | Retrosynthesis engine (pure algorithm, no LLM) | `github.com/wzhstat/SimpRetro4Learn` · forked to `github.com/xiaolixl/SimpRetro4Learn` |
| **retro-learn** | AI-agent skill layer + helper scripts (this repo) | `github.com/xiaolixl/retro-learn` |

The engine is included as a **git submodule** so it stays independently versioned and citable.

---

## Installation

### Prerequisites

| Requirement | Notes |
|-------------|-------|
| Python 3.9+ | With pip |
| Git | For cloning with submodules |
| Dependencies | rdkit-pypi, rdchiral, numpy, matplotlib, tqdm |

### Clone with Submodule

```bash
git clone --recurse-submodules https://github.com/xiaolixl/retro-learn.git
cd retro-learn
```

If you already cloned without `--recurse-submodules`:

```bash
git submodule update --init --recursive
```

### Install Dependencies

**Windows (one-click)**:

```bash
cd retro_agent
.\install.bat
```

**Any platform (manual)**:

```bash
pip install numpy==1.24.1 matplotlib==3.7.2 tqdm==4.65.0 rdkit-pypi==2023.3.3 rdchiral==1.1.0
```

> **Note**: The engine's `requirements.txt` lists `rdkit==2023.3.3`, but the actual PyPI package name is `rdkit-pypi`.

### Verify

```bash
python -c "import rdkit; import rdchiral; print('OK')"
python retro_agent/run_retro.py -s "CC(=O)C=C(C)C" -k 3 -o test.json
```

---

## Workflow

### Phase 1: Parse User Request

From the user's natural language, extract:

| Parameter | Required | Default | How to Extract |
|-----------|----------|---------|---------------|
| `target_smiles` | Yes | — | Look for a SMILES string (contains C,N,O,=,#,(),[],@,/,\,.) in the user's message |
| `target_name` | If no SMILES | — | If user gives a name like "aspirin" or "mesityl oxide" instead of SMILES |
| `steps` | No | Depends on request type (see below) | "single-step" / "1-step" → 1; "3-step" → 3; unspecified → determined by request type |
| `preferred_reactants` | No | [] | SMILES of preferred starting materials if user mentions them |
| `top_k` | No | Depends on request type | If user specifies "top-2", "top-3" etc., use that |

**SMILES recognition**:
- Pattern: continuous string with characters like `C`, `c`, `N`, `O`, `=`, `#`, `(`, `)`, `[`, `]`, `@`, no spaces or Chinese characters
- Examples: `CC(=O)C=C(C)C`, `c1ccccc1`, `CC(=O)Oc1ccccc1C(=O)O`
- If user provides only a name, attempt to resolve it using an online database (PubChem, CIR). If resolution fails, ask the user to provide a SMILES string.

**If the user gives an ambiguous request** (no clear target), ask a follow-up question before proceeding.

### Phase 2: Execute Retrosynthesis

Run `retro_agent/run_retro.py` using the Bash tool.

The execution path depends on the request type. Three scenarios:

---

#### Scenario A — Single-Step Request

User asks for retrosynthesis without suggesting multiple steps.

```bash
cd <project_root>
python retro_agent/run_retro.py \
  -s "CC(=O)C=C(C)C" \
  -k 5 \
  -o retro_result.json
```

- Returns the **top 5** ranked single-step routes
- Each route contains: score, reaction template, reaction conditions, and reactant SMILES
- **After presenting results, ask the user**: "Would you like to continue the retrosynthesis from any of these routes? Please tell me which route number to explore further."
- If the user picks a route, run single-step retrosynthesis on that route's reactant(s) and continue iteratively

---

#### Scenario B — Multi-Step Request WITHOUT a specified starting material

User asks for multi-step retrosynthesis (e.g. "3-step" or just "multi-step") but does not provide a preferred starting material.

| User says | Default behavior |
|-----------|-----------------|
| "3-step retrosynthesis" | 3 steps, top-1 (per_step_top_k=1, beam_width=5) |
| "multi-step retrosynthesis" (no number) | 3 steps, top-1 |
| "2-step retrosynthesis, top-2 per step" | 2 steps, per_step_top_k=2 |
| "5-step retrosynthesis, top-3" | 5 steps, per_step_top_k=3 |

**Default command** (user did not specify steps or top_k):

```bash
cd <project_root>
python retro_agent/run_retro.py \
  -s "CC(=O)C=C(C)C" \
  --steps 3 \
  --per-step-top-k 1 \
  -o retro_result.json
```

**If the user specifies steps and/or top_k**, honor those values. For example:

```bash
# User says: "2-step retrosynthesis, top-2 per step"
python retro_agent/run_retro.py \
  -s "CC(=O)C=C(C)C" \
  --steps 2 \
  --per-step-top-k 2 \
  --beam-width 8 \
  -o retro_result.json
```

---

#### Scenario C — Multi-Step Request WITH a specified starting material

User asks for multi-step retrosynthesis and provides a preferred starting material (e.g. "from butane", "starting from toluene"), but does not specify the number of steps.

**Iterative expanding search**: Run single-step, check if the starting material is reached, expand if not — repeat up to 5 steps.

**Algorithm**:
1. Start with `current_steps = 1`
2. Run single-step retrosynthesis with `top_k=5`
3. Check if any route's reactants contain the preferred starting material
4. If **yes**: output the route(s) that hit the starting material. Done.
5. If **no**: run multi-step with `steps = current_steps + 1`, per_step_top_k=1, beam_width=5
6. Check the final leaf reactants of the best route for the starting material
7. If **yes**: output and stop
8. If **no**: increment current_steps, repeat up to 5
9. If the starting material is not reached within 5 steps: show the 5-step result and explain that the starting material was not found

**Implementation** — run sequentially, checking results between steps:

```bash
# Step 1: single-step, check for starting material
python retro_agent/run_retro.py -s "TARGET" -k 5 -o step1.json

# Step 2 (if not found): 2-step beam search
python retro_agent/run_retro.py -s "TARGET" --steps 2 --per-step-top-k 1 -o step2.json

# ... continue up to step 5 if needed
```

**If the user also specifies the number of steps**, use that directly instead of the iterative search:

```bash
# User says: "3-step retrosynthesis from toluene"
python retro_agent/run_retro.py \
  -s "C1(/C=C/C2=CC=CC=C2)=CC=CC=C1" \
  --steps 3 \
  --per-step-top-k 1 \
  -pr Cc1ccccc1 \
  -o retro_result.json
```

#### Full Parameter Reference

| Flag | Description | Default |
|------|-------------|---------|
| `-s` / `--smiles` | Target molecule SMILES (required) | — |
| `--steps` | Retrosynthesis depth | 1 |
| `-k` / `--top-k` | Routes to return (single-step) | 5 |
| `-db` / `--database` | Stock database name | `emol_under_0_carbons` |
| `-tpl` / `--template` | Template file name | `reaction_template.json` |
| `-cond` / `--condition` | Condition file name | `template_condition.json` |
| `-w` / `--weights` | 4 scoring weights | `0.1 0.2 0.5 0.0` |
| `-pr` / `--preferred-reactants` | Preferred reactant SMILES | (none) |
| `--beam-width` | Beam width (multi-step) | 5 |
| `--per-step-top-k` | Top-K per step (multi-step) | 5 |
| `-o` / `--output` | Output JSON path | `retro_result.json` |

### Phase 3: Generate Visualization

Run the visualizer to convert the JSON result into an HTML page:

```bash
cd <project_root>
python retro_agent/visualize.py retro_result.json -o retro_result_view.html
```

This produces an HTML page with:
- **Target molecule**: SVG structure drawing + SMILES + molecular weight
- **Forward-direction layout**: Reactions are displayed in **forward synthesis direction** (starting materials → arrow → product), not retrosynthetic (product ← arrow ← reactants)
- **Single-arrow format**: Each reaction step uses a single long arrow. Reagents, solvents, and temperature are written **above** the arrow; the reaction type tag is written **below** the arrow
- **Scaled SVG structures**: Molecular structure images are rendered at ~60% scale (max-width ~100px) for a compact layout
- **Single-step mode**: Route cards showing reactant structures → arrow → target molecule, with scores, conditions, stock status, and reaction type
- **Multi-step mode**: Step-by-step chain showing the full forward sequence from starting materials → intermediate(s) → target
- Best route highlighted in green
- **Reaction type labels**: Each arrow includes a color-coded tag classifying the reaction type (e.g. "Reductive debromination", "Retro Diels-Alder", "Elimination", "Oxidation", "Halogenation"). Classification is inferred from the SMARTS template pattern and reaction conditions by `visualize.py`'s built-in `classify_reaction()` function.

**Reaction type classification categories**:

| Tag Color | Category | Examples |
|-----------|----------|---------|
| Green | Reduction | Reductive debromination, Reductive dechlorination, Carbonyl reduction, Catalytic hydrogenation |
| Blue | Elimination | Elimination (dehydrohalogenation) |
| Red | Oxidation | Oxidation (alcohol → carbonyl) |
| Purple | Halogenation | Bromination, Chlorination, Halogenation |
| Yellow | Cycloaddition | Retro Diels-Alder |
| Gray | General | Retro-synthetic disconnection (fallback) |

### Phase 4: Build Custom Visualization Page with Embedded RDKit SVGs

For a polished, presentation-ready result, build a dedicated HTML page that embeds complete RDKit SVG molecular structures directly into the synthetic route diagram.

**Naming convention** for the output file:
- `build_{target_slug}_from_{source_slug}_view.html` — when both target and starting material are specified
- `build_{target_slug}_view.html` — when no starting material is specified
- Use a short slug of the molecule name or SMILES (e.g. "stilbene_from_toluene", "butene_from_butane")
- **Versioning**: if the exact filename already exists on disk, append `_v1`, `_v2` etc. **after** `view` (e.g. `build_stilbene_from_toluene_view_v1.html`). First-time generation omits the version number.

**Examples**:
- First run, toluene to stilbene: `build_stilbene_from_toluene_view.html`
- Second run, same molecules: `build_stilbene_from_toluene_v1_view.html`
- No starting material, first run: `build_mesityl_oxide_view.html`

**Workflow**:

1. Write a Python script (`build_*.py` matching the HTML filename) that:
   a. Reads existing SVG files from `molecule_svgs/` directory
   b. For any molecule without an existing SVG, generates one using RDKit (`Chem.MolFromSmiles` → `AllChem.Compute2DCoords` → `rdMolDraw2D.MolDraw2DSVG`)
   c. Embeds all SVG files **completely and unmodified** into the HTML page (only strip the `<?xml...?>` processing instruction — keep ALL bond paths, atom label paths, double-bond indicators, and white background rect)
   d. Displays the route in **forward direction** (starting materials → arrow → product)
   e. Uses **single-arrow** format for each reaction step, with conditions written **above** the arrow and reaction type tag **below**
   f. Scales SVG structures to ~60% (max-width ~100px)

2. Execute the script using Python 3.9:
   ```bash
   cd <project_root>
   "C:\Program Files\Python39\python.exe" build_stilbene_from_toluene_view.py
   ```

3. Display the generated HTML file using `preview_url`

4. After displaying, offer to clean up the build script by deleting it (unless it may be useful to keep for reference or reruns).

### Phase 5: Present Results

1. Display the HTML file using `preview_url`
2. Provide a text summary:
   - Target molecule (name if identifiable, SMILES, MW)
   - Number of recommended routes (single-step: top 3; multi-step: 1 best)
   - Best route's reaction conditions and score
   - Whether reactants are in stock
3. Add a note: results are heuristic suggestions, not experimentally validated

---

## Scoring System

The engine uses 4 weighted scoring dimensions:

| Dimension | Default Weight | Description |
|-----------|---------------|-------------|
| Complexity Decrease (cd_score) | 0.1 | How much simpler products are vs. reactants |
| Availability (availability_score) | 0.2 | Whether reactants exist in the stock database |
| Ring Disconnection (ring_disconnection_score) | 0.5 | Whether ring structures are successfully opened |
| Template Frequency | 0.0 | How selective the template is (currently disabled) |

Default weights `[0.1, 0.2, 0.5, 0.0]` prioritize ring disconnection — ideal for organic chemistry teaching.

---

## Multi-Step Beam Search

When `steps ≥ 2`, the planner uses beam search:

- Each round: expand the non-stock leaf with the most atoms, generate single-step candidates
- Retain top `beam_width` candidates ranked by: completion → preferred reactant match → stock ratio → score
- Final output: the single best cumulative route

---

## Error Handling

| Situation | Response |
|-----------|----------|
| Invalid SMILES | "The SMILES format appears invalid. Please provide a canonical SMILES string." |
| No template matches | "The template library does not cover this molecule. Try simplifying the structure." |
| rdkit/rdchiral import fails | "Dependencies not installed. Run `retro_agent/install.bat` first." |
| Multi-step finds no route | "No viable multi-step route found. Try reducing the step count or removing preferred reactant constraints." |
| Molecule name unresolvable | "Could not resolve this molecule name. Please provide the SMILES string directly." |

---

## Common User Queries

| User Says | Action |
|-----------|--------|
| "Run retrosynthesis on CC(=O)C=C(C)C" | Single-step, top-5 → ask user which route to continue from |
| "3-step retrosynthesis for mesityl oxide" | Resolve name → SMILES, multi-step 3 steps top-1 |
| "multi-step retrosynthesis from butane" | Iterative search up to 5 steps, check if butane is reached at each step |
| "2-step retrosynthesis, top-2 from ethanol" | 2 steps, per_step_top_k=2, preferred=CCO |
| "Analyze synthetic routes for this molecule" | Ask for SMILES or name and whether single or multi-step |
| "Prefer ethanol as starting material" | Set `-pr CCO` |
| "Continue from Route 2" | Run single-step on Route 2's reactant SMILES with top_k=5 |
| "Show the last result" | Read the most recent `retro_result.json` and visualize |

---

## Limitations

- Template engine is inherently single-step; multi-step is iterative beam search on top
- Results are **heuristic suggestions**, not experimentally validated protocols
- Target must be a valid molecule parseable by RDKit
- Template coverage is finite — complex/unusual structures may not have matches
- For best results, users should provide canonical SMILES directly

---

## Tips

1. **SMILES > names**: Direct SMILES input is more reliable than molecule names
2. **Start with single-step**: Verify basic results before attempting multi-step
3. **Adjust weights**: Increase availability weight for practical synthesis; keep ring disconnection high for teaching
4. **Stock databases**: `emol_under_0_carbons` = unrestricted; `emol_under_5_carbons` = ≤5 carbons
5. **Canonical SMILES**: Use `rdkit.Chem.MolToSmiles(rdkit.Chem.MolFromSmiles(user_input))` to canonicalize if needed
