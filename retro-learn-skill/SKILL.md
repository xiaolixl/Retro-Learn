---
name: retro-learn
description: "This skill should be used when the user asks about retrosynthesis, synthesis route planning, or how to synthesize a target molecule. Triggers include: retrosynthesis, synthesize SMILES, synthesis route, retrosynthetic analysis, target molecule, route planning, single-step retrosynthesis, multi-step retrosynthesis, synthesize from, synthesis of, how to synthesize, how to make, prepare molecule, synthetic pathway, 逆合成, 合成路线, 逆合成分析, 单步逆合成, 多步逆合成, 目标分子, 合成某个分子, 从…合成, 怎么合成, 合成路径, 制备分子, 正向合成. It runs the SimpRetro template-based retrosynthesis engine, supports single-step and multi-step (beam search) route planning with optional carbon-count constraints, and generates HTML visualizations with RDKit molecular structure drawings."
agent_created: true
---

# RetroLearn — SimpRetro Retrosynthesis Skill

Template-based retrosynthetic route planning for organic chemistry target molecules. Parse user's natural language → invoke SimpRetro engine → generate HTML visualization with RDKit molecular structures.

## Architecture

```
SimpRetro4Learn/        ← Engine (git submodule, data files only)
├── main.py, name2smiles.py, reaction_template.json, template_condition.json
├── emol_*.txt           ← Chemical stock databases
│
retro-learn-skill/       ← This skill package
├── SKILL.md
├── scripts/run_retro.py                ← Unified CLI (single + multi-step)
├── scripts/visualize.py                ← JSON → HTML with SVG + reaction type labels
├── scripts/tree_view.py                ← Tree-style visualization (alternative)
├── scripts/template_reaction_types.json ← LLM-preclassified template→type cache (55 entries)
├── engine/                             ← Engine adapter (retro_engine.py, route_planner.py)
└── references/          ← Detailed reference docs
    ├── build_template.md   ← Phase 4 build script template + CSS
    ├── engine_params.md    ← CLI params, scoring, stock DB, dedup, errors
    └── scenarios_detail.md ← Scenario A/B/C detailed steps & examples
```

**Deploy**: `SimpRetro4Learn/` and `retro-learn-skill/` must be sibling directories.

## Installation

```bash
pip install numpy==1.24.1 matplotlib==3.7.2 tqdm==4.65.0 rdkit-pypi==2023.3.3 rdchiral==1.1.0
```

> PyPI package is `rdkit-pypi`, not `rdkit`.

## Workflow

### Phase 1: Parse User Request

Extract from natural language:

| Parameter | Required | Default | How |
|-----------|----------|---------|-----|
| `target_smiles` | Yes | — | SMILES string in message, or resolve name via PubChem/CIR |
| `steps` | No | By scenario | "single-step"→1, "3-step"→3, unspecified→determined by scenario |
| `preferred_reactants` | No | [] | SMILES of preferred starting materials |
| `top_k` | No | 5 (single) / 3 (multi) | If user specifies "top-2", "top-3" |
| carbon constraint | No | By scenario | "≤4C" → `-db emol_under_4_carbons` |

SMILES pattern: continuous string with `C`, `c`, `N`, `O`, `=`, `#`, `(`, `)`, `[`, `]`, `@`, no spaces.

If ambiguous, ask before proceeding.

### Phase 2: Execute Retrosynthesis

Run `scripts/run_retro.py` via Bash. **`<project_root>`** = parent of both `SimpRetro4Learn/` and `retro-learn-skill/`.

**Scenario A** — Single-step (no multi-step keywords):
```bash
python retro-learn-skill/scripts/run_retro.py -s "SMILES" -k 5 -o retro_result.json
```
After results, ask user if they want to continue expanding any route.
Default stock DB: `emol_under_0_carbons` (unrestricted).

**Scenario B** — Multi-step without starting material or carbon constraint:
```bash
python retro-learn-skill/scripts/run_retro.py -s "SMILES" --steps 3 --per-step-top-k 3 --beam-width 8 -db emol_under_0_carbons -o retro_result.json
```
Default: 3 steps, unrestricted stock DB (`emol_under_0_carbons`). Honor user's step/top_k if specified.

**Scenario C** — Multi-step WITH starting material or carbon constraint:
```bash
# With starting material (NO carbon constraint) → unrestricted DB
python retro-learn-skill/scripts/run_retro.py -s "SMILES" --steps 5 --per-step-top-k 3 --beam-width 8 -pr "PREFERRED" -db emol_under_0_carbons -o retro_result.json

# With carbon constraint (NO starting material) → matching emol DB
python retro-learn-skill/scripts/run_retro.py -s "SMILES" --steps 3 --per-step-top-k 3 --beam-width 8 -db emol_under_4_carbons -o retro_result.json

# With BOTH starting material AND carbon constraint → matching emol DB + preferred
python retro-learn-skill/scripts/run_retro.py -s "SMILES" --steps 5 --per-step-top-k 3 --beam-width 8 -pr "PREFERRED" -db emol_under_4_carbons -o retro_result.json
```
Default: 5 steps. Direct beam search (no incremental trial).

**Carbon-count constraint mapping** (detect from "不超过X碳", "≤X carbon", "X碳以内", etc.):

| Max carbons | `-db` value |
|-------------|-------------|
| Unrestricted | `emol_under_0_carbons` |
| ≤3C | `emol_under_3_carbons` |
| ≤4C | `emol_under_4_carbons` |
| ≤6C | `emol_under_6_carbons` |

Full mapping: 1–6 carbons → `emol_under_{N}_carbons`.

**Database selection summary by scenario**:
- Scenario A (single-step): `emol_under_0_carbons` (unrestricted, default)
- Scenario B (multi-step, no constraints): `emol_under_0_carbons` (unrestricted)
- Scenario C (starting material, no carbon constraint): `emol_under_0_carbons` (unrestricted, starting material acts as stock constraint)
- Scenario C (carbon constraint only): matching `emol_under_{N}_carbons`
- Scenario C (both starting material + carbon constraint): matching `emol_under_{N}_carbons` + preferred reactant

**Direct beam search** (Scenario C): Run beam search directly with 5 steps (no incremental trial). Engine checks preferred reactant hit internally at each step. If unreached after 5 steps, use LLM knowledge for 1–2 routes.

**Route ordering rule** (Scenario C combined engine + LLM):

| Best SimpRetro score | Display order |
|---------------------|---------------|
| ≥ 0 | Engine routes first → LLM routes |
| < 0 or no engine routes | LLM routes first (max 2) → engine routes |

### Visual Language Convention (MANDATORY)

Follow the user's input language in all visual output:
- **English input** → English labels (reaction names, step titles, badges, legends)
- **Chinese input (中文输入)** → Chinese labels (中文标签、步骤标题、图例)

Applies to: Phase 3 `visualize.py`, Phase 4 custom build scripts, tree views, and any hand-crafted HTML visualizations produced for the user.

### Reaction Type Label Convention (MANDATORY)

Arrow labels: **reaction conditions above**, reaction type below. Follow these rules strictly:

1. **Only label known types**: oxidation, reduction, elimination, addition, and named reactions (e.g., Diels-Alder, Wittig, Grignard, Gilman). Never invent or guess.
2. **No parenthetical details**: use `Oxidation` not `Oxidation (alcohol to carbonyl)`. Keep labels clean and short.
3. **Unknown → omit**: if the reaction type cannot be confidently identified, leave the type label **blank** (show only conditions). **NEVER** display `Retro-synthetic disconnection` or any generic fallback label.

### Phase 3: Generate Visualization

```bash
# Default: forward-direction flowchart
python retro-learn-skill/scripts/visualize.py retro_result.json -o build_target_view.html
# With starting material
python retro-learn-skill/scripts/visualize.py retro_result.json -o build_target_from_source_view.html
# Tree view (only when user explicitly asks)
python retro-learn-skill/scripts/tree_view.py retro_result.json -o build_target_tree.html
```

Produces HTML with: target molecule SVG + SMILES + MW, forward-direction layout (reactants → arrow → product), single-arrow format (conditions above, reaction type below), ~60% scaled SVGs, color-coded reaction type labels.

**Multi-step rendering**: For multi-step routes, `visualize.py` now renders each forward step as a vertical stack of rows. Each step shows ALL reactants from `expanded_smiles` (with `+` separators for multi-reactant steps) → arrow with conditions → product. This correctly handles branched/converging routes where co-reactants are prepared by independent branches. Steps are labeled "Step 1", "Step 2", etc.

### Phase 4: Custom Visualization (use when needed)

Build a custom HTML page for situations where `visualize.py` output is insufficient: when you need annotated mechanism explanations, side-by-side route comparisons, explicit branch labels (Branch A/B), or Chinese/English manual labels. Use the standard template in `references/build_template.md` — includes `load_or_gen()`, `mb()`, `mb_target()`, `ar()` helpers and full CSS.

**Naming convention**:
- `build_{target_slug}_from_{source_slug}_view.html` — with starting material
- `build_{target_slug}_view.html` — no starting material
- Version: append `_v1`, `_v2` if file exists

**Target molecule display rule (MANDATORY)**: Target molecule uses `mb_target()` — same layout as `mb()` but with indigo (#4338CA) "target" pill label and indigo SVG border.

**Best route highlighting (MANDATORY)**: Highest-scoring route gets `engine` class + "Best Route" badge + "Best" tag pill.

### Phase 5: Present Results

1. Display HTML via `preview_url`
2. Text summary: target name/SMILES/MW, route count, best route conditions/score, stock status
3. Disclaimer: computational suggestions, not experimentally validated

## Reaction Type Classification

Reaction types are determined by **static template cache lookup** (no runtime LLM call). `run_retro.py` loads `scripts/template_reaction_types.json` (pre-classified by LLM in retro-agent), looks up each `reaction_template` in the cache, and injects `reaction_type` into the output JSON. `visualize.py` and `tree_view.py` consume `reaction_type` directly — no hardcoded SMARTS rules.

| Category | CSS Class | Color |
|----------|-----------|-------|
| Diels-Alder | `rxn-type-da` | Yellow |
| Esterification / Amidation / Acylation | `rxn-type-ester` | Pink |
| Coupling / Cross-coupling | `rxn-type-coupling` | Orange |
| Elimination / Dehydrohalogenation | `rxn-type-elim` | Blue |
| Reduction / Hydrogenation | `rxn-type-red` | Green |
| Oxidation | `rxn-type-ox` | Red |
| Halogenation / Bromination | `rxn-type-hal` | Purple |
| Hydrolysis | `rxn-type-hydro` | Teal |
| Substitution (SN1 / SN2) | `rxn-type-sub` | Steel-blue |
| Addition (Grignard / Wittig / Aldol) | `rxn-type-add` | Light-green |
| Alkylation | `rxn-type-alk` | Indigo |
| Unknown / Other | `rxn-type-default` | Gray |

To add new template classifications, update `template_reaction_types.json` (key = reaction_template SMARTS, value = type label) and regenerate via `retro-agent/classify_templates.py` if LLM classification is needed for new templates.

## Limitations

- Template engine is single-step; multi-step is iterative beam search
- Results are computational suggestions, not experimentally validated
- Template coverage is finite — complex/unusual structures may have no matches
- Prefer canonical SMILES input for best results
