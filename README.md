# RetroLearn — Retro-Synthesis AI Agent Skill & CLI Tools

**An AI-agent skill layer on top of the [SimpRetro4Learn](https://github.com/wzhstat/SimpRetro4Learn) retrosynthesis engine.**

This repository provides an AI agent-compatible SKILL (`retro-learn/SKILL.md`) and helper scripts (`retro_agent/`) that enable an AI agent (Claude Code, Codex, WorkBuddy, etc.) to perform template-based retrosynthetic route planning for organic chemistry target molecules using natural language.

---

## Repository Structure

```
Retro-Learn/                           ← This repository (xiaolixl/Retro-Learn)
├── README.md                          ← This file
├── .gitignore
├── .gitmodules
│
├── SimpRetro4Learn/       ← Engine (git submodule)
│   └── (Template matching + scoring + beam search)
│
├── retro_agent/                       ← Helper scripts
│   ├── run_retro.py                   ←   Unified CLI (single-step + multi-step)
│   ├── visualize.py                   ←   JSON → HTML with SVG molecules + reaction type labels
│   └── install.bat                    ←   One-click dependency installer (Windows)
│
└── retro-learn/                       ←   skill package
    └── SKILL.md                       ←   Skill definition file
```

---

## Engine Description

**Two-repo architecture:**

| Repository | Purpose | Utility |
|------------|---------|---------|
| [SimpRetro4Learn](https://github.com/wzhstat/SimpRetro4Learn) | Retrosynthesis engine (pure algorithm, no LLM) | Can be used as a stand-alone package without AI |
| [Retro-Learn](https://github.com/xiaolixl/Retro-Learn) | AI-agent skill layer + helper scripts (this repo) | Requires AI agent (WorkBuddy / Claude Code) |

---

## Quick Start

### Prerequisites

- **Python 3.9+** with pip
- **Git** (for cloning with submodules)

### Clone

```bash
git clone --recurse-submodules https://github.com/xiaolixl/retro-learn.git
cd retro-learn
```

> If you already cloned without `--recurse-submodules`, run:
> ```bash
> git submodule update --init --recursive
> ```

### Install Dependencies

**Windows (one-click):**

```powershell
cd retro_agent
.\install.bat
```

**Any platform (manual):**

```bash
pip install numpy==1.24.1 matplotlib==3.7.2 tqdm==4.65.0 rdkit-pypi==2023.3.3 rdchiral==1.1.0
```

> **Note**: The engine's `requirements.txt` lists `rdkit==2023.3.3`, but the actual PyPI package name is `rdkit-pypi`.

### Verify Installation

```bash
python -c "import rdkit; import rdchiral; print('OK')"
python retro_agent/run_retro.py -s "CC(=O)C=C(C)C" -k 3 -o test.json
```

---
## Usage

### 🤖 Use with AI Agent (natural language interface)

The `retro-learn` skill enables AI agents (WorkBuddy, Claude Code, or any LLM with tool-calling ability) to perform retrosynthesis analysis through natural language conversation.

**How it works:**

1. The user describes the target molecule and desired route in natural language (e.g. "Find a 3-step retrosynthesis route for CC(=O)C=C(C)C", or "Find a synthesis route for stillbene from tolune")
2. The AI agent parses the request — converting chemical names to SMILES, extracting the target SMILES, step count, and any preferred starting materials
3. The agent runs `run_retro.py` to execute the retrosynthesis calculation (calling the main.py under SimpRetro4Learn)
4. The agent generates an HTML visualization page with RDKit molecular structures
5. The agent displays the result and provides a educational interpretation


**Claude Code / other AI agents:**

Import the `RetroLearn_SKILL.md` content as a system instruction or custom instruction (you need to rename the SKILL.md). The AI agent will follow the defined workflow to invoke `run_retro.py` and `visualize.py` using its tool-calling capabilities. 

**WorkBuddy:**

Upload `retro-learn/SKILL.md` (packaged as `retro-learn.zip`) to WorkBuddy via **Add Skill → Upload Skill**. Once installed, the skill will automatically activate when the user mentions retrosynthesis-related queries.

### ⚙️Pure Algorithm (no LLM required)

Run single-step retrosynthesis:

```bash
cd retro-learn/
python retro_agent/run_retro.py -s "CC(=O)C=C(C)C" -k 5 -o result.json
python retro_agent/visualize.py result.json -o result_view.html
```

Multi-step retrosynthesis:

```bash
python retro_agent/run_retro.py -s "CC(=O)C=C(C)C" --steps 3 -o result.json
```

---
