#!/usr/bin/env python3
"""
SimpRetro unified CLI entry point.
Supports single-step and multi-step retrosynthesis.

Usage:
  python run_retro.py -s "CC(=O)C=C(C)C"                          # single-step, top 5
  python run_retro.py -s "CC(=O)C=C(C)C" --steps 3                # multi-step (beam search)
  python run_retro.py -s "CC(=O)C=C(C)C" -k 3 -o result.json      # custom top-k and output
"""

import argparse
import json
import os
import sys


def _resolve_engine_dir():
    """Find the SimpRetro engine directory.

    Search order:
    1. ENGINE_DIR environment variable
    2. Walk up from this script's directory to find SimpRetro4Learn/
       (handles both retro-learn-skill/scripts/ and retro_agent/ layouts)
    3. SimpRetro4Learn/ as a sibling of the current working directory
       Fallback for development or custom setups.
    """
    # 1. Explicit override
    env_dir = os.environ.get("ENGINE_DIR")
    if env_dir and os.path.isfile(os.path.join(env_dir, "main.py")):
        return os.path.abspath(env_dir)

    # 2. Walk up from script dir to find SimpRetro4Learn/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    search_dir = script_dir
    for _ in range(4):  # search up to 4 levels
        candidate = os.path.join(search_dir, "SimpRetro4Learn")
        if os.path.isdir(candidate) and os.path.isfile(os.path.join(candidate, "main.py")):
            return os.path.abspath(candidate)
        parent = os.path.dirname(search_dir)
        if parent == search_dir:
            break
        search_dir = parent

    # 3. Fallback: look in cwd
    candidate = os.path.join(os.getcwd(), "SimpRetro4Learn")
    if os.path.isdir(candidate) and os.path.isfile(os.path.join(candidate, "main.py")):
        return os.path.abspath(candidate)

    print("[ERROR] Cannot find SimpRetro engine directory.", file=sys.stderr)
    print("Expected SimpRetro4Learn/ with main.py as:", file=sys.stderr)
    print(f"  - sibling of the repository root", file=sys.stderr)
    print(f"  - sibling of current working directory", file=sys.stderr)
    print(f"  - or set ENGINE_DIR env variable", file=sys.stderr)
    sys.exit(1)


ENGINE_DATA_DIR = _resolve_engine_dir()  # SimpRetro4Learn/ submodule (data files only)

# Ensure engine package (retro-learn-skill/engine/) is importable
_script_dir = os.path.dirname(os.path.abspath(__file__))
_skill_root = os.path.dirname(_script_dir)  # retro-learn-skill/
if _skill_root not in sys.path:
    sys.path.insert(0, _skill_root)
# Also keep repo root on path (for backward compat)
_repo_root = os.path.dirname(ENGINE_DATA_DIR)
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from engine.retro_engine import (
    DEFAULT_CONDITION_FILE,
    DEFAULT_DATABASE,
    DEFAULT_TEMPLATE_FILE,
    DEFAULT_WEIGHTS,
    run_retrosynthesis,
)
from engine.route_planner import plan_retrosynthesis


def run_single_step(args):
    result = run_retrosynthesis(
        smiles=args.smiles,
        database_name=args.database,
        template_file=args.template,
        condition_file=args.condition,
        weights=args.weights,
        preferred_reactants=args.preferred_reactants,
        top_k=args.top_k,
        base_dir=ENGINE_DATA_DIR,
        show_progress=True,
    )
    return result


def run_multi_step(args):
    result = plan_retrosynthesis(
        smiles=args.smiles,
        steps=args.steps,
        database_name=args.database,
        template_file=args.template,
        condition_file=args.condition,
        weights=args.weights,
        preferred_reactants=args.preferred_reactants,
        base_dir=ENGINE_DATA_DIR,
        beam_width=args.beam_width,
        per_step_top_k=args.per_step_top_k,
    )
    return result


def main():
    parser = argparse.ArgumentParser(
        description="SimpRetro Retrosynthesis CLI — single-step and multi-step"
    )
    parser.add_argument("-s", "--smiles", required=True, help="Target molecule SMILES string")
    parser.add_argument("--steps", type=int, default=1, help="Number of retrosynthesis steps (default: 1, single-step)")
    parser.add_argument("-db", "--database", default=DEFAULT_DATABASE, help=f"Stock database name (default: {DEFAULT_DATABASE})")
    parser.add_argument("-tpl", "--template", default=DEFAULT_TEMPLATE_FILE, help="Reaction template file")
    parser.add_argument("-cond", "--condition", default=DEFAULT_CONDITION_FILE, help="Reaction condition file")
    parser.add_argument("-w", "--weights", type=float, nargs=4, default=list(DEFAULT_WEIGHTS), help="Four scoring weights")
    parser.add_argument("-k", "--top-k", type=int, default=5, help="Number of routes to return for single-step (default: 5)")
    parser.add_argument("-pr", "--preferred-reactants", nargs="*", default=[], help="Preferred reactant SMILES to prioritize")
    parser.add_argument("-o", "--output", default="retro_result.json", help="Output JSON file path")
    parser.add_argument("--beam-width", type=int, default=5, help="Beam width for multi-step search (default: 5)")
    parser.add_argument("--per-step-top-k", type=int, default=5, help="Top-K per step for multi-step search (default: 5)")
    args = parser.parse_args()

    print("=" * 50)
    print(f"Target: {args.smiles}")
    print(f"Mode:   {'multi-step' if args.steps >= 2 else 'single-step'} ({args.steps} step(s))")
    print("=" * 50)

    if args.steps >= 2:
        result = run_multi_step(args)
    else:
        result = run_single_step(args)

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    # Print summary
    data = result.get("data", {})
    mode = data.get("mode", "single_step")
    if mode == "single_step":
        routes = data.get("retrosynthesis_routes", [])
        print(f"\n[Done] {len(routes)} route(s) returned.")
        for r in routes[:3]:
            rank = r["route_rank"]
            score = r["score"]
            cond = ", ".join(r.get("reaction_condition", []))
            n_reactants = len(r.get("reactants", []))
            print(f"  Route {rank}: score={score:.4f}, {n_reactants} reactant(s), conditions: {cond}")
    else:
        route = data.get("recommended_route")
        if route:
            print(f"\n[Done] Best route: {route['actual_steps']} step(s), score={route['route_score']:.4f}")
            for step in route.get("steps", []):
                n = step["step_number"]
                smi = step["expanded_smiles"][:30]
                print(f"  Step {n}: expand {smi}..., score={step['step_score']:.4f}")
        else:
            print("\n[Done] No viable multi-step route found.")

    print(f"\nResult saved to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
