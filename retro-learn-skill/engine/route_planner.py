"""SimpRetro multi-step route planner — beam search."""

import os
import json
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .retro_engine import (
    run_retrosynthesis, get_engine_state, canonical_smiles,
    DEFAULT_DATABASE, DEFAULT_TEMPLATE_FILE, DEFAULT_CONDITION_FILE, DEFAULT_WEIGHTS,
)


def _count_atoms(smiles):
    """Get heavy atom count from SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumAtoms() if mol else 0


def _route_key(route_state):
    """Sorting key: prefer completed routes, then by score."""
    return (not route_state.get('completed', False), -route_state.get('score', 0))


def _all_leaf_reactants(route_state):
    """Get all leaf reactant SMILES from a route state."""
    return list(route_state.get('leaf_reactants', {}).keys())


def _check_preferred(target_smiles_list, preferred_set):
    """Check if any target matches preferred reactants."""
    if not preferred_set:
        return False
    for s in target_smiles_list:
        canonical = canonical_smiles(s)
        if canonical and canonical in preferred_set:
            return True
    return False


def plan_retrosynthesis(smiles, steps=2, database_name=None,
                         template_file=None, condition_file=None,
                         weights=None, preferred_reactants=None,
                         base_dir=None, beam_width=5, per_step_top_k=5):
    """Multi-step retrosynthesis via beam search."""
    if database_name is None:
        database_name = DEFAULT_DATABASE
    if template_file is None:
        template_file = DEFAULT_TEMPLATE_FILE
    if condition_file is None:
        condition_file = DEFAULT_CONDITION_FILE
    if weights is None:
        weights = list(DEFAULT_WEIGHTS)
    if preferred_reactants is None:
        preferred_reactants = []
    if base_dir is None:
        base_dir = os.path.dirname(os.path.abspath(__file__))

    preferred_set = set()
    for smi in preferred_reactants:
        cs = canonical_smiles(smi)
        if cs:
            preferred_set.add(cs)

    # Step 1: initial single-step
    result = run_retrosynthesis(
        smiles=smiles, database_name=database_name,
        template_file=template_file, condition_file=condition_file,
        weights=weights, preferred_reactants=preferred_reactants,
        top_k=per_step_top_k, base_dir=base_dir, show_progress=False,
    )

    if result.get("status") != "success":
        return {
            "status": "error",
            "message": f"Initial step failed: {result.get('message', 'unknown')}",
            "data": None
        }

    routes = result["data"].get("retrosynthesis_routes", [])
    if not routes:
        return {
            "status": "success",
            "message": "No routes found.",
            "data": {
                "target_molecule": result["data"]["target_molecule"],
                "mode": "multi_step",
                "recommended_route": None,
                "all_routes": []
            }
        }

    # Build initial beam states
    engine = get_engine_state()
    beam = []
    for route in routes:
        leaf_reactants = {}
        for r in route["reactants"]:
            leaf_reactants[r["smiles"]] = {
                "smiles": r["smiles"],
                "in_stock": r["in_stock"],
                "molecular_weight": r.get("molecular_weight")
            }

        leaf_list = list(leaf_reactants.keys())
        all_in_stock = all(leaf_reactants[s]["in_stock"] for s in leaf_list)
        all_preferred = _check_preferred(leaf_list, preferred_set) if preferred_set else False

        step_record = {
            "step_number": 1,
            "target_smiles": smiles,
            "expanded_smiles": leaf_list[0] if len(leaf_list) == 1 else ".".join(leaf_list),
            "reaction_template": route.get("reaction_template", ""),
            "reaction_condition": route.get("reaction_condition", []),
            "step_score": route.get("score", 0),
        }

        state = {
            "score": route.get("score", 0),
            "leaf_reactants": leaf_reactants,
            "steps_history": [step_record],
            "completed": all_in_stock or (preferred_set and all_preferred),
        }
        beam.append(state)

    # Sort beam
    beam.sort(key=_route_key)

    if steps <= 1:
        return _build_multi_result(smiles, result["data"]["target_molecule"], beam)

    # Iterative expansion
    for step_idx in range(2, steps + 1):
        new_beam = []
        for state in beam[:beam_width]:
            if state["completed"]:
                new_beam.append(state)
                continue

            leaf = state["leaf_reactants"]
            # Pick the non-stock leaf with most atoms
            non_stock = [(s, d) for s, d in leaf.items() if not d["in_stock"]]
            if not non_stock:
                state["completed"] = True
                new_beam.append(state)
                continue

            target_leaf = max(non_stock, key=lambda x: _count_atoms(x[0]))[0]

            # Run single-step on this leaf
            sub = run_retrosynthesis(
                smiles=target_leaf, database_name=database_name,
                template_file=template_file, condition_file=condition_file,
                weights=weights, preferred_reactants=preferred_reactants,
                top_k=per_step_top_k, base_dir=base_dir, show_progress=False,
            )

            if sub.get("status") != "success":
                # Leaf cannot be expanded — keep this state as a finished shorter route
                state["completed"] = True
                new_beam.append(state)
                continue

            sub_routes = sub["data"].get("retrosynthesis_routes", [])

            if not sub_routes:
                # No templates matched this leaf — keep as finished shorter route
                state["completed"] = True
                new_beam.append(state)
                continue

            for sr in sub_routes[:per_step_top_k]:
                new_leaf = dict(leaf)
                del new_leaf[target_leaf]
                for r in sr["reactants"]:
                    new_leaf[r["smiles"]] = {
                        "smiles": r["smiles"],
                        "in_stock": r["in_stock"],
                        "molecular_weight": r.get("molecular_weight")
                    }

                leaf_list = list(new_leaf.keys())
                all_preferred = _check_preferred(leaf_list, preferred_set) if preferred_set else False

                new_score = state["score"] + sr.get("score", 0)
                new_step = {
                    "step_number": step_idx,
                    "target_smiles": target_leaf,
                    "expanded_smiles": sr["reactants"][0]["smiles"] if len(sr["reactants"]) == 1
                        else ".".join(r["smiles"] for r in sr["reactants"]),
                    "reaction_template": sr.get("reaction_template", ""),
                    "reaction_condition": sr.get("reaction_condition", []),
                    "step_score": sr.get("score", 0),
                }

                new_beam.append({
                    "score": new_score,
                    "leaf_reactants": new_leaf,
                    "steps_history": state["steps_history"] + [new_step],
                    "completed": all(d["in_stock"] for d in new_leaf.values())
                                  or (preferred_set and all_preferred),
                })

        # Deduplicate: keep best state per leaf set
        dedup = {}
        for s in new_beam:
            key = tuple(sorted(s["leaf_reactants"].keys()))
            if key not in dedup or s["score"] > dedup[key]["score"]:
                dedup[key] = s
        beam = sorted(dedup.values(), key=_route_key)

    return _build_multi_result(smiles, result["data"]["target_molecule"], beam)


def _build_multi_result(target_smiles, target_mol_info, beam):
    """Build the final multi-step result structure."""
    if not beam:
        return {
            "status": "success",
            "message": "No viable multi-step route found.",
            "data": {
                "target_molecule": target_mol_info,
                "mode": "multi_step",
                "recommended_route": None,
                "all_routes": []
            }
        }

    # Check for preferred reactant match
    best = beam[0]
    final_leaves = [
        {"smiles": s, "in_stock": d["in_stock"], "molecular_weight": d["molecular_weight"]}
        for s, d in best["leaf_reactants"].items()
    ]

    recommended = {
        "actual_steps": len(best["steps_history"]),
        "route_score": round(best["score"], 4),
        "leaf_reactants": final_leaves,
        "steps": best["steps_history"],
    }

    all_routes = []
    for i, state in enumerate(beam[:10]):
        all_routes.append({
            "route_rank": i + 1,
            "route_score": round(state["score"], 4),
            "steps": len(state["steps_history"]),
            "leaf_reactants": [
                {"smiles": s, "in_stock": d["in_stock"]}
                for s, d in state["leaf_reactants"].items()
            ],
            "steps_history": state["steps_history"],
        })

    return {
        "status": "success",
        "message": f"Found {len(all_routes)} multi-step route(s).",
        "data": {
            "target_molecule": target_mol_info,
            "mode": "multi_step",
            "recommended_route": recommended,
            "all_routes": all_routes,
        }
    }
