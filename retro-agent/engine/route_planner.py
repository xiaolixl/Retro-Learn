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

    # Step 1: initial single-step — use wider top_k to capture all viable routes
    # (including low-score but chemically valid paths like Grignard coupling).
    # Subsequent expansion rounds use per_step_top_k for narrower branching.
    initial_top_k = max(per_step_top_k, 5)
    result = run_retrosynthesis(
        smiles=smiles, database_name=database_name,
        template_file=template_file, condition_file=condition_file,
        weights=weights, preferred_reactants=preferred_reactants,
        top_k=initial_top_k, base_dir=base_dir, show_progress=False,
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
        return _build_multi_result(smiles, result["data"]["target_molecule"], beam, requested_steps=steps)

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

            # Cycle detection: collect all molecules previously expanded as targets
            # in this route's history. If target_leaf was already expanded, skip.
            expanded_targets = set()
            for step in state["steps_history"]:
                ct = canonical_smiles(step.get("target_smiles", ""))
                if ct:
                    expanded_targets.add(ct)
            target_canon = canonical_smiles(target_leaf)
            if target_canon and target_canon in expanded_targets:
                # Cycle: this molecule was already retrosynthesized earlier
                state["completed"] = True
                new_beam.append(state)
                continue

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
                # Cycle detection: skip if any new intermediate re-introduces
                # a molecule that was already a target in this route's history
                cycle = False
                for r in sr["reactants"]:
                    rc = canonical_smiles(r["smiles"])
                    if rc and rc in expanded_targets:
                        cycle = True
                        break
                if cycle:
                    continue

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

    return _build_multi_result(smiles, result["data"]["target_molecule"], beam, requested_steps=steps)


def _build_multi_result(target_smiles, target_mol_info, beam, requested_steps=3):
    """Build multi-step results, dropping routes that are prefixes of longer routes.

    If route A's expanded_smiles signature is a strict prefix of route B's,
    route A is discarded — it carries no additional synthetic information
    beyond what B already provides.
    """
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

    def route_signature(history):
        """Extract expanded_smiles at each step as a route fingerprint."""
        return [step.get("expanded_smiles", "") for step in history]

    def is_prefix(short_sig, long_sig):
        """True if short_sig is a strict prefix of long_sig."""
        return len(short_sig) < len(long_sig) and short_sig == long_sig[:len(short_sig)]

    # --- Step 1: deduplicate by signature, keep best score per sig ---
    by_signature = {}
    for original_rank, state in enumerate(beam, start=1):
        history = state["steps_history"]
        raw_score = round(sum(step.get("step_score", 0) for step in history), 4)
        sig = tuple(route_signature(history))
        candidate = {
            "original_rank": original_rank,
            "route_score": raw_score,
            "steps": len(history),
            "leaf_reactants": list(state.get("leaf_reactants", {}).values()),
            "steps_history": history,
            "signature": list(sig),
        }
        existing = by_signature.get(sig)
        if existing is None or (
            raw_score,
            len(history),
            -original_rank,
        ) > (
            existing["route_score"],
            existing["steps"],
            -existing["original_rank"],
        ):
            by_signature[sig] = candidate

    # --- Step 2: prefix filtering — drop short routes that are
    #        just the first few steps of a longer route ---
    candidates = list(by_signature.values())
    kept = []
    for route in candidates:
        sig = route["signature"]
        if any(is_prefix(sig, other["signature"]) for other in candidates):
            continue
        kept.append(route)

    # --- Step 3: sort by score desc, steps desc, original rank asc ---
    kept.sort(key=lambda route: (
        -route["route_score"],
        -route["steps"],
        route["original_rank"],
    ))

    best = kept[0]
    recommended = {
        "actual_steps": best["steps"],
        "route_score": best["route_score"],
        "leaf_reactants": best["leaf_reactants"],
        "steps": best["steps_history"],
    }

    all_routes = []
    for rank, route in enumerate(kept, start=1):
        all_routes.append({
            "route_rank": rank,
            "route_score": route["route_score"],
            "steps": route["steps"],
            "leaf_reactants": route["leaf_reactants"],
            "steps_history": route["steps_history"],
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
