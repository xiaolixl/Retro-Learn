"""SimpRetro retrosynthesis engine — template matching + scoring."""

import os
import re
import json
import numpy as np
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from rdchiral.main import rdchiralRun
from rdchiral.initialization import rdchiralReaction, rdchiralReactants

# ── Defaults ──────────────────────────────────────────────
DEFAULT_TEMPLATE_FILE = "reaction_template.json"
DEFAULT_CONDITION_FILE = "template_condition.json"
DEFAULT_DATABASE = "emol_under_0_carbons"
DEFAULT_WEIGHTS = (0.1, 0.2, 0.5, 0.0)

# ── Scoring functions ─────────────────────────────────────
def CDScore(p_mol, r_mols):
    p_atom_count = p_mol.GetNumAtoms()
    n_r_mols = len(r_mols)
    if len(r_mols) == 1:
        return 0
    r_atom_count = [len([int(num[1:]) for num in re.findall(r':\d+', r_mol) if int(num[1:]) < 900]) for r_mol in r_mols]
    main_r = r_mols[np.argmax(r_atom_count)]
    if len(Chem.MolFromSmiles(main_r).GetAtoms()) >= p_atom_count:
        return 0
    MAE = 1 / n_r_mols * sum([abs(p_atom_count / n_r_mols - r_atom_count[i]) for i in range(n_r_mols)])
    return 1 / (1 + MAE) * p_atom_count

def canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))

def ASScore(p_mol, r_mol_dict, in_stock):
    p_atom_count = p_mol.GetNumAtoms()
    r_mols = list(r_mol_dict.keys())
    r_atom_count = [len([int(num[1:]) for num in re.findall(r':\d+', r_mol) if int(num[1:]) < 900]) for r_mol in r_mols]
    main_r = r_mols[np.argmax(r_atom_count)]
    asscore = 0
    for k, v in r_mol_dict.items():
        if v in in_stock:
            add = len([int(num[1:]) for num in re.findall(r':\d+', k) if int(num[1:]) < 900])
            if len(Chem.MolFromSmiles(main_r).GetAtoms()) < p_atom_count:
                asscore += add
            else:
                asscore += add if add > 2 else 0
        if ('Mg' in v or 'Li' in v or 'Zn' in v) and v not in in_stock:
            asscore -= 5
    return asscore

def RDScore(p_mol, r_mols):
    p_ring_count = p_mol.GetRingInfo().NumRings()
    r_rings_s = [r_mol.GetRingInfo().AtomRings() for r_mol in r_mols]
    r_ring_count = 0
    for r_rings, r_mol in zip(r_rings_s, r_mols):
        for r_ring in r_rings:
            mapnums = [r_mol.GetAtomWithIdx(i).GetAtomMapNum() for i in r_ring]
            symbols = [r_mol.GetAtomWithIdx(i).GetSymbol() for i in r_ring]
            if 'B' in symbols or 'Si' in symbols:
                continue
            if min(mapnums) < 900:
                r_ring_count += 1
    if p_ring_count > r_ring_count:
        return 1
    else:
        return 0

# ── Engine ────────────────────────────────────────────────
class _EngineState:
    """Hold loaded engine data (replaces global vars)."""
    def __init__(self):
        self.in_stock = None
        self.templates_raw = None
        self.template_list = None
        self.tpl_condition = None

_engine = _EngineState()

def _load_engine_data(database_name, template_file, condition_file,
                       preferred_reactants, base_dir):
    """Load stock DB, templates, and conditions into the engine state."""
    # Database
    if database_name in ('emol_under_0_carbons', 'unrestricted'):
        in_stock = set()
    else:
        db_path = os.path.join(base_dir, database_name + '.txt')
        if not os.path.exists(db_path):
            raise FileNotFoundError(f"Stock database not found: {db_path}")
        with open(db_path, 'r', encoding='utf-8') as f:
            in_stock = set(f.read().splitlines())

    if preferred_reactants:
        in_stock.update(preferred_reactants)

    in_stock = set(canonical_smiles(s) for s in in_stock if canonical_smiles(s) is not None)

    # Templates
    tpl_path = os.path.join(base_dir, template_file)
    if not os.path.exists(tpl_path):
        raise FileNotFoundError(f"Template file not found: {tpl_path}")
    with open(tpl_path, 'r', encoding='utf-8') as f:
        templates_raw = json.load(f)

    template_list = []
    for tpl in tqdm(templates_raw, desc="Init templates"):
        template_list.append(rdchiralReaction(tpl))

    # Conditions
    cond_path = os.path.join(base_dir, condition_file)
    if not os.path.exists(cond_path):
        raise FileNotFoundError(f"Condition file not found: {cond_path}")
    with open(cond_path, 'r', encoding='utf-8') as f:
        tpl_condition = json.load(f)

    _engine.in_stock = in_stock
    _engine.templates_raw = templates_raw
    _engine.template_list = template_list
    _engine.tpl_condition = tpl_condition

def _format_single_step(target_smiles, results):
    """Format single-step results into the standard JSON structure."""
    target_mol = Chem.MolFromSmiles(target_smiles)
    target_mw = rdMolDescriptors.CalcExactMolWt(target_mol) if target_mol else None

    output = {
        "target_molecule": {
            "smiles": target_smiles,
            "molecular_weight": round(target_mw, 2) if target_mw else None,
            "in_stock": target_smiles in _engine.in_stock
        },
        "retrosynthesis_routes": []
    }

    if not results:
        return output

    for rank, (reactants_smiles_str, data) in enumerate(results.items(), start=1):
        reactants_list = []
        for r_smiles in reactants_smiles_str.split('.'):
            r_mol = Chem.MolFromSmiles(r_smiles)
            r_mw = rdMolDescriptors.CalcExactMolWt(r_mol) if r_mol else None
            reactants_list.append({
                "smiles": r_smiles,
                "molecular_weight": round(r_mw, 2) if r_mw else None,
                "in_stock": r_smiles in _engine.in_stock
            })

        route = {
            "route_rank": rank,
            "score": round(data['Score'], 4),
            "reaction_template": data['Template'],
            "reaction_condition": data['Condition'],
            "all_reactants_in_stock": all(r['in_stock'] for r in reactants_list),
            "reactants": reactants_list
        }
        output["retrosynthesis_routes"].append(route)

    return output

def run_retrosynthesis(smiles, database_name=None, template_file=None,
                        condition_file=None, weights=None,
                        preferred_reactants=None, top_k=5,
                        base_dir=None, show_progress=True):
    """Run single-step retrosynthesis on a target molecule."""
    if database_name is None:
        database_name = DEFAULT_DATABASE
    if template_file is None:
        template_file = DEFAULT_TEMPLATE_FILE
    if condition_file is None:
        condition_file = DEFAULT_CONDITION_FILE
    if weights is None:
        weights = DEFAULT_WEIGHTS
    if preferred_reactants is None:
        preferred_reactants = []
    if base_dir is None:
        base_dir = os.path.dirname(os.path.abspath(__file__))

    # Validate
    if not Chem.MolFromSmiles(smiles):
        return {"status": "error", "message": f"Invalid SMILES: {smiles}", "data": None}

    # Load
    _load_engine_data(database_name, template_file, condition_file,
                       preferred_reactants, base_dir)

    w1, w2, w3, w4 = weights
    if show_progress:
        print(f"[Config] Weights: w1={w1}, w2={w2}, w3={w3}, w4={w4}")

    results = {}
    result_set = set()
    p_mol_rdchiral = rdchiralReactants(smiles)
    p_mol = Chem.MolFromSmiles(smiles)

    for idx, (template, template_raw) in enumerate(zip(_engine.template_list,
                                                         _engine.templates_raw)):
        mapped_curr_results = rdchiralRun(template, p_mol_rdchiral, keep_mapnums=True)
        for r in mapped_curr_results:
            canonical_r = canonical_smiles(r)
            canonical_r_dict = {r_: canonical_smiles(r_) for r_ in r.split('.')}

            if canonical_r in result_set:
                continue
            result_set.add(canonical_r)

            r_mols = [Chem.MolFromSmiles(r_) for r_ in r.split('.')]
            rdscore = RDScore(p_mol, r_mols)

            score = 1 * (w1 * CDScore(p_mol, r.split('.')) +
                         w2 * ASScore(p_mol, canonical_r_dict, _engine.in_stock) +
                         w3 * rdscore +
                         w4 * 1 / max(len(mapped_curr_results), 1))

            results[canonical_r] = {
                'Score': score,
                'Template': template_raw,
                'Template_id': idx,
                'Condition': _engine.tpl_condition.get(template_raw, None)
            }

            results = dict(sorted(results.items(),
                                   key=lambda x: x[1]['Score'], reverse=True))

    # Keep top-k
    results = dict(list(results.items())[:top_k])

    formatted = _format_single_step(smiles, results)
    return {
        "status": "success",
        "message": f"Retrosynthesis completed. {len(formatted['retrosynthesis_routes'])} route(s).",
        "data": {**formatted, "mode": "single_step"}
    }

# ── Helper: expose engine state for route_planner ─────────
def get_engine_state():
    """Return the loaded engine state (in_stock, templates, conditions)."""
    return _engine
